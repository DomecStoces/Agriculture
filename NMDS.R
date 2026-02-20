### Bray-Curtis Dissimilirity Matrix ###
### 1. Data Prep ###
df_nmds <- long_format %>%
  mutate(Site_ID = paste(Village, Locality, Trap, sep = "_"))

comm_wide <- df_nmds %>%
  group_by(Site_ID, Village, Locality, Treatment, Crop, Species) %>%
  summarise(Count = sum(Count, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Species, values_from = Count, values_fill = 0)

comm_matrix <- as.matrix(comm_wide %>% select(-Site_ID, -Village, -Locality, -Treatment, -Crop))
rownames(comm_matrix) <- comm_wide$Site_ID
metadata <- comm_wide %>% select(Site_ID, Village, Locality, Treatment, Crop)

# --- THE FIX: Square Root Transformation ---
# This is standard practice before Bray-Curtis for over-abundant count data
comm_sqrt <- sqrt(comm_matrix)


### 2. PERMDISP ###
# Use the sqrt-transformed data + Bray
bray_dist <- vegdist(comm_sqrt, method = "bray")

dispersion_treatment <- betadisper(bray_dist, metadata$Treatment)
permutest(dispersion_treatment, permutations = 999)
boxplot(dispersion_treatment)
tukey_dispersion <- TukeyHSD(dispersion_treatment)


### 3. PERMANOVA ###
set.seed(123)
# Use comm_sqrt + Bray
permanova_res <- adonis2(comm_sqrt ~ Village + Crop + Treatment, 
                         data = metadata, 
                         method = "bray", 
                         by = "margin",               
                         strata = metadata$Locality,
                         permutations = 999)
print(permanova_res)


### 4. NMDS ###
set.seed(123)
# Use comm_sqrt + Bray, and keep autotransform = FALSE so we have manual control
nmds_res <- metaMDS(comm_sqrt, distance = "bray", k = 2, trymax = 100, autotransform = FALSE) 

print(paste("NMDS Stress:", round(nmds_res$stress, 3)))


### 5. GGPLOT2 Visualization ###
nmds_coords <- as.data.frame(vegan::scores(nmds_res, display = "sites"))
nmds_coords <- cbind(nmds_coords, metadata)

d3<-ggplot(nmds_coords, aes(x = NMDS1, y = NMDS2, color = Treatment, shape = Crop)) +
  geom_point(size = 3.5, alpha = 0.8) +
  stat_ellipse(aes(group = Treatment), type = "t", linetype = 2, linewidth = 1, show.legend = FALSE) +
  theme_classic() +
  labs(subtitle = paste("Stress =", round(nmds_res$stress, 3)),
       x = "NMDS 1",
       y = "NMDS 2") +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11)
  )
d3

### 6. INDICATOR SPECIES ANALYSIS ###
table(metadata$Treatment)
# raw comm_matrix counts
indval_res <- multipatt(comm_matrix, 
                        metadata$Treatment, 
                        func = "IndVal.g", 
                        max.order = 1,
                        control = how(nperm = 999))
summary(indval_res, alpha = 0.05, indvalcomp = TRUE)