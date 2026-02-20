### Bray-Curtis Dissimilirity Matrix ###
df_nmds <- long_format %>%
  mutate(Site_ID = paste(Village, Locality, Trap, sep = "_"))

comm_wide <- df_nmds %>%
  group_by(Site_ID, Village, Locality, Treatment, Crop, Species) %>%
  summarise(Count = sum(Count, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Species, values_from = Count, values_fill = 0)

comm_matrix <- as.matrix(comm_wide %>% select(-Site_ID, -Village, -Locality, -Treatment, -Crop))
rownames(comm_matrix) <- comm_wide$Site_ID

# Save metadata
metadata <- comm_wide %>% select(Site_ID, Village, Locality, Treatment, Crop)

# PERMANOVA
set.seed(123)
permanova_res <- adonis2(comm_matrix ~ Village + Crop + Treatment, 
                         data = metadata, 
                         method = "bray", 
                         by = "margin",               
                         strata = metadata$Locality,
                         permutations = 999)
print(permanova_res)

set.seed(123)
nmds_res <- metaMDS(comm_matrix, distance = "bray", k = 2, trymax = 100)

print(paste("NMDS Stress:", round(nmds_res$stress, 3)))

# Visualization of NMDS 
nmds_coords <- as.data.frame(scores(nmds_res, display = "sites"))
nmds_coords <- cbind(nmds_coords, metadata)
ggplot(nmds_coords, aes(x = NMDS1, y = NMDS2, color = Treatment, shape = Crop)) +
  # Add confidence ellipses for the Treatments
  stat_ellipse(aes(group = Treatment), type = "t", linetype = 2, linewidth = 1) +
  # Add the site points
  geom_point(size = 3.5, alpha = 0.8) +
  # Apply nice theme and labels
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

# INDICATOR SPECIES ANALYSIS
indval_res <- multipatt(comm_matrix, 
                        metadata$Treatment, 
                        func = "IndVal.g", 
                        control = how(nperm = 999))
summary(indval_res)
