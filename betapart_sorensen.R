# 1. Convert community matrix to Presence/Absence (Binary)
comm_pa <- ifelse(comm_matrix > 0, 1, 0)

# 2. Partition using the Sørensen framework (as per your paper)
beta_part_sor <- beta.pair(comm_pa, index.family = "sorensen")

# Extract the three matrices:
dist_turnover <- beta_part_sor$beta.sim  # Simpson Turnover
dist_nested <- beta_part_sor$beta.sne    # Nestedness-resultant
dist_total <- beta_part_sor$beta.sor     # Total Sørensen dissimilarity

# 3. Calculate Distance to Centroids (Dispersion) for each matrix
disp_turnover <- betadisper(sqrt(dist_turnover), metadata$Treatment)
disp_nested <- betadisper(dist_nested, metadata$Treatment)
disp_total <- betadisper(dist_total, metadata$Treatment)

# Extract distances and create dataframes
df_turnover <- data.frame(Distance = disp_turnover$distances, 
                          Treatment = metadata$Treatment, 
                          Component = "Turnover")

df_nested <- data.frame(Distance = disp_nested$distances, 
                        Treatment = metadata$Treatment, 
                        Component = "Nestedness")

df_total <- data.frame(Distance = disp_total$distances, 
                       Treatment = metadata$Treatment, 
                       Component = "βdiversity")

# Combine them all together
plot_data_box <- bind_rows(df_nested, df_turnover, df_total)

# Lock in the factor levels so they plot in the exact order of your example image!
plot_data_box$Component <- factor(plot_data_box$Component, 
                                  levels = c("Nestedness", "Turnover", "βdiversity"))
bw_colors <- c("EKOLOGIE" = "white", 
               "KONVENCE" = "grey75", 
               "REGENERACE" = "grey40")
d4<-ggplot(plot_data_box, aes(x = Treatment, y = Distance, fill = Treatment)) +
  geom_boxplot(color = "black", outlier.shape = 16, outlier.size = 1.5, alpha = 0.9) +
  facet_wrap(~Component, scales = "free_y") +
  scale_fill_manual(values = bw_colors) +
  theme_bw(base_size = 15) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "plain", size = 14, color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    legend.position = "none" 
  ) +
  labs(x = "Treatment",
       y = "Distance to group centroid")
d4
# --- A. TOTAL BETA DIVERSITY (PERMANOVA) ---
set.seed(123)
permanova_total <- adonis2(dist_total ~ Village + Crop + Treatment, 
                           data = metadata, 
                           by = "margin", 
                           strata = metadata$Locality,
                           permutations = 999)
print("--- PERMANOVA: TOTAL BETA DIVERSITY ---")
print(permanova_total)

# --- B. PERMUTATION TESTS FOR COMPONENTS (Dispersion) ---
# This tests if the boxplots you just drew are significantly different in height
print("--- PERMUTATION TEST: TURNOVER DISPERSION ---")
permutest(disp_turnover, permutations = 999)

print("--- PERMUTATION TEST: NESTEDNESS DISPERSION ---")
permutest(disp_nested, permutations = 999)

# To get the exact pairwise differences between treatments (e.g., EKO vs KON):
TukeyHSD(disp_turnover)
TukeyHSD(disp_nested)
