# 1. Convert community matrix to Presence/Absence (Binary)
comm_pa <- ifelse(comm_matrix > 0, 1, 0)

# 2. Partition using the Sørensen framework
beta_part_sor <- beta.pair(comm_pa, index.family = "sorensen")

# Extract the three matrices:
dist_turnover <- beta_part_sor$beta.sim  # Simpson Turnover
dist_nested <- beta_part_sor$beta.sne    # Nestedness-resultant
dist_total <- beta_part_sor$beta.sor     # Total Sørensen dissimilarity

# 3. Calculate Distance to Centroids (Dispersion) for each matrix
disp_turnover <- betadisper(dist_turnover, metadata$Treatment)
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
ann_text <- data.frame(
  Treatment = c("EKOLOGIE", "KONVENCE", "REGENERACE", 
                "EKOLOGIE", "KONVENCE", "REGENERACE", 
                "EKOLOGIE", "KONVENCE", "REGENERACE"),
  Distance = 0.75,
  label = c("a", "a", "a",    # Nestedness (PERMDISP)
            "ab", "b", "a",   # Turnover (PERMDISP)
            "a", "b", "a"),   # Beta Diversity (PERMDISP)
  Component = factor(c("Nestedness", "Nestedness", "Nestedness", 
                       "Turnover", "Turnover", "Turnover", 
                       "βdiversity", "βdiversity", "βdiversity"),
                     levels = c("Nestedness", "Turnover", "βdiversity"))
)
d4<-ggplot(plot_data_box, aes(x = Treatment, y = Distance, fill = Treatment)) +
  geom_boxplot(color = "black", outlier.shape = 16, outlier.size = 1.5, alpha = 0.9) +
  # Add the labels
  geom_text(data = ann_text, aes(label = label), vjust = -0.5, size = 5, color = "black") +
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

### Pairwise comparisons of TOTAL BETA DIVERSITY ###
# Get all unique pairs of treatments
treatments <- unique(metadata$Treatment)
pairs <- combn(treatments, 2, simplify = FALSE)

# Create an empty dataframe to store results
pairwise_permanova_results <- data.frame(Pair = character(), P_Value = numeric(), stringsAsFactors = FALSE)
set.seed(123)
for (i in seq_along(pairs)) {
  pair <- pairs[[i]]
  # Subset metadata for just these two treatments
  meta_sub <- metadata %>% filter(Treatment %in% pair)
  
  # Subset the total Sørensen distance matrix for just these sites
  dist_sub <- as.dist(as.matrix(dist_total)[meta_sub$Site_ID, meta_sub$Site_ID])
  
  # Run adonis2 with your exact spatial strata
  res <- adonis2(dist_sub ~ Village + Crop + Treatment, 
                 data = meta_sub, 
                 by = "margin", 
                 strata = meta_sub$Locality, 
                 permutations = 999)
  
  # Save the p-value for the Treatment variable
  pairwise_permanova_results <- rbind(pairwise_permanova_results, data.frame(
    Pair = paste(pair[1], "vs", pair[2]),
    P_Value = res["Treatment", "Pr(>F)"]
  ))
}
# --- 3. APPLY FDR CORRECTION ---
pairwise_permanova_results$P_Adj_FDR <- p.adjust(pairwise_permanova_results$P_Value, method = "fdr")
print("--- PAIRWISE PERMANOVA: Composition Differences ---")
print(pairwise_permanova_results)

### Pairwise comparisons of TURNOVER and NESTEDNESS DISPERSION ###
TukeyHSD(disp_turnover)
TukeyHSD(disp_nested)
TukeyHSD(disp_total)
