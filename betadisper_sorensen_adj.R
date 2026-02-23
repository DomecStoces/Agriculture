# 1. Data Prep & Partitioning Beta Diversity
# Convert community matrix to presence/absence
comm_pa <- ifelse(comm_matrix > 0, 1, 0)

# Partition beta diversity (Sørensen framework)
beta_part_sor <- beta.pair(comm_pa, index.family = "sorensen")

# Extract components
dist_turnover <- beta_part_sor$beta.sim   # Turnover (Simpson)
dist_nested   <- beta_part_sor$beta.sne   # Nestedness-resultant
dist_total    <- beta_part_sor$beta.sor   # Total Sørensen

# 2. PERMDISP (Homogeneity of Dispersions)
# Test homogeneity using Lingoes correction (add = TRUE) and bias adjustment
disp_turnover <- betadisper(dist_turnover, metadata$Treatment, add = TRUE, bias.adjust = TRUE)
disp_nested   <- betadisper(dist_nested,   metadata$Treatment, add = TRUE, bias.adjust = TRUE)
disp_total    <- betadisper(dist_total,    metadata$Treatment, add = TRUE, bias.adjust = TRUE)

# Global tests
perm_turnover <- permutest(disp_turnover, permutations = 999)
perm_nested   <- permutest(disp_nested,   permutations = 999)
perm_total    <- permutest(disp_total,    permutations = 999)

print(perm_turnover)
print(perm_nested)
print(perm_total)

tukey_turnover <- TukeyHSD(disp_turnover)
tukey_nested   <- TukeyHSD(disp_nested)
tukey_total    <- TukeyHSD(disp_total)

print(tukey_turnover)
print(tukey_nested)
print(tukey_total)
# 3. Prepare Data for boxplots

# Extract distances into a combined dataframe
df_turnover <- data.frame(Distance = disp_turnover$distances, Treatment = metadata$Treatment, Component = "Turnover")
df_nested   <- data.frame(Distance = disp_nested$distances,   Treatment = metadata$Treatment, Component = "Nestedness")
df_total    <- data.frame(Distance = disp_total$distances,    Treatment = metadata$Treatment, Component = "βdiversity")

plot_data_box <- bind_rows(df_nested, df_turnover, df_total) %>%
  mutate(Component = factor(Component, levels = c("Nestedness", "Turnover", "βdiversity")))

# Helper function to automatically extract Tukey letters into a dataframe
get_sig_letters <- function(disp_obj, comp_name) {
  tuk <- TukeyHSD(disp_obj)
  pvals <- tuk$group[, "p adj"]
  names(pvals) <- rownames(tuk$group)
  cld <- multcompLetters(pvals)$Letters
  
  data.frame(
    Treatment = names(cld),
    label = as.character(cld),
    Component = factor(comp_name, levels = c("Nestedness", "Turnover", "βdiversity")),
    stringsAsFactors = FALSE
  )
}

# Auto-generate annotation text dataframe for ggplot
ann_text <- bind_rows(
  get_sig_letters(disp_nested, "Nestedness"),
  get_sig_letters(disp_turnover, "Turnover"),
  get_sig_letters(disp_total, "βdiversity")
)

# Set Y-axis position for letters dynamically (slightly above the max distance)
ann_text$Distance <- max(plot_data_box$Distance)  

# 4. Plotting Boxplots
bw_colors <- c("EKOLOGIE" = "white", "KONVENCE" = "grey75", "REGENERACE" = "grey40")

d4 <- ggplot(plot_data_box, aes(x = Treatment, y = Distance, fill = Treatment)) +
  geom_boxplot(color = "black", outlier.shape = 16, outlier.size = 1.5, alpha = 0.9) +
  geom_text(data = ann_text, aes(label = label), vjust = 0, size = 5, color = "black") +
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
  labs(x = "Treatment", y = "Distance to group centroid")

print(d4)

# 5. PERMANOVA (Global and Pairwise)
set.seed(123)

# Global tests (Lingoes correction applied via add = TRUE)
global_nested   <- adonis2(dist_nested ~ Village + Treatment, data = metadata, by = "margin", permutations = 999, add = TRUE)
global_turnover <- adonis2(dist_turnover ~ Village + Treatment, data = metadata, by = "margin", permutations = 999, add = TRUE)
global_total    <- adonis2(dist_total ~ Village + Treatment, data = metadata, by = "margin", permutations = 999, add = TRUE)

print(global_nested)
print(global_turnover)
print(global_total)

# Custom function for Pairwise PERMANOVA (with Lingoes correction)
pairwise_adonis_custom <- function(dist_obj, metadata) {
  treatments <- unique(metadata$Treatment)
  pairs <- combn(treatments, 2, simplify = FALSE)
  
  results <- data.frame(Pair = character(), F_Value = numeric(), R2 = numeric(), P_Value = numeric(), stringsAsFactors = FALSE)
  
  set.seed(123)
  for (pair in pairs) {
    meta_sub <- metadata %>% filter(Treatment %in% pair)
    dist_sub <- as.dist(as.matrix(dist_obj)[meta_sub$Site_ID, meta_sub$Site_ID])
    
    # Run PERMANOVA on subset with Lingoes correction
    res <- adonis2(dist_sub ~ Village + Treatment, data = meta_sub, by = "margin", permutations = 999, add = TRUE)
    
    results <- rbind(results, data.frame(
      Pair = paste(pair[1], "vs", pair[2]),
      F_Value = round(res["Treatment", "F"], 3),
      R2 = round(res["Treatment", "R2"], 3),
      P_Value = res["Treatment", "Pr(>F)"]
    ))
  }
  
  # Adjust p-values using Benjamini-Hochberg
  results$P_Adj_BH <- round(p.adjust(results$P_Value, method = "BH"), 4)
  return(results)
}

# Pairwise comparisons outputs
pw_nested   <- pairwise_adonis_custom(dist_nested, metadata)
pw_turnover <- pairwise_adonis_custom(dist_turnover, metadata)
pw_total    <- pairwise_adonis_custom(dist_total, metadata)

print("Pairwise Nestedness:"); print(pw_nested)
print("Pairwise Turnover:");   print(pw_turnover)
print("Pairwise Total Beta:"); print(pw_total)
