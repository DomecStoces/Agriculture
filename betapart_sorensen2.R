# 1. Convert community matrix to presence/absence
comm_pa <- ifelse(comm_matrix > 0, 1, 0)

# 2. Partition beta diversity (Sørensen)
beta_part_sor <- beta.pair(comm_pa, index.family = "sorensen")

# Extract components
dist_turnover <- beta_part_sor$beta.sim   # turnover
dist_nested <- beta_part_sor$beta.sne     # nestedness-resultant
dist_total <- beta_part_sor$beta.sor      # total Sørensen

# Square-root transform to make distances Euclidean-compatible for betadisper
dist_turnover_sqrt <- sqrt(dist_turnover)
dist_nested_sqrt   <- sqrt(dist_nested)
dist_total_sqrt    <- sqrt(dist_total)

# 3. Test homogeneity of dispersions (PERMDISP)
disp_turnover <- betadisper(dist_turnover_sqrt, metadata$Treatment)
disp_nested   <- betadisper(dist_nested_sqrt, metadata$Treatment)
disp_total    <- betadisper(dist_total_sqrt, metadata$Treatment)

# Global tests
perm_turnover <- permutest(disp_turnover, permutations = 999)
perm_nested   <- permutest(disp_nested, permutations = 999)
perm_total    <- permutest(disp_total, permutations = 999)

print(perm_turnover)
print(perm_nested)
print(perm_total)

TukeyHSD(disp_turnover)
TukeyHSD(disp_total)
# 4. Extract distances and prepare plot data
df_turnover <- data.frame(Distance = disp_turnover$distances,
                          Treatment = metadata$Treatment,
                          Component = "Turnover")
df_nested <- data.frame(Distance = disp_nested$distances,
                        Treatment = metadata$Treatment,
                        Component = "Nestedness")
df_total <- data.frame(Distance = disp_total$distances,
                       Treatment = metadata$Treatment,
                       Component = "βdiversity")

plot_data_box <- bind_rows(df_nested, df_turnover, df_total)
plot_data_box$Component <- factor(plot_data_box$Component, 
                                  levels = c("Nestedness", "Turnover", "βdiversity"))

# 5. Generate significance letters automatically (Tukey HSD)
get_sig_letters <- function(disp_obj) {
  # 1. Run Tukey
  tuk <- TukeyHSD(disp_obj)
  
  # 2. Extract p-values using $group (the default name betadisper gives it)
  pvals <- tuk$group[, "p adj"]
  names(pvals) <- rownames(tuk$group)
  
  # 3. Get the letters
  library(multcompView)
  cld <- multcompLetters(pvals)$Letters
  
  # 4. multcompLetters already outputs a named vector (e.g., EKOLOGIE = "a")
  # We just need to convert it directly into a dataframe!
  df <- data.frame(Treatment = names(cld),
                   label = as.character(cld),
                   stringsAsFactors = FALSE)
  
  return(df)
}

# 6. Plot boxplots with significance letters
ann_text <- data.frame(
  Treatment = rep(c("EKOLOGIE", "KONVENCE", "REGENERACE"), 3),
  Component = factor(
    rep(c("Nestedness", "Turnover", "βdiversity"), each = 3), 
    levels = c("Nestedness", "Turnover", "βdiversity")
  ),
  label = c(
    "a", "a", "a",  # Nestedness letters
    "a", "b", "c",  # Turnover letters
    "a", "b", "c"   # βdiversity letters
  ),
  Distance = 0.75   # Y-axis position for the letters. Adjust if needed!
)
bw_colors <- c("EKOLOGIE" = "white", 
               "KONVENCE" = "grey75", 
               "REGENERACE" = "grey40")

d4<-ggplot(plot_data_box, aes(x = Treatment, y = Distance, fill = Treatment)) +
  geom_boxplot(color = "black", outlier.shape = 16, outlier.size = 1.5, alpha = 0.9) +
  geom_text(data = ann_text, aes(label = label), vjust = -0.5, size = 5, color = "black") +
  facet_wrap(~Component, scales = "free_y") +
  scale_fill_manual(values = bw_colors) +
  theme_bw(base_size = 15) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face = "plain", size = 14, color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.title = element_text(color = "black"),
        legend.position = "none") +
  labs(x = "Treatment", y = "Distance to group centroid")
d4

# 7. PERMANOVA on total beta diversity
set.seed(123)
permanova_total <- adonis2(dist_total ~ Village + Treatment, 
                           data = metadata, 
                           by = "margin",
                           permutations = 999)
print(permanova_total)

treatments <- unique(metadata$Treatment)
pairs <- combn(treatments, 2, simplify = FALSE)

pairwise_permanova_results <- data.frame(Pair = character(),
                                         P_Value = numeric(),
                                         stringsAsFactors = FALSE)
set.seed(123)
for (i in seq_along(pairs)) {
  pair <- pairs[[i]]
  meta_sub <- metadata %>% filter(Treatment %in% pair)
  
  dist_sub <- as.dist(as.matrix(dist_total)[meta_sub$Site_ID, meta_sub$Site_ID])
  
  res <- adonis2(dist_sub ~ Village + Treatment,
                 data = meta_sub,
                 by = "margin",
                 permutations = 999)
  
  pairwise_permanova_results <- rbind(pairwise_permanova_results, 
                                      data.frame(Pair = paste(pair[1], "vs", pair[2]),
                                                 P_Value = res["Treatment", "Pr(>F)"]))
}

# 8. PERMANOVA on turnover and nestedness components
pairwise_adonis_custom <- function(dist_obj, metadata) {
  treatments <- unique(metadata$Treatment)
  pairs <- combn(treatments, 2, simplify = FALSE)
  results <- data.frame(Pair = character(), 
                        F_Value = numeric(),
                        R2 = numeric(),
                        P_Value = numeric(), 
                        stringsAsFactors = FALSE)
  # Loop through each pair
  set.seed(123)
  for (i in seq_along(pairs)) {
    pair <- pairs[[i]]
    meta_sub <- metadata %>% filter(Treatment %in% pair)
    dist_matrix <- as.matrix(dist_obj)
    dist_sub <- as.dist(dist_matrix[meta_sub$Site_ID, meta_sub$Site_ID])
    res <- adonis2(dist_sub ~ Village + Treatment, 
                   data = meta_sub, 
                   by = "margin", 
                   permutations = 999)
    # Extract stats specifically for the Treatment
    results <- rbind(results, data.frame(
      Pair = paste(pair[1], "vs", pair[2]),
      F_Value = round(res["Treatment", "F"], 3),
      R2 = round(res["Treatment", "R2"], 3),
      P_Value = res["Treatment", "Pr(>F)"]
    ))
  }
  results$P_Adj_BH <- round(p.adjust(results$P_Value, method = "BH"), 4)
  return(results)
}

# 9a) Result outputs for nestedness, turnover and total betadiversity => global tests
set.seed(123)
global_nested <- adonis2(dist_nested ~ Village + Treatment, 
                         data = metadata, 
                         by = "margin",
                         permutations = 999)
print(global_nested)

set.seed(123)
global_turnover <- adonis2(dist_turnover ~ Village + Treatment, 
                           data = metadata, 
                           by = "margin",
                           permutations = 999)
print(global_turnover)

# 9b) Result outputs for nestedness, turnover and total betadiversity => pairwise comparisons
pw_nested <- pairwise_adonis_custom(dist_nested, metadata)
print(pw_nested)

pw_turnover <- pairwise_adonis_custom(dist_turnover, metadata)
print(pw_turnover)

pw_total <- pairwise_adonis_custom(dist_total, metadata)
print(pw_total)


