# =========================================
# 1. Convert community matrix to presence/absence
# =========================================
comm_pa <- ifelse(comm_matrix > 0, 1, 0)

# =========================================
# 2. Partition beta diversity (Sørensen)
# =========================================
beta_part_sor <- beta.pair(comm_pa, index.family = "sorensen")

# Extract components
dist_turnover <- beta_part_sor$beta.sim   # turnover
dist_nested <- beta_part_sor$beta.sne     # nestedness-resultant
dist_total <- beta_part_sor$beta.sor      # total Sørensen

# Square-root transform to make distances Euclidean-compatible for betadisper
dist_turnover_sqrt <- sqrt(dist_turnover)
dist_nested_sqrt   <- sqrt(dist_nested)
dist_total_sqrt    <- sqrt(dist_total)

# =========================================
# 3. Test homogeneity of dispersions (PERMDISP)
# =========================================
disp_turnover <- betadisper(dist_turnover_sqrt, metadata$Treatment)
disp_nested   <- betadisper(dist_nested_sqrt, metadata$Treatment)
disp_total    <- betadisper(dist_total_sqrt, metadata$Treatment)

# Global tests
perm_turnover <- permutest(disp_turnover, permutations = 999)
perm_nested   <- permutest(disp_nested, permutations = 999)
perm_total    <- permutest(disp_total, permutations = 999)

# =========================================
# 4. Extract distances and prepare plot data
# =========================================
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

# =========================================
# 5. Generate significance letters automatically (Tukey HSD)
# =========================================
get_sig_letters <- function(disp_obj) {
  tuk <- TukeyHSD(disp_obj)
  
  # Convert TukeyHSD result to named vector
  pvals <- tuk$Treatment[, "p adj"]
  names(pvals) <- rownames(tuk$Treatment)   # ensure names are character
  
  library(multcompView)
  cld <- multcompLetters(pvals)$Letters
  
  # Extract individual treatments from letters
  # Get all unique treatment names
  trt_names <- unique(unlist(strsplit(names(cld), "-")))
  
  # Map letters to treatments
  df <- data.frame(Treatment = trt_names,
                   label = NA,
                   stringsAsFactors = FALSE)
  
  for(i in seq_along(df$Treatment)) {
    # find all Tukey pairs containing this treatment
    idx <- grep(df$Treatment[i], names(cld))
    # take the letter assigned in multcompLetters (just the first one)
    df$label[i] <- cld[idx[1]]
  }
  
  return(df)
}

# =========================================
# 6. Plot boxplots with significance letters
# =========================================
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
# =========================================
# 7. PERMANOVA on total beta diversity
# =========================================
set.seed(123)
permanova_total <- adonis2(dist_total ~ Village + Crop + Treatment, 
                           data = metadata, 
                           by = "margin", 
                           strata = metadata$Locality,
                           permutations = 999)
print("--- PERMANOVA: TOTAL BETA DIVERSITY ---")
print(permanova_total)

# =========================================
# 8. Pairwise PERMANOVA for Treatment
# =========================================
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
  
  res <- adonis2(dist_sub ~ Village + Crop + Treatment,
                 data = meta_sub,
                 by = "margin",
                 strata = meta_sub$Locality,
                 permutations = 999)
  
  pairwise_permanova_results <- rbind(pairwise_permanova_results, 
                                      data.frame(Pair = paste(pair[1], "vs", pair[2]),
                                                 P_Value = res["Treatment", "Pr(>F)"]))
}

# BH correction
pairwise_permanova_results$P_Adj_BH <- p.adjust(pairwise_permanova_results$P_Value, method = "BH")
print("--- PAIRWISE PERMANOVA: Composition Differences ---")
print(pairwise_permanova_results)

### Pairwise comparisons of TURNOVER and NESTEDNESS DISPERSION ###
TukeyHSD(disp_turnover)
TukeyHSD(disp_nested)
TukeyHSD(disp_total)
