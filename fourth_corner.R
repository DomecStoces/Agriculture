library(dplyr)
library(tidyr)
library(ade4)
library(ggplot2)

long_format <- long_format %>%
  mutate(SampleID = paste(Village, Locality, Trap, Date, sep = "_"))

# --- 3) L MATRIX: Sites x Species Abundances ----------------------------------
L_data <- long_format %>%
  group_by(SampleID, Species) %>%
  summarise(Count = sum(Count, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Species, values_from = Count, values_fill = 0)

L <- as.data.frame(L_data %>% select(-SampleID))
rownames(L) <- L_data$SampleID

# --- 4) R MATRIX: Sites x Environmental Variables -----------------------------
R_data <- long_format %>%
  group_by(SampleID) %>%
  slice(1) %>% # Keep one row per site
  ungroup() %>%
  # Notice we include 'Month' here instead of 'Date'
  select(SampleID, Village, Locality, Vole, Trap, Month, Treatment, Crop, 
         Vole_dist, Digestate_fertilization)

R <- as.data.frame(R_data %>% select(-SampleID))
rownames(R) <- R_data$SampleID

# Coerce factors and numerics appropriately for dudi.hillsmith
R$Village   <- as.factor(R$Village)
R$Locality  <- as.factor(R$Locality)
R$Vole      <- as.factor(R$Vole)
R$Trap      <- as.factor(R$Trap)
R$Month     <- as.factor(R$Month)       
R$Treatment <- as.factor(R$Treatment)
R$Crop      <- as.factor(R$Crop)
R$Vole_dist <- as.numeric(R$Vole_dist)
R$Digestate_fertilization <- as.numeric(R$Digestate_fertilization)

# --- 5) Q MATRIX: Species x Traits --------------------------------------------
Q_data <- long_format %>%
  group_by(Species) %>%
  slice(1) %>% 
  ungroup() %>%
  select(Species, Diet, Bioindication, Moisture, Size, Breeding, 
         Wings, Biogeographical.affinity, IUCN, ZVL)

Q_raw <- as.data.frame(Q_data %>% select(-Species))
rownames(Q_raw) <- Q_data$Species

# --- 6) Fuzzy Traits Handling & Dummy Conversion ------------------------------
# Create missing target columns if needed
ensure_cols <- function(df, cols) {
  missing <- setdiff(cols, colnames(df))
  for (m in missing) df[[m]] <- 0
  df
}

# Row-normalize a set of columns so rows sum to 1
row_norm <- function(df, cols) {
  ok <- intersect(cols, colnames(df))
  if (length(ok) > 1) {
    rs <- rowSums(df[, ok, drop = FALSE])
    rs[rs == 0] <- 1
    df[, ok] <- sweep(df[, ok, drop = FALSE], 1, rs, "/")
  }
  df
}

# Convert standard categorical traits into dummy variables (0/1)
cat_cols <- c("Diet", "Bioindication", "Moisture", "Breeding", "IUCN", "ZVL")
Q_num <- Q_raw

for (col in cat_cols) {
  lvls <- unique(Q_num[[col]])
  for (l in lvls) {
    col_name <- paste(col, l, sep = "_")
    Q_num[[col_name]] <- ifelse(Q_num[[col]] == l, 1, 0)
  }
  Q_num[[col]] <- NULL # remove the original categorical column
}

# Fuzzy split for "Wings" (B, M, and the composite B/M)
Q_num$Wings_M  <- ifelse(Q_num$Wings == "M", 1, 0)
Q_num$Wings_B  <- ifelse(Q_num$Wings == "B", 1, 0)
Q_num$Wings_BM <- ifelse(Q_num$Wings == "B/M", 1, 0)

# Split 'Wings_BM' composite equally into 'B' and 'M'
Q_num$Wings_B <- Q_num$Wings_B + 0.5 * Q_num$Wings_BM
Q_num$Wings_M <- Q_num$Wings_M + 0.5 * Q_num$Wings_BM
Q_num$Wings_BM <- NULL
Q_num$Wings    <- NULL

# Ensure everything left is numeric
Q_num[] <- lapply(Q_num, function(x) as.numeric(as.character(x)))
Q_num[is.na(Q_num)] <- 0

# Row-normalize categorical dummy groups so they represent proportions
Q_num <- row_norm(Q_num, grep("Diet_", names(Q_num), value = TRUE))
Q_num <- row_norm(Q_num, grep("Bioindication_", names(Q_num), value = TRUE))
Q_num <- row_norm(Q_num, grep("Moisture_", names(Q_num), value = TRUE))
Q_num <- row_norm(Q_num, grep("Breeding_", names(Q_num), value = TRUE))
Q_num <- row_norm(Q_num, c("Wings_B", "Wings_M"))
Q_num <- row_norm(Q_num, grep("IUCN_", names(Q_num), value = TRUE))
Q_num <- row_norm(Q_num, grep("ZVL_", names(Q_num), value = TRUE))

# --- 7) Align species & sites between matrices --------------------------------
# Just to be strictly safe
common_sp <- intersect(colnames(L), rownames(Q_num))
L <- L[, common_sp, drop = FALSE]
Q <- Q_num[common_sp, , drop = FALSE]

common_sites <- intersect(rownames(L), rownames(R))
L <- L[common_sites, , drop = FALSE]
R <- R[common_sites, , drop = FALSE]

# --- 8) RLQ Analysis ----------------------------------------------------------
# Species abundance
afcL <- dudi.coa(L, scannf = FALSE)
# Environment (All variables)
acpR <- dudi.hillsmith(R, row.w = afcL$lw, scannf = FALSE)   
# Traits 
acpQ <- dudi.pca(Q, row.w = afcL$cw, scale = TRUE, scannf = FALSE)

# Run global RLQ
rlq_res <- rlq(acpR, afcL, acpQ, scannf = FALSE)
plot(rlq_res)

# --- 9) Fourth-corner Analysis ------------------------------------------------
R_treat <- R[, "Treatment", drop = FALSE]
set.seed(123) 
fc_treat <- fourthcorner(
  tabR = R_treat, 
  tabL = L, 
  tabQ = Q, 
  modeltype = 6, 
  p.adjust.method.G = "none",
  p.adjust.method.D = "none",
  nrepet = 999
)

# --- 10) Plot and Save Outputs ------------------------------------------------
par(mfrow = c(1, 2), mar = c(4,4,2,1))

# Plot FC1
plot(fc_treat, alpha = 0.05, stat = "D2")
# Plot FC2 Biplot
plot(fc_treat, x.rlq = rlq_res, alpha = 0.05, stat = "D2", type = "biplot")

fc2 <- fourthcorner2(R_treat, L, Q, modeltype = 6, p.adjust.method.G = "none", nrepet = 999)
fc2$trRLQ
# Save figure to TIFF
tiff("fourth_corner.tiff", units = "in", width = 10, height = 6, res = 600)
plot(fc, alpha = 0.05, stat = "D2")
dev.off()