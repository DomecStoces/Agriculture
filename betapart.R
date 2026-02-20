library(betapart)

metadata$Treatment <- factor(metadata$Treatment)
metadata$Village   <- factor(metadata$Village)
metadata$Crop      <- factor(metadata$Crop)
# 1. Convert to Presence/Absence
# betapart strictly requires binary (1/0) data for these specific indices
comm_pa <- ifelse(comm_matrix > 0, 1, 0)

# 2. Calculate Beta Diversity Partitions
# Jaccard family calculates the Jaccard-based turnover 
# and nestedness-resultant fractions (Baselga, 2012 framework)
beta_partitions <- beta.pair(comm_pa, index.family = "jaccard")

# Extract the three distance matrices for easy use:
# 1. Turnover (Species replacement)
dist_turnover <- beta_partitions$beta.jtu 
# 2. Nestedness (Species loss/gain)
dist_nestedness <- beta_partitions$beta.jne 
# 3. Total Beta Diversity (Jaccard)
dist_total <- beta_partitions$beta.jac

### Testing the partitions (PERMANOVA) ###
#TURNOVER
set.seed(123)
permanova_turnover <- adonis2(dist_turnover ~ Village + Crop + Treatment,
                               data = metadata,
                               by = "margin",
                               strata = metadata$Locality,
                               permutations = 999)
print(permanova_turnover)
# NESTEDNESS
set.seed(123)
permanova_nestedness <- adonis2(dist_nestedness ~ Village + Crop + Treatment, 
                                data = metadata, 
                                by = "margin", 
                                strata = metadata$Locality,
                                permutations = 999)
print(permanova_nestedness)
