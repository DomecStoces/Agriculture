library(readxl)
library(dplyr)
library(tidyr)
library(vegan)
library(glmmTMB)
library(emmeans)
library(car)
library(ggplot2)

long_format <- read_excel("long_format.xlsx")
long_format <- long_format %>%
  mutate(
    Treatment = as.factor(Treatment),
    Crop = as.factor(Crop),
    Village = as.factor(Village),
    Locality = as.factor(Locality),
    Trap = as.factor(Trap),
    Species = as.factor(Species),
    Date = as.factor(Date)
  )
long_format <- long_format %>%
  filter(Crop != "Pšenice tvrdá")

filtered_data <- long_format %>%
  filter(Vole_dist == 1 | Digestate_fertilization == 1)

### 1. Total abundance and species richness ###

# Test for overdispersion
glm_pois <- glm(total_abundance ~ Treatment * Crop,
                family = poisson,
                data = abundance_data)
dispersion <- sum(residuals(glm_pois, type="pearson")^2) /
  df.residual(glm_pois)
dispersion

### 1a) Total abundance per trap ###
# Rank defficient for: TreatmentREGENERACE:CropPšenice, TreatmentKONVENCE:CropPšenice tvrdá, TreatmentREGENERACE:CropPšenice tvrdá

abundance_data <- long_format %>%
  group_by(Village, Locality, Trap, Treatment, Crop,Date) %>%
  summarise(
    total_abundance = sum(Count),
    richness = n_distinct(Species[Count > 0]),
    .groups = "drop"
  )

glm_abund <- glmmTMB(total_abundance ~ Village + Crop + Treatment + (1 | Locality) + (1|Date),
                      family = nbinom2(link="log"),
                      data = abundance_data)
summary(glm_abund)
Anova(glm_abund, type = "II")
emm_treatment <- emmeans(glm_abund, ~ Treatment, type = "response")
emm_df <- as.data.frame(emm_treatment)
head(emm_df)
pairwise_results <- pairs(emm_treatment)
print(pairwise_results)
d1<-ggplot(emm_df, aes(x = Treatment, y = response, color = Treatment)) +
  # Updated column names: ymin = asymp.LCL, ymax = asymp.UCL
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), 
                width = 0.2, linewidth = 1) +
  # Add the points for the estimated means
  geom_point(size = 4) +
  # Customize the labels
  labs(
       x = "Treatment",
       y = "Estimated Total Abundance") + scale_y_log10()+
  # Apply a clean theme
  theme_classic() +
  theme(
    legend.position = "none", 
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, hjust = 0.5)
  )
d1

### 1b) Species richness per trap ###
richness_data <- long_format %>%
  group_by(Village, Locality, Trap, Treatment, Crop,Date) %>%
  summarise(
    total_richness = n_distinct(Species[Count > 0]),
    .groups = "drop"
  )
head(richness_data)
glm_rich <- glmmTMB(total_richness ~ Village + Crop + Treatment + (1 | Locality) + (1|Date),
                      family = nbinom2(link="log"),
                      data = richness_data)
summary(glm_rich)
Anova(glm_rich, type = "II")
emm_treatment <- emmeans(glm_rich, ~ Treatment, type = "response")
emm_df <- as.data.frame(emm_treatment)
head(emm_df)
pairwise_results <- pairs(emm_treatment)
print(pairwise_results)
d2<-ggplot(emm_df, aes(x = Treatment, y = response, color = Treatment)) +
  # Updated column names: ymin = asymp.LCL, ymax = asymp.UCL
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), 
                width = 0.2, linewidth = 1) +
  # Add the points for the estimated means
  geom_point(size = 4) +
  # Customize the labels
  labs(
    x = "Treatment",
    y = "Estimated Species Richness") + scale_y_log10()+
  # Apply a clean theme
  theme_classic() +
  theme(
    legend.position = "none", 
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, hjust = 0.5)
  )
d2
