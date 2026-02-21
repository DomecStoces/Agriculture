library(readxl)
library(dplyr)
library(tidyr)
library(vegan)
library(labdsv)
library(ggrepel)
library(ggnewscale)
library(indicspecies)
library(glmmTMB)
library(emmeans)
library(car)
library(ggplot2)
library(stringr)

# Data preparation
long_format <- read_excel("long_format.xlsx") %>%
  filter(Crop != "Pšenice tvrdá") %>%
  mutate(
    Date_str  = as.character(Date),
    Month_Num = as.numeric(stringr::str_split_i(Date_str, "\\.", 2)),
    Month = factor(
      Month_Num,
      levels = 6:10,
      labels = c("June", "July", "August", "September", "October")
    )
  ) %>%
  mutate(across(c(Treatment, Crop, Village, Locality, Trap, Species), as.factor)) %>%
  droplevels() 

### 1a) Total abundance per trap ###
# Rank defficient for: TreatmentREGENERACE:CropPšenice, TreatmentKONVENCE:CropPšenice tvrdá, TreatmentREGENERACE:CropPšenice tvrdá
abundance_data <- long_format %>%
  group_by(Village, Locality, Trap, Treatment, Crop, Month) %>%
  summarise(
    total_abundance = sum(Count, na.rm = TRUE),
    richness = n_distinct(Species[Count > 0]),
    .groups = "drop"
  ) %>%
  droplevels()

glm_rich_pois <- glmmTMB(
  total_abundance ~ Village + Crop + Treatment +
    (1 | Locality) + (1 | Month),
  family = poisson(link = "log"),
  data = abundance_data
)
performance::check_overdispersion(glm_rich_pois)
# Overdispersion detected.

glm_abund <- glmmTMB(total_abundance ~ Village + Crop + Treatment + (1 | Locality) + (1|Month),
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
  labs(x = "Treatment",
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
  group_by(Village, Locality, Trap, Treatment, Crop,Month) %>%
  summarise(
    total_richness = n_distinct(Species[Count > 0]),
    .groups = "drop"
  )  %>%
  droplevels()
head(richness_data)

glm_rich_pois <- glmmTMB(
  total_richness ~ Village + Crop + Treatment +
    (1 | Locality) + (1 | Month),
  family = poisson(link = "log"),
  data = richness_data
)
performance::check_overdispersion(glm_rich_pois)
# Species richness per trap was analyzed using a Poisson GLMM (log link). 
# Overdispersion diagnostics indicated no evidence of overdispersion (dispersion ratio = 0.85).

glm_rich <- glmmTMB(total_richness ~ Village + Crop + Treatment + (1 | Locality) + (1|Month),
                      family = poisson(link="log"),
                      data = richness_data)
summary(glm_rich)
Anova(glm_rich, type = "II")
emm_treatment <- emmeans(glm_rich, ~ Treatment, type = "response")
emm_df <- as.data.frame(emm_treatment)
head(emm_df)
pairwise_results <- pairs(emm_treatment)
print(pairwise_results)
d2 <- ggplot(emm_df, aes(x = Treatment, y = rate, color = Treatment)) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL),
                width = 0.2, linewidth = 1) +
  geom_point(size = 4) +
  labs(
    x = "Treatment",
    y = "Estimated Species Richness"
  ) +
  scale_y_log10() +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, hjust = 0.5)
  )

d2
