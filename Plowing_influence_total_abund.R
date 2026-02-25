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
library(DHARMa)
library(betapart)
library(multcompView)

wheat <- read_excel("wheat.xlsx") %>% filter(Treatment != "EKOLOGIE") %>%
  mutate(
    Date_str  = as.character(Date),
    Month_Num = as.numeric(stringr::str_split_i(Date_str, "\\.", 2)),
    Month = factor(
      Month_Num,
      levels = 6:10,
      labels = c("June", "July", "August", "September", "October")
    )
  ) %>%
  mutate(across(c(Treatment, Crop, Village, Locality, Trap, Species, Plowing), as.factor)) %>%
  droplevels() 
abundance_data <- wheat %>%
  group_by(Village, Locality, Trap, Treatment, Crop, Month, Plowing) %>%
  summarise(
    total_abundance = sum(Count, na.rm = TRUE),
    richness = n_distinct(Species[Count > 0]),
    .groups = "drop"
  ) %>%
  droplevels()

glm_rich_pois <- glmmTMB(
  total_abundance ~ Village + Treatment + Plowing +
    (1 | Locality/Trap),
  family = poisson(link = "log"),
  data = abundance_data
)
performance::check_overdispersion(glm_rich_pois)
# Overdispersion detected.

glm_abund <- glmmTMB(total_abundance ~ Village + Treatment + Plowing + (1|Locality/Trap),
                     family = nbinom2(link="log"),
                     data = abundance_data)
summary(glm_abund)
Anova(glm_abund, type = "II")
sim_res <- simulateResiduals(glm_abund)
plot(sim_res)

plotResiduals(sim_res, form = abundance_data$Treatment)
plotResiduals(sim_res, form = abundance_data$Crop)
plotResiduals(sim_res, form = abundance_data$Village)

emm_treatment <- emmeans(glm_abund, ~ Treatment, type = "response")
emm_df <- as.data.frame(emm_treatment)
head(emm_df)
pairwise_results <- pairs(emm_treatment)
print(pairwise_results)


library(ggpubr)

abundance_data$Plowing <- factor(abundance_data$Plowing, levels = c("Před", "Po"))
abundance_data$Treatment <- factor(abundance_data$Treatment, levels = c("KONVENCE", "REGENERACE"))

d1<-ggplot(abundance_data, aes(x = Treatment, y = total_abundance, fill = Plowing)) +
  # Boxploty s černým okrajem a šedou/bílou výplní
  geom_boxplot(outlier.shape = NA, color = "black", alpha = 1, width = 0.7, position = position_dodge(0.8)) +
  
  # Body (jitter) v černé barvě s průhledností
  geom_jitter(aes(group = Plowing), 
              position = position_dodge(0.8), 
              size = 1.2, alpha = 0.3, color = "black") +
  
  # Nastavení barev na černobílou škálu (bílá vs. světle šedá)
  scale_fill_manual(values = c("Před" = "white", "Po" = "grey85")) +
  
  # Párová porovnání (statistika)
  # Metoda Wilcoxon je bezpečná pro ne-normální data, 
  # label = "p.signif" zobrazí hvězdičky, "p.format" zobrazí přesné číslo
  stat_compare_means(aes(group = Plowing), 
                     label = "p.format", 
                     method = "wilcox.test", 
                     label.y = max(abundance_data$total_abundance) * 1.05) +
  
  # Popisky a téma
  labs(
    x = "Typ zemědělského hospodaření",
    y = "Odhad celkové početnosti",
  ) +
  theme_classic() + # Čisté bílé pozadí bez mřížek
  theme(
    legend.position = "top",
    axis.line = element_line(color = "black"),
    axis.text = element_text(color = "black", size = 10),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold")
  )
d1
tiff("Početnost_orba.tiff", width = 8, height = 6, units = "in", res = 300, compression = "lzw")
print(d1)
dev.off()

richness_data <- wheat %>%
  group_by(Village, Locality, Trap, Treatment, Crop,Month, Plowing) %>%
  summarise(
    total_richness = n_distinct(Species[Count > 0]),
    .groups = "drop"
  )  %>%
  droplevels()
head(richness_data)

glm_rich_pois <- glmmTMB(
  total_richness ~ Village + Treatment + Plowing +
    (1 | Locality/Trap),
  family = poisson(link = "log"),
  data = richness_data
)
performance::check_overdispersion(glm_rich_pois)
# Species richness per trap was analyzed using a Poisson GLMM (log link). 
# Overdispersion diagnostics indicated no evidence of overdispersion (dispersion ratio = 0.85).

glm_rich <- glmmTMB(total_richness ~ Village + Treatment + Plowing + (1 | Locality/Trap) + (1|Month),
                    family = poisson(link="log"),
                    data = richness_data)

summary(glm_rich)
Anova(glm_rich, type = "II")
sim_res <- simulateResiduals(glm_rich)
plot(sim_res)

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