library(tidyverse)

# 1. review abundance indices trend  ---------------------------------------------------------------

original_AI <- read.csv("results/stratified.mean.indices/longfin.squid/original.indices.csv") %>%
  add_column(MODEL = "original") %>%
  rename(fit = STRATIFIED_MEAN_N, lwr = lo_CI_95, upr = up_CI_95) %>%
  select(YEAR, SEASON, MODEL, fit, lwr, upr)


impacted_AI <- read.csv("results/stratified.mean.indices/longfin.squid/impacted.indices.csv") %>%
  add_column(MODEL = "impacted") %>%
  rename(fit = STRATIFIED_MEAN_N, lwr = lo_CI_95, upr = up_CI_95) %>%
  select(YEAR, SEASON, MODEL, fit, lwr, upr)

AI_df <- rbind(original_AI, impacted_AI) 

design_based_p <- ggplot(AI_df, aes(x = YEAR)) +
  # geom_ribbon(data = random_abundance_df, aes(ymin = QUANT_0.025, ymax = QUANT_0.975), fill = "grey", alpha = 0.3) + # first plot the ribbon for the random estimate
  geom_errorbar(aes(ymin = lwr, ymax = upr, color = MODEL), width = 0.2, position = position_dodge(width = 0.5)) +
  geom_point(aes(y = fit, color = MODEL), position = position_dodge(width = 0.5)) +
  geom_line(aes(y = fit, color = MODEL)) +
  facet_wrap(.~SEASON, nrow = 2) +
  scale_color_manual(values = c("indianred2", "steelblue3"), labels = c("impacted", "original"), name = "Time Series") +
  ylab("biomass indices [g/tow]") +
  theme_minimal() +
  theme(legend.position = c(0.05, 0.95), # Place the legend in the top-left corner
        legend.justification = c(0, 1))

#  ------------------------------------------------------------------------------------------------ #


# 2. total biomass calculation  ---------------------------------------------------------------

# all index-based method results are copied from the squid biomass calculation sheets from
# E:\Stony Brook job\NYSERDA Offshore wind\assessment model\squid
# a summary is pasted to E:\Stony Brook job\NYSERDA Offshore wind\results\stock assessment\squid

# the outputs are pasted as summary sheets available in
# in E:\Stony Brook job\NYSERDA Offshore wind\results\stock assessment\squid
# a summary is pasted to E:\Stony Brook job\NYSERDA Offshore wind\results\stock assessment\squid
