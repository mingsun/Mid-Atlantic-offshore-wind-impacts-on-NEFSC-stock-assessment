library(tidyverse)

# 1. catch by tow data ----


full_tow_df <- read.csv("results/stock assessment/squid/tow data/full.tow.list.for.assessment.csv") %>%
  # mutate(STRATUM = as.character(STRATUM)) %>%
  group_by(YEAR, SEASON, STRATUM) %>%
  mutate(TOTAL_N_STATION = length(STATION))

write.csv(full_tow_df, "results/indices for assessment/squid/catch_by_tow.csv", row.names = FALSE)


# ----------------------------------------------------------------------------------- #


# 2. original abundance indices ----
  
  ## 2.1 mean numbers by stratum ---------------------------------------------------------------------------------------------
  # for squid, the index-based method uses biomass

mean_N_stratum_df <- full_tow_df %>%
  group_by(YEAR, SEASON, STRATUM) %>%
  mutate(MEAN_N_STRATUM = sum(CATCH_WT_CAL, na.rm = TRUE)/TOTAL_N_STATION) %>% # mean within a strata
  mutate(VAR_STRATUM = sum((CATCH_WT_CAL - MEAN_N_STRATUM)^2, na.rm = TRUE)/ (TOTAL_N_STATION - 1)) %>% # variance by stratum
  # mutate(VAR_STRATUM = var(EXPCATCHNUM)) %>% # variance by stratum, same as last line
  ungroup() 


## 2.2 stratified mean numbers ---------------------------------------------------------------------------------------------

stratified_mean_N_df <- mean_N_stratum_df %>%
  select(c(YEAR, SEASON, STRATUM, TOTAL_N_STATION, STRATUM_AREA, MEAN_N_STRATUM, VAR_STRATUM)) %>%
  # select(c(YEAR, SEASON, STRATUM, N.STATION, REL_WEIGHT, MEAN_N_STRATUM)) %>%
  distinct() %>% # downsize the data frame to a minimal without repetitive rows
  group_by(YEAR, SEASON) %>%
  mutate(REL_WEIGHT = STRATUM_AREA/sum(STRATUM_AREA, na.rm = TRUE)) %>%
  summarize(STRATIFIED_MEAN_N = weighted.mean(MEAN_N_STRATUM, w = REL_WEIGHT), # stratified mean
            STRATIFIED_VAR = sum(REL_WEIGHT^2 * VAR_STRATUM / TOTAL_N_STATION, na.rm = TRUE), # variance
            STRATIFIED_SE = sqrt(STRATIFIED_VAR)) %>%   # standard deviance
  mutate(up_CI_95 = qnorm(0.975, mean = STRATIFIED_MEAN_N, sd = STRATIFIED_SE),
         lo_CI_95 = qnorm(0.025, mean = STRATIFIED_MEAN_N, sd = STRATIFIED_SE),
         lo_CI_95 = ifelse(lo_CI_95 < 0, yes = 0, no = lo_CI_95),
         CV = STRATIFIED_SE / STRATIFIED_MEAN_N) %>%
  ungroup()


write.csv(stratified_mean_N_df, "results/indices for assessment/squid/original.indices.csv", row.names = FALSE)

# ----------------------------------------------------------------------------------- #



# 3. WEE indices ---------------------------------------------------------------------------------------------

## filter catch by excluding the impacted tow
impacted_cat_df <- full_tow_df %>%
  filter(OWF  == "OUTSIDE") %>%
  group_by(YEAR, SEASON, STRATUM) %>%
  mutate(TOTAL_N_STATION = length(STATION))


  ## 3.1 mean numbers by stratum ----
impacted_mean_N_stratum_df <- impacted_cat_df %>%
  group_by(YEAR, SEASON, STRATUM) %>%
  # mutate(MEAN_N_STRATUM = sum(BIOMASS, na.rm = TRUE)/TOTAL_N_STATION) %>% # mean within a strata
  # mutate(VAR_STRATUM = sum((BIOMASS - MEAN_N_STRATUM)^2, na.rm = TRUE)/ (TOTAL_N_STATION - 1)) %>% # variance by stratum
  mutate(MEAN_N_STRATUM = sum(CATCH_WT_CAL, na.rm = TRUE)/TOTAL_N_STATION) %>% # mean within a strata
  mutate(VAR_STRATUM = sum((CATCH_WT_CAL - MEAN_N_STRATUM)^2, na.rm = TRUE)/ (TOTAL_N_STATION - 1)) %>% # variance by stratum
  ungroup() 


  ## 3.2 stratified mean numbers ----
impacted_stratified_mean_N_df <- impacted_mean_N_stratum_df %>%
  select(c(YEAR, SEASON, STRATUM, TOTAL_N_STATION, STRATUM_AREA, MEAN_N_STRATUM, VAR_STRATUM)) %>%
  # select(c(YEAR, SEASON, STRATUM, N.STATION, REL_WEIGHT, MEAN_N_STRATUM)) %>%
  group_by(YEAR, SEASON) %>%
  mutate(REL_WEIGHT = STRATUM_AREA/sum(STRATUM_AREA, na.rm = TRUE)) %>%
  summarize(STRATIFIED_MEAN_N = weighted.mean(MEAN_N_STRATUM, w = REL_WEIGHT), # stratified mean
            STRATIFIED_VAR = sum(REL_WEIGHT^2 * VAR_STRATUM / TOTAL_N_STATION, na.rm = TRUE), # variance
            STRATIFIED_SE = sqrt(STRATIFIED_VAR)) %>%   # standard deviance
  mutate(up_CI_95 = qnorm(0.975, mean = STRATIFIED_MEAN_N, sd = STRATIFIED_SE),
         lo_CI_95 = qnorm(0.025, mean = STRATIFIED_MEAN_N, sd = STRATIFIED_SE),
         lo_CI_95 = ifelse(lo_CI_95 < 0, yes = 0, no = lo_CI_95),
         CV = STRATIFIED_SE / STRATIFIED_MEAN_N) %>%
  ungroup()


write.csv(impacted_stratified_mean_N_df, "results/indices for assessment/squid/WEE.indices.csv", row.names = FALSE)


# ----------------------------------------------------------------------------------- #




