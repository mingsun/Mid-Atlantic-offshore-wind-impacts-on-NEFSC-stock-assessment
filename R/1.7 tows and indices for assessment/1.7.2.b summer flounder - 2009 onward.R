library(tidyverse)


# the difference between this script and 1.7.2 is:
# for 2009 onward, when it is an individual abundance indices series, no calibration is needed for the vessel change
# but SWAN adjustment is needed as below (66th SAW)

  # SAM: The assessment model splits the Bigelow and Albatross survey indices. 
       # In preparation for application to the assessment model the Bigelow index most notably is scaled up for by-tow swept area, 
       # so that's the major difference that I assume you're noticing. It would be difficult to incorporate those modifications since they only cover 
       # the Bigelow years and I don't think those types of modifications would change your overall conclusions anyhow. 
       # There is some information in SAW 66 (66th Northeast Regional Stock Assessment Workshop (66th SAW) Assessment Report). 
       # If you search for "SWAN" (swept area numbers) you will find the relevant text.

# technical parameters: 


# 1. catch by tow data ----

  ## tow list information
full_tow_df <- read.csv("results/stock assessment/summer flounder/tow data/full.tow.list.for.assessment.csv") %>%
  filter(YEAR >= 2009) %>%  # Bigelow tow after 2009
  mutate(STRATUM = as.character(STRATUM)) %>%
  group_by(YEAR, SEASON, STRATUM) %>%
  mutate(TOTAL_N_STATION = length(STATION))


  ## load and apply the calibration factor at length
load("results/BTS_catch_at_len_by_tow.Rdata")

cat_df <- catch_at_len_df %>%
  filter(SVSPP == 103) %>%
  select(-c(CRUISE, SVSPP, CATCHSEX, CATCHSEX)) %>%
  mutate(YEAR = as.numeric(substr(CRUISE6, 1, 4))) %>%
  filter(YEAR >= 2009) %>%  # Bigelow tow after 2009
  mutate(STRATUM = ifelse(nchar(STRATUM) == 5, str_sub(STRATUM, 2, 5), STRATUM)) %>%
  mutate(ID = paste(CRUISE6, STRATUM, TOW, sep = "."),
         BLOCK = paste(YEAR, SEASON, STRATUM, sep = ".")) %>%
  group_by(YEAR, SEASON, BLOCK, CRUISE6, STRATUM, TOW, STATION, ID) %>%
  summarize(NUMBER = sum(EXPNUMLEN))


  ## add zero catch tows
cat_df <- cat_df %>%
  right_join(full_tow_df) %>% # right join to add zero catch tows
  mutate(NUMBER = ifelse(is.na(NUMBER), 0, NUMBER)) %>%
  select(YEAR, SEASON, STRATUM, BLOCK, everything()) %>%
  arrange(BLOCK)



# ----------------------------------------------------------------------------------- #





# 2. original abundance indices ----


  ## 2.1 mean numbers by stratum ---------------------------------------------------------------------------------------------
    # all methods below follow (https://noaa-edab.github.io/survdat/articles/calc_strat_mean.html)

mean_N_stratum_df <- cat_df %>%
  group_by(YEAR, SEASON, STRATUM) %>%
  mutate(MEAN_N_STRATUM = mean(NUMBER)) %>% # mean within a strata
  mutate(VAR_STRATUM = sum((NUMBER - MEAN_N_STRATUM)^2, na.rm = TRUE)/ (TOTAL_N_STATION - 1)) %>% # variance by stratum
  # mutate(VAR_STRATUM = var(EXPCATCHNUM)) %>% # variance by stratum, same as last line
  ungroup() 



  ## 2.2 stratified mean numbers ---------------------------------------------------------------------------------------------

stratified_mean_N_df <- mean_N_stratum_df %>%
  select(c(YEAR, SEASON, STRATUM, TOTAL_N_STATION, STRATUM_AREA, MEAN_N_STRATUM, VAR_STRATUM)) %>%
  distinct() %>% # downsize the data frame to a minimal without repetitive rows
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

write.csv(stratified_mean_N_df, "results/indices for assessment/summer flounder/original.indices.BIG.csv", row.names = FALSE)




  ## 2.3 catch efficiency  ---------------------------------------------------------------------------------------------

# catch efficiency are length-specific, so we calculate the number by length here
# compare with Table A50 of SW66





catch_efficiency_df <- read.csv("data/NOAA.stock.data/summer.flounder/STOCKEFF_SV_172735_UNIT_NONE_length_based_catch_efficiency.csv")

load("data/NOAA.stock.data/summer.flounder/fluke_N.W.RData")




# ----------------------------------------------------------------------------------- #





# 3. WEE indices ----

  ## filter catch by excluding the impacted tow
impacted_cat_df <- cat_df %>%
  filter(OWF  == "OUTSIDE") %>%
  group_by(YEAR, SEASON, STRATUM) %>%
  mutate(TOTAL_N_STATION = length(STATION))


  ## 3.1 mean numbers by stratum ----
impacted_mean_N_stratum_df <- impacted_cat_df %>%
  group_by(YEAR, SEASON, STRATUM) %>%
  mutate(MEAN_N_STRATUM = mean(NUMBER)) %>% # mean within a strata
  mutate(VAR_STRATUM = sum((NUMBER - MEAN_N_STRATUM)^2, na.rm = TRUE)/ (TOTAL_N_STATION - 1)) %>% # variance by stratum
  mutate(VAR = var(NUMBER)) %>% # variance by stratum
  ungroup() 


  ## 3.2 stratified mean numbers ----
impacted_stratified_mean_N_df <- impacted_mean_N_stratum_df %>%
  select(c(YEAR, SEASON, STRATUM, TOTAL_N_STATION, STRATUM_AREA, MEAN_N_STRATUM, VAR_STRATUM)) %>%
  distinct() %>% # downsize the data frame to a minimal without repetitive rows
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

write.csv(impacted_stratified_mean_N_df, "results/indices for assessment/summer flounder/WEE.indices.BIG.csv", row.names = FALSE)


# ----------------------------------------------------------------------------------- #


