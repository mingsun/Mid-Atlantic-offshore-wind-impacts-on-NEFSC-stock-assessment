library(tidyverse)
library(ASAPplots)


# 0. see the value of the indices value ----

data <- ReadASAP3DatFile("assessment model/summer flounder/ASAP input data/ASAP3_MTA2023_FINAL.DAT")
data$survey.names

# 2 is BTS ALB spring, 3 is BTS ALB fall, they are 1982-2008
# 25 is BTS BIG spring, 26 is BTS BIG fall, they are 2009-2022
# col.1 is year 1982-2022, 1982-2008 correspond to row 1-27, 2009-2022 correspond to row 28-41, 
# col.2 is sum value, col.3 is CV
# column 4-11 correspond to age 1-8,  col. 11 is sample size

ALB_spring <- data$dat$IAA_mats[[2]][1:27,c(2,4:11)]
ALB_fall <- data$dat$IAA_mats[[3]][1:27,c(2,4:11)]

BIG_spring <- data$dat$IAA_mats[[25]][28:41,c(2,4:11)] # 2020 not considered, will be set as unchanged in this script
BIG_fall <- data$dat$IAA_mats[[26]][28:41,c(2,4:11)] # 2017 and 2020 not considered, will be set as unchanged in this script


# ------------------------------------------------- #





# 2. ALB (1982 - 2008) ----



 ## 2.0 catch by tow data (original) ----

## tow list information
full_tow_df <- read.csv("results/stock assessment/summer flounder/tow data/full.tow.list.for.assessment.csv") %>%
  mutate(STRATUM = as.character(STRATUM)) %>%
  group_by(YEAR, SEASON, STRATUM) %>%
  mutate(TOTAL_N_STATION = length(STATION)) %>%
  filter(YEAR != 2023)


## Albatross IV decommission in 2008, the NEFSC transitioned to the NOAA Ship Henry B. Bigelow starting from spring 2009
## Bigelow.catch = Albatross.catch * calibration.factor
## summer flounder assessment uses all abundance-based indices, and separate Alb and Bigelow survey series, so no need to calibrate


## load and apply the calibration factor at length
load("results/BTS_catch_at_len_by_tow.Rdata")
caliFactor_df <- read.csv("data/NOAA.stock.data/summer.flounder/STOCKEFF_SV_172735_UNIT_NONE_length_based_calibration.csv")[,c(8:10)]

cat_df <- catch_at_len_df %>%
  filter(SVSPP == 103) %>%
  select(-c(CRUISE, SVSPP, CATCHSEX, CATCHSEX)) %>%
  mutate(YEAR = as.numeric(substr(CRUISE6, 1, 4))) %>%
  filter(YEAR >= 1982) %>% # summer flounder assessment uses survey data since 1982
  mutate(STRATUM = ifelse(nchar(STRATUM) == 5, str_sub(STRATUM, 2, 5), STRATUM)) %>%
  mutate(ID = paste(CRUISE6, STRATUM, TOW, sep = "."),
         BLOCK = paste(YEAR, SEASON, STRATUM, sep = ".")) %>%
  left_join(caliFactor_df) %>% # integrate with the calibration factors
  mutate(CALIBRATION_FACTOR = ifelse(YEAR > 2008, CALIBRATION_FACTOR, 1), # only apply calibration factor to >2008
         NUMBER_LEN = EXPNUMLEN / CALIBRATION_FACTOR) %>% # calculate the numbers
  group_by(YEAR, SEASON, BLOCK, CRUISE6, STRATUM, TOW, STATION, ID) %>%
  summarize(NUMBER = sum(NUMBER_LEN))


## add zero catch tows
cat_df <- cat_df %>%
  right_join(full_tow_df) %>% # right join to add zero catch tows
  mutate(NUMBER = ifelse(is.na(NUMBER), 0, NUMBER)) %>%
  select(YEAR, SEASON, STRATUM, BLOCK, everything()) %>%
  arrange(BLOCK)

remove(caliFactor_df, catch_at_len_df)



write.csv(cat_df, "results/indices for assessment/summer flounder/catch_by_tow_ALB.csv", row.names = FALSE)


  ## ---------------------------------- ##




 ## 2.1 full dataset ----

mean_N_stratum_df <- read.csv("results/indices for assessment/summer flounder/catch_by_tow_ALB.csv") %>%
  filter(YEAR <= 2008) %>%
  group_by(YEAR, SEASON, STRATUM) %>%
  mutate(MEAN_N_STRATUM = mean(NUMBER)) %>% # mean within a strata
  mutate(VAR_STRATUM = var(NUMBER), na.rm = TRUE) %>% # variance by stratum, same as last line
  ungroup() 

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

write.csv(stratified_mean_N_df, "results/indices for assessment/summer flounder/original.ALB.csv", row.names = FALSE)


 ## ------------------------ ##


rm(list = ls())




  ## 2.2 WEE dataset ----

WEE_cat_df <- read.csv("results/indices for assessment/summer flounder/catch_by_tow_ALB.csv")  %>%
  filter(OWF  == "OUTSIDE") %>%
  group_by(YEAR, SEASON, STRATUM) %>%
  mutate(TOTAL_N_STATION = length(STATION))


mean_N_stratum_df <- WEE_cat_df %>%
  filter(YEAR <= 2008) %>%
  group_by(YEAR, SEASON, STRATUM) %>%
  mutate(MEAN_N_STRATUM = mean(NUMBER)) %>% # mean within a strata
  mutate(VAR_STRATUM = var(NUMBER), na.rm = TRUE) %>% # variance by stratum, same as last line
  ungroup() 

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

write.csv(stratified_mean_N_df, "results/indices for assessment/summer flounder/WEE/WEE.ALB.csv", row.names = FALSE)


  ## extract WEE ratio relative to original
original.ALB_df <- read.csv("results/indices for assessment/summer flounder/original.ALB.csv")


ALB_WEE_SPRING_ratio <- subset(stratified_mean_N_df, SEASON == "SPRING")$STRATIFIED_MEAN_N/subset(original.ALB_df, SEASON == "SPRING")$STRATIFIED_MEAN_N
ALB_WEE_FALL_ratio <- subset(stratified_mean_N_df, SEASON == "FALL")$STRATIFIED_MEAN_N/subset(original.ALB_df, SEASON == "FALL")$STRATIFIED_MEAN_N


save(ALB_WEE_SPRING_ratio, ALB_WEE_FALL_ratio, file = "results/indices for assessment/summer flounder/WEE/WEE.ALB_ratio.Rdata")

rm(list =ls())



  ## ------------------------ ##




  ## 2.3 MIT.1 dataset ----


  # !!! the loaded catch_by_tow data are generated from R script 3.0

MIT.1_cat_df <- read.csv("results/indices for assessment/summer flounder/mitigation_1_WEA_distance/catch_by_tow_ALB.csv") %>%
  group_by(YEAR, SEASON, STRATUM) %>%
  mutate(TOTAL_N_STATION = length(STATION)) %>%
  filter(YEAR != 2023)

mean_N_stratum_df <- MIT.1_cat_df %>%
  filter(YEAR <= 2008) %>%
  group_by(YEAR, SEASON, STRATUM) %>%
  mutate(MEAN_N_STRATUM = mean(NUMBER)) %>% # mean within a strata
  mutate(VAR_STRATUM = var(NUMBER), na.rm = TRUE) %>% # variance by stratum, same as last line
  ungroup() 

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

write.csv(stratified_mean_N_df, "results/indices for assessment/summer flounder/mitigation_1_WEA_distance/MIT.1.ALB.csv", row.names = FALSE)


  ## extract WEE ratio relative to original
original.ALB_df <- read.csv("results/indices for assessment/summer flounder/original.ALB.csv")


ALB_MIT.1_SPRING_ratio <- subset(stratified_mean_N_df, SEASON == "SPRING")$STRATIFIED_MEAN_N/subset(original.ALB_df, SEASON == "SPRING")$STRATIFIED_MEAN_N
ALB_MIT.1_FALL_ratio <- subset(stratified_mean_N_df, SEASON == "FALL")$STRATIFIED_MEAN_N/subset(original.ALB_df, SEASON == "FALL")$STRATIFIED_MEAN_N


save(ALB_MIT.1_SPRING_ratio, ALB_MIT.1_FALL_ratio, file = "results/indices for assessment/summer flounder/mitigation_1_WEA_distance/MIT.1.ALB_ratio.Rdata")

rm(list =ls())

  ## ------------------------ ##




  ## 2.4 MIT.2 dataset ----


# !!! the loaded catch_by_tow data are generated from R script 4.1

MIT.2_cat_df <- read.csv("results/indices for assessment/summer flounder/mitigation_2_model-standardized/Full_ALB_RF_standardized_tow.csv") %>%
  group_by(YEAR, SEASON, STRATUM) %>%
  mutate(TOTAL_N_STATION = length(STATION)) %>%
  mutate(NUMBER = Final_NUMBER) %>% # !! important, need to replace the value 
  filter(YEAR != 2023)

mean_N_stratum_df <- MIT.2_cat_df %>%
  filter(YEAR <= 2008) %>%
  group_by(YEAR, SEASON, STRATUM) %>%
  mutate(MEAN_N_STRATUM = mean(NUMBER)) %>% # mean within a strata
  mutate(VAR_STRATUM = var(NUMBER), na.rm = TRUE) %>% # variance by stratum, same as last line
  ungroup() 

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

write.csv(stratified_mean_N_df, "results/indices for assessment/summer flounder/mitigation_2_model-standardized/MIT.2.ALB.csv", row.names = FALSE)


    ## extract MIT.2 ratio relative to original
original.ALB_df <- read.csv("results/indices for assessment/summer flounder/original.ALB.csv")


ALB_MIT.2_SPRING_ratio <- subset(stratified_mean_N_df, SEASON == "SPRING")$STRATIFIED_MEAN_N/subset(original.ALB_df, SEASON == "SPRING")$STRATIFIED_MEAN_N
ALB_MIT.2_FALL_ratio <- subset(stratified_mean_N_df, SEASON == "FALL")$STRATIFIED_MEAN_N/subset(original.ALB_df, SEASON == "FALL")$STRATIFIED_MEAN_N


save(ALB_MIT.2_SPRING_ratio, ALB_MIT.2_FALL_ratio, file = "results/indices for assessment/summer flounder/mitigation_2_model-standardized/MIT.2.ALB_ratio.Rdata")

rm(list =ls())

  ## ------------------------ ##



  ## 2.5 MIT.3 dataset ----

    # !!! the loaded indices data are generated from R script 5.1

MIT.3.ALB.fall.df <- read.csv("results/indices for assessment/summer flounder/mitigation_3_model-based indices/ALB_fall_Indices.csv") %>%
  add_column(SEASON = "FALL")

MIT.3.ALB.spring.df <- read.csv("results/indices for assessment/summer flounder/mitigation_3_model-based indices/ALB_spring_Indices.csv") %>%
  add_column(SEASON = "SPRING")

stratified_mean_N_df <- rbind(MIT.3.ALB.fall.df, MIT.3.ALB.spring.df)
remove(MIT.3.ALB.fall.df, MIT.3.ALB.spring.df)

write.csv(stratified_mean_N_df, "results/indices for assessment/summer flounder/mitigation_3_model-based indices/MIT.3.ALB.csv", row.names = FALSE)



    ## extract MIT.3 ratio relative to original
original.ALB_df <- read.csv("results/indices for assessment/summer flounder/original.ALB.csv")


ALB_MIT.3_SPRING_ratio <- subset(stratified_mean_N_df, SEASON == "SPRING")$fit/subset(original.ALB_df, SEASON == "SPRING")$STRATIFIED_MEAN_N
ALB_MIT.3_FALL_ratio <- subset(stratified_mean_N_df, SEASON == "FALL")$fit/subset(original.ALB_df, SEASON == "FALL")$STRATIFIED_MEAN_N


save(ALB_MIT.3_SPRING_ratio, ALB_MIT.3_FALL_ratio, file = "results/indices for assessment/summer flounder/mitigation_3_model-based indices/MIT.3.ALB_ratio.Rdata")

rm(list =ls())


  ## ------------------------ ##







# ----------------------------------------------------------------------------------- #




# 3. BIG (2009 - 2022) ----


  ## 3.0 catch by tow data (original) ----

## tow list information
full_tow_df <- read.csv("results/stock assessment/summer flounder/tow data/full.tow.list.for.assessment.csv") %>%
  filter(YEAR > 2008 & YEAR != 2023) %>% # summer flounder assessment uses survey data since 1982
  mutate(STRATUM = as.character(STRATUM)) %>%
  group_by(YEAR, SEASON, STRATUM) %>%
  mutate(TOTAL_N_STATION = length(STATION))


## here we need to do a SWAN following SW66 report P80
caliFactor_df <- read.csv("data/NOAA.stock.data/summer.flounder/rock_efficiencies.csv") [,c(1,2)] %>%
  rename(LENGTH = length,
         CALIBRATION_FACTOR = SSq_cst)


  ## apply the calibration factor at length
load("results/BTS_catch_at_len_by_tow.Rdata")

  ## number by length 
cat_df <- catch_at_len_df %>%
  filter(SVSPP == 103) %>%
  select(-c(CRUISE, SVSPP, CATCHSEX, CATCHSEX)) %>%
  mutate(YEAR = as.numeric(substr(CRUISE6, 1, 4))) %>%
  filter(YEAR > 2008) %>% # summer flounder assessment uses survey data since 1982
  mutate(STRATUM = ifelse(nchar(STRATUM) == 5, str_sub(STRATUM, 2, 5), STRATUM)) %>%
  mutate(ID = paste(CRUISE6, STRATUM, TOW, sep = "."),
         BLOCK = paste(YEAR, SEASON, STRATUM, sep = ".")) %>%
  left_join(caliFactor_df) %>% # integrate with the calibration factors
  mutate(NUMBER_LEN = EXPNUMLEN / CALIBRATION_FACTOR) %>% 
  group_by(YEAR, SEASON, BLOCK, CRUISE6, STRATUM, TOW, STATION, ID) %>%
  summarize(NUMBER = sum(NUMBER_LEN))



## add zero catch tows
cat_df <- cat_df %>%
  right_join(full_tow_df) %>% # right join to add zero catch tows
  mutate(NUMBER = ifelse(is.na(NUMBER), 0, NUMBER)) %>%
  select(YEAR, SEASON, STRATUM, BLOCK, everything()) %>%
  arrange(BLOCK)

remove(caliFactor_df, catch_at_len_df, full_tow_df)


write.csv(cat_df, "results/indices for assessment/summer flounder/catch_by_tow_BIG.csv", row.names = FALSE)

rm(list = ls())
  ## ---------------------------------- ##




  ## 3.1 full dataset ----

mean_N_stratum_df <- read.csv("results/indices for assessment/summer flounder/catch_by_tow_BIG.csv") %>%
  filter(YEAR > 2008) %>%
  group_by(YEAR, SEASON, STRATUM) %>%
  mutate(MEAN_N_STRATUM = mean(NUMBER)) %>% # mean within a strata
  mutate(VAR_STRATUM = var(NUMBER), na.rm = TRUE) %>% # variance by stratum, same as last line
  ungroup() 

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


    ## calculate Swpt Area Numbers adjust by swept area 
      # standard BIG area swept per tow is 0.00647931 sqnm = 22223.41 sqm (times 3429904)
      # survey coverage area: spring 27855, fall 17924

AREA_df <- data.frame(SEASON = c("SPRING", "FALL"), SURVEY.AREA = c(27855, 17924), SWEPT.AREA = 0.00647931)


SWAN <- stratified_mean_N_df %>%
  select(YEAR, SEASON, STRATIFIED_MEAN_N) %>%
  arrange(SEASON) %>%
  left_join(AREA_df) %>%
  mutate(SWAN = STRATIFIED_MEAN_N/SWEPT.AREA * SURVEY.AREA/1000)

write.csv(SWAN, "results/indices for assessment/summer flounder/original.BIG.csv", row.names = FALSE)

rm(list =ls())

  ## ---------------------------------- ##




  ## 3.2 WEE dataset ----

WEE_cat_df <- read.csv("results/indices for assessment/summer flounder/catch_by_tow_BIG.csv")  %>%
  filter(OWF  == "OUTSIDE") %>%
  group_by(YEAR, SEASON, STRATUM) %>%
  mutate(TOTAL_N_STATION = length(STATION)) %>%
  filter(YEAR != 2023)


mean_N_stratum_df <- WEE_cat_df %>%
  filter(YEAR > 2008) %>%
  group_by(YEAR, SEASON, STRATUM) %>%
  mutate(MEAN_N_STRATUM = mean(NUMBER)) %>% # mean within a strata
  mutate(VAR_STRATUM = var(NUMBER), na.rm = TRUE) %>% # variance by stratum, same as last line
  ungroup() 

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


    ## calculate Swpt Area Numbers adjust by swept area 
      # standard BIG area swept per tow is 0.00647931 sqnm = 22223.41 sqm (times 3429904)
      # survey coverage area: spring 27855, fall 17924

AREA_df <- data.frame(SEASON = c("SPRING", "FALL"), SURVEY.AREA = c(27855, 17924), SWEPT.AREA = 0.00647931)


SWAN <- stratified_mean_N_df %>%
  select(YEAR, SEASON, STRATIFIED_MEAN_N) %>%
  arrange(SEASON) %>%
  left_join(AREA_df) %>%
  mutate(SWAN = STRATIFIED_MEAN_N/SWEPT.AREA * SURVEY.AREA/1000)

write.csv(SWAN, "results/indices for assessment/summer flounder/WEE/WEE.BIG.csv", row.names = FALSE)



  ## extract WEE ratio relative to original
original.BIG_df <- read.csv("results/indices for assessment/summer flounder/original.BIG.csv")

BIG_WEE_SPRING_ratio <- subset(SWAN, SEASON == "SPRING")$SWAN/subset(original.BIG_df, SEASON == "SPRING")$SWAN
BIG_WEE_SPRING_ratio[12] = 1 # 2020 not changed


BIG_WEE_FALL_ratio <- subset(SWAN, SEASON == "FALL")$SWAN/subset(original.BIG_df, SEASON == "FALL")$SWAN
BIG_WEE_FALL_ratio[9] = 1 # 2017 not changed
BIG_WEE_FALL_ratio <- append(BIG_WEE_FALL_ratio, 1, after = 11) # 2020 not changed 


save(BIG_WEE_SPRING_ratio, BIG_WEE_FALL_ratio, file = "results/indices for assessment/summer flounder/WEE/WEE.BIG_ratio.Rdata")

rm(list =ls())

  ## ---------------------------------- ##



  ## 3.3 MIT.1 dataset ----


# !!! the loaded catch_by_tow data are generated from R script 3.0


MIT.1_cat_df <- read.csv("results/indices for assessment/summer flounder/mitigation_1_WEA_distance/catch_by_tow_BIG.csv")  %>%
  group_by(YEAR, SEASON, STRATUM) %>%
  mutate(TOTAL_N_STATION = length(STATION)) %>%
  filter(YEAR != 2023)


mean_N_stratum_df <- MIT.1_cat_df %>%
  filter(YEAR > 2008) %>%
  group_by(YEAR, SEASON, STRATUM) %>%
  mutate(MEAN_N_STRATUM = mean(NUMBER)) %>% # mean within a strata
  mutate(VAR_STRATUM = var(NUMBER), na.rm = TRUE) %>% # variance by stratum, same as last line
  ungroup() 

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


  ## calculate Swpt Area Numbers adjust by swept area 
    # standard BIG area swept per tow is 0.00647931 sqnm = 22223.41 sqm (times 3429904)
    # survey coverage area: spring 27855, fall 17924

AREA_df <- data.frame(SEASON = c("SPRING", "FALL"), SURVEY.AREA = c(27855, 17924), SWEPT.AREA = 0.00647931)


SWAN <- stratified_mean_N_df %>%
  select(YEAR, SEASON, STRATIFIED_MEAN_N) %>%
  arrange(SEASON) %>%
  left_join(AREA_df) %>%
  mutate(SWAN = STRATIFIED_MEAN_N/SWEPT.AREA * SURVEY.AREA/1000)

write.csv(SWAN, "results/indices for assessment/summer flounder/mitigation_1_WEA_distance/MIT.1.BIG.csv", row.names = FALSE)



## extract MIT.1 ratio relative to original
original.BIG_df <- read.csv("results/indices for assessment/summer flounder/original.BIG.csv")

BIG_MIT.1_SPRING_ratio <- subset(SWAN, SEASON == "SPRING")$SWAN/subset(original.BIG_df, SEASON == "SPRING")$SWAN
BIG_MIT.1_SPRING_ratio[12] = 1 # 2020 not changed


BIG_MIT.1_FALL_ratio <- subset(SWAN, SEASON == "FALL")$SWAN/subset(original.BIG_df, SEASON == "FALL")$SWAN
BIG_MIT.1_FALL_ratio[9] = 1 # 2017 not changed
BIG_MIT.1_FALL_ratio <- append(BIG_MIT.1_FALL_ratio, 1, after = 11) # 2020 not changed 


save(BIG_MIT.1_SPRING_ratio, BIG_MIT.1_FALL_ratio, file = "results/indices for assessment/summer flounder/mitigation_1_WEA_distance/MIT.1.BIG_ratio.Rdata")

rm(list =ls())

  ## ---------------------------------- ##


  ## 3.4 MIT.2 dataset ----


# !!! the loaded catch_by_tow data are generated from R script 4.1


MIT.2_cat_df <- read.csv("results/indices for assessment/summer flounder/mitigation_2_model-standardized/Full_BIG_RF_standardized_tow.csv")  %>%
  group_by(YEAR, SEASON, STRATUM) %>%
  mutate(TOTAL_N_STATION = length(STATION)) %>%
  mutate(NUMBER = Final_NUMBER) %>% # !! important, need to replace the value 
  filter(YEAR != 2023)


mean_N_stratum_df <- MIT.2_cat_df %>%
  filter(YEAR > 2008) %>%
  group_by(YEAR, SEASON, STRATUM) %>%
  mutate(MEAN_N_STRATUM = mean(NUMBER)) %>% # mean within a strata
  mutate(VAR_STRATUM = var(NUMBER), na.rm = TRUE) %>% # variance by stratum, same as last line
  ungroup() 

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


## calculate Swpt Area Numbers adjust by swept area 
# standard BIG area swept per tow is 0.00647931 sqnm = 22223.41 sqm (times 3429904)
# survey coverage area: spring 27855, fall 17924

AREA_df <- data.frame(SEASON = c("SPRING", "FALL"), SURVEY.AREA = c(27855, 17924), SWEPT.AREA = 0.00647931)


SWAN <- stratified_mean_N_df %>%
  select(YEAR, SEASON, STRATIFIED_MEAN_N) %>%
  arrange(SEASON) %>%
  left_join(AREA_df) %>%
  mutate(SWAN = STRATIFIED_MEAN_N/SWEPT.AREA * SURVEY.AREA/1000)

write.csv(SWAN, "results/indices for assessment/summer flounder/mitigation_2_model-standardized/MIT.2.BIG.csv", row.names = FALSE)



## extract MIT.2 ratio relative to original
original.BIG_df <- read.csv("results/indices for assessment/summer flounder/original.BIG.csv")

BIG_MIT.2_SPRING_ratio <- subset(SWAN, SEASON == "SPRING")$SWAN/subset(original.BIG_df, SEASON == "SPRING")$SWAN
BIG_MIT.2_SPRING_ratio[12] = 1 # 2020 not changed


BIG_MIT.2_FALL_ratio <- subset(SWAN, SEASON == "FALL")$SWAN/subset(original.BIG_df, SEASON == "FALL")$SWAN
BIG_MIT.2_FALL_ratio[9] = 1 # 2017 not changed
BIG_MIT.2_FALL_ratio <- append(BIG_MIT.2_FALL_ratio, 1, after = 11) # 2020 not changed 


save(BIG_MIT.2_SPRING_ratio, BIG_MIT.2_FALL_ratio, file = "results/indices for assessment/summer flounder/mitigation_2_model-standardized/MIT.2.BIG_ratio.Rdata")

rm(list =ls())

  ## ---------------------------------- ##



  ## 3.5 MIT.3 dataset ----


  # !!! the loaded indices data are generated from R script 5.1

MIT.3.BIG.fall.df <- read.csv("results/indices for assessment/summer flounder/mitigation_3_model-based indices/BIG_fall_Indices.csv") %>%
  add_column(SEASON = "FALL") %>%
  filter(YEAR != 2020)

MIT.3.BIG.spring.df <- read.csv("results/indices for assessment/summer flounder/mitigation_3_model-based indices/BIG_spring_Indices.csv") %>%
  add_column(SEASON = "SPRING")

stratified_mean_N_df <- rbind(MIT.3.BIG.fall.df, MIT.3.BIG.spring.df)
remove(MIT.3.BIG.fall.df, MIT.3.BIG.spring.df)


    ## calculate Swpt Area Numbers adjust by swept area 
      # standard BIG area swept per tow is 0.00647931 sqnm = 22223.41 sqm (times 3429904)
      # survey coverage area: spring 27855, fall 17924

AREA_df <- data.frame(SEASON = c("SPRING", "FALL"), SURVEY.AREA = c(27855, 17924), SWEPT.AREA = 0.00647931)


SWAN <- stratified_mean_N_df %>%
  select(YEAR, SEASON, fit) %>%
  arrange(SEASON) %>%
  left_join(AREA_df) %>%
  mutate(SWAN = fit/SWEPT.AREA * SURVEY.AREA/1000)

write.csv(SWAN, "results/indices for assessment/summer flounder/mitigation_3_model-based indices/MIT.3.BIG.csv", row.names = FALSE)



## extract MIT.3 ratio relative to original
original.BIG_df <- read.csv("results/indices for assessment/summer flounder/original.BIG.csv")

BIG_MIT.3_SPRING_ratio <- subset(SWAN, SEASON == "SPRING")$SWAN/subset(original.BIG_df, SEASON == "SPRING")$SWAN
BIG_MIT.3_SPRING_ratio[12] = 1 # 2020 not changed


BIG_MIT.3_FALL_ratio <- subset(SWAN, SEASON == "FALL")$SWAN/subset(original.BIG_df, SEASON == "FALL")$SWAN
BIG_MIT.3_FALL_ratio[9] = 1 # 2017 not changed
BIG_MIT.3_FALL_ratio <- append(BIG_MIT.3_FALL_ratio, 1, after = 11) # 2020 not changed 


save(BIG_MIT.3_SPRING_ratio, BIG_MIT.3_FALL_ratio, file = "results/indices for assessment/summer flounder/mitigation_3_model-based indices/MIT.3.BIG_ratio.Rdata")

rm(list =ls())


