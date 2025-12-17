library(tidyverse)

# ALK for 2009 - 2021
ALK_fall_21_df <- read.csv("assessment model/summer flounder/ALK/al_prop_fluk_fall_2009-2021.csv") %>%
  mutate(YEAR = as.numeric(substr(CRUISE6, 1, 4))) %>%
  select(YEAR, LENGTH, age0_prop_num:age20_prop_num) %>%
  rowwise() %>%
  mutate(age8_prop_num = sum(c_across(age8_prop_num:age20_prop_num), na.rm = TRUE)) %>%
  # mutate(total = sum(c_across(age0_prop_num:age8_prop_num))) %>%
  ungroup() %>%
  select(-c(age9_prop_num:age20_prop_num)) %>%
  add_column(SEASON = "FALL", .before = 1)

ALK_spring_21_df <- read.csv("assessment model/summer flounder/ALK/al_prop_fluk_spring_2009-2021.csv") %>%
  mutate(YEAR = as.numeric(substr(CRUISE6, 1, 4))) %>%
  select(YEAR, LENGTH, age0_prop_num:age20_prop_num) %>%
  rowwise() %>%
  mutate(age8_prop_num = sum(c_across(age8_prop_num:age20_prop_num), na.rm = TRUE)) %>%
  # mutate(total = sum(c_across(age0_prop_num:age8_prop_num))) %>%
  ungroup() %>%
  select(-c(age9_prop_num:age20_prop_num)) %>%
  add_column(SEASON = "SPRING", .before = 1)

ALK_21_df <- rbind(ALK_fall_21_df, ALK_spring_21_df)
write.csv(ALK_21_df, "assessment model/summer flounder/ALK/ALK_before_2021.csv", row.names = FALSE)

remove(ALK_fall_21_df, ALK_spring_21_df, ALK_21_df)

# ALK for 2022
ALK_fall_22_df <- read.csv("assessment model/summer flounder/ALK/al_prop_fluk_fall_2022.csv") %>%
  mutate(YEAR = as.numeric(substr(CRUISE6, 1, 4))) %>%
  select(YEAR, LENGTH, age0_prop_num:age20_prop_num) %>%
  rowwise() %>%
  mutate(age8_prop_num = sum(c_across(age8_prop_num:age20_prop_num), na.rm = TRUE)) %>%
  # mutate(total = sum(c_across(age0_prop_num:age8_prop_num))) %>%
  ungroup() %>%
  select(-c(age9_prop_num:age20_prop_num)) %>%
  add_column(SEASON = "FALL", .before = 1)

ALK_spring_22_df <- read.csv("assessment model/summer flounder/ALK/al_prop_fluk_spring_2022.csv") %>%
  mutate(YEAR = as.numeric(substr(CRUISE6, 1, 4))) %>%
  select(YEAR, LENGTH, age0_prop_num:age20_prop_num) %>%
  rowwise() %>%
  mutate(age8_prop_num = sum(c_across(age8_prop_num:age20_prop_num), na.rm = TRUE)) %>%
  # mutate(total = sum(c_across(age0_prop_num:age8_prop_num))) %>%
  ungroup() %>%
  select(-c(age9_prop_num:age20_prop_num)) %>%
  add_column(SEASON = "SPRING", .before = 1)

ALK_22_df <- rbind(ALK_fall_22_df, ALK_spring_22_df)
write.csv(ALK_22_df, "assessment model/summer flounder/ALK/ALK_2022.csv", row.names = FALSE)

remove(ALK_fall_22_df, ALK_spring_22_df, ALK_22_df)


### below are code for evaluating age effect in summer flounder for offshore wind

# 2.2.2 age-specific abundance numbers from BTS data ----

load("results/BTS_catch_at_len_by_tow.Rdata") # catch by length
caliFactor_df <- read.csv("data/NOAA.stock.data/summer.flounder/STOCKEFF_SV_172735_UNIT_NONE_length_based_calibration.csv")[,c(8:10)] # length-based calibration after 2008, not needed actually


# load ALK
ALK_21_df <- read.csv("assessment model/summer flounder/ALK/ALK_before_2021.csv")
ALK_22_df <- read.csv("assessment model/summer flounder/ALK/ALK_2022.csv")
ALK_df <- rbind(ALK_21_df, ALK_22_df) %>%
  rename_with(~ str_replace(., "age(\\d+)_prop_num", "AGE\\1"), starts_with("age"))

remove(ALK_21_df, ALK_22_df)


# original sample size by year by size
cat_original_df <- catch_at_len_df %>%
  filter(SVSPP == 103) %>%
  select(-c(CRUISE, SVSPP, CATCHSEX, CATCHSEX)) %>%
  mutate(YEAR = as.numeric(substr(CRUISE6, 1, 4))) %>%
  filter(YEAR >= 2009) %>% # summer flounder assessment uses survey data since 1982
  mutate(STRATUM = ifelse(nchar(STRATUM) == 5, str_sub(STRATUM, 2, 5), STRATUM)) %>%
  mutate(ID = paste(CRUISE6, STRATUM, TOW, sep = "."),
         BLOCK = paste(YEAR, SEASON, STRATUM, sep = ".")) %>%
  group_by(YEAR, SEASON, LENGTH) %>%
  summarize(NUMBER = sum(EXPNUMLEN)) %>%
  mutate(GROUP = ifelse(YEAR < 2022, 2021, 2022))
# left_join(caliFactor_df) %>% # integrate with the calibration factors
# mutate(NUMBER_LEN = EXPNUMLEN / CALIBRATION_FACTOR) %>% # calculate the numbers
# group_by(YEAR, SEASON, BLOCK, CRUISE6, STRATUM, TOW, STATION, ID, LENGTH) %>%
# summarize(NUMBER = sum(EXPNUMLEN))

cat_original_by_age_df <- cat_original_df %>%
  # filter(LENGTH >= min(ALK_df$LENGTH)) %>% 
  right_join(ALK_df)




# WEA sample size by year by size

full_tow_df <- read.csv("results/stratified.mean.indices/summer.flounder/tow.list.csv") # full tow list
BTS_overlay_df <- read.csv("results/BTS_tows_overlay_final.csv") # impacted tow list

WEA_tow_df <- full_tow_df %>%
  filter(!ID %in% unique(BTS_overlay_df$ID))

cat_WEA_df <- catch_at_len_df %>%
  filter(SVSPP == 103, !ID %in% unique(BTS_overlay_df$ID)) %>%
  select(-c(CRUISE, SVSPP, CATCHSEX, CATCHSEX)) %>%
  mutate(YEAR = as.numeric(substr(CRUISE6, 1, 4))) %>%
  filter(YEAR >= 2009) %>% # summer flounder assessment uses survey data since 1982
  mutate(STRATUM = ifelse(nchar(STRATUM) == 5, str_sub(STRATUM, 2, 5), STRATUM)) %>%
  mutate(ID = paste(CRUISE6, STRATUM, TOW, sep = "."),
         BLOCK = paste(YEAR, SEASON, STRATUM, sep = ".")) %>%
  left_join(caliFactor_df) %>% # integrate with the calibration factors
  mutate(NUMBER_LEN = EXPNUMLEN / CALIBRATION_FACTOR) %>% # calculate the numbers
  group_by(YEAR, SEASON, BLOCK, CRUISE6, STRATUM, TOW, STATION, ID, LENGTH) %>%
  summarize(NUMBER = sum(NUMBER_LEN)) %>%
  mutate(GROUP = ifelse(YEAR < 2022, 2021, 2022))
