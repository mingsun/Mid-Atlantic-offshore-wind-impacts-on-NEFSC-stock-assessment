library(tidyverse)

# this mitigation scenario will reallocate the lost survey effort to the tows that are closer to the WEAs
# this reallocation is handled for each year for the top ranking tows based on lost effort



# 1. summer flounder ----

  ## 1.1 generate tow list ----

 # load tow list with effort,WEA distance, and ID
tow_list_df <- read.csv("results/Tow distance from WEAs/SFL_Rank_Final.csv") %>%
  arrange(SEASON, YEAR, NEW_RANK)


  # prepare a empty df for mitigated tow list
mitigtated_tow_list_df <- data.frame()


  # loop to generate a full catch database with reallocated effort
# s = "SPRING"; i = 1982

for (s in unique(tow_list_df$SEASON)) {
  
  for (i in min(tow_list_df$YEAR):max(tow_list_df$YEAR)) {
    
    if(i == 2020 && s == "FALL") next() # 2020 has no fall survey so escape it 
    
    temp_tow_list <- subset(tow_list_df, SEASON == s & YEAR == i)
    
    # determine lost tows
    N_TOW_lost <- sum(temp_tow_list$NEW_RANK == 0)
    
    # pick the same amount of tows based on distance to WEA (closest)
    temp_tow_added <- temp_tow_list %>% 
      filter(NEW_RANK != 0) %>% 
      arrange(NEW_RANK) %>%  
      slice_head(n = N_TOW_lost)
    
    # combine them
    temp_tow_mitigtated_list <- temp_tow_list %>% 
      filter(NEW_RANK != 0) %>% 
      bind_rows(temp_tow_added)
      
    # combine to the full database
    mitigtated_tow_list_df <- rbind(mitigtated_tow_list_df, temp_tow_mitigtated_list)
    
    remove(temp_tow_list, N_TOW_lost, temp_tow_added, temp_tow_mitigtated_list)
  }
  
}; remove(i,s)




  # load catch data to generate the new catch database

cat_ALB_df <- read.csv("results/indices for assessment/summer flounder/catch_by_tow_ALB.csv") 
cat_BIG_df <- read.csv("results/indices for assessment/summer flounder/catch_by_tow_BIG.csv") 

mitigtated_tow_list_ALB_df <- subset(mitigtated_tow_list_df, YEAR <= 2008)
mitigtated_tow_list_BIG_df <- subset(mitigtated_tow_list_df, YEAR > 2008 & YEAR != 2023)



  # generate the catch database based on mitigation strategies

cat_mitigtated_ALB_df <- cat_ALB_df[match(mitigtated_tow_list_ALB_df$ID, cat_ALB_df$ID), ]
rownames(cat_mitigtated_ALB_df) <- 1:nrow(cat_mitigtated_ALB_df)

cat_mitigtated_BIG_df <- cat_BIG_df[match(mitigtated_tow_list_BIG_df$ID, cat_BIG_df$ID), ]
rownames(cat_mitigtated_BIG_df) <- 1:nrow(cat_mitigtated_BIG_df)


write.csv(mitigtated_tow_list_df, "results/indices for assessment/summer flounder/mitigation_1_WEA_distance/tow_list.csv", row.names = FALSE)
write.csv(mitigtated_tow_list_ALB_df, "results/indices for assessment/summer flounder/mitigation_1_WEA_distance/ALB_tow_list.csv", row.names = FALSE)
write.csv(mitigtated_tow_list_BIG_df, "results/indices for assessment/summer flounder/mitigation_1_WEA_distance/BIG_tow_list.csv", row.names = FALSE)


write.csv(cat_mitigtated_ALB_df,      "results/indices for assessment/summer flounder/mitigation_1_WEA_distance/catch_by_tow_ALB.csv", row.names = FALSE)
write.csv(cat_mitigtated_BIG_df,      "results/indices for assessment/summer flounder/mitigation_1_WEA_distance/catch_by_tow_BIG.csv", row.names = FALSE)

remove(cat_df, tow_list_df)






  ## 1.2 generate abundance indices ----


  ### 1.2.1 mean numbers by stratum ---------------------------------------------------------------------------------------------
# all methods below follow (https://noaa-edab.github.io/survdat/articles/calc_strat_mean.html)

mean_N_stratum_df <- cat_mitigtated_df %>%
  group_by(YEAR, SEASON, STRATUM) %>%
  mutate(MEAN_N_STRATUM = mean(NUMBER)) %>% # mean within a strata
  mutate(VAR_STRATUM = sum((NUMBER - MEAN_N_STRATUM)^2, na.rm = TRUE)/ (TOTAL_N_STATION - 1)) %>% # variance by stratum
  # mutate(VAR_STRATUM = var(EXPCATCHNUM)) %>% # variance by stratum, same as last line
  ungroup() 



  ### 1.2.2 stratified mean numbers ---------------------------------------------------------------------------------------------

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

write.csv(stratified_mean_N_df, "results/indices for assessment/summer flounder/mitigation.1.indices.csv", row.names = FALSE)



# ----------------------------------------------------- #











# 2. squid ----


## 2.1 generate tow list ----


  # load tow list with effort,WEA distance, and ID
tow_list_df <- read.csv("results/Tow distance from WEAs/LFS_Rank_Final.csv") %>%
  arrange(SEASON, YEAR, NEW_RANK)


  # prepare a empty df for mitigated tow list
mitigtated_tow_list_df <- data.frame()


# loop to generate a full catch database with reallocated effort
# s = "SPRING"; i = 1982

for (s in unique(tow_list_df$SEASON)) {
  
  for (i in min(tow_list_df$YEAR):max(tow_list_df$YEAR)) {
    
    if(i == 2020 && s == "FALL") next() # 2020 has no fall survey so escape it 
    
    temp_tow_list <- subset(tow_list_df, SEASON == s & YEAR == i)
    
    # determine lost tows
    N_TOW_lost <- sum(temp_tow_list$NEW_RANK == 0)
    
    # pick the same amount of tows based on distance to WEA (closest)
    temp_tow_added <- temp_tow_list %>% 
      filter(NEW_RANK != 0) %>% 
      arrange(NEW_RANK) %>%  
      slice_head(n = N_TOW_lost)
    
    # combine them
    temp_tow_mitigtated_list <- temp_tow_list %>% 
      filter(NEW_RANK != 0) %>% 
      bind_rows(temp_tow_added)
    
    # combine to the full database
    mitigtated_tow_list_df <- rbind(mitigtated_tow_list_df, temp_tow_mitigtated_list)
    
    remove(temp_tow_list, N_TOW_lost, temp_tow_added, temp_tow_mitigtated_list)
  }
  
}; remove(i,s)




# load catch data to generate the new catch database
cat_df <- tow_list_df  %>% # the tow list already has catch weight info
  group_by(YEAR, SEASON, STRATUM) %>%
  mutate(TOTAL_N_STATION = length(STATION)) %>% # this is important because number of tows changed due to tow filtering
  ungroup() 



# generate the catch data base based on mitigation strategies
cat_mitigtated_df <- cat_df[match(mitigtated_tow_list_df$ID, cat_df$ID), ]
rownames(cat_mitigtated_df) <- 1:nrow(cat_mitigtated_df)


write.csv(mitigtated_tow_list_df, "results/indices for assessment/squid/mitigation_1_WEA_distance/tow_list.csv", row.names = FALSE)
write.csv(cat_mitigtated_df,      "results/indices for assessment/squid/mitigation_1_WEA_distance/catch_by_tow.csv", row.names = FALSE)

remove(cat_df, tow_list_df)



  ## 2.2 generate abundance indices ----




    ### 2.2.1 mean numbers by stratum ---------------------------------------------------------------------------------------------
      # all methods below follow (https://noaa-edab.github.io/survdat/articles/calc_strat_mean.html)

mean_N_stratum_df <- cat_mitigtated_df %>%
  group_by(YEAR, SEASON, STRATUM) %>%
  mutate(MEAN_N_STRATUM = mean(CATCH_WT_CAL, na.rm = TRUE)) %>% # mean within a strata
  mutate(VAR_STRATUM = sum((CATCH_WT_CAL - MEAN_N_STRATUM)^2, na.rm = TRUE)/ (TOTAL_N_STATION - 1)) %>% # variance by stratum
  ungroup() 



    ### 2.2.2 stratified mean numbers ---------------------------------------------------------------------------------------------

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

write.csv(stratified_mean_N_df, "results/indices for assessment/squid/mitigation.1.indices.csv", row.names = FALSE)




# ----------------------------------------------------- #






# 3. surfclam ----

  ## 3.1 generate tow list ----

# load tow list with effort,WEA distance, and ID
tow_list_df <- read.csv("results/Tow distance from WEAs/ASC_Rank_Final.csv") %>%
  arrange(YEAR, NEW_RANK)


# prepare a empty df for mitigated tow list
mitigtated_tow_list_df <- data.frame()


# loop to generate a full catch database with reallocated effort
# s = "SPRING"; i = 1982


for (i in min(tow_list_df$YEAR):max(tow_list_df$YEAR)) {
  
  # if(i == 2020 && s == "FALL") next() # 2020 has no fall survey so escape it 
  
  temp_tow_list <- subset(tow_list_df, YEAR == i)
  
  # determine lost tows
  N_TOW_lost <- sum(temp_tow_list$NEW_RANK == 0)
  
  # pick the same amount of tows based on distance to WEA (closest)
  temp_tow_added <- temp_tow_list %>% 
    filter(NEW_RANK != 0) %>% 
    arrange(NEW_RANK) %>%  
    slice_head(n = N_TOW_lost)
  
  # combine them
  temp_tow_mitigtated_list <- temp_tow_list %>% 
    filter(NEW_RANK != 0) %>% 
    bind_rows(temp_tow_added)
  
  # combine to the full database
  mitigtated_tow_list_df <- rbind(mitigtated_tow_list_df, temp_tow_mitigtated_list)
  
  remove(temp_tow_list, N_TOW_lost, temp_tow_added, temp_tow_mitigtated_list)
}

  
remove(i)




# load catch data to generate the new catch database
cat_df <- read.csv("results/stratified.mean.indices/surfclam/total.catch.by.tow.csv",)

# generate the catch data base based on mitigation strategies
cat_mitigtated_df <- cat_df[match(mitigtated_tow_list_df$ID, cat_df$ID), ]
rownames(cat_mitigtated_df) <- 1:nrow(cat_mitigtated_df)


write.csv(mitigtated_tow_list_df, "results/indices for assessment/surfclam/mitigation_1_WEA_distance/tow_list.csv", row.names = FALSE)
write.csv(cat_mitigtated_df,      "results/indices for assessment/surfclam/mitigation_1_WEA_distance/catch_by_tow.csv", row.names = FALSE)

remove(cat_df, tow_list_df)


## 3.2 generate abundance indices ----


  ### 3.2.1 calculate mean biomass by stratum ---------------------------------------------------------------------------------------------

mean_N_stratum_df <- cat_mitigtated_df %>%
  group_by(YEAR, REGION, STRATUM) %>%
  mutate(MEAN_N_STRATUM = mean(NPERTOW)) %>% # mean within a strata
  mutate(VAR_STRATUM = sum((NPERTOW - MEAN_N_STRATUM)^2, na.rm = TRUE)/ (TOTAL_N_STATION - 1)) %>% # variance by stratum
  mutate(SE_STRATUM = sqrt(VAR_STRATUM)/sqrt(TOTAL_N_STATION)) %>%
  # mutate(VAR_STRATUM = var(EXPCATCHNUM)) %>% # variance by stratum, same as last line
  ungroup() %>% 
  select(c(YEAR, REGION, STRATUM, TOTAL_N_STATION, WEIGHT, MEAN_N_STRATUM, VAR_STRATUM)) %>%
  distinct() %>%
  arrange(YEAR, REGION)


  ### 3.2.2 calculate stratified mean biomass ---------------------------------------------------------------------------------------------

stratified_mean_N_df <- mean_N_stratum_df %>%
  # select(c(YEAR, REGION, STRATUM, TOTAL_N_STATION, WEIGHT, MEAN_N_STRATUM, VAR_STRATUM)) %>%
  group_by(YEAR, REGION) %>%
  summarize(STRATIFIED_MEAN_N = weighted.mean(MEAN_N_STRATUM, w = WEIGHT), # stratified mean
            STRATIFIED_VAR = sum(WEIGHT^2 * VAR_STRATUM / TOTAL_N_STATION, na.rm = TRUE), # variance
            STRATIFIED_SE = sqrt(STRATIFIED_VAR)) %>%   # standard deviance
  mutate(up_CI_95 = qnorm(0.975, mean = STRATIFIED_MEAN_N, sd = STRATIFIED_SE),
         lo_CI_95 = qnorm(0.025, mean = STRATIFIED_MEAN_N, sd = STRATIFIED_SE),
         lo_CI_95 = ifelse(lo_CI_95 < 0, yes = 0, no = lo_CI_95),
         CV = STRATIFIED_SE / STRATIFIED_MEAN_N) %>%
  ungroup()

write.csv(stratified_mean_B_df, "results/indices for assessment/surfclam/mitigation.1.indices.csv", row.names = FALSE)




# ----------------------------------------------------- #







# 4. quahog ----

  ## 4.1 generate tow list ----

    # load tow list with effort,WEA distance, and ID
tow_list_df <- read.csv("results/Tow distance from WEAs/OQ_Rank_Final.csv") %>%
  arrange(YEAR, NEW_RANK)


# prepare a empty df for mitigated tow list
mitigtated_tow_list_df <- data.frame()


# loop to generate a full catch database with reallocated effort
# s = "SPRING"; i = 1982


for (i in min(tow_list_df$YEAR):max(tow_list_df$YEAR)) {
  
  # if(i == 2020 && s == "FALL") next() # 2020 has no fall survey so escape it 
  
  temp_tow_list <- subset(tow_list_df, YEAR == i)
  
  # determine lost tows
  N_TOW_lost <- sum(temp_tow_list$NEW_RANK == 0)
  
  # pick the same amount of tows based on distance to WEA (closest)
  temp_tow_added <- temp_tow_list %>% 
    filter(NEW_RANK != 0) %>% 
    arrange(NEW_RANK) %>%  
    slice_head(n = N_TOW_lost)
  
  # combine them
  temp_tow_mitigtated_list <- temp_tow_list %>% 
    filter(NEW_RANK != 0) %>% 
    bind_rows(temp_tow_added)
  
  # combine to the full database
  mitigtated_tow_list_df <- rbind(mitigtated_tow_list_df, temp_tow_mitigtated_list)
  
  remove(temp_tow_list, N_TOW_lost, temp_tow_added, temp_tow_mitigtated_list)
}


remove(i)




# load catch data to generate the new catch database
cat_df <- read.csv("results/stratified.mean.indices/quahog//total.catch.by.tow.csv",)

# generate the catch data base based on mitigation strategies
cat_mitigtated_df <- cat_df[match(mitigtated_tow_list_df$ID, cat_df$ID), ]
rownames(cat_mitigtated_df) <- 1:nrow(cat_mitigtated_df)


write.csv(mitigtated_tow_list_df, "results/indices for assessment/quahog/mitigation_1_WEA_distance/tow_list.csv", row.names = FALSE)
write.csv(cat_mitigtated_df,      "results/indices for assessment/quahog/mitigation_1_WEA_distance/catch_by_tow.csv", row.names = FALSE)

remove(cat_df, tow_list_df)


  ## 4.2 generate abundance indices ----


    ### 4.2.1 calculate mean biomass by stratum ---------------------------------------------------------------------------------------------

mean_N_stratum_df <- cat_mitigtated_df %>%
  group_by(YEAR, REGION, STRATUM) %>%
  mutate(MEAN_N_STRATUM = mean(NPERTOW)) %>% # mean within a strata
  mutate(VAR_STRATUM = sum((NPERTOW - MEAN_N_STRATUM)^2, na.rm = TRUE)/ (TOTAL_N_STATION - 1)) %>% # variance by stratum
  mutate(SE_STRATUM = sqrt(VAR_STRATUM)/sqrt(TOTAL_N_STATION)) %>%
  # mutate(VAR_STRATUM = var(EXPCATCHNUM)) %>% # variance by stratum, same as last line
  ungroup() %>% 
  select(c(YEAR, REGION, STRATUM, TOTAL_N_STATION, WEIGHT, MEAN_N_STRATUM, VAR_STRATUM)) %>%
  distinct() %>%
  arrange(YEAR, REGION)


    ### 4.2.2 calculate stratified mean biomass ---------------------------------------------------------------------------------------------

stratified_mean_N_df <- mean_N_stratum_df %>%
  # select(c(YEAR, REGION, STRATUM, TOTAL_N_STATION, WEIGHT, MEAN_N_STRATUM, VAR_STRATUM)) %>%
  group_by(YEAR, REGION) %>%
  summarize(STRATIFIED_MEAN_N = weighted.mean(MEAN_N_STRATUM, w = WEIGHT), # stratified mean
            STRATIFIED_VAR = sum(WEIGHT^2 * VAR_STRATUM / TOTAL_N_STATION, na.rm = TRUE), # variance
            STRATIFIED_SE = sqrt(STRATIFIED_VAR)) %>%   # standard deviance
  mutate(up_CI_95 = qnorm(0.975, mean = STRATIFIED_MEAN_N, sd = STRATIFIED_SE),
         lo_CI_95 = qnorm(0.025, mean = STRATIFIED_MEAN_N, sd = STRATIFIED_SE),
         lo_CI_95 = ifelse(lo_CI_95 < 0, yes = 0, no = lo_CI_95),
         CV = STRATIFIED_SE / STRATIFIED_MEAN_N) %>%
  ungroup()

write.csv(stratified_mean_N_df, "results/indices for assessment/quahog/mitigation.1.indices.csv", row.names = FALSE)

# ----------------------------------------------------- #






