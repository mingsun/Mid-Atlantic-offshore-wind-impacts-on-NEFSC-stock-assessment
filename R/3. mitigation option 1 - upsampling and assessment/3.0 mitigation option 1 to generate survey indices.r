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
# added tows are flagged as MIT_ADDED = 1
# s = "SPRING"; i = 1982

for (s in unique(tow_list_df$SEASON)) {
  
  for (i in min(tow_list_df$YEAR):max(tow_list_df$YEAR)) {
    
    if(i == 2020 && s == "FALL") next() # 2020 has no fall survey so escape it 
    
    temp_tow_list <- subset(tow_list_df, SEASON == s & YEAR == i)
    
    # determine lost tows
    N_TOW_lost <- sum(temp_tow_list$NEW_RANK == 0)
    
    # retain all accessible tows
    temp_tow_retained <- temp_tow_list %>%
      filter(NEW_RANK != 0) %>%
      mutate(MIT_ADDED = 0)
    
    # pick the same amount of tows based on distance to WEA (closest)
    temp_tow_added <- temp_tow_list %>% 
      filter(NEW_RANK != 0) %>% 
      arrange(NEW_RANK) %>%  
      slice_head(n = N_TOW_lost) %>%
      mutate(MIT_ADDED = 1)
    
    # combine retained observations and added pseudo-observations
    temp_tow_mitigtated_list <- bind_rows(temp_tow_retained, temp_tow_added)
    
    # combine to the full database
    mitigtated_tow_list_df <- rbind(mitigtated_tow_list_df, temp_tow_mitigtated_list)
    
    remove(temp_tow_list, N_TOW_lost, temp_tow_retained, temp_tow_added, temp_tow_mitigtated_list)
  }
  
}; remove(i,s)




  # load catch data to generate the new catch database

cat_ALB_df <- read.csv("results/indices for assessment/summer flounder/catch_by_tow_ALB.csv") 
cat_BIG_df <- read.csv("results/indices for assessment/summer flounder/catch_by_tow_BIG.csv") 

mitigtated_tow_list_ALB_df <- subset(mitigtated_tow_list_df, YEAR <= 2008)
mitigtated_tow_list_BIG_df <- subset(mitigtated_tow_list_df, YEAR > 2008 & YEAR != 2023)


  # load year x season x stratum variance estimates

variance_ALB_df <- read.csv("results/indices for assessment/summer flounder/stratum_year.ALB.csv") %>%
  select(YEAR, SEASON, STRATUM, VAR_STRATUM) %>%
  distinct()

variance_BIG_df <- read.csv("results/indices for assessment/summer flounder/stratum_year.BIG.csv") %>%
  select(YEAR, SEASON, STRATUM, VAR_STRATUM) %>%
  distinct()


  # separate ALB and BIG periods

mitigtated_tow_list_ALB_df <- subset(mitigtated_tow_list_df, YEAR <= 2008)
mitigtated_tow_list_BIG_df <- subset(mitigtated_tow_list_df, YEAR > 2008 & YEAR != 2023)


  # generate the catch database based on mitigation strategies

cat_mitigtated_ALB_df <- cat_ALB_df[match(mitigtated_tow_list_ALB_df$ID, cat_ALB_df$ID), ]
rownames(cat_mitigtated_ALB_df) <- 1:nrow(cat_mitigtated_ALB_df)

cat_mitigtated_BIG_df <- cat_BIG_df[match(mitigtated_tow_list_BIG_df$ID, cat_BIG_df$ID), ]
rownames(cat_mitigtated_BIG_df) <- 1:nrow(cat_mitigtated_BIG_df)


  # carry the indicator of whether each row is an original observation
  # or an added pseudo observation

cat_mitigtated_ALB_df$MIT_ADDED <- mitigtated_tow_list_ALB_df$MIT_ADDED
cat_mitigtated_BIG_df$MIT_ADDED <- mitigtated_tow_list_BIG_df$MIT_ADDED


## 1.2 add stratum-level error to the pseudo observation ----

set.seed(1234)

# ALB

cat_mitigtated_ALB_df <- cat_mitigtated_ALB_df %>%
  left_join(variance_ALB_df,by = c("YEAR", "SEASON", "STRATUM")) %>% 
  mutate(VAR_STRATUM = ifelse(is.na(VAR_STRATUM), 0, VAR_STRATUM))

cat_mitigtated_ALB_df <- cat_mitigtated_ALB_df %>%
  mutate(MIT_ERROR = ifelse(
    MIT_ADDED == 1,
    rnorm(
      n(),
      mean = 0,
      sd = sqrt(VAR_STRATUM)
    ),
    0
  ),
  # add error only to replacement observations and constrain abundance to be non-negative
  NUMBER = pmax(0, NUMBER + MIT_ERROR)
  )


# BIG

cat_mitigtated_BIG_df <- cat_mitigtated_BIG_df %>%
  left_join(variance_BIG_df, by = c("YEAR", "SEASON", "STRATUM")) %>% 
  mutate(VAR_STRATUM = ifelse(is.na(VAR_STRATUM), 0, VAR_STRATUM))

cat_mitigtated_BIG_df <- cat_mitigtated_BIG_df %>%
  mutate(MIT_ERROR = ifelse(
      MIT_ADDED == 1,
      rnorm(
        n(),
        mean = 0,
        sd = sqrt(VAR_STRATUM)
      ),
      0
    ),
    # add error only to replacement observations and constrain abundance to be non-negative
    NUMBER = pmax(0, NUMBER + MIT_ERROR)
  )


# remove temporary variables so final catch databases retain original format

cat_mitigtated_ALB_df <- cat_mitigtated_ALB_df %>%
  select(-MIT_ADDED, -VAR_STRATUM, -MIT_ERROR)

cat_mitigtated_BIG_df <- cat_mitigtated_BIG_df %>%
  select(-MIT_ADDED, -VAR_STRATUM, -MIT_ERROR)



## 1.3 save the tow information ----

write.csv(mitigtated_tow_list_df, "results/indices for assessment/summer flounder/mitigation_1_WEA_distance/tow_list.csv", row.names = FALSE)
write.csv(mitigtated_tow_list_ALB_df, "results/indices for assessment/summer flounder/mitigation_1_WEA_distance/ALB_tow_list.csv", row.names = FALSE)
write.csv(mitigtated_tow_list_BIG_df, "results/indices for assessment/summer flounder/mitigation_1_WEA_distance/BIG_tow_list.csv", row.names = FALSE)


write.csv(cat_mitigtated_ALB_df,      "results/indices for assessment/summer flounder/mitigation_1_WEA_distance/catch_by_tow_ALB.csv", row.names = FALSE)
write.csv(cat_mitigtated_BIG_df,      "results/indices for assessment/summer flounder/mitigation_1_WEA_distance/catch_by_tow_BIG.csv", row.names = FALSE)

remove(cat_df, tow_list_df)


# ----------------------------------------------------- #


## 1.4 generate abundance indices ----


  # !!!! go to script 1.7.2 for calculating MIT.1 indices


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
    
    # retain all accessible tows
    temp_tow_retained <- temp_tow_list %>%
      filter(NEW_RANK != 0) %>%
      mutate(MIT_ADDED = 0)
    
    # pick the same amount of tows based on distance to WEA (closest)
    temp_tow_added <- temp_tow_list %>% 
      filter(NEW_RANK != 0) %>% 
      arrange(NEW_RANK) %>%  
      slice_head(n = N_TOW_lost) %>%
      mutate(MIT_ADDED = 1)
    
    # combine retained observations and added pseudo-observations
    temp_tow_mitigtated_list <- bind_rows(temp_tow_retained, temp_tow_added)
    
    # combine to the full database
    mitigtated_tow_list_df <- rbind(mitigtated_tow_list_df, temp_tow_mitigtated_list)
    
    remove(temp_tow_list, N_TOW_lost, temp_tow_retained, temp_tow_added, temp_tow_mitigtated_list)
  }
  
}; remove(i,s)




# load catch data to generate the new catch database
cat_df <- tow_list_df  %>% # the tow list already has catch weight info
  group_by(YEAR, SEASON, STRATUM) %>%
  mutate(TOTAL_N_STATION = length(STATION)) %>% # this is important because number of tows changed due to tow filtering
  ungroup() 


# load year x season x stratum variance estimates
variance_df <- read.csv("results/indices for assessment/squid/stratum_year.csv") %>%
  select(YEAR, SEASON, STRATUM, VAR_STRATUM) %>%
  distinct()


# generate the catch data base based on mitigation strategies
cat_mitigtated_df <- cat_df[match(mitigtated_tow_list_df$ID, cat_df$ID), ]
rownames(cat_mitigtated_df) <- 1:nrow(cat_mitigtated_df)
cat_mitigtated_df$MIT_ADDED <- mitigtated_tow_list_df$MIT_ADDED



##  2.2 add stratum-level error to the pseudo observations ----

cat_mitigtated_df <- cat_mitigtated_df %>%
  left_join(variance_df, by = c("YEAR", "SEASON", "STRATUM")) %>%
  mutate(
    # if stratum variance is unavailable, assign variance = 0
    VAR_STRATUM = ifelse(
      is.na(VAR_STRATUM),
      0,
      VAR_STRATUM
    ),
    # initialize error as zero for all observations
    MIT_ERROR = 0
  )


# fixed seed for reproducibility
set.seed(1234)


# identify only the added pseudo-observations
added_rows <- which(cat_mitigtated_df$MIT_ADDED == 1)


# generate one error draw for each added tow using its corresponding
# year x season x stratum variance
cat_mitigtated_df$MIT_ERROR[added_rows] <- rnorm(length(added_rows), mean = 0,
                                                 sd = sqrt(cat_mitigtated_df$VAR_STRATUM[added_rows])
)


# adjust catch weight only for added pseudo-observations
# and constrain the minimum value to zero
cat_mitigtated_df$CATCH_WT_CAL <- pmax(0, cat_mitigtated_df$CATCH_WT_CAL + cat_mitigtated_df$MIT_ERROR)

# remove temporary variables so final catch database retains original format
cat_mitigtated_df <- cat_mitigtated_df %>%
  select(-MIT_ADDED, -VAR_STRATUM, -MIT_ERROR)


  ## 2.3 save the tow information ----

write.csv(mitigtated_tow_list_df, "results/indices for assessment/squid/mitigation_1_WEA_distance/tow_list.csv", row.names = FALSE)
write.csv(cat_mitigtated_df,      "results/indices for assessment/squid/mitigation_1_WEA_distance/catch_by_tow.csv", row.names = FALSE)

remove(cat_df, tow_list_df)



  ## 2.4 generate abundance indices ----


    ### 2.4.1 mean numbers by stratum ---------------------------------------------------------------------------------------------
      # all methods below follow (https://noaa-edab.github.io/survdat/articles/calc_strat_mean.html)

mean_N_stratum_df <- cat_mitigtated_df %>%
  group_by(YEAR, SEASON, STRATUM) %>%
  mutate(MEAN_N_STRATUM = mean(CATCH_WT_CAL, na.rm = TRUE)) %>% # mean within a strata
  mutate(VAR_STRATUM = sum((CATCH_WT_CAL - MEAN_N_STRATUM)^2, na.rm = TRUE)/ (TOTAL_N_STATION - 1)) %>% # variance by stratum
  ungroup() 



    ### 2.4.2 stratified mean numbers ---------------------------------------------------------------------------------------------

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
  
  temp_tow_list <- subset(tow_list_df, YEAR == i)
  
  # determine lost tows
  N_TOW_lost <- sum(temp_tow_list$NEW_RANK == 0)
  
    # retain accessible tows
  temp_tow_retained <- temp_tow_list %>%
    filter(NEW_RANK != 0) %>%
    mutate(MIT_ADDED = 0)
  
  # pick the same amount of tows based on distance to WEA (closest)
  temp_tow_added <- temp_tow_list %>%
    filter(NEW_RANK != 0) %>%
    arrange(NEW_RANK) %>%
    slice_head(n = N_TOW_lost) %>%
    mutate(MIT_ADDED = 1)
  
  # combine retained and added tows
  temp_tow_mitigtated_list <- bind_rows(temp_tow_retained, temp_tow_added)
  
  
  # combine to the full database
  mitigtated_tow_list_df <- bind_rows(mitigtated_tow_list_df, temp_tow_mitigtated_list)
  
  remove(temp_tow_list, N_TOW_lost, temp_tow_retained, temp_tow_added, temp_tow_mitigtated_list)
}

remove(i)



  ## 3.2 add stratum-level observation error ----

variance_RDtrendS_df <- read.csv("results/indices for assessment/surfclam/stratum_year.RDtrendS.csv") %>%
  select(YEAR, STRATUM, VAR_STRATUM) %>%
  distinct() %>%
  rename(VAR_RDtrendS = VAR_STRATUM)


variance_RDscaleS_df <- read.csv("results/indices for assessment/surfclam/stratum_year.RDscaleS.csv") %>%
  select(YEAR, STRATUM, VAR_STRATUM) %>%
  distinct() %>%
  rename(VAR_RDscaleS = VAR_STRATUM)


variance_MCDS_df <- read.csv("results/indices for assessment/surfclam/stratum_year.MCDS.csv") %>%
  select(YEAR, STRATUM, VAR_STRATUM) %>%
  distinct() %>%
  rename(VAR_MCDS = VAR_STRATUM)


# attach the three variance estimates to the mitigation tow list

mitigtated_tow_list_df <- mitigtated_tow_list_df %>%
  left_join(variance_RDtrendS_df, by = c("YEAR", "STRATUM")) %>%
  left_join(variance_RDscaleS_df, by = c("YEAR", "STRATUM")) %>%
  left_join(variance_MCDS_df,by = c("YEAR", "STRATUM")) %>%
  mutate(
    # assign variance = 0 where a year x stratum variance is unavailable
    VAR_RDtrendS = ifelse(is.na(VAR_RDtrendS), 0, VAR_RDtrendS),
    VAR_RDscaleS = ifelse(is.na(VAR_RDscaleS), 0, VAR_RDscaleS),
    VAR_MCDS     = ifelse(is.na(VAR_MCDS),     0, VAR_MCDS),
    # initialize mitigation errors
    MIT_ERROR_RDtrendS = 0,
    MIT_ERROR_RDscaleS = 0,
    MIT_ERROR_MCDS     = 0
  )


# fixed seed for reproducibility
set.seed(1234)

# identify added pseudo-observations
added_rows <- which(mitigtated_tow_list_df$MIT_ADDED == 1)


# generate separate observation errors for each abundance time series

mitigtated_tow_list_df$MIT_ERROR_RDtrendS[added_rows] <- rnorm(
  length(added_rows),
  mean = 0,
  sd = sqrt(mitigtated_tow_list_df$VAR_RDtrendS[added_rows])
)

mitigtated_tow_list_df$MIT_ERROR_RDscaleS[added_rows] <- rnorm(
  length(added_rows),
  mean = 0,
  sd = sqrt(mitigtated_tow_list_df$VAR_RDscaleS[added_rows])
)

mitigtated_tow_list_df$MIT_ERROR_MCDS[added_rows] <- rnorm(
  length(added_rows),
  mean = 0,
  sd = sqrt(mitigtated_tow_list_df$VAR_MCDS[added_rows])
)


# generate catch database for RDscaleS and MCDS

cat_df <- read.csv("results/stratified.mean.indices/surfclam/total.catch.by.tow.csv")


# generate catch database based on mitigation strategy
cat_mitigtated_df <- cat_df[match(mitigtated_tow_list_df$ID, cat_df$ID), ]

rownames(cat_mitigtated_df) <- 1:nrow(cat_mitigtated_df)


# carry mitigation information into catch database
cat_mitigtated_df$MIT_ADDED <- mitigtated_tow_list_df$MIT_ADDED
cat_mitigtated_df$MIT_ERROR_RDscaleS <- mitigtated_tow_list_df$MIT_ERROR_RDscaleS
cat_mitigtated_df$MIT_ERROR_MCDS <- mitigtated_tow_list_df$MIT_ERROR_MCDS



# apply the appropriate error according to abundance time series
# RDtrendS is handled separately downstream using tow_list.csv

cat_mitigtated_df <- cat_mitigtated_df %>%
  mutate(
    NPERTOW = case_when(
      MIT_ADDED == 1 & REGION == "RDscaleS" ~ pmax(0, NPERTOW + MIT_ERROR_RDscaleS),
      MIT_ADDED == 1 & REGION == "MCDS" ~ pmax(0, NPERTOW + MIT_ERROR_MCDS),
      TRUE ~ NPERTOW
    )
  )


## 3.3 save the tow information ----

write.csv(mitigtated_tow_list_df, "results/indices for assessment/surfclam/mitigation_1_WEA_distance/tow_list.csv", row.names = FALSE)

cat_mitigtated_df <- cat_mitigtated_df %>%
  select(-MIT_ADDED, -MIT_ERROR_RDscaleS, -MIT_ERROR_MCDS)

write.csv(cat_mitigtated_df,"results/indices for assessment/surfclam/mitigation_1_WEA_distance/catch_by_tow.csv", row.names = FALSE)

remove(cat_df, tow_list_df)


## 3.4 generate abundance indices ----


# !!!! go to script 1.7.4 for calculating MIT.1 indices

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
  
  temp_tow_list <- subset(tow_list_df, YEAR == i)
  
  # determine lost tows
  N_TOW_lost <- sum(temp_tow_list$NEW_RANK == 0)
  
  # retain accessible tows
  temp_tow_retained <- temp_tow_list %>%
    filter(NEW_RANK != 0) %>%
    mutate(MIT_ADDED = 0)
  
  # pick the same amount of tows based on distance to WEA (closest)
  temp_tow_added <- temp_tow_list %>%
    filter(NEW_RANK != 0) %>%
    arrange(NEW_RANK) %>%
    slice_head(n = N_TOW_lost) %>%
    mutate(MIT_ADDED = 1)
  
  # combine retained and added tows
  temp_tow_mitigtated_list <- bind_rows(temp_tow_retained, temp_tow_added)
  
  
  # combine to the full database
  mitigtated_tow_list_df <- bind_rows(mitigtated_tow_list_df, temp_tow_mitigtated_list)
  
  remove(temp_tow_list, N_TOW_lost, temp_tow_retained, temp_tow_added, temp_tow_mitigtated_list)
}

remove(i)



## 4.2 add stratum-level observation error ----

variance_RDtrendS_df <- read.csv("results/indices for assessment/quahog/stratum_year.RDtrendS.csv") %>%
  select(YEAR, STRATUM, VAR_STRATUM) %>%
  distinct() %>%
  rename(VAR_RDtrendS = VAR_STRATUM)

variance_RDscaleS_df <- read.csv("results/indices for assessment/quahog/stratum_year.RDscaleS.csv") %>%
  select(YEAR, STRATUM, VAR_STRATUM) %>%
  distinct() %>%
  rename(VAR_RDscaleS = VAR_STRATUM)

variance_MCDS_df <- read.csv("results/indices for assessment/quahog/stratum_year.MCDS.csv") %>%
  select(YEAR, STRATUM, VAR_STRATUM) %>%
  distinct() %>%
  rename(VAR_MCDS = VAR_STRATUM)


# attach the three variance estimates to the mitigation tow list

mitigtated_tow_list_df <- mitigtated_tow_list_df %>%
  left_join(variance_RDtrendS_df, by = c("YEAR", "STRATUM")) %>%
  left_join(variance_RDscaleS_df, by = c("YEAR", "STRATUM")) %>%
  left_join(variance_MCDS_df,by = c("YEAR", "STRATUM")) %>%
  mutate(
    # assign variance = 0 where a year x stratum variance is unavailable
    VAR_RDtrendS = ifelse(is.na(VAR_RDtrendS), 0, VAR_RDtrendS),
    VAR_RDscaleS = ifelse(is.na(VAR_RDscaleS), 0, VAR_RDscaleS),
    VAR_MCDS     = ifelse(is.na(VAR_MCDS),     0, VAR_MCDS),
    # initialize mitigation errors
    MIT_ERROR_RDtrendS = 0,
    MIT_ERROR_RDscaleS = 0,
    MIT_ERROR_MCDS     = 0
  )


# fixed seed for reproducibility
set.seed(1234)

# identify added pseudo-observations
added_rows <- which(mitigtated_tow_list_df$MIT_ADDED == 1)


# generate separate observation errors for each abundance time series

mitigtated_tow_list_df$MIT_ERROR_RDtrendS[added_rows] <- rnorm(
  length(added_rows),
  mean = 0,
  sd = sqrt(mitigtated_tow_list_df$VAR_RDtrendS[added_rows])
)

mitigtated_tow_list_df$MIT_ERROR_RDscaleS[added_rows] <- rnorm(
  length(added_rows),
  mean = 0,
  sd = sqrt(mitigtated_tow_list_df$VAR_RDscaleS[added_rows])
)

mitigtated_tow_list_df$MIT_ERROR_MCDS[added_rows] <- rnorm(
  length(added_rows),
  mean = 0,
  sd = sqrt(mitigtated_tow_list_df$VAR_MCDS[added_rows])
)


# generate catch database for RDscaleS and MCDS

cat_df <- read.csv("results/stratified.mean.indices/quahog/total.catch.by.tow.csv")


# generate catch database based on mitigation strategy
cat_mitigtated_df <- cat_df[match(mitigtated_tow_list_df$ID, cat_df$ID), ]

rownames(cat_mitigtated_df) <- 1:nrow(cat_mitigtated_df)


# carry mitigation information into catch database
cat_mitigtated_df$MIT_ADDED <- mitigtated_tow_list_df$MIT_ADDED
cat_mitigtated_df$MIT_ERROR_RDscaleS <- mitigtated_tow_list_df$MIT_ERROR_RDscaleS
cat_mitigtated_df$MIT_ERROR_MCDS <- mitigtated_tow_list_df$MIT_ERROR_MCDS



# apply the appropriate error according to abundance time series
# RDtrendS is handled separately downstream using tow_list.csv

cat_mitigtated_df <- cat_mitigtated_df %>%
  mutate(
    NPERTOW = case_when(
      MIT_ADDED == 1 & REGION == "RDscaleS" ~ pmax(0, NPERTOW + MIT_ERROR_RDscaleS),
      MIT_ADDED == 1 & REGION == "MCDS" ~ pmax(0, NPERTOW + MIT_ERROR_MCDS),
      TRUE ~ NPERTOW
    )
  )


## 4.3 save the tow information ----

write.csv(mitigtated_tow_list_df, "results/indices for assessment/quahog/mitigation_1_WEA_distance/tow_list.csv", row.names = FALSE)

cat_mitigtated_df <- cat_mitigtated_df %>%
  select(-MIT_ADDED, -MIT_ERROR_RDscaleS, -MIT_ERROR_MCDS)

write.csv(cat_mitigtated_df,"results/indices for assessment/quahog/mitigation_1_WEA_distance/catch_by_tow.csv", row.names = FALSE)

remove(cat_df, tow_list_df)


## 4.4 generate abundance indices ----


# !!!! go to script 1.7.5 for calculating MIT.1 indices

# ----------------------------------------------------- #







# 5. not used quahog ----

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






