library(tidyverse)

# RDtrendS: total N stratifed mean, 1982-2011, assuming numbers per m2 consistent
# RDscaleS: total N straitifed mean times the area size, 1997-2011
# MCDS: similar to RDscaleS but based on sensor tow distance (DISTANCEDETAIL = SENSOR), 2012-2022



stra_area_df <- read.csv("results/stratified.mean.indices/surfclam/strata.area.csv") %>% # strata area as mean
  filter(REGION == "SVAtoSNE") %>%
  mutate(WEIGHT = AREASQNM/sum(AREASQNM))

AS_QQ_overlay_df <- read.csv("results/AS_OQ_tows_overlay_final.csv") %>%
  mutate(ID.temp = paste(CRUISE6, STATION, sep = "."))



# 1. RDtrends ----

tow_RDtrendS_df <- read.csv("data/NOAA.stock.data/surfclam/TOW_DATA.CSV") %>% 
  select(REGNAM, CRUISE6, YR, STATION, STRATUM, TOW, TOTALN, ASWEPT_M2, DISTANCEDETAIL) %>%
  rename(YEAR = YR) %>%
  filter(YEAR <= 2011, REGNAM == "SVAtoSNE") %>%
  filter(nchar(TOW) != 9) %>% # remove the repeated tows, looks like the data are already filled in to the original record
  mutate(ID = paste(CRUISE6, STRATUM, STATION, sep = ".")) %>%
  distinct() %>% 
  group_by(YEAR, STRATUM) %>%
  mutate(TOTAL_N_STATION = length(STATION)) %>%
  ungroup()


  ## 1.1 full dataset ----
RDtrendS_full_df <- tow_RDtrendS_df %>% 
  mutate(N_per_M2 = TOTALN) %>% # average number per m2
  group_by(YEAR, STRATUM) %>%
  summarise(mean_N_per_M2 = mean(N_per_M2, na.rm = TRUE), # average number per m2 by strata
            VAR_STRATUM = var(N_per_M2, na.rm = TRUE),
            SD_STRATUM = sqrt(VAR_STRATUM),
            TOTAL_N_STATION = n(),
            .groups = "drop") %>%
  left_join(stra_area_df) %>% # add the area weighting
  group_by(YEAR) %>%
  summarise(VALUE = weighted.mean(mean_N_per_M2, w = WEIGHT), # stratified mean
            VAR = sum((WEIGHT^2) * (VAR_STRATUM / TOTAL_N_STATION)), # stratified variance
            SE = sqrt(VAR),
            CV = SE/VALUE,
            log_SE = sqrt(log(1 + CV^2))) 

write.csv(RDtrendS_full_df, "results/indices for assessment/surfclam/original.RDtrendS.csv", row.names = FALSE)



  ## 1.2 WEE dataset ----

RDtrendS_WEE_df <- tow_RDtrendS_df %>% 
  mutate(ID.temp = paste(CRUISE6, STATION, sep = ".")) %>% # using a new id here because the overlay ID have different stratum coding system
  filter(!ID.temp %in% unique(AS_QQ_overlay_df$ID.temp)) %>% 
  mutate(N_per_M2 = TOTALN) %>% # average number per m2
  group_by(YEAR, STRATUM) %>% 
  summarise(mean_N_per_M2 = mean(N_per_M2, na.rm = TRUE), # average number per m2 by strata
            VAR_STRATUM = var(N_per_M2, na.rm = TRUE),
            SD_STRATUM = sqrt(VAR_STRATUM),
            TOTAL_N_STATION = n(),
            .groups = "drop") %>%
  left_join(stra_area_df) %>% # add the area weighting
  group_by(YEAR) %>%
  summarise(VALUE = weighted.mean(mean_N_per_M2, w = WEIGHT), # stratified mean
            VAR = sum((WEIGHT^2) * (VAR_STRATUM / TOTAL_N_STATION)), # stratified variance
            SE = sqrt(VAR),
            CV = SE/VALUE,
            log_SE = sqrt(log(1 + CV^2))) 

write.csv(RDtrendS_WEE_df, "results/indices for assessment/surfclam/WEE/WEE.RDtrendS.csv", row.names = FALSE)


  ## extract WEE ratio relative to original
RDtrendS_WEE_VALUE_ratio <- RDtrendS_WEE_df$VALUE/RDtrendS_full_df$VALUE
RDtrendS_WEE_STDERR_ratio <- RDtrendS_WEE_df$log_SE/RDtrendS_full_df$log_SE

save(RDtrendS_WEE_VALUE_ratio, RDtrendS_WEE_STDERR_ratio, file = "results/indices for assessment/surfclam/WEE/WEE.RDtrendS_ratio.Rdata")


rm(list = setdiff(ls(), c("stra_area_df", "AS_QQ_overlay_df", "RDtrendS_full_df", "tow_RDtrendS_df")))



  ## 1.3 MIT.1 dataset ----

mitigtated_tow_list_df <- read.csv("results/indices for assessment/surfclam/mitigation_1_WEA_distance/tow_list.csv") %>%
  filter(YEAR <= 2011) 

RDtrendS_MIT.1_df <- tow_RDtrendS_df[match(mitigtated_tow_list_df$ID, tow_RDtrendS_df$ID), ] # upsample the original tow list

RDtrendS_MIT.1_df <- RDtrendS_MIT.1_df %>% 
  mutate(N_per_M2 = TOTALN) %>% # average number per m2
  group_by(YEAR, STRATUM) %>% 
  summarise(mean_N_per_M2 = mean(N_per_M2, na.rm = TRUE), # average number per m2 by strata
            VAR_STRATUM = var(N_per_M2, na.rm = TRUE),
            SD_STRATUM = sqrt(VAR_STRATUM),
            TOTAL_N_STATION = n(),
            .groups = "drop") %>%
  left_join(stra_area_df) %>% # add the area weighting
  group_by(YEAR) %>%
  summarise(VALUE = weighted.mean(mean_N_per_M2, w = WEIGHT), # stratified mean
            VAR = sum((WEIGHT^2) * (VAR_STRATUM / TOTAL_N_STATION)), # stratified variance
            SE = sqrt(VAR),
            CV = SE/VALUE,
            log_SE = sqrt(log(1 + CV^2))) 


write.csv(RDtrendS_MIT.1_df, "results/indices for assessment/surfclam/mitigation_1_WEA_distance/MIT.1.RDtrendS.csv", row.names = FALSE)


## extract MIT.1 ratio relative to original
RDtrendS_MIT.1_VALUE_ratio <- RDtrendS_MIT.1_df$VALUE/RDtrendS_full_df$VALUE
RDtrendS_MIT.1_STDERR_ratio <- RDtrendS_MIT.1_df$log_SE/RDtrendS_full_df$log_SE

save(RDtrendS_MIT.1_VALUE_ratio, RDtrendS_MIT.1_STDERR_ratio, file = "results/indices for assessment/surfclam/mitigation_1_WEA_distance/MIT.1.RDtrendS_ratio.Rdata")



rm(list = setdiff(ls(), c("stra_area_df", "AS_QQ_overlay_df", "RDtrendS_full_df")))



  ## 1.4 MIT.2 dataset ----

RDtrendS_MIT.2_df <- read.csv("results/indices for assessment/surfclam/mitigation_2_model-standardized/Full_RF_standardized_tow.csv") %>%
  filter(YEAR <= 2011)  %>% # the list of tow for MIT.2
  mutate(TOTALN = Final_NPERTOW * ASWEPT_M2) # calculate TOTALN based on model standardized NPERTOW

RDtrendS_MIT.2_df <- RDtrendS_MIT.2_df %>% 
  mutate(N_per_M2 = TOTALN) %>% # average number per m2
  group_by(YEAR, STRATUM) %>% 
  summarise(mean_N_per_M2 = mean(N_per_M2, na.rm = TRUE), # average number per m2 by strata
            VAR_STRATUM = var(N_per_M2, na.rm = TRUE),
            SD_STRATUM = sqrt(VAR_STRATUM),
            TOTAL_N_STATION = n(),
            .groups = "drop") %>%
  left_join(stra_area_df) %>% # add the area weighting
  group_by(YEAR) %>%
  summarise(VALUE = weighted.mean(mean_N_per_M2, w = WEIGHT), # stratified mean
            VAR = sum((WEIGHT^2) * (VAR_STRATUM / TOTAL_N_STATION)), # stratified variance
            SE = sqrt(VAR),
            CV = SE/VALUE,
            log_SE = sqrt(log(1 + CV^2))) 


write.csv(RDtrendS_MIT.2_df, "results/indices for assessment/surfclam/mitigation_2_model-standardized/MIT.2.RDtrendS.csv", row.names = FALSE)


  ## extract MIT.2 ratio relative to original
RDtrendS_MIT.2_VALUE_ratio <- RDtrendS_MIT.2_df$VALUE/RDtrendS_full_df$VALUE
RDtrendS_MIT.2_STDERR_ratio <- RDtrendS_MIT.2_df$log_SE/RDtrendS_full_df$log_SE

save(RDtrendS_MIT.2_VALUE_ratio, RDtrendS_MIT.2_STDERR_ratio, file = "results/indices for assessment/surfclam/mitigation_2_model-standardized/MIT.2.RDtrendS_ratio.Rdata")



rm(list = setdiff(ls(), c("stra_area_df", "AS_QQ_overlay_df", "RDtrendS_full_df")))




  ## 1.5 MIT.3 dataset ----

RDtrendS_MIT.3_df <- read.csv("results/indices for assessment/surfclam/mitigation_3_model-based indices/VAST_RDtrendS/final model/VAST_RDtrendS_Indices.csv") %>%
  rename(VALUE = fit) %>%
  mutate(YEAR = c(1982, 1983, 1984, 1986, 1989, 1992, 1994, 1997, 1999, 2002, 2005, 2008, 2011)) %>%
  mutate(SE = (upr - lwr)/(2 * 1.96),
         VAR = SE^2,
         CV = SE/VALUE,
         log_SE = sqrt(log(1 + CV^2))) %>%
  select(YEAR, VALUE, VAR, SE, CV, log_SE)
  

write.csv(RDtrendS_MIT.3_df, "results/indices for assessment/surfclam/mitigation_3_model-based indices/MIT.3.RDtrendS.csv", row.names = FALSE)


## extract MIT.3 ratio relative to original
RDtrendS_MIT.3_VALUE_ratio <- RDtrendS_MIT.3_df$VALUE/RDtrendS_full_df$VALUE
RDtrendS_MIT.3_STDERR_ratio <- RDtrendS_MIT.3_df$log_SE/RDtrendS_full_df$log_SE

save(RDtrendS_MIT.3_VALUE_ratio, RDtrendS_MIT.3_STDERR_ratio, file = "results/indices for assessment/surfclam/mitigation_3_model-based indices/MIT.3.RDtrendS_ratio.Rdata")



rm(list = setdiff(ls(), c("stra_area_df", "AS_QQ_overlay_df", "RDtrendS_full_df")))



# ----------------------------------------------------------------------------------------------------- #





# 2. RDscaleS ----

catch_by_tow_df <- read.csv("results/stratified.mean.indices/surfclam/total.catch.by.tow.csv")%>%
  filter(YEAR <= 2011 & YEAR >= 1997, REGION == "SVAtoSNE")

  ## 2.1 full dataset ----

  # mean abundance by stratum
full_mean_N_stratum_df <- catch_by_tow_df %>%
  group_by(YEAR, REGION, STRATUM) %>%
  mutate(MEAN_N_STRATUM = mean(NPERTOW)) %>% # mean within a strata
  mutate(VAR_STRATUM = var(NPERTOW)) %>% # variance by stratum
  ungroup() %>% 
  select(c(YEAR, REGION, STRATUM, TOTAL_N_STATION, WEIGHT, MEAN_N_STRATUM, VAR_STRATUM)) %>%
  distinct() %>%
  arrange(YEAR, REGION)


  # stratified mean
full_stratified_mean_N_df <- full_mean_N_stratum_df %>%
  group_by(YEAR, REGION) %>%
  summarize(STRATIFIED_MEAN_N = weighted.mean(MEAN_N_STRATUM, w = WEIGHT), # stratified mean
            VAR = sum(WEIGHT^2 * VAR_STRATUM / TOTAL_N_STATION, na.rm = TRUE), # variance
            SE = sqrt(VAR)) %>%   # standard deviance
  mutate(CV = SE / STRATIFIED_MEAN_N) %>%
  ungroup()


  # adjust to total abundance based on population area
RDscaleS_full_df <- catch_by_tow_df %>%
  select(YEAR, STRATUM, AREASQNM) %>%
  distinct()  %>%
  group_by(YEAR) %>%
  summarize(total.AREASQM = sum(AREASQNM * 3429904))  %>% # convert sqnm to square meter
  ungroup() %>%
  right_join(full_stratified_mean_N_df) %>%
  mutate(VALUE = round(total.AREASQM * STRATIFIED_MEAN_N/1000, -1),
         SE = round(total.AREASQM * SE/1000, -1),
         log_SE = sqrt(log(1 + CV^2))) %>%
  select(YEAR, VALUE, VAR, SE, CV, log_SE)


write.csv(RDscaleS_full_df, "results/indices for assessment/surfclam/original.RDscaleS.csv", row.names = FALSE)



rm(list = setdiff(ls(), c("stra_area_df", "AS_QQ_overlay_df", "catch_by_tow_df", "RDscaleS_full_df")))




  ## 2.2 WEE dataset ----

catch_by_tow_WEE_df <- catch_by_tow_df %>% 
  mutate(ID.temp = paste(CRUISE6, STATION, sep = ".")) %>% 
  filter(!ID.temp %in% unique(AS_QQ_overlay_df$ID.temp))

# mean abundance by stratum
WEE_mean_N_stratum_df <- catch_by_tow_WEE_df %>%
  group_by(YEAR, REGION, STRATUM) %>%
  mutate(MEAN_N_STRATUM = mean(NPERTOW)) %>% # mean within a strata
  mutate(VAR_STRATUM = var(NPERTOW)) %>% # variance by stratum
  ungroup() %>% 
  select(c(YEAR, REGION, STRATUM, TOTAL_N_STATION, WEIGHT, MEAN_N_STRATUM, VAR_STRATUM)) %>%
  distinct() %>%
  arrange(YEAR, REGION)


# stratified mean
WEE_stratified_mean_N_df <- WEE_mean_N_stratum_df %>%
  group_by(YEAR, REGION) %>%
  summarize(STRATIFIED_MEAN_N = weighted.mean(MEAN_N_STRATUM, w = WEIGHT), # stratified mean
            VAR = sum(WEIGHT^2 * VAR_STRATUM / TOTAL_N_STATION, na.rm = TRUE), # variance
            SE = sqrt(VAR)) %>%   # standard deviance
  mutate(CV = SE / STRATIFIED_MEAN_N) %>%
  ungroup()



# adjust to total abundance based on population area
RDscaleS_WEE_df <- catch_by_tow_df %>%
  select(YEAR, STRATUM, AREASQNM) %>%
  distinct()  %>%
  group_by(YEAR) %>%
  summarize(total.AREASQM = sum(AREASQNM * 3429904))  %>% # convert sqnm to square meter
  ungroup() %>%
  right_join(WEE_stratified_mean_N_df) %>%
  mutate(VALUE = round(total.AREASQM * STRATIFIED_MEAN_N/1000, -1),
         SE = round(total.AREASQM * SE/1000, -1),
         log_SE = sqrt(log(1 + CV^2))) %>%
  select(YEAR, VALUE, VAR, SE, CV, log_SE)


write.csv(RDscaleS_WEE_df, "results/indices for assessment/surfclam/WEE/WEE.RDscaleS.csv", row.names = FALSE)


  ## extract WEE ratio relative to original
RDscaleS_WEE_VALUE_ratio <- RDscaleS_WEE_df$VALUE/RDscaleS_full_df$VALUE
RDscaleS_WEE_STDERR_ratio <- RDscaleS_WEE_df$log_SE/RDscaleS_full_df$log_SE

save(RDscaleS_WEE_VALUE_ratio, RDscaleS_WEE_STDERR_ratio, file = "results/indices for assessment/surfclam/WEE/WEE.RDscaleS_ratio.Rdata")



rm(list = setdiff(ls(), c("stra_area_df", "AS_QQ_overlay_df", "catch_by_tow_df", "RDscaleS_full_df")))



  ## 2.3 MIT.1 dataset ----

mitigtated_tow_list_df <- read.csv("results/indices for assessment/surfclam/mitigation_1_WEA_distance/tow_list.csv") %>%
  filter(YEAR <= 2011 & YEAR >= 1997) 


catch_by_tow_MIT.1_df <- catch_by_tow_df[match(mitigtated_tow_list_df$ID, catch_by_tow_df$ID), ] # up sample the original tow list


# mean abundance by stratum
MIT.1_mean_N_stratum_df <- catch_by_tow_MIT.1_df %>%
  group_by(YEAR, REGION, STRATUM) %>%
  mutate(MEAN_N_STRATUM = mean(NPERTOW)) %>% # mean within a strata
  mutate(VAR_STRATUM = var(NPERTOW)) %>% # variance by stratum
  ungroup() %>% 
  select(c(YEAR, REGION, STRATUM, TOTAL_N_STATION, WEIGHT, MEAN_N_STRATUM, VAR_STRATUM)) %>%
  distinct() %>%
  arrange(YEAR, REGION)


# stratified mean
MIT.1_stratified_mean_N_df <- MIT.1_mean_N_stratum_df %>%
  group_by(YEAR, REGION) %>%
  summarize(STRATIFIED_MEAN_N = weighted.mean(MEAN_N_STRATUM, w = WEIGHT), # stratified mean
            VAR = sum(WEIGHT^2 * VAR_STRATUM / TOTAL_N_STATION, na.rm = TRUE), # variance
            SE = sqrt(VAR)) %>%   # standard deviance
  mutate(CV = SE / STRATIFIED_MEAN_N) %>%
  ungroup()


# adjust to total abundance based on population area
RDscaleS_MIT.1_df <- catch_by_tow_df %>%
  select(YEAR, STRATUM, AREASQNM) %>%
  distinct()  %>%
  group_by(YEAR) %>%
  summarize(total.AREASQM = sum(AREASQNM * 3429904))  %>% # convert sqnm to square meter
  ungroup() %>%
  right_join(MIT.1_stratified_mean_N_df) %>%
  mutate(VALUE = round(total.AREASQM * STRATIFIED_MEAN_N/1000, -1),
         log_SE = sqrt(log(1 + CV^2))) %>%
  select(YEAR, VALUE, VAR, SE, CV, log_SE)


write.csv(RDscaleS_MIT.1_df, "results/indices for assessment/surfclam/mitigation_1_WEA_distance/MIT.1.RDscaleS.csv", row.names = FALSE)


  ## extract MIT.1 ratio relative to original
RDscaleS_MIT.1_VALUE_ratio <- RDscaleS_MIT.1_df$VALUE/RDscaleS_full_df$VALUE
RDscaleS_MIT.1_STDERR_ratio <- RDscaleS_MIT.1_df$log_SE/RDscaleS_full_df$log_SE

save(RDscaleS_MIT.1_VALUE_ratio, RDscaleS_MIT.1_STDERR_ratio, file = "results/indices for assessment/surfclam/mitigation_1_WEA_distance/MIT.1.RDscaleS_ratio.Rdata")


rm(list = setdiff(ls(), c("stra_area_df", "AS_QQ_overlay_df", "catch_by_tow_df", "RDscaleS_full_df")))



  ## 2.4 MIT.2 dataset ----

catch_by_tow_MIT.2_df <- read.csv("results/indices for assessment/surfclam/mitigation_2_model-standardized/Full_RF_standardized_tow.csv") %>%
  filter(YEAR <= 2011 & YEAR >= 1997)  # the list of tow for MIT.2


# mean abundance by stratum
MIT.2_mean_N_stratum_df <- catch_by_tow_MIT.2_df %>%
  group_by(YEAR, STRATUM) %>%
  mutate(MEAN_N_STRATUM = mean(Final_NPERTOW)) %>% # mean within a strata
  mutate(VAR_STRATUM = var(Final_NPERTOW)) %>% # variance by stratum
  ungroup() %>% 
  select(c(YEAR, STRATUM, TOTAL_N_STATION, WEIGHT, MEAN_N_STRATUM, VAR_STRATUM)) %>%
  distinct() %>%
  arrange(YEAR)


# stratified mean
MIT.2_stratified_mean_N_df <- MIT.2_mean_N_stratum_df %>%
  group_by(YEAR) %>%
  summarize(STRATIFIED_MEAN_N = weighted.mean(MEAN_N_STRATUM, w = WEIGHT), # stratified mean
            VAR = sum(WEIGHT^2 * VAR_STRATUM / TOTAL_N_STATION, na.rm = TRUE), # variance
            SE = sqrt(VAR)) %>%   # standard deviance
  mutate(CV = SE / STRATIFIED_MEAN_N) %>%
  ungroup()


# adjust to total abundance based on population area
RDscaleS_MIT.2_df <- catch_by_tow_df %>%
  select(YEAR, STRATUM, AREASQNM) %>%
  distinct()  %>%
  group_by(YEAR) %>%
  summarize(total.AREASQM = sum(AREASQNM * 3429904))  %>% # convert sqnm to square meter
  ungroup() %>%
  right_join(MIT.2_stratified_mean_N_df) %>%
  mutate(VALUE = round(total.AREASQM * STRATIFIED_MEAN_N/1000, -1),
         log_SE = sqrt(log(1 + CV^2))) %>%
  select(YEAR, VALUE, VAR, SE, CV, log_SE)


write.csv(RDscaleS_MIT.2_df, "results/indices for assessment/surfclam/mitigation_2_model-standardized//MIT.2.RDscaleS.csv", row.names = FALSE)


## extract MIT.2 ratio relative to original
RDscaleS_MIT.2_VALUE_ratio <- RDscaleS_MIT.2_df$VALUE/RDscaleS_full_df$VALUE
RDscaleS_MIT.2_STDERR_ratio <- RDscaleS_MIT.2_df$log_SE/RDscaleS_full_df$log_SE

save(RDscaleS_MIT.2_VALUE_ratio, RDscaleS_MIT.2_STDERR_ratio, file = "results/indices for assessment/surfclam/mitigation_2_model-standardized/MIT.2.RDscaleS_ratio.Rdata")


rm(list = setdiff(ls(), c("stra_area_df", "AS_QQ_overlay_df", "catch_by_tow_df", "RDscaleS_full_df")))




  ## 2.5 MIT.3 dataset ----

catch_by_tow_MIT.3_df <- read.csv("results/stratified.mean.indices/surfclam/total.catch.by.tow.csv") %>%
  filter(YEAR <= 2011 & YEAR >= 1997)  # the list of tow for MIT.3


# stratified mean
MIT.3_stratified_mean_N_df <- read.csv("results/indices for assessment/surfclam/mitigation_3_model-based indices/VAST_RDscaleS/final model/VAST_RDscaleS_Indices.csv") %>%
  filter(YEAR <= 2011 & YEAR >= 1997) %>%
  rename(STRATIFIED_MEAN_N = fit) %>%
  mutate(SE = (upr - lwr)/(2 * 1.96),
       VAR = SE^2,
       CV = SE/STRATIFIED_MEAN_N,
       log_SE = sqrt(log(1 + CV^2))) %>%
  select(YEAR, STRATIFIED_MEAN_N, VAR, SE, CV, log_SE)


# adjust to total abundance based on population area
RDscaleS_MIT.3_df <- catch_by_tow_MIT.3_df %>%
  select(YEAR, STRATUM, AREASQNM) %>%
  distinct()  %>%
  group_by(YEAR) %>%
  summarize(total.AREASQM = sum(AREASQNM * 3429904))  %>% # convert sqnm to square meter
  ungroup() %>%
  right_join(MIT.3_stratified_mean_N_df) %>%
  mutate(VALUE = round(total.AREASQM * STRATIFIED_MEAN_N/1000, -1),
         log_SE = sqrt(log(1 + CV^2))) %>%
  select(YEAR, VALUE, VAR, SE, CV, log_SE)


write.csv(RDscaleS_MIT.3_df, "results/indices for assessment/surfclam/mitigation_3_model-based indices/MIT.3.RDscaleS.csv", row.names = FALSE)


## extract MIT.3 ratio relative to original
RDscaleS_MIT.3_VALUE_ratio <- RDscaleS_MIT.3_df$VALUE/RDscaleS_full_df$VALUE
RDscaleS_MIT.3_STDERR_ratio <- RDscaleS_MIT.3_df$log_SE/RDscaleS_full_df$log_SE

save(RDscaleS_MIT.3_VALUE_ratio, RDscaleS_MIT.3_STDERR_ratio, file = "results/indices for assessment/surfclam/mitigation_3_model-based indices/MIT.3.RDscaleS_ratio.Rdata")


rm(list = setdiff(ls(), c("stra_area_df", "AS_QQ_overlay_df", "catch_by_tow_df", "RDscaleS_full_df")))




# ----------------------------------------------------------------------------------------------------- #





# 3. MCDS ----


tow_MCDS_df <- read.csv("results/stratified.mean.indices/surfclam/total.catch.by.tow.csv")%>%
  filter(YEAR %in% c(2012, 2015, 2018, 2022), REGION == "SVAtoSNE")

  ## 3.1 full dataset ----

    # mean abundance by stratum
MCDS_mean_N_stratum_df <- tow_MCDS_df %>%
  group_by(YEAR, REGION, STRATUM) %>%
  mutate(MEAN_N_STRATUM = mean(NPERTOW)) %>% # mean within a strata
  mutate(VAR_STRATUM = var(NPERTOW)) %>% # variance by stratum
  ungroup() %>% 
  select(c(YEAR, REGION, STRATUM, TOTAL_N_STATION, WEIGHT, MEAN_N_STRATUM, VAR_STRATUM)) %>%
  distinct() %>%
  arrange(YEAR, REGION)


  # stratified mean
MCDS_stratified_mean_N_df <- MCDS_mean_N_stratum_df %>%
  group_by(YEAR, REGION) %>%
  summarize(STRATIFIED_MEAN_N = weighted.mean(MEAN_N_STRATUM, w = WEIGHT), # stratified mean
            VAR = sum(WEIGHT^2 * VAR_STRATUM / TOTAL_N_STATION, na.rm = TRUE), # variance
            SE = sqrt(VAR)) %>%   # standard deviance
  mutate(CV = SE / STRATIFIED_MEAN_N) %>%
  ungroup()



# adjust to total abundance based on population area
MCDS_full_df <- tow_MCDS_df %>%
  select(YEAR, STRATUM, AREASQNM) %>%
  distinct()  %>%
  group_by(YEAR) %>%
  summarize(total.AREASQM = sum(AREASQNM * 3429904))  %>% # convert sqnm to square meter
  ungroup() %>%
  right_join(MCDS_stratified_mean_N_df) %>%
  mutate(VALUE = round(total.AREASQM * STRATIFIED_MEAN_N/1000, -1),
         SE = round(total.AREASQM * SE/1000, -1),
         log_SE = sqrt(log(1 + CV^2))) %>%
  select(YEAR, VALUE, VAR, SE, CV, log_SE)


write.csv(MCDS_full_df, "results/indices for assessment/surfclam/original.MCDS.csv", row.names = FALSE)

rm(list = setdiff(ls(), c("stra_area_df", "AS_QQ_overlay_df", "tow_MCDS_df", "MCDS_full_df")))



  ## 3.2 WEE dataset ----

tow_MCDS_WEE_df <- tow_MCDS_df %>% 
  mutate(ID.temp = paste(CRUISE6, STATION, sep = ".")) %>% 
  filter(!ID.temp %in% unique(AS_QQ_overlay_df$ID.temp))

# mean abundance by stratum
MCDS_WEE_mean_N_stratum_df <- tow_MCDS_WEE_df %>%
  group_by(YEAR, REGION, STRATUM) %>%
  mutate(MEAN_N_STRATUM = mean(NPERTOW)) %>% # mean within a strata
  mutate(VAR_STRATUM = var(NPERTOW)) %>% # variance by stratum
  ungroup() %>% 
  select(c(YEAR, REGION, STRATUM, TOTAL_N_STATION, WEIGHT, MEAN_N_STRATUM, VAR_STRATUM)) %>%
  distinct() %>%
  arrange(YEAR, REGION)


# stratified mean
MCDS_WEE_stratified_mean_N_df <- MCDS_WEE_mean_N_stratum_df %>%
  group_by(YEAR, REGION) %>%
  summarize(STRATIFIED_MEAN_N = weighted.mean(MEAN_N_STRATUM, w = WEIGHT), # stratified mean
            VAR = sum(WEIGHT^2 * VAR_STRATUM / TOTAL_N_STATION, na.rm = TRUE), # variance
            SE = sqrt(VAR)) %>%   # standard deviance
  mutate(CV = SE / STRATIFIED_MEAN_N) %>%
  ungroup()


# adjust to total abundance based on population area
MCDS_WEE_df <- tow_MCDS_WEE_df %>%
  select(YEAR, STRATUM, AREASQNM) %>%
  distinct()  %>%
  group_by(YEAR) %>%
  summarize(total.AREASQM = sum(AREASQNM * 3429904))  %>% # convert sqnm to square meter
  ungroup() %>%
  right_join(MCDS_WEE_stratified_mean_N_df) %>%
  mutate(VALUE = round(total.AREASQM * STRATIFIED_MEAN_N/1000, -1),
         SE = round(total.AREASQM * SE/1000, -1),
         log_SE = sqrt(log(1 + CV^2))) %>%
  select(YEAR, VALUE, VAR, SE, CV, log_SE)

write.csv(MCDS_WEE_df, "results/indices for assessment/surfclam/WEE/WEE.MCDS.csv", row.names = FALSE)


  ## extract WEE ratio relative to original
MCDS_WEE_VALUE_ratio <- MCDS_WEE_df$VALUE/MCDS_full_df$VALUE
MCDS_WEE_STDERR_ratio <- MCDS_WEE_df$log_SE/MCDS_full_df$log_SE

save(MCDS_WEE_VALUE_ratio, MCDS_WEE_STDERR_ratio, file = "results/indices for assessment/surfclam/WEE/WEE.MCDS_ratio.Rdata")


rm(list = setdiff(ls(), c("stra_area_df", "AS_QQ_overlay_df", "tow_MCDS_df", "MCDS_full_df")))



  ## 3.3 MIT.1 dataset ----

mitigtated_tow_list_df <- read.csv("results/indices for assessment/surfclam/mitigation_1_WEA_distance/tow_list.csv") %>%
  filter(YEAR %in% c(2012, 2015, 2018, 2022))


tow_MCDS_MIT.1_df <- tow_MCDS_df[match(mitigtated_tow_list_df$ID, tow_MCDS_df$ID), ] # up sample the original tow list



# mean abundance by stratum
MCDS_MIT.1_mean_N_stratum_df <- tow_MCDS_MIT.1_df %>%
  group_by(YEAR, REGION, STRATUM) %>%
  mutate(MEAN_N_STRATUM = mean(NPERTOW)) %>% # mean within a strata
  mutate(VAR_STRATUM = var(NPERTOW)) %>% # variance by stratum
  ungroup() %>% 
  select(c(YEAR, REGION, STRATUM, TOTAL_N_STATION, WEIGHT, MEAN_N_STRATUM, VAR_STRATUM)) %>%
  distinct() %>%
  arrange(YEAR, REGION)


# stratified mean
MCDS_MIT.1_stratified_mean_N_df <- MCDS_MIT.1_mean_N_stratum_df %>%
  group_by(YEAR, REGION) %>%
  summarize(STRATIFIED_MEAN_N = weighted.mean(MEAN_N_STRATUM, w = WEIGHT), # stratified mean
            VAR = sum(WEIGHT^2 * VAR_STRATUM / TOTAL_N_STATION, na.rm = TRUE), # variance
            SE = sqrt(VAR)) %>%   # standard deviance
  mutate(CV = SE / STRATIFIED_MEAN_N) %>%
  ungroup()


# adjust to total abundance based on population area
MCDS_MIT.1_df <- tow_MCDS_MIT.1_df %>%
  select(YEAR, STRATUM, AREASQNM) %>%
  distinct()  %>%
  group_by(YEAR) %>%
  summarize(total.AREASQM = sum(AREASQNM * 3429904))  %>% # convert sqnm to square meter
  ungroup() %>%
  right_join(MCDS_MIT.1_stratified_mean_N_df) %>%
  mutate(VALUE = round(total.AREASQM * STRATIFIED_MEAN_N/1000, -1),
         log_SE = sqrt(log(1 + CV^2))) %>%
  select(YEAR, VALUE, VAR, SE, CV, log_SE)

write.csv(MCDS_MIT.1_df, "results/indices for assessment/surfclam/mitigation_1_WEA_distance/MIT.1.MCDS.csv", row.names = FALSE)


  ## extract MIT.1 ratio relative to original
MCDS_MIT.1_VALUE_ratio <- MCDS_MIT.1_df$VALUE/MCDS_full_df$VALUE
MCDS_MIT.1_STDERR_ratio <- MCDS_MIT.1_df$log_SE/MCDS_full_df$log_SE

save(MCDS_MIT.1_VALUE_ratio, MCDS_MIT.1_STDERR_ratio, file = "results/indices for assessment/surfclam/mitigation_1_WEA_distance/MIT.1.MCDS_ratio.Rdata")



rm(list = setdiff(ls(), c("stra_area_df", "AS_QQ_overlay_df", "tow_MCDS_df", "MCDS_full_df")))



  ## 3.4 MIT.2 dataset ----


tow_MCDS_MIT.2_df <- read.csv("results/indices for assessment/surfclam/mitigation_2_model-standardized/Full_RF_standardized_tow.csv") %>%
  filter(YEAR %in% c(2012, 2015, 2018, 2022)) # load the tow list for MIT.2



# mean abundance by stratum
MCDS_MIT.2_mean_N_stratum_df <- tow_MCDS_MIT.2_df %>%
  group_by(YEAR, STRATUM) %>%
  mutate(MEAN_N_STRATUM = mean(Final_NPERTOW)) %>% # mean within a strata
  mutate(VAR_STRATUM = var(Final_NPERTOW)) %>% # variance by stratum
  ungroup() %>% 
  select(c(YEAR, STRATUM, TOTAL_N_STATION, WEIGHT, MEAN_N_STRATUM, VAR_STRATUM)) %>%
  distinct() %>%
  arrange(YEAR)


# stratified mean
MCDS_MIT.2_stratified_mean_N_df <- MCDS_MIT.2_mean_N_stratum_df %>%
  group_by(YEAR, ) %>%
  summarize(STRATIFIED_MEAN_N = weighted.mean(MEAN_N_STRATUM, w = WEIGHT), # stratified mean
            VAR = sum(WEIGHT^2 * VAR_STRATUM / TOTAL_N_STATION, na.rm = TRUE), # variance
            SE = sqrt(VAR)) %>%   # standard deviance
  mutate(CV = SE / STRATIFIED_MEAN_N) %>%
  ungroup()


# adjust to total abundance based on population area
MCDS_MIT.2_df <- tow_MCDS_MIT.2_df %>%
  left_join(select(stra_area_df, -WEIGHT)) %>%
  select(YEAR, STRATUM, AREASQNM) %>%
  distinct()  %>%
  group_by(YEAR) %>%
  summarize(total.AREASQM = sum(AREASQNM * 3429904))  %>% # convert sqnm to square meter
  ungroup() %>%
  right_join(MCDS_MIT.2_stratified_mean_N_df) %>%
  mutate(VALUE = round(total.AREASQM * STRATIFIED_MEAN_N/1000, -1),
         log_SE = sqrt(log(1 + CV^2))) %>%
  select(YEAR, VALUE, VAR, SE, CV, log_SE)

write.csv(MCDS_MIT.2_df, "results/indices for assessment/surfclam/mitigation_2_model-standardized//MIT.2.MCDS.csv", row.names = FALSE)


## extract MIT.2 ratio relative to original
MCDS_MIT.2_VALUE_ratio <- MCDS_MIT.2_df$VALUE/MCDS_full_df$VALUE
MCDS_MIT.2_STDERR_ratio <- MCDS_MIT.2_df$log_SE/MCDS_full_df$log_SE

save(MCDS_MIT.2_VALUE_ratio, MCDS_MIT.2_STDERR_ratio, file = "results/indices for assessment/surfclam/mitigation_2_model-standardized/MIT.2.MCDS_ratio.Rdata")



rm(list = setdiff(ls(), c("stra_area_df", "AS_QQ_overlay_df", "tow_MCDS_df", "MCDS_full_df")))



  ## 3.5 MIT.3 dataset ----

catch_by_tow_MIT.3_df <- read.csv("results/stratified.mean.indices/surfclam/total.catch.by.tow.csv") %>%
  filter(YEAR %in% c(2012, 2015, 2018, 2022))  # the list of tow for MIT.3


# stratified mean
MCDS_MIT.3_stratified_mean_N_df <- read.csv("results/indices for assessment/surfclam/mitigation_3_model-based indices/VAST_MCDS/final model/VAST_MCDS_Indices.csv") %>%
  filter(YEAR %in% c(2012, 2015, 2018, 2022)) %>%
  rename(STRATIFIED_MEAN_N = fit) %>%
  mutate(SE = (upr - lwr)/(2 * 1.96),
         VAR = SE^2,
         CV = SE/STRATIFIED_MEAN_N,
         log_SE = sqrt(log(1 + CV^2))) %>%
  select(YEAR, STRATIFIED_MEAN_N, VAR, SE, CV, log_SE)


# adjust to total abundance based on population area
MCDS_MIT.3_df <- catch_by_tow_MIT.3_df %>%
  select(YEAR, STRATUM, AREASQNM) %>%
  distinct()  %>%
  group_by(YEAR) %>%
  summarize(total.AREASQM = sum(AREASQNM * 3429904))  %>% # convert sqnm to square meter
  ungroup() %>%
  right_join(MCDS_MIT.3_stratified_mean_N_df) %>%
  mutate(VALUE = round(total.AREASQM * STRATIFIED_MEAN_N/1000, -1),
         log_SE = sqrt(log(1 + CV^2))) %>%
  select(YEAR, VALUE, VAR, SE, CV, log_SE)


write.csv(MCDS_MIT.3_df, "results/indices for assessment/surfclam/mitigation_3_model-based indices/MIT.3.MCDS.csv", row.names = FALSE)


## extract MIT.3 ratio relative to original
MCDS_MIT.3_VALUE_ratio <- MCDS_MIT.3_df$VALUE/MCDS_full_df$VALUE
MCDS_MIT.3_STDERR_ratio <- MCDS_MIT.3_df$log_SE/MCDS_full_df$log_SE

save(MCDS_MIT.3_VALUE_ratio, MCDS_MIT.3_STDERR_ratio, file = "results/indices for assessment/surfclam/mitigation_3_model-based indices/MIT.3.MCDS_ratio.Rdata")


rm(list = setdiff(ls(), c("stra_area_df", "AS_QQ_overlay_df", "catch_by_tow_df", "MCDS_full_df")))





# ----------------------------------------------------------------------------------------------------- #


