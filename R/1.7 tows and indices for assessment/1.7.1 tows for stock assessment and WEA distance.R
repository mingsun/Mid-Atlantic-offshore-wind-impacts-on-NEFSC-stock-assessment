library(tidyverse)

# A. BTS ----

  ## station list ----
station.spring.df <- read.csv("data/22561_NEFSCSpringFisheriesIndependentBottomTrawlData/22561_UNION_FSCS_SVSTA.csv") 
station.fall.df <- read.csv("data/22560_NEFSCFallFisheriesIndependentBottomTrawlData/22560_UNION_FSCS_SVSTA.csv") %>%
  mutate(STRATUM = substr(STRATUM, 2, 5))

station.spring.df$SEASON <- "SPRING"
station.fall.df$SEASON <- "FALL"

BTS_tows_list_df <- rbind(station.spring.df, station.fall.df)  %>%
  select(CRUISE6, STRATUM, TOW, STATION, AREA, SVVESSEL, SVGEAR, EST_YEAR, DECDEG_BEGLAT, DECDEG_BEGLON, SEASON, SHG, TOGA) %>%
  mutate(ID = paste(CRUISE6, STRATUM, TOW, sep = ".")) %>%
  rename(YEAR = EST_YEAR)
  
remove(station.spring.df, station.fall.df)


  ## strata area ----

stra_area_df <- read.csv("data/22560_NEFSCFallFisheriesIndependentBottomTrawlData/SVDBS_SupportTables/SVDBS_SVMSTRATA.csv") %>%
  select(stratum, stratum_area) %>%
  mutate(stratum = as.character(stratum)) %>%
  rename(STRATUM = stratum, STRATUM_AREA = stratum_area) %>%
  distinct(.keep_all = TRUE)


 ## add inside/outside OWF feature
BTS_overlay_df <- read.csv("results/BTS_tows_overlay_final.csv")

BTS_tows_list_df <- BTS_tows_list_df %>%
  mutate(OWF = ifelse(ID %in% unique(BTS_overlay_df$ID), "INSIDE", "OUTSIDE"))  %>%
  left_join(stra_area_df)   %>%
  left_join(stra_area_df) # add strata area

# ----------------------------------------------------------------------------- #



  ## 1. summer flounder ----

    # create the list of tow used by season
fall_SF_stra_list <- c(01010, 01050, 01090, 01610, 01650, 01690, 01730, 03010, 03020, 03030, 03040, 03050, 03060, 03070, 
                        03080, 03090, 03100, 03110, 03120, 03130, 03140, 03150, 03160, 03170, 03180, 03190, 03200, 03210, 
                        03220, 03230, 03240, 03250, 03260, 03270, 03280, 03290, 03300, 03310, 03320, 03330, 03340, 03350, 
                        03360, 03370, 03380, 03390, 03400, 03410, 03420, 03430, 03440, 03450, 03460, 03470, 03480, 03490, 
                        03500, 03510, 03520, 03530, 03540, 03550, 03560, 03570, 03580, 03590, 03600, 03610) # 68
fall_SF_stra_list <- substr(fall_SF_stra_list, 1,4) # remove the first 0 as strata in other excel files can be broken 

spring_SF_stra_list <- c(01010, 01020, 01030, 01040, 01050, 01060, 01070, 01080, 01090, 01100, 01110, 01120, 01610, 01620, 
                          01630, 01640, 01650, 01660, 01670, 01680, 01690, 01700, 01710, 01720, 01730, 01740, 01750, 01760) # 28
spring_SF_stra_list <- substr(spring_SF_stra_list, 1,4) # remove the first 0 as strata in other excel files can be broken 

              # the list is acquired from SAW 66 table in appendix

  # generate tow list by season

SF_fall_tow_df <- BTS_tows_list_df  %>%
  filter(YEAR >= 1982 & SEASON == "FALL") %>% # summer flounder assessment uses survey data since 1982
  filter(STRATUM %in% fall_SF_stra_list)  

SF_spring_tow_df <- BTS_tows_list_df  %>%
  filter(YEAR >= 1982 & SEASON == "SPRING") %>% # summer flounder assessment uses survey data since 1982
  filter(STRATUM %in% spring_SF_stra_list) 


  # combine 
SF_tow_list <- rbind(SF_fall_tow_df, SF_spring_tow_df) %>%
  distinct(ID, .keep_all = TRUE) %>%
  filter((SHG <= 136 & YEAR <= 2008) |
         (TOGA <= 1330  & YEAR >=2009)) # tow evaluation criteria



write.csv(SF_tow_list, "results/stock assessment/summer flounder/tow data/full.tow.list.for.assessment.csv", row.names = FALSE)

remove(fall_SF_stra_list, spring_SF_stra_list, SF_fall_tow_df, SF_spring_tow_df, SF_tow_list, BTS_overlay_df)


# ----------------------------------------------------------------------------- #




  ## 2. squid ----

  # the list of tows are directly provided by NOAA Jessica


# 8617
LS_tow_list <- read.csv("results/stratified.mean.indices/longfin.squid/calibrated.total.biomass.by.tow.csv") %>%
  select(YEAR, SEASON, STRATUM, ID, SHG, TOGA, CATCH_WT_CAL) %>%
  filter(YEAR >= 1976) %>% # longfin squid assessment uses survey data since 1976
  mutate(STRATUM = as.character(STRATUM)) %>%
  distinct() %>%
  left_join(BTS_tows_list_df)

write.csv(LS_tow_list, "results/stock assessment/squid/tow data/full.tow.list.for.assessment.csv", row.names = FALSE)


LS_tow_list[duplicated(LS_tow_list$ID) | duplicated(LS_tow_list$ID, fromLast = TRUE), ]

# ----------------------------------------------------------------------------- #



# 3. surfclam ----

full_tow_df <- read.csv("results/stratified.mean.indices/surfclam/total.catch.by.tow.csv") # biomass/abundance both included

AS_QQ_overlay_df <- read.csv("results/AS_OQ_tows_overlay_final.csv") %>%
  mutate(ID.temp = paste(CRUISE6, STATION, sep = "."))

full_tow_df <- full_tow_df %>% 
  mutate(ID.temp = paste(CRUISE6, STATION, sep = ".")) %>% # add this because the WEA overlay strata are old BTS strata
  mutate(OWF = ifelse(ID.temp %in% unique(AS_QQ_overlay_df$ID.temp), "INSIDE", "OUTSIDE")) %>% 
  select(-c(ID.temp)) %>%
  filter(REGION == "SVAtoSNE", !is.na(LON)) 

write.csv(full_tow_df, "results/stock assessment/surfclam/tow data/full.info.tow.list.for.assessment.csv", row.names = FALSE)

remove(full_tow_df, AS_QQ_overlay_df)




#  4. quahog ----

full_tow_df <- read.csv("results/stratified.mean.indices/quahog/total.catch.by.tow.csv") # biomass/abundance both included

AS_QQ_overlay_df <- read.csv("results/AS_OQ_tows_overlay_final.csv") %>%
  mutate(ID.temp = paste(CRUISE6, STATION, sep = "."))

full_tow_df <- full_tow_df %>% 
  mutate(ID.temp = paste(CRUISE6, STATION, sep = ".")) %>% # add this because the WEA overlay strata are old BTS strata
  mutate(OWF = ifelse(ID.temp %in% unique(AS_QQ_overlay_df$ID.temp), "INSIDE", "OUTSIDE")) %>% 
  select(-c(ID.temp)) %>%
  filter(REGION == "SVAtoSNE", !is.na(LON)) 

write.csv(full_tow_df, "results/stock assessment/quahog/tow data/full.info.tow.list.for.assessment.csv", row.names = FALSE)


