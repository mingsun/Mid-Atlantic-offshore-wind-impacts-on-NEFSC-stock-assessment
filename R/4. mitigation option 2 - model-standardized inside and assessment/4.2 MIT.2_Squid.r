library(tidyverse)
library(randomForest)


# this mitigation scenario will use model to standardized tow level abundance to replace the original tow abundance
# note that all tows (inside and outside OWF) are standardized using full dataset
# according to the first paper, the best model is Random Forest for all scenario



# 1. generate full dataset with env.var ----

cat_df <- read.csv("results/indices for assessment/squid/catch_by_tow.csv")

# integrate env.var

station.spring.df <- read.csv("data/22561_NEFSCSpringFisheriesIndependentBottomTrawlData/22561_UNION_FSCS_SVSTA.csv") 
station.spring.df$SEASON <- "SPRING"

station.fall.df <- read.csv("data/22560_NEFSCFallFisheriesIndependentBottomTrawlData/22560_UNION_FSCS_SVSTA.csv") 
station.fall.df$SEASON <- "FALL"

station.df <- rbind(station.spring.df, station.fall.df)

station.df <- station.df %>%
  mutate(STRATUM = ifelse(nchar(STRATUM) == 5, str_sub(STRATUM, 2, 5), STRATUM),
         ID = paste(CRUISE6, STRATUM, TOW, sep = "."),
         BLOCK = paste(EST_YEAR, SEASON, STRATUM, sep = ".")) %>%
  select(ID, SURFTEMP, BOTTEMP, AVGDEPTH) %>%
  filter(ID %in% unique(cat_df$ID))

full_df <- cat_df %>%
  mutate(STRATUM = as.character(STRATUM)) %>%
  left_join(station.df, by = "ID") %>%
  rename(Year = YEAR, Lat = DECDEG_BEGLAT, Lon = DECDEG_BEGLON) %>% # the exact spelling is required by VAST 
  add_column(AreaSwept_km2 = 0.024, # 24000 m2 from doi:10.1093/icesjms/fsv166 and BTS protocol
             Vessel = NA, Pass = 0) %>% # required by the model, don't know why  
  filter(!if_any(c(Lat, Lon), is.na))

remove(station.spring.df, station.fall.df, station.df, cat_df)


full_RF_df <- full_df %>%   # convert YEAR into a factor
  mutate(Year = as.factor(Year)) %>%
  filter(!if_any(c(CATCH_WT_CAL), is.na)) %>%
  filter(Year != 2023)

# ---------------------------------------------------------------------------  #




# 2. training model build up ---- 

#     ## determine the number of trees to grow (n.tree) and number of variables runadomly samples at each split (mtry)
# set.seed(1) # for reproducibility
# train_index <- sample(1:nrow(full_df), 0.8 * nrow(full_df)) # 80% for training
# train_data <- full_df[train_index, ]
# valid_data <- full_df[-train_index, ]
# 
#   ### Perform cross-validation to determine optimal mtry
# rf_tune <- tuneRF(x = train_data[,c("YEAR", "LAT", "BOTTEMP", "AVGDEPTH")],
#                   y = train_data$BIOMASS,
#                   ntree = 1000, # change n.tree manually 1000, 1500
#                   mtryStart = 2,
#                   stepFactor = 1, improve = 0.01, trace = TRUE, plot = TRUE)
# 
# remove(train_index, train_data, valid_data, rf_tune)


# ---------------------------------------------------------------------------  #





# 3. fall model ----

fall_df <- subset(full_RF_df, SEASON == "FALL")


  ## 3.1 Kriging the missing BOTTEM by Year ----

library(gstat)

no_BOTTEMP_df <- data.frame() # an empty df to store all interpolated tows

for (y in unique(fall_df$Year)) {
  
  annual_df <- subset(fall_df, Year == y)
  
  if (sum(is.na(annual_df$BOTTEMP)) == 0) next
  
  # get the df to develop Kriging
  with_BOTTEMP_df <- subset(annual_df, !is.na(BOTTEMP))
  sp::coordinates(with_BOTTEMP_df) <-  ~Lon + Lat
  
  # develop a variogram model
  temp.vgm <- variogram(BOTTEMP ~ Lon + Lat, data = with_BOTTEMP_df, width = 0.01)
  # temp.fit <- fit.variogram(temp.vgm, vgm(c( "Nug")))
  temp.fit <- fit.variogram(temp.vgm, vgm(c("Nug", "Sph","Mat", "Pen" , "Ste", "Cir" , "Bes", "Exc")))
  
  # create the spatial scale to be interpolated
  missing_BT_coord <- subset(annual_df, is.na(BOTTEMP))[, c("Lat", "Lon")]
  sp::coordinates(missing_BT_coord) <-  ~Lon + Lat
  
  # perform kriging interpolation for the tows without BOTTEMP
  temp.kriged <- krige(BOTTEMP ~ Lon + Lat, with_BOTTEMP_df, missing_BT_coord, model = temp.fit)
  
  # store the interpolated value
  missing_BT_df <- subset(annual_df, is.na(BOTTEMP))
  missing_BT_df$BOTTEMP <- temp.kriged$var1.pred
  
  no_BOTTEMP_df <- rbind(no_BOTTEMP_df, missing_BT_df)
}

remove(y, annual_df, with_BOTTEMP_df, temp.vgm, temp.fit, missing_BT_coord, temp.kriged, missing_BT_df)


# generate the final df for RF
fall_RF_df <- rbind(subset(fall_df, !is.na(BOTTEMP)),
                    no_BOTTEMP_df) %>%
              arrange(ID)

# last touch to fix the missing AVGDEPTH for ID = 201105.1050.1
fall_RF_df[is.na(fall_RF_df$AVGDEPTH),]$AVGDEPTH <- 37 # a value manually extracted from the same ID with that information


  ## 3.2 run Random Forest ----

set.seed(2) # ensure reproducibility
RF_fall <- randomForest(CATCH_WT_CAL ~ Year + Lat + BOTTEMP + AVGDEPTH, data = subset(fall_RF_df, OWF == "OUTSIDE"), mtry = 2, ntree = 1000)


  ## 3.3  predict the results and save ----
fall_RF_df <- fall_RF_df %>%
  mutate(PREDICTED_BIOMASS =  predict(RF_fall, newdata = ., type = "response")) %>%
  mutate(Final_CATCH_WT_CAL = case_when(OWF == "INSIDE" ~ PREDICTED_BIOMASS,
                                        OWF == "OUTSIDE" ~  CATCH_WT_CAL))
         
write.csv(fall_RF_df, "results/indices for assessment/squid/mitigation_2_model-standardized/Fall_RF_standardized_tow.csv", row.names = FALSE)


rm(list = setdiff(ls(), c("full_RF_df")))

# ---------------------------------------------------------------------------  #





# 4. spring model ----


spring_df <- subset(full_RF_df, SEASON == "SPRING")


## 3.1 Kriging the missing BOTTEM by Year ----

library(gstat)

no_BOTTEMP_df <- data.frame() # an empty df to store all interpolated tows

for (y in unique(spring_df$Year)) {
  
  annual_df <- subset(spring_df, Year == y)
  
  if (sum(is.na(annual_df$BOTTEMP)) == 0) next
  
  # get the df to develop Kriging
  with_BOTTEMP_df <- subset(annual_df, !is.na(BOTTEMP))
  sp::coordinates(with_BOTTEMP_df) <-  ~Lon + Lat
  
  # develop a variogram model
  temp.vgm <- variogram(BOTTEMP ~ Lon + Lat, data = with_BOTTEMP_df, width = 0.01)
  temp.fit <- fit.variogram(temp.vgm, vgm(c( "Nug")))
  # temp.fit <- fit.variogram(temp.vgm, vgm(c("Nug", "Sph","Mat", "Pen" , "Ste", "Cir" , "Bes", "Exc")))
  
  # create the spatial scale to be interpolated
  missing_BT_coord <- subset(annual_df, is.na(BOTTEMP))[, c("Lat", "Lon")]
  sp::coordinates(missing_BT_coord) <-  ~Lon + Lat
  
  # perform kriging interpolation for the tows without BOTTEMP
  temp.kriged <- krige(BOTTEMP ~ Lon + Lat, with_BOTTEMP_df, missing_BT_coord, model = temp.fit)
  
  # store the interpolated value
  missing_BT_df <- subset(annual_df, is.na(BOTTEMP))
  missing_BT_df$BOTTEMP <- temp.kriged$var1.pred
  
  no_BOTTEMP_df <- rbind(no_BOTTEMP_df, missing_BT_df)
}

remove(y, annual_df, with_BOTTEMP_df, temp.vgm, temp.fit, missing_BT_coord, temp.kriged, missing_BT_df)


# generate the final df for RF
spring_RF_df <- rbind(subset(spring_df, !is.na(BOTTEMP)),
                    no_BOTTEMP_df) %>%
  arrange(ID)

# last touch to fix the missing AVGDEPTH by ID 
spring_RF_df[is.na(spring_RF_df$AVGDEPTH),]$ID

spring_RF_df[spring_RF_df$ID == "200803.3260.2" ,]$AVGDEPTH <- 22 # same as 200803.3260.1
spring_RF_df[spring_RF_df$ID == "201302.3610.4" ,]$AVGDEPTH <- 47 # same as 201302.3610.2
spring_RF_df[spring_RF_df$ID == "201702.1020.2" ,]$AVGDEPTH <- 58 # same as the other 201702.1020.2
spring_RF_df[spring_RF_df$ID == "201702.3110.3" ,]$AVGDEPTH <- 23 # same as the other 201702.3110.3
spring_RF_df[spring_RF_df$ID == "201702.3380.1" ,]$AVGDEPTH <- mean(c(31,23)) # the average of strata 3350 and 3410
spring_RF_df[spring_RF_df$ID == "201802.1640.3" ,]$AVGDEPTH <- 183 # same as 201802.1640.2




## 3.2 run Random Forest ----

set.seed(2) # ensure reproducibility
RF_spring <- randomForest(CATCH_WT_CAL ~ Year + Lat + BOTTEMP + AVGDEPTH, data = subset(spring_RF_df, OWF == "OUTSIDE"), mtry = 2, ntree = 1000)


## 3.3  predict the results and save ----
spring_RF_df <- spring_RF_df %>%
  mutate(PREDICTED_BIOMASS =  predict(RF_spring, newdata = ., type = "response")) %>%
  mutate(Final_CATCH_WT_CAL = case_when(OWF == "INSIDE" ~ PREDICTED_BIOMASS,
                                        OWF == "OUTSIDE" ~  CATCH_WT_CAL))

write.csv(spring_RF_df, "results/indices for assessment/squid/mitigation_2_model-standardized/spring_RF_standardized_tow.csv", row.names = FALSE)




# -------------------------------------------------------------------------------------------------------- #

rm(list = ls())



# 5. Design-based abundance indices ----

fall_RF_df <- read.csv("results/indices for assessment/squid/mitigation_2_model-standardized/fall_RF_standardized_tow.csv")
spring_RF_df <- read.csv("results/indices for assessment/squid/mitigation_2_model-standardized/spring_RF_standardized_tow.csv")

RF_df <- rbind(fall_RF_df, spring_RF_df) %>%
  select(Year, SEASON, STRATUM, TOW, ID, TOTAL_N_STATION, STRATUM_AREA, Final_CATCH_WT_CAL)



  ## 5.1 mean numbers by stratum ---------------------------------------------------------------------------------------------

mean_N_stratum_df <- RF_df %>%
  group_by(Year, SEASON, STRATUM) %>%
  mutate(MEAN_N_STRATUM = sum(Final_CATCH_WT_CAL, na.rm = TRUE)/TOTAL_N_STATION) %>% # mean within a strata
  mutate(VAR_STRATUM = sum((Final_CATCH_WT_CAL - MEAN_N_STRATUM)^2, na.rm = TRUE)/ (TOTAL_N_STATION - 1)) %>% # variance by stratum
  # mutate(VAR_STRATUM = var(EXPCATCHNUM)) %>% # variance by stratum, same as last line
  ungroup() 


  ## 5.2 stratified mean numbers ---------------------------------------------------------------------------------------------

stratified_mean_N_df <- mean_N_stratum_df %>%
  select(c(Year, SEASON, STRATUM, TOTAL_N_STATION, STRATUM_AREA, MEAN_N_STRATUM, VAR_STRATUM)) %>%
  # select(c(YEAR, SEASON, STRATUM, N.STATION, REL_WEIGHT, MEAN_N_STRATUM)) %>%
  distinct() %>% # downsize the data frame to a minimal without repetitive rows
  group_by(Year, SEASON) %>%
  mutate(REL_WEIGHT = STRATUM_AREA/sum(STRATUM_AREA, na.rm = TRUE)) %>%
  summarize(STRATIFIED_MEAN_N = weighted.mean(MEAN_N_STRATUM, w = REL_WEIGHT), # stratified mean
            STRATIFIED_VAR = sum(REL_WEIGHT^2 * VAR_STRATUM / TOTAL_N_STATION, na.rm = TRUE), # variance
            STRATIFIED_SE = sqrt(STRATIFIED_VAR)) %>%   # standard deviance
  mutate(up_CI_95 = qnorm(0.975, mean = STRATIFIED_MEAN_N, sd = STRATIFIED_SE),
         lo_CI_95 = qnorm(0.025, mean = STRATIFIED_MEAN_N, sd = STRATIFIED_SE),
         lo_CI_95 = ifelse(lo_CI_95 < 0, yes = 0, no = lo_CI_95),
         CV = STRATIFIED_SE / STRATIFIED_MEAN_N) %>%
  ungroup()


write.csv(stratified_mean_N_df, "results/indices for assessment/squid/mitigation_2_model-standardized/MIT.2.stratified.mean.indices.csv", row.names = FALSE)

# Then paste the calculated abundance indices to the worksheet in
# E:\Stony Brook job\NYSERDA Offshore wind\assessment model\squid\MIT.2 model

# -------------------------------------------------------------------------------------------------------- #

