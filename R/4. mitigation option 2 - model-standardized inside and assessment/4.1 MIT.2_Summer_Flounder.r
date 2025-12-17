library(tidyverse)
library(randomForest)
library(gstat)


# this mitigation scenario will use model to standardized tow level abundance to replace the original tow abundance
# note that all tows (inside and outside OWF) are standardized using full dataset
# according to the first paper, the best model is Random Forest for all scenario


# 1. ALB (1982 - 2008) ----


  ## 1.1 full dataset with env.var ----

cat_ALB_df <- read.csv("results/indices for assessment/summer flounder/catch_by_tow_ALB.csv") %>%
  filter(YEAR <= 2008)


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
  filter(ID %in% unique(cat_ALB_df$ID))


full_ALB_df <- cat_ALB_df %>%
  mutate(STRATUM = as.character(STRATUM)) %>%
  left_join(station.df, by = "ID") %>%
  rename(LAT = DECDEG_BEGLAT, LON = DECDEG_BEGLON)


  # manually check if location data are missing
sum(is.na(full_ALB_df$LON))

  # manually fix missing AVGDEPTH
sum(is.na(full_ALB_df$AVGDEPTH))
full_ALB_df[is.na(full_ALB_df$AVGDEPTH),]$AVGDEPTH <- 19 # same as ID=199804.3350.1



full_ALB_RF_df <- full_ALB_df %>%   # convert YEAR into a factor
  mutate(Year = as.factor(YEAR)) %>%
  filter(!if_any(c(NUMBER), is.na))

remove(station.spring.df, station.fall.df, station.df, cat_ALB_df, full_ALB_df)



## ----------------------------------- ##


  ## 1.2 data availability ----

    ## checking BOTTEMP data availability by year 

ggplot(full_ALB_RF_df) +
  geom_point(aes(x = LON, y = LAT)) +
  facet_wrap(.~YEAR)


temp_data_avail_df <- full_ALB_RF_df %>%
  summarize(N.TOW.MISSING = sum(is.na(BOTTEMP)), 
            N.TOTAL = length(BOTTEMP), 
            PROP.NO = N.TOW.MISSING/N.TOTAL,
            .by = c(YEAR, SEASON))


# save it and manually edit the INPUT.YEAR and then read back
# write.csv(temp_data_avail_df, "results/indices for assessment/summer flounder/mitigation_2_model-standardized/BOTTEMP.data.avail.csv", row.names = FALSE)
temp_data_avail_df <- read.csv("results/indices for assessment/summer flounder/mitigation_2_model-standardized/BOTTEMP.data.avail.csv")


  ## ----------------------------------- ##




  ## 1.3 Kriging the missing BOTTEM by Year ----


no_TEMP_df <- data.frame() # an empty df to store all interpolated tows

for (y in unique(full_ALB_RF_df$YEAR)) {
  
  for (s in unique(full_ALB_RF_df$SEASON)) {
    
    annual_df <- subset(full_ALB_RF_df, Year == y & SEASON == s)
    INPUT.YEAR <- subset(temp_data_avail_df, YEAR == y & SEASON == s)$INPUT.YEAR
    
    if (sum(is.na(annual_df$BOTTEMP)) == 0) next
    
    # get the df to develop Kriging
    with_TEMP_df <- subset(subset(full_ALB_RF_df, Year == INPUT.YEAR), !is.na(BOTTEMP))
    sp::coordinates(with_TEMP_df) <-  ~LON + LAT
    
    # develop a variogram model
    temp.vgm <- variogram(BOTTEMP ~ LON + LAT, data = with_TEMP_df, width = 0.01)
    # temp.fit <- fit.variogram(temp.vgm, vgm(c( "Nug")))
    temp.fit <- fit.variogram(temp.vgm, vgm(c("Nug", "Sph","Mat", "Pen" , "Ste", "Cir" , "Bes", "Exc")))
    
    
    # create the spatial scale to be interpolated
    missing_TEMP_coord <- subset(annual_df, is.na(BOTTEMP))[, c("LAT", "LON")]
    sp::coordinates(missing_TEMP_coord) <-  ~LON + LAT
    
    # perform kriging interpolation for the tows without BOTTEMP
    temp.kriged <- krige(BOTTEMP ~ LON + LAT, with_TEMP_df, missing_TEMP_coord, model = temp.fit)
    
    # store the interpolated value
    missing_TEMP_df <- subset(annual_df, is.na(BOTTEMP))
    missing_TEMP_df$BOTTEMP <- temp.kriged$var1.pred
    
    no_TEMP_df <- rbind(no_TEMP_df, missing_TEMP_df)
  
  }
}
  
remove(y, s, annual_df, INPUT.YEAR, with_TEMP_df, temp.vgm, temp.fit, missing_TEMP_coord, temp.kriged, missing_TEMP_df)

  
  # generate the final df for RF
full_ALB_RF_df <- rbind(subset(full_ALB_RF_df, !is.na(BOTTEMP)),
                    no_TEMP_df) %>%
              arrange(ID)

  ## ----------------------------------- ##



  ## 1.4 fit spring model ---- 

set.seed(2) # ensure reproducibility
RF_ALB_spring <- randomForest(NUMBER ~ YEAR + LAT + BOTTEMP + AVGDEPTH, data = subset(full_ALB_RF_df, SEASON == "SPRING" & OWF == "OUTSIDE"), mtry = 2, ntree = 1000)

  ##  predict the results and save
full_ALB_RF_spring_df <- full_ALB_RF_df %>%
  filter(SEASON == "SPRING") %>%
  mutate(PREDICTED_NUMBER =  predict(RF_ALB_spring, newdata = ., type = "response")) %>%
  mutate(Final_NUMBER = case_when(OWF == "INSIDE" ~ PREDICTED_NUMBER,
                                  OWF == "OUTSIDE" ~  NUMBER))

  ## ----------------------------------- ##


  ## 1.5 fit fall model ---- 

set.seed(2) # ensure reproducibility
RF_ALB_fall <- randomForest(NUMBER ~ YEAR + LAT + BOTTEMP + AVGDEPTH, data = subset(full_ALB_RF_df, SEASON == "FALL" & OWF == "OUTSIDE"), mtry = 2, ntree = 1000)

##  predict the results and save
full_ALB_RF_fall_df <- full_ALB_RF_df %>%
  filter(SEASON == "FALL") %>%
  mutate(PREDICTED_NUMBER =  predict(RF_ALB_fall, newdata = ., type = "response")) %>%
  mutate(Final_NUMBER = case_when(OWF == "INSIDE" ~ PREDICTED_NUMBER,
                                  OWF == "OUTSIDE" ~  NUMBER))

  ## ----------------------------------- ##

  
  ## 1.6 combine and save

full_ALB_RF_df <- rbind(full_ALB_RF_spring_df, full_ALB_RF_fall_df)

write.csv(full_ALB_RF_df, "results/indices for assessment/summer flounder/mitigation_2_model-standardized/Full_ALB_RF_standardized_tow.csv", row.names = FALSE)


rm(list = ls())


  ## ----------------------------------- ##




# 2. BIG (1982 - 2008) ----


  ## 2.1 full dataset with env.var ----

cat_BIG_df <- read.csv("results/indices for assessment/summer flounder/catch_by_tow_BIG.csv") 


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
  filter(ID %in% unique(cat_BIG_df$ID))


full_BIG_df <- cat_BIG_df %>%
  mutate(STRATUM = as.character(STRATUM)) %>%
  left_join(station.df, by = "ID") %>%
  rename(LAT = DECDEG_BEGLAT, LON = DECDEG_BEGLON)


# manually check if location data are missing
sum(is.na(full_BIG_df$LON))

# manually fix missing AVGDEPTH
sum(is.na(full_BIG_df$AVGDEPTH))
full_BIG_df[is.na(full_BIG_df$AVGDEPTH),]$AVGDEPTH <- 19 # same as ID=199804.3350.1



full_BIG_RF_df <- full_BIG_df %>%   # convert YEAR into a factor
  mutate(Year = as.factor(YEAR)) %>%
  filter(!if_any(c(NUMBER), is.na))

remove(station.spring.df, station.fall.df, station.df, cat_BIG_df, full_BIG_df)



## ----------------------------------- ##


  ## 2.2 data availability ----

## checking BOTTEMP data availability by year 

ggplot(full_BIG_RF_df) +
  geom_point(aes(x = LON, y = LAT)) +
  facet_wrap(.~YEAR)


temp_data_avail_df <- full_BIG_RF_df %>%
  summarize(N.TOW.MISSING = sum(is.na(BOTTEMP)), 
            N.TOTAL = length(BOTTEMP), 
            PROP.NO = N.TOW.MISSING/N.TOTAL,
            .by = c(YEAR, SEASON))


# save it and manually edit the INPUT.YEAR and then read back
# write.csv(temp_data_avail_df, "results/indices for assessment/summer flounder/mitigation_2_model-standardized/BOTTEMP.data.avail_BIG.csv", row.names = FALSE)
temp_data_avail_df <- read.csv("results/indices for assessment/summer flounder/mitigation_2_model-standardized/BOTTEMP.data.avail_BIG.csv")


## ----------------------------------- ##




  ## 2.3 Kriging the missing BOTTEM by Year ----


no_TEMP_df <- data.frame() # an empty df to store all interpolated tows

for (y in unique(full_BIG_RF_df$YEAR)) {
  
  for (s in unique(full_BIG_RF_df$SEASON)) {
    
    annual_df <- subset(full_BIG_RF_df, Year == y & SEASON == s)
    INPUT.YEAR <- subset(temp_data_avail_df, YEAR == y & SEASON == s)$INPUT.YEAR
    
    if (sum(is.na(annual_df$BOTTEMP)) == 0) next
    
    # get the df to develop Kriging
    with_TEMP_df <- subset(subset(full_BIG_RF_df, Year == INPUT.YEAR), !is.na(BOTTEMP))
    sp::coordinates(with_TEMP_df) <-  ~LON + LAT
    
    # develop a variogram model
    temp.vgm <- variogram(BOTTEMP ~ LON + LAT, data = with_TEMP_df, width = 0.01)
    # temp.fit <- fit.variogram(temp.vgm, vgm(c( "Nug")))
    temp.fit <- fit.variogram(temp.vgm, vgm(c("Nug", "Sph","Mat", "Pen" , "Ste", "Cir" , "Bes", "Exc")))
    
    
    # create the spatial scale to be interpolated
    missing_TEMP_coord <- subset(annual_df, is.na(BOTTEMP))[, c("LAT", "LON")]
    sp::coordinates(missing_TEMP_coord) <-  ~LON + LAT
    
    # perform kriging interpolation for the tows without BOTTEMP
    temp.kriged <- krige(BOTTEMP ~ LON + LAT, with_TEMP_df, missing_TEMP_coord, model = temp.fit)
    
    # store the interpolated value
    missing_TEMP_df <- subset(annual_df, is.na(BOTTEMP))
    missing_TEMP_df$BOTTEMP <- temp.kriged$var1.pred
    
    no_TEMP_df <- rbind(no_TEMP_df, missing_TEMP_df)
    
  }
}

remove(y, s, annual_df, INPUT.YEAR, with_TEMP_df, temp.vgm, temp.fit, missing_TEMP_coord, temp.kriged, missing_TEMP_df)


# generate the final df for RF
full_BIG_RF_df <- rbind(subset(full_BIG_RF_df, !is.na(BOTTEMP)),
                        no_TEMP_df) %>%
                  arrange(ID)

## ----------------------------------- ##



  ## 2.4 fit spring model ---- 

set.seed(2) # ensure reproducibility
RF_BIG_spring <- randomForest(NUMBER ~ YEAR + LAT + BOTTEMP + AVGDEPTH, data = subset(full_BIG_RF_df, SEASON == "SPRING" & OWF == "OUTSIDE"), mtry = 2, ntree = 1000)

##  predict the results and save
full_BIG_RF_spring_df <- full_BIG_RF_df %>%
  filter(SEASON == "SPRING") %>%
  mutate(PREDICTED_NUMBER =  predict(RF_BIG_spring, newdata = ., type = "response")) %>%
  mutate(Final_NUMBER = case_when(OWF == "INSIDE" ~ PREDICTED_NUMBER,
                                  OWF == "OUTSIDE" ~  NUMBER))

## ----------------------------------- ##


  ## 2.5 fit fall model ---- 

set.seed(2) # ensure reproducibility
RF_BIG_fall <- randomForest(NUMBER ~ YEAR + LAT + BOTTEMP + AVGDEPTH, data = subset(full_BIG_RF_df, SEASON == "FALL" & OWF == "OUTSIDE"), mtry = 2, ntree = 1000)

##  predict the results and save
full_BIG_RF_fall_df <- full_BIG_RF_df %>%
  filter(SEASON == "FALL") %>%
  mutate(PREDICTED_NUMBER =  predict(RF_BIG_fall, newdata = ., type = "response")) %>%
  mutate(Final_NUMBER = case_when(OWF == "INSIDE" ~ PREDICTED_NUMBER,
                                  OWF == "OUTSIDE" ~  NUMBER))

## ----------------------------------- ##


  ## 2.6 combine and save

full_BIG_RF_df <- rbind(full_BIG_RF_spring_df, full_BIG_RF_fall_df)

write.csv(full_BIG_RF_df, "results/indices for assessment/summer flounder/mitigation_2_model-standardized/Full_BIG_RF_standardized_tow.csv", row.names = FALSE)


rm(list = ls())


## ----------------------------------- ##
