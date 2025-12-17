library(tidyverse)
library(splines)
library(Metrics)
library(VAST)


# this mitigation scenario will use model-based indices to replace the original indices
# since there are different indices used, we are fitting new models for abundance here
# note the these model-based indices will only use WEE dataset
# let's use VAST for abundance here


# 1. RDtrends (using TOTALN)  ----


  ## 1.1 generate full dataset  ----

AS_QQ_overlay_df <- read.csv("results/AS_OQ_tows_overlay_final.csv") %>%
  mutate(ID.temp = paste(CRUISE6, STATION, sep = "."))

tow_RDtrendS_df <- read.csv("data/NOAA.stock.data/surfclam/TOW_DATA.CSV") %>% 
  mutate(ID.temp = paste(CRUISE6, STATION, sep = ".")) %>% # using a new id here because the overlay ID have different stratum coding system
  filter(!ID.temp %in% unique(AS_QQ_overlay_df$ID.temp)) %>% 
  select(REGNAM, CRUISE6, YR, STATION, STRATUM, TOW, TOTALN, ASWEPT_M2, DISTANCEDETAIL, LON, LAT) %>%
  rename(YEAR = YR) %>%
  filter(YEAR <= 2011, REGNAM == "SVAtoSNE") %>%
  filter(nchar(TOW) != 9) %>% # remove the repeated tows, looks like the data are already filled in to the original record
  mutate(ID = paste(CRUISE6, STRATUM, STATION, sep = ".")) %>%
  distinct() %>% 
  group_by(YEAR, STRATUM) %>%
  mutate(TOTAL_N_STATION = length(STATION)) %>%
  ungroup() %>%
  filter(!is.na(TOTALN)) %>%
  mutate(LON = -LON)


## --------------------------------------------- ##



  ## 1.2 VAST model prep ----

RDtrendS_VAST_df <- tow_RDtrendS_df %>%  
  mutate(AreaSwept_km2 = ASWEPT_M2/(10^6)) %>%
  add_column(Vessel = NA, Pass = 0) %>% # required by the model, don't know why  
  rename(Lat = LAT, Lon = LON) %>% # the exact spelling is required by VAST
  mutate(Year = as.numeric(factor(YEAR)))

# quick visual

ggplot(RDtrendS_VAST_df) +
  geom_point(aes(x = Lon, y = Lat)) +
  facet_wrap(.~YEAR)

## --------------------------------------------- ##



  ## 1.3 prepare extrapolation area ----

extrapolation_region <- make_extrapolation_info(Region = "northwest_atlantic",  # use the embedded dataset
                                                strata.limits = data.frame(STRATA = "EPU"),
                                                epu_to_use = c("All", "Georges_Bank", "Mid_Atlantic_Bight", "Scotian_Shelf", "Gulf_of_Maine", "Other")[c(2,3)])

colnames(extrapolation_region$Data_Extrap)[4] <- "Area_km2"

ggplot(subset(extrapolation_region$Data_Extrap, Include == 1)) +
  geom_point(aes(x = Lon, y = Lat))
# plot(extrapolation_region)

## --------------------------------------------- ##




  ## 1.4 creating setting ----

settings <-  make_settings(n_x = c(100, 250, 500, 1000, 2000)[4], # number of knots, the more knots the longer it runs
                           purpose = "index2", # to calculate index
                           Region = "user",
                           # when region is user, then input_grid is needed be specified in the fit_model function as the extrapolation we generated earlier
                           # the setting for temporal modelign is: FieldConfig = c("Omega1" = 0, "Epsilon1" = 0, "Omega2" = 1, "Epsilon2" = 1),
                           FieldConfig = c("Omega1" = 1, "Epsilon1" = 1, "Omega2" = 1, "Epsilon2" = 1),
                           # specify various options for turning on/off spatial and spatial-temporal variation factors in the two linear predictors
                           # omega: spatial variation ; epsilon: spatiotemporal variation
                           # details see here: https://rdrr.io/github/James-Thorson/VAST/man/make_data.html
                           # !!! turned off first predictor for longfin squid
                           RhoConfig = c("Beta1" = 0, "Beta2" = 0, "Epsilon1" = 0, "Epsilon2" = 0),
                           # specify whether either intercepts (Beta1 and Beta2) or spatio-temporal variation (Epsilon1 and Epsilon2) is structured among time intervals
                           # untouched default setting: all is zero, so all are fixed effect
                           OverdispersionConfig= c(0,0),
                           # vessel/targeting effects for encounter probability and positive catch rate
                           # untouched default setting (0,0) means no effects
                           ObsModel = c(2,1),
                           # link functions for encounter probability and positive catch rate are available
                           # untouched default setting here
                           bias.correct = TRUE # to save time set as false, epsilon bias-correction
)

## --------------------------------------------- ##



  ## 1.5 run model with loop ----


scenario_RDtrendS_df <- expand.grid(error.dist = c(1,2,4,9), link.func = c(0,1), AIC = NA, RMSE = NA) ##details see "make_data" description


# c(4,1) with env.variable is the optimal model, but env messed up the indices, so we ditched the env.var
# error - delta-lognormal
# link - Poisson link  delta model 

# i = 7 # optimal model


for (i in 1:nrow(scenario_RDtrendS_df)) {
  
  # change the observation model setting, 
  settings$ObsModel <- c(scenario_RDtrendS_df$error.dist[i], scenario_RDtrendS_df$link.func[i])
  
  
  # run model with tryCatch to skip errors when scenario does not converge
  
  VAST_RDtrendS <- tryCatch(
    # attempt to fit the model
    {
      VAST_RDtrendS <- fit_model(settings = settings,
                                      input_grid = extrapolation_region$Data_Extrap, 
                                      Lat_i = RDtrendS_VAST_df$Lat,
                                      Lon_i = RDtrendS_VAST_df$Lon,
                                      t_i = as.numeric(RDtrendS_VAST_df$Year), 
                                      b_i = RDtrendS_VAST_df$TOTALN,
                                      a_i = as_units(RDtrendS_VAST_df$AreaSwept_km2, "km^2"),
                                      working_dir = "results/indices for assessment/surfclam/mitigation_3_model-based indices/VAST_RDtrendS")
      VAST_RDtrendS # return the fit model if successful
    },
    
    # return NULL if error
    error = function(e) { NULL }
  )
  
  # skip to the next iteration if NULL is returned
  if (is.null(VAST_RDtrendS)) { next }
  if (length(VAST_RDtrendS$Report) == 1) {
    if (VAST_RDtrendS$Report == "Model is not converged") { next }
  }
  
  scenario_RDtrendS_df$AIC[i] = VAST_RDtrendS$parameter_estimates$AIC
  
  save(VAST_RDtrendS, file = paste0("results/indices for assessment/surfclam/mitigation_3_model-based indices/VAST_RDtrendS/VAST_RDtrendS_",
                                    scenario_RDtrendS_df$error.dist[i], ".",scenario_RDtrendS_df$link.func[i], "_model.Rdata"))
  
  write.csv(scenario_RDtrendS_df, "results/indices for assessment/surfclam/mitigation_3_model-based indices/VAST_RDtrendS/scenario_RDtrendS.csv", row.names = FALSE)
  
}



## --------------------------------------------- ##

remove(list=ls())


  ## 1.6. extract abundance indices ----

load("results/indices for assessment/surfclam/mitigation_3_model-based indices/VAST_RDtrendS/final model/VAST_RDtrendS_model.Rdata")


VAST_RDtrendS_data <- plot(VAST_RDtrendS, working_dir = "results/indices for assessment/surfclam/mitigation_3_model-based indices/VAST_RDtrendS/final model/figure/")

VAST_RDtrendS_df <- data.frame(YEAR = VAST_RDtrendS_data$Index$Table$Time,
                           Season = NA,
                           MODEL = "SVAtoSNE",
                           fit = VAST_RDtrendS_data$Index$Table$Estimate, 
                           lwr = VAST_RDtrendS_data$Index$Table$Estimate - 1.96 * VAST_RDtrendS_data$Index$Table$`Std. Error for Estimate`, 
                           upr = VAST_RDtrendS_data$Index$Table$Estimate + 1.96 * VAST_RDtrendS_data$Index$Table$`Std. Error for Estimate`) 



write.csv(VAST_RDtrendS_df, "results/indices for assessment/surfclam/mitigation_3_model-based indices/VAST_RDtrendS/final model/VAST_RDtrendS_Indices.csv", row.names = FALSE)

## --------------------------------------------- ##



# -------------------------------------------------------------------------------------------------------- #




remove(list=ls())



# 2. RDscales (using NPERTOW)  ----



  ## 2.1 generate full dataset  ----

AS_QQ_overlay_df <- read.csv("results/AS_OQ_tows_overlay_final.csv") %>%
  mutate(ID.temp = paste(CRUISE6, STATION, sep = "."))

tow_RDscaleS_df <- read.csv("results/stratified.mean.indices/surfclam/total.catch.by.tow.csv")%>%
  filter(YEAR %in% c(1997, 1999, 2002, 2005, 2008, 2011), REGION == "SVAtoSNE") %>% 
  mutate(ID.temp = paste(CRUISE6, STATION, sep = ".")) %>% 
  filter(!ID.temp %in% unique(AS_QQ_overlay_df$ID.temp))

remove(AS_QQ_overlay_df)


## --------------------------------------------- ##



  ## 2.2 VAST model prep ----

RDscaleS_VAST_df <- tow_RDscaleS_df %>%  
  mutate(AreaSwept_km2 = ASWEPT_M2/(10^6)) %>%
  add_column(Vessel = NA, Pass = 0) %>% # required by the model, don't know why  
  rename(Lat = LAT, Lon = LON) %>% # the exact spelling is required by VAST
  mutate(Year = as.numeric(factor(YEAR)))

# quick visual

ggplot(RDscaleS_VAST_df) +
  geom_point(aes(x = Lon, y = Lat)) +
  facet_wrap(.~YEAR)

## --------------------------------------------- ##




  ## 2.3 prepare extrapolation area ----

extrapolation_region <- make_extrapolation_info(Region = "northwest_atlantic",  # use the embedded dataset
                                                strata.limits = data.frame(STRATA = "EPU"),
                                                epu_to_use = c("All", "Georges_Bank", "Mid_Atlantic_Bight", "Scotian_Shelf", "Gulf_of_Maine", "Other")[c(2,3)])

colnames(extrapolation_region$Data_Extrap)[4] <- "Area_km2"

ggplot(subset(extrapolation_region$Data_Extrap, Include == 1)) +
  geom_point(aes(x = Lon, y = Lat))
# plot(extrapolation_region)

## --------------------------------------------- ##




  ## 2.4 creating setting ----

settings <-  make_settings(n_x = c(100, 250, 500, 1000, 2000)[4], # number of knots, the more knots the longer it runs
                           purpose = "index2", # to calculate index
                           Region = "user",
                           # when region is user, then input_grid is needed be specified in the fit_model function as the extrapolation we generated earlier
                           # the setting for temporal modelign is: FieldConfig = c("Omega1" = 0, "Epsilon1" = 0, "Omega2" = 1, "Epsilon2" = 1),
                           FieldConfig = c("Omega1" = 1, "Epsilon1" = 1, "Omega2" = 1, "Epsilon2" = 1),
                           # specify various options for turning on/off spatial and spatial-temporal variation factors in the two linear predictors
                           # omega: spatial variation ; epsilon: spatiotemporal variation
                           # details see here: https://rdrr.io/github/James-Thorson/VAST/man/make_data.html
                           # !!! turned off first predictor for longfin squid
                           RhoConfig = c("Beta1" = 0, "Beta2" = 0, "Epsilon1" = 0, "Epsilon2" = 0),
                           # specify whether either intercepts (Beta1 and Beta2) or spatio-temporal variation (Epsilon1 and Epsilon2) is structured among time intervals
                           # untouched default setting: all is zero, so all are fixed effect
                           OverdispersionConfig= c(0,0),
                           # vessel/targeting effects for encounter probability and positive catch rate
                           # untouched default setting (0,0) means no effects
                           ObsModel = c(2,1),
                           # link functions for encounter probability and positive catch rate are available
                           # untouched default setting here
                           bias.correct = TRUE # to save time set as false, epsilon bias-correction
)

## --------------------------------------------- ##



  ## 2.5 run model with loop ----


scenario_RDscaleS_df <- expand.grid(error.dist = c(1,2,4,9), link.func = c(0,1), AIC = NA, RMSE = NA) ##details see "make_data" description


# c(9,1) with env.variable is the optimal model, but env messed up the indices, so we ditched the env.var
# error - delta-generalized gamma
# link - Poisson link  delta model 

# i = 8 # optimal model


for (i in 1:nrow(scenario_RDscaleS_df)) {
  
  # change the observation model setting, 
  settings$ObsModel <- c(scenario_RDscaleS_df$error.dist[i], scenario_RDscaleS_df$link.func[i])
  
  
  # run model with tryCatch to skip errors when scenario does not converge
  
  VAST_RDscaleS <- tryCatch(
    # attempt to fit the model
    {
      VAST_RDscaleS <- fit_model(settings = settings,
                                 input_grid = extrapolation_region$Data_Extrap, 
                                 Lat_i = RDscaleS_VAST_df$Lat,
                                 Lon_i = RDscaleS_VAST_df$Lon,
                                 t_i = as.numeric(RDscaleS_VAST_df$Year), 
                                 b_i = RDscaleS_VAST_df$NPERTOW,
                                 a_i = as_units(RDscaleS_VAST_df$AreaSwept_km2, "km^2"),
                                 working_dir = "results/indices for assessment/surfclam/mitigation_3_model-based indices/VAST_RDscaleS")
      VAST_RDscaleS # return the fit model if successful
    },
    
    # return NULL if error
    error = function(e) { NULL }
  )
  
  # skip to the next iteration if NULL is returned
  if (is.null(VAST_RDscaleS)) { next }
  if (length(VAST_RDscaleS$Report) == 1) {
    if (VAST_RDscaleS$Report == "Model is not converged") { next }
  }
  
  scenario_RDscaleS_df$AIC[i] = VAST_RDscaleS$parameter_estimates$AIC
  
  save(VAST_RDscaleS, file = paste0("results/indices for assessment/surfclam/mitigation_3_model-based indices/VAST_RDscaleS/VAST_RDscaleS_",
                                    scenario_RDscaleS_df$error.dist[i], ".",scenario_RDscaleS_df$link.func[i], "_model.Rdata"))
  
  write.csv(scenario_RDscaleS_df, "results/indices for assessment/surfclam/mitigation_3_model-based indices/VAST_RDscaleS/scenario_RDscaleS.csv", row.names = FALSE)
  
}



## --------------------------------------------- ##

remove(list=ls())


  ## 2.6. extract abundance indices ----

load("results/indices for assessment/surfclam/mitigation_3_model-based indices/VAST_RDscaleS/final model/VAST_RDscaleS_model.Rdata")


VAST_RDscaleS_data <- plot(VAST_RDscaleS, working_dir = "results/indices for assessment/surfclam/mitigation_3_model-based indices/VAST_RDscaleS/final model/figure/")

VAST_RDscaleS_df <- data.frame(YEAR = VAST_RDscaleS_data$Index$Table$Time,
                               Season = NA,
                               MODEL = "SVAtoSNE",
                               fit = VAST_RDscaleS_data$Index$Table$Estimate, 
                               lwr = VAST_RDscaleS_data$Index$Table$Estimate - 1.96 * VAST_RDscaleS_data$Index$Table$`Std. Error for Estimate`, 
                               upr = VAST_RDscaleS_data$Index$Table$Estimate + 1.96 * VAST_RDscaleS_data$Index$Table$`Std. Error for Estimate`) 



write.csv(VAST_RDscaleS_df, "results/indices for assessment/surfclam/mitigation_3_model-based indices/VAST_RDscaleS/final model/VAST_RDscaleS_Indices.csv", row.names = FALSE)





# -------------------------------------------------------------------------------------------------------- #


remove(list = ls())



# 3. MCDS (using NPERTOW)  ----



  ## 3.1 generate full dataset  ----

AS_QQ_overlay_df <- read.csv("results/AS_OQ_tows_overlay_final.csv") %>%
  mutate(ID.temp = paste(CRUISE6, STATION, sep = "."))

tow_MCDS_df <- read.csv("results/stratified.mean.indices/surfclam/total.catch.by.tow.csv")%>%
  filter(YEAR %in% c(2012, 2015, 2018, 2022), REGION == "SVAtoSNE") %>% 
  mutate(ID.temp = paste(CRUISE6, STATION, sep = ".")) %>% 
  filter(!ID.temp %in% unique(AS_QQ_overlay_df$ID.temp))

remove(AS_QQ_overlay_df)

tow_MCDS_df[tow_MCDS_df$ID == "201560.5S.140" ,]$LON <- mean(c(-74.9333, -73.7823, -75.4405, -74.2491, -73.7816, -74.3187)) # fill the missing Lon with average of station 140


## --------------------------------------------- ##



  ## 3.2 VAST model prep ----

MCDS_VAST_df <- tow_MCDS_df %>%  
  mutate(AreaSwept_km2 = ASWEPT_M2/(10^6)) %>%
  add_column(Vessel = NA, Pass = 0) %>% # required by the model, don't know why  
  rename(Lat = LAT, Lon = LON) %>% # the exact spelling is required by VAST
  mutate(Year = as.numeric(factor(YEAR)))

# quick visual

ggplot(MCDS_VAST_df) +
  geom_point(aes(x = Lon, y = Lat)) +
  facet_wrap(.~YEAR)

## --------------------------------------------- ##




  ## 3.3 prepare extrapolation area ----

extrapolation_region <- make_extrapolation_info(Region = "northwest_atlantic",  # use the embedded dataset
                                                strata.limits = data.frame(STRATA = "EPU"),
                                                epu_to_use = c("All", "Georges_Bank", "Mid_Atlantic_Bight", "Scotian_Shelf", "Gulf_of_Maine", "Other")[c(2,3)])

colnames(extrapolation_region$Data_Extrap)[4] <- "Area_km2"

ggplot(subset(extrapolation_region$Data_Extrap, Include == 1)) +
  geom_point(aes(x = Lon, y = Lat))
# plot(extrapolation_region)

## --------------------------------------------- ##




  ## 3.4 creating setting ----

settings <-  make_settings(n_x = c(100, 250, 500, 1000, 2000)[4], # number of knots, the more knots the longer it runs
                           purpose = "index2", # to calculate index
                           Region = "user",
                           # when region is user, then input_grid is needed be specified in the fit_model function as the extrapolation we generated earlier
                           # the setting for temporal modelign is: FieldConfig = c("Omega1" = 0, "Epsilon1" = 0, "Omega2" = 1, "Epsilon2" = 1),
                           FieldConfig = c("Omega1" = 1, "Epsilon1" = 1, "Omega2" = 1, "Epsilon2" = 1),
                           # specify various options for turning on/off spatial and spatial-temporal variation factors in the two linear predictors
                           # omega: spatial variation ; epsilon: spatiotemporal variation
                           # details see here: https://rdrr.io/github/James-Thorson/VAST/man/make_data.html
                           # !!! turned off first predictor for longfin squid
                           RhoConfig = c("Beta1" = 0, "Beta2" = 0, "Epsilon1" = 0, "Epsilon2" = 0),
                           # specify whether either intercepts (Beta1 and Beta2) or spatio-temporal variation (Epsilon1 and Epsilon2) is structured among time intervals
                           # untouched default setting: all is zero, so all are fixed effect
                           OverdispersionConfig= c(0,0),
                           # vessel/targeting effects for encounter probability and positive catch rate
                           # untouched default setting (0,0) means no effects
                           ObsModel = c(2,1),
                           # link functions for encounter probability and positive catch rate are available
                           # untouched default setting here
                           bias.correct = TRUE # to save time set as false, epsilon bias-correction
)

## --------------------------------------------- ##



  ## 3.5 run model with loop ----


scenario_MCDS_df <- expand.grid(error.dist = c(1,2,4,9), link.func = c(0,1), AIC = NA, RMSE = NA) ##details see "make_data" description


# c(2,0) with env.variable is the optimal model, but env messed up the indices, so we ditched the env.var
# error - delta gamma
# link - Conventional delta-model using logit-link for encounter probability and log-link for positive catch rates

# i = 2 # optimal model


for (i in 1:nrow(scenario_MCDS_df)) {
  
  # change the observation model setting, 
  settings$ObsModel <- c(scenario_MCDS_df$error.dist[i], scenario_MCDS_df$link.func[i])
  
  
  # run model with tryCatch to skip errors when scenario does not converge
  
  VAST_MCDS <- tryCatch(
    # attempt to fit the model
    {
      VAST_MCDS <- fit_model(settings = settings,
                                 input_grid = extrapolation_region$Data_Extrap, 
                                 Lat_i = MCDS_VAST_df$Lat,
                                 Lon_i = MCDS_VAST_df$Lon,
                                 t_i = as.numeric(MCDS_VAST_df$Year), 
                                 b_i = MCDS_VAST_df$NPERTOW,
                                 a_i = as_units(MCDS_VAST_df$AreaSwept_km2, "km^2"),
                                 working_dir = "results/indices for assessment/surfclam/mitigation_3_model-based indices/VAST_MCDS")
      VAST_MCDS # return the fit model if successful
    },
    
    # return NULL if error
    error = function(e) { NULL }
  )
  
  # skip to the next iteration if NULL is returned
  if (is.null(VAST_MCDS)) { next }
  if (length(VAST_MCDS$Report) == 1) {
    if (VAST_MCDS$Report == "Model is not converged") { next }
  }
  
  scenario_MCDS_df$AIC[i] = VAST_MCDS$parameter_estimates$AIC
  
  save(VAST_MCDS, file = paste0("results/indices for assessment/surfclam/mitigation_3_model-based indices/VAST_MCDS/VAST_MCDS_",
                                    scenario_MCDS_df$error.dist[i], ".",scenario_MCDS_df$link.func[i], "_model.Rdata"))
  
  write.csv(scenario_MCDS_df, "results/indices for assessment/surfclam/mitigation_3_model-based indices/VAST_MCDS/scenario_MCDS.csv", row.names = FALSE)
  
}



## --------------------------------------------- ##

remove(list=ls())


  ## 3.6. extract abundance indices ----

load("results/indices for assessment/surfclam/mitigation_3_model-based indices/VAST_MCDS/final model/VAST_MCDS_model.Rdata")


VAST_MCDS_data <- plot(VAST_MCDS, working_dir = "results/indices for assessment/surfclam/mitigation_3_model-based indices/VAST_MCDS/final model/figure/")

VAST_MCDS_df <- data.frame(YEAR = VAST_MCDS_data$Index$Table$Time,
                               Season = NA,
                               MODEL = "SVAtoSNE",
                               fit = VAST_MCDS_data$Index$Table$Estimate, 
                               lwr = VAST_MCDS_data$Index$Table$Estimate - 1.96 * VAST_MCDS_data$Index$Table$`Std. Error for Estimate`, 
                               upr = VAST_MCDS_data$Index$Table$Estimate + 1.96 * VAST_MCDS_data$Index$Table$`Std. Error for Estimate`) 



write.csv(VAST_MCDS_df, "results/indices for assessment/surfclam/mitigation_3_model-based indices/VAST_MCDS/final model/VAST_MCDS_Indices.csv", row.names = FALSE)





# -------------------------------------------------------------------------------------------------------- #


remove(list = ls())


# 4. reformat to stratified mean value ---------------------------------------------------------------------------------------------


 ## 4.1 RDtrendS ----

VAST_RDtrendS_df <- read.csv("results/indices for assessment/surfclam/mitigation_3_model-based indices/VAST_RDtrendS/final model/VAST_RDtrendS_Indices.csv") %>%
  select(-c(MODEL, Season)) %>%
  mutate(YEAR = c(1982, 1983, 1984, 1986, 1989, 1992, 1994, 1997, 1999, 2002, 2005, 2008, 2011)) 


# extract the total area size from VAST to estimate density

load("results/indices for assessment/surfclam/mitigation_3_model-based indices/VAST_RDtrendS/final model/VAST_RDtrendS_model.Rdata")

VAST_extra_area <- as.numeric(sum(VAST_RDtrendS$extrapolation_list$Area_km2_x)) * 10^6

VAST_RDtrendS_df <- VAST_RDtrendS_df  %>%
  mutate(fit = fit/VAST_extra_area * 1000,
         lwr = lwr/VAST_extra_area * 1000,
         upr = upr/VAST_extra_area * 1000)

write.csv(VAST_RDtrendS_df, "results/indices for assessment/surfclam/mitigation_3_model-based indices/VAST_RDtrendS/final model/VAST_RDtrendS_Indices.csv", row.names = FALSE)

## --------------------------------------------- ##



  ## 4.2 RDscaleS ----

VAST_RDscaleS_df <- read.csv("results/indices for assessment/surfclam/mitigation_3_model-based indices/VAST_RDscaleS/final model/VAST_RDscaleS_Indices.csv") %>%
  select(-c(MODEL, Season)) %>%
  mutate(YEAR = c(1997, 1999, 2002, 2005, 2008, 2011)) 


# extract the total area size from VAST to estimate density

load("results/indices for assessment/surfclam/mitigation_3_model-based indices/VAST_RDscaleS/final model/VAST_RDscaleS_model.Rdata")

VAST_extra_area <- as.numeric(sum(VAST_RDscaleS$extrapolation_list$Area_km2_x)) * 10^6

VAST_RDscaleS_df <- VAST_RDscaleS_df  %>%
  mutate(fit = fit/VAST_extra_area * 1000,lwr = lwr/VAST_extra_area * 1000,
         upr = upr/VAST_extra_area * 1000)

write.csv(VAST_RDscaleS_df, "results/indices for assessment/surfclam/mitigation_3_model-based indices/VAST_RDscaleS/final model/VAST_RDscaleS_Indices.csv", row.names = FALSE)


## --------------------------------------------- ##


## 4.3 MCDS ----

VAST_MCDS_df <- read.csv("results/indices for assessment/surfclam/mitigation_3_model-based indices/VAST_MCDS/final model/VAST_MCDS_Indices.csv") %>%
  select(-c(MODEL, Season)) %>%
  mutate(YEAR = c(2012, 2015, 2018, 2022)) 


# extract the total area size from VAST to estimate density

load("results/indices for assessment/surfclam/mitigation_3_model-based indices/VAST_MCDS/final model/VAST_MCDS_model.Rdata")

VAST_extra_area <- as.numeric(sum(VAST_MCDS$extrapolation_list$Area_km2_x)) * 10^6

VAST_MCDS_df <- VAST_MCDS_df  %>%
  mutate(fit = fit/VAST_extra_area * 1000,lwr = lwr/VAST_extra_area * 1000,
         upr = upr/VAST_extra_area * 1000)

write.csv(VAST_MCDS_df, "results/indices for assessment/surfclam/mitigation_3_model-based indices/VAST_MCDS/final model/VAST_MCDS_Indices.csv", row.names = FALSE)

# -------------------------------------------------------------------------------------------------------- #



