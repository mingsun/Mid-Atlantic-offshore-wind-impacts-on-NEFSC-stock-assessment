library(tidyverse)
library(VAST)
library(splines)
library(Metrics)

# this mitigation scenario will use model-based indices to replace the original indices
# since there are different indices used, new models will need to be fitted here
# note the these model-based indices will only use WEE dataset
# according to the first paper, the best model is VAST for all scenario


# 1. ALB (1982-2008) ----


  ## 1.1 generate full dataset with env.var ----

  cat_df <- read.csv("results/indices for assessment/summer flounder/catch_by_tow_ALB.csv")


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

  ## ---------------------------------------  ##



  ## 1.2 VAST model prep ----
  
full_VAST_df <- full_df %>%   # convert YEAR into a factor
  filter(OWF == "OUTSIDE") 

  ## ---------------------------------------  ##



  ## 1.3 Fall model ----

fall_VAST_df <- subset(full_VAST_df, SEASON == "FALL")


  ### 1.3.1 prepare extrapolation area ----

extrapolation_region <- make_extrapolation_info(Region = "northwest_atlantic",  # use the embedded dataset
                                                strata.limits = data.frame(STRATA = "EPU"),
                                                epu_to_use = c("All", "Georges_Bank", "Mid_Atlantic_Bight", "Scotian_Shelf", "Gulf_of_Maine", "Other")[c(2,3)],
                                                # if use "Georges_Bank", "Mid_Atlantic_Bight", there will be 436 rows missing for summer flounder
                                                # all outside OWF, not affecting the gap analysis
                                                DirPath = "results/indices for assessment/summer flounder/mitigation_3_model-based indices/VAST run/")

colnames(extrapolation_region$Data_Extrap)[4] <- "Area_km2"

ggplot(subset(extrapolation_region$Data_Extrap, Include == 1)) +
  geom_point(aes(x = Lon, y = Lat))

  ### ---------------------------------------  ###


      ### 1.3.2 prepare settings and covariate & formula ----

## create the setting file
settings <-  make_settings(n_x = c(100, 250, 500, 1000, 2000)[4], # number of knots, the more knots the longer it runs
                           purpose = "index2", # to calculate index
                           Region = "user", 
                           # when region is user, then input_grid is needed be specified in the fit_model function as the extrapolation we generated earlier
                           # FieldConfig = c("Omega1" = 1, "Epsilon1" = 1, "Omega2" = 0, "Epsilon2" = 0), 
                           FieldConfig = c("Omega1" = 1, "Epsilon1" = 1, "Omega2" = 1, "Epsilon2" = 1),
                           # specify various options for turning on/off spatial and spatial-temporal variation factors in the two linear predictors
                           # omega: spatial variation ; epsilon: spatiotemporal variation
                           # details see here: https://rdrr.io/github/James-Thorson/VAST/man/make_data.html
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





  ### ---------------------------------------  ###


    ### 1.3.3 loop over multiple error/link options to identify the best model ----

# only error dist = c(1,2,4,9) have model fitted
scenario_df <- expand.grid(error.dist = c(1,2,4,9), link.func = c(0:1), AIC = NA, RMSE = NA, finished  = NA) ##details see "make_data" description

# c(9,1) with env.variable is the optimal model, but env messed up the indices, so we ditched the env.var
# error - delta-generalized gamma
# link - Poisson link  delta model 

i = 8 # optimal model

for (i in 1:nrow(scenario_df)) {
  scenario_df$finished[i] = "yes"
  write.csv(scenario_df, "results/indices for assessment/summer flounder/mitigation_3_model-based indices/VAST run/scenario.csv", row.names = FALSE)
  
  # change the observation model setting, 
  settings$ObsModel <- c(scenario_df$error.dist[i], scenario_df$link.func[i])
  
  # run model with tryCatch to skip errors when scenario does not converge
  
  VAST_temp_fall <- tryCatch(
    # attempt to fit the model
    {
      VAST_temp_fall <- fit_model(settings = settings,
                                  input_grid = extrapolation_region$Data_Extrap, 
                                  Lat_i = fall_VAST_df$Lat,
                                  Lon_i = fall_VAST_df$Lon,
                                  t_i = as.numeric(fall_VAST_df$Year), 
                                  b_i = as_units(fall_VAST_df$NUMBER, "count"),
                                  a_i = as_units(fall_VAST_df$AreaSwept_km2, "km^2"),
                                  working_dir = "results/indices for assessment/summer flounder/mitigation_3_model-based indices/VAST run")
      VAST_temp_fall # return the fit model if successful
    },
    
    # return NULL if error
    error = function(e) { NULL }
  )
  
  # skip to the next iteration if NULL is returned
  if (is.null(VAST_temp_fall)) { next }
  if (length(VAST_temp_fall$Report) == 1) {
    if (VAST_temp_fall$Report == "Model is not converged") { next }
  }
  
  scenario_df$AIC[i] = VAST_temp_fall$parameter_estimates$AIC
  scenario_df$finished[i] = "yes"
  
  save(VAST_temp_fall, file = paste0("results/indices for assessment/summer flounder/mitigation_3_model-based indices/VAST run/VAST_fall_",
                                     scenario_df$error.dist[i], ".",scenario_df$link.func[i], "_model.Rdata"))
  
  write.csv(scenario_df, "results/indices for assessment/summer flounder/mitigation_3_model-based indices/VAST run/scenario.csv", row.names = FALSE)
  
}

## ---------------------------------------  ##


rm(list = setdiff(ls(), "full_VAST_df"))



  ## 1.4 Spring model ----

spring_VAST_df <- subset(full_VAST_df, SEASON == "SPRING")


    ### 1.4.1 prepare extrapolation area ----

extrapolation_region <- make_extrapolation_info(Region = "northwest_atlantic",  # use the embedded dataset
                                                strata.limits = data.frame(STRATA = "EPU"),
                                                epu_to_use = c("All", "Georges_Bank", "Mid_Atlantic_Bight", "Scotian_Shelf", "Gulf_of_Maine", "Other")[c(2,3)],
                                                # if use "Georges_Bank", "Mid_Atlantic_Bight", there will be 436 rows missing for summer flounder
                                                # all outside OWF, not affecting the gap analysis
                                                DirPath = "results/indices for assessment/summer flounder/mitigation_3_model-based indices/VAST run/")

colnames(extrapolation_region$Data_Extrap)[4] <- "Area_km2"

ggplot(subset(extrapolation_region$Data_Extrap, Include == 1)) +
  geom_point(aes(x = Lon, y = Lat))

    ### ---------------------------------------  ###


    ### 1.4.2 prepare settings and covariate & formula ----

## create the setting file
settings <-  make_settings(n_x = c(100, 250, 500, 1000, 2000)[4], # number of knots, the more knots the longer it runs
                           purpose = "index2", # to calculate index
                           Region = "user", 
                           # when region is user, then input_grid is needed be specified in the fit_model function as the extrapolation we generated earlier
                           # FieldConfig = c("Omega1" = 1, "Epsilon1" = 1, "Omega2" = 0, "Epsilon2" = 0), 
                           FieldConfig = c("Omega1" = 1, "Epsilon1" = 1, "Omega2" = 1, "Epsilon2" = 1),
                           # specify various options for turning on/off spatial and spatial-temporal variation factors in the two linear predictors
                           # omega: spatial variation ; epsilon: spatiotemporal variation
                           # details see here: https://rdrr.io/github/James-Thorson/VAST/man/make_data.html
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





    ### ---------------------------------------  ###


    ### 1.4.3 loop over multiple error/link options to identify the best model ----

# only error dist = c(1,2,4,9) have model fitted
scenario_df <- expand.grid(error.dist = c(1,2,4,9), link.func = c(0:1), AIC = NA, RMSE = NA, finished  = NA) ##details see "make_data" description

# c(4,1) with env.variable is the optimal model, but env messed up the indices, so we ditched the env.var
# error - delta-lognormal
# link - Poisson link  delta model

i = 7 # optimal model

for (i in 1:nrow(scenario_df)) {
  scenario_df$finished[i] = "yes"
  write.csv(scenario_df, "results/indices for assessment/summer flounder/mitigation_3_model-based indices/VAST run/scenario.csv", row.names = FALSE)
  
  # change the observation model setting, 
  settings$ObsModel <- c(scenario_df$error.dist[i], scenario_df$link.func[i])
  
  # run model with tryCatch to skip errors when scenario does not converge
  
  VAST_temp_spring <- tryCatch(
    # attempt to fit the model
    {
      VAST_temp_spring <- fit_model(settings = settings,
                                  input_grid = extrapolation_region$Data_Extrap, 
                                  Lat_i = spring_VAST_df$Lat,
                                  Lon_i = spring_VAST_df$Lon,
                                  t_i = as.numeric(spring_VAST_df$Year), 
                                  b_i = as_units(spring_VAST_df$NUMBER, "count"),
                                  a_i = as_units(spring_VAST_df$AreaSwept_km2, "km^2"),
                                  working_dir = "results/indices for assessment/summer flounder/mitigation_3_model-based indices/VAST run")
      VAST_temp_spring # return the fit model if successful
    },
    
    # return NULL if error
    error = function(e) { NULL }
  )
  
  # skip to the next iteration if NULL is returned
  if (is.null(VAST_temp_spring)) { next }
  if (length(VAST_temp_spring$Report) == 1) {
    if (VAST_temp_spring$Report == "Model is not converged") { next }
  }
  
  scenario_df$AIC[i] = VAST_temp_spring$parameter_estimates$AIC
  scenario_df$finished[i] = "yes"
  
  save(VAST_temp_spring, file = paste0("results/indices for assessment/summer flounder/mitigation_3_model-based indices/VAST run/VAST_spring_",
                                     scenario_df$error.dist[i], ".",scenario_df$link.func[i], "_model.Rdata"))
  
  write.csv(scenario_df, "results/indices for assessment/summer flounder/mitigation_3_model-based indices/VAST run/scenario.csv", row.names = FALSE)
  
}



# ---------------------------------------------------------------------------  #

rm(list = ls())



# 2. BIG (2009-2022) ----


  ## 2.1 generate full dataset with env.var ----

cat_df <- read.csv("results/indices for assessment/summer flounder/catch_by_tow_BIG.csv")


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

  ## ---------------------------------------  ##



  ## 2.2 VAST model prep ----

full_VAST_df <- full_df %>%   # convert YEAR into a factor
  filter(OWF == "OUTSIDE") 

## ---------------------------------------  ##



  ## 2.3 Fall model ----

fall_VAST_df <- subset(full_VAST_df, SEASON == "FALL")


    ### 2.3.1 prepare extrapolation area ----

extrapolation_region <- make_extrapolation_info(Region = "northwest_atlantic",  # use the embedded dataset
                                                strata.limits = data.frame(STRATA = "EPU"),
                                                epu_to_use = c("All", "Georges_Bank", "Mid_Atlantic_Bight", "Scotian_Shelf", "Gulf_of_Maine", "Other")[c(2,3)],
                                                # if use "Georges_Bank", "Mid_Atlantic_Bight", there will be 436 rows missing for summer flounder
                                                # all outside OWF, not affecting the gap analysis
                                                DirPath = "results/indices for assessment/summer flounder/mitigation_3_model-based indices/VAST run/")

colnames(extrapolation_region$Data_Extrap)[4] <- "Area_km2"

ggplot(subset(extrapolation_region$Data_Extrap, Include == 1)) +
  geom_point(aes(x = Lon, y = Lat))

    ### ---------------------------------------  ###


    ### 2.3.2 prepare settings and covariate & formula ----

      ## create the setting file
settings <-  make_settings(n_x = c(100, 250, 500, 1000, 2000)[4], # number of knots, the more knots the longer it runs
                           purpose = "index2", # to calculate index
                           Region = "user", 
                           # when region is user, then input_grid is needed be specified in the fit_model function as the extrapolation we generated earlier
                           # FieldConfig = c("Omega1" = 1, "Epsilon1" = 1, "Omega2" = 0, "Epsilon2" = 0), 
                           FieldConfig = c("Omega1" = 1, "Epsilon1" = 1, "Omega2" = 1, "Epsilon2" = 1),
                           # specify various options for turning on/off spatial and spatial-temporal variation factors in the two linear predictors
                           # omega: spatial variation ; epsilon: spatiotemporal variation
                           # details see here: https://rdrr.io/github/James-Thorson/VAST/man/make_data.html
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




    ### ---------------------------------------  ###


    ### 2.3.3 loop over multiple error/link options to identify the best model ----

# only error dist = c(1,2,4,9) have model fitted
scenario_df <- expand.grid(error.dist = c(1,2,4,9), link.func = c(0:1), AIC = NA, RMSE = NA, finished  = NA) ##details see "make_data" description

# c(4,1) with env.variable is the optimal model, but env messed up the indices, so we ditched the env.var
# error - delta-lognormal
# link - Poisson link  delta model

i = 7 # optimal model

for (i in 1:nrow(scenario_df)) {
  scenario_df$finished[i] = "yes"
  write.csv(scenario_df, "results/indices for assessment/summer flounder/mitigation_3_model-based indices/VAST run/scenario.csv", row.names = FALSE)
  
  # change the observation model setting, 
  settings$ObsModel <- c(scenario_df$error.dist[i], scenario_df$link.func[i])
  
  # run model with tryCatch to skip errors when scenario does not converge
  
  VAST_temp_fall <- tryCatch(
    # attempt to fit the model
    {
      VAST_temp_fall <- fit_model(settings = settings,
                                  input_grid = extrapolation_region$Data_Extrap, 
                                  Lat_i = fall_VAST_df$Lat,
                                  Lon_i = fall_VAST_df$Lon,
                                  t_i = as.numeric(fall_VAST_df$Year), 
                                  b_i = as_units(fall_VAST_df$NUMBER, "count"),
                                  a_i = as_units(fall_VAST_df$AreaSwept_km2, "km^2"),
                                  working_dir = "results/indices for assessment/summer flounder/mitigation_3_model-based indices/VAST run")
      VAST_temp_fall # return the fit model if successful
    },
    
    # return NULL if error
    error = function(e) { NULL }
  )
  
  # skip to the next iteration if NULL is returned
  if (is.null(VAST_temp_fall)) { next }
  if (length(VAST_temp_fall$Report) == 1) {
    if (VAST_temp_fall$Report == "Model is not converged") { next }
  }
  
  scenario_df$AIC[i] = VAST_temp_fall$parameter_estimates$AIC
  scenario_df$finished[i] = "yes"
  
  save(VAST_temp_fall, file = paste0("results/indices for assessment/summer flounder/mitigation_3_model-based indices/VAST run/VAST_fall_",
                                     scenario_df$error.dist[i], ".",scenario_df$link.func[i], "_model.Rdata"))
  
  write.csv(scenario_df, "results/indices for assessment/summer flounder/mitigation_3_model-based indices/VAST run/scenario.csv", row.names = FALSE)
  
}

  ## ---------------------------------------  ##


rm(list = setdiff(ls(), "full_VAST_df"))



  ## 2.4 Spring model ----

spring_VAST_df <- subset(full_VAST_df, SEASON == "SPRING")


    ### 2.4.1 prepare extrapolation area ----

extrapolation_region <- make_extrapolation_info(Region = "northwest_atlantic",  # use the embedded dataset
                                                strata.limits = data.frame(STRATA = "EPU"),
                                                epu_to_use = c("All", "Georges_Bank", "Mid_Atlantic_Bight", "Scotian_Shelf", "Gulf_of_Maine", "Other")[c(2,3)],
                                                # if use "Georges_Bank", "Mid_Atlantic_Bight", there will be 436 rows missing for summer flounder
                                                # all outside OWF, not affecting the gap analysis
                                                DirPath = "results/indices for assessment/summer flounder/mitigation_3_model-based indices/VAST run/")

colnames(extrapolation_region$Data_Extrap)[4] <- "Area_km2"

ggplot(subset(extrapolation_region$Data_Extrap, Include == 1)) +
  geom_point(aes(x = Lon, y = Lat))

    ### ---------------------------------------  ###


    ### 2.4.2 prepare settings and covariate & formula ----

      ## create the setting file
settings <-  make_settings(n_x = c(100, 250, 500, 1000, 2000)[4], # number of knots, the more knots the longer it runs
                           purpose = "index2", # to calculate index
                           Region = "user", 
                           # when region is user, then input_grid is needed be specified in the fit_model function as the extrapolation we generated earlier
                           # FieldConfig = c("Omega1" = 1, "Epsilon1" = 1, "Omega2" = 0, "Epsilon2" = 0), 
                           FieldConfig = c("Omega1" = 1, "Epsilon1" = 1, "Omega2" = 1, "Epsilon2" = 1),
                           # specify various options for turning on/off spatial and spatial-temporal variation factors in the two linear predictors
                           # omega: spatial variation ; epsilon: spatiotemporal variation
                           # details see here: https://rdrr.io/github/James-Thorson/VAST/man/make_data.html
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

    ### ---------------------------------------  ###







    ### 2.4.3 loop over multiple error/link options to identify the best model ----

# only error dist = c(1,2,4,9) have model fitted
scenario_df <- expand.grid(error.dist = c(1,2,4,9), link.func = c(0:1), AIC = NA, RMSE = NA, finished  = NA) ##details see "make_data" description

# c(9,0) with env.variable is the optimal model, but env messed up the indices, so we ditched the env.var
# error - delta-generalized gamma
# link - Conventional delta-model using logit-link for encounter probability and log-link for positive catch rates

i = 4 # optimal model

for (i in 1:nrow(scenario_df)) {
  scenario_df$finished[i] = "yes"
  write.csv(scenario_df, "results/indices for assessment/summer flounder/mitigation_3_model-based indices/VAST run/scenario.csv", row.names = FALSE)
  
  # change the observation model setting, 
  settings$ObsModel <- c(scenario_df$error.dist[i], scenario_df$link.func[i])
  
  # run model with tryCatch to skip errors when scenario does not converge
  
  VAST_temp_spring <- tryCatch(
    # attempt to fit the model
    {
      VAST_temp_spring <- fit_model(settings = settings,
                                    input_grid = extrapolation_region$Data_Extrap, 
                                    Lat_i = spring_VAST_df$Lat,
                                    Lon_i = spring_VAST_df$Lon,
                                    t_i = as.numeric(spring_VAST_df$Year), 
                                    b_i = as_units(spring_VAST_df$NUMBER, "count"),
                                    a_i = as_units(spring_VAST_df$AreaSwept_km2, "km^2"),
                                    working_dir = "results/indices for assessment/summer flounder/mitigation_3_model-based indices/VAST run")
      VAST_temp_spring # return the fit model if successful
    },
    
    # return NULL if error
    error = function(e) { NULL }
  )
  
  # skip to the next iteration if NULL is returned
  if (is.null(VAST_temp_spring)) { next }
  if (length(VAST_temp_spring$Report) == 1) {
    if (VAST_temp_spring$Report == "Model is not converged") { next }
  }
  
  scenario_df$AIC[i] = VAST_temp_spring$parameter_estimates$AIC
  scenario_df$finished[i] = "yes"
  
  save(VAST_temp_spring, file = paste0("results/indices for assessment/summer flounder/mitigation_3_model-based indices/VAST run/VAST_spring_",
                                       scenario_df$error.dist[i], ".",scenario_df$link.func[i], "_model.Rdata"))
  
  write.csv(scenario_df, "results/indices for assessment/summer flounder/mitigation_3_model-based indices/VAST run/scenario.csv", row.names = FALSE)
  
}





# ---------------------------------------------------------------------------  #




# 3. extract abundance indices ----


  ## 3.1 ALB fall model ----

load("results/indices for assessment/summer flounder/mitigation_3_model-based indices/VAST run/final model/VAST_fall_ALB_model.Rdata")

VAST_fall_data <- plot(VAST_temp_fall, working_dir = "results/indices for assessment/summer flounder/mitigation_3_model-based indices/VAST run/final model/figure/")

VAST_temp_fall_df <- data.frame(YEAR = VAST_fall_data$Index$Table$Time,
                                Season = "FALL",
                                MODEL = "VAST",
                                fit = VAST_fall_data$Index$Table$Estimate, 
                                lwr = VAST_fall_data$Index$Table$Estimate - 1.96 * VAST_fall_data$Index$Table$`Std. Error for Estimate`, 
                                upr = VAST_fall_data$Index$Table$Estimate + 1.96 * VAST_fall_data$Index$Table$`Std. Error for Estimate`) 

write.csv(VAST_temp_fall_df, "results/indices for assessment/summer flounder/mitigation_3_model-based indices/VAST run/final model/VAST_fall_ALB_Indices.csv", row.names = FALSE)


## ---------- ##



  ## 3.2 ALB spring model ----

load("results/indices for assessment/summer flounder/mitigation_3_model-based indices/VAST run/final model/VAST_spring_ALB_model.Rdata")

VAST_spring_data <- plot(VAST_temp_spring, working_dir = "results/indices for assessment/summer flounder/mitigation_3_model-based indices/VAST run/final model/figure/")

VAST_temp_spring_df <- data.frame(YEAR = VAST_spring_data$Index$Table$Time,
                                Season = "SPRING",
                                MODEL = "VAST",
                                fit = VAST_spring_data$Index$Table$Estimate, 
                                lwr = VAST_spring_data$Index$Table$Estimate - 1.96 * VAST_spring_data$Index$Table$`Std. Error for Estimate`, 
                                upr = VAST_spring_data$Index$Table$Estimate + 1.96 * VAST_spring_data$Index$Table$`Std. Error for Estimate`) 

write.csv(VAST_temp_spring_df, "results/indices for assessment/summer flounder/mitigation_3_model-based indices/VAST run/final model/VAST_spring_ALB_Indices.csv", row.names = FALSE)

  ## ---------- ##



  ## 3.3 BIG fall model ----
load("results/indices for assessment/summer flounder/mitigation_3_model-based indices/VAST run/final model/VAST_fall_BIG_model.Rdata")

VAST_fall_data <- plot(VAST_temp_fall, working_dir = "results/indices for assessment/summer flounder/mitigation_3_model-based indices/VAST run/final model/figure/")

VAST_temp_fall_df <- data.frame(YEAR = VAST_fall_data$Index$Table$Time,
                                Season = "FALL",
                                MODEL = "VAST",
                                fit = VAST_fall_data$Index$Table$Estimate, 
                                lwr = VAST_fall_data$Index$Table$Estimate - 1.96 * VAST_fall_data$Index$Table$`Std. Error for Estimate`, 
                                upr = VAST_fall_data$Index$Table$Estimate + 1.96 * VAST_fall_data$Index$Table$`Std. Error for Estimate`) 

write.csv(VAST_temp_fall_df, "results/indices for assessment/summer flounder/mitigation_3_model-based indices/VAST run/final model/VAST_fall_BIG_Indices.csv", row.names = FALSE)

  ## ---------- ##




  ## 3.4 BIG spring model ----

load("results/indices for assessment/summer flounder/mitigation_3_model-based indices/VAST run/final model/VAST_spring_BIG_model.Rdata")

VAST_spring_data <- plot(VAST_temp_spring, working_dir = "results/indices for assessment/summer flounder/mitigation_3_model-based indices/VAST run/final model/figure/")

VAST_temp_spring_df <- data.frame(YEAR = VAST_spring_data$Index$Table$Time,
                                  Season = "SPRING",
                                  MODEL = "VAST",
                                  fit = VAST_spring_data$Index$Table$Estimate, 
                                  lwr = VAST_spring_data$Index$Table$Estimate - 1.96 * VAST_spring_data$Index$Table$`Std. Error for Estimate`, 
                                  upr = VAST_spring_data$Index$Table$Estimate + 1.96 * VAST_spring_data$Index$Table$`Std. Error for Estimate`) 

write.csv(VAST_temp_spring_df, "results/indices for assessment/summer flounder/mitigation_3_model-based indices/VAST run/final model/VAST_spring_BIG_Indices.csv", row.names = FALSE)

  ## ---------- ##


# ---------------------------------------------------------------------------  #





# 4. reformat to the area based on density ----


  ## 4.1 ALB fall ----
  
ALB_fall_df <- read.csv("results/indices for assessment/summer flounder/mitigation_3_model-based indices/VAST run/final model/VAST_fall_ALB_Indices.csv") %>%
  select(-c(MODEL, Season)) %>%
  filter(YEAR <= 2008) 


# extract the total area size from VAST to estimate density

load("results/indices for assessment/summer flounder/mitigation_3_model-based indices/VAST run/final model/VAST_fall_ALB_model.Rdata")

VAST_extra_area <- as.numeric(sum(VAST_temp_fall$extrapolation_list$Area_km2_x)) * 10^6

ALB_fall_df <- ALB_fall_df  %>%
  mutate(fit = fit/VAST_extra_area * 24000, # multiplied with the area covered by each tow - 24000 m2 
         lwr = lwr/VAST_extra_area * 24000,
         upr = upr/VAST_extra_area * 24000)

write.csv(ALB_fall_df, "results/indices for assessment/summer flounder/mitigation_3_model-based indices/ALB_fall_Indices.csv", row.names = FALSE)


    ## -------------------------- ##
    


  ## 4.2 ALB spring ----

ALB_spring_df <- read.csv("results/indices for assessment/summer flounder/mitigation_3_model-based indices/VAST run/final model/VAST_spring_ALB_Indices.csv") %>%
  select(-c(MODEL, Season)) %>%
  filter(YEAR <= 2008) 


# extract the total area size from VAST to estimate density

load("results/indices for assessment/summer flounder/mitigation_3_model-based indices/VAST run/final model/VAST_spring_ALB_model.Rdata")

VAST_extra_area <- as.numeric(sum(VAST_temp_spring$extrapolation_list$Area_km2_x)) * 10^6

ALB_spring_df <- ALB_spring_df  %>%
  mutate(fit = fit/VAST_extra_area * 24000, # multiplied with the area covered by each tow - 24000 m2 
         lwr = lwr/VAST_extra_area * 24000,
         upr = upr/VAST_extra_area * 24000)

write.csv(ALB_spring_df, "results/indices for assessment/summer flounder/mitigation_3_model-based indices/ALB_spring_Indices.csv", row.names = FALSE)

  ## -------------------------- ##


  ## 4.3 BIG fall ----

BIG_fall_df <- read.csv("results/indices for assessment/summer flounder/mitigation_3_model-based indices/VAST run/final model/VAST_fall_BIG_Indices.csv") %>%
  select(-c(MODEL, Season))


# extract the total area size from VAST to estimate density

load("results/indices for assessment/summer flounder/mitigation_3_model-based indices/VAST run/final model/VAST_fall_BIG_model.Rdata")

VAST_extra_area <- as.numeric(sum(VAST_temp_fall$extrapolation_list$Area_km2_x)) * 10^6

BIG_fall_df <- BIG_fall_df  %>%
  mutate(fit = fit/VAST_extra_area * 24000, # multiplied with the area covered by each tow - 24000 m2 
         lwr = lwr/VAST_extra_area * 24000,
         upr = upr/VAST_extra_area * 24000)

write.csv(BIG_fall_df, "results/indices for assessment/summer flounder/mitigation_3_model-based indices/BIG_fall_Indices.csv", row.names = FALSE)


  ## -------------------------- ##


  ## 4.4 BIG spring ----

BIG_spring_df <- read.csv("results/indices for assessment/summer flounder/mitigation_3_model-based indices/VAST run/final model/VAST_spring_BIG_Indices.csv") %>%
  select(-c(MODEL, Season)) 


# extract the total area size from VAST to estimate density

load("results/indices for assessment/summer flounder/mitigation_3_model-based indices/VAST run/final model/VAST_spring_BIG_model.Rdata")

VAST_extra_area <- as.numeric(sum(VAST_temp_spring$extrapolation_list$Area_km2_x)) * 10^6

BIG_spring_df <- BIG_spring_df  %>%
  mutate(fit = fit/VAST_extra_area * 24000, # multiplied with the area covered by each tow - 24000 m2 
         lwr = lwr/VAST_extra_area * 24000,
         upr = upr/VAST_extra_area * 24000)

write.csv(BIG_spring_df, "results/indices for assessment/summer flounder/mitigation_3_model-based indices/BIG_spring_Indices.csv", row.names = FALSE)

  ## -------------------------- ##



# 6. apply VAST indices to the ASAP input ----

# for the indices it requires several steps
# 1st need to get density, for design-based indices it is stratified mean
# 2nd, for bigelow indices after 2009, there is a SWAN adjustment as below

# SAM: The assessment model splits the Bigelow and Albatross survey indices. 
  # In preparation for application to the assessment model the Bigelow index most notably is scaled up for by-tow swept area, 
  # so that's the major difference that I assume you're noticing. It would be difficult to incorporate those modifications since they only cover 
  # the Bigelow years and I don't think those types of modifications would change your overall conclusions anyhow. 
  # There is some information in SAW 66 (66th Northeast Regional Stock Assessment Workshop (66th SAW) Assessment Report). 
  #  If you search for "SWAN" (swept area numbers) you will find the relevant text.


  # ASAP original inputs

data <- ReadASAP3DatFile("assessment model/summer flounder/ASAP input data/ASAP3_MTA2023_FINAL.DAT")
data$survey.names

# 2 is BTS ALB spring, 3 is BTS ALB fall, they are 1982-2008
# 25 is BTS BIG spring, 26 is BTS BIG fall, they are 2009-2022
# col.1 is year 1982-2022, 1982-2008 correspond to row 1-27, 2009-2022 correspond to row 28-41, 
# col.2 is sum value, col.3 is CV
# column 4-11 correspond to age 1-8,  col. 11 is sample size

ALB_spring <- data$dat$IAA_mats[[2]][1:27,c(2,4:11)]
ALB_fall <- data$dat$IAA_mats[[3]][1:27,c(2,4:11)]

BIG_spring <- data$dat$IAA_mats[[25]][28:41,c(2,4:11)]
BIG_fall <- data$dat$IAA_mats[[26]][28:41,c(2,4:11)]


  # area size

load("results/indices for assessment/summer flounder/mitigation_3_model-based indices/VAST run/final model/VAST_fall_model.Rdata")
fall_area <- sum(VAST_temp_fall$extrapolation_list$Area_km2_x)
remove(VAST_temp_fall)

load("results/indices for assessment/summer flounder/mitigation_3_model-based indices/VAST run/final model/VAST_spring_model.Rdata")
spring_area <- sum(VAST_temp_spring$extrapolation_list$Area_km2_x)
remove(VAST_temp_spring)
  

  # VAST indices

VAST_fall_df <- read.csv("results/indices for assessment/summer flounder/mitigation_3_model-based indices/VAST run/final model/VAST_fall_Indices.csv")
VAST_spring_df <- read.csv("results/indices for assessment/summer flounder/mitigation_3_model-based indices/VAST run/final model/VAST_spring_Indices.csv")


## 6.1 ALB fall  ----

# 1982 - 2008

ALB_fall

VAST_fall_df$fit/151225.1


312,000 

## --------------------------------------- ##



## 6.2 ALB spring ----

# 1982 - 2008

## --------------------------------------- ##



## 6.3 BIG fall ----

# 2009 - 2022

## --------------------------------------- ##



## 6.4 BIG spring ----

# 2009 - 2022

## --------------------------------------- ##




# ---------------------------------------------------------------------------  #





# 25 is# 2 is BTS ALB spring, 3 is BTS ALB fall, they are 1982-2008
# 25 is BTS BIG spring, 26 is BTS BIG fall, they are 2009-2022
# col.1 is year 1982-2022, 1982-2008 correspond to row 1-27, 2009-2022 correspond to row 28-41, 
# col.2 is sum value, col.3 is CV
# column 4-11 correspond to age 1-8,  col. 11 is sample size