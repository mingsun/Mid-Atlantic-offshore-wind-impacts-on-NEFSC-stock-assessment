library(tidyverse)
library(VAST)
library(splines)
library(Metrics)

# this mitigation scenario will use model-based indices to replace the original indices
# since there are different indices used, new models will need to be fitted here
# note the these model-based indices will only use WEE dataset
# according to the first paper, the best model is VAST for all scenario



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

# ---------------------------------------------------------------------------  #




# 2. VAST model prep ----

full_VAST_df <- full_df %>%   # convert YEAR into a factor
  filter(OWF == "OUTSIDE") 



# ---------------------------------------------------------------------------  #



# 3. Fall model ----
fall_VAST_df <- subset(full_VAST_df, SEASON == "FALL")

  ## 3.1 prepare extrapolation area ----

extrapolation_region <- make_extrapolation_info(Region = "other",
                                                strata.limits = data.frame(STRATA = "All_areas"),
                                                observations_LL = fall_VAST_df[,c('Lat','Lon')],
                                                maximum_distance_from_sample = 15)

# write.csv(extrapolation_region$Data_Extrap, file = "results/indices for assessment/squid/mitigation_2_model-based indices/VAST run/extra_list.csv", row.names = FALSE)

ggplot(subset(extrapolation_region$Data_Extrap, Include == 1)) +
  geom_point(aes(x = Lon, y = Lat))


   ## 3.2 prepare settings and covariate & formula ----

#     ### create the setting file
# settings <-  make_settings(n_x = c(100, 250, 500, 1000, 2000)[4], # number of knots, the more knots the longer it runs
#                            purpose = "index2", # to calculate index
#                            Region = "user",
#                            # when region is user, then input_grid is needed be specified in the fit_model function as the extrapolation we generated earlier
#                            # the setting for temporal modelign is: FieldConfig = c("Omega1" = 0, "Epsilon1" = 0, "Omega2" = 1, "Epsilon2" = 1), 
#                            FieldConfig = c("Omega1" = 1, "Epsilon1" = 1, "Omega2" = 1, "Epsilon2" = 1),
#                            # specify various options for turning on/off spatial and spatial-temporal variation factors in the two linear predictors
#                            # omega: spatial variation ; epsilon: spatiotemporal variation
#                            # details see here: https://rdrr.io/github/James-Thorson/VAST/man/make_data.html
#                            # !!! turned off first predictor for longfin squid
#                            RhoConfig = c("Beta1" = 0, "Beta2" = 0, "Epsilon1" = 0, "Epsilon2" = 0), 
#                            # specify whether either intercepts (Beta1 and Beta2) or spatio-temporal variation (Epsilon1 and Epsilon2) is structured among time intervals
#                            # untouched default setting: all is zero, so all are fixed effect
#                            OverdispersionConfig= c(0,0),
#                            # vessel/targeting effects for encounter probability and positive catch rate
#                            # untouched default setting (0,0) means no effects
#                            ObsModel = c(2,1),
#                            # link functions for encounter probability and positive catch rate are available
#                            # untouched default setting here
#                            bias.correct = TRUE # to save time set as false, epsilon bias-correction
# )

  ## 3.3 loop over multiple error/link options to identify the best model ----


# setting for spatiotemporal model: scenario_df <- expand.grid(error.dist = c(2,4,9), link.func = c(3,4), AIC = NA, RMSE = NA) ##details see "make_data" description
scenario_df <- expand.grid(error.dist = c(1,2,4,9), link.func = c(0,1), AIC = NA, RMSE = NA) ##details see "make_data" description

# c(2,0) with env.variable is the optimal model, but env messed up the indices, so we ditched the env.var
# error - delta-gamma
# link - Conventional delta-model using logit-link for encounter probability and log-link for positive catch rates 

i = 2 # optimal model

for (i in 5:nrow(scenario_df)) {
  
  # change the observation model setting, 
  settings$ObsModel <- c(scenario_df$error.dist[i], scenario_df$link.func[i])
  # settings$FieldConfig = c("Omega1" = 1, "Epsilon1" = 1, "Omega2" = 1, "Epsilon2" = 1) # turn off sth as suggested by a failed run
  
  # run model with tryCatch to skip errors when scenario does not converge
  
  VAST_temp_fall <- tryCatch(
    # attempt to fit the model
    {
      VAST_temp_fall <- fit_model(settings = settings,
                                  input_grid = extrapolation_region$Data_Extrap,
                                  Lat_i = fall_VAST_df$Lat,
                                  Lon_i = fall_VAST_df$Lon,
                                  t_i = as.numeric(fall_VAST_df$Year), 
                                  b_i = fall_VAST_df$CATCH_WT_CAL,
                                  a_i = as_units(fall_VAST_df$AreaSwept_km2, "km^2"),
                                  # max_cells = 10000, # avoid automatic grid reduction 
                                  # covariate_data = full_VAST_df, # use full_VAST_df here to avoid missing year (2020) for fall
                                  # X1_formula = p1_formula,
                                  # X2_formula = p2_formula,
                                  working_dir = "results/indices for assessment/squid/mitigation_3_model-based indices/VAST run")
      VAST_temp_fall # return the fit model if successful
    },
    
    # return NULL if error
    error = function(e) { NULL }
  )
  
  # skip to the next iteration if NULL is returned
  if (is.null(VAST_temp_fall)) { next }
  if (length(VAST_temp_fall$Report) == 1) {
    if (VAST_temp_fall$Report != "Model is not converged") { next }
  }
  
  scenario_df$AIC[i] = VAST_temp_fall$parameter_estimates$AIC
  
  save(VAST_temp_fall, file = paste0("results/indices for assessment/squid/mitigation_3_model-based indices/VAST run/VAST_fall_",
                                     scenario_df$error.dist[i], ".",scenario_df$link.func[i], "_model.Rdata"))
  
  # save(VAST_temp_fall, full_VAST_df, file = paste0("results/temporal model/longfin.squid/models objects/VAST_temp_fall.Rdata"))
  # load(file = paste0("results/temporal model/longfin.squid/models objects/VAST_temp_fall.Rdata"))
  
  write.csv(scenario_df, "results/indices for assessment/squid/mitigation_3_model-based indices/VAST run/scenario.csv", row.names = FALSE)
  
}


# ---------------------------------------------------------------------------  #


# 4. spring model ----

spring_VAST_df <- subset(full_VAST_df, SEASON == "SPRING")

  ## 4.1 prepare extrapolation area ----

# the other one generates dots all over the map so ditched it
extrapolation_region <- make_extrapolation_info(Region = "northwest_atlantic",  # use the embedded dataset
                                                strata.limits = data.frame(STRATA = "EPU"),
                                                epu_to_use = c("All", "Georges_Bank", "Mid_Atlantic_Bight", "Scotian_Shelf", "Gulf_of_Maine", "Other")[c(2,3)],
                                                # if use "Georges_Bank", "Mid_Atlantic_Bight", there will be 436 rows missing for summer flounder
                                                # all outside OWF, not affecting the gap analysis
                                                DirPath = "results/spatial model/longfin.squid/models objects/VAST/")
colnames(extrapolation_region$Data_Extrap)[4] <- "Area_km2"

# extrapolation_region <- make_extrapolation_info(Region = "other",
#                                                 strata.limits = data.frame(STRATA = "All_areas"),
#                                                 observations_LL = spring_VAST_df[,c('Lat','Lon')],
#                                                 maximum_distance_from_sample = 15)


ggplot(subset(extrapolation_region$Data_Extrap, Include == 1)) +
  geom_point(aes(x = Lon, y = Lat))


  ## 4.2 prepare settings and covariate & formula ----

    ### create the setting file
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



  ## 4.3 loop over multiple error/link options to identify the best model ----

# setting for spatiotemporal model: scenario_df <- expand.grid(error.dist = c(2,4,9), link.func = c(3,4), AIC = NA, RMSE = NA) ##details see "make_data" description
scenario_df <- expand.grid(error.dist = c(1,2,4,9), link.func = c(0,1), AIC = NA, RMSE = NA) ##details see "make_data" description

# c(9,0) with env.variable is the optimal model, but env messed up the indices, so we ditched the env.var
# error - Delta-Generalized Gamma
# link - Conventional delta-model using logit-link for encounter probability and log-link for positive catch rates 

i = 4 # optimal model

for (i in 5:nrow(scenario_df)) {
  
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
                                    b_i = spring_VAST_df$CATCH_WT_CAL,
                                    a_i = as_units(spring_VAST_df$AreaSwept_km2, "km^2"),
                                    # covariate_data = full_VAST_df, # use full_VAST_df here to avoid missing year (2020) for fall
                                    # X1_formula = p1_formula,
                                    # X2_formula = p2_formula,
                                    working_dir = "results/indices for assessment/squid/mitigation_3_model-based indices/VAST run")
      VAST_temp_spring # return the fit model if successful
    },
    
    # return NULL if error
    error = function(e) { NULL }
  )
  
  # skip to the next iteration if NULL is returned
  if (is.null(VAST_temp_spring)) { next }
  if (length(VAST_temp_spring$Report) == 1) {
    if (VAST_temp_spring$Report != "Model is not converged") { next }
  }
  
  scenario_df$AIC[i] = VAST_temp_spring$parameter_estimates$AIC
  
  save(VAST_temp_spring, file = paste0("results/indices for assessment/squid/mitigation_3_model-based indices/VAST run/VAST_spring_",
                                       scenario_df$error.dist[i], ".",scenario_df$link.func[i], "_model.Rdata"))
  
  # save(VAST_temp_spring, file = paste0("results/temporal model/longfin.squid/models objects/VAST_temp_spring.Rdata"))
  
  # load(file = paste0("results/temporal model/longfin.squid/models objects/VAST_temp_spring.Rdata"))
  
  write.csv(scenario_df, "results/indices for assessment/squid/mitigation_3_model-based indices/VAST run/scenario.csv", row.names = FALSE)
  
}


# ---------------------------------------------------------------------------  #



# 5. extract abundance indices ----


  ## 5.1 fall model ----

load("results/indices for assessment/squid/mitigation_3_model-based indices/VAST run/final model/VAST_fall_model.Rdata")


VAST_fall_data <- plot(VAST_temp_fall, working_dir = "results/indices for assessment/squid/mitigation_2_model-based indices/VAST run/final model/figure/")

VAST_temp_fall_df <- data.frame(YEAR = VAST_fall_data$Index$Table$Time,
                                Season = "FALL",
                                MODEL = "VAST",
                                fit = VAST_fall_data$Index$Table$Estimate, 
                                lwr = VAST_fall_data$Index$Table$Estimate - 1.96 * VAST_fall_data$Index$Table$`Std. Error for Estimate`, 
                                upr = VAST_fall_data$Index$Table$Estimate + 1.96 * VAST_fall_data$Index$Table$`Std. Error for Estimate`) 

write.csv(VAST_temp_fall_df, "results/indices for assessment/squid/mitigation_3_model-based indices/VAST run/final model/VAST_fall_Indices.csv", row.names = FALSE)



  ## 5.2 spring model ----

load("results/indices for assessment/squid/mitigation_23_model-based indices/VAST run/final model/VAST_spring_model.Rdata")

VAST_spring_data <- plot(VAST_temp_spring, working_dir = "results/indices for assessment/squid/mitigation_2_model-based indices/VAST run/final model/figure/")

VAST_temp_spring_df <- data.frame(YEAR = VAST_spring_data$Index$Table$Time,
                                  Season = "SPRING",
                                  MODEL = "VAST",
                                  fit = VAST_spring_data$Index$Table$Estimate, 
                                  lwr = VAST_spring_data$Index$Table$Estimate - 1.96 * VAST_spring_data$Index$Table$`Std. Error for Estimate`, 
                                  upr = VAST_spring_data$Index$Table$Estimate + 1.96 * VAST_spring_data$Index$Table$`Std. Error for Estimate`) 

write.csv(VAST_temp_spring_df, "results/indices for assessment/squid/mitigation_3_model-based indices/VAST run/final model/VAST_spring_Indices.csv", row.names = FALSE)


# ---------------------------------------------------------------------------  #



# 6. apply VAST indices to the biomass index method ----

  # use the total biomass as the input to the spreadsheet 
  # E:\Stony Brook job\NYSERDA Offshore wind\assessment model\squid\MIT.2 model





