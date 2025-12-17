library(tidyverse)
library(randomForest)



# this mitigation scenario will use model to standardized tow level abundance to replace the original tow abundance
# note that all tows (inside and outside OWF) are standardized using full dataset
# according to the first paper, the best model is Random Forest for all scenario


# 1. generate full dataset with env.var ----

catch_by_tow_df <- read.csv("results/stratified.mean.indices/surfclam/total.catch.by.tow.csv") %>%
  select(ID, REGION, WEIGHT)

full_df <- read.csv("results/stratified.mean.indices/surfclam/tow.list.csv") %>%
  mutate(ID.temp = paste(CRUISE6, STATION, sep = ".")) %>%
  left_join(catch_by_tow_df) 

AS_QQ_overlay_df <- read.csv("results/AS_OQ_tows_overlay_final.csv") %>%
  mutate(ID.temp = paste(CRUISE6, STATION, sep = "."))

full_df <- full_df %>% 
  mutate(OWF = ifelse(ID.temp %in% unique(AS_QQ_overlay_df$ID.temp), "INSIDE", "OUTSIDE")) %>% 
  select(-c(ID.temp)) %>%
  rename(BIOMASS = KGPERTOW)  %>%
  filter(REGION == "SVAtoSNE")

full_df[is.na(full_df$LON),]$LON <-  mean(c(-74.3187, -73.7816)) # manually fix the only station without LON

full_RF_df <- full_df %>%   # convert YEAR into a factor
  mutate(Year = as.factor(YEAR)) %>%
  filter(!if_any(c(NPERTOW), is.na))

remove(catch_by_tow_df, AS_QQ_overlay_df)

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



# 3. fit model ---------------------------------------------------------------------------------------------



  ## 3.0 checking temp data availability by year 

ggplot(full_RF_df) +
  geom_point(aes(x = LON, y = LAT)) +
  facet_wrap(.~YEAR)


temp_data_avail_df <- full_RF_df %>%
  mutate(MONTH = substr(CRUISE6, 5, 6)) %>%
  summarize(N.TOW.MISSING = sum(is.na(TEMP)), 
            N.TOTAL = length(TEMP), 
            PROP.NO = N.TOW.MISSING/N.TOTAL,.by = c(YEAR, MONTH))


# save it and manually edti the INPUT.YEAR and then read back
# write.csv(temp_data_avail_df, "results/indices for assessment/surfclam/mitigation_2_model-standardized/temp.data.avail.csv", row.names = FALSE)
temp_data_avail_df <- read.csv("results/indices for assessment/surfclam/mitigation_2_model-standardized/temp.data.avail.csv")




## 3.1 Kriging the missing BOTTEM by Year ----

library(gstat)

no_TEMP_df <- data.frame() # an empty df to store all interpolated tows

for (y in unique(full_RF_df$Year)) {
  
  annual_df <- subset(full_RF_df, Year == y)
  INPUT.YEAR <- subset(temp_data_avail_df, YEAR == y)$INPUT.YEAR
  
  if (sum(is.na(annual_df$TEMP)) == 0) next
  
  # get the df to develop Kriging
  with_TEMP_df <- subset(subset(full_RF_df, Year == INPUT.YEAR), !is.na(TEMP))
  sp::coordinates(with_TEMP_df) <-  ~LON + LAT
  
  # develop a variogram model
  temp.vgm <- variogram(TEMP ~ LON + LAT, data = with_TEMP_df, width = 0.01)
  # temp.fit <- fit.variogram(temp.vgm, vgm(c( "Nug")))
  temp.fit <- fit.variogram(temp.vgm, vgm(c("Nug", "Sph","Mat", "Pen" , "Ste", "Cir" , "Bes", "Exc")))
  

  # create the spatial scale to be interpolated
  missing_TEMP_coord <- subset(annual_df, is.na(TEMP))[, c("LAT", "LON")]
  sp::coordinates(missing_TEMP_coord) <-  ~LON + LAT
  
  # perform kriging interpolation for the tows without BOTTEMP
  temp.kriged <- krige(TEMP ~ LON + LAT, with_TEMP_df, missing_TEMP_coord, model = temp.fit)
  
  # store the interpolated value
  missing_TEMP_df <- subset(annual_df, is.na(TEMP))
  missing_TEMP_df$TEMP <- temp.kriged$var1.pred
  
  no_TEMP_df <- rbind(no_TEMP_df, missing_TEMP_df)
}

remove(y, annual_df, INPUT.YEAR, with_TEMP_df, temp.vgm, temp.fit, missing_TEMP_coord, temp.kriged, missing_TEMP_df)


# generate the final df for RF
full_RF_df <- rbind(subset(full_RF_df, !is.na(TEMP)),
                    no_TEMP_df) %>%
              arrange(ID)




    ## 3.2 run Random Forest ----

set.seed(2) # ensure reproducibility
RF_full <- randomForest(NPERTOW ~ Year + LAT + TEMP + DEPTH, data = subset(full_RF_df, OWF == "OUTSIDE"), mtry = 2, ntree = 1000)


    ## 3.3  predict the results and save ----
full_RF_df <- full_RF_df %>%
  mutate(PREDICTED_NPERTOW =  predict(RF_full, newdata = ., type = "response")) %>%
  mutate(Final_NPERTOW = case_when(OWF == "INSIDE" ~ PREDICTED_NPERTOW,
                                        OWF == "OUTSIDE" ~  NPERTOW))

write.csv(full_RF_df, "results/indices for assessment/surfclam/mitigation_2_model-standardized/Full_RF_standardized_tow.csv", row.names = FALSE)


rm(list = ls())



# -------------------------------------------------------------------------------------------------------- #


# 4. Design-based abundance indices ----

RF_df <- read.csv("results/indices for assessment/surfclam/mitigation_2_model-standardized/Full_RF_standardized_tow.csv") %>%
  select(YEAR, STRATUM, TOW, ID, TOTAL_N_STATION, WEIGHT, Final_NPERTOW)



  ## 4.1 mean numbers by stratum ---------------------------------------------------------------------------------------------

mean_N_stratum_df <- RF_df %>%
  group_by(YEAR, STRATUM) %>%
  mutate(MEAN_N_STRATUM = mean(Final_NPERTOW)) %>% # mean within a strata
  mutate(VAR_STRATUM = var(Final_NPERTOW)) %>% 
  ungroup() %>% 
  select(c(YEAR, STRATUM, TOTAL_N_STATION, WEIGHT, MEAN_N_STRATUM, VAR_STRATUM)) %>%
  distinct() %>% # downsize the data frame to a minimal without repetitive rows
  arrange(YEAR)


  ## 4.2 stratified mean numbers ---------------------------------------------------------------------------------------------

stratified_mean_N_df <- mean_N_stratum_df %>%
  group_by(YEAR) %>%
  summarize(STRATIFIED_MEAN_N = weighted.mean(MEAN_N_STRATUM, w = WEIGHT), # stratified mean
            STRATIFIED_VAR = sum(WEIGHT^2 * VAR_STRATUM / TOTAL_N_STATION, na.rm = TRUE), # variance
            STRATIFIED_SE = sqrt(STRATIFIED_VAR)) %>%   # standard deviance
  mutate(up_CI_95 = qnorm(0.975, mean = STRATIFIED_MEAN_N, sd = STRATIFIED_SE),
         lo_CI_95 = qnorm(0.025, mean = STRATIFIED_MEAN_N, sd = STRATIFIED_SE),
         lo_CI_95 = ifelse(lo_CI_95 < 0, yes = 0, no = lo_CI_95),
         CV = STRATIFIED_SE / STRATIFIED_MEAN_N) %>%
  ungroup()


write.csv(stratified_mean_N_df, "results/indices for assessment/surfclam/mitigation_2_model-standardized/MIT.2.stratified.mean.indices.csv", row.names = FALSE)


# -------------------------------------------------------------------------------------------------------- #


