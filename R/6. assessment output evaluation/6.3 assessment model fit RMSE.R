library(tidyverse)
library(ASAPplots)
library(r4ss)


# RMSE for input data: total catch ,impacted surveys, and total surveys


# 1. summer flounder ----

  ## 1.1 base model ----

asap.name <- "ASAP3_MTA2023_FINAL"
asap <- dget("assessment model/summer flounder/ASAP/base model/ASAP3_MTA2023_FINAL.rdat") 
od <- "results/stock assessment/summer flounder/model fit/"


PlotRMSEtable(asap.name, asap, save.plots = FALSE, od, plotf = "jpg")

read.csv("results/stock assessment/summer flounder/model fit/RMSE.Table.ASAP3_MTA2023_FINAL.csv") %>%
  rename("Quantity" = asap.name) %>%
  add_column(Scenario = "Base Model") %>%
  write.csv("results/stock assessment/summer flounder/model fit/RMSE_base.csv", row.names = FALSE)


# 25 is BTS BIG spring, 26 is BTS BIG fall, they are 2009-2022
RMSE_base_df <- read.csv("results/stock assessment/summer flounder/model fit/RMSE_base.csv") %>%
  filter(Quantity %in% c("catch.tot", "ind2", "ind3" , "ind25", "ind26", "ind.total"))

remove(asap.name, asap, od)




  ## 1.2 WEE_design-based_indices_model assessment ----

asap.name <- "ASAP3_MTA2023_WEE_BTS"
asap <- dget("assessment model/summer flounder/ASAP/WEE_design-based_indices_model/ASAP3_MTA2023_WEE_BTS.RDAT") 
od <- "results/stock assessment/summer flounder/model fit/"


PlotRMSEtable(asap.name, asap, save.plots = FALSE, od, plotf = "jpg")

read.csv("results/stock assessment/summer flounder/model fit/RMSE.Table.ASAP3_MTA2023_WEE_BTS.csv") %>%
  rename("Quantity" = asap.name) %>%
  add_column(Scenario = "WEE Model") %>%
  write.csv("results/stock assessment/summer flounder/model fit/RMSE_WEE_BTS.csv", row.names = FALSE)


# 25 is BTS BIG spring, 26 is BTS BIG fall, they are 2009-2022
RMSE_WEE_df <- read.csv("results/stock assessment/summer flounder/model fit/RMSE_WEE_BTS.csv") %>%
  filter(Quantity %in% c("catch.tot", "ind02", "ind03" , "ind25", "ind26", "ind.total"))

remove(asap.name, asap, od)



    ## 1.3 MIT.1 model ----

asap.name <- "ASAP3_MTA2023_MIT.1_BTS"
asap <- dget("assessment model/summer flounder/ASAP/MIT.1 model/ASAP3_MTA2023_MIT.1_BTS.RDAT") 
od <- "results/stock assessment/summer flounder/model fit/"


PlotRMSEtable(asap.name, asap, save.plots = FALSE, od, plotf = "jpg")

read.csv("results/stock assessment/summer flounder/model fit/RMSE.Table.ASAP3_MTA2023_MIT.1_BTS.csv") %>%
  rename("Quantity" = asap.name) %>%
  add_column(Scenario = "MIT.1 Model") %>%
  write.csv("results/stock assessment/summer flounder/model fit/RMSE_MIT.1.csv", row.names = FALSE)


# 25 is BTS BIG spring, 26 is BTS BIG fall, they are 2009-2022
RMSE_MIT.1_df <- read.csv("results/stock assessment/summer flounder/model fit/RMSE_MIT.1.csv") %>%
  filter(Quantity %in% c("catch.tot", "ind2", "ind3" , "ind25", "ind26", "ind.total"))

remove(asap.name, asap, od)




    ## 1.6 combine and save ----

RMSE_combine_df <- rbind(RMSE_base_df, RMSE_WEE_df, RMSE_MIT.1_df) %>%
  add_column(Species = "Summer Flounder") %>%
  write.csv("results/stock assessment/summer flounder/model fit/RMSE_combined.csv", row.names = FALSE)





# ------------------------------------------------------------------------ #




# 2. squid ----

# For squid with index-based method, we cannot do RMSE
# we calculate the final indices cv


base_df <- read.csv("results/stock assessment/squid/index-based summary_base.csv")
WEE_df <-  read.csv("results/stock assessment/squid/index-based summary_WEE.csv")
MIT.1_df <- read.csv("results/stock assessment/squid/index-based summary_MIT.1.csv")

index_df <- rbind(base_df, WEE_df, MIT.1_df) %>%
  filter(Year >= 1987 & !is.na(Annualized.Exploitation.Indices.final))  %>%
  select(Scenario, Year, Annual.Catch, Annualized.Exploitation.Indices.final) %>%
  # mutate(Ratio = Annual.Catch/Annualized.Exploitation.Indices.final) %>%
  group_by(Scenario) %>%
  summarize(CV = sd(Annualized.Exploitation.Indices.final)/mean(Annualized.Exploitation.Indices.final)) %>%
  add_column(Species = "Longfin Squid") %>%
  write.csv("results/stock assessment/squid/model fit/stability.csv", row.names = FALSE)

remove(index_df)


# ------------------------------------------------------------------------ #





# 3. surfclam ----

  ## 3.1 base model ----

assessment_results <- SS_output("assessment model/surfclam/base model/", verbose = FALSE)


  # RMSE by survey indices
RMSE_ind_survey_df <- assessment_results$cpue %>%
  mutate(residual = log(Obs) - log(Exp)) %>% # log is needed to be consistent with ASAp
  group_by(Fleet_name) %>%
  summarise(N = n(), RMSE = sqrt(mean(residual^2, na.rm = TRUE))) %>%
  ungroup()


  # RMSE for all surveys combined
RMSE_all_survey_df <- assessment_results$cpue %>%
  mutate(residual = log(Obs) - log(Exp)) %>%
  summarise(Fleet_name = "Survey.total", N = n(), RMSE = sqrt(mean(residual^2, na.rm = TRUE))) 
  

RMSE_df <- rbind(RMSE_ind_survey_df, RMSE_all_survey_df) %>%
  add_column(Scenario = "Base Model")

write.csv(RMSE_df, "results/stock assessment/surfclam/model fit/RMSE_base.csv", row.names = FALSE)

remove(assessment_results, RMSE_ind_survey_df, RMSE_all_survey_df, RMSE_df)




  ## 3.2 WEE model ----

assessment_results <- SS_output("assessment model/surfclam/WEE_design-based_indices_model/", verbose = FALSE)


# RMSE by survey indices
RMSE_ind_survey_df <- assessment_results$cpue %>%
  mutate(residual = log(Obs) - log(Exp)) %>% # log is needed to be consistent with ASAp
  group_by(Fleet_name) %>%
  summarise(N = n(), RMSE = sqrt(mean(residual^2, na.rm = TRUE))) %>%
  ungroup()


# RMSE for all surveys combined
RMSE_all_survey_df <- assessment_results$cpue %>%
  mutate(residual = log(Obs) - log(Exp)) %>%
  summarise(Fleet_name = "Survey.total", N = n(), RMSE = sqrt(mean(residual^2, na.rm = TRUE))) 


RMSE_df <- rbind(RMSE_ind_survey_df, RMSE_all_survey_df) %>%
  add_column(Scenario = "WEE Model")

write.csv(RMSE_df, "results/stock assessment/surfclam/model fit/RMSE_WEE.csv", row.names = FALSE)

remove(assessment_results, RMSE_ind_survey_df, RMSE_all_survey_df, RMSE_df)



    ## 3.3 MIT.1 model ----

assessment_results <- SS_output("assessment model/surfclam/MIT.1 model/", verbose = FALSE)


# RMSE by survey indices
RMSE_ind_survey_df <- assessment_results$cpue %>%
  mutate(residual = log(Obs) - log(Exp)) %>% # log is needed to be consistent with ASAp
  group_by(Fleet_name) %>%
  summarise(N = n(), RMSE = sqrt(mean(residual^2, na.rm = TRUE))) %>%
  ungroup()


# RMSE for all surveys combined
RMSE_all_survey_df <- assessment_results$cpue %>%
  mutate(residual = log(Obs) - log(Exp)) %>%
  summarise(Fleet_name = "Survey.total", N = n(), RMSE = sqrt(mean(residual^2, na.rm = TRUE))) 


RMSE_df <- rbind(RMSE_ind_survey_df, RMSE_all_survey_df) %>%
  add_column(Scenario = "MIT.1 Model")

write.csv(RMSE_df, "results/stock assessment/surfclam/model fit/RMSE_MIT.1.csv", row.names = FALSE)

remove(assessment_results, RMSE_ind_survey_df, RMSE_all_survey_df, RMSE_df)



    ## 3.6 combine and save ----

RMSE_base_df <- read.csv("results/stock assessment/surfclam/model fit/RMSE_base.csv")
RMSE_WEE_df <- read.csv("results/stock assessment/surfclam/model fit/RMSE_WEE.csv")
RMSE_MIT.1_df <- read.csv("results/stock assessment/surfclam/model fit/RMSE_MIT.1.csv")


RMSE_combine_df <- rbind(RMSE_base_df, RMSE_WEE_df, RMSE_MIT.1_df) %>%
  add_column(Species = "Surfclam") %>%
  write.csv("results/stock assessment/surfclam/model fit/RMSE_combined.csv", row.names = FALSE)


# ------------------------------------------------------------------------ #







# 4. quahog ----


  ## 4.1 base model ----

assessment_results <- SS_output("assessment model/quahog/base model/", verbose = FALSE)


# RMSE by survey indices
RMSE_ind_survey_df <- assessment_results$cpue %>%
  mutate(residual = log(Obs) - log(Exp)) %>% # log is needed to be consistent with ASAp
  group_by(Fleet_name) %>%
  summarise(N = n(), RMSE = sqrt(mean(residual^2, na.rm = TRUE))) %>%
  ungroup()


# RMSE for all surveys combined
RMSE_all_survey_df <- assessment_results$cpue %>%
  mutate(residual = log(Obs) - log(Exp)) %>%
  summarise(Fleet_name = "Survey.total", N = n(), RMSE = sqrt(mean(residual^2, na.rm = TRUE))) 


RMSE_df <- rbind(RMSE_ind_survey_df, RMSE_all_survey_df) %>%
  add_column(Scenario = "Base Model")

write.csv(RMSE_df, "results/stock assessment/quahog/model fit/RMSE_base.csv", row.names = FALSE)

remove(assessment_results, RMSE_ind_survey_df, RMSE_all_survey_df, RMSE_df)




    ## 4.2 WEE model ----

assessment_results <- SS_output("assessment model/quahog/WEE_design-based_indices_model/", verbose = FALSE)


# RMSE by survey indices
RMSE_ind_survey_df <- assessment_results$cpue %>%
  mutate(residual = log(Obs) - log(Exp)) %>% # log is needed to be consistent with ASAp
  group_by(Fleet_name) %>%
  summarise(N = n(), RMSE = sqrt(mean(residual^2, na.rm = TRUE))) %>%
  ungroup()


# RMSE for all surveys combined
RMSE_all_survey_df <- assessment_results$cpue %>%
  mutate(residual = log(Obs) - log(Exp)) %>%
  summarise(Fleet_name = "Survey.total", N = n(), RMSE = sqrt(mean(residual^2, na.rm = TRUE))) 


RMSE_df <- rbind(RMSE_ind_survey_df, RMSE_all_survey_df) %>%
  add_column(Scenario = "WEE Model")

write.csv(RMSE_df, "results/stock assessment/quahog/model fit/RMSE_WEE.csv", row.names = FALSE)

remove(assessment_results, RMSE_ind_survey_df, RMSE_all_survey_df, RMSE_df)




    ## 3.3 MIT.1 model ----

assessment_results <- SS_output("assessment model/quahog/MIT.1 model/", verbose = FALSE)


# RMSE by survey indices
RMSE_ind_survey_df <- assessment_results$cpue %>%
  mutate(residual = log(Obs) - log(Exp)) %>% # log is needed to be consistent with ASAp
  group_by(Fleet_name) %>%
  summarise(N = n(), RMSE = sqrt(mean(residual^2, na.rm = TRUE))) %>%
  ungroup()


# RMSE for all surveys combined
RMSE_all_survey_df <- assessment_results$cpue %>%
  mutate(residual = log(Obs) - log(Exp)) %>%
  summarise(Fleet_name = "Survey.total", N = n(), RMSE = sqrt(mean(residual^2, na.rm = TRUE))) 


RMSE_df <- rbind(RMSE_ind_survey_df, RMSE_all_survey_df) %>%
  add_column(Scenario = "MIT.1 Model")

write.csv(RMSE_df, "results/stock assessment/quahog/model fit/RMSE_MIT.1.csv", row.names = FALSE)

remove(assessment_results, RMSE_ind_survey_df, RMSE_all_survey_df, RMSE_df)



## 3.6 combine and save ----

RMSE_base_df <- read.csv("results/stock assessment/quahog/model fit/RMSE_base.csv")
RMSE_WEE_df <- read.csv("results/stock assessment/quahog/model fit/RMSE_WEE.csv")
RMSE_MIT.1_df <- read.csv("results/stock assessment/quahog/model fit/RMSE_MIT.1.csv")


RMSE_combine_df <- rbind(RMSE_base_df, RMSE_WEE_df, RMSE_MIT.1_df) %>%
  add_column(Species = "Quahog") %>%
  write.csv("results/stock assessment/quahog/model fit/RMSE_combined.csv", row.names = FALSE)


# ------------------------------------------------------------------------ #
