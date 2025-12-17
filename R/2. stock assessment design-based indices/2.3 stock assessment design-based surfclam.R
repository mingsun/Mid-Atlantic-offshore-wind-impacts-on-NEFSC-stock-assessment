library(r4ss)
library(tidyverse)


# 1. run SS3 with original full dataset ----

  ## 1.1 run model ----
run_dir <- "assessment model/surfclam/base model/"

r4ss::run(dir = run_dir) # skip any folders that already contain a "Report.sso" file




  ## 1.2 extract outputs ----

# create a list of quantities for the outputs
assessment_results <- SS_output(run_dir, verbose = FALSE)
# SS_plots(assessment_results)

# see the names of items available
sort(names(assessment_results))





  ## 1.3 timeseries data and reference points ----

    ### 1.3.1 SSB ----

SSB_df <- assessment_results$derived_quants %>%
  filter(grepl("^SSB_\\d{4}$", Label))

  # SSB-threshold is 1/4 of the maximum historical SSB 
SSB.MSY <- 1/2 * max(SSB_df$Value)
SSB.threshold <- 1/2 * SSB.MSY

SSB_df$Ratio <- SSB_df$Value/SSB.threshold
SSB_df$Ratio.SD <- SSB_df$StdDev/SSB.threshold
SSB_df$Scenario <- "Base Model"

write.csv(SSB_df, "results/stock assessment/surfclam/SSB_original.csv")

    ### 1.3.2 F ----

  # F time series, F is combined with number-weighted mean
F_df <- assessment_results$timeseries %>%
  select(Yr, Area, 
         F1 = `F:_1`, F2 = `F:_2`,
         N1 = `SmryNum_SX:1_GP:1`, 
         N2 = `SmryNum_SX:1_GP:2`) %>%
  group_by(Yr) %>%
  summarize(
    F1 = max(F1, na.rm = TRUE),
    F2 = max(F2, na.rm = TRUE),
    N1 = max(N1, na.rm = TRUE),
    N2 = max(N2, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(F_weighted = (F1 * N1 + F2 * N2) / (N1 + N2))


  # grab the F std from the derived quants
Fstd_df <- assessment_results$derived_quants %>%
  filter(grepl("^F_\\d{4}$", Label)) %>%
  mutate(Yr = as.numeric(substr(Label, 3, 6)))  %>%
  select(Yr, Value, StdDev)


  # combine F with Fstd
F_df <- F_df %>%
  left_join(Fstd_df) %>%
  mutate(F_weighted_std = F_weighted/Value*StdDev)  %>%
  select(-c(Value, StdDev))



  # F threshold is defined based on an algorithm provided by Dan

source("R/functions/2.3.1 functions for surfclam stock assessment.r")

GetFref(rlst = assessment_results, rhoF = rho)

F_df$Ratio <- F_df$F_weighted/0.1526799
F_df$Ratio.SD <- F_df$F_weighted_std/0.1526799
F_df$Scenario <- "Base Model"

write.csv(F_df, "results/stock assessment/surfclam/F_original.csv")



  ### 1.3.3 catch limits ----

ABC_df <- assessment_results$derived_quants %>%
  filter(grepl("ForeCatch_\\d{4}$", Label)) %>%
  add_column(Scenario = "Base Model")
  
write.csv(ABC_df, "results/stock assessment/surfclam/ABC_original.csv", row.names = FALSE)



# ------------------------------------------------------------------------------------------ #




# 2. SS3 input adjustment for WEE  ----

data <- SS_readdat("assessment model/surfclam/base model/data.ss")
data$Nsurveys # 6 surveys, RD means research dredge, MCD means modified commercial dredge
data$fleetnames # need to play with RDtrendS, RDscaleS, and MCDs
data$CPUEinfo # units are all numbers according to data.ss, 


  ## 2.1 RDtrendS: index 3, numbers per m2, 1982-2011 ----

load("results/indices for assessment/surfclam/WEE/WEE.RDtrendS_ratio.Rdata")

RDtrendS <- data$CPUE %>%
  filter(index == 3) %>%
  mutate(obs = obs * RDtrendS_WEE_VALUE_ratio,
         se_log = se_log * RDtrendS_WEE_STDERR_ratio)

data$CPUE[data$CPUE$index == 3, ] <- RDtrendS



  ## 2.2 RDscaleS: index 4, numbers per tow using the more precise sensor tow distances, 1997-2011 ----

load("results/indices for assessment/surfclam/WEE/WEE.RDscaleS_ratio.Rdata")

RDscaleS <- data$CPUE %>%
  filter(index == 4) %>%
  mutate(obs = obs * RDscaleS_WEE_VALUE_ratio,
         se_log = se_log * RDscaleS_WEE_STDERR_ratio)

data$CPUE[data$CPUE$index == 4, ] <- RDscaleS



  ## 2.3 MCDS: index 4, numbers per tow using the more precise sensor tow distances, 1997-2011 ----

load("results/indices for assessment/surfclam/WEE/WEE.MCDS_ratio.Rdata")

MCDS <- data$CPUE %>%
  filter(index == 5) %>%
  mutate(obs = obs * MCDS_WEE_VALUE_ratio,
         se_log = se_log * MCDS_WEE_STDERR_ratio)

data$CPUE[data$CPUE$index == 5, ] <- MCDS


  ## 2.4 save the WEE input into a data.ss ----


SS_writedat(datlist = data, outfile = "assessment model/surfclam/WEE_design-based_indices_model/data.ss", overwrite = TRUE)



# ------------------------------------------------------------------------------------------ #






# 3. run SS3 with the WEE input ----

  ## 3.1 run model ----

run_dir_WEE <- "assessment model/surfclam/WEE_design-based_indices_model/"

r4ss::run(dir = run_dir_WEE) # skip any folders that already contain a "Report.sso" file




  ## 3.2 extract outputs ----

    # create a list of quantities for the outputs
assessment_results_WEE <- SS_output(run_dir_WEE, verbose = FALSE)
    # SS_plots(assessment_results)

    # see the names of items available
sort(names(assessment_results_WEE))





  ## 3.3 timeseries data and reference points ----

    ### 3.3.1 SSB ----

SSB_WEE_df <- assessment_results_WEE$derived_quants %>%
  filter(grepl("^SSB_\\d{4}$", Label))

# SSB-threshold is 1/4 of the maximum historical SSB 
SSB_WEE.MSY <- 1/2 * max(SSB_WEE_df$Value)
SSB_WEE.threshold <- 1/2 * SSB_WEE.MSY

SSB_WEE_df$Ratio <- SSB_WEE_df$Value/SSB_WEE.threshold
SSB_WEE_df$Ratio.SD <- SSB_WEE_df$StdDev/SSB_WEE.threshold

SSB_WEE_df$Scenario <- "WEE Model"

write.csv(SSB_WEE_df, "results/stock assessment/surfclam/SSB_WEE.csv")


    ### 3.3.2 F ----

# F time series, F is combined with number-weighted mean
F_WEE_df <- assessment_results_WEE$timeseries %>%
  select(Yr, Area, 
         F1 = `F:_1`, F2 = `F:_2`,
         N1 = `SmryNum_SX:1_GP:1`, 
         N2 = `SmryNum_SX:1_GP:2`) %>%
  group_by(Yr) %>%
  summarize(
    F1 = max(F1, na.rm = TRUE),
    F2 = max(F2, na.rm = TRUE),
    N1 = max(N1, na.rm = TRUE),
    N2 = max(N2, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(F_weighted = (F1 * N1 + F2 * N2) / (N1 + N2))

# grab the F std from the derived quants
Fstd_WEE_df <- assessment_results_WEE$derived_quants %>%
  filter(grepl("^F_\\d{4}$", Label)) %>%
  mutate(Yr = as.numeric(substr(Label, 3, 6)))  %>%
  select(Yr, Value, StdDev)


# combine F with Fstd
F_WEE_df <- F_WEE_df %>%
  left_join(Fstd_WEE_df) %>%
  mutate(F_weighted_std = F_weighted/Value*StdDev)  %>%
  select(-c(Value, StdDev))


# F threshold is defined based on an algorithm provided by Dan

source("R/functions/2.3.1 functions for surfclam stock assessment.r")

GetFref(rlst = assessment_results_WEE, rhoF = rho)

F_WEE_df$Ratio <- F_WEE_df$F_weighted/0.1526799
F_WEE_df$Ratio.SD <- F_WEE_df$F_weighted_std/0.1526799
F_WEE_df$Scenario <- "WEE Model"


write.csv(F_WEE_df, "results/stock assessment/surfclam/F_WEE.csv")


    ### 3.3.3 catch limits ----

ABC_df <- assessment_results_WEE$derived_quants %>%
  filter(grepl("ForeCatch_\\d{4}$", Label)) %>%
  add_column(Scenario = "WEE Model")

write.csv(ABC_df, "results/stock assessment/surfclam/ABC_WEE.csv", row.names = FALSE)



# ------------------------------------------------------------------------------------------ #










