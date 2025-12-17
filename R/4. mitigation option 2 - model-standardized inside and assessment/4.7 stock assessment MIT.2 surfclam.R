library(r4ss)
library(tidyverse)


# for surfclam, all adjusted indices for the scenarios are generated from the R script 
# E:\Stony Brook job\NYSERDA Offshore wind\R\1.7 tows and indices for assessment\1.7.4 surfclam (ratio)


# 1. SS3 input adjustment for MIT.2  ----

data <- SS_readdat("assessment model/surfclam/base model/data.ss")
data$Nsurveys # 6 surveys, RD means research dredge, MCD means modified commercial dredge
data$fleetnames # need to play with RDtrendS, RDscaleS, and MCDs
data$CPUEinfo # units are all numbers according to data.ss, 


## 2.1 RDtrendS: index 3, numbers per m2, 1982-2011 ----

load("results/indices for assessment/surfclam/mitigation_2_model-standardized/MIT.2.RDtrendS_ratio.Rdata")

RDtrendS <- data$CPUE %>%
  filter(index == 3) %>%
  mutate(obs = obs * RDtrendS_MIT.2_VALUE_ratio,
         se_log = se_log * RDtrendS_MIT.2_STDERR_ratio)

data$CPUE[data$CPUE$index == 3, ] <- RDtrendS



## 2.2 RDscaleS: index 4, numbers per tow using the more precise sensor tow distances, 1997-2011 ----

load("results/indices for assessment/surfclam/mitigation_2_model-standardized/MIT.2.RDscaleS_ratio.Rdata")

RDscaleS <- data$CPUE %>%
  filter(index == 4) %>%
  mutate(obs = obs * RDscaleS_MIT.2_VALUE_ratio,
         se_log = se_log * RDscaleS_MIT.2_STDERR_ratio)

data$CPUE[data$CPUE$index == 4, ] <- RDscaleS



## 2.3 MCDS: index 4, numbers per tow using the more precise sensor tow distances, 1997-2011 ----

load("results/indices for assessment/surfclam/mitigation_2_model-standardized/MIT.2.MCDS_ratio.Rdata")

MCDS <- data$CPUE %>%
  filter(index == 5) %>%
  mutate(obs = obs * MCDS_MIT.2_VALUE_ratio,
         se_log = se_log * MCDS_MIT.2_STDERR_ratio)

data$CPUE[data$CPUE$index == 5, ] <- MCDS


## 2.4 save the MIT.2 input into a data.ss ----


SS_writedat(datlist = data, outfile = "assessment model/surfclam/MIT.2 model/data.ss", overwrite = TRUE)




# ------------------------------------------------------------------------------------------ #






# 2. run SS3 with the WEE input ----

  ## 2.1 run model ----

run_dir_MIT.2 <- "assessment model/surfclam/MIT.2 model/"

r4ss::run(dir = run_dir_MIT.2) # skip any folders that already contain a "Report.sso" file




  ## 2.2 extract outputs ----

    # create a list of quantities for the outputs
assessment_results_MIT.2 <- SS_output(run_dir_MIT.2, verbose = FALSE)
    # SS_plots(assessment_results)

    # see the names of items available
sort(names(assessment_results_MIT.2))





  ## 2.3 timeseries data and reference points ----

    ### 2.3.1 SSB ----

SSB_MIT.2_df <- assessment_results_MIT.2$derived_quants %>%
  filter(grepl("^SSB_\\d{4}$", Label))

# SSB-threshold is 1/4 of the maximum historical SSB 
SSB_MIT.2.MSY <- 1/2 * max(SSB_MIT.2_df$Value)
SSB_MIT.2.threshold <- 1/2 * SSB_MIT.2.MSY

SSB_MIT.2_df$Ratio <- SSB_MIT.2_df$Value/SSB_MIT.2.threshold
SSB_MIT.2_df$Ratio.SD <- SSB_MIT.2_df$StdDev/SSB_MIT.2.threshold

SSB_MIT.2_df$Scenario <- "MIT.2 Model"

write.csv(SSB_MIT.2_df, "results/stock assessment/surfclam/SSB_MIT.2.csv")


    ### 2.3.2 F ----

# F time series, F is combined with number-weighted mean
F_MIT.2_df <- assessment_results_MIT.2$timeseries %>%
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
Fstd_MIT.2_df <- assessment_results_MIT.2$derived_quants %>%
  filter(grepl("^F_\\d{4}$", Label)) %>%
  mutate(Yr = as.numeric(substr(Label, 3, 6)))  %>%
  select(Yr, Value, StdDev)


# combine F with Fstd
F_MIT.2_df <- F_MIT.2_df %>%
  left_join(Fstd_MIT.2_df) %>%
  mutate(F_weighted_std = F_weighted/Value*StdDev)  %>%
  select(-c(Value, StdDev))


# F threshold is defined based on an algorithm provided by Dan

source("R/functions/2.3.1 functions for surfclam stock assessment.r")

GetFref(rlst = assessment_results_MIT.2, rhoF = rho)

F_MIT.2_df$Ratio <- F_MIT.2_df$F_weighted/0.1526799
F_MIT.2_df$Ratio.SD <- F_MIT.2_df$F_weighted_std/0.1526799
F_MIT.2_df$Scenario <- "MIT.2 Model"


write.csv(F_MIT.2_df, "results/stock assessment/surfclam/F_MIT.2.csv")


    ### 2.3.3 catch limits ----

ABC_df <- assessment_results_MIT.2$derived_quants %>%
  filter(grepl("ForeCatch_\\d{4}$", Label)) %>%
  add_column(Scenario = "MIT.2 Model")

write.csv(ABC_df, "results/stock assessment/surfclam/ABC_MIT.2.csv", row.names = FALSE)



# ------------------------------------------------------------------------------------------ #










