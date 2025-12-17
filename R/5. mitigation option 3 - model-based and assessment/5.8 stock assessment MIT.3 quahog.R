library(r4ss)
library(tidyverse)


# for quahog, all adjusted indices for the scenarios are generated from the R script 
# E:\Stony Brook job\NYSERDA Offshore wind\R\1.7 tows and indices for assessment\1.7.4 quahog (ratio)


# 1. SS3 input adjustment for MIT.3  ----

data <- SS_readdat("assessment model/quahog/base model/Wquahog1.dat")
data$Nsurveys # 6 surveys, RD means research dredge, MCD means modified commercial dredge
data$fleetnames # need to play with RDtrendS, RDscaleS, and MCDs
data$CPUEinfo # units are all numbers according to data.ss, 


  ## 2.1 RDtrendS: index 3, numbers per m2, 1982-2011 ----

load("results/indices for assessment/quahog/mitigation_3_model-based indices/MIT.3.RDtrendS_ratio.Rdata")

RDtrendS <- data$CPUE %>%
  filter(index == 3) %>%
  mutate(obs = obs * RDtrendS_MIT.3_VALUE_ratio,
         se_log = se_log * RDtrendS_MIT.3_STDERR_ratio)

data$CPUE[data$CPUE$index == 3, ] <- RDtrendS



  ## 2.2 RDscaleS: index 4, numbers per tow using the more precise sensor tow distances, 1997-2011 ----

load("results/indices for assessment/quahog/mitigation_3_model-based indices/MIT.3.RDscaleS_ratio.Rdata")

RDscaleS <- data$CPUE %>%
  filter(index == 4) %>%
  mutate(obs = obs * RDscaleS_MIT.3_VALUE_ratio,
         se_log = se_log * RDscaleS_MIT.3_STDERR_ratio)

data$CPUE[data$CPUE$index == 4, ] <- RDscaleS



  ## 2.3 MCDS: index 4, numbers per tow using the more precise sensor tow distances, 1997-2011 ----

load("results/indices for assessment/quahog/mitigation_3_model-based indices/MIT.3.MCDS_ratio.Rdata")

MCDS <- data$CPUE %>%
  filter(index == 5) %>%
  mutate(obs = obs * MCDS_MIT.3_VALUE_ratio,
         se_log = se_log * MCDS_MIT.3_STDERR_ratio)

data$CPUE[data$CPUE$index == 5, ] <- MCDS


  ## 2.4 save the MIT.3 input into a data.ss ----


SS_writedat(datlist = data, outfile = "assessment model/quahog/MIT.3 model/Wquahog1.dat", overwrite = TRUE)




# ------------------------------------------------------------------------------------------ #






# 2. run SS3 with the WEE input ----

  ## 2.1 run model ----

run_dir_MIT.3 <- "assessment model/quahog/MIT.3 model/"

r4ss::run(dir = run_dir_MIT.3) # skip any folders that already contain a "Report.sso" file




  ## 2.2 extract outputs ----

    # create a list of quantities for the outputs
assessment_results_MIT.3 <- SS_output(run_dir_MIT.3, verbose = FALSE)
    # SS_plots(assessment_results)

    # see the names of items available
sort(names(assessment_results_MIT.3))





  ## 2.3 timeseries data and reference points ----

    ### 2.3.1 SSB ----

SSB_MIT.3_df <- assessment_results_MIT.3$derived_quants %>%
  filter(grepl("^SSB_\\d{4}$", Label))

# SSB-threshold is 1/4 of the maximum historical SSB 
SSB_MIT.3.MSY <- 1/2 * max(SSB_MIT.3_df$Value)
SSB_MIT.3.threshold <- 1/2 * SSB_MIT.3.MSY

SSB_MIT.3_df$Ratio <- SSB_MIT.3_df$Value/SSB_MIT.3.threshold
SSB_MIT.3_df$Ratio.SD <- SSB_MIT.3_df$StdDev/SSB_MIT.3.threshold

SSB_MIT.3_df$Scenario <- "MIT.3 Model"

write.csv(SSB_MIT.3_df, "results/stock assessment/quahog/SSB_MIT.3.csv")


    ### 2.3.2 F ----

# F time series, F is combined with number-weighted mean
F_MIT.3_df <- assessment_results_MIT.3$timeseries %>%
  filter(Yr <= 2020) %>%
  select(Yr, Area, 
         F1 = `F:_1`, F2 = `F:_2`,
         N1 = `SmryNum_SX:1_GP:1`) %>%
  mutate(F = if_else(Area == 1, F1, F2),
         N = N1) %>%
  select(Yr, Area, F, N) %>%
  pivot_wider(names_from = Area, values_from = c(F, N), names_prefix = "") %>%
  rename(F1 = `F_1`, F2 = `F_2`, N1 = `N_1`, N2 = `N_2`) %>%
  mutate(F_weighted = (F1 * N1 + F2 * N2) / (N1 + N2))

# grab the F std from the derived quants
Fstd_MIT.3_df <- assessment_results_MIT.3$derived_quants %>%
  filter(grepl("^F_\\d{4}$", Label)) %>%
  mutate(Yr = as.numeric(substr(Label, 3, 6)))  %>%
  select(Yr, Value, StdDev)


# combine F with Fstd
F_MIT.3_df <- F_MIT.3_df %>%
  left_join(Fstd_MIT.3_df) %>%
  mutate(F_weighted_std = F_weighted/Value*StdDev)  %>%
  select(-c(Value, StdDev))


# FMSY is F.Threshold is 0.019

F_MIT.3_df$Ratio <- F_MIT.3_df$F_weighted/0.019
F_MIT.3_df$Ratio.SD <- F_MIT.3_df$F_weighted_std/0.019

F_MIT.3_df$Scenario <- "MIT.3 Model"

write.csv(F_MIT.3_df, "results/stock assessment/quahog/F_MIT.3.csv")


    ### 2.3.3 catch limits ----

ABC_df <- assessment_results_MIT.3$derived_quants %>%
  filter(grepl("ForeCatch_\\d{4}$", Label)) %>%
  add_column(Scenario = "MIT.3 Model")

write.csv(ABC_df, "results/stock assessment/quahog/ABC_MIT.3.csv", row.names = FALSE)



# ------------------------------------------------------------------------------------------ #










