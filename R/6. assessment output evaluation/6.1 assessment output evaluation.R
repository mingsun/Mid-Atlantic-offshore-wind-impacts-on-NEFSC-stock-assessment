library(tidyverse)
library(ASAPplots)
library(r4ss)

# 1. load results ----

  ## summer flounder ASAP ----

SF_ASAP <- list.files("results/stock assessment/summer flounder/", pattern = "^ASAP_summary.*\\.csv$", full.names = TRUE) %>%
  lapply(read.csv) %>%
  bind_rows()



SF.ASAP.base <- read.csv("results/stock assessment/summer flounder/ASAP_summary_ASAP3_MTA2023_FINAL.csv")
SF.ASAP.WEE <- read.csv("results/stock assessment/summer flounder/ASAP_summary_ASAP3_MTA2023_WEE_BTS.csv")



SF.ASAP.MIT.1 <- read.csv("results/stock assessment/summer flounder/ASAP_summary_ASAP3_MTA2023_MIT.1_BTS.csv")
SF.ASAP.MIT.2 <- read.csv("results/stock assessment/summer flounder/ASAP_summary_ASAP3_MTA2023_MIT.2_BTS.csv")
SF.ASAP.MIT.3 <- read.csv("results/stock assessment/summer flounder/ASAP_summary_ASAP3_MTA2023_MIT.3_BTS.csv")


  ## squid index-based ---- 
LS.INDEX.base <- read.csv("results/stock assessment/squid/index-based summary_base.csv")
LS.INDEX.WEE <- read.csv("results/stock assessment/squid/index-based summary_WEE.csv")
LS.INDEX.MIT.1 <- read.csv("results/stock assessment/squid/index-based summary_MIT.1.csv")
LS.INDEX.MIT.2 <- read.csv("results/stock assessment/squid/index-based summary_MIT.2.csv")
LS.INDEX.MIT.3 <- read.csv("results/stock assessment/squid/index-based summary_MIT.3.csv")


  ## surfclam SS ----
SSB_AS_base_df <- read.csv("results/stock assessment/surfclam/SSB_original.csv")
SSB_AS_WEE_df <- read.csv("results/stock assessment/surfclam/SSB_WEE.csv")
SSB_AS_MIT.1_df <- read.csv("results/stock assessment/surfclam/SSB_MIT.1.csv")
SSB_AS_MIT.1_df <- read.csv("results/stock assessment/surfclam/SSB_MIT.1.csv")
SSB_AS_MIT.1_df <- read.csv("results/stock assessment/surfclam/SSB_MIT.1.csv")

F_AS_base_df <- read.csv("results/stock assessment/surfclam/F_original.csv")
F_AS_WEE_df <- read.csv("results/stock assessment/surfclam/F_WEE.csv")
F_AS_MIT.1_df <- read.csv("results/stock assessment/surfclam/F_MIT.1.csv")
F_AS_MIT.1_df <- read.csv("results/stock assessment/surfclam/F_MIT.1.csv")
F_AS_MIT.1_df <- read.csv("results/stock assessment/surfclam/F_MIT.1.csv")

ABC_AS_base_df <-  read.csv("results/stock assessment/surfclam/ABC_original.csv")
ABC_AS_WEE_df <-  read.csv("results/stock assessment/surfclam/ABC_WEE.csv")
ABC_AS_MIT.1_df <- read.csv("results/stock assessment/surfclam/ABC_MIT.1.csv")
ABC_AS_MIT.1_df <- read.csv("results/stock assessment/surfclam/ABC_MIT.1.csv")
ABC_AS_MIT.1_df <- read.csv("results/stock assessment/surfclam/ABC_MIT.1.csv")


  ## quahog SS ----
SSB_OQ_base_df <- read.csv("results/stock assessment/quahog/SSB_original.csv")
SSB_OQ_WEE_df <- read.csv("results/stock assessment/quahog/SSB_WEE.csv")
SSB_OQ_MIT.1_df <- read.csv("results/stock assessment/quahog/SSB_MIT.1.csv")

  
F_OQ_base_df <- read.csv("results/stock assessment/quahog/F_original.csv")
F_OQ_WEE_df <- read.csv("results/stock assessment/quahog/F_WEE.csv")
F_OQ_MIT.1_df <- read.csv("results/stock assessment/quahog/F_MIT.1.csv")

ABC_OQ_base_df <-  read.csv("results/stock assessment/quahog/ABC_original.csv")
ABC_OQ_WEE_df <-  read.csv("results/stock assessment/quahog/ABC_WEE.csv")
ABC_OQ_MIT.1_df <-  read.csv("results/stock assessment/quahog/ABC_MIT.1.csv")




# ------------------------------------------------------------------------ #





# 2. SSB ----

  ## 2.1 summer flounder ----

SSB_SF_df <- rbind(SF.ASAP.base, SF.ASAP.WEE, SF.ASAP.MIT.1) %>%
  select(Scenario, Year, SSB, SSB_95_lo, SSB_95_hi)

write.csv(SSB_SF_df, "results/stock assessment/summer flounder/SSB_combined_df.csv", row.names = FALSE)

SSB_SF_plot <- ggplot(SSB_SF_df, aes(x = Year)) +
  geom_errorbar(aes(ymin = SSB_95_lo, ymax = SSB_95_hi, color = Scenario), width = 0.2, position = position_dodge(width = 0.5)) +
  geom_point(aes(y = SSB, color = Scenario), position = position_dodge(width = 0.5)) +
  geom_line(aes(y = SSB, color = Scenario)) +
  # scale_color_manual(values = c("steelblue3", "indianred2"), labels = c("Base Model", "WEE Model"), name = "Time Series") +
  ylab("SSB [mt]") +
  theme(legend.position = c(0.05, 0.95), # Place the legend in the top-left corner
        legend.justification = c(0, 1))

remove(SSB_SF_plot)



  ## 2.2 squid ----

SSB_LS_df <- rbind(LS.INDEX.base, LS.INDEX.WEE, LS.INDEX.MIT.1)
write.csv(SSB_LS_df, "results/stock assessment/squid/B_Indices_combined_df.csv", row.names = FALSE)


SSB_LS_plot <- ggplot(SSB_LS_df, aes(x = Year)) +
  geom_point(aes(y = Annualized.B, color = Scenario), position = position_dodge(width = 0.5)) +
  geom_line(aes(y = Annualized.B, color = Scenario)) +
  # scale_color_manual(values = c("steelblue3", "indianred2"), labels = c("Base Model", "WEE Model"), name = "Time Series") +
  ylab("Annualized Biomass [mt]") +
  theme(legend.position = c(0.05, 0.95), # Place the legend in the top-left corner
        legend.justification = c(0, 1))

remove(SSB_LS_df, SSB_LS_plot)



  ## 2.3 surfclam ----

SSB_AS_df <- rbind(SSB_AS_base_df, SSB_AS_WEE_df, SSB_AS_MIT.1_df) %>%
  mutate(Year = as.numeric(substr(Label, 5, 8))) %>%
  select(Scenario, Year, Value) %>%
  filter(Year <= 2023)

write.csv(SSB_AS_df, "results/stock assessment/surfclam/SSB_combined_df.csv", row.names = FALSE)

  
SSB_AS_plot <- ggplot(SSB_AS_df, aes(x = Year)) +
  geom_point(aes(y = Value, color = Scenario), position = position_dodge(width = 0.5)) +
  geom_line(aes(y = Value, color = Scenario)) +
  # scale_color_manual(values = c("steelblue3", "indianred2"), labels = c("Base Model", "WEE Model"), name = "Time Series") +
  ylab("Annualized Biomass [mt]") +
  theme(legend.position = c(0.05, 0.95), # Place the legend in the top-left corner
        legend.justification = c(0, 1))

remove(SSB_AS_df, SSB_AS_base_df, SSB_AS_WEE_df, SSB_AS_plot)


  ## 2.4 quahog ----

SSB_OQ_df <- rbind(SSB_OQ_base_df, SSB_OQ_WEE_df, SSB_OQ_MIT.1_df) %>%
  mutate(Year = as.numeric(substr(Label, 5, 8))) %>%
  select(Scenario, Year, Value) %>%
  filter(Year <= 2019)

write.csv(SSB_OQ_df, "results/stock assessment/quahog/SSB_combined_df.csv", row.names = FALSE)


SSB_OQ_plot <- ggplot(SSB_OQ_df, aes(x = Year)) +
  geom_point(aes(y = Value, color = Scenario), position = position_dodge(width = 0.5)) +
  geom_line(aes(y = Value, color = Scenario)) +
  # scale_color_manual(values = c("steelblue3", "indianred2"), labels = c("Base Model", "WEE Model"), name = "Time Series") +
  ylab("Annualized Biomass [mt]") +
  theme(legend.position = c(0.05, 0.95), # Place the legend in the top-left corner
        legend.justification = c(0, 1))

remove(SSB_OQ_df, SSB_OQ_base_df, SSB_OQ_WEE_df, SSB_OQ_plot)





# ------------------------------------------------------------------------ #




# 3. F ----

  ## 3.1 summer flounder ----

F_SF_df <- rbind(SF.ASAP.base, SF.ASAP.WEE, SF.ASAP.MIT.1) %>%
  select(Scenario, Year, Freport, Freport_95_lo, Freport_95_hi) 

write.csv(F_SF_df, "results/stock assessment/summer flounder/F_combined_df.csv", row.names = FALSE)

F_SF_plot <- ggplot(F_SF_df, aes(x = Year)) +
  geom_errorbar(aes(ymin = Freport_95_lo, ymax = Freport_95_hi, color = Scenario), width = 0.2, position = position_dodge(width = 0.5)) +
  geom_point(aes(y = Freport, color = Scenario), position = position_dodge(width = 0.5)) +
  geom_line(aes(y = Freport, color = Scenario)) +
  # scale_color_manual(values = c("steelblue3", "indianred2"), labels = c("Base Model", "WEA_BTS Model"), name = "Time Series") +
  ylab("Fishing Mortality") +
  theme(legend.position = c(0.05, 0.95), # Place the legend in the top-left corner
        legend.justification = c(0, 1))

remove( F_SF_plot)



  ## 3.2 squid ----

F_LS_df <- rbind(LS.INDEX.base, LS.INDEX.WEE, LS.INDEX.MIT.1)

# csv already saved

F_LS_plot <- ggplot(F_LS_df, aes(x = Year)) +
  geom_point(aes(y = Annualized.Exploitation.Indices, color = Scenario), position = position_dodge(width = 0.5)) +
  geom_line(aes(y = Annualized.Exploitation.Indices, color = Scenario)) +
  # scale_color_manual(values = c("steelblue3", "indianred2"), labels = c("Base Model", "WEE Model"), name = "Time Series") +
  ylab("Annualized Exploitation Indices") +
  theme(legend.position = c(0.05, 0.95), # Place the legend in the top-left corner
        legend.justification = c(0, 1))

remove(F_LS_df, F_LS_plot)



  ## 3.3 surfclam ----

F_AS_df <- rbind(F_AS_base_df, F_AS_WEE_df, F_AS_MIT.1_df) %>%
  rename(Year = Yr) %>%
  select(Scenario, Year, F_weighted)  %>%
  filter(Year <= 2023)

write.csv(F_AS_df, "results/stock assessment/surfclam/F_combined_df.csv", row.names = FALSE)

F_AS_plot <- ggplot(F_AS_df, aes(x = Year)) +
  geom_point(aes(y = F_weighted, color = Scenario), position = position_dodge(width = 0.5)) +
  geom_line(aes(y = F_weighted, color = Scenario)) +
  # scale_color_manual(values = c("steelblue3", "indianred2"), labels = c("Base Model", "WEE Model"), name = "Time Series") +
  ylab("Annualized Biomass [mt]") +
  theme(legend.position = c(0.05, 0.95), # Place the legend in the top-left corner
        legend.justification = c(0, 1))

remove(F_AS_df, F_AS_base_df, F_AS_WEE_df, F_AS_plot)



  ## 3.4 quahog ----

F_OQ_df <- rbind(F_OQ_base_df, F_OQ_WEE_df, F_OQ_MIT.1_df) %>%
  rename(Year = Yr) %>%
  select(Scenario, Year, F_weighted) %>%
  filter(Year <= 2019)

write.csv(F_OQ_df, "results/stock assessment/quahog/F_combined_df.csv", row.names = FALSE)

F_OQ_plot <- ggplot(F_OQ_df, aes(x = Year)) +
  geom_point(aes(y = F_weighted, color = Scenario), position = position_dodge(width = 0.5)) +
  geom_line(aes(y = F_weighted, color = Scenario)) +
  # scale_color_manual(values = c("steelblue3", "indianred2"), labels = c("Base Model", "WEE Model"), name = "Time Series") +
  ylab("Annualized Biomass [mt]") +
  theme(legend.position = c(0.05, 0.95), # Place the legend in the top-left corner
        legend.justification = c(0, 1))

remove(F_OQ_df, F_OQ_base_df, F_OQ_WEE_df, F_OQ_plot)




# ------------------------------------------------------------------------ #





# 4. Stock Status ----

  ## 4.1 summer flounder ----

Status_SF_df <- merge(F_SF_df, SSB_SF_df) %>%
  filter(Year == 2022) %>%
  add_column(FMSY = 0.451, SSBMSY = 49561) %>% 
  # translate from wide to long
  pivot_longer(cols = c(Freport, SSB), names_to = "Metric", values_to = "Value") %>%
  mutate(Lower_CI = case_when(Metric == "Freport" ~ Freport_95_lo, Metric == "SSB" ~ SSB_95_lo),
         Upper_CI = case_when(Metric == "Freport" ~ Freport_95_hi, Metric == "SSB" ~ SSB_95_hi),
         MSY_Ref  = case_when(Metric == "Freport" ~ FMSY,Metric == "SSB" ~ SSBMSY)) %>%
  select(Scenario, Year, Metric, Value, Lower_CI, Upper_CI, MSY_Ref) %>%
  # convert to relative value
  mutate(Value = Value/MSY_Ref, Lower_CI = Lower_CI/MSY_Ref, Upper_CI = Upper_CI/MSY_Ref)
  
Status_SF_plot <- ggplot(Status_SF_df, aes(x = Scenario)) +
  geom_errorbar(aes(ymin = Lower_CI, ymax = Upper_CI), width = 0.2, linetype = 2) +
  geom_point(aes(y = Value)) +
  geom_hline(yintercept = 1) +
  facet_wrap(.~ Metric, scale = "free_y") +
  # scale_color_manual(values = c("steelblue3", "indianred2"), labels = c("Base Model", "WEA_BTS Model"), name = "Time Series") +
  ylab("Stock Status") 

remove(Status_SF_df, Status_SF_plot)


  ## 4.2 squid ----

Status_LS_df <- rbind(LS.INDEX.base, LS.INDEX.WEE, LS.INDEX.MIT.1)

Status_LS_plot <- ggplot(Status_LS_df, aes(x = Year)) +
  geom_point(aes(y = Annualized.Exploitation.Indices.final, color = Scenario), position = position_dodge(width = 0.5)) +
  geom_line(aes(y = Annualized.Exploitation.Indices.final, color = Scenario)) +
  # scale_color_manual(values = c("steelblue3", "indianred2"), labels = c("Base Model", "WEE Model"), name = "Time Series") +
  ylab("Status Indices") +
  theme(legend.position = c(0.05, 0.95), # Place the legend in the top-left corner
        legend.justification = c(0, 1))

remove(Status_LS_df, Status_LS_plot)


  ## 4.3 surfclam ----

SSB_status_AS_df <- rbind(SSB_AS_base_df, SSB_AS_WEE_df, SSB_AS_MIT.1_df) %>%
  mutate(Year = as.numeric(substr(Label, 5, 8))) %>%
  select(Scenario, Year, Ratio, Ratio.SD)  %>%
  add_column(Metric = "SSB") 

F_status_AS_df <- rbind(F_AS_base_df, F_AS_WEE_df, F_AS_MIT.1_df) %>%
  rename(Year = Yr) %>%
  select(Scenario, Year, Ratio, Ratio.SD) %>%
  add_column(Metric = "F")

Status_AS_df <- rbind(SSB_status_AS_df, F_status_AS_df) %>%
  filter(Year == 2023) %>%
  mutate(Lower_CI = Ratio - 1.96 * Ratio.SD,
         Upper_CI = Ratio + 1.96 * Ratio.SD)


# remove(SSB_status_AS_df, F_status_AS_df)



Status_AS_plot <- ggplot(Status_AS_df, aes(x = Scenario)) +
  geom_errorbar(aes(ymin = Lower_CI, ymax = Upper_CI), width = 0.2, linetype = 2) +
  geom_point(aes(y = Ratio)) +
  geom_hline(yintercept = 1) +
  facet_wrap(.~ Metric, scale = "free_y") +
  # scale_color_manual(values = c("steelblue3", "indianred2"), labels = c("Base Model", "WEA_BTS Model"), name = "Time Series") +
  ylab("Stock Status") 

remove(Status_AS_df, Status_AS_plot)


  ## 4.4 quahog ----

SSB_status_OQ_df <- rbind(SSB_OQ_base_df, SSB_OQ_WEE_df, SSB_OQ_MIT.1_df) %>%
  mutate(Year = as.numeric(substr(Label, 5, 8))) %>%
  select(Scenario, Year, Ratio, Ratio.SD)  %>%
  add_column(Metric = "SSB")

F_status_OQ_df <- rbind(F_OQ_base_df, F_OQ_WEE_df, F_OQ_MIT.1_df) %>%
  rename(Year = Yr) %>%
  select(Scenario, Year, Ratio, Ratio.SD) %>%
  add_column(Metric = "F")

Status_OQ_df <- rbind(SSB_status_OQ_df, F_status_OQ_df) %>%
  filter(Year == 2019) %>%
  mutate(Lower_CI = Ratio - 1.96 * Ratio.SD,
         Upper_CI = Ratio + 1.96 * Ratio.SD)


remove(SSB_status_OQ_df, F_status_OQ_df)

Status_OQ_plot <- ggplot(Status_OQ_df, aes(x = Scenario)) +
  geom_errorbar(aes(ymin = Lower_CI, ymax = Upper_CI), width = 0.2, linetype = 2) +
  geom_point(aes(y = Ratio)) +
  geom_hline(yintercept = 1) +
  facet_wrap(.~ Metric, scale = "free_y") +
  # scale_color_manual(values = c("steelblue3", "indianred2"), labels = c("Base Model", "WEA_BTS Model"), name = "Time Series") +
  ylab("Stock Status") 

remove(Status_OQ_df, Status_OQ_plot)









# ------------------------------------------------------------------------ #






# 5. Catch limits ----

# this will be the instant catch advice for the year following the recenet assessment

  ## 5.1 summer flounder ---- 
    # values see the "notes" in E:\Stony Brook job\NYSERDA Offshore wind\assessment model\summer flounder





  ## 5.2 squid ----
    # they are pasted here from the spreadsheet E:\Stony Brook job\NYSERDA Offshore wind\assessment model\squid\base model

ABC_LS_df <- rbind(data.frame(Scenario = c("Base Model",  Year = 2023, Value = 22.315, Lower_CI = NA, Upper_CI = NA, Metric = "ABC")),
                   data.frame(Scenario = c("WEE Model",   Year = 2023, Value = 22.943, Lower_CI = NA, Upper_CI = NA, Metric = "ABC")),
                   data.frame(Scenario = c("MIT.1 Model", Year = 2023, Value = 23.742, Lower_CI = NA, Upper_CI = NA, Metric = "ABC")))
  




  ## 5.3 surfclam ----


ABC_AS_df <- rbind(ABC_AS_base_df, ABC_AS_WEE_df, ABC_AS_MIT.1_df) %>%
  mutate(Year = as.numeric(substr(Label, 11, 14)),
         Lower_CI = Value - 1.96 * StdDev,
         Upper_CI = Value + 1.96 * StdDev) %>%
  select(Scenario, Year, Value, Lower_CI, Upper_CI)  %>%
  filter(Year == 2024) %>%
  add_column(Metric = "ABC")


ABC_AS_plot <- ggplot(ABC_AS_df, aes(x = Scenario)) +
  geom_errorbar(aes(ymin = Lower_CI, ymax = Upper_CI), width = 0.2, linetype = 2) +
  geom_point(aes(y = Value)) +
  # scale_color_manual(values = c("steelblue3", "indianred2"), labels = c("Base Model", "WEA_BTS Model"), name = "Time Series") +
  ylab("ABC") 

remove(ABC_AS_df, ABC_AS_plot)



  ## 5.4 quahog ----

ABC_OQ_df <- rbind(ABC_OQ_base_df, ABC_OQ_WEE_df, ABC_OQ_MIT.1_df) %>%
  mutate(Year = as.numeric(substr(Label, 11, 14)),
         Lower_CI = Value - 1.96 * StdDev,
         Upper_CI = Value + 1.96 * StdDev) %>%
  select(Scenario, Year, Value, Lower_CI, Upper_CI)  %>%
  filter(Year == 2020) %>%
  add_column(Metric = "ABC")


ABC_OQ_plot <- ggplot(ABC_OQ_df, aes(x = Scenario)) +
  geom_errorbar(aes(ymin = Lower_CI, ymax = Upper_CI), width = 0.2, linetype = 2) +
  geom_point(aes(y = Value)) +
  # scale_color_manual(values = c("steelblue3", "indianred2"), labels = c("Base Model", "WEA_BTS Model"), name = "Time Series") +
  ylab("ABC") 

remove(ABC_OQ_df, ABC_OQ_plot)


# ------------------------------------------------------------------------ #














