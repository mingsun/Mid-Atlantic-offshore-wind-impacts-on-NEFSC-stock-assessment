library(tidyverse)
library(ASAPplots)
library(r4ss)

# retrospective error by species

# 1. load results ----

# ------------------------------------------------------------------------ #




# 1. summer flounder ----

  ## 1.1 base model ----

wd <- "assessment model/summer flounder/ASAP/base model/retrospective run/"
asap.name <- "ASAP3_MTA2023_FINAL_000" # need to use the suffix 000
asap <- dget("assessment model/summer flounder/ASAP/base model/retrospective run/ASAP3_MTA2023_FINAL_000.rdat") # need to use the suffix 000
od <- "results/stock assessment/summer flounder/retrospective error/"

PlotRetroWrapper(wd, asap.name, asap, save.plots = FALSE, od, plotf = "jpg") # note that the Mohn's rho is a five year average

read.csv("results/stock assessment/summer flounder/retrospective error/Retro.rho.values_ASAP3_MTA2023_FINAL_000.csv") %>%
  add_column(Scenario = "Base Model") %>%
  write.csv("results/stock assessment/summer flounder/retrospective error/ASAP_retro_error_final.csv", row.names = FALSE)

remove(wd, asap.name, asap, od)



  ## 1.2 WEE_design-based_indices_model assessment ----

wd <- "assessment model/summer flounder/ASAP/WEE_design-based_indices_model/retrospective run/"
asap.name <- "ASAP3_MTA2023_WEE_BTS_000" # need to use the suffix 000
asap <- dget("assessment model/summer flounder/ASAP/WEE_design-based_indices_model/retrospective run/ASAP3_MTA2023_WEE_BTS_000.rdat") # need to use the suffix 000
od <- "results/stock assessment/summer flounder/retrospective error/"

PlotRetroWrapper(wd, asap.name, asap, save.plots = FALSE, od, plotf = "jpg")

read.csv("results/stock assessment/summer flounder/retrospective error/Retro.rho.values_ASAP3_MTA2023_WEE_BTS_000.csv") %>%
  add_column(Scenario = "WEE Model") %>%
  write.csv("results/stock assessment/summer flounder/retrospective error/ASAP_retro_error_WEE_BTS.csv", row.names = FALSE)

remove(wd, asap.name, asap, od)




## 1.3 MIT.1 assessment ----
# get the variable needed

wd <- "assessment model/summer flounder/ASAP/MIT.1 model/retrospective run/"
asap.name <- "ASAP3_MTA2023_MIT.1_BTS_000"
asap <- dget("assessment model/summer flounder/ASAP/MIT.1 model/retrospective run/ASAP3_MTA2023_MIT.1_BTS_000.rdat") 
od <- "results/stock assessment/summer flounder/retrospective error/"


PlotRetroWrapper(wd, asap.name, asap, save.plots = FALSE, od, plotf = "jpg")

read.csv("results/stock assessment/summer flounder/retrospective error/Retro.rho.values_ASAP3_MTA2023_MIT.1_BTS_000.csv") %>%
  add_column(Scenario = "MIT.1 Model") %>%
  write.csv("results/stock assessment/summer flounder/retrospective error/ASAP_retro_error_MIT.1_BTS.csv", row.names = FALSE)

remove(wd, asap.name, asap, od)




  ## 1.6 plot ----

base_df <- read.csv("results/stock assessment/summer flounder/retrospective error/ASAP_retro_error_final.csv")
WEE_df <- read.csv("results/stock assessment/summer flounder/retrospective error/ASAP_retro_error_WEE_BTS.csv")
MIT.1_df <- read.csv("results/stock assessment/summer flounder/retrospective error/ASAP_retro_error_MIT.1_BTS.csv")

plot_df <- rbind(base_df, WEE_df, MIT.1_df) %>%
  filter(X == "Mohn.rho") %>%
  select(f.rho, ssb.rho, Scenario) %>%
  pivot_longer(cols = c(f.rho, ssb.rho), names_to = "Metric", values_to = "Value") %>%
  add_column(Species = "Summer Flounder")

write.csv(plot_df, "results/stock assessment/summer flounder/retrospective error/retro_error_combined.csv", row.names = FALSE)

rho_SF_plot <- ggplot(plot_df, aes(x = Scenario)) +
  geom_bar(aes(y = Value), stat = "identity") +
  geom_hline(yintercept = 0) +
  facet_wrap(.~ Metric) +
  ylab("5-year Mohn's rho") 
  

# plot_df <- rbind(base_df, WEE_df, MIT.1_df) %>%
#   filter(X != "Mohn.rho") %>%
#   select(f.rho, ssb.rho, Scenario) %>%
#   pivot_longer(cols = c(f.rho, ssb.rho), names_to = "Metric", values_to = "Value")
# 
# rho_SF_plot <- ggplot(plot_df, aes(x = Scenario, y = Value)) +
#   geom_bar(stat = "summary", fun = mean) +
#   geom_point() +
#   geom_hline(yintercept = 0) +
#   facet_wrap(.~ Metric) +
#   ylab("5-year Mohn's rho") 



remove(base_df, WEE_df, plot_df, rho_SF_plot)



# ------------------------------------------------------------------------ #







# 2. squid ----

# For squid with index-based method, we cannot do retrospective error
# we calculate the relationship between catch and final indices

source("R/functions/4.2.1 squid retro function.R")

  ## 2.1 base model ----

index_df <- read.csv("results/stock assessment/squid/index-based summary_base.csv") %>%
  filter(Year >= 1987 & !is.na(Annualized.Exploitation.Indices.final))

results <- peel_regression_analysis(index_df, max_peel = 5) %>%
  add_column(Scenario = "Base Model") %>%
  write.csv("results/stock assessment/squid/retrospective error/retro_error_base.csv", row.names = FALSE)

remove(index_df, results)


  ## 2.2 WEE model ----

index_df <- read.csv("results/stock assessment/squid/index-based summary_WEE.csv") %>%
  filter(Year >= 1987 & !is.na(Annualized.Exploitation.Indices.final))

results <- peel_regression_analysis(index_df, max_peel = 5) %>%
  add_column(Scenario = "WEE Model") %>%
  write.csv("results/stock assessment/squid/retrospective error/retro_error_WEE.csv", row.names = FALSE)

remove(index_df, results)


  ## 2.3 MIT.1 model ----

index_df <- read.csv("results/stock assessment/squid/index-based summary_MIT.1.csv") %>%
  filter(Year >= 1987 & !is.na(Annualized.Exploitation.Indices.final))

results <- peel_regression_analysis(index_df, max_peel = 5) %>%
  add_column(Scenario = "MIT.1 Model") %>%
  write.csv("results/stock assessment/squid/retrospective error/retro_error_MIT.1.csv", row.names = FALSE)

remove(index_df, results)


  ## 2.6 plot ----

base_df <- read.csv("results/stock assessment/squid/retrospective error/retro_error_base.csv")
WEE_df <- read.csv("results/stock assessment/squid/retrospective error/retro_error_WEE.csv")
MIT.1_df <- read.csv("results/stock assessment/squid/retrospective error/retro_error_MIT.1.csv")

plot_df <- rbind(base_df, WEE_df, MIT.1_df) %>%
  filter(!is.na(Peel)) %>%
  group_by(Scenario) %>%
  summarise(Mohn.rho = mean(RelativeDeviation)) %>%
  add_column(Metric = "f.rho")%>%
  add_column(Species = "Longfin Squid")

write.csv(plot_df, "results/stock assessment/squid/retrospective error/retro_error_combined.csv", row.names = FALSE)


rho_SF_plot <- ggplot(plot_df, aes(x = Scenario)) +
  geom_bar(aes(y = Value), stat = "identity") +
  geom_hline(yintercept = 0) +
  facet_wrap(.~ Metric) +
  ylab("5-year Mohn's rho") 

remove(base_df, WEE_df, plot_df, rho_SF_plot)





# ------------------------------------------------------------------------ #




# 3. surfclam ----

  ## 3.1 base model ----

  # do retro fitting 

retro(dir = "assessment model/surfclam/base model/retrospective run/",
      oldsubdir = "",
      newsubdir = "retrospectives",
      subdirstart = "retro",
      years = 0:-5)

  # evaluate retro results

retroModels <- SSgetoutput(dirvec = file.path("assessment model/surfclam/base model/retrospective run/", "retrospectives", paste("retro", 0:-5, sep = "")))

retroSummary <- SSsummarize(retroModels)
endyrvec <- retroSummary[["endyrs"]] + 0:-5
Mohnrho_list <- SSmohnsrho(retroSummary, endyrvec, startyr = 1963, verbose = TRUE)

  # generate a dataframe below in a format consistent with the ASAP output
Mohnrho_base_df <- data.frame(X = "Mohn.rho", f.rho = Mohnrho_list$F/5, ssb.rho = Mohnrho_list$SSB/5, recr.rho = Mohnrho_list$Rec/5, jan1b.rho = Mohnrho_list$Bratio,
           explb.rho = NA, stockn.rho = NA, Age.1 = NA, Age.2 = NA, Age.3 = NA, Age.4 = NA, Age.5 = NA, Age.6 = NA, Age.7 = NA, Age.8 = NA,
           Scenario = "Base Model")

write.csv(Mohnrho_base_df, "results/stock assessment/surfclam/retrospective error/SS_retro_error_base.csv", row.names = FALSE)

remove(retroModels, retroSummary, endyrvec, Mohnrho_list, Mohnrho_base_df)



  ## 3.2 WEE model ----

# do retro fitting 

retro(dir = "assessment model/surfclam/WEE_design-based_indices_model/retrospective run/",
      oldsubdir = "",
      newsubdir = "retrospectives",
      subdirstart = "retro",
      years = 0:-5)

# evaluate retro results

retroModels <- SSgetoutput(dirvec = file.path("assessment model/surfclam/WEE_design-based_indices_model/retrospective run/", "retrospectives", paste("retro", 0:-5, sep = "")))

retroSummary <- SSsummarize(retroModels)
endyrvec <- retroSummary[["endyrs"]] + 0:-5
Mohnrho_list <- SSmohnsrho(retroSummary, endyrvec, startyr = 1963, verbose = TRUE)

# generate a dataframe below in a format consistent with the ASAP output
Mohnrho_WEE_df <- data.frame(X = "Mohn.rho", f.rho = Mohnrho_list$F/5, ssb.rho = Mohnrho_list$SSB/5, recr.rho = Mohnrho_list$Rec/5, jan1b.rho = Mohnrho_list$Bratio,
                              explb.rho = NA, stockn.rho = NA, Age.1 = NA, Age.2 = NA, Age.3 = NA, Age.4 = NA, Age.5 = NA, Age.6 = NA, Age.7 = NA, Age.8 = NA,
                              Scenario = "WEE Model")

write.csv(Mohnrho_WEE_df, "results/stock assessment/surfclam/retrospective error/SS_retro_error_WEE.csv", row.names = FALSE)

remove(retroModels, retroSummary, endyrvec, Mohnrho_list, Mohnrho_WEE_df)



  ## 3.3 MIT.1 model ----

# do retro fitting 

retro(dir = "assessment model/surfclam/MIT.1 model/retrospective run/",
      oldsubdir = "",
      newsubdir = "retrospectives",
      subdirstart = "retro",
      years = 0:-5)

# evaluate retro results

retroModels <- SSgetoutput(dirvec = file.path("assessment model/surfclam/MIT.1 model/retrospective run/", "retrospectives", paste("retro", 0:-5, sep = "")))

retroSummary <- SSsummarize(retroModels)
endyrvec <- retroSummary[["endyrs"]] + 0:-5
Mohnrho_list <- SSmohnsrho(retroSummary, endyrvec, startyr = 1963, verbose = TRUE)

# generate a dataframe below in a format consistent with the ASAP output
Mohnrho_MIT.1_df <- data.frame(X = "Mohn.rho", f.rho = Mohnrho_list$F/5, ssb.rho = Mohnrho_list$SSB/5, recr.rho = Mohnrho_list$Rec/5, jan1b.rho = Mohnrho_list$Bratio,
                             explb.rho = NA, stockn.rho = NA, Age.1 = NA, Age.2 = NA, Age.3 = NA, Age.4 = NA, Age.5 = NA, Age.6 = NA, Age.7 = NA, Age.8 = NA,
                             Scenario = "MIT.1 Model")

write.csv(Mohnrho_MIT.1_df, "results/stock assessment/surfclam/retrospective error/SS_retro_error_MIT.1.csv", row.names = FALSE)



  ## 3.3 MIT.2 model ----

# do retro fitting 

retro(dir = "assessment model/surfclam/MIT.2 model/retrospective run/",
      oldsubdir = "",
      newsubdir = "retrospectives",
      subdirstart = "retro",
      years = 0:-5)

# evaluate retro results

retroModels <- SSgetoutput(dirvec = file.path("assessment model/surfclam/MIT.2 model/retrospective run/", "retrospectives", paste("retro", 0:-5, sep = "")))

retroSummary <- SSsummarize(retroModels)
endyrvec <- retroSummary[["endyrs"]] + 0:-5
Mohnrho_list <- SSmohnsrho(retroSummary, endyrvec, startyr = 1963, verbose = TRUE)

# generate a dataframe below in a format consistent with the ASAP output
Mohnrho_MIT.2_df <- data.frame(X = "Mohn.rho", f.rho = Mohnrho_list$F/5, ssb.rho = Mohnrho_list$SSB/5, recr.rho = Mohnrho_list$Rec/5, jan1b.rho = Mohnrho_list$Bratio,
                               explb.rho = NA, stockn.rho = NA, Age.1 = NA, Age.2 = NA, Age.3 = NA, Age.4 = NA, Age.5 = NA, Age.6 = NA, Age.7 = NA, Age.8 = NA,
                               Scenario = "MIT.2 Model")

write.csv(Mohnrho_MIT.2_df, "results/stock assessment/surfclam/retrospective error/SS_retro_error_MIT.2.csv", row.names = FALSE)



  ## 3.4 MIT.3 model ----

# do retro fitting 

retro(dir = "assessment model/surfclam/MIT.3 model/retrospective run/",
      oldsubdir = "",
      newsubdir = "retrospectives",
      subdirstart = "retro",
      years = 0:-5)

# evaluate retro results

retroModels <- SSgetoutput(dirvec = file.path("assessment model/surfclam/MIT.3 model/retrospective run/", "retrospectives", paste("retro", 0:-5, sep = "")))

retroSummary <- SSsummarize(retroModels)
endyrvec <- retroSummary[["endyrs"]] + 0:-5
Mohnrho_list <- SSmohnsrho(retroSummary, endyrvec, startyr = 1963, verbose = TRUE)

# generate a dataframe below in a format consistent with the ASAP output
Mohnrho_MIT.3_df <- data.frame(X = "Mohn.rho", f.rho = Mohnrho_list$F/5, ssb.rho = Mohnrho_list$SSB/5, recr.rho = Mohnrho_list$Rec/5, jan1b.rho = Mohnrho_list$Bratio,
                               explb.rho = NA, stockn.rho = NA, Age.1 = NA, Age.2 = NA, Age.3 = NA, Age.4 = NA, Age.5 = NA, Age.6 = NA, Age.7 = NA, Age.8 = NA,
                               Scenario = "MIT.3 Model")

write.csv(Mohnrho_MIT.3_df, "results/stock assessment/surfclam/retrospective error/SS_retro_error_MIT.3.csv", row.names = FALSE)




  ## 3.6 plot ----

base_df <- read.csv("results/stock assessment/surfclam/retrospective error/SS_retro_error_base.csv")
WEE_df <- read.csv("results/stock assessment/surfclam/retrospective error/SS_retro_error_WEE.csv")
MIT.1_df <- read.csv("results/stock assessment/surfclam/retrospective error/SS_retro_error_MIT.1.csv")

plot_df <- rbind(base_df, WEE_df, MIT.1_df) %>%
  filter(X == "Mohn.rho") %>%
  select(f.rho, ssb.rho, Scenario) %>%
  pivot_longer(cols = c(f.rho, ssb.rho), names_to = "Metric", values_to = "Value") %>%
  add_column(Species = "Surfclam")

write.csv(plot_df, "results/stock assessment/surfclam/retrospective error/retro_error_combined.csv", row.names = FALSE)

rho_SF_plot <- ggplot(plot_df, aes(x = Scenario)) +
  geom_bar(aes(y = Value), stat = "identity") +
  geom_hline(yintercept = 0) +
  facet_wrap(.~ Metric) +
  ylab("5-year Mohn's rho") 


remove(base_df, WEE_df, plot_df, rho_SF_plot)

 





# ------------------------------------------------------------------------ #




# 4. quahog ----

  ## 4.1 base model ----

  # do retro fitting 

retro(dir = "assessment model/quahog/base model/retrospective run/",
      oldsubdir = "",
      newsubdir = "retrospectives",
      subdirstart = "retro",
      years = 0:-5)

# evaluate retro results

retroModels <- SSgetoutput(dirvec = file.path("assessment model/quahog/base model/retrospective run/", "retrospectives", paste("retro", 0:-5, sep = "")))

retroSummary <- SSsummarize(retroModels)
endyrvec <- retroSummary[["endyrs"]] + 0:-5
Mohnrho_list <- SSmohnsrho(retroSummary, endyrvec, startyr = 1982, verbose = TRUE)

# generate a dataframe below in a format consistent with the ASAP output
Mohnrho_base_df <- data.frame(X = "Mohn.rho", f.rho = Mohnrho_list$F/5, ssb.rho = Mohnrho_list$SSB/5, recr.rho = Mohnrho_list$Rec/5, jan1b.rho = Mohnrho_list$Bratio,
                              explb.rho = NA, stockn.rho = NA, Age.1 = NA, Age.2 = NA, Age.3 = NA, Age.4 = NA, Age.5 = NA, Age.6 = NA, Age.7 = NA, Age.8 = NA,
                              Scenario = "Base Model")

write.csv(Mohnrho_base_df, "results/stock assessment/quahog/retrospective error/SS_retro_error_base.csv", row.names = FALSE)

remove(retroModels, retroSummary, endyrvec, Mohnrho_list, Mohnrho_base_df)




  ## 4.2 WEE model ----

  # do retro fitting 

retro(dir = "assessment model/quahog/WEE_design-based_indices_model/retrospective run/",
      oldsubdir = "",
      newsubdir = "retrospectives",
      subdirstart = "retro",
      years = 0:-5)

  # evaluate retro results

retroModels <- SSgetoutput(dirvec = file.path("assessment model/quahog/WEE_design-based_indices_model/retrospective run/", "retrospectives", paste("retro", 0:-5, sep = "")))

retroSummary <- SSsummarize(retroModels)
endyrvec <- retroSummary[["endyrs"]] + 0:-5
Mohnrho_list <- SSmohnsrho(retroSummary, endyrvec, startyr = 1982, verbose = TRUE)

  # generate a dataframe below in a format consistent with the ASAP output
Mohnrho_WEE_df <- data.frame(X = "Mohn.rho", f.rho = Mohnrho_list$F/5, ssb.rho = Mohnrho_list$SSB/5, recr.rho = Mohnrho_list$Rec/5, jan1b.rho = Mohnrho_list$Bratio,
                             explb.rho = NA, stockn.rho = NA, Age.1 = NA, Age.2 = NA, Age.3 = NA, Age.4 = NA, Age.5 = NA, Age.6 = NA, Age.7 = NA, Age.8 = NA,
                             Scenario = "WEE Model")

write.csv(Mohnrho_WEE_df, "results/stock assessment/quahog/retrospective error/SS_retro_error_WEE.csv", row.names = FALSE)

remove(retroModels, retroSummary, endyrvec, Mohnrho_list, Mohnrho_WEE_df)




    ## 4.3 MIT.1 model ----

      # do retro fitting 

retro(dir = "assessment model/quahog/MIT.1 model/retrospective run/",
      oldsubdir = "",
      newsubdir = "retrospectives",
      subdirstart = "retro",
      years = 0:-5)

      # evaluate retro results

retroModels <- SSgetoutput(dirvec = file.path("assessment model/quahog/MIT.1 model/retrospective run/", "retrospectives", paste("retro", 0:-5, sep = "")))

retroSummary <- SSsummarize(retroModels)
endyrvec <- retroSummary[["endyrs"]] + 0:-5
Mohnrho_list <- SSmohnsrho(retroSummary, endyrvec, startyr = 1982, verbose = TRUE)

# generate a dataframe below in a format consistent with the ASAP output
Mohnrho_MIT.1_df <- data.frame(X = "Mohn.rho", f.rho = Mohnrho_list$F/5, ssb.rho = Mohnrho_list$SSB/5, recr.rho = Mohnrho_list$Rec/5, jan1b.rho = Mohnrho_list$Bratio,
                             explb.rho = NA, stockn.rho = NA, Age.1 = NA, Age.2 = NA, Age.3 = NA, Age.4 = NA, Age.5 = NA, Age.6 = NA, Age.7 = NA, Age.8 = NA,
                             Scenario = "MIT.1 Model")

write.csv(Mohnrho_MIT.1_df, "results/stock assessment/quahog/retrospective error/SS_retro_error_MIT.1.csv", row.names = FALSE)

remove(retroModels, retroSummary, endyrvec, Mohnrho_list, Mohnrho_WEE_df)



  ## 4.4 MIT.2 model ----

# do retro fitting 

retro(dir = "assessment model/quahog/MIT.2 model/retrospective run/",
      oldsubdir = "",
      newsubdir = "retrospectives",
      subdirstart = "retro",
      years = 0:-5)

# evaluate retro results

retroModels <- SSgetoutput(dirvec = file.path("assessment model/quahog/MIT.2 model/retrospective run/", "retrospectives", paste("retro", 0:-5, sep = "")))

retroSummary <- SSsummarize(retroModels)
endyrvec <- retroSummary[["endyrs"]] + 0:-5
Mohnrho_list <- SSmohnsrho(retroSummary, endyrvec, startyr = 1982, verbose = TRUE)

# generate a dataframe below in a format consistent with the ASAP output
Mohnrho_MIT.2_df <- data.frame(X = "Mohn.rho", f.rho = Mohnrho_list$F/5, ssb.rho = Mohnrho_list$SSB/5, recr.rho = Mohnrho_list$Rec/5, jan1b.rho = Mohnrho_list$Bratio,
                               explb.rho = NA, stockn.rho = NA, Age.1 = NA, Age.2 = NA, Age.3 = NA, Age.4 = NA, Age.5 = NA, Age.6 = NA, Age.7 = NA, Age.8 = NA,
                               Scenario = "MIT.2 Model")

write.csv(Mohnrho_MIT.2_df, "results/stock assessment/quahog/retrospective error/SS_retro_error_MIT.2.csv", row.names = FALSE)

remove(retroModels, retroSummary, endyrvec, Mohnrho_list, Mohnrho_WEE_df)



  ## 3.4 MIT.3 model ----

# do retro fitting 

retro(dir = "assessment model/quahog/MIT.3 model/retrospective run/",
      oldsubdir = "",
      newsubdir = "retrospectives",
      subdirstart = "retro",
      years = 0:-5)

# evaluate retro results

retroModels <- SSgetoutput(dirvec = file.path("assessment model/quahog/MIT.3 model/retrospective run/", "retrospectives", paste("retro", 0:-5, sep = "")))

retroSummary <- SSsummarize(retroModels)
endyrvec <- retroSummary[["endyrs"]] + 0:-5
Mohnrho_list <- SSmohnsrho(retroSummary, endyrvec, startyr = 1982, verbose = TRUE)

# generate a dataframe below in a format consistent with the ASAP output
Mohnrho_MIT.3_df <- data.frame(X = "Mohn.rho", f.rho = Mohnrho_list$F/5, ssb.rho = Mohnrho_list$SSB/5, recr.rho = Mohnrho_list$Rec/5, jan1b.rho = Mohnrho_list$Bratio,
                               explb.rho = NA, stockn.rho = NA, Age.1 = NA, Age.2 = NA, Age.3 = NA, Age.4 = NA, Age.5 = NA, Age.6 = NA, Age.7 = NA, Age.8 = NA,
                               Scenario = "MIT.3 Model")

write.csv(Mohnrho_MIT.3_df, "results/stock assessment/quahog/retrospective error/SS_retro_error_MIT.3.csv", row.names = FALSE)




    ## 4.6 plot ----

base_df <- read.csv("results/stock assessment/quahog/retrospective error/SS_retro_error_base.csv")
WEE_df <- read.csv("results/stock assessment/quahog/retrospective error/SS_retro_error_WEE.csv")
MIT.1_df <- read.csv("results/stock assessment/quahog/retrospective error/SS_retro_error_MIT.1.csv")

plot_df <- rbind(base_df, WEE_df, MIT.1_df) %>%
  filter(X == "Mohn.rho") %>%
  select(f.rho, ssb.rho, Scenario) %>%
  pivot_longer(cols = c(f.rho, ssb.rho), names_to = "Metric", values_to = "Value") %>%
  add_column(Species = "Quahog")

write.csv(plot_df, "results/stock assessment/quahog/retrospective error/retro_error_combined.csv", row.names = FALSE)

rho_SF_plot <- ggplot(plot_df, aes(x = Scenario)) +
  geom_bar(aes(y = Value), stat = "identity") +
  geom_hline(yintercept = 0) +
  facet_wrap(.~ Metric) +
  ylab("5-year Mohn's rho") 


remove(base_df, WEE_df, plot_df, rho_SF_plot)


















# ------------------------------------------------------------------------ #




