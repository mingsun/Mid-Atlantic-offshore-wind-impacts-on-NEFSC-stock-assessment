library(tidyverse)
library(ASAPplots)


# 1. review abundance indices trend  ---------------------------------------------------------------

original_AI <- read.csv("results/indices for assessment/summer flounder/original.indices.csv") %>%
  add_column(MODEL = "original") %>%
  rename(fit = STRATIFIED_MEAN_N, lwr = lo_CI_95, upr = up_CI_95) %>%
  select(YEAR, SEASON, MODEL, fit, lwr, upr)


WEE_AI <- read.csv("results/indices for assessment/summer flounder/WEE.indices.csv") %>%
  add_column(MODEL = "WEE") %>%
  rename(fit = STRATIFIED_MEAN_N, lwr = lo_CI_95, upr = up_CI_95) %>%
  select(YEAR, SEASON, MODEL, fit, lwr, upr)

AI_df <- rbind(original_AI, WEE_AI) 


design_based_p <- ggplot(AI_df, aes(x = YEAR)) +
  # geom_ribbon(data = random_abundance_df, aes(ymin = QUANT_0.025, ymax = QUANT_0.975), fill = "grey", alpha = 0.3) + # first plot the ribbon for the random estimate
  geom_errorbar(aes(ymin = lwr, ymax = upr, color = MODEL), width = 0.2, position = position_dodge(width = 0.5)) +
  geom_point(aes(y = fit, color = MODEL), position = position_dodge(width = 0.5)) +
  geom_line(aes(y = fit, color = MODEL)) +
  facet_wrap(.~SEASON, nrow = 2) +
  scale_color_manual(values = c("indianred2", "steelblue3"), labels = c("impacted", "original"), name = "Time Series") +
  ylab("biomass indices [g/tow]") +
  theme_minimal() +
  theme(legend.position = c(0.05, 0.95), # Place the legend in the top-left corner
        legend.justification = c(0, 1))

remove(original_AI, WEE_AI, design_based_p)

#  ------------------------------------------------------------------------------------------------ #


# 2. total abundance loss  ---------------------------------------------------------------

# in this section, I applied the ratio between the original_AI and impacted_AI to the original NEFSC indices used in assessment

  ## 2.1 extract the ratio ----

load("results/indices for assessment/summer flounder/WEE/WEE.ALB_ratio.Rdata")
load("results/indices for assessment/summer flounder/WEE/WEE.BIG_ratio.Rdata")



  ## 2.2 extract the abundance indices from the ASAP model ----

# ASAP assessment input value 

data <- ReadASAP3DatFile("assessment model/summer flounder/ASAP input data/ASAP3_MTA2023_FINAL.DAT")
data$survey.names

# a test
# WriteASAP3DatFile("assessment model/summer flounder/ASAP input data/ASAP3_MTA2023_test_output.DAT", data, 
#                   header.text = "Summer Flounder 2023: test output \n# \n# \n#") 
# test <- ReadASAP3DatFile("assessment model/summer flounder/ASAP input data/ASAP3_MTA2023_test_output.DAT")

# 2 is BTS ALB spring, 3 is BTS ALB fall, they are 1982-2008
# 25 is BTS BIG spring, 26 is BTS BIG fall, they are 2009-2022
# col.1 is year 1982-2022, 1982-2008 correspond to row 1-27, 2009-2022 correspond to row 28-41, 
# col.2 is sum value, col.3 is CV
# column 4-11 correspond to age 1-8,  col. 11 is sample size

ALB_spring <- data$dat$IAA_mats[[2]][1:27,c(2,4:11)]
ALB_fall <- data$dat$IAA_mats[[3]][1:27,c(2,4:11)]

BIG_spring <- data$dat$IAA_mats[[25]][28:41,c(2,4:11)]
BIG_fall <- data$dat$IAA_mats[[26]][28:41,c(2,4:11)]


  ## 2.3 apply the ratio to estimate Wind Energy Excluded indices (WEE) ----

# note: be careful with the round and -999 values, need to be very consistent with the original indices format 

ALB_spring_WEE <- round(ALB_spring * ALB_WEE_SPRING_ratio, 2)
ALB_fall_WEE <-   round(ALB_fall *   ALB_WEE_FALL_ratio, 2)

BIG_spring_WEE <- round(BIG_spring * BIG_WEE_SPRING_ratio)
BIG_fall_WEE <-   round(BIG_fall *   BIG_WEE_FALL_ratio)


#  ------------------------------------------------------------------------------------------------ #




# 3. generate ASAP input with WEE abundance indices  -----------------------------------------------------------------------

# ALB spring (1982-2008)
data$dat$IAA_mats[[2]][1:27,c(2,4:11)] <- ALB_spring_WEE

# ALB fall (1982-2008)
data$dat$IAA_mats[[3]][1:27,c(2,4:11)] <- ALB_fall_WEE

# BIG spring (2009-2022)
data$dat$IAA_mats[[25]][28:41,c(2,4:11)] <- BIG_spring_WEE

# BIG fall (2009-2022)
data$dat$IAA_mats[[26]][28:41,c(2,4:11)] <- BIG_fall_WEE

WriteASAP3DatFile("assessment model/summer flounder/ASAP input data/ASAP3_MTA2023_WEE_BTS.DAT", data, 
                  header.text = "Summer Flounder 2023: BTS indices using WEE dataset (Sun 2024) \n \n \n ") 



#  ------------------------------------------------------------------------------------------------ #




# 4. extract assessment output  -----------------------------------------------------------------------


  ## 4.1 based model assessment ----

  # get the variable needed

wd <- "assessment model/summer flounder/ASAP/base model/"
asap.name <- "ASAP3_MTA2023_FINAL"
asap <- dget("assessment model/summer flounder/ASAP/base model/ASAP3_MTA2023_FINAL.rdat") 
gn <- GrabNames(wd, asap.name, asap) # names of the catch fleet and survey
fleet.names <- gn$fleet.names
index.names <- gn$index.names

  # reads in auxiliary files with values from the *.std, *.par, and *.cor
a1 <- GrabAuxFiles(wd, asap.name, asap = asap, fleet.names, index.names)

  # generate the estimate output, SSB and F
SummarizeASAP(asap, a1, od = "results/stock assessment/summer flounder/")

output <- read.csv("results/stock assessment/summer flounder/ASAP_summary_ASAP3_MTA2023_FINAL.csv") %>%
  add_column(Scenario = "Base Model") %>%
  write.csv("results/stock assessment/summer flounder/ASAP_summary_ASAP3_MTA2023_FINAL.csv", row.names = FALSE)





  ## 4.2 WEE_design-based_indices_model assessment ----

  # get the variable needed

wd.WEE <- "assessment model/summer flounder/ASAP/WEE_design-based_indices_model/"
asap.name.WEE <- "ASAP3_MTA2023_WEE_BTS"
asap.WEE <- dget("assessment model/summer flounder/ASAP/WEE_design-based_indices_model/ASAP3_MTA2023_WEE_BTS.RDAT") 
gn.WEE <- GrabNames(wd.WEE, asap.name.WEE, asap.WEE) # names of the catch fleet and survey
fleet.names.WEE <- gn.WEE$fleet.names
index.names.WEE <- gn.WEE$index.names

# reads in auxiliary files with values from the *.std, *.par, and *.cor
a1.WEE <- GrabAuxFiles(wd.WEE, asap.name.WEE, asap = asap.WEE, fleet.names.WEE, index.names.WEE)

# generate the output 
SummarizeASAP(asap.WEE, a1.WEE, od = "results/stock assessment/summer flounder/")

output <- read.csv("results/stock assessment/summer flounder/ASAP_summary_ASAP3_MTA2023_WEE_BTS.csv") %>%
  add_column(Scenario = "WEE Model") %>%
  write.csv("results/stock assessment/summer flounder/ASAP_summary_ASAP3_MTA2023_WEE_BTS.csv", row.names = FALSE)




#  ------------------------------------------------------------------------------------------------ #










