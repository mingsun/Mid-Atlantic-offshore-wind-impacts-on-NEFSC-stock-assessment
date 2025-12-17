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

MIT_1_AI <- read.csv("results/indices for assessment/summer flounder/mitigation.1.indices.csv") %>%
  add_column(MODEL = "MIT.1") %>%
  rename(fit = STRATIFIED_MEAN_N, lwr = lo_CI_95, upr = up_CI_95) %>%
  select(YEAR, SEASON, MODEL, fit, lwr, upr)

AI_df <- rbind(original_AI, WEE_AI, MIT_1_AI)
remove(original_AI, WEE_AI, MIT_1_AI)


ggplot(AI_df, aes(x = YEAR)) +
  # geom_ribbon(data = random_abundance_df, aes(ymin = QUANT_0.025, ymax = QUANT_0.975), fill = "grey", alpha = 0.3) + # first plot the ribbon for the random estimate
  geom_errorbar(aes(ymin = lwr, ymax = upr, color = MODEL), width = 0.2, position = position_dodge(width = 0.5)) +
  geom_point(aes(y = fit, color = MODEL), position = position_dodge(width = 0.5)) +
  geom_line(aes(y = fit, color = MODEL)) +
  facet_wrap(.~SEASON, nrow = 2) +
  # scale_color_manual(values = c("indianred2", "steelblue3"), labels = c("impacted", "original"), name = "Time Series") +
  ylab("biomass indices [g/tow]") +
  theme_minimal() +
  theme(legend.position = c(0.05, 0.95), # Place the legend in the top-left corner
        legend.justification = c(0, 1))

#  ------------------------------------------------------------------------------------------------ #






# 2. mitigation 3 ----


# in this section, I applied the ratio between the original_AI and impacted_AI to the original NEFSC indices used in assessment

  ## 2.1 extract the ratio ----

load("results/indices for assessment/summer flounder/mitigation_3_model-based indices/MIT.3.ALB_ratio.Rdata")
load("results/indices for assessment/summer flounder/mitigation_3_model-based indices/MIT.3.BIG_ratio.Rdata")



  ## 2.2 extract the abundance indices from the ASAP model ----

# ASAP assessment input value 

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



    ## 2.3 apply the ratio to estimate Wind Energy Excluded indices (MIT.3) ----

# note: be careful with the round and -999 values, need to be very consistent with the original indices format 

ALB_spring_MIT.3 <- round(ALB_spring * ALB_MIT.3_SPRING_ratio, 2)
ALB_fall_MIT.3 <-   round(ALB_fall *   ALB_MIT.3_FALL_ratio, 2)

BIG_spring_MIT.3 <- round(BIG_spring * BIG_MIT.3_SPRING_ratio)
BIG_spring_MIT.3[c(12),] <- -999

BIG_fall_MIT.3 <-   round(BIG_fall *   BIG_MIT.3_FALL_ratio)
BIG_fall_MIT.3[c(9, 12),] <- -999


#  ------------------------------------------------------------------------------------------------ #


  # 3. generate ASAP input with WEE abundance indices  -----------------------------------------------------------------------

# ALB spring (1982-2008)
data$dat$IAA_mats[[2]][1:27,c(2,4:11)] <- ALB_spring_MIT.3

# ALB fall (1982-2008)
data$dat$IAA_mats[[3]][1:27,c(2,4:11)] <- ALB_fall_MIT.3

# BIG spring (2009-2022)
data$dat$IAA_mats[[25]][28:41,c(2,4:11)] <- BIG_spring_MIT.3

# BIG fall (2009-2022)
data$dat$IAA_mats[[26]][28:41,c(2,4:11)] <- BIG_fall_MIT.3

WriteASAP3DatFile("assessment model/summer flounder/ASAP input data/ASAP3_MTA2023_MIT.3_BTS.DAT", data, 
                  header.text = "Summer Flounder 2023: BTS indices using MIT.3 dataset (Sun 2024) \n \n \n ") 



#  ------------------------------------------------------------------------------------------------ #




# 4. extract assessment output  -----------------------------------------------------------------------


# get the variable needed


wd <- "assessment model/summer flounder/ASAP/MIT.3 model/"
asap.name <- "ASAP3_MTA2023_MIT.3_BTS"
asap <- dget("assessment model/summer flounder/ASAP/MIT.3 model/ASAP3_MTA2023_MIT.3_BTS.RDAT") 
gn <- GrabNames(wd, asap.name, asap) # names of the catch fleet and survey
fleet.names <- gn$fleet.names
index.names <- gn$index.names

# reads in auxiliary files with values from the *.std, *.par, and *.cor
a1 <- GrabAuxFiles(wd, asap.name, asap = asap, fleet.names, index.names)

# generate the estimate output, SSB and F
SummarizeASAP(asap, a1, od = "results/stock assessment/summer flounder/")

output <- read.csv("results/stock assessment/summer flounder/ASAP_summary_ASAP3_MTA2023_MIT.3_BTS.csv") %>%
  add_column(Scenario = "MIT.3 Model") %>%
  write.csv("results/stock assessment/summer flounder/ASAP_summary_ASAP3_MTA2023_MIT.3_BTS.csv", row.names = FALSE)




#  ------------------------------------------------------------------------------------------------ #