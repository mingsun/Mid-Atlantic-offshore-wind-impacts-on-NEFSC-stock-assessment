library(tidyverse)
library(ASAPplots)
library(r4ss)
source("R/functions/4.4.3 ASAP get parameter function.r")

# let's compare the Negative Log-Likelihoods (NLLs), the lower the value, the better the fit 
# the higher the individual proportion, the larger effort or "penalty weighting" is used to fit that type of data

# we also compare AIC here, defined as AIC = 2k−2ln(L), where K is the number of estimated parameters in the model and ln(L) is the log-likelihood (NLL above)


# 1. summer flounder ----

  ## 1.1 base model ----

asap <- dget("assessment model/summer flounder/ASAP/base model/ASAP3_MTA2023_FINAL.rdat") 

Total <- asap$like$lk.total
Survey <- asap$like$lk.index.fit.total
Param.N <- sum(get_asap_parameters(asap)$Count)
AIC <- 2 * Param.N - 2 * (-Total)  # negate Total to get logLik

base_df <- data.frame(Scenario = "Base Model", Species = "Summer Flounder", 
                      Total.NLL = Total, Survey.NLL = Survey,
                      Param.N = Param.N, AIC = AIC)

remove(asap, Total, Survey, Param.N, AIC)



  ## 1.2 WEE_design-based_indices_model assessment ----

asap <- dget("assessment model/summer flounder/ASAP/WEE_design-based_indices_model/ASAP3_MTA2023_WEE_BTS.RDAT") 

Total <- asap$like$lk.total
Survey <- asap$like$lk.index.fit.total
Param.N <- sum(get_asap_parameters(asap)$Count)
AIC <- 2 * Param.N - 2 * (-Total)  # negate Total to get logLik

WEE_df <- data.frame(Scenario = "WEE Model", Species = "Summer Flounder", 
                     Total.NLL = Total, Survey.NLL = Survey,
                     Param.N = Param.N, AIC = AIC)

remove(asap, Total, Survey, Param.N, AIC)



  ## 1.3 MIT.1 assessment ----

asap <- dget("assessment model/summer flounder/ASAP/MIT.1 model/ASAP3_MTA2023_MIT.1_BTS.RDAT") 

Total <- asap$like$lk.total
Survey <- asap$like$lk.index.fit.total
Param.N <- sum(get_asap_parameters(asap)$Count)
AIC <- 2 * Param.N - 2 * (-Total)  # negate Total to get logLik

MIT.1_df <- data.frame(Scenario = "MIT.1 Model", Species = "Summer Flounder", 
                           Total.NLL = Total, Survey.NLL = Survey,
                           Param.N = Param.N, AIC = AIC)

remove(asap, Total, Survey, Param.N, AIC)



  ## 1.6 combine ----

SF_df <- rbind(base_df, WEE_df, MIT.1_df)
remove(base_df, WEE_df, MIT.1_df)

write.csv(SF_df, "results/stock assessment/summer flounder/model fit/NLL_AIC.csv", row.names = FALSE)

# ------------------------------------------------------------------------ #


# 2. squid ----

# For squid with index-based method, we cannot do likelihood component
# we calculate the SNR based on the loess fit

source("R/functions/4.4.2 squid SNR function.R")


base_df <- read.csv("results/stock assessment/squid/index-based summary_base.csv")
WEE_df <-  read.csv("results/stock assessment/squid/index-based summary_WEE.csv")

index_df <- rbind(base_df, WEE_df) %>%
  filter(Year >= 1987 & !is.na(Annualized.Exploitation.Indices.final))  %>%
  select(Scenario, Year, Annualized.Exploitation.Indices.final) %>%
  group_by(Scenario) %>%
  group_modify(~ {
    snr <- calculate_snr(df = .x, index_col = "Annualized.Exploitation.Indices.final", year_col = "Year",
      span = 0.5, method = "variance")
    tibble(SNR = snr)
  }) %>%  
  ungroup() %>%
  write.csv("results/stock assessment/squid/model fit/SNR.csv", row.names = FALSE)

remove(index_df, base_df, WEE_df)





# ------------------------------------------------------------------------ #



# 3. surfclam ----


  ## 3.1 base model ----

assessment_results <- SS_output("assessment model/surfclam/base model/", verbose = FALSE)

Total <- assessment_results$likelihoods_used$values[1]
Survey <- assessment_results$likelihoods_used$values[4]
Param.N <- sum(assessment_results$parameters$Active_Cnt > 0, na.rm = TRUE)
AIC <- 2 * Param.N - 2 * (-Total)  # negate Total to get logLik

base_df <- data.frame(Scenario = "Base Model", Species = "Surfclam", 
                      Total.NLL = Total, Survey.NLL = Survey,
                      Param.N = Param.N, AIC = AIC)

remove(assessment_results, Total, Survey, Param.N, AIC)


  ## 3.2 WEE_design-based_indices_model assessment ----

assessment_results <- SS_output("assessment model/surfclam/WEE_design-based_indices_model/", verbose = FALSE)

Total <- assessment_results$likelihoods_used$values[1]
Survey <- assessment_results$likelihoods_used$values[4]
Param.N <- sum(assessment_results$parameters$Active_Cnt > 0, na.rm = TRUE)
AIC <- 2 * Param.N - 2 * (-Total)  # negate Total to get logLik

WEE_df <- data.frame(Scenario = "WEE Model", Species = "Surfclam", 
                      Total.NLL = Total, Survey.NLL = Survey,
                      Param.N = Param.N, AIC = AIC)

remove(assessment_results, Total, Survey, Param.N, AIC)


  ## 3.3 MIT.1 assessment ----

assessment_results <- SS_output("assessment model/surfclam/MIT.1 model/", verbose = FALSE)

Total <- assessment_results$likelihoods_used$values[1]
Survey <- assessment_results$likelihoods_used$values[4]
Param.N <- sum(assessment_results$parameters$Active_Cnt > 0, na.rm = TRUE)
AIC <- 2 * Param.N - 2 * (-Total)  # negate Total to get logLik

MIT.1_df <- data.frame(Scenario = "MIT.1 Model", Species = "Surfclam", 
                     Total.NLL = Total, Survey.NLL = Survey,
                     Param.N = Param.N, AIC = AIC)

remove(assessment_results, Total, Survey, Param.N, AIC)



  ## 3.6 combine ----

AS_df <- rbind(base_df, WEE_df, MIT.1_df) 
remove(base_df, WEE_df, MIT.1_df)

write.csv(AS_df, "results/stock assessment/surfclam/model fit/NLL_AIC.csv", row.names = FALSE)


# ------------------------------------------------------------------------ #



# 4. quahog ----


  ## 4.1 base model ----

assessment_results <- SS_output("assessment model/quahog/base model/", verbose = FALSE)

Total <- assessment_results$likelihoods_used$values[1]
Survey <- assessment_results$likelihoods_used$values[4]
Param.N <- sum(assessment_results$parameters$Active_Cnt > 0, na.rm = TRUE)
AIC <- 2 * Param.N - 2 * (-Total)  # negate Total to get logLik

base_df <- data.frame(Scenario = "Base Model", Species = "Quahog", 
                      Total.NLL = Total, Survey.NLL = Survey,
                      Param.N = Param.N, AIC = AIC)

remove(assessment_results, Total, Survey, Param.N, AIC)


  ## 4.2 WEE_design-based_indices_model assessment ----

assessment_results <- SS_output("assessment model/quahog/WEE_design-based_indices_model/", verbose = FALSE)

Total <- assessment_results$likelihoods_used$values[1]
Survey <- assessment_results$likelihoods_used$values[4]
Param.N <- sum(assessment_results$parameters$Active_Cnt > 0, na.rm = TRUE)
AIC <- 2 * Param.N - 2 * (-Total)  # negate Total to get logLik

WEE_df <- data.frame(Scenario = "WEE Model", Species = "Quahog", 
                     Total.NLL = Total, Survey.NLL = Survey,
                     Param.N = Param.N, AIC = AIC)

remove(assessment_results, Total, Survey, Param.N, AIC)


  ## 4.3 MIT.1 model assessment ----

assessment_results <- SS_output("assessment model/quahog/MIT.1 model/", verbose = FALSE)

Total <- assessment_results$likelihoods_used$values[1]
Survey <- assessment_results$likelihoods_used$values[4]
Param.N <- sum(assessment_results$parameters$Active_Cnt > 0, na.rm = TRUE)
AIC <- 2 * Param.N - 2 * (-Total)  # negate Total to get logLik

MIT.1_df <- data.frame(Scenario = "MIT.1 Model", Species = "Quahog", 
                     Total.NLL = Total, Survey.NLL = Survey,
                     Param.N = Param.N, AIC = AIC)

remove(assessment_results, Total, Survey, Param.N, AIC)


  ## 4.6 combine ----

OQ_df <- rbind(base_df, WEE_df, MIT.1_df) 
remove(base_df, WEE_df, MIT.1_df)

write.csv(OQ_df, "results/stock assessment/quahog/model fit/NLL_AIC.csv", row.names = FALSE)


# ------------------------------------------------------------------------ #





# 5 plot ----

NLL_OQ_df <- rbind(NLL_base_df, NLL_WEE_df) %>%
  add_column(Species = "Quahog")  %>%
  # calculate other proportion
  mutate(lk.other = TOTAL - Survey) %>%
  # select(Scenario, Species, lk.catch.total, lk.index.fit.total, lk.other) %>%
  # wide to long
  pivot_longer(cols = c(lk.catch.total, lk.index.fit.total, lk.other), names_to = "Component", values_to = "Value") 


ggplot(NLL_OQ_df, aes(x = Scenario, y = values, fill = Item)) +
  geom_bar(stat = "identity") +  # stacked barplot
  geom_text(aes(label = round(Value, 1)),   # <--- Add values inside bars
            position = position_stack(vjust = 0.5),   # Centered within each block
            size = 3) +
  geom_text(aes(x = Scenario, y = lk.total, label = round(lk.total, 1)),
            vjust = -0.5,
            size = 3.5) +  # Total label above bar
  facet_wrap(~Species) +
  ylab("Negative Log-Likelihood (NLL)") +
  theme_minimal() +
  ggtitle("Likelihood Components by Scenario") +
  scale_fill_brewer(palette = "Set2") +
  theme(axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5))





NLL_SF_df <- rbind( NLL_base_df, NLL_WEE_df) %>%
  add_column(Species = "Summer Flounder")  %>%
  # calculate other proportion
  mutate(lk.other = lk.total - lk.catch.total - lk.index.fit.total) %>%
  # select(Scenario, Species, lk.catch.total, lk.index.fit.total, lk.other) %>%
  # wide to long
  pivot_longer(cols = c(lk.catch.total, lk.index.fit.total, lk.other), names_to = "Component", values_to = "Value") 


ggplot(NLL_SF_df, aes(x = Scenario, y = Value, fill = Component)) +
  geom_bar(stat = "identity") +  # stacked barplot
  geom_text(aes(label = round(Value, 1)),   # <--- Add values inside bars
            position = position_stack(vjust = 0.5),   # Centered within each block
            size = 3) +
  geom_text(aes(x = Scenario, y = lk.total, label = round(lk.total, 1)),
            vjust = -0.5,
            size = 3.5) +  # Total label above bar
  facet_wrap(~Species) +
  ylab("Negative Log-Likelihood (NLL)") +
  theme_minimal() +
  ggtitle("Likelihood Components by Scenario") +
  scale_fill_brewer(palette = "Set2") +
  theme(axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5))



remove(asap, NLL_base_df, NLL_WEE_df)

