# a function to calculate the total number of parameters 

get_asap_parameters <- function(asap) {
  nyears <- asap$parms$nyears
  nfleets <- asap$parms$nfleets
  nindices <- asap$parms$nindices
  nselparms <- asap$like.additional$nselparms
  nindexselparms <- asap$like.additional$nindexselparms
  
  phases <- asap$control.parms$phases
  
  df <- data.frame(
    Parameter = c("Selectivity", "Index Selectivity", 
                  "Recruitment Deviations", "Fmult Deviations", 
                  "Fmult Year 1", "Q Year 1", 
                  "Q Deviations", "N Year 1", 
                  "SR Steepness", "SR Scaler"),
    
    Phase = c(NA, NA, 
              phases$phase.recruit.devs, phases$phase.Fmult.devs, 
              phases$phase.Fmult.year1, phases$phase.q.year1, 
              phases$phase.q.devs, phases$phase.N.year1.devs, 
              phases$phase.steepness, phases$phase.SR.scaler),
    
    Included = c(TRUE, TRUE, 
                 phases$phase.recruit.devs >= 1, phases$phase.Fmult.devs >= 1, 
                 phases$phase.Fmult.year1 >= 1, phases$phase.q.year1 >= 1, 
                 phases$phase.q.devs >= 1, phases$phase.N.year1.devs >= 1, 
                 phases$phase.steepness >= 1, phases$phase.SR.scaler >= 1),
    
    Count = c(nselparms, nindexselparms,
              ifelse(phases$phase.recruit.devs >= 1, nyears, 0),
              ifelse(phases$phase.Fmult.devs >= 1, nyears * nfleets, 0),
              ifelse(phases$phase.Fmult.year1 >= 1, nfleets, 0),
              ifelse(phases$phase.q.year1 >= 1, nindices, 0),
              ifelse(phases$phase.q.devs >= 1, nyears * nindices, 0),
              ifelse(phases$phase.N.year1.devs >= 1, 1, 0),
              ifelse(phases$phase.steepness >= 1, 1, 0),
              ifelse(phases$phase.SR.scaler >= 1, 1, 0)
    ),
    stringsAsFactors = FALSE
  )
  
  return(df)
}
