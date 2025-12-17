ExtractSpecifiedIndexAgeCompResids <- function(asap, index.names, od) {
  
  ages <- seq(1, asap$parms$nages)
  nages <- length(ages)
  years <- seq(asap$parms$styr, asap$parms$endyr)  # FULL RANGE DEFINED ONCE
  nyrs <- asap$parms$nyears
  
  all_residuals <- list()
  
  for (i in 1:asap$parms$nindices) {
    
    this_index_name <- index.names[i]
    
    # Only process if index name is not NA
    if (!is.na(this_index_name)) {
      
      acomp.obs1 <- as.data.frame(asap$index.comp.mats[[2*i - 1]])
      acomp.obs <- acomp.obs1
      acomp.pred <- as.data.frame(asap$index.comp.mats[[2*i]])
      
      index.yrs <- which(asap$index.Neff.init[i, ] > 0)  # POSITIONS with positive Neff
      
      if (length(index.yrs) > 0) {
        
        s.age <- asap$control.parms$index.sel.start.age[i]
        e.age <- asap$control.parms$index.sel.end.age[i]
        
        for (j in 1:length(index.yrs)) {
          
          year_idx <- index.yrs[j]            # POSITION
          year_val <- years[year_idx]          # ACTUAL year value from full year vector
          
          # Raw residuals: obs - pred
          resids <- as.numeric(acomp.obs[year_idx, ] - acomp.pred[year_idx, ])
          
          # Standard deviation
          tmp.sd <- as.numeric(sqrt(acomp.pred[year_idx, ] * (1 - acomp.pred[year_idx, ]) / asap$index.Neff.init[i, year_idx]))
          
          # Pearson standardized residuals
          pearson_resids <- resids / tmp.sd
          
          # Only store within selected age range
          ages_used <- s.age:e.age
          
          for (a in ages_used) {
            all_residuals[[length(all_residuals) + 1]] <- data.frame(
              Index = this_index_name,
              Year = year_val,
              Age = a,
              Residual = pearson_resids[a]
            )
          }
          
        } # end year loop
      } # end if positive Neff
    } # end if non-NA index
  } # end index loop
  
  # Combine into a dataframe
  residuals_df <- do.call(rbind, all_residuals)
  
  # Save to CSV
  output_file <- file.path(od, "Specified_Index_AgeComp_Pearson_Residuals.csv")
  write.csv(residuals_df, output_file, row.names = FALSE)
  
  return(residuals_df)
}



output_file <- file.path(paste0(od, "Index_AgeComp_Pearson_Residuals_", asap.name, ".csv"))

