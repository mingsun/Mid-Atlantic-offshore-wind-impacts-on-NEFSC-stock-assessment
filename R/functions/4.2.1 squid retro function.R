# peel regression analyses for squid "retrospective error"
# Difference between the predicted value from the peel regression and the full regression at the peeled year, 
# standardized by the full regression prediction:
# (predicted(peel) - predicted(full)) / predicted(full)



peel_regression_analysis <- function(data, response = "Annualized.Exploitation.Indices.final",
                                     predictor = "Annual.Catch", year_col = "Year", max_peel = 5) {
  # Clean and sort
  data <- data[!is.na(data[[response]]), ]
  data <- data[order(data[[year_col]]), ]
  
  # Unique valid years
  unique_years <- sort(unique(data[[year_col]]))
  terminal_year <- max(unique_years)
  
  # Last max_peel years BEFORE terminal year
  comparison_years <- tail(unique_years[unique_years < terminal_year], max_peel)
  
  # Full model and prediction
  full_model <- lm(as.formula(paste(response, "~", predictor)), data = data)
  full_pred_df <- data.frame(
    Year = data[[year_col]],
    FullPred = predict(full_model, newdata = data)
  )
  
  # Peel results
  results <- data.frame(
    Peel = seq_len(length(comparison_years)),
    FinalYear = comparison_years,
    Slope = NA,
    Intercept = NA,
    R2 = NA,
    N = NA,
    RelativeDeviation = NA
  )
  
  for (i in seq_along(comparison_years)) {
    peel_to_year <- comparison_years[i]
    df_peel <- data[data[[year_col]] < peel_to_year, ]
    
    if (nrow(df_peel) < 3) next
    
    fit <- lm(as.formula(paste(response, "~", predictor)), data = df_peel)
    
    final_row <- data[data[[year_col]] == peel_to_year, ]
    pred_peel <- predict(fit, newdata = final_row)
    pred_full <- full_pred_df$FullPred[full_pred_df$Year == peel_to_year]
    
    results$Slope[i]     <- coef(fit)[2]
    results$Intercept[i] <- coef(fit)[1]
    results$R2[i]        <- summary(fit)$r.squared
    results$N[i]         <- nrow(df_peel)
    results$RelativeDeviation[i] <- (pred_peel - pred_full) / pred_full
  }
  
  # Add full-model terminal year reference row
  ref_row <- data.frame(
    Peel = NA,
    FinalYear = terminal_year,
    Slope = coef(full_model)[2],
    Intercept = coef(full_model)[1],
    R2 = summary(full_model)$r.squared,
    N = nrow(data),
    RelativeDeviation = 0
  )
  
  # Combine and fix row names
  results <- rbind(results, ref_row)
  row.names(results) <- NULL  
  
  return(results)
}
