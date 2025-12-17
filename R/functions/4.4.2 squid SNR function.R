# this function is used to calculate Signal-to-Noise Ratio (SNR),  

calculate_snr <- function(df, index_col, year_col = "Year", span = 0.5, method = "variance") {
  # Ensure index_col and year_col exist
  if (!all(c(index_col, year_col) %in% names(df))) {
    stop("Specified column names not found in the data frame.")
  }
  
  # Remove NA values
  temp <- df[!is.na(df[[index_col]]), c(year_col, index_col)]
  colnames(temp) <- c("Year", "Index")
  
  # Fit LOESS smoother
  loess_fit <- loess(Index ~ Year, data = temp, span = span)
  temp$Smoothed <- predict(loess_fit)
  temp$Residual <- temp$Index - temp$Smoothed
  
  # Compute SNR
  if (method == "variance") {
    snr <- var(temp$Smoothed, na.rm = TRUE) / var(temp$Residual, na.rm = TRUE)
  } else if (method == "mean") {
    snr <- mean(temp$Smoothed, na.rm = TRUE) / sd(temp$Residual, na.rm = TRUE)
  } else {
    stop("Method must be either 'variance' or 'mean'.")
  }
  
  return(snr)
}