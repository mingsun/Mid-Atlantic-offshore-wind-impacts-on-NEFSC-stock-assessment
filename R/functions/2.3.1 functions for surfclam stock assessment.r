# pool.cv used to calcualte mean CVd ----

pool.cv <- function(x) {return(sqrt(mean(na.omit(x)^2)))}



# calculate Fthreshold ----

rho = 4.136303

GetFref = function(rlst, rhoF){
  
  #Find the F reference point
  yrs  <-  unique(rlst$timeseries$Yr)
  WF <- rlst$derived_quants[match(paste("F_", yrs, sep = ""), rlst$derived_quants$Label), 1:3]
  WF$Yr <-  yrs
  
  #The Fref is supposed to be based on the southern area
  FrefS <-  mean(rlst$timeseries$`F:_1`[rlst$timeseries$Yr%in%c(1982:2015) & rlst$timeseries$Area==1])
  FrefW <-  FrefS * rhoF
  WF$CV <-  WF$StdDev/WF$Value
  FMSY <-  data.frame("Value" = FrefW, "CV" = pool.cv(WF$CV[WF$Yr%in%(1982:2015)]))
  return(FMSY)
}



