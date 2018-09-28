FindSS <- function(vec){
  for (i in 1:(length(vec) - 10)) {
    #From start, search for 3 points each 5 unit time apart where delta < 10^-3
    thisSlope <- abs(vec[i+1] - vec[i])
    nextSlope <- abs(vec[i+5] - vec[i+4])
    lastSlope <- abs(vec[i+10] - vec[i+9])
    
    if (thisSlope < 10^-3 & nextSlope < 10^-3 & lastSlope < 10^-3) {
      #Verification: make sure no change after plateau: extrapolate
      #Extrapolate by thisSlope and calculate value at end of integration time
      ExtrEndTimeValue = vec[i] + lastSlope*(length(vec) - i)
      if ((ExtrEndTimeValue < (vec[length(vec)]) + 0.05)& (ExtrEndTimeValue > (vec[length(vec)]) - 0.05)) {
        return(i)
      }
      #There is a plateau in the middle, but in the end density changes again: find next plateau and verify agin
      
    }
  }
  #Til the end no steady state found
  return(NaN)
}