df <- subset.data.frame(Input, Input$SAMPLE_ID == 'WEHI-1251426-001')
load("debug.RData")
library(plyr)
library(drc)

#Inputs
fixedUpperOn = FALSE
fixedUpperVal = c('Upper Limit:(Intercept)'=110)
fixedLowerOn = FALSE
fixedLowerVal = c('Lower Limit:(Intercept)'=-10)
fixedSlopeOn = FALSE
fixedSlopeVal = c('Slope Limit:(Intercept)'='')
minLower = -Inf
maxLower = 10
minUpper = 50
maxUpper = Inf
maxInflection = Inf
minInflection = 0
slopeLimits = '-Inf,0'
fitForEach = 'replicate'
EDLevel = 50
EDType
linSlopeConfLimit = 2

#Inputs
Input
fixedUpperVal
fixedLowerVal
fixedSlopeVal
minLower
maxLower
minUpper
maxUpper
maxInflection
minInflection
slopeLimits
fitForEach
EDLevel
EDType
linSlopeConfLimit

#Convert inputs
fixedLower <- as.numeric(fixedLowerVal)
fixedLowerOn <- if (is.na(fixedLower)) FALSE else TRUE
fixedUpper <- as.numeric(fixedUpperVal)
fixedUpperOn <- if (is.na(fixedUpper)) FALSE else TRUE
fixedSlope <- as.numeric(fixedSlopeVal)
fixedSlopeOn <- if (is.na(fixedSlope)) FALSE else TRUE
minLower <- as.numeric(minLower)
maxLower <- as.numeric(maxLower)
minUpper <- as.numeric(minUpper)
maxUpper <- as.numeric(maxUpper)
maxInflection  <- as.numeric(maxInflection)
minInflection <- as.numeric(minInflection)
maxSlope <- as.numeric(substring(slopeLimits,regexpr(',',slopeLimits)[1]+1))
minSlope <- as.numeric(substr(slopeLimits,1,regexpr(',',slopeLimits)[1]-1))

#Return fits for each sample, or each sample on each plate
groupBy <- if (fitForEach == 'Sample') ~ SAMPLE_ID + EXPERIMENT_ID + PROTOCOL_ID else ~ SAMPLE_ID + EXPERIMENT_ID + SAMPLE_PLATE_ID + PROTOCOL_ID

#Remove excluded points from calculation
sampleData <- subset.data.frame(Input, (Input$KNOCKOUT != "Y" | is.na(Input$KNOCKOUT)))
#Remove anything not marked as a sample
sampleData <- subset.data.frame(sampleData, sampleData$SAMPTYPE == "S")

#Seperate ctrl wells
highCtrls <- subset.data.frame(Input, Input$SAMPTYPE == 'H')
lowCtrls <- subset.data.frame(Input, Input$SAMPTYPE == 'L')
AverageHighCtrl <- mean(highCtrls$Normalized)
AverageLowCtrl <- mean(lowCtrls$Normalized)

#Change fixed values to NA if they're not being used
if(fixedUpperOn == FALSE) fixedUpperVal = NA 
if(fixedLowerOn == FALSE) fixedLowerVal = NA 
if(fixedSlopeOn == FALSE) fixedSlopeVal = NA              

#Set limits for fitting, skip limit for any fixed values
lowerlimit <- c(if (fixedSlopeOn) NULL else minSlope, #Slope
                if (fixedLowerOn) NULL else minLower, #Lower
                if (fixedUpperOn) NULL else minUpper, #Upper
                minInflection) # Inflection point
upperlimit <- c(if (fixedSlopeOn) NULL else maxSlope, #Slope
                if (fixedLowerOn) NULL else maxLower, #Lower
                if (fixedUpperOn) NULL else maxUpper, #Upper
                maxInflection) # Inflection point

#Function to perform fit and reformat parameters
fit4pl <- function (df) {
  #Store per-group parameters
  print(df$SAMPLE_ID[1])
  groups <- c(UNIQUE_PROP_ID = paste(unique(df$UNIQUE_PROP_ID), collapse = ", "))
  groups <- c(groups, SAMPLE_PLATE_ID = paste(unique(df$SAMPLE_PLATE_ID), collapse = ", "))
  groups <- c(groups, "Max Conc" = max(df$CONC), "Min Conc" = min(df$CONC))
  
  #Calculate an Log(x) ~ Normalized linear fit
  df$Log10CONC <- log10(df$CONC)
  lin_fit<-nls(Normalized ~ yIntercept + slope * Log10CONC,
           data = df,
           start = list(slope= 0, yIntercept = 0))
  lin_params <- summary(lin_fit)$parameters
  lin_flat <- c(lin_params)
  names(lin_flat) <- c(outer(paste("LinFit", rownames(lin_params)), colnames(lin_params), paste, sep = ":"))
 
  # Calculate a logistic fit
  tryCatch({
    fit <- drm(df$Normalized~df$CONC, 
               data = sampleData, 
               fct = LL.4(fixed=c(fixedSlope, fixedLower, fixedUpper, NA), 
                          names = c("Slope", "Lower Limit", "Upper Limit", "ed")), 
               control = drmc(errorm = FALSE, warnVal = -1, noMessage = TRUE),
               lowerl = lowerlimit,
               upperl = upperlimit
              )
    #reformat paramters so all types of fits are consistent, and add ED## calculation
    parms <- c(fit$parmMat)
    if (length(parms) > 0) {
      names(parms) <- unlist(fit$parNames[1])
      if (fixedSlopeOn) parms <- append(parms, fixedSlopeVal,0)
      if (fixedLowerOn) parms <- append(parms, fixedLowerVal,1)
      if (fixedUpperOn) parms <- append(parms, fixedUpperVal,2)
      ed <- ED(fit,EDLevel, interval = "delta", type = EDType)
      edflat <- c(ed)
      names(edflat) <- c(outer("ED", colnames(ed), paste, sep = ":"))
      parms <- c(parms, edflat)

    } else {
      #If fit fails to converge, insert NaNs in place of fit parameters
      parms <- rep(NaN, 8)
    }
    
  #Calculate the ED alpha, check if it's too big, small or not there
  EDLowLimit <- max(as.numeric(groups['Min Conc']), minInflection)
  EDHighLimit <- min(as.numeric(groups['Max Conc']), maxInflection)
  #if it's not there or doesn't make sense
  #check if the values are closer to the high ctrl or low ctrl
  if (is.na(parms[5]) | 
      (slopeLimits == "-Inf, 0" & lin_params['slope','Estimate'] < 0) |
      (slopeLimits == "0, Inf" & lin_params['slope','Estimate'] > 0) |
      (lin_params['slope','t value'] < linSlopeConfLimit)) {
    meanNormalized <- mean(df$Normalized)
    if ((meanNormalized - AverageHighCtrl) > (AverageLowCtrl - meanNormalized)) {
      ED_Alpha <- paste("<", EDLowLimit)
    } else {
      ED_Alpha <- paste(">", EDHighLimit)
    }
  } 
  # If the ED is outside the conc or fit limits flag as >= or <=
  else if (parms[5] >= EDHighLimit) 
    {ED_Alpha <- paste(">=", EDHighLimit)}
  else if (parms[5] <= EDLowLimit) 
  {ED_Alpha <- paste("<=", EDLowLimit)}
  # Otherwise report the ED value as is
  else {ED_Alpha <- as.character(parms[5])}
  })
  
  highestPoint <- max(df$Normalized)
  lowestPoint <- min(df$Normalized)
  EDPoint <- ((EDLevel / 100) * (highestPoint - lowestPoint)) + lowestPoint
  pointsAbove <- subset.data.frame(df, df$Normalized > EDPoint)
  pointsBelow <- subset.data.frame(df, df$Normalized < EDPoint)
  pointAbove <- subset.data.frame(df, df$Normalized == min(pointsAbove$Normalized))
  pointBelow <- subset.data.frame(df, df$Normalized == min(pointsBelow$Normalized))
  midConc <- mean(mean(pointBelow$CONC), mean(pointAbove$CONC))
  
  return(c(groups,parms,ED_Alpha, lin_flat, "midConc" = midConc))
}

# Group data by SAMPLE_ID and apply fit function to each group
out <- ddply(sampleData, groupBy, function(t) fit4pl(t))


names(out) <- c("SAMPLE_ID",
                "EXPERIMENT_ID",
                "PROTOCOL_ID",
                "UNIQUE_PROP_ID",
                "SAMPLE_PLATE_ID",
                "Max Conc",
                "Min Conc",
                "Slope",
                "Lower Limit",
                "Upper Limit",
                "Inflection Point",
                "ED Estimate",
                "ED Std Error",
                "ED Lower CI",
                "ED Upper CI",
                "ED Alpha",
                "LinFit slope:Estimate",
                "LinFit yIntercept:Estimate",
                "LinFit slope:Std. Error",    
                "LinFit yIntercept:Std. Error",
                "LinFit slope:t value",       
                "LinFit yIntercept:t value",
                "LinFit slope:Pr(>|t|)",  
                "LinFit yIntercept:Pr(>|t|)",
                "Simple ED Conc") 
out <- cbind(out,
            "fixedLower" = fixedLower,
            "fixedUpper" = fixedUpper,
            "fixedSlope" = fixedSlope,
            "minLower" = minLower,
            "maxLower" = maxLower,
            "minUpper" = minUpper,
            "maxUpper" = maxUpper,
            "maxInflection" = maxInflection,
            "minInflection" = minInflection,
            "maxSlope" = maxSlope,
            "minSlope" = minSlope,
            "EDLevel" = EDLevel,
            "EDType" = EDType)
           
#Convert specified columns to numeric
char_cols = c(1,3,4,5,16,38);    
out[,-char_cols] = apply(out[,-char_cols], 2, function(x) as.numeric(as.character(x)))
                  
#Function to round all numeric columns in a data frame
round_df <- function(df, digits) {
  nums <- vapply(df, is.numeric, FUN.VALUE = logical(1))
  
  df[,nums] <- round(df[,nums], digits = digits)
  
  (df)
}
out <- round_df(out,30)
