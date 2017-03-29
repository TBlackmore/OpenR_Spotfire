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
AverageHighCtrl <- mean(highCtrls$RESPONSE)
AverageLowCtrl <- mean(lowCtrls$RESPONSE)

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
  print(df$SAMPLE_ID[1])
  groups <- c(UNIQUE_PROP_ID = paste(unique(df$UNIQUE_PROP_ID), collapse = ", "))
  groups <- c(groups, SAMPLE_PLATE_ID = paste(unique(df$SAMPLE_PLATE_ID), collapse = ", "))
  groups <- c(groups, "Max Conc" = max(df$CONC), "Min Conc" = min(df$CONC))

  tryCatch({
    fit <- drm(df$RESPONSE~df$CONC, 
               data = sampleData, 
               fct = LL.4(fixed=c(fixedSlope, fixedLower, fixedUpper, NA), 
                          names = c("Slope", "Lower Limit", "Upper Limit", "ed")), 
               control = drmc(errorm = FALSE, warnVal = -1, noMessage = TRUE),
               lowerl = lowerlimit,
               upperl = upperlimit
              )
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
      parms <- rep(NaN, 8)
    }
  #Calculate the ED alpha, check if it's too big, small or not there
  if (is.na(parms[5])) {
    #if it's not there check if the values are closer to the high ctrl or low ctrl
    meanResponse <- mean(df$RESPONSE)
    if ((meanResponse - AverageHighCtrl) > (AverageLowCtrl - meanResponse)) {
      ED_Alpha <- paste(">", groups['Max Conc'])
    } else {
      ED_Alpha <- paste("<", groups['Min Conc'])
    }
  } 
  else if (parms[5] >= maxInflection) 
    {ED_Alpha <- paste(">", groups['Max Conc'])}
  else if (parms[5] <= minInflection) 
    {ED_Alpha <- paste("<", groups['Min Conc'])}
  })
  
  return(c(groups,parms,ED_Alpha))
}

# Group data by SAMPLE_ID and apply fit function to each group
out <- ddply(sampleData, groupBy, function(t) fit4pl(t))


names(out) <- c("SAMPLE_ID", "EXPERIMENT_ID", "PROTOCOL_ID", "UNIQUE_PROP_ID", "SAMPLE_PLATE_ID", "Max Conc", "Min Conc", "Slope", "Lower Limit", "Upper Limit", "Inflection Point", "ED Estimate", "ED Std Error", "ED Lower CI", "ED Upper CI", "ED Alpha") 
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
           
#Replace NaN with NA
#out[is.na(out)] <- NA 

#Convert specified columns to numeric
cols = c(2, 6, 7,8,9,10,11,12,13,14,15);    
out[,cols] = apply(out[,cols], 2, function(x) as.numeric(as.character(x)))
                  
#Function to round all numeric columns in a data frame
round_df <- function(df, digits) {
  nums <- vapply(df, is.numeric, FUN.VALUE = logical(1))
  
  df[,nums] <- round(df[,nums], digits = digits)
  
  (df)
}
out <- round_df(out,30)
