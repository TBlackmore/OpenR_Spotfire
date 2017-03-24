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
fitForEach = 'Sample'
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

#Calculate fit limits based on filtered data
#minUpper <- -Inf
#maxUpper <- max(sampleData$RESPONSE) * 2
#minLower <- min(sampleData$RESPONSE) - (maxUpper /2)
#maxLower <- Inf

if(fixedUpperOn == FALSE) fixedUpperVal = NA 
if(fixedLowerOn == FALSE) fixedLowerVal = NA 
if(fixedSlopeOn == FALSE) fixedSlopeVal = NA              

#Put some limit on inflection point max to speed calculations
#maxInflection <- max(sampleData$CONC) * 10

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
  fit <- drm(df$RESPONSE~df$CONC, 
      data = sampleData, 
      fct = LL.4(fixed=c(fixedSlopeVal, fixedLowerVal, fixedUpperVal, NA), 
            names = c("Slope", "Lower Limit", "Upper Limit", "ed")), 
      control = drmc(errorm = FALSE, warnVal = -1, noMessage = TRUE),
      lowerl = lowerlimit,
      upperl = upperlimit
      )
  parms <- c(fit$parmMat)
  names(parms) <- unlist(fit$parNames[1])
  #add fixed parameters to parms to keep number of columns consistent
  if (fixedSlopeOn) parms <- append(parms, fixedSlopeVal,0)
  if (fixedLowerOn) parms <- append(parms, fixedLowerVal,1)
  if (fixedUpperOn) parms <- append(parms, fixedUpperVal,2)
  ed <- ED(fit,EDLevel, interval = "delta", type = EDType)
  edflat <- c(ed)
  names(edflat) <- c(outer("ED", colnames(ed), paste, sep = ":"))
  parms <- c(parms, edflat)
  return(parms)
}


# Group data by SAMPLE_ID and apply fit function to each group
out <- ddply(sampleData, groupBy,
             function(t) tryCatch(fit4pl(t), 
                                  error=function(x) {
                                    rep(NaN,8)
                                  })
)


names(out) <- c("SAMPLE_ID", "EXPERIMENT_ID", "PROTOCOL_ID", "Slope", "Lower Limit", "Upper Limit", "Inflection Point", "ED Estimate", "ED Std Error", "ED Lower CI", "ED Upper CI")
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
out[is.na(out)] <- NA            

#Function to round all numeric columns in a data frame
round_df <- function(df, digits) {
  nums <- vapply(df, is.numeric, FUN.VALUE = logical(1))
  
  df[,nums] <- round(df[,nums], digits = digits)
  
  (df)
}

out <- round_df(out,30)


#tf <- daply(out, ~SAMPLE_ID, function(x) nchar(x$SAMPLE_ID)) > 10
#out[c(which(tf)),]

help(ddply)
fit4pl(tmp)

### test with generated data
# variables
n <- 22 # number of data points
t <- 10^(seq(log10(0.1),log10(100),,22))
min <- 0
max <- 100
hill <- -.7
X50 <- 20
c.unif <- runif(n)
c.norm <- rnorm(n) *2


# generate data and calculate "y"
set.seed(1)
y1 <- min + ((max - min) / (1 + 10^(hill*(log10(X50)- log10(t)))))  + c.norm # uniform error

fit <- drm(y1~t, 
           fct = LL.4(fixed=c(NA, NA, NA, NA), 
          names = c("Slope", "Lower Limit", "Upper Limit", "ed")), 
           control = drmc(errorm = FALSE, warnVal = -1, noMessage = TRUE))
ED(fit,50)
plot(fit)

tmp <- subset.data.frame(Input, Input$SAMPLE_ID == "WEHI-1213557-001")
fit <- drm(tmp$RESPONSE~tmp$CONC, 
           fct = LL.5(fixed=c(NA, 1, NA, NA, 1), 
                      names = c("Slope", "Lower Limit", "Upper Limit", "ed", "skew")), 
           control = drmc(errorm = TRUE, warnVal = 1, noMessage = FALSE))

lowerlimit <- c(-Inf,  -Inf, 0)
upperlimit <- c(Inf,  maxUpper, maxInflection)

fit <- drm(RESPONSE~CONC, 
           data = tmp, 
           fct = LL.5(fixed=c(NA, 0, NA, NA, 1), 
                      names = c("Slope", "Lower Limit", "Upper Limit", "ed", "Skew")), 
           control = drmc(errorm = FALSE, warnVal = -1, noMessage = TRUE),
           lowerl = lowerlimit,
           upperl = upperlimit
)
fit

out <- ddply(sampleData, "SAMPLE_ID", fit4pl)
out <- dlply(sampleData,"SAMPLE_ID",
               function(t) tryCatch(fit4pl(t), 
                                    error=function(x) {
                                      cat("error occurred for:\n", 
                                          "\n...skipping this ticker\n")
                                    })
)
