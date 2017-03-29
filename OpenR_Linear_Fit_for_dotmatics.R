library(plyr)

#Inputs
Input
fitForEach

linfit <- function (df) {
  groups <- c(UNIQUE_PROP_ID = paste(unique(df$UNIQUE_PROP_ID), collapse = ", "))
  groups <- c(groups, SAMPLE_PLATE_ID = paste(unique(df$SAMPLE_PLATE_ID), collapse = ", "))
  
  df$Log10CONC <- log10(df$CONC)
  fit<-nls(RESPONSE ~ yIntercept + slope * Log10CONC,
      data = df,
      start = list(slope= 0, yIntercept = 0))
  params <- summary(fit)$parameters
  flat <- c(params)
  
  names(flat) <- c(outer(paste("LinFit", rownames(params)), colnames(params), paste, sep = ":"))
  return(c(groups,flat))
  }

sampleData <- subset.data.frame(Input, Input$SAMPTYPE == "S")
sampleData <- subset.data.frame(sampleData, (sampleData$KNOCKOUT != "Y" | is.na(sampleData$KNOCKOUT)))

#Return fits for each sample, or each sample on each plate
groupBy <- if (fitForEach == 'Sample') ~ SAMPLE_ID + EXPERIMENT_ID + PROTOCOL_ID else ~ SAMPLE_ID + EXPERIMENT_ID + SAMPLE_PLATE_ID + PROTOCOL_ID

out <- ddply(sampleData, groupBy, linfit)

#Convert specified columns to numeric
cols = c(2, 6, 7,8,9,10,11,12,13);    
out[,cols] = apply(out[,cols], 2, function(x) as.numeric(as.character(x)))

#Function to round all numeric columns in a data frame
round_df <- function(df, digits) {
  nums <- vapply(df, is.numeric, FUN.VALUE = logical(1))
  
  df[,nums] <- round(df[,nums], digits = digits)
  
  (df)
}

out <- round_df(out,30)
