library(plyr)
linfit <- function (df) {
  df$Log10CONC <- log10(df$CONC)
  fit<-nls(RESPONSE ~ yIntercept + slope * Log10CONC,
      data = df,
      start = list(slope= 0, yIntercept = 0))
  params <- summary(fit)$parameters
  flat <- c(params)
  names(flat) <- c(outer(paste("LinFit", rownames(params)), colnames(params), paste, sep = ":"))
  flat
  }

sampleData <- subset.data.frame(Input, Input$SAMPTYPE == "S")
sampleData <- subset.data.frame(sampleData, (sampleData$KNOCKOUT != "Y" | is.na(sampleData$KNOCKOUT)))

out <- ddply(sampleData, "SAMPLE_ID", linfit)

#Function to round all numeric columns in a data frame
round_df <- function(df, digits) {
  nums <- vapply(df, is.numeric, FUN.VALUE = logical(1))
  
  df[,nums] <- round(df[,nums], digits = digits)
  
  (df)
}

out <- round_df(out,30)
