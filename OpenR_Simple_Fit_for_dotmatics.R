#Remove excluded points from calculation
sampleData <- subset.data.frame(Input, (Input$KNOCKOUT != "Y" | is.na(Input$KNOCKOUT)))
sampleData <- subset.data.frame(sampleData, sampleData$SAMPTYPE == "S")

simpleX50 <- function(df) {
  max <- max(df$RESPONSE)
  min <- min(df$RESPONSE)
  mid <- mean(c(max,min))
  
  topHalf <- subset.data.frame(df, df$RESPONSE >= mid)
  bottomHalf <- subset.data.frame(df, df$RESPONSE <= mid)
  
  midUpper <- subset.data.frame(df, df$RESPONSE == min(topHalf$RESPONSE))
  midLower <- subset.data.frame(df, df$RESPONSE == max(bottomHalf$RESPONSE))
  
  midUpperConc <- midUpper$CONC
  midLowerConc <- midLower$CONC
  midUpperResponse <- midUpper$RESPONSE
  midLowerResponse <- midLower$RESPONSE
 
  xdif <- midUpperConc - midLowerConc
  ydif <- midUpperResponse - midLowerResponse
  
  m <- ydif / xdif
  midConc <- midLowerConc + ((mid - midLowerResponse) / m)
  
  
  final <- c("midUpperConc" = midUpperConc,
             "midLowerConc" = midLowerConc,
             "midUpperResponse" = midUpperResponse,
             "midLowerResponse" = midLowerResponse,
             "midResponse" = mid,
             "midConc" =midConc)
}

sampleOnly <- ~ SAMPLE_ID + EXPERIMENT_ID + PROTOCOL_ID
sampleAndPlate <- ~ SAMPLE_ID + EXPERIMENT_ID + SAMPLE_PLATE_ID + PROTOCOL_ID

outSimple <- ddply(sampleData,sampleOnly,simpleX50)
