# Based on script written by Rick Reneau.

#setwd("your directory here")

decimal <- 4
file_list <- list.files(getwd(),"^1")
date_list <- substr(file_list,10,19)
subj_list <- substr(file_list,1,8)
subj_list <- subj_list[!duplicated(subj_list)]
start_date <- date_list[1]
end_date <- date_list[length(date_list)]
DaysVec <- seq(as.Date(start_date), as.Date(end_date), by="days")
final_df <- as.data.frame(matrix(nrow = length(subj_list)+1,ncol = length(DaysVec)+1))
final_df[1,1] <- "subject_ID"
final_df[1,2:ncol(final_df)] <- as.character(DaysVec)
final_df[2:nrow(final_df),1] <- subj_list
final_df[1,1] <- "subject_ID"

for(a in 1:length(file_list)){
  date1 <- date_list[a]
  subj <- substr(file_list[a],1,8)
  for(g in 1:ncol(final_df)){
    if(grepl(date1,final_df[1,g])){
      colindex <- g
      break
    }
  }
  for(h in 1:nrow(final_df)){
    if(grepl(subj,final_df[h,1])){
      rowindex <- h
      break
    }
  }
  ftable <- read.csv(file_list[a]) # Read in file for participant (single day's file).
  if(nrow(ftable)==1){
    final_df[rowindex,colindex] <- NA
    next
  }
  pointtable <- ftable[,2:3]
  library(fields)
  
  # Filter.
  minimum_distance <- .05
  largest_angle <- .4
  b <- 2
  while(b<(nrow(ftable))){
    p12 <- rdist.earth(pointtable[b-1,],pointtable[b,])
    p23 <- rdist.earth(pointtable[b,],pointtable[b+1,])
    p13 <- rdist.earth(pointtable[b-1,],pointtable[b+1,])
    angle <- (p12^2+p23^2-p13^2)/(2*p12*p23)
    if(is.na(angle)){
      b <- b+1
      next()
    }
    if(angle > 1)
      angle <- 1
    if(angle < -1)
      angle <- -1
    if(acos(angle)<largest_angle){
      ftable <- ftable[-(b),]
      pointtable <- pointtable[-(b),]
    }
    else{
      b <- b+1
    }
  }
  c <- 1
  while(c<(nrow(ftable)-1)){
    if(rdist.earth(pointtable[c,],pointtable[c+1,])<minimum_distance){
      ftable <- ftable[-(c+1),]
      pointtable <- pointtable[-(c+1),]
    }
    else{
      c <- c+1
    }
  }
  LonVec <- round(as.numeric(pointtable[,2]), digits = decimal)
  LatVec <- round(as.numeric(pointtable[,1]), digits = decimal)
  DateVec <- substr(ftable[,1],1,10)
  TimeVec <- substr(ftable[,1],12,19)
  DateTimeVec <- substr(ftable[,1],1,19)
  positionMatrix = matrix(
    c(LonVec,LatVec),
    nrow=length(LonVec),
    ncol=2)

    
    # Initialize day matrix.
    dmsize = 1440
    dayMatrix = matrix(
      nrow=dmsize,
      ncol=2)
    
    # Find data points with given date.
    dataLocation <- character(2000) # Preallocating for max possible number of trackpoints in a day.
    f <- 1
    for (t in 1:length(DateVec)){
      if(DateVec[t]==date1){
        dataLocation[f] <- t
        f <- f +1
      }
      if(DateVec[t]>date1)
        break
    }
    dataLocation <- dataLocation[!dataLocation %in% ""]
    
    beginDate <- as.POSIXct(as.character(date1))
    
    # Place coordinates in proper timeline (each row of a matrix represents a minute of the timeframe of collected data).
    if(length(dataLocation) >= 1)
    {
      
      # Check for daylight savings.
       if(difftime(as.POSIXct(DateTimeVec[as.numeric(dataLocation[length(dataLocation)])],format = "%Y-%m-%dT%H:%M:%S"),
                  as.POSIXct(as.character(date1)),units = "mins")>1440){
        dmsize <- 1500
        dayMatrix = matrix(
          nrow=dmsize,
          ncol=2
        )
      }
      
      startPoint <- 1
      for(u in 1:length(dataLocation)){
        tTest <- as.POSIXct(DateTimeVec[as.numeric(dataLocation[u])],format = "%Y-%m-%dT%H:%M:%S")
        timeInMinutes <- floor(as.numeric(difftime(tTest,beginDate,units = "mins")))+1
        
        if(startPoint==1){
          for(v in 1:(timeInMinutes-1)){
            dayMatrix[v,] <- positionMatrix[as.numeric(dataLocation[u]),]
          }
        }
        
        if(u==length(dataLocation)&&startPoint<=dmsize){
          for(w in startPoint:dmsize){
            dayMatrix[w,] <- positionMatrix[as.numeric(dataLocation[u]),]
          }
        }
        
        if(timeInMinutes-startPoint>0&&u!=1&&u!=length(dataLocation)){
          
          previous_lat <- as.numeric(positionMatrix[as.numeric(dataLocation[u])-1,1])
          current_lat <- as.numeric(positionMatrix[as.numeric(dataLocation[u]),1])
          previous_lon <- as.numeric(positionMatrix[as.numeric(dataLocation[u])-1,2])
          current_lon <- as.numeric(positionMatrix[as.numeric(dataLocation[u]),2])
          previous_time <- startPoint-1
          current_time <- timeInMinutes
          avg_lat_velocity <- (current_lat-previous_lat)/(current_time-previous_time)
          avg_lon_velocity <- (current_lon-previous_lon)/(current_time-previous_time)
          
          for(z in startPoint:(timeInMinutes-1)){
            dayMatrix[z,1] <- round((avg_lat_velocity*(z-previous_time)+previous_lat), digits = decimal)
            dayMatrix[z,2] <- round((avg_lon_velocity*(z-previous_time)+previous_lon), digits = decimal)
          }
        }
        
        dayMatrix[timeInMinutes,] <- positionMatrix[as.numeric(dataLocation[u]),]
        startPoint <- timeInMinutes+1  
      } 
    }
    
    # Convert matrix to dataframe.
    df <- as.data.frame(dayMatrix,row.names=NULL)
    colnames(df) <- c("lon","lat")
    
    # Calculate roaming entropy for that date and place it in RE vector.
    library(plyr)
    uniquePositions <- ddply(df,.(lon,lat),nrow)
    uniquePositions$p <- uniquePositions$V1/dmsize
    # Take natural log and divide by log(648000000).
    uniquePositions$plogp <- uniquePositions$p*log(uniquePositions$p)
    sampleRE = (-sum(uniquePositions$plogp))/log(648000000)
    final_df[rowindex,colindex] <- sampleRE
}
setwd("/Users/nataliesaragosaharris/Desktop")
write.csv(final_df,file = "MSD_Entropy_Filtered.csv", row.names = FALSE)