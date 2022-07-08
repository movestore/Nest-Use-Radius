library('move')
library('foreach')
library('geosphere')
library('ggmap')
library('ggplot2')
library('sp')
library('rgeos')
library('fields')

rFunction <- function(data,radii=500,selName=NULL,trackVar=NULL)
{
  Sys.setenv(tz="UTC")
  
  # create circles data frame from the centers data frame
  make_circles <- function(centers, radius, nPoints = 100){
    # centers: the data frame of centers with ID
    # radius: radius measured in kilometer
    #
    meanLat <- mean(centers$location.lat)
    # length per longitude changes with lattitude, so need correction
    radiusLon <- radius/1000 /111 / cos(meanLat/57.3) 
    radiusLat <- radius/1000 / 111
    circleDF <- data.frame(ID = rep(centers$trackId, each = nPoints))
    angle <- seq(0,2*pi,length.out = nPoints)
    
    circleDF$lon <- unlist(lapply(centers$location.long, function(x) x + radiusLon * cos(angle)))
    circleDF$lat <- unlist(lapply(centers$location.lat, function(x) x + radiusLat * sin(angle)))
    return(circleDF)
  }
  
  
  if (is.null(selName) | (selName %in% namesIndiv(data))==FALSE) 
  {
    logger.info("Your property name (e.g. 'nest') is not set or does not exist in the data. This will lead to an error.")
    result <- NULL
  } else
  {
    ix <- which(namesIndiv(data)==selName)
    selT <- data[[ix]]
    data <- data[[-ix]]
    
    selT@data[,c("location.long","location.lat")] <- coordinates(selT)
    
    data.split <- move::split(data)
  
    radiuss <- sort(as.numeric(trimws(strsplit(as.character(radii),",")[[1]])),decreasing=FALSE)
    n.rad <- length(radiuss)
    ids <- namesIndiv(data)[(which(namesIndiv(data) %in% selT@data[,trackVar]))] # names of tracks in data for which nests are available
    n.ids <- length(ids)

    rad_table <- data.frame("track"=rep(c(ids,"mean","sd"),each=n.rad),"radius"=rep(radiuss,times=n.ids+2),"n.loc"=numeric((n.ids+2)*n.rad),"prop.locs"=numeric((n.ids+2)*n.rad),"prop.dur"=numeric((n.ids+2)*n.rad))
    
    for (i in seq(along=data.split))
    {
      datai <- data.split[[i]]
      nest.long <- coordinates(selT)[selT@data[,trackVar]==namesIndiv(datai),1] #nest locations are automatically the locations in the element called selVar
      nest.lat <- coordinates(selT)[selT@data[,trackVar]==namesIndiv(datai),2]
      
      #dur <- c(timeLag(datai,units="hours"),NA)
      dur <- datai@data$timelag # from TimeLag App
      dist.nest <- distVincentyEllipsoid(p1=c(nest.long,nest.lat),p2=coordinates(datai)) #metres
      n.loc <- prop.loc <- prop.dur <- numeric(n.rad)
      
      for (j in seq(along=radiuss))
      {
        ix <- which(dist.nest<radiuss[j])
        n.loc[j] <- length(ix)
        prop.loc[j] <- n.loc[j]/length(datai)
        prop.dur[j] <- sum(dur[ix],na.rm=TRUE)/sum(dur,na.rm=TRUE)
      }
      
      rad_table[which(rad_table$track==namesIndiv(datai)),3:5] <- data.frame(n.loc,prop.loc,prop.dur)
    }
    
    for (k in seq(along=radiuss))
    {
      ixk <- which(rad_table$radius==radiuss[k] & rad_table$track %in% ids)
      rad_table[which(rad_table$track=="mean" & rad_table$radius==radiuss[k]),3] <- length(rad_table$n.loc[ixk])
      rad_table[which(rad_table$track=="mean" & rad_table$radius==radiuss[k]),4] <- mean(rad_table$prop.loc[ixk])
      rad_table[which(rad_table$track=="mean" & rad_table$radius==radiuss[k]),5] <- mean(rad_table$prop.dur[ixk])
      
      rad_table[which(rad_table$track=="sd" & rad_table$radius==radiuss[k]),3] <- length(rad_table$n.loc[ixk])
      rad_table[which(rad_table$track=="sd" & rad_table$radius==radiuss[k]),4] <- sd(rad_table$prop.loc[ixk])
      rad_table[which(rad_table$track=="sd" & rad_table$radius==radiuss[k]),5] <- sd(rad_table$prop.dur[ixk])
    }
    
    write.csv(rad_table,paste0(Sys.getenv(x = "APP_ARTIFACTS_DIR", "/tmp/"), "Radius_NestUse.csv"),row.names=FALSE)
    
    selT.df <- as.data.frame(selT)
    
    pdf(paste0(Sys.getenv(x = "APP_ARTIFACTS_DIR", "/tmp/"), "Tracks_withRadii_onMap.pdf"),width=12,height=8)
    lapply(data.split, function(z){
      
      map <- get_map(bbox(extent(z)),source="osm",force=TRUE,zoom=10)
      if (any(names(z)=="location.long"))
      {
        datai.df <- as.data.frame(moveStack(z))
      } else 
      {
          datai.df <- data.frame(coordinates(z),as.data.frame(moveStack(z)))
          names(datai.df)[1:2] <- c("location.long","location.lat")
      }
      
      selTi.df <- selT.df[selT.df[,trackVar]==namesIndiv(z),]
      circs.i <- make_circles(selTi.df,radiuss[1])
      for (i in seq(along=radiuss)[-1]) circs.i <- rbind(circs.i,make_circles(selTi.df,radiuss[i]))
      
      mymap <- ggmap(map) +
          geom_path(data=datai.df, 
                    aes(x=location.long, y=location.lat, col=trackId),show.legend=FALSE) +
          geom_point(data=selTi.df,
                     aes(x=location.long, y=location.lat),
                     size=3,color="red",alpha=0.8,shape=17) +
        geom_polygon(data = circs.i, aes(lon, lat, group = ID), color = "blue", alpha = 0) +
          theme(legend.justification = "top") +
          labs(x="Longitude", y="Latitude") +
        scale_fill_manual(name="Track", values=tim.colors(length(namesIndiv(data))),aesthetics=c("colour","fill"))
        mymap
      })
      dev.off()
  }
  
  result <- data #return full data set
  return(result)
}

  
  
  
  
  
  
  
  
  
  
