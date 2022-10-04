library('move')
library('foreach')
library('geosphere')
library('ggmap')
library('ggplot2')
library('sp')
library('rgeos')
library('fields')
library('grid')
library('gridExtra')

rFunction <- function(data,radii=500,selName=NULL,trackVar=NULL,gap_adapt=FALSE)
{
  Sys.setenv(tz="UTC")
  
  if (gap_adapt==TRUE) TL <- "timelag2" else TL <- "timelag"
  
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
    
    data.split.sel <- data.split[which(names(data.split) %in% ids)]
    notin <- names(data.split[which((names(data.split) %in% ids)==FALSE)])
    if (length(notin)>0) logger.info (paste("The tracks",paste(notin,collapse=", "),"do not relate to any detected nesting attempt. They will not be analysed further, but are included in the output data."))

    rad_table <- data.frame("track"=rep(c(ids,"mean","sd"),each=n.rad),"radius"=rep(radiuss,times=n.ids+2),"n.loc"=numeric((n.ids+2)*n.rad),"prop.locs"=numeric((n.ids+2)*n.rad),"prop.dur"=numeric((n.ids+2)*n.rad))
    avg.table <- data.frame("trackId"=c(ids,"all"),"n.pts"=numeric(n.ids+1),"mean.pts.dist"=numeric(n.ids+1),"sd.pts.dist"=numeric(n.ids+1),"mean.dur.dist"=numeric(n.ids+1),"sd.dur.dist"=numeric(n.ids+1))
    
    pdf(paste0(Sys.getenv(x = "APP_ARTIFACTS_DIR", "/tmp/"), "Histograms_Dist2Nest.pdf"),width=12,height=8)
    out_sel <- dists_all <- dur_all <- numeric()
    for (i in seq(along=data.split.sel))
    {
      datai <- data.split.sel[[i]]
      seli <- which(selT@data[,trackVar]==namesIndiv(datai))
      if (length(seli)>1) 
        {
        logger.info(paste("The track",namesIndiv(datai),"is related to more than one detected nesting site. Here the first nest candidate is used. Go back and adapt your nesting file (delete row of incorrect duplicate nesting sites) if this does not work for you."))
        out_sel <- c(out_sel,seli[-1])
        }
      

      nest.long <- coordinates(selT)[selT@data[,trackVar]==namesIndiv(datai),1][1] #nest locations are automatically the locations in the element called selVar
      nest.lat <- coordinates(selT)[selT@data[,trackVar]==namesIndiv(datai),2][1]
        
      dur <- datai@data[,TL] # from TimeLag App
      dist.nest <- distVincentyEllipsoid(p1=c(nest.long,nest.lat),p2=coordinates(datai)) #metres
        
      if (any(!is.na(dist.nest)))
      {
        min_dist <- min(dist.nest,na.rm=TRUE)
        if (min_dist<0) min_dist <- c(min_dist,0)
        max_dist <- max(dist.nest,na.rm=TRUE)
        hist(dist.nest,xlim=c(quantile(dist.nest,probs=0.01,na.rm=TRUE),quantile(dist.nest,probs=0.99,na.rm=TRUE)),breaks=c(min_dist,radiuss,max_dist),main=paste("Histogramme of", namesIndiv(datai)),xlab="distance to nest",freq=FALSE,col="blue",ylab="Probability density")
      }
      
      dists_all <- c(dists_all,dist.nest)
      dur_all <- c(dur_all,dur)
      
      avg.table$n.pts[i] <- length(dist.nest)
      avg.table$mean.pts.dist[i] <- mean(dist.nest,na.rm=TRUE)
      avg.table$sd.pts.dist[i] <- sd(dist.nest,na.rm=TRUE)
      mu_i <- sum(dist.nest*dur,na.rm=TRUE)/sum(dur,na.rm=TRUE)
      avg.table$mean.dur.dist[i] <- mu_i
      avg.table$sd.dur.dist[i] <- sqrt(sum((dist.nest-mu_i)*(dist.nest-mu_i)*dur,na.rm=TRUE)/sum(dur,na.rm=TRUE)) #sqrt(weighted variance)
        
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
    
    if (any(!is.na(dists_all)))
    {
      min_distA <- min(dists_all,na.rm=TRUE)
      if(min_distA<0) min_distA <- c(min_distA,0)
      max_distA <- max(dists_all,na.rm=TRUE)
      hist(dist.nest,xlim=c(quantile(dists_all,probs=0.01,na.rm=TRUE),quantile(dists_all,probs=0.99,na.rm=TRUE)),breaks=c(min_distA,radiuss,max_distA),main=paste("Histogramme of all tracks' locations"),xlab="distance to nest of respective individual",freq=FALSE,col="red",ylab="Probability density")
    }
    dev.off()
    
    avg.table$n.pts[n.ids+1] <- length(dists_all)
    avg.table$mean.pts.dist[n.ids+1] <- mean(dists_all,na.rm=TRUE)
    avg.table$sd.pts.dist[n.ids+1] <- sd(dists_all,na.rm=TRUE)
    mu_n <- sum(dists_all*dur_all,na.rm=TRUE)/sum(dur_all,na.rm=TRUE)
    avg.table$mean.dur.dist[n.ids+1] <- mu_n
    avg.table$sd.dur.dist[n.ids+1] <- sqrt(sum((dists_all-mu_n)*(dists_all-mu_n)*dur_all,na.rm=TRUE)/sum(dur_all,na.rm=TRUE)) #sqrt(weighted variance)
    
    for (k in seq(along=radiuss))
    {
      ixk <- which(rad_table$radius==radiuss[k] & rad_table$track %in% ids)
      rad_table[which(rad_table$track=="mean" & rad_table$radius==radiuss[k]),3] <- length(rad_table$n.loc[ixk])
      rad_table[which(rad_table$track=="mean" & rad_table$radius==radiuss[k]),4] <- mean(rad_table$prop.loc[ixk],na.rm=TRUE)
      rad_table[which(rad_table$track=="mean" & rad_table$radius==radiuss[k]),5] <- mean(rad_table$prop.dur[ixk],na.rm=TRUE)
      
      rad_table[which(rad_table$track=="sd" & rad_table$radius==radiuss[k]),3] <- length(rad_table$n.loc[ixk])
      rad_table[which(rad_table$track=="sd" & rad_table$radius==radiuss[k]),4] <- sd(rad_table$prop.loc[ixk],na.rm=TRUE)
      rad_table[which(rad_table$track=="sd" & rad_table$radius==radiuss[k]),5] <- sd(rad_table$prop.dur[ixk],na.rm=TRUE)
    }
    
    write.csv(avg.table,paste0(Sys.getenv(x = "APP_ARTIFACTS_DIR", "/tmp/"), "Avg_NestDists.csv"),row.names=FALSE)
    write.csv(rad_table,paste0(Sys.getenv(x = "APP_ARTIFACTS_DIR", "/tmp/"), "Radius_NestUseProps.csv"),row.names=FALSE)
    
    if (length(out_sel)>0) #only plot tracks with first detected nest per track
    {
      selT.df <- as.data.frame(selT[-out_sel])
      data.split.nn <- data.split.sel[-out_sel]
    } else 
    {
        selT.df <- as.data.frame(selT)
        data.split.nn <- data.split.sel 
    }
    
    mymaps <- lapply(data.split.nn,function(datai) 
    {
       map <- get_map(bbox(extent(datai)),source="osm",force=TRUE,zoom=10)
      if (any(names(datai)=="location.long"))
      {
        datai.df <- as.data.frame(moveStack(datai))
      } else 
      {
          datai.df <- data.frame(coordinates(datai),as.data.frame(moveStack(z)))
          names(datai.df)[1:2] <- c("location.long","location.lat")
      }
      
      selTi.df <- selT.df[selT.df[,trackVar]==namesIndiv(datai),]
      circs.i <- make_circles(selTi.df,radiuss[1])
      logger.info(circs.i[1,])
      for (ii in seq(along=radiuss)[-1]) circs.i <- rbind(circs.i,make_circles(selTi.df,radiuss[ii]))
      
      mymap <- ggmap(map) +
          geom_path(data=datai.df, 
                    aes(x=location.long, y=location.lat, col=trackId),show.legend=FALSE) +
          geom_point(data=datai.df, 
                  aes(x=location.long, y=location.lat, col=trackId),color="orange",size=1) +
          geom_point(data=selTi.df,
                     aes(x=location.long, y=location.lat),
                     size=3,color="red",alpha=0.8,shape=17) +
        geom_polygon(data = circs.i, aes(lon, lat, group = ID), color = "blue", alpha = 0) +
          theme(legend.justification = "top") +
          labs(x="Longitude", y="Latitude",title=namesIndiv(datai)) +
        scale_fill_manual(name="Track", values=tim.colors(length(namesIndiv(data))),aesthetics=c("colour","fill"))
        mymap
    })
    
    mymaps_p  <- marrangeGrob(mymaps, nrow = 1, ncol = 1)
    ggsave(paste0(Sys.getenv(x = "APP_ARTIFACTS_DIR", "/tmp/"),"Tracks_withRadii_onMap.pdf"), plot = mymaps_p, width = 21, height = 29.7, units = "cm")
   
    #pdf(paste0(Sys.getenv(x = "APP_ARTIFACTS_DIR", "/tmp/"), "Tracks_withRadii_onMap.pdf"),width=12,height=8)
    #mymaps
    #dev.off()
  }
  
  result <- data #return full data set
  return(result)
}

  
  
  
  
  
  
  
  
  
  
