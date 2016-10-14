#####################################################################
# OIB footprint testing sandbox		 							    #
# JKing 13/10/16													#
# Input: 															#
# Output: 															#
#####################################################################

#Required libraries
require(gstat) 
require(rgdal)
require(sp)
require(rgeos)
require(ggplot2)
require(raster)

#####################################################################
# FUNCTIONS															#
#####################################################################

#Gives barring between from two points
bearingDeg <- function(yy,xx) {
  return(90 - (180/pi)*atan2(yy, xx))
}

#Degrees to radians
deg2rad <- function(deg) {(deg * pi) / (180)}

#Creates a maximum extent polygon of a set of spatial points
extentPoly <- function(bbox, offset, GCS="+proj=longlat +ellps=WGS84", PCS=NULL){
  rownames(bbox) <- c("lon","lat")
  colnames(bbox) <- c('min','max')
  bbox[1,1] = bbox[1,1] - offset
  bbox[2,1] = bbox[2,1] - offset
  bbox[1,2] = bbox[1,2] + offset
  bbox[2,2] = bbox[2,2] + offset
  # clockwise, 5 points to close it
  bboxMat <- rbind( c(bbox['lon','min'],bbox['lat','min']), c(bbox['lon','min'],
    bbox['lat','max']), c(bbox['lon','max'],bbox['lat','max']), c(bbox['lon','max'],
	bbox['lat','min']), c(bbox['lon','min'],bbox['lat','min']) ) 
  bboxSP <- SpatialPolygons( list(Polygons(list(Polygon(bboxMat)),"bbox")))

  if(!is.null(PCS)){
    proj4string(bboxSP) <-CRS(PCS)
  } else {
    proj4string(bboxSP) <- CRS(GCS)
  }
  
  return(bboxSP)
}

#########################################
# CONSTANTS								#
#########################################
setwd("C:/Users/KingJ/Documents/R") #Set this directory to where the script is located
snowfile <- "April13LiDARSite_Magna_ArcReady.csv"
fpAlong <- 5 # in meters alongtrack radar footprint
fpAcross <- 10 # in meters acrosstrack radar footprint; 
numLines <- 100  #must be even; number of flight lines; 
lengthLine <- 200 #in meters; length of each flight line; 
GCS <- "+proj=longlat +ellps=WGS84" #CRS geographic coordniate system WGS84; 
PCS <- "+proj=utm +zone=16 +north +units=m +ellps=WGS84" #CRS projected coordniate system UTM16N/WGS84

#########################################
# PROCESS DATA  						#
#########################################
#Load raw magnaprobe data
snow <-read.csv(file=snowfile,head=TRUE,sep=",")
coordinates(snow) <- ~ Longitude+Latitude

#Identify and remove GPS duplicates 
zd <- zerodist(snow, zero = 0.0)
snow <- snow[-zd[,2], ]

#Identify and remove maxed out measurements
snow <- snow[snow$DepthCm!=120,]

#Set coordinate system and project to UTM 16N
proj4string(snow) <- CRS(GCS)
snow <- spTransform(snow, CRS(PCS))

hist(snow$DepthCm, breaks=20,freq = FALSE, main="Magnaprobe distrobution", xlab="Snow depth (cm)", ylab="PDF")
mean(snow$DepthCm)
min(snow$DepthCm)
max(snow$DepthCm)

#Create random ponits within the observation area, match them as pairs, and build pseudo flight lines
mpExtent <- extentPoly(bbox(snow),-10, GCS,PCS)
randomP <- spsample(mpExtent, numLines*2,"random")
randomDF <- data.frame(ID=1:length(randomP))
randomDF$lineId[sample(1:nrow(randomDF), nrow(randomDF), FALSE)] <- rep(1:(numLines), rep(2, numLines))
randomPDF <- SpatialPointsDataFrame(randomP,randomDF)
fpList <- lapply(split(randomPDF, randomPDF$lineId), function(x) Lines(list(Line(coordinates(x))), x$lineId[1L])) #fp is short for flight path
fpLength <- unlist(lapply(1:length(fpList), function(x) {sqrt((coordinates(fpList[[x]])[[1]][1] - coordinates(fpList[[x]])[[1]][2])^2 + (coordinates(fpList[[x]])[[1]][3] - coordinates(fpList[[x]])[[1]][4])^2)}))
fpDirection <- unlist(lapply(1:length(fpList), function(x) {bearingDeg(coordinates(fpList[[x]])[[1]][1] - coordinates(fpList[[x]])[[1]][2],coordinates(fpList[[x]])[[1]][3] - coordinates(fpList[[x]])[[1]][4])}))

#Create spatial lines
fpLines <- SpatialLines(fpList)
fpData <- data.frame(line = unique(randomPDF$lineId), length=fpLength)
rownames(fpData)<-fpData$line
fpLines <- SpatialLinesDataFrame(fpLines, fpData)
proj4string(fpLines) <- CRS(PCS)

#Plot data
mpPlot.df = snow@data
mpPlot.df$x = coordinates(snow)[,1]
mpPlot.df$y = coordinates(snow)[,2]
dev.new()
ggplot(mpPlot.df, aes(x,y))+
  geom_point(aes(colour = mpPlot.df$DepthCm), size = 1)+
  scale_colour_gradient(low = "green", high = "red", limit=c(min( mpPlot.df$DepthCm),max( mpPlot.df$DepthCm)), name="Snow depth (cm)")+
  labs(title = "Magnaprobe FYI Snow Grid", x="Northing (m)", y="Easting(m)")

#This extends the randomly generated lines to a min length specified in the constants
npDf <- data.frame(lineId = rep(0, length(fpList)*2), x= rep(0, length(fpList)), y = rep(0, length(fpList)))
for (a in 1:length(fpList)){
  nx1 <- coordinates(fpList[[a]])[[1]][1]+((lengthLine - fpLength[a])/2)*cos(deg2rad(fpDirection[a]))
  ny1 <- coordinates(fpList[[a]])[[1]][3]+((lengthLine - fpLength[a])/2)*sin(deg2rad(fpDirection[a]))

  nx2 <- coordinates(fpList[[a]])[[1]][2]-((lengthLine - fpLength[a])/2)*cos(deg2rad(fpDirection[a]))
  ny2 <- coordinates(fpList[[a]])[[1]][4]-((lengthLine - fpLength[a])/2)*sin(deg2rad(fpDirection[a]))
	
  npDf[(a-1)*2+1,] = as.numeric(c(fpList[[a]]@ID,nx1,ny1))
  npDf[(a-1)*2+2,] = as.numeric(c(fpList[[a]]@ID,nx2,ny2))
}

coordinates(npDf) <- ~x+y
nfpList <- lapply(split(npDf, npDf$lineId), function(x) Lines(list(Line(coordinates(x))), x$lineId[1L]))
nfpLines <-SpatialLines(nfpList)
nfpData <- data.frame(line = unique(npDf$lineId), length=lengthLine)
rownames(fpData)<-nfpData$line
nfpLines <- SpatialLinesDataFrame(nfpLines, nfpData)
proj4string(nfpLines) <- CRS(PCS)

#Plot the lines over the extent box.
dev.new()
plot(mpExtent)
lines(nfpLines)

#Generate the snow radar footprints
#TODO: Write a function that can mimic variations in filght atitude.
numFp <- lengthLine/fpAlong
radarPoint <- lapply(1:length(nfpLines), function(x) {spsample(nfpLines[x,], numFp,type="regular")})
radarBuffer <- lapply(1:length(radarPoint), function(x) gBuffer(radarPoint[[x]], width = fpAlong/2, byid=TRUE,capStyle="ROUND", quadsegs=10))
radarSeg <- lapply(1:length(radarBuffer), function(x) gIntersection(nfpLines[x,],radarBuffer[[x]], byid=TRUE))
radarFootprint <- lapply(1:length(radarSeg), function(x) gBuffer(radarSeg[[x]], width=fpAcross/2, byid=TRUE, capStyle="FLAT")) #TODO, adjust width dynamicaly based on topography
radarFootprint <- lapply(1:length(radarSeg), function(x) radarFootprint[[x]][-c(1,length(radarFootprint[[x]]))]) #last one might be the wrong size, remove incase
radarDf <- lapply(1:length(radarFootprint), function(x) data.frame(Line = rep(x,length(radarFootprint[[x]])))) #add barring of the flight line here
radarFootprint <- do.call(bind, radarFootprint)
radarDf <- do.call(bind, radarDf)
radarFootprint <- SpatialPolygonsDataFrame(radarFootprint,radarDf)
radarPointMerge <- lapply(1:length(radarPoint), function(x) radarPoint[[x]][-c(1,length(radarPoint[[x]]))]) #last one might be the wrong size, remove incase
radarPointMerge <- do.call(bind, radarPointMerge)

#Display the generated footprints along with the in situ observations
ggplot(mpPlot.df, aes(x,y))+
  scale_colour_gradient(low = "green", high = "red", limit=c(min( mpPlot.df$DepthCm),max( mpPlot.df$DepthCm)), name="Snow depth (cm)")+
  geom_polygon(data=radarFootprint, aes(x=long, y=lat,group=id), alpha=1,colour="darkgrey", fill=NA, size=0.1)+	
  geom_point(aes(colour = mpPlot.df$DepthCm), size = 1)+
  labs(title = "Synthetic snowradar footprints", x="Northing (m)", y="Easting(m)")

#Extract point observations within each footprint and derive statistics		 
mpPolylist <- over(radarFootprint, snow, returnList = TRUE, fn=NULL) #get point list to do precentile
resultMp <- data.frame(count = rep(NA,length(mpPolylist)),mean = rep(NA,length(mpPolylist)), sd = rep(NA,length(mpPolylist)))
for ( i in 1:nrow(resultMp)){
	 mpPoints.df = as.data.frame(mpPolylist[[i]])
		if(nrow( mpPoints.df)>0){
			resultMp$mean[i] = mean(mpPoints.df$DepthCm)
			resultMp$count[i] = nrow(mpPoints.df)
			resultMp$sd[i] = sd(mpPoints.df$DepthCm)
		}
	}

#Extract NN in situ measurement for each footprint centroid, this emulates the CRESIS lat long phase center
distNN <- gDistance(radarPointMerge, snow,byid=TRUE)
NN <- apply(distNN , 2, function(X) rownames(distNN )[order(X)][1])
resultMp$nnDist <- round(apply(distNN , 2, function(X) min(X)),2)
resultMp$nnVal <- snow$DepthCm[as.numeric(NN)]

#Clean up analysis; remove footprints with less than 10 in situ meausrements
resultMp <- resultMp[complete.cases(resultMp),]
resultMp <- resultMp[resultMp$count>10,]

#Plot relationship between nn and integrated sampling
plot(resultMp$nnVal, resultMp$mean,
  xlab="NN depth (cm)", ylab="Footprint depth (cm)",
  xlim=c(0, 120), ylim=c(0, 120))
  abline(0,1)
rmseMp = sqrt(sum((resultMp$nnVal-resultMp$mean)^2)/nrow(resultMp))
biasMp = mean(resultMp$nnVal-resultMp$mean)
corMp = cor(resultMp$nnVal,resultMp$mean)
text(10, 105, paste0("R = ",round(corMp,2)), cex = .8) 
text(10, 110, paste0("Bias = ",round(biasMp,2)), cex = .8) 
text(10, 115, paste0("RMSE = ",round(rmseMp,2)), cex = .8)  
  
