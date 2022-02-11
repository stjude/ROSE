#!/usr/bin/env Rscript
#============================================================================
#==============SUPER-ENHANCER CALLING AND PLOTTING FUNCTIONS=================
#============================================================================

#This function calculates the cutoff by sliding a diagonal line and finding where it is tangential (or as close as possible)
calculate_cutoff <- function(inputVector, drawPlot=TRUE,...){
 	inputVector <- sort(inputVector)
	inputVector[inputVector<0]<-0 #set those regions with more control than ranking equal to zero
	slope <- (max(inputVector)-min(inputVector))/length(inputVector) #This is the slope of the line we want to slide. This is the diagonal.
	xPt <- floor(optimize(numPts_below_line,lower=1,upper=length(inputVector),myVector= inputVector,slope=slope)$minimum) #Find the x-axis point where a line passing through that point has the minimum number of points below it. (ie. tangent)
	y_cutoff <- inputVector[xPt] #The y-value at this x point. This is our cutoff.
	
	if(drawPlot){  #if TRUE, draw the plot
		plot(1:length(inputVector), inputVector,type="l",...)
		b <- y_cutoff-(slope* xPt)
		abline(v= xPt,h= y_cutoff,lty=2,col=8)
		points(xPt,y_cutoff,pch=16,cex=0.9,col=2)
		abline(coef=c(b,slope),col=2)
		title(paste("x=",xPt,"\ny=",signif(y_cutoff,3),"\nFold over Median=",signif(y_cutoff/median(inputVector),3),"x\nFold over Mean=",signif(y_cutoff/mean(inputVector),3),"x",sep=""))
		axis(1,sum(inputVector==0),sum(inputVector==0),col.axis="pink",col="pink") #Number of regions with zero signal
	}
	return(list(absolute=y_cutoff,overMedian=y_cutoff/median(inputVector),overMean=y_cutoff/mean(inputVector)))
}

#this is an accessory function, that determines the number of points below a diagnoal passing through [x,yPt]
numPts_below_line <- function(myVector,slope,x){
	yPt <- myVector[x]
	b <- yPt-(slope*x)
	xPts <- 1:length(myVector)
	return(sum(myVector<=(xPts*slope+b)))
}


convert_stitched_to_bed <- function(inputStitched,trackName,trackDescription,outputFile,splitSuper=TRUE,score=c(),superRows=c(),baseColor="0,0,0",superColor="255,0,0"){
	outMatrix <- matrix(data="",ncol=4+ifelse(length(score)==nrow(inputStitched),1,0),nrow=nrow(inputStitched))
	
	outMatrix[,1] <- as.character(inputStitched$CHROM)
	outMatrix[,2] <- as.character(inputStitched$START)
	outMatrix[,3] <- as.character(inputStitched$STOP)
	outMatrix[,4] <- as.character(inputStitched$REGION_ID)
	if(length(score)==nrow(inputStitched)){
		score <- rank(score,ties.method="first")
		score <- length(score)-score+1  #Stupid rank only does smallest to largest. 
		outMatrix[,5] <- as.character(score)
	}
	trackDescription <- paste(trackDescription,"\nCreated on ",format(Sys.time(), "%b %d %Y"),collapse="",sep="")
	trackDescription <- gsub("\n","\t", trackDescription)
	tName <- gsub(" ","_",trackName)
	cat('track name="', tName,'" description="', trackDescription,'" itemRGB=On color=',baseColor,"\n",sep="",file=outputFile)
	write.table(file= outputFile,outMatrix,sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE,append=TRUE)
	if(splitSuper==TRUE){
		cat("\ntrack name=\"Super_", tName,'" description="Super ', trackDescription,'" itemRGB=On color=', superColor,"\n",sep="",file=outputFile,append=TRUE)
		write.table(file= outputFile,outMatrix[superRows,],sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE,append=TRUE)
	}
}

writeSuperEnhancer_table <- function(superEnhancer,description,outputFile,additionalData=NA){
	description <- paste("#",description,"\nCreated on ",format(Sys.time(), "%b %d %Y"),collapse="",sep="")
	description <- gsub("\n","\n#",description)
	cat(description,"\n",file=outputFile)
	if(is.matrix(additionalData)){
		if(nrow(additionalData)!=nrow(superEnhancer)){
			warning("Additional data does not have the same number of rows as the number of super stitched peaks.\n--->>> ADDITIONAL DATA NOT INCLUDED <<<---\n")
		}else{
			superEnhancer <- cbind(superEnhancer,additionalData)
			superEnhancer = superEnhancer[order(superEnhancer$stitchedPeakRank),]
			
		}
	}
	write.table(file=outputFile,superEnhancer,sep="\t",quote=FALSE,row.names=FALSE,append=TRUE)
}



#============================================================================
#===================SUPER-ENHANCER CALLING AND PLOTTING======================
#============================================================================



args <- commandArgs()

print('THESE ARE THE ARGUMENTS')
print(args)

#ARGS
outFolder = args[6] #3
enhancerFile = args[7] #4
enhancerName = args[8] #5
wceName = args[9] #6



#Read enhancer regions with closestGene columns
stitched_regions <- read.delim(file= enhancerFile,sep="\t")


#perform WCE subtraction. Using pipeline table to match samples to proper background. 
rankBy_factor = colnames(stitched_regions)[7]

prefix = unlist(strsplit(rankBy_factor,'_'))[1]

if(wceName == 'NONE'){
	
	
	rankBy_vector = as.numeric(stitched_regions[,7])
	
}else{
	wceName = colnames(stitched_regions)[8]
	print('HERE IS THE WCE NAME')
	print(wceName)
	
	rankBy_vector = as.numeric(stitched_regions[,7])-as.numeric(stitched_regions[,8])
}	
	

#SETTING NEGATIVE VALUES IN THE rankBy_vector to 0

rankBy_vector[rankBy_vector < 0] <- 0


#FIGURING OUT THE CUTOFF

cutoff_options <- calculate_cutoff(rankBy_vector, drawPlot=FALSE,xlab=paste(rankBy_factor,' Stitched peaks'),ylab=paste(rankByFactor,' Signal','- ',wceName),lwd=2,col=4)


#These are the super-enhancers
superEnhancerRows <- which(rankBy_vector> cutoff_options$absolute)
typicalEnhancers = setdiff(1:nrow(stitched_regions),superEnhancerRows)
enhancerDescription <- paste(enhancerName," Stitched Peaks\nCreated from ", enhancerFile,"\nRanked by ",rankBy_factor,"\nUsing cutoff of ",cutoff_options$absolute," for Super-Stitched Peaks",sep="",collapse="")


#MAKING HOCKEY STICK PLOT
plotFileName = paste(outFolder,enhancerName,'_Plot_points.png',sep='')
png(filename=plotFileName,height=600,width=600)
signalOrder = order(rankBy_vector,decreasing=TRUE)
if(wceName == 'NONE'){
	plot(length(rankBy_vector):1,rankBy_vector[signalOrder], col='red',xlab=paste(rankBy_factor,' Stitched peaks'),ylab=paste(rankBy_factor,' Signal'),pch=19,cex=2)	
	
}else{
	plot(length(rankBy_vector):1,rankBy_vector[signalOrder], col='red',xlab=paste(rankBy_factor,' Stitched peaks'),ylab=paste(rankBy_factor,' Signal','- ',wceName),pch=19,cex=2)
}
abline(h=cutoff_options$absolute,col='grey',lty=2)
abline(v=length(rankBy_vector)-length(superEnhancerRows),col='grey',lty=2)
lines(length(rankBy_vector):1,rankBy_vector[signalOrder],lwd=4, col='red')
text(0,0.8*max(rankBy_vector),paste(' Cutoff used: ',cutoff_options$absolute,'\n','Super-Stitched peaks identified: ',length(superEnhancerRows)),pos=4)

dev.off()




#Writing a bed file
bedFileName = paste(outFolder,enhancerName,'_Stitched_withSuper.bed',sep='')
convert_stitched_to_bed(stitched_regions,paste(rankBy_factor,"Stitched"), enhancerDescription,bedFileName,score=rankBy_vector,splitSuper=TRUE,superRows= superEnhancerRows,baseColor="0,0,0",superColor="255,0,0")



#This matrix is just the super_enhancers
true_super_enhancers <- stitched_regions[superEnhancerRows,]

additionalTableData <- matrix(data=NA,ncol=2,nrow=nrow(stitched_regions))
colnames(additionalTableData) <- c("stitchedPeakRank","isSuper")
additionalTableData[,1] <- nrow(stitched_regions)-rank(rankBy_vector,ties.method="first")+1
additionalTableData[,2] <- 0
additionalTableData[superEnhancerRows,2] <- 1


#Writing enhancer and super-enhancer tables with enhancers ranked and super status annotated
enhancerTableFile = paste(outFolder,enhancerName,'_AllStitched.table.txt',sep='')
writeSuperEnhancer_table(stitched_regions, enhancerDescription,enhancerTableFile, additionalData= additionalTableData)

superTableFile = paste(outFolder,enhancerName,'_SuperStitched.table.txt',sep='')
writeSuperEnhancer_table(true_super_enhancers, enhancerDescription,superTableFile, additionalData= additionalTableData[superEnhancerRows,])
