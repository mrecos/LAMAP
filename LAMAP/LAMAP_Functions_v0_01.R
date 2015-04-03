### LAMAP v0.01 PRE-ALPHA
### 
### TODO: results function; exists, but make sure it works alright and includes the right metrics
### TODO: data preperation
### TODO: test to see if scale(varaibles) and using SD for shift is a good idea

# Required packages
require(np)
require(raster)
require(rgeos)
require(ROCR)
# end packages 

roundingFactor <- 1 # significant digits of variabel measures before computation of CDFs
siteNeighborhood <- 10 # number of nearest sites for conditional probability of Loi
distKernel <- "exponential" # one of exponential or uniform (default)
shiftValues <- matrix(c(1,2), nrow=1) # 1xp matrix, p = variables, row values are shift value
colnames(shiftValues) <- c("X1", "X2")  # colnames are varaible names, must match trainingDf

# step 1, send training data to createSiteCDF(.) to return list of CDFs, one for each site
allSitesCDF <- createSiteCDF(trainingDf, roundingFactor)
save(allSitesCDF,file="./allSitesCDF.RData") # save	
# step 2, send allSitesCDF and other objects to buildModel(.), returns prediction raster
# set shiftValues as matrix of shift per variable
predRast <- buildModel(allSitesCDF, vectPoints, r_stack, rasterMask, distKernel="uniform", roundingFactor, shiftValues)
# step 3, get prediction results, polys can be train or test sites depending on intention
results <- lamapResults(predRaster, polys, sampleMultiplier = 10)

# function to check for zero-variance in variable measurements at a given site (cat)
constantCheck <- function(trainingDf, colsSelect){
  siteNames <- unique(trainingDf$cat) # unique sites
  for(i in seq_along(siteNames)){
    d <- trainingDf[which(trainingDf$cat == siteNames[i]),]
    uniqueLength <- apply(d[,colsSelect],2,function(x) length(unique(x)))
    for(j in length(uniqueLength)){
      if(1 %in% uniqueLength[j]){
        value <- trainingDf[which(trainingDf$cat == siteNames[i]),colsSelect[j]][1]
        trainingDf[which(trainingDf$cat == siteNames[i]),colsSelect[j]][1] <- value + 0.1 # arbitrary adjustment
      }
    }
  }
  return(trainingDf)
}

createSiteCDF <- function(trainingDf, roundingFactor){
	colsSelect <- colnames(trainingDf)[!(colnames(trainingDf) %in% c("cat","x","y"))] #hardcode these???
	trainingDf[,colsSelect] <- round(trainingDf[,colsSelect],roundingFactor) # reduce percision of data 
	trainingDf <- constantCheck(trainingDf, colsSelect) # check for constant values on sites, add 0.1 if so
	siteNames <- unique(trainingDf$cat)
	allSitesCDF <- list() # was npu_pdf_bw, changed name to allSitesCDF
	for (j in seq_along(siteNames)) { # For each site, filter data, create CDF, add to list
		siteData <- trainingDf[which(trainingDf$cat == siteNames[j]), colsSelect]
		CDF <- np::npudistbw(siteData) # create non-parametric cumulative distribution function (CDF)
		allSitesCDF[[siteNames[j]]] <- list(siteNames[j],CDF,siteData)
	}
	# allSitesCDF <- allSitesCDF[!sapply(allSitesCDF,is.null)] ## test this a bit, not bad to have NULL checker
return(allSitesCDF)
}

buildModel <- function(allSitesCDF, vectPoints, r_stack, rasterMask, distKernel="uniform", roundingFactor, shiftValues, siteNeighborhood) {
	print(paste('Process began at',Sys.time()))
	print('Calculating probabilities...')
	modelRastDf <- as.data.frame(rasterMask)
	numcells <- length(modelRastDf[,1])
	modelRastDf['index'] <- (seq(1,numcells,1)) # index is cell number to assign value to
	modelRastDf <- cbind(modelRastDf, coordinates(rasterMask))
	modelRastDf['probability'] <- 0
	modelRastDf <- modelRastDf[,-1] # remove "layer" pesky column  
	modelRastDfAll <- modelRastDf   # A copy of modelRastDf with all cells to attach results to
	modelRastDf <- modelRastDf[complete.cases(modelRastDf),]   ### remove NA cells
	# test this ----------------------------------------------
	colsSelect <- names(r_stack) ### unless passed in, this is not checked against allSitesCDF colSelect, it will fail if not equal !!! TEST !!!
	# end test -----------------------------------------------
	pb <- txtProgressBar(min = 0, max = length(modelRastDf$index), style = 3)
	for (i in 1:nrow(modelRastDf)) { 
		originContain <- modelRastDf[i,] # isolate cell in modelRast for which to compute probability
		rastValues <- NULL
		for(k in seq_along(colsSelect)){
		  xyValue <-  r_stack[[k]][originContain[,"index"]]
		  rastValues <- cbind(rastValues, xyValue)
		}
		colnames(rastValues) <- colsSelect
		originContain <- data.frame(originContain,rastValues)
		originContain[,colsSelect] <- round(originContain[,colsSelect],roundingFactor) # rounding
		# FACTOR CONTROL GOES HERE --------------
		modelRastDf[i,'probability'] <- calcProb(originContain, vectPoints, allSitesCDF, distKernel, colsSelect, shiftValues)
		setTxtProgressBar(pb, i)
	}
	close(pb)
	print(paste('Process complete at',Sys.time()))
	# dm_text <- paste(c(as.character(quantile(modelRastDf$probability, seq(0,1,0.1)))), collapse=',')
	# twitter_update(dm_text)
	predRaster <- createProbRaster(modelRastDfAll, modelRastDf)
	return(predRaster)
}

createProbRaster <- function(modelRastDfAll, modelRastDf){
	predRasterdf <- merge(modelRastDfAll, modelRastDf, by="index", all = TRUE, stringsAsFactors = FALSE)
	predRasterdf <- predRasterdf[,c("index", "x.x", "y.x", "probability.y")]
	colnames(predRasterdf) <- c("index", "x", "y", "probability")
	predRaster <- raster(rasterMask)
	predRaster[] <- predRasterdf$probability
return(predRaster)
}

calcProb <- function(originContain, vectPoints, allSitesCDF, distKernel, colsSelect, shiftValues, siteNeighborhood) {
	siteDistance <- data.frame(raster::pointDistance(originContain[c('x','y')], vectPoints[,c('x','y')], lonlat=FALSE)) #set up for UTM, add switch for Lat/Lon
	siteDistance$cat <- vectPoints$cat
	colnames(siteDistance) <- c("distance","cat")
	# import colsSelect from buildModel(.), check to make sure this works...
	# colsSelect <- colnames(originContain)[!(colnames(originContain) %in% c("cat","x","y","index","probability"))] # original line
	originMin <- originContain
	originMax <- originContain
	for (j in colsSelect) {
		colValue <- originContain[,j]
		# shift <- NULL
		# shift <- ifelse(j == "X1", 300, 1.0) # cludgy hard code... ### NEEDS ATTENTION!!! pass in as matrix of shift values
		originMin[,j] <- colValue - shiftValues[1,j] 
		originMax[,j] <- colValue + shiftValues[1,j] 
	}
	predCDF <- data.frame(cat = NA, probability = NA, stringsAsFactors = FALSE)
	for (i in seq_along(allSitesCDF)) { # switched to npudist
		sitesBWS <- allSitesCDF[[i]][[2]]
		fitMin <- np::npudist(bws=sitesBWS,edat=originMin[colsSelect])
		probMin <- fitted(fitMin)
		fitMax <- np::npudist(bws=sitesBWS,edat=originMax[colsSelect])
		probMax <- fitted(fitMax)
		prob <- probMax - probMin # distance between min and max is prob of intersecting range of site CDF
		predCDF[i,] <- c(allSitesCDF[[i]][[1]],round(prob,5))
	}
# 	predCDF <- predCDF[-which(apply(predCDF,1,function(x)all(is.na(x)))),] # NA Checker, make somthing that works
	loiSitePred <- merge(siteDistance, predCDF, by="cat", stringsAsFactors = FALSE)
	loiSitePred$probability <- as.numeric(loiSitePred$probability) # hack b/c it keeps going to chr
	loiSitePred <- loiSitePred[which(loiSitePred['probability'] > 0),] # remove zero prob 
	loiSitePred <- loiSitePred[order(loiSitePred$distance),] # order sites, nearest to furthest
	loiSitePred$weight <- getDistanceWeights(loiSitePred, distKernel) # uniform or exponential
	loiSitePred$weightedProb <-  round(loiSitePred$probability * loiSitePred$weight,3)
	if (length(loiSitePred[,'cat']) > siteNeighborhood){  
		loiSitePred <- loiSitePred[c(1:siteNeighborhood),]
	}
	if (length(loiSitePred[,'cat']) > 1) {
		p <- 1-prod(1-loiSitePred$weightedProb) # computes Union of probabilities; inclusion/exclusion principle
		return(round(p,3))
	} 
	else if (length(loiSitePred[,'cat']) == 0) { # in case of zero sites
		return(0)
	} 
	else { return(round(loiSitePred[,'probability'],3)) } # in case of 1 site
}

getDistanceWeights <- function(loiSitePred, distKernel){
	if( distKernel == "uniform"){
		distanceWeigths <- rep(1, nrow(loiSitePred))
	} else if( distKernel == "exponential"){
		distanceWeigths <- round((1-(1/loiSitePred$distance^(1/seq_along(loiSitePred$distance)))),3) # inverse decay weight 
	}
return(distanceWeigths)
}

lamapResults <- function(predRaster, polys, sampleMultiplier = 10){
  polyValue <- extract(predRaster, polys)
  names(polyValue) <- polys$SITENO
  # Site specific aggregate
  polyMean <- do.call(rbind, lapply(polyValue,mean))
  polyMedian <- do.call(rbind, lapply(polyValue,median))
  polyMax <- do.call(rbind, lapply(polyValue,max))
  polyMin <- do.call(rbind, lapply(polyValue,min))
  polyMetrics <- round(cbind(polyMean, polyMedian, polyMax, polyMin),3)
  colnames(polyMetrics) <- c("Mean", "Median", "Max", "Min")
  polyMetrics <- data.frame(polyMetrics, row.names = seq(1,nrow(polyMetrics)))
  polyMetrics$cat <- polys$SITENO
  # all Sites vs. background
  siteValues <- unlist(polyValue)
  backValues <- sampleRandom(predRaster, (nrow(siteValues) * sampleMultiplier), na.rm=TRUE)
  # auc 
  pred <- c(as.numeric(siteValues), backValues)
  actual <- c(rep(1,length(siteValues)), rep(0,length(backValues)))
  predROC <- ROCR::prediction(pred, actual) #predicted, and actaul as 1/0
  perf <- ROCR::performance(predROC, measure = "tpr", x.measure = "fpr") 
  auc <- ROCR::performance(predROC, "auc")
  auc <- round((auc@y.values[[1]]),3)
  return(list(auc, polyMetrics))
}

