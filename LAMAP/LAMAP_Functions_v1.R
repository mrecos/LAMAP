### LAMAP 
### Trying to figure out of npudist or npudens.  I think dens, but so far all zero propbabilities (densities)

constantCheck <- function(training_df, cols_select){
  groups <- unique(training_df$cat)
  for(i in seq_along(groups)){
    d <- training_df[which(training_df$cat == groups[i]),]
    uniqueLength <- apply(d[,cols_select],2,function(x) length(unique(x)))
    for(j in length(uniqueLength)){
      if(1 %in% uniqueLength[j]){
        value <- training_df[which(training_df$cat == groups[i]),cols_select[j]][1]
        training_df[which(training_df$cat == groups[i]),cols_select[j]][1] <- value + 0.1
      }
    }
  }
  return(training_df)
}

require(np)
require(raster)
require(rgeos)

# training_df <- backup
# training_df[,1:2] <- scale(training_df[,1:2]) # might be worth a try, but need to do to raster as well
rounding_factor <- 1
predRast <- buildModel(NULL, vect_pts, training_df, r_stack, MASK, rounding_factor)

buildModel <- function(outrastname,vect_pts,training_df,r_stack, MASK, distkernel="uniform", rounding_factor) {
	print(paste('Process began at',Sys.time()))
	rast_list <- r_stack
	rasts <- names(r_stack)
	##prepare grid for holding probabilities
	print('Extracting raster MASK as probability grid from GRASS')
# 	modelRast <- raster(MASK)
# 	projection(modelRast) <- CRS('proj4string') # need to complete
	cols_select <- colnames(training_df)[!(colnames(training_df) %in% c("cat","x","y"))]
  training_df[,cols_select] <- round(training_df[,cols_select],rounding_factor) # reduce percision of data 
  training_df <- constantCheck(training_df, cols_select) # check for constant values on sites, add 0.1 if so
  siteNames <- unique(training_df$cat)
	npu_pdf_bw <- list() # empty list
	for (j in seq_along(siteNames)) { # For each site
    siteData <- training_df[which(training_df$cat == siteNames[j]), cols_select]
		tmp <- np::npudistbw(siteData)
		npu_pdf_bw[[siteNames[j]]] <- list(siteNames[j],tmp,siteData)
	}
# 	npu_pdf_bw <- npu_pdf_bw[!sapply(npu_pdf_bw,is.null)] ## tried this with NULL examples, but could not make it do much; perhaps bc the list structure is crazy
	save(npu_pdf_bw,file="./npu_pdf_bw.RData") # save
	print('Calculating probabilities...')
  modelRast_df <- as.data.frame(MASK)
  numcells <- length(modelRast_df[,1])
  modelRast_df['index'] <- (seq(1,numcells,1)) # index is cell number to assign value to
  modelRast_df <- cbind(modelRast_df, coordinates(MASK))
  modelRast_df['probability'] <- 0
  modelRast_df <- modelRast_df[,-1] # remove "layer" pesky column  
  modelRast_dfAll <- modelRast_df   # A copy of modelRast_df with all cells to attach results to
  modelRast_df <- modelRast_df[complete.cases(modelRast_df),]   ### THIS GETS RID OF NA Cells!!!
  pb <- txtProgressBar(min = 0, max = length(modelRast_df$index), style = 3)

print(paste('Process began at',Sys.time()))
pb <- txtProgressBar(min = 0, max = 2500, style = 3)  
for (i in 1:nrow(modelRast_df)) { # for (i in 1:2500) { #   
		##isolate cell in modelRast for which to compute probability
		origin_contain <- modelRast_df[i,]
		origin_xycoords <- matrix(unlist(c(origin_contain['x'],origin_contain['y'])), nrow=1)
		##set up names for origin dataframe that will be passed to calcProb()
		##sample data rasters at origin_xycoords, fill origin dataframe, and rename columns to match training data dataframes
    rastValues <- NULL
    for(k in seq_along(cols_select)){
      value <-  r_stack[[k]][origin_contain[,"index"]] # extract(r_stack[[k]],origin_xycoords)
      rastValues <- cbind(rastValues, value)
    }
    colnames(rastValues) <- cols_select
		origin_contain <- data.frame(origin_contain,rastValues)
    ## ROUNDING
    origin_contain[,cols_select] <- round(origin_contain[,cols_select],rounding_factor)
    ## FACTOR CONTROL GOES HERE -- I took it out for testing
    modelRast_df[i,'probability'] <- calcProb(origin_contain,vect_pts,npu_pdf_bw,distkernel="Uniform")
    setTxtProgressBar(pb, i)
	}
  close(pb)
print(paste('Process end at',Sys.time()))

dm_text <- paste(c(as.character(quantile(modelRast_df$probability, seq(0,1,0.1)))), collapse=',')
twitter_update(dm_text)

predRasterdf <- merge(modelRast_dfAll, modelRast_df, by="index", all = TRUE, stringsAsFactors = FALSE)
predRasterdf <- predRasterdf[,c("index", "x.x", "y.x", "probability.y")]
colnames(predRasterdf) <- c("index", "x", "y", "probability")
predRaster <- raster(MASK)
predRaster[] <- predRasterdf$probability
writeRaster(predRaster, filename = "../output/raster/R91RS3_test1.tif", format="GTiff", datatype="FLT4S")
#   plot(predRaster)
#   plot(polys, add=TRUE) 
results <- lamapResults(predRaster, polys, 10)
print("Writing raster map...")
return(list(predRaster, results))
}

calcProb <- function(origin_contain,vect_pts,npu_pdf_bw,distkernel) {
  dist_ <- data.frame(raster::pointDistance(origin_contain[c('x','y')],vect_pts[,c('x','y')], lonlat=FALSE)) # works
	dist_$cat <- vect_pts$cat
  colnames(dist_) <- c("dist","cat")
	cols_select <- colnames(origin_contain)[!(colnames(origin_contain) %in% c("cat","x","y","index","probability"))]
	origin_min <- origin_contain
	origin_max <- origin_contain
	for (j in cols_select) {
    colValue <- origin_contain[,j]
    shift <- NULL
    shift <- ifelse(j == "X1", 300, 1.0) # cludgy hard code...
		origin_min[,j] <- colValue - shift 
		origin_max[,j] <- colValue + shift 
	}
	pdf_pred <- data.frame(cat = NA, probability = NA, stringsAsFactors = FALSE)
	for (i in seq_along(npu_pdf_bw)) { # switched to npudist
	  sitesBWS <- npu_pdf_bw[[i]][[2]]
    fit_min <- np::npudist(bws=sitesBWS,edat=origin_min[cols_select])
		prob_min <- fitted(fit_min)
    fit_max <- np::npudist(bws=sitesBWS,edat=origin_max[cols_select])
		prob_max <- fitted(fit_max)
		prob <- prob_max - prob_min
		pdf_pred[i,] <- c(npu_pdf_bw[[i]][[1]],round(prob,8))
	}
# 	pdf_pred <- pdf_pred[-which(apply(pdf_pred,1,function(x)all(is.na(x)))),] # breaks things
	tmp <- merge(dist_,pdf_pred,by="cat",stringsAsFactors = FALSE)
  tmp$probability <- as.numeric(tmp$probability) # hack b/c it keeps going to chr
  tmp <- tmp[which(tmp['probability'] > 0),]
	tmp <- tmp[order(tmp$dist),]
  tmp$weight <- round((1-(1/tmp$dist^(1/seq_along(tmp$dist)))),3) # i came up with this, for inverse decay 
  tmp$weightedProb <-  round(tmp$probability * tmp$weight,3)
	if (length(tmp[,'cat']) > 10){  
		tmp <- tmp[c(1:10),]
	}
	if (length(tmp[,'cat']) > 1) {
#     p <- 1-prod(1-tmp$probability) # computes Union, unweighted
	  p <- 1-prod(1-tmp$weightedProb)
    return(round(p,3))
	} # in case of zero sites
	else if (length(tmp[,'cat']) == 0) {
		return(0.0)
	} # in case of 1 site
	else {return(round(tmp[,'probability'],3))}
}

lamapResults <- function(predRaster, polys, rSampleSize = 10){
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
  backValues <- sampleRandom(predRaster, rSampleSize, na.rm=TRUE)
  # auc 
  pred <- c(as.numeric(siteValues), backValues)
  actual <- c(rep(1,length(siteValues)), rep(0,length(backValues)))
#   predValues <- cbind(pred, actual)
  require(ROCR)
  predROC <- prediction(pred, actual) #predicted, and actaul as 1/0
  perf <- performance(predROC, measure = "tpr", x.measure = "fpr") 
  auc <- performance(predROC, "auc")
  auc <- round((auc@y.values[[1]]),3)
  return(list(auc, polyMetrics))
}

