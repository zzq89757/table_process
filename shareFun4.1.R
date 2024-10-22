#### File for store all functions shared by variety test analysis.
#### modified and coded by Guoqi Yao 2017.11.24

#### Function to do anova across locations, return anov summary information.
## return summary of anova and means caculated with NA replaced by least square sum estimation.
hz.aov_tot <- function(df, column="variety"){
    # df: a breeding yield test dataframe with head "location, variety, rep1, rep2, ....	
    library(reshape)
    df.2 <- melt(df,id=c("location", "variety"))
    colnames(df.2)[3:4]<-c("rep", "y")
    ##ANOVA
    aov.res <- aov(y ~ location * variety + rep %in% location, data=df.2)
    l.aov.res.sum <-summary(aov.res)
    ##Calulate number of yield data for each variety or location. 
    if (column == "variety") {
	df.3 <- cast(df.2, location + rep ~ variety ,value = "y")
    }else{
	df.3 <- cast(df.2, variety + rep ~ location, value = "y")
    }
    r <- !is.na(df.3) 
    r <- r * 1
    obj_col <- colSums(r, na.rm=T)
    obj_col <- obj_col[-c(1, 2)]
    ##calculate mean yield with NA replaced by least square sum estimation.
    df_repNA <- df
    if(nMissdata(df_repNA) > 0) df_repNA <- lstSqEstMissdata(df_repNA)
    df_repNA_melt <- melt(df_repNA, id=c("location", "variety"))
    colnames(df_repNA_melt)[3:4] <- c("rep", "y")
    if (column == "variety") {
	df_repNA.3 <- cast(df_repNA_melt, location + rep ~ variety ,value = "y")
    }else{
	df_repNA.3 <- cast(df_repNA_melt, variety + rep ~ location, value = "y")
    }
    mean_col <- colMeans(df_repNA.3)
    meanAndCnt <- cbind(mean_col, obj_col)
    meanAndCnt <- as.data.frame(meanAndCnt)
    meanAndCnt <- meanAndCnt[order(meanAndCnt$mean_col, decreasing=T), ]
    return(list(aov_res=l.aov.res.sum, meanAndCnt=meanAndCnt, dataRepNA=df_repNA ))
}


#### Function to do local anova, return anov summary information.
hz.aovLoc <-function(df)
{
# df: a breeding yield test dataframe with head of "location, variety, rep1, rep2, ....	
	library(reshape)
	df.2 <- melt(df[,-1], id=c("variety"))
	colnames(df.2)[2:3] <- c("rep","y")
	aov.res<-aov(y ~ variety + rep, data=df.2)
	## Caculate varieties' mean yield.
	if(nMissdata(df) > 0){
	    df <- lstSqEstMissdata(df)
	    df.2 <- melt(df[,-1], id=c("variety"))
	    colnames(df.2)[2:3] <- c("rep","y")
	}
	n.meanYield <- mean(df.2$y)
	df.aov.res.sum <-summary(aov.res)[[1]]
	return(list(aov.result=df.aov.res.sum, meanYield = n.meanYield))
}


####Function to estimate missing data with Least Square
## Estimate missing data from linear model. Recive a dataframe, return a new dataframe 
##    with miss data replaced.
lstSqEstMissdata <- function(df_fieldTest){
    #df - a yield test dataframe across all locations with structure of " locaion, variety, rep1yield, rep2yield ... 
    library(reshape2)
    uni_loc=unique(df_fieldTest$location)
    for (n in 1:length(uni_loc)){
	df <- df_fieldTest[df_fieldTest$location == uni_loc[n],]
	names_col <- names(df)
	names_loc <- df[,1]
	names_var <- df[,2]
	#rownames(df) <- names_var
	df_nonLoc <- df[,-1]
	df_nonLoc_melt <- melt(df_nonLoc, id="variety")
	names(df_nonLoc_melt) <- c("variety", "rep", "yield")
	result_lm <- lm(yield ~ variety + rep, data=df_nonLoc_melt)
	#get the index for all missing data
	mt_df <- df[,-c(1,2)]
	r <- is.na(mt_df)
	x <- r * 1
	ind_miss <- which(r == 1)
	if(length(ind_miss) > 0){ 
	    colInd_miss <- ceiling(ind_miss / nrow(r))
	    rowInd_miss <- ind_miss %% nrow(r) 
	    rowInd_miss[rowInd_miss == 0] <- nrow(r)
	    for(i in 1:length(ind_miss)){
		est_missdata <- result_lm$coefficients[1]
		row_index <- rowInd_miss[i]
		if(row_index > 1) est_missdata <- est_missdata + result_lm$coefficients[row_index]
		col_index <- colInd_miss[i]
		if (col_index > 1) est_missdata <- est_missdata + result_lm$coefficients[col_index - 1 + nrow(r)]
		mt_df[row_index, col_index] <- est_missdata
	    }
	    df <- cbind.data.frame(names_loc, names_var, mt_df)
	    names(df) <- names_col
	}
	if (n == 1) new_dfFieldTest <- df
	else new_dfFieldTest <- rbind.data.frame(new_dfFieldTest, df)
    }
    return (new_dfFieldTest)
}


####Function to count missing data number 
nMissdata <- function(df_fieldTest){
    #df - a yield test dataframe across all locations with structure of " locaion, variety, rep1yield, rep2yield ... 
    mt_df <- df_fieldTest[,-c(1,2)]
    r <- is.na(mt_df)
    x <- r * 1
    return(sum(x))
}
####Function convert Zero to NA
ZeroToNA <- function(df){
    #df - a yield test dataframe across all locations with structure of " locaion, variety, rep1yield, rep2yield ... 
    df[-c(1, 2)][df[-c(1, 2)]==0] <- NA
    return(df)
}

####Function convert joint yield test data with multirepeat to meanyiead table across variey * enviroment.
####2017.12.7
meanYd.VarLoc <- function(yieldData)
#yieldData: a dataframe with head and structure of "location, variety, rep1, rep2. Missindata should been coded with NA.
#return value: a meanyield dataframe with head of "variety, envirnoment1, enviroment2, ...
{
    require(reshape2)
    ## estimate missingdata.  
    if(nMissdata(yieldData) > 0) yieldData <- lstSqEstMissdata(yieldData)
    ## caculate average yield across each repeat in a location.
    yieldData_mat <- yieldData[,-c(1:2)]
    repMeans <- rowMeans(yieldData_mat)
    rowMeans_df <- data.frame(location=yieldData$location, variety=yieldData$variety, repMeanYd=repMeans)
    meanYdVarLoc <- dcast(rowMeans_df, variety~location)
    rownames(meanYdVarLoc) <- meanYdVarLoc$variety
    meanYdVarLoc$variety <- NULL
    return(meanYdVarLoc)
}
##########END######################


