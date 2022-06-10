
#' Subset SRT object based on domain labels of interest
#' @param simsrt A SRTsim object
#' @param sel_label    A vector of selected domain labels used for the data generation
#' @return Returns a spatialExperiment-based object 
#' 
#' @export 

subsetSRT <- function(simsrt,sel_label=NULL){
    if(is.null(sel_label)){
        cat("## No specific label provided in sel_label, original object returned\n")
    }else{
        locfile    <- simsrt@refcolData
        countfile  <- simsrt@refCounts
    
        subLoc     <- locfile[locfile$label %in% sel_label ,]
        subCount   <- countfile[,rownames(subLoc)]

        simsrt@refcolData <- subLoc
        simsrt@refCounts  <- subCount 
        cat("## ",nrow(locfile)-nrow(subLoc),"spatial locations are removed \n")
    }

	return(simsrt)
}
