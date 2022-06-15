
#' Subset SRT object based on domain labels of interest
#' @param simsrt A SRTsim object
#' @param sel_label    A vector of selected domain labels used for the data generation
#' @return Returns a spatialExperiment-based object 
#' 
#' @export 
#' @examples
#'
#' ## Create a simSRT object
#' toySRT  <- createSRT(count_in=toyData$toyCount,loc_in = toyData$toyInfo)
#' set.seed(1)
#' ## Only Keep the Spatial Locations labelled as "A" in the reference data
#' subtoySRT <- subsetSRT(toySRT,sel_label = "A")
#'

subsetSRT <- function(simsrt,sel_label=NULL){
    if(is.null(sel_label)){
        message("## No specific label provided in sel_label, original object returned")
    }else{
        locfile    <- simsrt@refcolData
        countfile  <- simsrt@refCounts
    
        subLoc     <- locfile[locfile$label %in% sel_label ,]
        subCount   <- countfile[,rownames(subLoc)]

        simsrt@refcolData <- subLoc
        simsrt@refCounts  <- subCount 
        message("## ",nrow(locfile)-nrow(subLoc),"spatial locations are removed")
    }

	return(simsrt)
}
