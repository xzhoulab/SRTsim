
#' Create simSRT object
#' @param count_in     A gene expression count \code{matrix}
#' @param loc_in       A location \code{dataframe} with colnames x,y,label
#' @param refID    A \code{character} reference sample identifier. Default = \code{ref1}.
#' @return Returns a spatialExperiment-based object 
#' @importFrom methods new
#' @export 
#' @examples
#'
#' ## Create a simSRT object
#' toySRT  <- createSRT(count_in=toyData$toyCount,loc_in = toyData$toyInfo)
#'
#' ## Explore the object
#' toySRT 


createSRT <- function(count_in, loc_in,refID="ref1"){
  	
    ## check the location index in two files
    if(!identical(colnames(count_in),rownames(loc_in))){
        stop("colnames in the count file is different from the rownames in the location file")
    }

    ## check the location index in two files
    if(ncol(loc_in)!=3){
        stop("location file should be organized in three columns: x,y and label!")
    }

    ## in case the columns are not named in this way
    colnames(loc_in) <- c("x","y","label")

	# srt <- SpatialExperiment(
	#     assay = list(counts=count_in), 
	#     colData = loc_in, 
	#     spatialCoordsNames = c("x", "y"),
	#     sample_id=refID
 #    )

    simsrt <- new(
        Class = "simSRT",
        refCounts = count_in,
        # simCounts=matrix(3,3,3),
        refcolData = as(loc_in,"DataFrame"),
        refID=refID
    )

	return(simsrt)
}
