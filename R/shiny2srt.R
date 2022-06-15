#' Create a SRTsim object from reference-free shinyoutput
#' @param shinyOutput A list of Shiny Output. Including a simCount, simInfo,simcountParam,simLocParam
#' @return Returns a SRTsim object with user-specified parameters stored in metaParam slot.
#' @importFrom S4Vectors SimpleList
#' @importFrom methods as
#' @export
#' @examples
#'
#' shinySRT <- Shiny2SRT(toyShiny)
#'
#' ## Explore the new SRT object
#' shinySRT@metaParam
#' shinySRT@simCounts[1:3,1:3]
#' shinySRT@simcolData

Shiny2SRT <- function(shinyOutput){
    simsrt <- new(
        Class = "simSRT",
        simCounts = as(as.matrix(shinyOutput$simCount),"sparseMatrix"),
        # simCounts=matrix(3,3,3),
        simcolData = as(shinyOutput$simInfo,"DataFrame"),
        metaParam = list(shinySRTParam = S4Vectors::SimpleList(
        	simcountParam = shinyOutput$simcountParam,
        	simLocParam= shinyOutput$simLocParam),
        	simType = "shinySRT")
    )
    return(simsrt)
}



