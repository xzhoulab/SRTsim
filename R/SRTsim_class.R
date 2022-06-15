#' @importFrom Matrix colSums rowSums colMeans rowMeans
#' @importClassesFrom Matrix dgCMatrix
#' @importClassesFrom S4Vectors DataFrame SimpleList

setClassUnion(name = 'AnyMatrix', members = c("matrix", "dgCMatrix","NULL"))
setClassUnion(name = 'AnyDF', members = c("DataFrame", "NULL"))
setClassUnion(name = 'OptionalList', members = c('NULL', 'list'))
setClassUnion(name = 'OptionalCharacter', members = c('NULL', 'character'))

setClass("simSRT",
    representation(
        refCounts="AnyMatrix",             # Data -- e.g., list of matrices
        refcolData="DataFrame",            # columns and their annotations
        refrowData="DataFrame",            # columns and their annotations
        simCounts="AnyMatrix", 
        simcolData="AnyDF",
        simrowData="AnyDF",
        EstParam="OptionalList", 
        metaParam="OptionalList", 
        refID="OptionalCharacter",
        elementMetadata="DataFrame"
    ),
    prototype(
        refcolData=new("DFrame"),
        refrowData=new("DFrame"),
        elementMetadata=new("DFrame")
    )
)


#' @importFrom S4Vectors coolcat
setMethod("show", "simSRT",
    function(object)
{
    cat("class:", class(object), "\n")
    cat("Reference name:", object@refID, "\n")
    cat("Reference dim:", dim(object@refCounts), "\n")

    obj_refcounts <- object@refCounts
    
    ## rownames()
    rownames <- rownames(obj_refcounts)
    if (!is.null(rownames)) coolcat("Reference rownames(%d): %s\n", rownames)
    else cat("Reference rownames: NULL\n")

    # ## rowData`()
    # coolcat("rowData names(%d): %s\n", names(rowData(object, use.names=FALSE)))

    ## colnames()
    colnames <- colnames(obj_refcounts)
    if (!is.null(colnames)) coolcat("Reference colnames(%d): %s\n", colnames)
    else cat("Reference colnames: NULL\n")

    ## colData()
    coolcat("refcolData names(%d): %s\n", names(refcolData(object)))
})


#' Access reference count matrix
#' @param x SRTsim object
#' @export
#' @return Returns a reference count matrix
#' @examples
#'
#' ## Create a simSRT object
#' toySRT  <- createSRT(count_in=toyData$toyCount,loc_in = toyData$toyInfo)
#' set.seed(1)
#' ## Estimate model parameters for data generation
#' toySRT <- srtsim_fit(toySRT,sim_schem="tissue")
#' ## Generate synthetic data with estimated parameters
#' toySRT <- srtsim_count(toySRT)
#' refCounts(toySRT)[1:3,1:3]
#'
setGeneric("refCounts", function(x) x@refCounts)

#' Access reference colData 
#' @param x SRTsim object
#' @export
#' @return Returns the colData of reference data
#' @examples
#'
#' ## Create a simSRT object
#' toySRT  <- createSRT(count_in=toyData$toyCount,loc_in = toyData$toyInfo)
#' set.seed(1)
#' ## Estimate model parameters for data generation
#' toySRT <- srtsim_fit(toySRT,sim_schem="tissue")
#' ## Generate synthetic data with estimated parameters
#' toySRT <- srtsim_count(toySRT)
#' refcolData(toySRT)
#'
setGeneric("refcolData", function(x) x@refcolData)

#' Access reference rowData 
#' @param x SRTsim object
#' @export
#' @return Returns the rowData of reference data
#' @examples
#'
#' ## Create a simSRT object
#' toySRT  <- createSRT(count_in=toyData$toyCount,loc_in = toyData$toyInfo)
#' set.seed(1)
#' ## Estimate model parameters for data generation
#' toySRT <- srtsim_fit(toySRT,sim_schem="tissue")
#' ## Generate synthetic data with estimated parameters
#' toySRT <- srtsim_count(toySRT)
#' refrowData(toySRT)
#'
setGeneric("refrowData", function(x) x@refrowData)


#' Access synthetic count matrix 
#' @param x SRTsim object
#' @export
#' @return Returns a synthetic count matrix
#' @examples
#'
#' ## Create a simSRT object
#' toySRT  <- createSRT(count_in=toyData$toyCount,loc_in = toyData$toyInfo)
#' set.seed(1)
#' ## Estimate model parameters for data generation
#' toySRT <- srtsim_fit(toySRT,sim_schem="tissue")
#' ## Generate synthetic data with estimated parameters
#' toySRT <- srtsim_count(toySRT)
#' simCounts(toySRT)[1:3,1:3]
#'
setGeneric("simCounts", function(x) x@simCounts)

#' Access synthetic colData 
#' @param x SRTsim object
#' @export
#' @return Returns the colData of synthetic data
#' @examples
#'
#' ## Create a simSRT object
#' toySRT  <- createSRT(count_in=toyData$toyCount,loc_in = toyData$toyInfo)
#' set.seed(1)
#' ## Estimate model parameters for data generation
#' toySRT <- srtsim_fit(toySRT,sim_schem="tissue")
#' ## Generate synthetic data with estimated parameters
#' toySRT <- srtsim_count(toySRT)
#' simcolData(toySRT)
#'
setGeneric("simcolData", function(x) x@simcolData)

#' Access synthetic rowData 
#' @param x SRTsim object
#' @export
#' @return Returns the rowData of synthetic data
#' @examples
#'
#' ## Create a simSRT object
#' toySRT  <- createSRT(count_in=toyData$toyCount,loc_in = toyData$toyInfo)
#' set.seed(1)
#' ## Estimate model parameters for data generation
#' toySRT <- srtsim_fit(toySRT,sim_schem="tissue")
#' ## Generate synthetic data with estimated parameters
#' toySRT <- srtsim_count(toySRT)
#' simrowData(toySRT)
#'
setGeneric("simrowData", function(x) x@simrowData)


#' Access Model Fitting Parameters 
#' @param x SRTsim object
#' @export
#' @return Returns a list of estimated parameters by fitting models
#' @examples
#'
#' ## Create a simSRT object
#' toySRT  <- createSRT(count_in=toyData$toyCount,loc_in = toyData$toyInfo)
#' set.seed(1)
#' ## Estimate model parameters for data generation
#' toySRT <- srtsim_fit(toySRT,sim_schem="tissue")
#' ## Generate synthetic data with estimated parameters
#' toySRT <- srtsim_count(toySRT)
#' EstParam(toySRT)
#'
setGeneric("EstParam", function(x) x@EstParam)

#' Access User-Specified Parameters
#' @param x SRTsim object
#' @export
#' @return Returns a list of user-specified parameters
#' @examples
#'
#' ## Create a simSRT object
#' toySRT  <- createSRT(count_in=toyData$toyCount,loc_in = toyData$toyInfo)
#' set.seed(1)
#' ## Estimate model parameters for data generation
#' toySRT <- srtsim_fit(toySRT,sim_schem="tissue")
#' ## Generate synthetic data with estimated parameters
#' toySRT <- srtsim_count(toySRT)
#' metaParam(toySRT)
#'
setGeneric("metaParam", function(x) x@metaParam)

