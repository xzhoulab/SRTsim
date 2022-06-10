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

setMethod("show", "simSRT",
    function(object)
{
    cat("class:", class(object), "\n")
    cat("Reference name:", object@refID, "\n")
    cat("Reference dim:", dim(object@refCounts), "\n")

    obj_refcounts <- object@refCounts
    
    ## metadata()
    # expt <- names(metadata(object))
    # if (is.null(expt))
    #     expt <- character(length(metadata(object)))
    # coolcat("metadata(%d): %s\n", expt)

    # ## assays()
    # nms <- assayNames(object)
    # if (is.null(nms))
    #     nms <- character(length(assays(obj_assay, withDimnames=FALSE)))
    # coolcat("assays(%d): %s\n", nms)

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

# setGeneric("colData", function(x, ...) standardGeneric("colData"))

# ## Fix old DataFrame instances on-the-fly.
# setMethod("colData", "SRT",
#     function(x, ...) updateObject(x@colData, check=FALSE)
# )

# setGeneric("counts", function(x, ...) standardGeneric("counts"))
# setMethod("counts", "SRT",
#     function(x, ...) updateObject(x@obj_refcounts, check=FALSE)
# )


#' Access reference count matrix
#' @param x SRTsim object
#' @export
setGeneric("refCounts", function(x) x@refCounts)

#' Access reference colData 
#' @param x SRTsim object
#' @export
setGeneric("refcolData", function(x) x@refcolData)

#' Access reference rowData 
#' @param x SRTsim object
#' @export
setGeneric("refrowData", function(x) x@refrowData)


#' Access synthetic count matrix 
#' @rdname simCounts
#' @param x SRTsim object
#' @export
setGeneric("simCounts", function(x) x@simCounts)

#' Access synthetic colData 
#' @param x SRTsim object
#' @export
setGeneric("simcolData", function(x) x@simcolData)

#' Access synthetic rowData 
#' @param x SRTsim object
#' @export

setGeneric("simrowData", function(x) x@simrowData)


#' Access Model Fitting Parameters 
#' @param x SRTsim object
#' @export

setGeneric("EstParam", function(x) x@EstParam)

#' Access User-Specified Parameters
#' @param x SRTsim object
#' @export
setGeneric("metaParam", function(x) x@metaParam)

