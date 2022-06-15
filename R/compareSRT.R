
#' Summarize metrics for reference data and synthetic data
#' @param simsrt A SRTsim object
#' @return Returns an object with summarized metrics
#' @importFrom S4Vectors DataFrame
#' @export
#' @examples
#'
#' ## Create a simSRT object
#' toySRT  <- createSRT(count_in=toyData$toyCount,loc_in = toyData$toyInfo)
#' set.seed(1)
#' ## Estimate model parameters for data generation
#' toySRT <- srtsim_fit(toySRT,sim_schem="tissue")
#' ## Generate synthetic data with estimated parameters
#' toySRT <- srtsim_count(toySRT)
#' ## Compute metrics 
#' toySRT   <- compareSRT(toySRT)


compareSRT <- function(simsrt){
	tmp_ref <- DataFrame(get_stats_loc(simsrt@refCounts,group="Reference"))
	simsrt@refcolData <- cbind(simsrt@refcolData,tmp_ref)
	simsrt@refrowData <- DataFrame(get_stats_gene(simsrt@refCounts,group="Reference"))
	rm(tmp_ref)

	tmp_sim <- DataFrame(get_stats_loc(simsrt@simCounts,group="SRTsim"))
	simsrt@simcolData <- cbind(simsrt@simcolData,tmp_sim)
	simsrt@simrowData <- DataFrame(get_stats_gene(simsrt@simCounts,group="SRTsim"))

	rm(tmp_sim)
	simsrt@metaParam$compared <- TRUE
	return(simsrt)
}


#' Summarize gene-wise summary metrics 
#' @param mat  A count matrix
#' @param group  A group label
#' @param log_trans  A logical constant indicating whether to log transform the gene mean and variance
#' @return Returns a n by 5 dataframe with location metrics
#' @importFrom Matrix colSums rowSums
get_stats_gene <- function(mat, group, log_trans = TRUE){
  mean  <- Matrix::rowMeans(mat)
  var   <- apply(mat,1,var)
  cv    <- sqrt(var)/mean
  zero_gene <- Matrix::rowSums(mat < 1e-5)/ncol(mat)

  if(log_trans){
    mean 	<- log10(mean + 1)
    var 	<- log10(var + 1)
  }
  
  rowStat <- cbind.data.frame(sampleLabel=group,GeneMean=mean,GeneVar=var,GeneCV=cv,GeneZeroProp=zero_gene)
  return(rowStat)
}


#' Summarize location-wise summary metrics 
#' @param mat  A count matrix
#' @param group  A group label
#' @param log_trans  A logical constant indicating whether to log transform the libsize
#' @return Returns a n by 3 dataframe with location metrics
#' @importFrom Matrix colSums

get_stats_loc <- function(mat, group, log_trans = TRUE){
  zero_cell <- Matrix::colSums(mat < 1e-5)/nrow(mat)
  libsize   <- Matrix::colSums(mat)
  if(log_trans){
    libsize <- log10(libsize + 1)
  }
  colStat <- cbind.data.frame(sampleLabel=group,LocZeroProp = zero_cell, LocLibSize = libsize)
  return(colStat)
}
