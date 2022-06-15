#' Data used for creating vignettes 
#'
#' A data list containing the a gene expression matrix and a location dataframe
#'
#' @format A data list
#' \describe{
#'   \item{example_count}{A sparse matrix with 80 rows and 3611 columns.}
#'   \item{example_loc}{A data frame with 3611 rows and 6 columns.}
#' }
#' @source {created based on a SpatialLIBD data (SampleID: 151673) to serve as an example}
#' @examples 
#' data(exampleLIBD)       #Lazy loading. Data becomes visible as soon as called
"exampleLIBD"




#' A toyExample to showcase reference-based simulations 
#'
#' A data list containing the a gene expression matrix and a location dataframe
#'
#' @format A data list
#' \describe{
#'   \item{toyCount}{A sparse matrix with 100 rows and 251 columns.}
#'   \item{toyInfo}{A data frame with 251 rows and 3 columns.}
#' }
#' @source {created based on a ST Human Breast Cancer data to serve as an example}
#' @examples 
#' data(toyData)       #Lazy loading. Data becomes visible as soon as called
"toyData"



#' A toyExample to showcase reference-free simulations
#'
#' A list of shiny output 
#'
#' @format A list of shiny output 
#' \describe{
#'   \item{simCount}{A data frame with 150 rows and 980 columns.}
#'   \item{simInfo}{A data frame with 980 rows and 4 columns: x, y, group, foldchange}
#'   \item{simcountParam}{A list of user-specified parameters for count generation}
#'   \item{simLocParam}{A list of user-specified parameters for pattern design}
#' }
#' @source {created based using the SRTsim_shiny()}
#' @examples 
#' data(toyShiny)       #Lazy loading. Data becomes visible as soon as called
"toyShiny"









