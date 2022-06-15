#' Run the SRTsim Shiny Application
#' @importFrom shiny shinyApp runApp
#' @importFrom ggplot2 theme_set theme_bw
#' @return A list that contains a count matrix, a location dataframe, and all parameter specifications.
#' @export
#' @examples
#'
#'\dontrun{
#'  ## Will Load an Interactive Session
#' shinyOutput <- SRTsim_shiny()
#'}

SRTsim_shiny <- function(){
  # vv = NULL
  theme_set(theme_bw())
  out <- runApp(shinyApp(ui,server))
  return(out)
}
