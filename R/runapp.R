#' Run the SRTsim Shiny Application
#' @importFrom shiny shinyApp
#' @return A list that contains a count matrix, a location dataframe, and all parameter specifications.
#' @export
SRTsim_shiny <- function(){
  # vv = NULL
  theme_set(theme_bw())
  out <- runApp(shinyApp(ui,server))
  return(out)
}
