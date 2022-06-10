#' @import shiny
#' @import shinydashboard
ui <- shinydashboard::dashboardPage(
    title="Spatially Resolved Transcriptomics Simulator",
    shinydashboard::dashboardHeader(
        title = span('SRTsim',style = "font-size: 30px; font-family:'Copperplate'; font-weight: bold"), 
        tags$li(
            class = "dropdown",
            tags$style(".main-header {max-height: 50px ;}"),
            tags$style(".main-header .logo {height: 50px;}"),
            tags$style(".sidebar-toggle {height: 50px; color: #FFFFFF !important}"),
            tags$style(".navbar {min-height:50px !important}"),
            tags$style(type = 'text/css',".badge{min-height: 20px;}")
            # ,
            # tags$script("document.getElementsByClassName('sidebar-toggle')[0].style.visibility = 'hidden';") ## hide the sidebar-toggle
        )
    ),
    sidebar,
    body
) ## end of ui (dashboardPage)
