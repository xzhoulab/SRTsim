
jscode <- "Shiny.addCustomMessageHandler('closeWindow', function(m) {window.close();});"

sidebar <- dashboardSidebar(
  tags$style(".left-side, .main-sidebar {padding-top: 65px}"),
  tags$head(tags$script(HTML(jscode))),
  sidebarMenu(      
    menuItem("Spatial Pattern Design",tabName = "dashboard",icon = icon("pencil-ruler")),
    menuItem("Count Data Generation",tabName = "simCounts", icon = icon("laptop-code")),

    fluidRow(
      column(width=12,
        column(width=4,actionButton("exit", "Exit",icon=icon("sign-out-alt"))),
       
        # column(width=4,actionButton("close", "Close window",icon=icon("sign-out-alt"))),
        column(width=4,actionButton("reload", "Restart",icon=icon("redo")))
      )
    )
  ) ## end of the sidebarMenu
)## end of the dashboardSidebar
