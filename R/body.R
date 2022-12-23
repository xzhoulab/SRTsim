#' @import shiny
#' @import shinydashboard
#' @importFrom plotly plotlyOutput
#' @importFrom shinyBS bsModal
#' @importFrom dashboardthemes shinyDashboardThemeDIY
body <- shinydashboard::dashboardBody(
    tags$head(
        tags$style(
            HTML(
                "pre { white-space: pre-wrap; word-break: keep-all; word-wrap:break-word;} 
                #countbrush {
                color: red;
                background: yellow;
                font-family: 'Times New Roman', Times, serif;
                font-size: 12px;
                font-style: italic;
                }
                #brush3 {
                color: red;
                background: yellow;
                font-family: 'Times New Roman', Times, serif;
                font-size: 14px;
                font-weight: bold;
                }
                #brush4 {
                color: blue;
                background: orange;
                font-family: Arial, Helvetica, sans-serif;
                font-size: 14px;
                }
                "
            ),
        )
    ), ## end of the tags$head

dashboardthemes::shinyDashboardThemeDIY(
  ### general
  appFontFamily = "Garamond"
  ,appFontColor = "rgb(45,45,45)"
  # ,primaryFontColor = "rgb(15,15,15)"
    ,primaryFontColor = "rgb(255,255,255)"
  ,infoFontColor = "rgb(15,15,15)"
  # ,successFontColor = "rgb(15,15,15)"
   ,successFontColor = "rgb(255,255,255)"
  ,warningFontColor = "rgb(255,255,255)"
  # ,warningFontColor = "rgb(45,45,45)"
  ,dangerFontColor = "rgb(15,15,15)"
  ,bodyBackColor = "rgb(255,255,255)" ## background of the output

  ### header

  # ,logoBackColor = "#ff6633"
  # ,headerButtonBackColor = "#ff6633"
  # ,headerButtonBackColorHover = "#ff6633"
  # ,headerButtonIconColorHover = "#ff6633"
  # ,headerBackColor = "#ff6633"
  # ,headerBoxShadowColor = "#ff6633"

  ,logoBackColor = "rgb(76,76,255)"
  ,headerButtonBackColor = "rgb(76,76,255)"
  ,headerButtonBackColorHover = "rgb(76,76,255)"
  ,headerButtonIconColorHover = "rgb(76,76,255)"
  ,headerBackColor ="rgb(76,76,255)"
  ,headerBoxShadowColor = "rgb(76,76,255)"
  ,headerButtonIconColor = "rgb(220,220,220)"
  ,headerBoxShadowSize = "0px 0px 0px"

  ### sidebar
  ,sidebarBackColor = "rgb(240,240,240)"
  # ,sidebarBackColor = "rgb(255,240,240)"
  ,sidebarPadding = 10

  ,sidebarMenuBackColor = "transparent"
  ,sidebarMenuPadding = 0
  ,sidebarMenuBorderRadius = 0

  ,sidebarShadowRadius = "0px 0px 0px"
  ,sidebarShadowColor = "#dfdfdf"
  ,sidebarUserTextColor = "rgb(115,115,115)"
  ,sidebarSearchBackColor = "rgb(240,240,240)"
  ,sidebarSearchIconColor = "rgb(100,100,100)"
  ,sidebarSearchBorderColor = "rgb(220,220,220)"

  ,sidebarTabTextColor = "rgb(100,100,100)"
  ,sidebarTabTextSize = 14
  ,sidebarTabBorderStyle = "none"
  ,sidebarTabBorderColor = "none"
  ,sidebarTabBorderWidth = 0

  ,sidebarTabBackColorSelected = "rgb(230,230,230)"
  ,sidebarTabTextColorSelected = "rgb(0,0,0)"
  ,sidebarTabRadiusSelected = "0px"

  ,sidebarTabBackColorHover = "rgb(245,245,245)"
  ,sidebarTabTextColorHover = "rgb(0,0,0)"
  ,sidebarTabBorderStyleHover = "none solid none none"
  ,sidebarTabBorderColorHover = "rgb(200,200,200)"
  ,sidebarTabBorderWidthHover = 4
  ,sidebarTabRadiusHover = "0px"


  # ,boxBackColor = "rgb(248,248,248)"
  ,boxBackColor = "rgb(255,255,255)"
  ,boxBorderRadius = 5
  ,boxShadowSize = "none"
  ,boxShadowColor = ""
  ,boxTitleSize = 18
  ,boxDefaultColor = "rgb(225,225,225)"
  # ,boxPrimaryColor = "rgb(95,155,213)"
  ,boxPrimaryColor = "#006c69"
  ,boxInfoColor = "rgb(180,180,180)"
    ,boxSuccessColor = "#194d7f"
  # ,boxSuccessColor = "#339bff"
  # ,boxSuccessColor = "rgb(102,102,185)"
  # ,boxSuccessColor = "rgb(112,173,71)"
  # ,boxWarningColor = "rgb(237,125,49)"
  ,boxWarningColor = "rgb(0,39,76)"
  ,boxDangerColor = "rgb(232,76,34)"

  ,tabBoxTabColor = "rgb(248,248,248)"
  ,tabBoxTabTextSize = 14
  ,tabBoxTabTextColor = "rgb(100,100,100)"
  ,tabBoxTabTextColorSelected = "rgb(45,45,45)"
  ,tabBoxBackColor = "rgb(248,248,248)"
  ,tabBoxHighlightColor = "rgb(200,200,200)"
  ,tabBoxBorderRadius = 1

  ### inputs
  ,buttonBackColor = "rgb(215,215,215)"
  ,buttonTextColor = "rgb(45,45,45)"
  ,buttonBorderColor = "rgb(150,150,150)"
  ,buttonBorderRadius = 5

  ,buttonBackColorHover = "rgb(190,190,190)"
  ,buttonTextColorHover = "rgb(0,0,0)"
  ,buttonBorderColorHover = "rgb(150,150,150)"

  ,textboxBackColor = "rgb(255,255,255)"
  ,textboxBorderColor = "rgb(118,118,118)"
  ,textboxBorderRadius = 5
  ,textboxBackColorSelect = "rgb(245,245,245)"
  ,textboxBorderColorSelect = "rgb(108,108,108)"

  ### tables
  # ,tableBackColor = "rgb(248,248,248)"
  ,tableBackColor = "rgb(255,255,255)"
  ,tableBorderColor = "rgb(238,238,238)"
  ,tableBorderTopSize = 1
  ,tableBorderRowSize = 1
),



    tabItems(
        tabItem("dashboard",
            fluidRow(
                column(width=8, 
                    box(title="Visualization",status="success",solidHeader=TRUE,width = NULL,
                    # align="middle",plotly::plotlyOutput('plot1',inline=FALSE),height=600),
                     # align="middle",plotly::plotlyOutput('plot1',inline=TRUE,width="auto",height="auto")
                        align="middle",plotly::plotlyOutput('plot1')
                     ),
                    br(),
                    fluidRow(
                    column(width=6,box(title="Interactive Selection Info",width=NULL,status="primary",solidHeader=TRUE,verbatimTextOutput('brush'))),
                    column(width=6,box(title="Group Summary", width=NULL,status="primary",solidHeader=TRUE,column(width=12,align="center",tableOutput('summary_table')))),

                    ),
                    verbatimTextOutput('brush3')
                ),

                column(width=4,
                    box(
                        width=NULL,solidHeader=TRUE,status="warning",
                        title="Choose CSV File",
                        height=140,
                        div(style = "margin-top:-20px"),
                        fluidRow(width = 10,
                          column(width=5,checkboxInput("header", "Header", TRUE)),
                          column(width=5,checkboxInput("rownames", "Row Names", TRUE))
                        ),
                        div(style = "margin-top:-50px"),
                        fileInput("datafile", label=NULL,
                          accept = c(
                            "text/csv",
                            "text/comma-separated-values,text/plain",
                            ".csv")
                        ),
                    ),

                    box(width=NULL,solidHeader=TRUE,status="warning",title="Shape and Layout",
                        selectInput('shape', 'Shape', c("Square","Circle","User Define Shape","User Define Spots")),
                        fluidRow(
                          column(width=6,selectInput('arrange', 'Spots Layout', c("Random","Grid"))),  
                          column(width=6, numericInput('sid', 'Reproducible Seed', 1))
                        ),
                        sliderInput("numloc", "Expected Number of Spots:",min = 100, max = 10000, value = 1000, step = 100)
                    ),

                    box(width=NULL,solidHeader=TRUE,status="warning",title="Group Assignment",
                        fluidRow(column(6,textInput('NewGroup', label = 'New Group ID')),
                        column(6,numericInput('fc', 'Fold Change in Mean', 1))),
                        # fluidRow(
                        #     column(width=12,
                        #         column(width=6,actionButton('Change',label=HTML('Confirm <br/> Assignment'),icon=icon("check-circle"))),
                        #         column(width=6,actionButton("loc_pop", label="View Location File", icon=icon("table")))
                        #     )


                            fluidRow(
                                column(width=8,offset=3,actionButton('Change',label='Confirm Assignment',icon=icon("check-circle"),width='92%',style = "margin-bottom: 2px;")
                            ),
                            fluidRow(
                                column(width=8,offset=3,actionButton("loc_pop", label="View Location File", icon=icon("table"),width='90%'))
                            )
                        )
                    ),

                    box(solidHeader=TRUE,status="warning",title="Figure Adjustment",width = NULL,
                        sliderInput("ptsize","Point Size",min=0.5,max=10,value=1),
                        sliderInput("textsize","Text Size",min=5,max=20,value=10)
                    )
                )
            ), # end of the fluidRow 

            bsModal("modalExample", "Location Data Table", "loc_pop", size = "large",
                dataTableOutput("distTable")
                # downloadButton("downloadPopLoc", "Download Location File")
            ),
            ## for now, nothing would come 
            verbatimTextOutput('brush4')
        ), ## end of tabItem dashboard

        tabItem("simCounts",
            fluidRow(
                width=12,
                column(width=3,
                    box(solidHeader=TRUE,status="warning",width=NULL,
                        title="Simulation Setting",
                        numericInput("numHighSig", "Number of Higher Signal Genes", 50, min = 0),
                        numericInput("numLowSig", "Number of Lower Signal Genes", 50, min = 0),
                        numericInput("numNoise", "Number of Noise Genes", 50, min = 0),
                        numericInput("count_sid", "Reproducible Seed", 1),

                        fluidRow(
                            column(width=4,numericInput("zeroProp", "Zero%", 0,min=0,max=1)),  
                            column(width=4,numericInput("disper_para", "Dispersion", 0,min=0)),
                            column(width=4,numericInput("meanCount_para", "Mean", 0,min=0))
                        ),

                        fluidRow(
                            column(width=8,offset=2,actionButton("countGenerate",'Generate New Data',icon=icon("spinner")),style = "margin-bottom: 2px;")
                        ),
                        fluidRow(
                            column(width=8,offset=2,actionButton("count_pop", "View CountData File", icon=icon("table")))
                        )
                    ),
                    verbatimTextOutput('countbrush'),
                    box(solidHeader=TRUE,status="warning",width=NULL,
                        title="Display Option",
                        checkboxInput("dosignal", "Show Expression Pattern For Signal Genes", value = T),
                        checkboxInput("donoise", "Show Expression Pattern For Noise Genes", value = T),

                        fluidRow(
                            column(width=6,
                                numericInput("sigidx", "Signal Gene Index", 1, min = 1, max = 10000,step=1)
                            ),
                            column(width=6,
                                numericInput("noidx", "Noise Gene Index", 1, min = 1, max = 10000,step=1)
                            )
                        )
                    ),

                    box(solidHeader=TRUE,status="warning",width=NULL,
                        title="Figure Adjustment",
                        sliderInput("ptsizeCount","Point Size",min=0.5,max=10,value=3),
                        sliderInput("textsizeCount","Text Size",min=10,max=30,value=15)
                    )
                ),

                column(width=9,
                    fluidRow(
                        column(width=6,
                            box(solidHeader=TRUE,status="success",width=NULL,
                                title="Spots Spatial Pattern",
                                plotOutput('SpotPlot'))
                        ),
                        column(width=6,
                            box(solidHeader=TRUE,status="primary",width=NULL,
                                title="Interactive Simulation Info",
                                verbatimTextOutput('countInfobrush')),
                            box(solidHeader=TRUE,status="primary",width=NULL,
                                title="Group Parameter Summary",
                                tableOutput('summary_count_table'))
                        )
                    ),
                ),

                hr(),
                column(width=9,
                    box(solidHeader=TRUE,status="success",width=NULL,
                        title="Spatial Expression Pattern",        
                        plotOutput('ExpressionPlot')
                    )
                ),

                bsModal("modalCount", "Count Data Table", "count_pop", 
                    size = "large",
                    dataTableOutput("countTable")
                    # downloadButton("downloadPopCount", "Download CountData File")
                )
            )
        )## end of tabItem simCounts
  )## end of the tabItems
)## end of dashboardBody

