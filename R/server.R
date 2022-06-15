#' @importFrom plotly renderPlotly ggplotly layout event_data 
#' @importFrom ggpubr ggarrange
#' @importFrom viridis scale_color_viridis
#' @importFrom utils read.csv write.csv
#' @importFrom stats aggregate
#' @importFrom magrittr %>%
#' @import ggplot2
server <- function(input, output, session){
    
    ## global reactiveValues list
    globalVarList <- reactiveValues()

    filedata <- reactive({
        inFile <- input$datafile
        if(is.null(inFile)) return(NULL)
        read.csv(
            inFile$datapath, 
            header = input$header,
            row.names=as.numeric(input$rownames)
        )
    })

    newLocdata <- reactive({
        get_new_coordinates2(
          ns=input$numloc, 
          bg_shape=input$shape,
          isid=input$sid,
          lay_out=input$arrange,
          preLoc = filedata()
        )
    })

    observeEvent(input$reload, {
        session$reload()
    })


    output$plot1 <- plotly::renderPlotly({
        if(is.null(globalVarList$vv)){
            loclist <- newLocdata()
            g1      <- loclist$simLoc
            globalVarList$locdf <- g1
            globalVarList$simLocParam <- loclist$simLocParam
        }else{
            g1  <- globalVarList$locdf 
        }

        gg_scatter <- ggplot(data=g1,aes(x=x, y=y,color=group)) + 
                        geom_point(size=input$ptsize) + 
                        xlab("Coordinate 1") + 
                        ylab("Coordinate 2")+ 
                        theme(legend.title=element_blank(),
                                text = element_text(size=input$textsize))+ 
                        coord_fixed()

        plotly::ggplotly(gg_scatter, height = 500, width= 1.1*500*(diff(range(g1$x))/diff(range(g1$y)))) %>% plotly::layout(dragmode = 'select',legend = list(title=list(text='<b> Group </b>',side="top"),orientation = "v", x=1.05))
    })

    ## initial table
    output$summary_table <- renderTable({
            g1 <- globalVarList$locdf
            summary_df <- cbind.data.frame(aggregate(g1$foldchange~g1$group,FUN=length),
                                    aggregate(g1$foldchange~g1$group,FUN=mean)[,2])
            colnames(summary_df) <- c("Group","NumSpots","Fold Change")
            summary_df$PropSpots <- round(summary_df$NumSpots/sum(summary_df$NumSpots),3)
            return(summary_df[,c("Group","NumSpots","PropSpots","Fold Change")])
        },
        digits=3
    )

    output$brush <- renderText({
        if(is.null(filedata())&(input$shape %in% c("User Define Shape","User Define Spots"))){
            "!!! Warning: No dataframe provided, simulate using the square backgound shape \n"
        }
        # g1 <- newLocdata()$simLoc
        g1 <- globalVarList$locdf
        group <- g1[,3]
        names(group) <- rownames(g1)
        d  <- plotly::event_data('plotly_selected')
        if (is.null(d)){
            "Select points (i.e., box/lasso) to define new group"
        }else{
            dd <- round(cbind(d[[3]],d[[4]]),3)
            vv <- group[which(round(g1[,1],3) %in% dd[,1] & round(g1[,2],3) %in% dd[,2])]
            globalVarList$vv <- vv
            paste0("Number of points selected: ",length(vv)," (",100*round(length(vv)/length(group),3),"%)")
        } 

    })

    output$distTable <- renderDataTable({
            locdf <- globalVarList$locdf
            display_locdf <- cbind.data.frame(LocID=rownames(locdf),locdf)
            display_locdf
        }, 
        options = list(pageLength=10)
    )

    output$SpotPlot <- renderPlot({
        locdf <- globalVarList$locdf
        gg_scatter_spot <- ggplot(data=locdf,aes(x=x, y=y,color=group)) + 
                            geom_point(size=input$ptsizeCount) + 
                            xlab("Coordinate 1") + 
                            ylab("Coordinate 2")+ 
                            theme(legend.title=element_blank(),
                                text = element_text(size=input$textsizeCount)
                            )+ coord_fixed()
        ptlist <- list(gg_scatter_spot)
        ggpubr::ggarrange(plotlist=ptlist,ncol=length(ptlist))
    })

    observeEvent(input$Change > 0, {
        if (!is.null(globalVarList$vv)){
            ## update the location file in the global list
            vv      <- globalVarList$vv
            locdf   <- globalVarList$locdf
            locdf[which(row.names(locdf) %in% names(vv)),]$group <- input$NewGroup
            locdf[which(row.names(locdf) %in% names(vv)),]$foldchange <- input$fc
            globalVarList$locdf <- locdf

            ## check if there is multiple effects in the same group
            output$brush3 <- renderText({
                 locdf <- globalVarList$locdf
                check_multi <- aggregate(locdf$foldchange~locdf$group,FUN=unique)
                if(length(grep(",",check_multi[,2]))!=0){
                    paste0("Group ", check_multi[,1][grep(",",check_multi[,2])]," has more than one effect size, consider redefining the effect size for the group")
                }
            })

            output$SpotPlot <- renderPlot({
                locdf <- globalVarList$locdf
                gg_scatter_spot <- ggplot(data=locdf,aes(x=x, y=y,color=group)) + 
                                    geom_point(size=input$ptsizeCount) + 
                                    xlab("Coordinate 1") + ylab("Coordinate 2")+ 
                                    theme(legend.title=element_blank(),
                                            text = element_text(size=input$textsizeCount))+
                                    coord_fixed()
                ptlist <- list(gg_scatter_spot)
                ggpubr::ggarrange(plotlist=ptlist,ncol=length(ptlist))
            })
        } ## end of the is.null if
    })  

    ## Note global loc file is used here
    newCountdata <- reactive({
        pattern_count_func2(
            pattern_in = globalVarList$locdf,
            numHighSignal = input$numHighSig,
            numLowSignal = input$numLowSig,
            numBG = input$numNoise,
            disper_in = input$disper_para,
            zero_in = input$zeroProp,
            mu_in = input$meanCount_para,
            isid = input$count_sid
        )
    })

    ## without >0, won't run before click
    observeEvent(input$countGenerate,{
        locdf <- globalVarList$locdf
        output$countbrush <- NULL
        if(input$zeroProp < 0){
            output$countbrush <- renderText({
                paste0("The minimum zero proportion allowed is 0!")
            })
        }else if(input$zeroProp> 1){
            output$countbrush <- renderText({
                paste0("The max zero proportion allowed is 1!")
            })
        }

        if((input$numHighSig + input$numLowSig + input$numNoise)== 0){
            output$countbrush <- renderText({
                paste0("Total number of genes to simulate is zero, please change the simulation setting")
            })
        }else{
            output$countInfobrush <- renderText({
                paste0("Baseline Mu: ",input$meanCount_para,"\n",
                    "Baseline Overdispersion: ",input$disper_para,"\n",
                    "Additional Zero Proportion: ",input$zeroProp,"\n",
                    "Number of Locations: ", nrow(locdf),"\n",
                    "Number of Higher Signal Genes Simulated: ", input$numHighSig,"\n",
                    "Number of Lower Signal Genes Simulated: ", input$numLowSig,"\n",
                    "Number of Noise Genes Simulated: ", input$numNoise,"\n")
            })
            
            countlist       <- newCountdata()
            countdf         <- countlist$ctdf
            simcountParam   <- countlist$simcountParam

            globalVarList$countdf <- countdf
            globalVarList$simcountParam <- simcountParam

            ## has to be after the countdf, otherwise, the plot is not updated with seed change.
            newPlotdata <- reactive({
                combined_df     <- cbind.data.frame(locdf,apply(countdf,1,relative_func))
                selected_gene   <- c(paste0("signal",input$sigidx),paste0("noise",input$noidx))
                npdf            <- combined_df[,c("x","y","group",selected_gene)]
                names(npdf)     <- c("x","y","group","signal_gene","noise_gene")
                return(npdf)
            })

            ## count summary table
            output$summary_count_table <- renderTable({
                    loc_mu_df <- locdf
                    summary_df <- cbind.data.frame(aggregate(loc_mu_df$foldchange~loc_mu_df$group,FUN=length),
                            aggregate(loc_mu_df$foldchange~loc_mu_df$group,FUN=mean)[,2])
                    summary_df$Mu <- input$meanCount_para*summary_df[,3]
                    colnames(summary_df) <- c("Group","NumSpots","Fold Change","Mu")
                    summary_df$PropSpots <- round(summary_df$NumSpots/sum(summary_df$NumSpots),3)
                    summary_df$Dispersion    <- input$disper_para
                    return(summary_df[,c("Group","NumSpots","PropSpots","Mu","Dispersion")])
                },digits=3)

            ## update pop display of the countdata after the group assignment
            output$countTable <- renderDataTable({
                display_count_df <- cbind.data.frame(GeneID=rownames(countdf),countdf)
                display_count_df[,1:10]
                }, options = list(pageLength=10)
            )

            output$ExpressionPlot <- renderPlot({
                pltdf <- newPlotdata()
                if(input$dosignal){
                    gg_scatter_signal <- ggplot(data=pltdf,aes(x=x, y=y,color=signal_gene)) + 
                                            geom_point(size=input$ptsizeCount)+ 
                                            viridis::scale_color_viridis(option="D",direction=-1) + 
                                            xlab("Coordinate 1") + ylab("Coordinate 2")+ 
                                            ggtitle(paste0("Signal Gene #",input$sigidx)) +
                                            theme(legend.title=element_blank(),
                                                text = element_text(size=input$textsizeCount))+ coord_fixed()
                }else{
                    gg_scatter_signal <- NULL
                }

                if(input$donoise){
                    gg_scatter_noise <- ggplot(data=pltdf,aes(x=x, y=y,color=noise_gene)) + 
                                          geom_point(size=input$ptsizeCount)+ 
                                          viridis::scale_color_viridis(option="D",direction=-1) + 
                                          xlab("Coordinate 1") + ylab("Coordinate 2")+ 
                                          ggtitle(paste0("Noise Gene #",input$noidx)) + 
                                          theme(legend.title=element_blank(),text = element_text(size=input$textsizeCount)) + coord_fixed()
                }else{
                    gg_scatter_noise <- NULL
                }
                ptlist    <- list(gg_scatter_signal,gg_scatter_noise)
                to_delete <- !sapply(ptlist,is.null)
                ptlist    <- ptlist[to_delete] 
                if(length(ptlist)==0){
                    return(NULL)
                }else if(length(ptlist)==2){
                    ggpubr::ggarrange(plotlist=ptlist,ncol=length(ptlist), common.legend = TRUE, legend="right")
                }else{
                    ggpubr::ggarrange(plotlist=ptlist,ncol=length(ptlist))
                }
            }) ## end of ExpressionPlot
        }## end of gene number else
    })## end of countGenerate  

    ## Quit the app and get back to the R
    observeEvent(input$exit,{
        countdf         <- globalVarList$countdf 
        simcountParam   <- globalVarList$simcountParam 

        locdf       <- globalVarList$locdf
        simLocParam <- globalVarList$simLocParam

        session$sendCustomMessage(type = "closeWindow", message = "message")
        stopApp(list(simCount = countdf,simInfo = locdf,simcountParam = simcountParam,simLocParam=simLocParam))
    })
} ## end of server
