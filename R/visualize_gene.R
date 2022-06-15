#' Visualize expression pattern for the gene of interest in reference data and synthetic data
#' @param simsrt A SRTsim object
#' @param plotgn A gene name selected for visualization
#' @param ptsizeCount Specification of point size. Default is 2. 
#' @param textsizeCount Specification of axis font size. Default is 12.
#' @param rev_y Logical indicating whether to reverse the y axis. Default is FALSE. Useful for Visualize the LIBD data.
#' @param virOption Specification of \code{option} in the \code{scale_color_viridis}. Default is "D". User can choose a letter from 'A' to 'H'.
#' @param virDirection Specification of \code{direction} in the \code{scale_color_viridis}. Default is "-1".  User can choose '1' or '-1'.
#' @return Returns two expression plots for the gene of interest
#' @importFrom ggpubr ggarrange
#' @importFrom viridis scale_color_viridis
#' @import ggplot2
#' @export
#' @examples
#' 
#' ## Create a simSRT object
#' toySRT  <- createSRT(count_in=toyData$toyCount,loc_in = toyData$toyInfo)
#' set.seed(1)
#' ## Create New Locations within Profile
#' toySRT2 <- srtsim_newlocs(toySRT,new_loc_num=1000)
#'
#' ## Estimate model parameters for data generation
#' toySRT2 <- srtsim_fit(toySRT2,sim_schem="tissue")
#'
#' ## Generate synthetic data with estimated parameters
#' toySRT2     <- srtsim_count(toySRT2,rrr=1)
#'
#' ## compare the expression pattern of HLA-B in synthetic data and reference data
#' visualize_gene(simsrt=toySRT2,plotgn = "HLA-B")


visualize_gene <- function(simsrt,plotgn=NULL,ptsizeCount=2,textsizeCount= 12,rev_y = FALSE,virOption="D",virDirection=-1){
    ref_pd <- cbind.data.frame(simsrt@refcolData[,c("x","y")],relative_func(simsrt@refCounts[plotgn,]))
    sim_pd <- cbind.data.frame(simsrt@simcolData[,c("x","y")],relative_func(simsrt@simCounts[plotgn,]))
    colnames(ref_pd) <- colnames(sim_pd) <- c("x","y","sel_gene")

    gg_scatter_real <- ggplot(data=ref_pd,aes(x=x, y=y,color=sel_gene)) + 
            geom_point(size=ptsizeCount)+
            viridis::scale_color_viridis(name="Relative Expression",option=virOption,direction=virDirection)+ 
            theme_bw()+ 
            theme(
                legend.text=element_text(size=textsizeCount),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                axis.line = element_blank(),
                axis.text = element_blank(),
                axis.ticks = element_blank(),
                text = element_text(size=textsizeCount),
                plot.title = element_text(size = textsizeCount, hjust = 0.5))+ 
            ggtitle("Reference")+
            coord_fixed()+
            guides(colour = guide_colourbar(title.position="top", title.hjust = 0.5,barwidth = 10))


    gg_scatter_sim <- ggplot(data=sim_pd,aes(x=x, y=y,color=sel_gene)) + 
            geom_point(size=ptsizeCount)+
            viridis::scale_color_viridis(name="Relative Expression",option=virOption,direction=virDirection) + 
            theme_bw()+ 
            theme(
                legend.text=element_text(size=textsizeCount),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                axis.line = element_blank(),
                axis.text = element_blank(),
                axis.ticks = element_blank(),
                text = element_text(size=textsizeCount),
                plot.title = element_text(size = textsizeCount, hjust = 0.5))+ 
            ggtitle("SRTsim")+
            coord_fixed()+
            guides(colour = guide_colourbar(title.position="top", title.hjust = 0.5,barwidth = 10))

    if(rev_y){
        combined_plot <- ggpubr::ggarrange(gg_scatter_real+ylim(rev(range(ref_pd$y))),gg_scatter_sim+ylim(rev(range(sim_pd$y))),common.legend = TRUE,legend="bottom") 
    }else{
        combined_plot <- ggpubr::ggarrange(gg_scatter_real,gg_scatter_sim,common.legend = TRUE,legend="bottom") 
    }

    plot2show   <- ggpubr::annotate_figure(combined_plot, top = ggpubr::text_grob(plotgn, face = "bold", size = textsizeCount+2))

    return(plot2show)
}


# ## also defined in the shiny utilies_func
# relative_func <- function(rx){
#   rexpr   = (rx-min(rx))/(max(rx)-min(rx))
#   return(rexpr)
# }


