#' Visualize summarized metrics for reference data and synthetic data
#' @param simsrt A SRTsim object
#' @param metric_type Specification of metrics to be plotted. Default value is 'all', which will plot all six metrics: including four gene-wise metrics and two location-wise metrics. "genewise" will produce violin plots for all four gene-wise metrics; "locwise" will produce violin plots for all two location-wise metrics; "GeneMean", "GeneVar", "GeneCV", "GeneZeroProp", "LocZeroProp", and "LocLibSize" will produce single violin plot for the corresponding metric. 
#' @param colorpalette Specification of color palette to be passed to \code{palette} in the \code{scale_fill_brewer}. Default is "Set3"
#' @param axistextsize Specification of axis font size. Default is 12. 
#' @return Returns a list of ggplots
#' @importFrom ggpubr ggarrange
#' @export
visualize_metrics <- function(simsrt,metric_type=c("all","genewise","locwise","GeneMean","GeneVar","GeneCV","GeneZeroProp","LocZeroProp","LocLibSize"),colorpalette="Set3",axistextsize=12){

	if(is.null(simsrt@metaParam$compared)){
		stop("Please run compareSRT prior to visualize the metrics!")
	}

	metric_type <- match.arg(metric_type)

	# plotlist 	<- list()
	allmetriclist <- c("GeneMean","GeneVar","GeneCV","GeneZeroProp","LocZeroProp","LocLibSize")
	if(metric_type=="all"){
		plotlist <- base::lapply(allmetriclist,function(x){visualize_metric_single(simsrt,selmetric=x)})
	}else if(metric_type=="genewise"){
		plotlist <- base::lapply(allmetriclist[1:4],function(x){visualize_metric_single(simsrt,selmetric=x)})
	}else if(metric_type=="locwise"){
		plotlist <- base::lapply(allmetriclist[5:6],function(x){visualize_metric_single(simsrt,selmetric=x)})
	}else{
		plotlist <- base::lapply(metric_type,function(x){visualize_metric_single(simsrt,selmetric=x)})
	}

	return(ggpubr::ggarrange(plotlist=plotlist))
}


#' Extracted summarized metrics for reference data and synthetic data
#' @param simsrt A SRTsim object
#' @param metric Specification of metrics to be plotted. 
#' @return Returns a dataframe for ggplot 

get_metrics_pd <- function(simsrt,metric="GeneMean"){
	if(metric %in% c("GeneMean","GeneVar","GeneCV","GeneZeroProp")){
		refpd <- cbind.data.frame(value=simsrt@refrowData[,metric],GeneName=rownames(simsrt@refrowData),group=simsrt@refrowData$sampleLabel)
		simpd <- cbind.data.frame(value=simsrt@simrowData[,metric],GeneName=rownames(simsrt@simrowData),group=simsrt@simrowData$sampleLabel)
		pd 	<- rbind(refpd,simpd)
	}else{
		refpd <- cbind.data.frame(value=simsrt@refcolData[,metric],LocID=rownames(simsrt@refcolData),group=simsrt@refcolData$sampleLabel)
		simpd <- cbind.data.frame(value=simsrt@simcolData[,metric],LocID=rownames(simsrt@simcolData),group=simsrt@simcolData$sampleLabel)
		pd 	<- rbind(refpd,simpd)
	}
	return(pd)
}


#' Extracted summarized metrics for reference data and synthetic data
#' @param simsrt A SRTsim object
#' @param selmetric Specification of metrics to be plotted. 
#' @param colorpalette Specification of color palette to be passed to \code{palette} in the \code{scale_fill_brewer}. Default is "Set3"
#' @param axistextsize Specification of axis font size. Default is 12. 
#' @return Returns a ggplot object

#' @noRd
#' @keywords internal


visualize_metric_single <- function(simsrt,selmetric="GeneMean",colorpalette="Set3",axistextsize=12){
	gpd <- get_metrics_pd(simsrt,metric = selmetric)

	metriclist <-  c("Gene Mean", "Gene Variance", "Gene CV",
                    "Gene Zero Prop.", "Loc Zero Prop.", "Loc Library Size")
	names(metriclist) <- c("GeneMean","GeneVar","GeneCV","GeneZeroProp","LocZeroProp","LocLibSize")

	stats_plot <- ggplot(gpd, aes(x = group, y = value,fill=group)) +
		  geom_violin(scale = 'width', trim = TRUE) +
		  scale_fill_brewer(palette = colorpalette)+
		  theme_bw()+
		  theme(panel.grid.major = element_blank(),
		        panel.grid.minor = element_blank(),
		        axis.text = element_text(size = axistextsize),
		        legend.position="none",
		        plot.title = element_text(size = axistextsize+3, face = "bold",hjust = 0.5)) +
		  ggtitle(metriclist[selmetric])+
		  xlab("") + ylab("")
	return(stats_plot)
}
