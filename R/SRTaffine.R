#' Generate Data with Estimated Parameters For A New Designed Pattern
#' @param simsrt A SRTsim object with estimated parameters from fitting step
#' @param reflabel A character vector specifying labels for reference regions
#' @param targetlabel A character vector specifying labels for target regions
#' @param breaktie A character string specifying how ties are treated.  Same as the "tie.method" in rank function
#' @param nn_func A character string specifying how the psedo-count to be generated. options include 'mean','median' and 'ransam'.
#' @param nn_num A integer of nearest neighbors, default is 5.
#' @param local_sid  A numberic seed used locally for the affine transformation. Default is NULL. 
#' @param numCores A number of cores to use 
#' @return Returns a SRTsim object with a newly generated count matrix
#' 
#' @export
#' @importFrom matrixStats rowMedians
#' @importFrom methods as
#' @importFrom parallel mclapply
#' @examples
#' 
#' ## Prepare Data From LIBD Sample
#' subinfo <- exampleLIBD$info[,c("imagecol","imagerow","layer")]
#' colnames(subinfo) <- c("x","y","label")
#' gns 	<- c("ENSG00000168314","ENSG00000183036", "ENSG00000132639" )
#' 
#' ## Create a simSRT Object with Three Genes For a Fast Example
#' simSRT1 <- createSRT(count_in= exampleLIBD$count[gns,],loc_in =subinfo)
#' 
#' ## Estimate model parameters for data generation: domain-specific 
#' simSRT1 <- srtsim_fit(simSRT1,sim_schem="domain")
#' 
#' ## Define New Layer Structures
#' simSRT1@refcolData$target_label <- "NL1"
#' simSRT1@refcolData$target_label[simSRT1@refcolData$label %in% paste0("Layer",4:5)] <- "NL2"
#' simSRT1@refcolData$target_label[simSRT1@refcolData$label %in% c("Layer6","WM")] <- "NL3"
#' 
#' ## Perform Data Generation for New Defined Layer Structures
#' ## Reference: WM --> NL3, Layer5--> NL2, Layer3 --> NL1
#' simSRT1 <- srtsim_count_affine(simSRT1,
#' 								reflabel=c("Layer3","Layer5","WM"),
#' 								targetlabel=c("NL1","NL2","NL3"),
#' 								nn_func="ransam"
#' 								)
#' 
#' ## Visualize the Expression Pattern for Gene of Interest
#' visualize_gene(simsrt=simSRT1,plotgn = "ENSG00000168314",rev_y=TRUE,ptsizeCount=1)


srtsim_count_affine <- function(simsrt,
								reflabel, 
								targetlabel,  
								breaktie="random",
								nn_func = c('mean','median',"ransam"),
								nn_num = 5,
								local_sid=NULL,
								numCores=1){

	nn_func     	<- match.arg(nn_func)

	if(length(reflabel)!=length(targetlabel)){
		stop("reflabel should have same length as targetlabel!")
	}

	newrawinfo <- refcolData(simsrt)[,c("x","y","target_label")]
	colnames(newrawinfo) <- c("x","y","label")
	newinfo 	<- newrawinfo[newrawinfo$label %in% targetlabel,]

	all_affine_data <- list()
	for(inl in 1:length(targetlabel)){
		all_affine_data[[inl]]	<- srtsim_count_affine_single(simsrt,
											orig_label=reflabel[inl],
											new_label=targetlabel[inl],
											breaktie=breaktie,
											nn_func = nn_func,
											nn_num = nn_num,
											local_sid=local_sid,
											numCores=numCores
										)
	}

	all_affine_df  			<- do.call(cbind,all_affine_data)
	rm(all_affine_data)
	final_df 				<- all_affine_df[,rownames(newinfo)]
	# simcolData(simsrt) 		<- newinfo
	# simCounts(simsrt) 		<- as(final_df,"sparseMatrix")

	simsrt@simcolData		<- newinfo
	simsrt@simCounts 		<- as(final_df,"sparseMatrix")

	return(simsrt)
}

#' Generate Count Data with Estimated Parameters For a Single Region
#' @param simsrt A SRTsim object with estimated parameters from fitting step
#' @param orig_label A character string specifying label for reference region
#' @param new_label A character string specifying label for target region
#' @param breaktie A character string specifying how ties are treated.  Same as the "tie.method" in rank function
#' @param nn_func A character string specifying how the psedo-count to be generated. options include 'mean','median' and 'ransam'.
#' @param nn_num A integer of nearest neighbors, default is 5.
#' @param local_sid  A numberic seed used locally for the affine transformation. Default is NULL. 
#' @param numCores A number of cores to use 
#' @return Returns a newly generated count matrix
#' @noRd
#' @keywords internal
#' @importFrom Matrix rowMeans
#' @importFrom stats rnbinom rbinom
#' @importFrom pdist pdist

srtsim_count_affine_single <- function(simsrt,
								orig_label, 
								new_label,  
								breaktie="random",
								nn_func = c('mean','median',"ransam"),
								nn_num = 5,
								local_sid=NULL,
								numCores=1){

    nn_func     	<- match.arg(nn_func)
    # param_res   	<- srtsim_fitres$fit_result

	param_res       <- simsrt@EstParam
    model_params 	<- param_res[[orig_label]]

    ## set the imagerow and imagecol outside, here, we only need the relative locations
    ## for the visualization purpose, we can always replace the location matrix
	original_info 	<- simsrt@refcolData

    ## check the location index in two files
    if(ncol(original_info)!=4){
        stop("location file should be organized in four columns: x,y,label, and target_label!")
    }

    ## in case the columns are not named in this way
    colnames(original_info) <- c("x","y","label","target_label")

    ## find the transmat, and create the projection for reference
	old_slt_idx   	<- which(original_info$label%in%orig_label)
	new_slt_idx   	<- which(original_info$target_label%in%new_label)

	ref_trans_xy 	<- transform_mat_func(info_df= as.data.frame(original_info),
											idx1 = old_slt_idx,
											idx2 = new_slt_idx,
											local_sid=local_sid)

	new_xy        	<- original_info[new_slt_idx,c("x","y")]

    ## typically, if label is well-defined outside, then the unorder should be same as the refdata
	## well defined, means only one orig_label
	refdata 		<- refCounts(simsrt)[,rownames(ref_trans_xy)]

	p1 <- length(model_params$gene_sel1)
	p2 <- length(model_params$gene_sel2)

	## for the rank, need to make sure the realdata supplied is same as the number of genes being fitted
	if(p2 >0){
		clean_refdata <- refdata[-model_params$gene_sel2,]
	}else{
		clean_refdata <- refdata
	}

	## calculate pairwise distance 
	dd_matrix   	<- t(apply(new_xy,1,function(x){as.matrix(pdist(ref_trans_xy,x))}))
	nn_mat          <- t(apply(dd_matrix,1,function(x){order(x)[1:min(50,ncol(dd_matrix))]}))

	## find the nearest neighbor
	if(nn_num==1){
		data_for_rank <- clean_refdata[,nn_mat[,1]]
	}else{
		if(nn_func=="mean"){
            data_for_rank <- apply(nn_mat[,1:nn_num],1,function(x){rowMeans(clean_refdata[,x])})
        }else if(nn_func=="median"){
            data_for_rank <- apply(nn_mat[,1:nn_num],1,function(x){matrixStats::rowMedians(clean_refdata[,x])})
        }else if(nn_func=="ransam"){
            # set.seed(1)
            ransam_idx    <- apply(nn_mat[,1:nn_num],1,function(xx){sample(xx,1)})
            data_for_rank <- clean_refdata[,ransam_idx]
        }else{
            stop("The nearest neighbor impute function has to be mean, median or ransam")
        }
	}   

	
	newNumLoc = length(new_slt_idx)
	## initialize all zeros, and for those nozero values, generate from the negative binomial 
	result <- matrix(0, nrow = p1 + p2, ncol = newNumLoc)

	result_parallel <- mclapply(1:p1, function(iter){
	        xrank     <- rank(data_for_rank[iter,],ties.method=breaktie)
	        param     <- model_params$marginal_param1[iter, ]
	        simRawExpr <- rbinom(newNumLoc, 1, 1-param[1,1]) * rnbinom(newNumLoc, size = param[1,2], mu = param[1,3])
	        return(sort(simRawExpr)[xrank])
	      },
	      mc.set.seed = FALSE,
	      mc.cores=numCores)
	result31 <- do.call(rbind,result_parallel)

	all_zero_indx <- which(apply(result31,1,sum)==0)

	if(length(all_zero_indx)>1){
		result31[all_zero_indx,] <- t(apply(result31[all_zero_indx,],1,function(x){x[sample(1:length(x),max(round(newNumLoc/ncol(refdata)),1))]<-1;return(x)}))
		result[model_params$gene_sel1, ] <- result31
	}else if(length(all_zero_indx)==1){
		tmpgene <- result31[all_zero_indx,] 
		tmpgene[sample(1:length(tmpgene),1)]<-1
		result31[all_zero_indx,]  <- tmpgene
		result[model_params$gene_sel1, ] <- result31
	}else{
		result[model_params$gene_sel1, ] <- result31
	}

	## location order should be consistent with that in the info
	colnames(result) <- rownames(new_xy)
	rownames(result) <- rownames(refdata)
	return(result)
}



#' Generate Data with Estimated Parameters For A New Designed Pattern
#' @param info_df A dataframe with four columns: x, y, orig_label and new_label
#' @param idx1 A numeric vector specifying the index for the orig_label
#' @param idx2 A numeric vector specifying the index for the new_label
#' @param local_sid  A numberic seed used locally for the affine transformation. Default is NULL. 
#' @param trans_type A character string specifying the transformation type used in the \code{computeTransform}
#' @return Returns a dataframe with transformed location information
#' @noRd
#' @keywords internal
#' @importFrom sf st_as_sf
#' @importFrom concaveman concaveman
#' @importFrom Morpho computeTransform applyTransform
#' @importFrom magrittr %>%


transform_mat_func <- function(info_df,idx1,idx2,local_sid=NULL,trans_type="affine"){

	loc1 <- info_df[idx1,]
	loc2 <- info_df[idx2,]

	pnts1    <- loc1[,c("x", "y")] %>% st_as_sf(coords = c("x", "y"))
	polygon1 <- concaveman(pnts1,2.0) 
	poly_coords1 = as.data.frame(as.matrix(polygon1$polygons[[1]]))

	pnts2    <- loc2[,c("x", "y")] %>% st_as_sf(coords = c("x", "y"))
	polygon2 <- concaveman(pnts2,2.0) 
	poly_coords2 = as.data.frame(as.matrix(polygon2$polygons[[1]]))

	n1 <- nrow(poly_coords1)
	n2 <- nrow(poly_coords2)
	k  <- min(n1,n2)

	## set seed locally
	## https://stackoverflow.com/questions/14324096/setting-seed-locally-not-globally-in-r
	if(!is.null(local_sid)){
		old <- .Random.seed
		set.seed(local_sid)
		moving_data <- as.matrix(poly_coords1[sort(sample(n1,k)),])
		fixed_data  <- as.matrix(poly_coords2[sort(sample(n2,k)),])
		.Random.seed <- old
	}else{
		moving_data <- as.matrix(poly_coords1[sort(sample(n1,k)),])
		fixed_data  <- as.matrix(poly_coords2[sort(sample(n2,k)),])
	}

	trafo <- computeTransform(
	  x=fixed_data,
	  y=moving_data,
	  type = trans_type,
	  reflection = FALSE,
	  lambda = 1e-08,
	  weights = NULL,
	  centerweight = FALSE,
	  threads = 1
	)

	trans_loc_mat <- applyTransform(as.matrix(loc1[,1:2]),trafo)
	trans_loc_df  <- as.data.frame(trans_loc_mat[,1:2])
	colnames(trans_loc_df)  <- c('newx','newy')
	return(trans_loc_df)
}