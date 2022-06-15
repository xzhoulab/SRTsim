#' Generate Data with Estimated Parameters 
#' @param simsrt A object with estimated parameters from fitting step
#' @param breaktie A character string specifying how ties are treated.  Same as the "tie.method" in rank function
#' @param total_count_new The (expected) total number of reads or UMIs in the simulated count matrix.
#' @param total_count_old The total number of reads or UMIs in the original count matrix.
#' @param rrr The ratio applies to the gene-specific mean estimate, used for the fixing average sequencing depth simulation. Default is null. Its specification will override the specification of total_count_new and total_count_old. 
#' @param nn_num A integer of nearest neighbors, default is 5.
#' @param nn_func A character string specifying how the psedo-count to be generated. options include 'mean','median' and 'ransam'.
#' @param numCores The number of cores to use 
#' @param verbose Whether to show running information for srtsim_count
#' @return Returns a SRTsim object with a newly generated count matrix
#' @importFrom pdist pdist
#' @importFrom parallel mclapply
#' @export
#'
#' @examples
#'
#' ## Create a simSRT object
#' toySRT  <- createSRT(count_in=toyData$toyCount,loc_in = toyData$toyInfo)
#' set.seed(1)
#' ## Estimate model parameters for data generation
#' toySRT <- srtsim_fit(toySRT,sim_schem="tissue")
#' ## Generate synthetic data with estimated parameters
#' toySRT <- srtsim_count(toySRT)
#'
#' ## Explore the synthetic count matrix
#' simCounts(toySRT)[1:3,1:3]

srtsim_count <- function(simsrt,
                        breaktie = "random",
                        total_count_new=NULL,
                        total_count_old=NULL,
                        rrr = NULL,
                        nn_num = 5,
                        nn_func = c('mean','median',"ransam"),
                        numCores = 1,
                        verbose =FALSE
                                ){

    nn_func         <- match.arg(nn_func)

    # srtsim_fitres   <- srt@metadata


    all_nn_mat      <- simsrt@metaParam$simLocParam$simref_nmat
    # loclist         <- list(srtsim_fitres$newlocfile,srtsim_fitres$orglocfile)
    if(is.null(simcolData(simsrt))) simsrt@simcolData <- refcolData(simsrt)


    if(simsrt@metaParam$fitParam$sim_scheme=="tissue"){
        realdatals <- list(all=refCounts(simsrt))
    }else{
        sel_label <- unique(refcolData(simsrt)$label)
        countfile <- refCounts(simsrt)
        realdatals <- list()
        for(ilabel in sel_label){
            realdatals[[ilabel]]   <- countfile[,which(refcolData(simsrt)$label %in% ilabel)]
        }
        rm(countfile)
    }

    oldnum_loc      <- nrow(refcolData(simsrt))
    newnum_loc      <- nrow(simcolData(simsrt))

    param_res       <- simsrt@EstParam

    if(is.null(total_count_old)){
        total_count_old  <- sum(sapply(param_res,function(x){x$n_read}))
    }

    if(is.null(total_count_new)){
        total_count_new <- sum(sapply(param_res,function(x){x$n_read}))
    }

    if(is.null(rrr)){
        r <- (total_count_new / newnum_loc) / (total_count_old / oldnum_loc)
    }else{
        r <- rrr
    }
    
    if(verbose){cat(paste0("The ratio between the seqdepth per location is: ",r,"\n"))}

    ## consider the case without new locations
    if(oldnum_loc == newnum_loc){
        if(length(param_res)==1){
            rawsimlist <- list(srtsim_count_single(model_params=param_res[[1]],
                                                    realdata=realdatals[[1]],
                                                    rr = r, breaktie=breaktie, 
                                                    numCores=numCores))
        }else{
            rawsimlist  <- sapply(1:length(param_res),function(pp){
                    srtsim_count_single(model_params=param_res[[pp]],
                                        realdata=realdatals[[pp]],
                                        rr = r,
                                        breaktie=breaktie,
                                        numCores=numCores)
                })

        }
        rawcount    <- do.call(cbind,rawsimlist)      
        outcount    <- rawcount[,rownames(simcolData(simsrt))]

    }else{
        if(length(param_res)==1){
            # if(ncol(realdatals[[1]])!=nrow(all_nn_mat)){
            #     stop("The number of locations in the count files are different from the location files !")
            # }
            outcount <- srtsim_simulate_count_rank_diff(model_params=param_res[[1]],
                                                    realdata = realdatals[[1]],
                                                    rr = r,
                                                    nn_num = nn_num,
                                                    nn_mat = all_nn_mat,
                                                    nn_func = nn_func ,
                                                    numCores = numCores,
                                                    breaktie = breaktie
                                                )
        }else{
            ## here, nn_mat should be recalculated
            rawsimlist <- list()
            sel_label = names(param_res)

            if(length(setdiff(sel_label,unique(simcolData(simsrt)$label)))!=0){
                stop(paste0(setdiff(sel_label,unique(simcolData(simsrt)$label))," is not in the location file label!"))
            }

            loclist <- list(simcolData(simsrt),refcolData(simsrt))
            for(ilabel in sel_label){
                if(verbose){cat("Simulate count for ",ilabel,"...\n")}
                sub_new_loc <- loclist[[1]][loclist[[1]]$label == ilabel,]
                sub_old_loc <- loclist[[2]][loclist[[2]]$label == ilabel,]

                reference_real <- realdatals[[ilabel]][,rownames(sub_old_loc)]

                ## calculate pairwise distance 
                new_old_dmat    <- t(apply(sub_new_loc[,c("x","y")],1,function(x){as.matrix(pdist(sub_old_loc[,c("x","y")],x))}))

                ## save nn_k50 matrix for further label assignment and data simulation.
                nn_k50      <- t(apply(new_old_dmat,1,function(x){order(x)[1:50]}))
                rawsimlist[[ilabel]]    <- srtsim_simulate_count_rank_diff(
                                                model_params=param_res[[ilabel]],
                                                realdata = realdatals[[ilabel]],
                                                rr = r,
                                                nn_num = nn_num,
                                                nn_mat = nn_k50,
                                                nn_func = nn_func ,
                                                numCores = numCores,
                                                breaktie = breaktie
                                            )
            }
            rawcount    <- do.call(cbind,rawsimlist)      
            outcount    <- rawcount[,rownames(simcolData(simsrt))]   
        }
    }

    simsrt@simCounts <- as(outcount,"sparseMatrix")
    return(simsrt)
}     




#' Generate Data with Estimated Parameters with Same Sample Size
#' @param model_params Estimated parameters for selected models
#' @param realdata a count matrix as reference for ranking
#' @param rr The ratio applies to the gene-specific mean estimate, used for the fixing average sequencing depth simulation. Default is 1. 
#' @param breaktie A character string specifying how ties are treated.  Same as the "tie.method" in rank function
#' @param numCores The number of cores to use 
#' @return Returns a count matrix
#' @importFrom stats rnbinom rbinom
#' @noRd
#' @keywords internal

srtsim_count_single <- function(model_params,
                                realdata = NULL,
                                rr = 1,
                                breaktie = "random",
                                numCores=1
                               ){

    p1 <- length(model_params$gene_sel1)
    p2 <- length(model_params$gene_sel2)

    gene_names <- rep(NA,p1+p2)
    gene_names[model_params$gene_sel1] <- names(model_params$gene_sel1)
    gene_names[model_params$gene_sel2] <- names(model_params$gene_sel2)

    ## for the rank, need to make sure the realdata supplied is same as the number of genes being fitted
    if(p2 >0){
        realdata <- realdata[-model_params$gene_sel2,]
    }

    # print("error1")
    ## to make the algorithm faster
    ## maybe a problem when the data is too large
    realdata    <- as.matrix(realdata)
    numloc      <- ncol(realdata)

    result  <- matrix(0, nrow = p1 + p2, ncol = numloc)
    if(p1 > 0){
        if(numCores==1){
            simMat  <- t(sapply(1:p1, function(iter){
              xrank     <- rank(realdata[iter,],ties.method=breaktie)
              param     <- model_params$marginal_param1[iter, ]
              simRawExpr <- rbinom(numloc, 1, 1-param[1,1]) * rnbinom(numloc, size = param[1,2], mu = rr*param[1,3])
              return(sort(simRawExpr)[xrank])
            }
            ))
        }else{
            result_parallel <- mclapply(1:p1, function(iter){
              xrank     <- rank(realdata[iter,],ties.method=breaktie)
              param     <- model_params$marginal_param1[iter, ]
              simRawExpr <- rbinom(numloc, 1, 1-param[1,1]) * rnbinom(numloc, size = param[1,2], mu = rr*param[1,3])
              return(sort(simRawExpr)[xrank])
            },mc.cores=numCores,mc.set.seed = FALSE)
            simMat <- do.call(rbind,result_parallel)
        }
    }

    # print("error2")

    all_zero_indx <- which(apply(simMat,1,sum)==0)

    if(length(all_zero_indx)>1){
        simMat[all_zero_indx,] <- t(apply(simMat[all_zero_indx,],1,function(x){x[sample(1:length(x),1)]<-1;return(x)}))
        result[model_params$gene_sel1, ] <- simMat
    }else if(length(all_zero_indx)==1){
        tmpgene <- simMat[all_zero_indx,] 
        tmpgene[sample(1:length(tmpgene),1)]<-1
        simMat[all_zero_indx,]  <- tmpgene
        result[model_params$gene_sel1, ] <- simMat
    }else{
        result[model_params$gene_sel1, ] <- simMat
    }

    ## need to retain the colnames in the realdata, necessary for data combination
    colnames(result) <- colnames(realdata)
    rownames(result) <- gene_names
    return(result)
}






#' Generate Data with Estimated Parameters 
#' @param model_params A object with estimated parameters from fitting step
#' @param realdata a count matrix used to construct pseudo-count reference matrix for ranking
#' @param rr The ratio applies to the gene-specific mean estimate, used for the fixing average sequencing depth simulation. Default is 1. 
#' @param nn_num A integer of nearest neighbors, default is 5.
#' @param nn_mat A matrix containing neighborhood information for each new constructed location. 
#' @param nn_func A character string specifying how the psedo-count to be generated. options include 'mean','median' and 'ransam'.
#' @param breaktie A character string specifying how ties are treated.  Same as the "tie.method" in rank function
#' @param numCores The number of cores to use 
#' @return Returns a count matrix
#' 
#' @noRd
#' @keywords internal


srtsim_simulate_count_rank_diff <- function(model_params, 
                              realdata = NULL,
                              rr = 1,
                              nn_num = 5,
                              nn_mat = NULL,
                              nn_func = c('mean','median',"ransam"),
                              breaktie = "random",
                              numCores = 1){

    ##============================
    ## need to rethink the nn_mat
    ## previous was calculated based on the all locations
    ## here, we may want to calculate one based on the same type

    nn_func   <- match.arg(nn_func)
    p1        <- length(model_params$gene_sel1)
    p2        <- length(model_params$gene_sel2)

    numloc = nrow(nn_mat)

    gene_names <- rep(NA,p1+p2)
    gene_names[model_params$gene_sel1] <- names(model_params$gene_sel1)
    gene_names[model_params$gene_sel2] <- names(model_params$gene_sel2)

    if(p2 > 0){
        realdata <- realdata[-model_params$gene_sel2,]
    }


    if(nn_num==1){
        data_for_rank <- realdata[,nn_mat[,1]]
    }else{
        if(nn_func=="mean"){
            data_for_rank <- apply(nn_mat[,1:nn_num],1,function(x){rowMeans(realdata[,x])})
        }else if(nn_func=="median"){
            data_for_rank <- apply(nn_mat[,1:nn_num],1,function(x){matrixStats::rowMedians(realdata[,x])})
        }else if(nn_func=="ransam"){
            # set.seed(1)
            ransam_idx    <- apply(nn_mat[,1:nn_num],1,function(xx){sample(xx,1)})
            data_for_rank <- realdata[,ransam_idx]
        }else{
            stop("The nearest neighbor impute function has to be mean, median or ransam")
        }
    }

    rownames(data_for_rank) <- names(model_params$gene_sel1)
    colnames(data_for_rank) <- rownames(nn_mat)

    ## initialize all zeros, and for those nonzero vectors, generate from the ZINB 
    result <- matrix(0, nrow = p1 + p2, ncol = numloc)

    if(p1 > 0){
        if(numCores==1){
            simMat  <- t(sapply(1:p1, function(iter){
                xrank     <- rank(data_for_rank[iter,],ties.method=breaktie)
                param     <- model_params$marginal_param1[iter, ]
                # simRawExpr <- rbinom(numloc, 1, 1-param[1]) * rnbinom(numloc, size = param[2], mu = rr*param[3])
                simRawExpr <- rbinom(numloc, 1, 1-param[1,1]) * rnbinom(numloc, size = param[1,2], mu = rr*param[1,3])
                return(sort(simRawExpr)[xrank])
            }))
        }else{
            result_parallel <- mclapply(1:p1, function(iter){
              xrank     <- rank(data_for_rank[iter,],ties.method=breaktie)
              param     <- model_params$marginal_param1[iter, ]
              # simRawExpr <- rbinom(numloc, 1, 1-param[1]) * rnbinom(numloc, size = param[2], mu = rr*param[3])
              simRawExpr <- rbinom(numloc, 1, 1-param[1,1]) * rnbinom(numloc, size = param[1,2], mu = rr*param[1,3])
              return(sort(simRawExpr)[xrank])
            },mc.cores=numCores,mc.set.seed = FALSE)
            simMat <- do.call(rbind,result_parallel)
        }
    }

    all_zero_indx <- which(apply(simMat,1,sum)==0)

    if(length(all_zero_indx)>1){
        simMat[all_zero_indx,] <- t(apply(simMat[all_zero_indx,],1,function(x){x[sample(1:length(x),max(round(numloc/ncol(realdata)),1))]<-1;return(x)}))
        result[model_params$gene_sel1, ] <- simMat
    }else if(length(all_zero_indx)==1){
        tmpgene <- simMat[all_zero_indx,] 
        tmpgene[sample(1:length(tmpgene),1)]<-1
        simMat[all_zero_indx,]  <- tmpgene
        result[model_params$gene_sel1, ] <- simMat
    }else{
        result[model_params$gene_sel1, ] <- simMat
    }

    colnames(result) <- rownames(nn_mat)
    rownames(result) <- gene_names
    return(result)
}

