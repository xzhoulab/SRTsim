
#' Fit the marginal distributions for each row of a count matrix
#' @param simsrt    A SRTsim object
#' @param marginal   Specification of the types of marginal distribution.Default value is 'auto_choose' which chooses between ZINB, NB, ZIP, and Poisson by a likelihood ratio test (lrt),AIC and whether there is underdispersion.'zinb' will fit the ZINB model. If there is underdispersion, it will fit the Poisson model. If there is no zero at all or an error occurs, it will fit an NB model instead.'nb' fits the NB model and chooses between NB and Poisson depending on whether there is underdispersion. 'poisson' simply fits the Poisson model.'zip' fits the ZIP model and chooses between ZIP and Poisson by a likelihood ratio test
#' @param sim_scheme   a character string specifying simulation scheme. "tissue" stands for tissue-based 
#'                      simulation; "domain" stands for domain-specific simulation. Default is "tissue".
#' @param min_nonzero_num  The minimum number of non-zero values required for a gene to be fitted. Default is 2.  
#' @param maxiter      The number of iterations for the model-fitting. Default is 500. 
#' @return Returns an object with estimated parameters
#' @importClassesFrom S4Vectors SimpleList
#' @export 
#'
#' @examples
#'
#' ## Create a simSRT object
#' toySRT  <- createSRT(count_in=toyData$toyCount,loc_in = toyData$toyInfo)
#' set.seed(1)
#' ## Estimate model parameters for data generation
#' toySRT <- srtsim_fit(toySRT,sim_schem="tissue")
#'

srtsim_fit <- function(simsrt,
                        marginal = c("auto_choose", "zinb", "nb", "poisson","zip"),
                        sim_scheme = c("tissue","domain"),
                        min_nonzero_num = 2,maxiter =500){

    countfile       <- refCounts(simsrt)
    locfile         <- refcolData(simsrt)[,c('x','y','label')]

    marginal    <- match.arg(marginal)
    sim_scheme  <- match.arg(sim_scheme)

    fit_result  <- list()

    if(sim_scheme == "tissue"){
        fit_result[['tissue']]  <- fit_single(x=countfile,
                                        marginal=marginal,
                                        min_nonzero_num=min_nonzero_num,
                                        maxiter=maxiter)
    }else if(sim_scheme == "domain"){
        sel_label <- unique(locfile$label)
        for(ilabel in sel_label){
            sel_count   <- countfile[,which(locfile$label %in% ilabel)]
            fit_result[[ilabel]] <- fit_single(x=sel_count,
                                                min_nonzero_num=min_nonzero_num,
                                                marginal=marginal,
                                                maxiter=maxiter)
        }

    }else{
        stop("sim_scheme has to be tissue or domain!")
    }

    simsrt@EstParam           <- fit_result
    simsrt@metaParam$fitParam <- S4Vectors::SimpleList(marginal=marginal,
                                    sim_scheme = sim_scheme,
                                    min_nonzero_num = min_nonzero_num,
                                    maxiter=maxiter
                                    )

    simsrt@metaParam$simType <- "RefSRT" 

    return(simsrt)
}


#' Fit the marginal distributions for a single label
#' @param x     A matrix of gene expression count for the locations in the selected domain label
#' @param marginal     Specification of the types of marginal distribution.Default value is 'auto_choose' which chooses between ZINB, NB, ZIP, and Poisson by a likelihood ratio test (lrt),AIC and whether there is underdispersion.'zinb' will fit the ZINB model. If there is underdispersion, it will fit the Poisson model. If there is no zero at all or an error occurs, it will fit an NB model instead.'nb' fits the NB model and chooses between NB and Poisson depending on whether there is underdispersion. 'poisson' simply fits the Poisson model.'zip' fits the ZIP model and chooses between ZIP and Poisson by a likelihood ratio test
#' @param min_nonzero_num  The minimum number of non-zero values required for a gene to be fitted. Default is 2.  
#' @param maxiter      The number of iterations for the model-fitting. Default is 500. 
#' @return Returns a list with estimated parameters for all the genes
#' @importClassesFrom S4Vectors SimpleList
#' @importFrom Matrix rowSums
#' @importFrom methods as
#' @noRd
#' @keywords internal

fit_single <- function(x, 
                marginal = c('auto_choose', 'zinb', 'nb', 'poisson','zip'), 
                min_nonzero_num = 2,
                maxiter =500){

    marginal <- match.arg(marginal)
    n <- ncol(x)
    p <- nrow(x)

    gene_zero_prop <- 1- Matrix::rowSums(x > 0) / n

    gene_sel1 <- which(gene_zero_prop < 1.0 - min_nonzero_num/n)
    gene_sel2 <- (1:p)[-gene_sel1]

    gene_names <- rownames(x)
    names(gene_sel1) <- gene_names[gene_sel1]
    names(gene_sel2) <- gene_names[gene_sel2]

    if(ncol(x)<10000){
        if(length(gene_sel1) > 0){
            marginal_result1 <- getparams(x[gene_sel1, ],marginal=marginal,maxiter=maxiter)
        }else{
            marginal_result1 = NULL
        }
    }else{
        ## split into small chunks to save the memory and improve the efficiency
        if(length(gene_sel1) > 0){
            xx <- x[gene_sel1, ]
            numChunks <- nrow(xx)%/%100
            marginal_result1 <- c()
            for(ichunk in 1:numChunks){
                if(ichunk==numChunks){
                    chunk_res <- getparams(xx[(1+100*(ichunk-1)):nrow(xx), ],marginal=marginal,maxiter=maxiter)
                }else{
                    chunk_res <- getparams(xx[(1+100*(ichunk-1)):(100*ichunk), ],marginal=marginal,maxiter=maxiter)
                }
                marginal_result1 <- rbind(marginal_result1,chunk_res)
            }
            rm(xx)
        }else{
            marginal_result1 = NULL
        }
    }

    marginal_result1 <- as(marginal_result1,"DFrame")

    return(S4Vectors::SimpleList(marginal_param1 = marginal_result1,
                gene_sel1 = gene_sel1, gene_sel2 = gene_sel2,
                min_nonzero_num = min_nonzero_num, 
                n_cell = n, n_read = sum(x)))
}


#' Fit the marginal distributions for single gene
#' @param x A matrix of gene expression count for the locations in the selected domain label
#' @param marginal Specification of the types of marginal distribution.Default value is 'auto_choose' which chooses between ZINB, NB, ZIP, and Poisson by a likelihood ratio test (lrt),AIC and whether there is underdispersion.'zinb' will fit the ZINB model. If there is underdispersion, it will fit the Poisson model. If there is no zero at all or an error occurs, it will fit an NB model instead.'nb' fits the NB model and chooses between NB and Poisson depending on whether there is underdispersion. 'poisson' simply fits the Poisson model.'zip' fits the ZIP model and chooses between ZIP and Poisson by a likelihood ratio test.
#' @param pval_cutoff  Cutoff of p-value of the lrt that determines whether there is zero inflation. Default is 0.05
#' @param maxiter      The number of iterations for the model-fitting. Default is 100 
#' @return Returns a n by 3 matrix with estimated parameters for all the genes
#' @importFrom stats pchisq
#' @importFrom dplyr mutate case_when
#' @importFrom magrittr %<>% 
#' @noRd
#' @keywords internal

getparams <- function(x, marginal = c('auto_choose', 'zinb', 'nb', 'poisson','zip'),pval_cutoff = 0.05,maxiter=100){
    p <- nrow(x)
    n <- ncol(x)
## parameter estimation
## parallel possible here, which will increase the efficiency
    marginal <- match.arg(marginal)
    if(marginal == 'auto_choose'){
        params <- t(apply(x, 1, function(gene){
            m <- mean(gene)
            v <- stats::var(gene)
            if(m >= v){
                c(0.0, Inf, m)
            }else{
                mle_NB <- fit_nb_optim(gene,maxiter=maxiter)
                # if(min(gene) > 0)
                ## if the empirical zero proportion is smaller than expected zero proportion 
                if(sum(gene==0)/length(gene) < (mle_NB[1]/(mle_NB[1]+mle_NB[2]))^{mle_NB[1]})
                    c(0.0, mle_NB[1], mle_NB[2])
                else
                    tryCatch({
                        mle_ZINB <- fit_zinb_optim(gene,maxiter=maxiter)
                        chisq_val <- 2 * (mle_ZINB[4] - mle_NB[3])
                        pvalue <- as.numeric(1 - pchisq(chisq_val, 1))
                        if(pvalue < pval_cutoff)
                            c(mle_ZINB[1], mle_ZINB[2], mle_ZINB[3])
                        else
                            tryCatch({
                                mle_ZIP   <- fit_zip_optim(gene,maxiter=maxiter)
                                aic_diff  <- 2 * (mle_ZIP[3]-mle_NB[3])
                                if(aic_diff > 2) 
                                    ## Model Selection and Multimodel Inference: A Practical Information-Theoretic Approach
                                    ## page 70, chapter 2.6, 4-7 considerable less
                                  c(mle_ZIP[1], Inf, mle_ZIP[2])
                                else
                                  c(0.0, mle_NB[1], mle_NB[2])
                            },
                            error = function(cond){
                              c(0.0, mle_NB[1], mle_NB[2])
                            })},
                        error = function(cond){
                        c(0.0, mle_NB[1], mle_NB[2])
                    })
            }
        }))
    }else if(marginal == 'zinb'){
        params <- t(apply(x, 1, function(gene){
            m <- mean(gene)
            v <- stats::var(gene)
            if(m >= v){
                c(0.0, Inf, m)
            }else{
                if(min(gene) > 0){
                  mle_NB <- fit_nb_optim(gene,maxiter=maxiter)
                  c(0.0, mle_NB[1], mle_NB[2])
                }
                else
                  tryCatch({
                    mle_ZINB <- fit_zinb_optim(gene,maxiter=maxiter)
                    c(mle_ZINB[1], mle_ZINB[2], mle_ZINB[3])
                  },
                  error = function(cond){
                    mle_NB <- fit_nb_optim(gene,maxiter=maxiter)
                    c(0.0, mle_NB[1], mle_NB[2])
                })
            }
        }))
    }else if(marginal == 'nb'){
        params <- t(apply(x, 1, function(gene){
          m <- mean(gene)
          v <- stats::var(gene)
          if(m >= v){
            c(0.0, Inf, m)
          }else{
            mle_NB <- fit_nb_optim(gene,maxiter=maxiter)
            c(0.0, mle_NB[1], mle_NB[2])
          }
        }))
    }else if(marginal == 'poisson'){
        params <- t(apply(x, 1, function(gene){
            c(0.0, Inf, mean(gene))
        }))
    }else if(marginal == 'zip'){
        params <- t(apply(x, 1, function(gene){
          m <- mean(gene)
          # v <- var(gene)
          mle_Poisson <- fit_pos_optim(gene,maxiter=maxiter)
          tryCatch({
            mle_ZIP   <- fit_zip_optim(gene,maxiter=maxiter)
            chisq_val <- 2 * (mle_ZIP[3] - mle_Poisson[2])
            pvalue <- as.numeric(1 - pchisq(chisq_val, 1))
            if(pvalue < pval_cutoff)
              c(mle_ZIP[1], Inf, mle_ZIP[2])
            else
              c(0.0, Inf, m)
          },
          error = function(cond){
            c(0.0, Inf, m)})
        }))
    }


    colnames(params) = c("pi0","theta","mu")
    params <- as.data.frame(params)
    params %<>% mutate(model_selected = case_when(
        pi0 == 0 & theta==Inf ~ "Poisson",
        pi0 == 0 & theta!=Inf ~ "NB",
        pi0 != 0 & theta!=Inf ~ "ZINB",
        pi0 != 0 & theta==Inf ~ "ZIP"
    ))

    # return(list(params = params))
    return(params)
}

