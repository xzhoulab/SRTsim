
#' log-likelihood for two-parameter negative binomial
#' @param param A vector of two parameters in the negative binomial distribution
#' @param y A vector of count values to be fitted
#' @return Returns a log-likelihood value
#' 
#' @noRd
#' @keywords internal

nb_loglik <- function(param,y){
    ## this version is two parameter estimation,slower than the fix mu estimation.
    ## but the result is same as the fitdistr

    # th = param[1]
    # mu = param[2]
    th = max(param[1],1e-10)
    mu = max(param[2],1e-10)

    return(-sum((lgamma(th + y) - lgamma(th) - lgamma(y + 1) + th * log(th) + y * 
        log(mu + (y == 0)) - (th + y) * log(th + mu))))
}


#' log-likelihood for two-parameter negative binomial but fix the mu
#' @param th A dispersion parameter to be estimated in the negative binomial distribution
#' @param mu A mean parameter being fixed in the estimation
#' @param y A vector of count values to be fitted
#' @return Returns a log-likelihood value
#' 
#' @noRd
#' @keywords internal

nb_loglik_mom <- function(th,mu, y){

    return(-sum((lgamma(th + y) - lgamma(th) - lgamma(y + 1) + th * log(th) + y * 
        log(mu + (y == 0)) - (th + y) * log(th + mu))))
}


#' fitting data with nb through optim function
#' @param x A vector of count values to be fitted
#' @param maxiter number of iteration
#' @return Returns a vector with dispersion theta, mean parameter mu, loglikelihood value llk, convergence
#' @importFrom MASS glm.nb
#' @importFrom stats logLik optim
#' @noRd
#' @keywords internal

fit_nb_optim <- function(x,maxiter= 500){
    m = mean(x)
    v = stats::var(x)

    if(v>m){
    size = m^2/(v-m)
    }else{
    size = 100
    }


    if(length(x)<1000){
        param2est <- c(size,m)
        ## the BFGS method is more stable than Nelder-Mead? need to figure out which is fast and which is better
        ## optim may give error in some case of 10X, similar problem observed in the MASS (fitdistr)
        fitres <- tryCatch({
                    opt_out <- optim(param2est,nb_loglik,y=x,control = list(maxit = maxiter),method="BFGS")
                    c(opt_out$par,-opt_out$value,1-opt_out$convergence)
                  },
                  error = function(cond){
                    # library(MASS)
                    glmfit <- glm.nb(x~1)
                    c(glmfit$theta, exp(glmfit$coefficients),as.numeric(logLik(glmfit)),min(glmfit$converged,is.null(glmfit$th.warn))) 
                  })
    }else{
        ## fix the mu through the moment matching to facilitate the estimation
        fitres <- tryCatch({
              opt_out <- optim(size,nb_loglik_mom,mu=m,y=x,control = list(maxit = maxiter), method="Brent",lower=0,upper=1000)
              c(opt_out$par,m,-opt_out$value,1-opt_out$convergence)
            },
            error = function(cond){
              # library(MASS)
              glmfit <- glm.nb(x~1)
              c(glmfit$theta, exp(glmfit$coefficients),as.numeric(logLik(glmfit)),min(glmfit$converged,is.null(glmfit$th.warn))) 
            })
    }

    names(fitres) <- c("theta","mu","llk","convergence")
    return(fitres)
}



#' log-likelihood for poisson distribution
#' @param mu the mean parameter in the poisson distribution
#' @param y A vector of count values to be fitted
#' @return Returns a log-likelihood value
#' 
#' @noRd
#' @keywords internal

pos_loglik <- function(mu,y){
	# th = param[1]
	# mu = param[2]
	return(-sum(y*log(mu) - mu - lgamma(y + 1)))
}

#' fitting data with poisson through optim function
#' @param x A vector of count values to be fitted
#' @param maxiter number of iteration
#' @return Returns a vector with mean paramter lambda, loglikelihood value llk, convergence
#' @importFrom stats optim

fit_pos_optim <- function(x,maxiter= 500){
  ## need to improve for the sparsematrix
	m = mean(x)

	opt_out <- optim(m,pos_loglik,y=x,control = list(maxit = maxiter), method="Brent",lower=0,upper=max(x))
	fitres <- c(opt_out$par,-opt_out$value,1-opt_out$convergence)
	names(fitres) <- c("lambda","llk","convergence")
	return(fitres)
}


#' log-likelihood for zero inflated poisson distribution
#' @param lam the mean parameter in the poisson distribution
#' @param y A vector of count values to be fitted
#' @return Returns a log-likelihood value
#'
#' @noRd
#' @keywords internal

zip_loglik <- function(lam,y){
    zeroidx <- which(y==0)
    nzero   <- length(zeroidx)
    numSam  <- length(y)

    if(nzero>0){
        nz_y <- y[-zeroidx]
    }else{
        nz_y <- y
    }

    pp <- (nzero/numSam-exp(-lam))/(1-exp(-lam))

    ll <- nzero*log(pp + (1-pp)*exp(-lam)) + sum(log(1-pp) + nz_y*log(lam) - lam -lgamma(nz_y + 1))

    return(-(ll))
}



#' fitting data with zip through optim function
#' @param x A vector of count values to be fitted
#' @param maxiter number of iteration
#' @return Returns a vector with zero proportion p0, mean parameter mu, loglikelihood value llk, convergence, and if zip
#' @importFrom stats optim
#' @noRd
#' @keywords internal

fit_zip_optim <- function(x,maxiter= 500){

    ## need to improve for the sparsematrix
    zeroidx <- which(x==0)
    n0    <- length(zeroidx)
    n     <- length(x)
    m     <- mean(x)

    opt_out <- optim(m,zip_loglik,y=x,control = list(maxit = maxiter), method="Brent",lower=min(-log(n0/n),0),upper=max(x))

    pp_est  <- (n0/n-exp(-opt_out$par))/(1-exp(-opt_out$par))

    if(pp_est<=0|pp_est==1|opt_out$convergence==1){
        opt_out   <- fit_pos_optim(x,maxiter= maxiter)
        fitres    <- c(0.0,opt_out,0)
    }else{
        fitres  <- c(pp_est,opt_out$par,-opt_out$value,1-opt_out$convergence,1)
    }
    names(fitres) <- c("p0","mu","llk","convergence","zip")
    return(fitres)
}



#' log-likelihood for zero inflated negative binomial distribution
#' @param par_init A vector of three parameters in the zero inflated negative binomial distribution
#' @param y A vector of count values to be fitted
#' @return Returns a log-likelihood value
#'
#' @noRd
#' @keywords internal


zinb_loglik <- function(par_init,y){

    th      <- par_init[1]
    # mu      <- par_init[2]
    mu      <- max(par_init[2],1e-10)
    rr      <- par_init[3]

    zeroidx <- which(y==0)
    nzero   <- length(zeroidx)
    numSam  <- length(y)

    if(nzero>0){
        nz_y <- y[-zeroidx]
    }else{
        nz_y <- y
    }

    # pp <- (nzero/numSam-(th/(mu+th))^th)/(1-(th/(mu+th))^th)
    pp  <- exp(rr)/(1+exp(rr))

    ll1 <- nzero*log(pp + (1-pp)*(th/(mu+th))^th) 
    ll2 <- sum(log(1-pp)+(lgamma(th + nz_y) - lgamma(th) - lgamma(nz_y + 1) + th * log(th) + nz_y * 
        log(mu) - (th + nz_y) * log(th + mu)))
    ll  <- ll1 + ll2

    return(-ll)
}


#' fitting data with zip through optim function
#' @param x A vector of count values to be fitted
#' @param maxiter number of iteration
#' @return Returns a vector with zero proportion p0, dispersion parameter theta, mean parameter mu, loglikelihood value llk, convergence, and if zinb
#' @importFrom stats optim
#' @noRd
#' @keywords internal

fit_zinb_optim <- function(x,maxiter= 500){

    ## need to improve for the sparsematrix

    zeroidx <- which(x==0)
    n0    <- length(zeroidx)
    n     <- length(x)
    m     <- mean(x)
    v     <- stats::var(x)
    if(v>m){
        size = m^2/(v-m)
    }else{
        size = 1
    }

    if(n0!=0){
        param_init <- c(size,m,log(n0/(n-n0)))

        opt_out <- optim(param_init,zinb_loglik,y=x,control = list(maxit = maxiter), method="Nelder-Mead")
        pp_est <- exp(opt_out$par[3])/(1+exp(opt_out$par[3]))

        # if not converge or the proportion is problematic, then switch to NB
        if(pp_est<=0|pp_est==1|opt_out$convergence==1){
          opt_out   <- fit_nb_optim(x,maxiter= maxiter)
          fitres    <- c(0.0,opt_out,0)
        }else{
          fitres    <- c(pp_est,opt_out$par[1],opt_out$par[2],-opt_out$value,1-opt_out$convergence,1)
        }
    }else{
        opt_out   <- fit_nb_optim(x,maxiter= maxiter)
        fitres    <- c(0.0,opt_out,0)
    }

    names(fitres) <- c("p0","theta","mu","llk","convergence","zinb")
    return(fitres)
}











