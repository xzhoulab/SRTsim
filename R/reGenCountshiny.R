#' ReSimulate Count Data with Parameters Specification from Shiny
#' @param shinyOutput A list of Shiny Output. Including a simCount, simInfo,simcountParam,simLocParam
#' @param NewSeed A new seed for data generation. Useful when multiple replicates are needed.
#' @return Returns a Count DataFrame 
#' @importFrom stats rnbinom rbinom
#' 
#' @export
#' @examples
#'
#' ## Re-generate Count Data based on ShinyOutput Parameters, should be same as simCount in ShinyOutput
#' cMat <- reGenCountshiny(toyShiny)
#'
#' ## Generate Count Data with A New Seed based on ShinyOutput Parameters 
#' cMat2 <- reGenCountshiny(toyShiny,NewSeed=2)
#' 
#' ## Comparison across the output
#' toyShiny$simCount[1:3,1:3]
#' cMat[1:3,1:3]
#' cMat2[1:3,1:3]
#'

reGenCountshiny <- function(shinyOutput,NewSeed=NULL){
  
  if(is.list(shinyOutput)){
    ShinyParam <- shinyOutput$simcountParam
  }else{
    ShinyParam <- shinyOutput@metaParam$shinySRTParam$simcountParam
  }

  pattern_in    <- ShinyParam$locfile
  numHighSignal <- ShinyParam$numHighSignal
  numLowSignal  <- ShinyParam$numLowSignal
  numBG         <- ShinyParam$numNoise
  disper_in     <- ShinyParam$Dispersion
  zero_in       <- ShinyParam$zeroProp
  mu_in         <- ShinyParam$mu
  isid          <- ShinyParam$seed

  if(!is.null(NewSeed)){isid <- NewSeed}

  set.seed(isid)

  numSpots <- nrow(pattern_in)
  numGenes <- numHighSignal + numLowSignal + numBG

  if(numHighSignal!=0){
    Highsignal_df <- t(sapply(1:numHighSignal,function(x){
      rbinom(n=numSpots, 1, 1-zero_in) * rnbinom(n=numSpots, size = 1/disper_in, mu = pattern_in$foldchange*mu_in)
    }))

    rownames(Highsignal_df) <- c(paste0("signal",1:numHighSignal))
  }else{
    Highsignal_df <- NULL
  }

  if(numLowSignal!=0){
    Lowsignal_df <- t(sapply(1:numLowSignal,function(x){
      rbinom(n=numSpots, 1, 1-zero_in) * rnbinom(n=numSpots, size = 1/disper_in, mu = (1/pattern_in$foldchange)*mu_in)
    }))

    rownames(Lowsignal_df) <- c(paste0("signal",(numHighSignal+1):(numHighSignal+numLowSignal)))
  }else{
    Lowsignal_df <- NULL
  }

  if(numBG!=0){
    noise_df    <- t(sapply(1:numBG,function(x){
      rbinom(n=numSpots, 1, 1-zero_in) * rnbinom(n=numSpots,size=1/disper_in,mu=mu_in)})
    )

    rownames(noise_df) <- c(paste0("noise",1:numBG))
  }else{
    noise_df <- NULL
  }

  count_out <- rbind.data.frame(Highsignal_df,Lowsignal_df,noise_df)

  colnames(count_out) <- rownames(pattern_in)

  return(count_out)
}






