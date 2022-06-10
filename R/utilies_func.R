#' @importFrom spatstat.random rpoispp
get_new_coordinates <- function(ns,bg_shape,isid,lay_out,preLoc=NULL){
  set.seed(isid)
  if(bg_shape=="Circle"){
    radius  <- 1
    if(lay_out=="Random"){
      pp      <- spatstat.random::rpoispp(ns/(pi*radius^2),win=disc(radius,centre=rep(radius,2)))
      simLoc  <- cbind.data.frame(x=pp[['x']],y=pp[['y']],group=rep("A",pp[["n"]]))
    }else{
      shp     <- st_as_sf(disc(radius,centre=rep(radius,2)))
      nspots  <- area(st_bbox(shp))/area(shp)*ns ## approxmiate number of spots generated to the one required
      grid    <- shp %>% 
        st_make_grid(n=round(sqrt(nspots)), what = "centers") %>% # grid of points
        st_intersection(shp)                               # only within the polygon

      simLoc  <- as.data.frame(st_coordinates(grid))
      simLoc$group <- rep("A",nrow(simLoc))
    }
    
  }else if(bg_shape=="Square"){
    if(lay_out=="Random"){
      pp      <- spatstat.random::rpoispp(ns)
      simLoc  <- cbind.data.frame(x=pp[['x']],y=pp[['y']],group=rep("A",pp[["n"]]))
    }else{
      shp     <- st_as_sf(owin(xrange=c(0,1), yrange=c(0,1)))
      nspots  <- area(st_bbox(shp))/area(shp)*ns ## approxmiate number of spots generated to the one required
      grid    <- shp %>% 
        st_make_grid(n=round(sqrt(nspots)), what = "centers") %>% # grid of points
        st_intersection(shp)                               # only within the polygon

      simLoc  <- as.data.frame(st_coordinates(grid))
      simLoc$group <- rep("A",nrow(simLoc))
    }
  }else if(bg_shape=="User Define Shape"){
    if(!is.null(preLoc)){
      # if(length(grep("group",colnames(preLoc)))==0){
      # if(simNew){
        pnts    <- preLoc[,1:2] %>% st_as_sf(coords = c("x", "y"))
        polygon <- concaveman(pnts,2.0) 
        poly_coords = as.data.frame(as.matrix(polygon$polygons[[1]]))
        colnames(poly_coords) <- c("x","y")
        Pl      <- Polygon(poly_coords)
        if(lay_out=="Random"){
          pts     <- spsample(Pl,n=ns,"random")
          simLoc  <- data.frame(x=pts$x,y=pts$y,group=rep("A",nrow(pts@coords)))
        }else{
          shp     <- polygon
          nspots  <- area(st_bbox(shp))/area(shp)*ns ## approxmiate number of spots generated to the one required
          grid    <- shp %>% 
            st_make_grid(n=round(sqrt(nspots)), what = "centers") %>% # grid of points
            st_intersection(shp)                               # only within the polygon

          simLoc  <- as.data.frame(st_coordinates(grid))
          simLoc$group <- rep("A",nrow(simLoc))
        }
    }else{
      # warnings("Please provide a dataframe with spots location for the boundary estimation")
      print("No dataframe provided, simulate using the square backgound shape")
      pp        <- spatstat.random::rpoispp(ns)
      simLoc    <- cbind.data.frame(x=pp[['x']],y=pp[['y']],group=rep("A",pp[["n"]]))
    }
  }else if(bg_shape=="User Define Spots"){
      if(!is.null(preLoc)){
          if(length(grep("group",colnames(preLoc)))!=0){
            simLoc  <- data.frame(x=preLoc[,1],y=preLoc[,2],group=preLoc[,grep("group",colnames(preLoc))])
          }else{
            simLoc <- data.frame(x=preLoc[,1],y=preLoc[,2])
            simLoc$group <- rep("A",nrow(simLoc))
          }
   
      }else{
        # warnings("Please provide a dataframe with spots location for the boundary estimation")
        print("No dataframe provided, simulate using the square backgound shape")
        pp        <- spatstat.random::rpoispp(ns)
        simLoc    <- cbind.data.frame(x=pp[['x']],y=pp[['y']],group=rep("A",pp[["n"]]))
      }
  }

  colnames(simLoc) <- c("x","y","group")
  rownames(simLoc) <- paste0("Loc",1:nrow(simLoc))
  simLoc$foldchange <- rep(1,nrow(simLoc))
  # simPara       <- param_tech_func(technique)
  # simLoc$mu0    <- rep(simPara$mu,nrow(simLoc))
  # simLoc$theta0 <- rep(simPara$theta,nrow(simLoc))
  # simLoc$technique <- rep(technique,nrow(simLoc))
  return(simLoc)
}

#' @importFrom spatstat.random rpoispp
get_new_coordinates2 <- function(ns,bg_shape,isid,lay_out,preLoc=NULL){
  set.seed(isid)
  if(bg_shape=="Circle"){
    radius  <- 1
    if(lay_out=="Random"){
      pp      <- spatstat.random::rpoispp(ns/(pi*radius^2),win=disc(radius,centre=rep(radius,2)))
      simLoc  <- cbind.data.frame(x=pp[['x']],y=pp[['y']],group=rep("A",pp[["n"]]))
    }else{
      shp     <- st_as_sf(disc(radius,centre=rep(radius,2)))
      nspots  <- area(st_bbox(shp))/area(shp)*ns ## approxmiate number of spots generated to the one required
      grid    <- shp %>% 
        st_make_grid(n=round(sqrt(nspots)), what = "centers") %>% # grid of points
        st_intersection(shp)                               # only within the polygon

      simLoc  <- as.data.frame(st_coordinates(grid))
      simLoc$group <- rep("A",nrow(simLoc))
    }
    
  }else if(bg_shape=="Square"){
    if(lay_out=="Random"){
      pp      <- spatstat.random::rpoispp(ns)
      simLoc  <- cbind.data.frame(x=pp[['x']],y=pp[['y']],group=rep("A",pp[["n"]]))
    }else{
      shp     <- st_as_sf(owin(xrange=c(0,1), yrange=c(0,1)))
      nspots  <- area(st_bbox(shp))/area(shp)*ns ## approxmiate number of spots generated to the one required
      grid    <- shp %>% 
        st_make_grid(n=round(sqrt(nspots)), what = "centers") %>% # grid of points
        st_intersection(shp)                               # only within the polygon

      simLoc  <- as.data.frame(st_coordinates(grid))
      simLoc$group <- rep("A",nrow(simLoc))
    }
  }else if(bg_shape=="User Define Shape"){
    if(!is.null(preLoc)){
      # if(length(grep("group",colnames(preLoc)))==0){
      # if(simNew){
        pnts    <- preLoc[,1:2] %>% st_as_sf(coords = c("x", "y"))
        polygon <- concaveman(pnts,2.0) 
        poly_coords = as.data.frame(as.matrix(polygon$polygons[[1]]))
        colnames(poly_coords) <- c("x","y")
        Pl      <- Polygon(poly_coords)
        if(lay_out=="Random"){
          pts     <- spsample(Pl,n=ns,"random")
          simLoc  <- data.frame(x=pts$x,y=pts$y,group=rep("A",nrow(pts@coords)))
        }else{
          shp     <- polygon
          nspots  <- area(st_bbox(shp))/area(shp)*ns ## approxmiate number of spots generated to the one required
          grid    <- shp %>% 
            st_make_grid(n=round(sqrt(nspots)), what = "centers") %>% # grid of points
            st_intersection(shp)                               # only within the polygon

          simLoc  <- as.data.frame(st_coordinates(grid))
          simLoc$group <- rep("A",nrow(simLoc))
        }
    }else{
      # warnings("Please provide a dataframe with spots location for the boundary estimation")
      print("No dataframe provided, simulate using the square backgound shape")
      pp        <- spatstat.random::rpoispp(ns)
      simLoc    <- cbind.data.frame(x=pp[['x']],y=pp[['y']],group=rep("A",pp[["n"]]))
    }
  }else if(bg_shape=="User Define Spots"){
      if(!is.null(preLoc)){
          if(length(grep("group",colnames(preLoc)))!=0){
            simLoc  <- data.frame(x=preLoc[,1],y=preLoc[,2],group=preLoc[,grep("group",colnames(preLoc))])
          }else{
            simLoc <- data.frame(x=preLoc[,1],y=preLoc[,2])
            simLoc$group <- rep("A",nrow(simLoc))
          }
   
      }else{
        # warnings("Please provide a dataframe with spots location for the boundary estimation")
        print("No dataframe provided, simulate using the square backgound shape")
        pp        <- spatstat.random::rpoispp(ns)
        simLoc    <- cbind.data.frame(x=pp[['x']],y=pp[['y']],group=rep("A",pp[["n"]]))
      }
  }

  colnames(simLoc) <- c("x","y","group")
  rownames(simLoc) <- paste0("Loc",1:nrow(simLoc))
  simLoc$foldchange <- rep(1,nrow(simLoc))

  simLocParam <- list(numloc = ns,
                      shape = bg_shape,
                      seed = isid,
                      lay_out=lay_out,
                      inputLoc = preLoc)
  return(list(simLoc=simLoc,simLocParam=simLocParam))
}



pattern_count_func <- function(pattern_in,numSignal=100,numBG=0,disper_in=0,zero_in=0,mu_in=0.5,verbose=FALSE,isid=1){
  
  set.seed(isid)

  numSpots <- nrow(pattern_in)
  numGenes <- numSignal + numBG

  if(verbose) {message("Simulating Signal Genes...")}
  if(numSignal!=0){

    signal_df <- t(sapply(1:numSignal,function(x){
      rbinom(n=numSpots, 1, 1-zero_in) * rnbinom(n=numSpots, size = 1/disper_in, mu = pattern_in$foldchange*mu_in)
    }))

    rownames(signal_df) <- c(paste0("signal",1:numSignal))

  }else{
    signal_df <- NULL
  }

  if(numBG!=0){
    if (verbose) {message("Simulating Noise Genes...")}
    noise_df    <- t(sapply(1:numBG,function(x){
      rbinom(n=numSpots, 1, 1-zero_in) * rnbinom(n=numSpots,size=1/disper_in,mu=mu_in)})
    )

    rownames(noise_df) <- c(paste0("noise",1:numBG))
    count_out   <- rbind.data.frame(signal_df,noise_df)
    # rownames(count_out) <- c(paste0("signal",1:numSignal),paste0("noise",1:numBG))
  }else{
    count_out <- signal_df
    # rownames(count_out) <- c(paste0("signal",1:numSignal))
  }
  colnames(count_out) <- rownames(pattern_in)
  return(count_out)
}


pattern_count_func2 <- function(pattern_in,numHighSignal=100,numLowSignal=0,numBG=0,disper_in=0,zero_in=0,mu_in=0.5,isid=1){
  
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

    # rownames(Lowsignal_df) <- c(paste0("lowsignal",1:numLowSignal))
    rownames(Lowsignal_df) <- c(paste0("signal",(numHighSignal+1):(numHighSignal+numLowSignal)))
  }else{
    Lowsignal_df <- NULL
  }

  if(numBG!=0){
    noise_df    <- t(sapply(1:numBG,function(x){
      rbinom(n=numSpots, 1, 1-zero_in) * rnbinom(n=numSpots,size=1/disper_in,mu=mu_in)})
    )

    rownames(noise_df) <- c(paste0("noise",1:numBG))
    # count_out   <- rbind.data.frame(signal_df,noise_df)
    # rownames(count_out) <- c(paste0("signal",1:numSignal),paste0("noise",1:numBG))
  }else{
    noise_df <- NULL
    # count_out <- signal_df
    # rownames(count_out) <- c(paste0("signal",1:numSignal))
  }

  count_out <- rbind.data.frame(Highsignal_df,Lowsignal_df,noise_df)

  colnames(count_out) <- rownames(pattern_in)

  simcountParam <- list(locfile=pattern_in,
                        numHighSignal=numHighSignal,
                        numLowSignal=numLowSignal,
                        numNoise = numBG,
                        Dispersion = disper_in,
                        zeroProp = zero_in,
                        mu = mu_in,
                        seed = isid
                        )
  return(list(ctdf=count_out,simcountParam=simcountParam))
}



relative_func <- function(rx){
  rexpr   = (rx-min(rx))/(max(rx)-min(rx))
  return(rexpr)
}


