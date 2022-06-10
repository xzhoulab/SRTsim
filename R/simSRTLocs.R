

#' Fit the marginal distributions for each row of a count matrix
#' @param simsrt    A SRTsim object
#' @param new_loc_num  A integer specifying the number of spatial locations in the synthetic data
#' @param voting_nn    A integer of nearest neighbors used in label assignment for new generated locations. Default is 3.  
#' @param loc_lay_out  a character string specifying arrangement of new generated spatial locations. Default is "grid"
#' @return Returns a object with estimated parameters
#' 
#' @export 

srtsim_newlocs <- function(simsrt,
                            new_loc_num = NULL,
                            loc_lay_out = c("grid","random"),
                            voting_nn = 3){


    locfile <- as.data.frame(simsrt@refcolData[,c("x","y","label")])
    ## may consider put this into simulate functions
    ## as one time model fit is enough for multiple sample size simulation
    if(!is.null(new_loc_num)){
        loc_lay_out     <- match.arg(loc_lay_out)
        new_loc_df      <- simNewLocs(newN = new_loc_num,
                                        lay_out=loc_lay_out, 
                                        preLoc=locfile)

        ## calculate pairwise distance 
        new_old_dmat    <- t(apply(new_loc_df,1,function(x){as.matrix(pdist(locfile[,c("x","y")],x))}))

        ## save nn_k50 matrix for further label assignment and data simulation.
        nn_k50          <- t(apply(new_old_dmat,1,function(x){order(x)[1:50]}))

        ## create labels for the new locations
        if(voting_nn==1){
            voting_knn_res    <- nn_k50[,1]
            new_loc_df$label  <- locfile$label[voting_knn_res]
        }else{
            voting_knn_res    <- nn_k50[,1:voting_nn]
            new_loc_df$label  <- apply(voting_knn_res,1,function(ii){vote_func(locfile$label[ii])})
        }  
        rownames(new_loc_df) <- paste0("loc",1:nrow(new_loc_df))

        if(loc_lay_out=="grid"){
            new_loc_df$x_grid <- convert_grid(new_loc_df$x)
            new_loc_df$y_grid <- convert_grid(new_loc_df$y)
        }

        new_loc_df  <- as(new_loc_df,"DFrame")
    }else{
        new_loc_df  <- as(locfile,"DFrame")
        nn_k50      <- NULL
    }


    simsrt@metaParam$simLocParam <- SimpleList(voting_nn=voting_nn,
                                        loc_lay_out = loc_lay_out,
                                        new_loc_num = new_loc_num,
                                        simref_nmat = nn_k50
                                        )

    simsrt@simcolData <- new_loc_df

    return(simsrt)
}




#' Majority voting for the labels
#' @param x A vector of labels from selected neighbors 
#' @return  Returns a label with highest frequency
#' @noRd
#' @keywords internal

vote_func <- function(x){
    tabx <- table(x)
    return(names(which.max(rank(tabx,ties.method="random"))))
}


#' Fit the marginal distributions for single gene
#' @param newN     A integer specifying the number of spatial locations in the synthetic data
#' @param lay_out  A character string specifying arrangement of new generated spatial locations. Default is "grid"
#' @param preLoc   A data frame of shape n by 3 that x, y coodinates and domain label
#' @return Returns a n by 2 dataframe with newly generated spatial locations
#' @importFrom spatstat.geom area
#' @noRd
#' @keywords internal

simNewLocs <- function(newN,lay_out=c("grid","random"),preLoc){
    lay_out <- match.arg(lay_out)
    pnts    <- preLoc[,1:2] %>% st_as_sf(coords = c("x", "y"))
    polygon <- concaveman(pnts,2.0) 
    poly_coords <- as.data.frame(as.matrix(polygon$polygons[[1]]))
    colnames(poly_coords) <- c("x","y")
    Pl      <- Polygon(poly_coords)
    if(lay_out=="random"){
        pts     <- spsample(Pl,n=newN,"random")
        simLoc  <- data.frame(x=pts$x,y=pts$y)
    }else{
        shp     <- polygon
        nspots  <- area(st_bbox(shp))/area(shp)*newN ## approxmiate number of spots generated to the one required
        grid    <- shp %>% 
        st_make_grid(n=round(sqrt(nspots)), what = "centers") %>% # grid of points
        st_intersection(shp)                               # only within the polygon
        simLoc  <- as.data.frame(st_coordinates(grid))
    }
    colnames(simLoc) <- c("x","y")
    return(simLoc)
}


#' Convert continuous coordinate into integer, essential for BayesSpace to determine the neighborhood info
#' @param x A numeric vector of continuous coordinate
#' @return Returns a numeric vector oof integer coordinate
#' @export
convert_grid <- function(x){
    uniqx <- unique(x)
    uniqx_order <- uniqx[order(uniqx)]
    return(as.numeric(factor(x,levels=uniqx_order)))
}














