
#' Generate Data with Cell-Cell Interaction Under Reference-Free Mode
#' @param zero_prop_in    A number specifying zero proportion for the count model, default is 0
#' @param disper_in    A number specifying dispersion for the count model, default is Inf. Same as the size parameter in rnbinom.
#' @param mu_in    A number specifying mean for background count model, default is 1
#' @param numGene    An integer specifying the number of genes in the synthetic data, default is 1000
#' @param location_in  A dataframe with x, y, and region_label
#' @param region_cell_map  A dataframe specifying the cell type proportion in each region. Row: region,Column: cell type. 
#' @param fc 	A number specifying effect size for ligand-receptor pairs that mediate the cel-cell communication, default is 3
#' @param LR_in 	A dataframe specifying ligand and receptor pairs, containing four columns: protein_a, protein_b, celltypeA, and celltype B
#' @param sim_seed 	A number for reproducible purpose
#' @param numKNN 	A number specifying number of nearest neighbors with elevated gene expressin levels, default is 4
#' @param numSingleCellType 	A number specifying number of spots in the background pool. Gene expression count are then sampled from this background pool.
#' @return Returns a SRTsim object with a newly generated count matrix and correspoding parameters
#' @importFrom FNN get.knn
#' @importFrom stats na.omit
#' @export 

srtsim_cci_free = function(
	# param_in = c(0,1,1),
	zero_prop_in = 0,
	disper_in = Inf,
	mu_in = 1, 
	numGene = 1000,
	location_in,
	region_cell_map,
	fc=3,
	LR_in,
	sim_seed = 1,
	numKNN = 4,
	numSingleCellType = 2000){	


	param_in <- c(zero_prop_in,disper_in,mu_in)

	set.seed(sim_seed)
	numLoc 	= nrow(location_in)
	region_label = location_in$region_label

	numEntry = numSingleCellType*numGene # 5000 gene x 2000 cells to be sampled for each cell type

	# assign cell types to regions
	print("Generate Gene Expression Data...")
	celltype_count_in = list()
	for(celltype in colnames(region_cell_map)){
		count_simu = rbinom(numEntry, 1, 1-param_in[1]) * rnbinom(numEntry, size = param_in[2], mu = param_in[3])
		celltype_count_in[[celltype]] = matrix(count_simu, numGene, numSingleCellType)
		rm(celltype)
	}

	# assign cell types to regions
	print("Generate Cell Types...")
	ztrue 	<- region_label
	c_simu 	<- rep(NA, length(ztrue))

	for(z in unique(region_label)){
		zi_idx <- ztrue == z # zi_idx is index for region z
		c_simu[zi_idx] <- sample(colnames(region_cell_map), sum(zi_idx), prob = region_cell_map[z,], replace = T)
	}

	# assign count data
	sim_cnt <- array(NA, c(numGene, numLoc))
	for(ct in unique(c_simu)){
		c_size <- sum(c_simu == ct)  # true number of cells in cell type c
		if(dim(celltype_count_in[[ct]])[2]<c_size){
			stop("Targeted cell number of ", ct, " is greater than prespecified cell number in background pool. Increase the numSingleCellType value")
		}else{
			cells_select <- sample(1:dim(celltype_count_in[[ct]])[2], c_size, replace = F)
		}
		
		# sample same number of group c cells in real data from generated cells
		sim_cnt[, c_simu == ct] <- as.matrix(celltype_count_in[[ct]][, cells_select])
		# for positions of original cell type c cells, assign generated counts of group c
	}
	colnames(sim_cnt) <- paste0("Cell" ,1:numLoc)
	rownames(sim_cnt) <- paste0("Gene" ,1:numGene)

	knn.list = get.knn(as.matrix(location_in[,c("x","y")]), k = numKNN)[[1]]
	LR_names = unique(unlist(LR_in[,1:2]))
	simu_count = sim_cnt
	rownames(simu_count)[1:length(LR_names)] = LR_names

	# assign cell types to regions
	print("Construct Communications...")
	for(celltype_A in unique(LR_in$celltypeA)){
		print(celltype_A)
		cellid_A = which(c_simu==celltype_A)
		for(cellid in cellid_A){
			nearest_cells = knn.list[cellid,]
			for(nearid in nearest_cells){
				nearest_celltypeB = c_simu[nearid]
				tmp_LR = LR_in[which(LR_in$celltypeA == celltype_A & LR_in$celltypeB %in% nearest_celltypeB),]
				gene_ind=na.omit(match(c(unlist(tmp_LR[,1:2])),rownames(simu_count)))
				if(length(gene_ind)>0)
				simu_count[gene_ind,c(cellid,nearid)] = fc*sim_cnt[gene_ind,c(cellid,nearid)]
			}
		}
	}

	simInfo <- cbind.data.frame(location_in,celltype=c_simu)
  simsrt 	<- new(
    Class = "simSRT",
    simCounts = as(as.matrix(simu_count),"sparseMatrix"),
    simcolData = as(simInfo,"DataFrame"),
    metaParam = list(
    	simSeed = sim_seed,
    	simCCIKNN= numKNN,
    	simParam = list(background_param = param_in, effect_size = fc),
    	simType = "CCI")
  )
  return(simsrt)
}


#' Generate Data with Cell-Cell Interaction Under Reference-Based Mode
#' @param EstParam 	A list of estimated parameters from srtsim_fit function, EstParam slot if the simSRT object.
#' @param numGene 	An integer specifying the number of genes in the synthetic data, default is 1000
#' @param location_in  A dataframe with x, y, and region_label
#' @param region_cell_map  A dataframe specifying the cell type proportion in each region. Row: region,Column: cell type. 
#' @param fc 	A number specifying effect size for ligand-receptor pairs that mediate the cel-cell communication, default is 3
#' @param LR_in 	A dataframe specifying ligand and receptor pairs, containing four columns: protein_a, protein_b, celltypeA, and celltype B
#' @param sim_seed 	A number for reproducible purpose
#' @param numKNN 	A number specifying number of nearest neighbors with elevated gene expressin levels, default is 4
#' @param numSingleCellType 	A number specifying number of spots in the background pool. Gene expression count are then sampled from this background pool.
#' @return Returns a SRTsim object with a newly generated count matrix and correspoding parameters
#' @importFrom FNN get.knn
#' @importFrom stats na.omit
#' @export 

srtsim_cci_ref = function(
	EstParam = NULL,
	numGene = 1000,
	location_in,
	region_cell_map,
	fc=3,
	LR_in,
	sim_seed = 1,
	numKNN = 4,
	numSingleCellType = 2000){	

	set.seed(sim_seed)
	if(is.null(EstParam)){stop("EstParam is NULL, consider srtsim_cci_free for reference free simulation!")}
	if(!is.list(EstParam)){stop("EstParam has to be a list!!")}
	if(length(EstParam)< ncol(region_cell_map)){
		warning("The length of EstParam is less than targeted celltype number.")
		idx1 				<- sample(1:length(EstParam),ncol(region_cell_map),replace=TRUE)
		newEstParam <- EstParam[idx1]
		names(newEstParam) <- paste0("Celltype",1:ncol(region_cell_map))
	}else if(length(EstParam)> ncol(region_cell_map)){
		idx1 <- sample(1:length(EstParam),ncol(region_cell_map),replace=FALSE)
		newEstParam <- EstParam[idx1]
		names(newEstParam) <- paste0("Celltype",1:ncol(region_cell_map))
		print(paste0(paste(names(EstParam)[idx1],collapse=",")," are selected for data generation"))
	}else if(length(EstParam)== ncol(region_cell_map)){
		newEstParam <- EstParam
		names(newEstParam) <- paste0("Celltype",1:ncol(region_cell_map))
		print(paste0(names(EstParam)," are renamed to ",paste0("Celltype",1:ncol(region_cell_map))))
	}


	numLoc 	= nrow(location_in)
	region_label = location_in$region_label

	# numEntry = numSingleCellType*numGene # 5000 gene x 2000 cells to be sampled for each cell type

	# assign cell types to regions
	print("Generate Gene Expression Data...")
	celltype_count_in = list()
	for(celltype in colnames(region_cell_map)){
		celltype_param 	<- newEstParam[[celltype]]$marginal_param1
		idx_gene 				<- sample(1:nrow(celltype_param),numGene,replace=TRUE)
		tmp_parm 				<- celltype_param[idx_gene,1:3] ## ignore the model

		count_simu 			<- t(apply(tmp_parm,1,function(param_in){rbinom(numSingleCellType, 1, 1-param_in[1]) * rnbinom(numSingleCellType, size = param_in[2], mu = param_in[3])}))

		celltype_count_in[[celltype]] = count_simu
		rm(celltype,count_simu)
	}

	# assign cell types to regions
	print("Generate Cell Types...")
	ztrue 	<- region_label
	c_simu 	<- rep(NA, length(ztrue))

	for(z in unique(region_label)){
		zi_idx <- ztrue == z # zi_idx is index for region z
		c_simu[zi_idx] <- sample(colnames(region_cell_map), sum(zi_idx), prob = region_cell_map[z,], replace = T)
	}

	# assign count data
	sim_cnt <- array(NA, c(numGene, numLoc))
	for(ct in unique(c_simu)){
		c_size <- sum(c_simu == ct)  # true number of cells in cell type c
		if(dim(celltype_count_in[[ct]])[2]<c_size){
			stop("Targeted cell number of ", ct, " is greater than prespecified cell number in background pool. Increase the numSingleCellType value")
		}else{
			cells_select <- sample(1:dim(celltype_count_in[[ct]])[2], c_size, replace = F)
		}
		
		# sample same number of group c cells in real data from generated cells
		sim_cnt[, c_simu == ct] <- as.matrix(celltype_count_in[[ct]][, cells_select])
		# for positions of original cell type c cells, assign generated counts of group c
	}
	colnames(sim_cnt) <- paste0("Cell" ,1:numLoc)
	rownames(sim_cnt) <- paste0("Gene" ,1:numGene)

	knn.list = get.knn(as.matrix(location_in[,c("x","y")]), k = numKNN)[[1]]
	LR_names = unique(unlist(LR_in[,1:2]))
	simu_count = sim_cnt
	# LR_names = unique(unlist(LR_dat))
	rownames(simu_count)[1:length(LR_names)] = LR_names
	# mat_foldchange_record = matrix(0, length(LR_names), dim(simu_count)[2])

	# assign cell types to regions
	print("Construct Communications...")
	for(celltype_A in unique(LR_in$celltypeA)){
		print(celltype_A)
		cellid_A = which(c_simu==celltype_A)
		for(cellid in cellid_A){
			nearest_cells = knn.list[cellid,]
			for(nearid in nearest_cells){
				nearest_celltypeB = c_simu[nearid]
				tmp_LR = LR_in[which(LR_in$celltypeA == celltype_A & LR_in$celltypeB %in% nearest_celltypeB),]
				gene_ind=na.omit(match(c(unlist(tmp_LR[,1:2])),rownames(simu_count)))
				if(length(gene_ind)>0)
				simu_count[gene_ind,c(cellid,nearid)] = fc*sim_cnt[gene_ind,c(cellid,nearid)]
			}
		}
	}

	simInfo <- cbind.data.frame(location_in,celltype=c_simu)
  simsrt 	<- new(
    Class = "simSRT",
    simCounts = as(as.matrix(simu_count),"sparseMatrix"),
    simcolData = as(simInfo,"DataFrame"),
    metaParam = list(
    	simSeed = sim_seed,
    	simCCIKNN= numKNN,
    	simParam = list(background_param = newEstParam, effect_size = fc),
    	simType = "CCI_REF")
  )
  return(simsrt)
}
