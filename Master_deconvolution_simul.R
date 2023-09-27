transformation = "none"
deconv_type = "sc"
normalization_scC = "TMM"
normalization_scT = "TMM"
method = "SCDC"
#install.packages("xbioc")
number_cells = round(as.numeric(100), digits = -2) #has to be multiple of 100
to_remove = "none"
num_cores = min(as.numeric(1),parallel::detectCores()-1)


### Helper functions + CIBERSORT external code
source('./helper_functions.R')
# source('CIBERSORT.R')

#-------------------------------------------------------
### Read data and metadata
#data = readRDS(list.files(path = "./example", pattern = "rds", full.names = TRUE)) # row: gene expression, column: cellID
#full_phenoData = read.table(list.files(path = "./example", pattern = "phenoData", full.names = TRUE), header=TRUE)

data=readRDS("baron_scRNA_cellAllType_simu1.rds")
full_phenoData=read.table("baron_phenoData_simu1.txt")

#-------------------------------------------------------
### QC 
require(dplyr); require(Matrix)

# First: cells with library size, mitochondrial or ribosomal content further than three MAD away were discarded
filterCells <- function(filterParam){
	cellsToRemove <- which(filterParam > median(filterParam) + 3 * mad(filterParam) | filterParam < median(filterParam) - 3 * mad(filterParam) )
	cellsToRemove
}

libSizes <- colSums(data)
gene_names <- rownames(data)

mtID <- grepl("^MT-|_MT-", gene_names, ignore.case = TRUE)
rbID <- grepl("^RPL|^RPS|_RPL|_RPS", gene_names, ignore.case = TRUE)

mtPercent <- colSums(data[mtID, ])/libSizes
rbPercent <- colSums(data[rbID, ])/libSizes

lapply(list(libSizes = libSizes, mtPercent = mtPercent, rbPercent = rbPercent), filterCells) %>% 
	unlist() %>% 
	unique() -> cellsToRemove

if(length(cellsToRemove) != 0){
	data <- data[,-cellsToRemove]
	full_phenoData <- full_phenoData[-cellsToRemove,]
} #dim(single_4): 20125 x 6042, dim(single_4'): 20125 x 6035

# Keep only "detectable" genes: at least 5% of cells (regardless of the group) have a read/UMI count different from 0
keep <- which(Matrix::rowSums(data > 0) >= round(0.05 * ncol(data)))
data = data[keep,] #dim(single_4):8761 x 6042

data_simu2=readRDS("baron_scRNA_cellAllType_simu2.rds")
full_phenoData_simu2=read.table("baron_phenoData_simu2.txt")

data_simu2 <- data_simu2[,-cellsToRemove]
full_phenoData_simu2 <- full_phenoData_simu2[-cellsToRemove,]
data_simu2 = data_simu2[keep,]


#-------------------------------------------------------
### Data split into training/test 
set.seed(24)
require(limma); require(dplyr); require(pheatmap)

colnames(full_phenoData) <- c("cellID", "cellType", "sampleID")
colnames(full_phenoData_simu2) <- c("cellID", "cellType", "sampleID")

original_cell_names = colnames(data) #length(original_cell_names)=6042
colnames(data) <- as.character(full_phenoData$cellType[match(colnames(data),full_phenoData$cellID)])  #for now let the colnames of data be cellType

original_cell_names_simu2 = colnames(data_simu2)

# Keep Cell Types  (CTs) with >= 50 cells after QC
cell_counts = table(colnames(data))
to_keep = names(cell_counts)[cell_counts >= 50]
pData <- full_phenoData[full_phenoData$cellType %in% to_keep,]
to_keep = which(colnames(data) %in% to_keep)   
data <- data[,to_keep]
original_cell_names <- original_cell_names[to_keep]


# Generate phenodata for reference matrix C
pDataC = pData

train <- data
test <- data_simu2

# "write.table" & "saveRDS" statements are optional, for users willing to avoid generation of matrix C every time:    
# write.table(pDataC, file = paste(dataset,"phenoDataC",sep="_"),row.names=FALSE,col.names=TRUE,sep="\t",quote=FALSE)

train_cellID = train
colnames(train_cellID) = original_cell_names
# saveRDS(object = train_cellID, file = paste(dataset,"qc_filtered_train.rds",sep="_")) #It has to contain cellID as colnames, not cellType (for scRNA-seq methods)
# saveRDS(object = test, file = paste(dataset,"qc_filtered_test.rds",sep="_"))

# reference matrix (C) + refProfiles.var from TRAINING dataset
cellType <- colnames(train)
group = list()
for(i in unique(cellType)){ 
	group[[i]] <- which(cellType %in% i)
}
C = lapply(group,function(x) Matrix::rowMeans(train[,x])) #C should be made with the mean (not sum) to agree with the way markers were selected #for cell type 1 and gene 1, mean(cell1+...+cell100+)
C = round(do.call(cbind.data.frame, C)) #the row of the C is gene, the column of C is cell type, if #gene=100, #cellType=5, then the dimension of C is 100*5
# write.table(C, file = "C",row.names=TRUE,col.names=TRUE,sep="\t",quote=FALSE,)

refProfiles.var = lapply(group,function(x) train[,x])
refProfiles.var = lapply(refProfiles.var, function(x) matrixStats::rowSds(Matrix::as.matrix(x)))
refProfiles.var = round(do.call(cbind.data.frame, refProfiles.var)) #the row of the refProfiles.var is gene, the column of refProfiles.var is cell type, if #gene=100, #cellType=5, then the dimension of refProfiles.var is 100*5
rownames(refProfiles.var) <- rownames(train)
# write.table(refProfiles.var, "refProfiles.var", quote=FALSE,row.names=TRUE,col.names=TRUE,sep="\t")

#-------------------------------------------------------
#Normalization of "train" followed by marker selection 

#for marker selection, keep genes where at least 30% of cells within a cell type have a read/UMI count different from 0
cellType = colnames(train) 
keep <- sapply(unique(cellType), function(x) {
	CT_hits = which(cellType %in% x)
 	size = ceiling(0.3*length(CT_hits))
 	Matrix::rowSums(train[,CT_hits,drop=FALSE] != 0) >= size
})
train = train[Matrix::rowSums(keep) > 0,]
train2 = Normalization(train)

# INITIAL CONTRASTS for marker selection WITHOUT taking correlated CT into account 
#[compare one group with average expression of all other groups]
annotation = factor(colnames(train2))
design <- model.matrix(~0+annotation)
colnames(design) <- unlist(lapply(strsplit(colnames(design),"annotation"), function(x) x[2]))
cont.matrix <- matrix((-1/ncol(design)),nrow=ncol(design),ncol=ncol(design))
colnames(cont.matrix) <- colnames(design)
diag(cont.matrix) <- (ncol(design)-1)/ncol(design)

v <- limma::voom(train2, design=design, plot=FALSE) 
fit <- limma::lmFit(v, design)
fit2 <- limma::contrasts.fit(fit, cont.matrix)
fit2 <- limma::eBayes(fit2, trend=TRUE)

markers = marker.fc(fit2, log2.threshold = log2(2))

#-------------------------------------------------------
### Generation of 1000 pseudo-bulk mixtures (T) (on test data)
cellType <- colnames(test)
colnames(test) <- original_cell_names_simu2

generator <- Generator(sce = test, phenoData = full_phenoData, Num.mixtures = 1000, pool.size = number_cells)
T <- generator[["T"]]
P <- generator[["P"]]

## Note, C is collected by train data, T is generated by test data

#-------------------------------------------------------
### Transformation, scaling/normalization, marker selection for bulk deconvolution methods and deconvolution:
if(deconv_type == "bulk"){

	T = Transformation(T, transformation)
	C = Transformation(C, transformation)

	T = Scaling(T, normalization)
	C = Scaling(C, normalization)

	# marker selection (on training data) 
	marker_distrib = marker_strategies(markers, marker_strategy, C)

	#If a cell type is removed, only meaningful mixtures where that CT was present (proportion < 0) are kept:
	if(to_remove != "none"){

		T <- T[,P[to_remove,] != 0]
		C <- C[, colnames(C) %in% rownames(P) & (!colnames(C) %in% to_remove)]
		P <- P[!rownames(P) %in% to_remove, colnames(T)]
		refProfiles.var = refProfiles.var[,colnames(refProfiles.var) %in% rownames(P) & (!colnames(refProfiles.var) %in% to_remove)]
		marker_distrib <- marker_distrib[marker_distrib$CT %in% rownames(P) & (marker_distrib$CT != to_remove),]
		
	}

	RESULTS = Deconvolution(T = T, C = C, method = method, P = P, elem = to_remove, marker_distrib = marker_distrib, refProfiles.var = refProfiles.var) 

} else if (deconv_type == "sc"){

	T = Transformation(T, transformation)
	C = Transformation(train_cellID, transformation)

	T = Scaling(T, normalization_scT)
	C = Scaling(C, normalization_scC)

	#If a cell type is removed, only meaningful mixtures where that CT was present (proportion < 0) are kept:
	if(to_remove != "none"){

		T <- T[,P[to_remove,] != 0]
		C <- C[,pDataC$cellType != to_remove]
		P <- P[!rownames(P) %in% to_remove, colnames(T)]
		pDataC <- pDataC[pDataC$cellType != to_remove,]

	}

	RESULTS = Deconvolution(T = T, C = C, method = method, phenoDataC = pDataC, P = P, elem = to_remove, refProfiles.var = refProfiles.var) 

}

RESULTS = RESULTS %>% dplyr::summarise(RMSE = sqrt(mean((observed_values-expected_values)^2)) %>% round(.,4), 
									   Pearson=cor(observed_values,expected_values) %>% round(.,4))
print(RESULTS)
# RMSE Pearson
# 0.0491  0.9785
