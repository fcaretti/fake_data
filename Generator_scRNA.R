############################ simulation 1 #########################################

library(SPARSim)

cellType1_dat <- readRDS("baron_scRNA_cellType1.rds")
cellType1_dat_norm <- scran_normalization(cellType1_dat)

cellType1_dat_colnames <- colnames(cellType1_dat) 
cellType1_dat_human1Index <- which(grepl("human1", cellType1_dat_colnames))
cellType1_dat_human2Index <- which(grepl("human2", cellType1_dat_colnames))
cellType1_dat_human3Index <- which(grepl("human3", cellType1_dat_colnames))
cellType1_dat_human4Index <- which(grepl("human4", cellType1_dat_colnames))

# Get column index for each experimental condition
cellType1_human1_column_index <- c(1:cellType1_dat_human1Index[[length(cellType1_dat_human1Index)]]) 
cellType1_human2_column_index <- c(cellType1_dat_human2Index[[1]]:cellType1_dat_human2Index[[length(cellType1_dat_human2Index)]]) 
cellType1_human3_column_index <- c(cellType1_dat_human3Index[[1]]:cellType1_dat_human3Index[[length(cellType1_dat_human3Index)]]) 
cellType1_human4_column_index <- c(cellType1_dat_human4Index[[1]]:cellType1_dat_human4Index[[length(cellType1_dat_human4Index)]]) 

# Create conditions param
cellType1_dat_conditions <- list(cond_A = cellType1_human1_column_index, 
                        cond_B = cellType1_human2_column_index, 
                        cond_C = cellType1_human3_column_index,
                        cond_D = cellType1_human4_column_index)

# Create SPARSim simulation parameter through the estimation from an existing count matrix
cellType1_dat_param <- SPARSim_estimate_parameter_from_data(raw_data = cellType1_dat, 
                                                   norm_data = cellType1_dat_norm, 
                                                   conditions = cellType1_dat_conditions)

# Run SPARSim simulation using the just created simulation parameter
cellType1_dat_simul <- SPARSim_simulation(dataset_parameter = cellType1_dat_param)

#Save, the cellType1_data_simul have five elemets: count_matrix, gene_matrix,
#abundanc_matrix, variability_matrix, batch_factors_matrix
saveRDS(cellType1_dat_simul, "baron_scRNA_cellType1_simu1.rds")

################# cell Type 2  ##########################################
cellType2_dat <- readRDS("baron_scRNA_cellType2.rds")

cellType2_dat_norm <- scran_normalization(cellType2_dat)

cellType2_dat_colnames <- colnames(cellType2_dat) 
cellType2_dat_human1Index <- which(grepl("human1", cellType2_dat_colnames))
cellType2_dat_human2Index <- which(grepl("human2", cellType2_dat_colnames))
cellType2_dat_human3Index <- which(grepl("human3", cellType2_dat_colnames))
cellType2_dat_human4Index <- which(grepl("human4", cellType2_dat_colnames))

# Get column index for each experimental condition
cellType2_human1_column_index <- c(1:cellType2_dat_human1Index[[length(cellType2_dat_human1Index)]]) 
cellType2_human2_column_index <- c(cellType2_dat_human2Index[[1]]:cellType2_dat_human2Index[[length(cellType2_dat_human2Index)]]) 
cellType2_human3_column_index <- c(cellType2_dat_human3Index[[1]]:cellType2_dat_human3Index[[length(cellType2_dat_human3Index)]]) 
cellType2_human4_column_index <- c(cellType2_dat_human4Index[[1]]:cellType2_dat_human4Index[[length(cellType2_dat_human4Index)]]) 

# Create conditions param
cellType2_dat_conditions <- list(cond_A = cellType2_human1_column_index, 
                                 cond_B = cellType2_human2_column_index, 
                                 cond_C = cellType2_human3_column_index,
                                 cond_D = cellType2_human4_column_index)

# Create SPARSim simulation parameter through the estimation from an existing count matrix
cellType2_dat_param <- SPARSim_estimate_parameter_from_data(raw_data = cellType2_dat, 
                                                            norm_data = cellType2_dat_norm, 
                                                            conditions = cellType2_dat_conditions)

# Run SPARSim simulation using the just created simulation parameter
cellType2_dat_simul <- SPARSim_simulation(dataset_parameter = cellType2_dat_param)

#Save
saveRDS(cellType2_dat_simul, "baron_scRNA_cellType2_simu1.rds")


############## cell Type 3  ##########################################

cellType3_dat <- readRDS("baron_scRNA_cellType3.rds")

cellType3_dat_norm <- scran_normalization(cellType3_dat)

cellType3_dat_colnames <- colnames(cellType3_dat) 
cellType3_dat_human1Index <- which(grepl("human1", cellType3_dat_colnames))
cellType3_dat_human2Index <- which(grepl("human2", cellType3_dat_colnames))
cellType3_dat_human3Index <- which(grepl("human3", cellType3_dat_colnames))
cellType3_dat_human4Index <- which(grepl("human4", cellType3_dat_colnames))

# Get column index for each experimental condition
cellType3_human1_column_index <- c(1:cellType3_dat_human1Index[[length(cellType3_dat_human1Index)]]) 
cellType3_human2_column_index <- c(cellType3_dat_human2Index[[1]]:cellType3_dat_human2Index[[length(cellType3_dat_human2Index)]]) 
cellType3_human3_column_index <- c(cellType3_dat_human3Index[[1]]:cellType3_dat_human3Index[[length(cellType3_dat_human3Index)]]) 
cellType3_human4_column_index <- c(cellType3_dat_human4Index[[1]]:cellType3_dat_human4Index[[length(cellType3_dat_human4Index)]]) 

# Create conditions param
cellType3_dat_conditions <- list(cond_A = cellType3_human1_column_index, 
                                 cond_B = cellType3_human2_column_index, 
                                 cond_C = cellType3_human3_column_index,
                                 cond_D = cellType3_human4_column_index)

# Create SPARSim simulation parameter through the estimation from an existing count matrix
cellType3_dat_param <- SPARSim_estimate_parameter_from_data(raw_data = cellType3_dat, 
                                                            norm_data = cellType3_dat_norm, 
                                                            conditions = cellType3_dat_conditions)

# Run SPARSim simulation using the just created simulation parameter
cellType3_dat_simul <- SPARSim_simulation(dataset_parameter = cellType3_dat_param)

#Save
saveRDS(cellType3_dat_simul, "baron_scRNA_cellType3_simu1.rds")



############## cell Type 4  ##########################################

cellType4_dat <- readRDS("baron_scRNA_cellType4.rds")

cellType4_dat_norm <- scran_normalization(cellType4_dat)

cellType4_dat_colnames <- colnames(cellType4_dat) 
cellType4_dat_human1Index <- which(grepl("human1", cellType4_dat_colnames))
cellType4_dat_human2Index <- which(grepl("human2", cellType4_dat_colnames))
cellType4_dat_human3Index <- which(grepl("human3", cellType4_dat_colnames))
cellType4_dat_human4Index <- which(grepl("human4", cellType4_dat_colnames))

# Get column index for each experimental condition
cellType4_human1_column_index <- c(1:cellType4_dat_human1Index[[length(cellType4_dat_human1Index)]]) 
cellType4_human2_column_index <- c(cellType4_dat_human2Index[[1]]:cellType4_dat_human2Index[[length(cellType4_dat_human2Index)]]) 
cellType4_human3_column_index <- c(cellType4_dat_human3Index[[1]]:cellType4_dat_human3Index[[length(cellType4_dat_human3Index)]]) 
cellType4_human4_column_index <- c(cellType4_dat_human4Index[[1]]:cellType4_dat_human4Index[[length(cellType4_dat_human4Index)]]) 

# Create conditions param
cellType4_dat_conditions <- list(cond_A = cellType4_human1_column_index, 
                                 cond_B = cellType4_human2_column_index, 
                                 cond_C = cellType4_human3_column_index,
                                 cond_D = cellType4_human4_column_index)

# Create SPARSim simulation parameter through the estimation from an existing count matrix
cellType4_dat_param <- SPARSim_estimate_parameter_from_data(raw_data = cellType4_dat, 
                                                            norm_data = cellType4_dat_norm, 
                                                            conditions = cellType4_dat_conditions)

# Run SPARSim simulation using the just created simulation parameter
cellType4_dat_simul <- SPARSim_simulation(dataset_parameter = cellType4_dat_param)

#Save
saveRDS(cellType4_dat_simul, "baron_scRNA_cellType4_simu1.rds")

#combine all the cell type together to make a new matrix
cellAllType_dat_simul <- cbind(cellType1_dat_simul$count_matrix, cellType2_dat_simul$count_matrix, cellType3_dat_simul$count_matrix, cellType4_dat_simul$count_matrix)
saveRDS(cellAllType_dat_simul, "baron_scRNA_cellAllType_simu1.rds")

#length(colnames(cellType1_dat_simul$count_matrix)) #2326
#1:2326 cell_type_1

#length(colnames(cellType2_dat_simul$count_matrix)) #2525 #2327+2525=4852
#2327:4851 cell_type_2

#length(colnames(cellType3_dat_simul$count_matrix)) #1077 #4852+1077=5929
#4852:5928 cell_type_3

#length(colnames(cellType4_dat_simul$count_matrix)) #601
#5929:length(colnames(cellAllType_dat_simul)) cell_type_4

# Create an empty character vector of length 20
pheno_vector2 <- character(length(colnames(cellAllType_dat_simul)))
pheno_vector2[1:2326] <- rep("cell_type_1", 2326)
pheno_vector2[2327:4851] <- rep("cell_type_2", 2525)
pheno_vector2[4852:5928] <- rep("cell_type_3", 1077)
pheno_vector2[5929:length(colnames(cellAllType_dat_simul))] <- rep("cell_type_4", 601)

pheno_vector1 <- colnames(cellAllType_dat_simul)
pheno_vector3 <- sub("_.*", "", pheno_vector1)

full_phenoData_simul <- data.frame(cellID = pheno_vector1, cellType = pheno_vector2, sampleID = pheno_vector3)

# Save the data frame as a tab-separated text file
write.table(full_phenoData_simul, "baron_phenoData_simu1.txt", sep = "\t", row.names = FALSE, col.names = FALSE)



############################ simulation 2 #############################################################


library(SPARSim)

cellType1_dat <- readRDS("baron_scRNA_cellType1.rds")
cellType1_dat_norm <- scran_normalization(cellType1_dat)

cellType1_dat_colnames <- colnames(cellType1_dat) 
cellType1_dat_human1Index <- which(grepl("human1", cellType1_dat_colnames))
cellType1_dat_human2Index <- which(grepl("human2", cellType1_dat_colnames))
cellType1_dat_human3Index <- which(grepl("human3", cellType1_dat_colnames))
cellType1_dat_human4Index <- which(grepl("human4", cellType1_dat_colnames))

# Get column index for each experimental condition
cellType1_human1_column_index <- c(1:cellType1_dat_human1Index[[length(cellType1_dat_human1Index)]]) 
cellType1_human2_column_index <- c(cellType1_dat_human2Index[[1]]:cellType1_dat_human2Index[[length(cellType1_dat_human2Index)]]) 
cellType1_human3_column_index <- c(cellType1_dat_human3Index[[1]]:cellType1_dat_human3Index[[length(cellType1_dat_human3Index)]]) 
cellType1_human4_column_index <- c(cellType1_dat_human4Index[[1]]:cellType1_dat_human4Index[[length(cellType1_dat_human4Index)]]) 

# Create conditions param
cellType1_dat_conditions <- list(cond_A = cellType1_human1_column_index, 
                                 cond_B = cellType1_human2_column_index, 
                                 cond_C = cellType1_human3_column_index,
                                 cond_D = cellType1_human4_column_index)

# Create SPARSim simulation parameter through the estimation from an existing count matrix
cellType1_dat_param <- SPARSim_estimate_parameter_from_data(raw_data = cellType1_dat, 
                                                            norm_data = cellType1_dat_norm, 
                                                            conditions = cellType1_dat_conditions)

# Run SPARSim simulation using the just created simulation parameter
cellType1_dat_simul <- SPARSim_simulation(dataset_parameter = cellType1_dat_param)

#Save, the cellType1_data_simul have five elemets: count_matrix, gene_matrix,
#abundanc_matrix, variability_matrix, batch_factors_matrix
saveRDS(cellType1_dat_simul, "baron_scRNA_cellType1_simu2.rds")

################# cell Type 2  ##########################################
cellType2_dat <- readRDS("baron_scRNA_cellType2.rds")

cellType2_dat_norm <- scran_normalization(cellType2_dat)

cellType2_dat_colnames <- colnames(cellType2_dat) 
cellType2_dat_human1Index <- which(grepl("human1", cellType2_dat_colnames))
cellType2_dat_human2Index <- which(grepl("human2", cellType2_dat_colnames))
cellType2_dat_human3Index <- which(grepl("human3", cellType2_dat_colnames))
cellType2_dat_human4Index <- which(grepl("human4", cellType2_dat_colnames))

# Get column index for each experimental condition
cellType2_human1_column_index <- c(1:cellType2_dat_human1Index[[length(cellType2_dat_human1Index)]]) 
cellType2_human2_column_index <- c(cellType2_dat_human2Index[[1]]:cellType2_dat_human2Index[[length(cellType2_dat_human2Index)]]) 
cellType2_human3_column_index <- c(cellType2_dat_human3Index[[1]]:cellType2_dat_human3Index[[length(cellType2_dat_human3Index)]]) 
cellType2_human4_column_index <- c(cellType2_dat_human4Index[[1]]:cellType2_dat_human4Index[[length(cellType2_dat_human4Index)]]) 

# Create conditions param
cellType2_dat_conditions <- list(cond_A = cellType2_human1_column_index, 
                                 cond_B = cellType2_human2_column_index, 
                                 cond_C = cellType2_human3_column_index,
                                 cond_D = cellType2_human4_column_index)

# Create SPARSim simulation parameter through the estimation from an existing count matrix
cellType2_dat_param <- SPARSim_estimate_parameter_from_data(raw_data = cellType2_dat, 
                                                            norm_data = cellType2_dat_norm, 
                                                            conditions = cellType2_dat_conditions)

# Run SPARSim simulation using the just created simulation parameter
cellType2_dat_simul <- SPARSim_simulation(dataset_parameter = cellType2_dat_param)

#Save
saveRDS(cellType2_dat_simul, "baron_scRNA_cellType2_simu2.rds")


############## cell Type 3  ##########################################

cellType3_dat <- readRDS("baron_scRNA_cellType3.rds")

cellType3_dat_norm <- scran_normalization(cellType3_dat)

cellType3_dat_colnames <- colnames(cellType3_dat) 
cellType3_dat_human1Index <- which(grepl("human1", cellType3_dat_colnames))
cellType3_dat_human2Index <- which(grepl("human2", cellType3_dat_colnames))
cellType3_dat_human3Index <- which(grepl("human3", cellType3_dat_colnames))
cellType3_dat_human4Index <- which(grepl("human4", cellType3_dat_colnames))

# Get column index for each experimental condition
cellType3_human1_column_index <- c(1:cellType3_dat_human1Index[[length(cellType3_dat_human1Index)]]) 
cellType3_human2_column_index <- c(cellType3_dat_human2Index[[1]]:cellType3_dat_human2Index[[length(cellType3_dat_human2Index)]]) 
cellType3_human3_column_index <- c(cellType3_dat_human3Index[[1]]:cellType3_dat_human3Index[[length(cellType3_dat_human3Index)]]) 
cellType3_human4_column_index <- c(cellType3_dat_human4Index[[1]]:cellType3_dat_human4Index[[length(cellType3_dat_human4Index)]]) 

# Create conditions param
cellType3_dat_conditions <- list(cond_A = cellType3_human1_column_index, 
                                 cond_B = cellType3_human2_column_index, 
                                 cond_C = cellType3_human3_column_index,
                                 cond_D = cellType3_human4_column_index)

# Create SPARSim simulation parameter through the estimation from an existing count matrix
cellType3_dat_param <- SPARSim_estimate_parameter_from_data(raw_data = cellType3_dat, 
                                                            norm_data = cellType3_dat_norm, 
                                                            conditions = cellType3_dat_conditions)

# Run SPARSim simulation using the just created simulation parameter
cellType3_dat_simul <- SPARSim_simulation(dataset_parameter = cellType3_dat_param)

#Save
saveRDS(cellType3_dat_simul, "baron_scRNA_cellType3_simu2.rds")



############## cell Type 4  ##########################################

cellType4_dat <- readRDS("baron_scRNA_cellType4.rds")

cellType4_dat_norm <- scran_normalization(cellType4_dat)

cellType4_dat_colnames <- colnames(cellType4_dat) 
cellType4_dat_human1Index <- which(grepl("human1", cellType4_dat_colnames))
cellType4_dat_human2Index <- which(grepl("human2", cellType4_dat_colnames))
cellType4_dat_human3Index <- which(grepl("human3", cellType4_dat_colnames))
cellType4_dat_human4Index <- which(grepl("human4", cellType4_dat_colnames))

# Get column index for each experimental condition
cellType4_human1_column_index <- c(1:cellType4_dat_human1Index[[length(cellType4_dat_human1Index)]]) 
cellType4_human2_column_index <- c(cellType4_dat_human2Index[[1]]:cellType4_dat_human2Index[[length(cellType4_dat_human2Index)]]) 
cellType4_human3_column_index <- c(cellType4_dat_human3Index[[1]]:cellType4_dat_human3Index[[length(cellType4_dat_human3Index)]]) 
cellType4_human4_column_index <- c(cellType4_dat_human4Index[[1]]:cellType4_dat_human4Index[[length(cellType4_dat_human4Index)]]) 

# Create conditions param
cellType4_dat_conditions <- list(cond_A = cellType4_human1_column_index, 
                                 cond_B = cellType4_human2_column_index, 
                                 cond_C = cellType4_human3_column_index,
                                 cond_D = cellType4_human4_column_index)

# Create SPARSim simulation parameter through the estimation from an existing count matrix
cellType4_dat_param <- SPARSim_estimate_parameter_from_data(raw_data = cellType4_dat, 
                                                            norm_data = cellType4_dat_norm, 
                                                            conditions = cellType4_dat_conditions)

# Run SPARSim simulation using the just created simulation parameter
cellType4_dat_simul <- SPARSim_simulation(dataset_parameter = cellType4_dat_param)

#Save
saveRDS(cellType4_dat_simul, "baron_scRNA_cellType4_simu2.rds")

#combine all the cell type together to make a new matrix
cellAllType_dat_simul <- cbind(cellType1_dat_simul$count_matrix, cellType2_dat_simul$count_matrix, cellType3_dat_simul$count_matrix, cellType4_dat_simul$count_matrix)
saveRDS(cellAllType_dat_simul, "baron_scRNA_cellAllType_simu2.rds")

#length(colnames(cellType1_dat_simul$count_matrix)) #2326
#1:2326 cell_type_1

#length(colnames(cellType2_dat_simul$count_matrix)) #2525 #2327+2525=4852
#2327:4851 cell_type_2

#length(colnames(cellType3_dat_simul$count_matrix)) #1077 #4852+1077=5929
#4852:5928 cell_type_3

#length(colnames(cellType4_dat_simul$count_matrix)) #601
#5929:length(colnames(cellAllType_dat_simul)) cell_type_4

# Create an empty character vector of length 20
pheno_vector2 <- character(length(colnames(cellAllType_dat_simul)))
pheno_vector2[1:2326] <- rep("cell_type_1", 2326)
pheno_vector2[2327:4851] <- rep("cell_type_2", 2525)
pheno_vector2[4852:5928] <- rep("cell_type_3", 1077)
pheno_vector2[5929:length(colnames(cellAllType_dat_simul))] <- rep("cell_type_4", 601)

pheno_vector1 <- colnames(cellAllType_dat_simul)
pheno_vector3 <- sub("_.*", "", pheno_vector1)

full_phenoData_simul <- data.frame(cellID = pheno_vector1, cellType = pheno_vector2, sampleID = pheno_vector3)

# Save the data frame as a tab-separated text file
write.table(full_phenoData_simul, "baron_phenoData_simu2.txt", sep = "\t", row.names = FALSE, col.names = FALSE)







#######################################################################################################

library(SPARSim)

grp = "PFC_L2_3"

data.dir.github = "./PFC_data/"

#The final 5 genes are MT-RNR2, MT-TM, MT-RNR1, MT-ND4L and MT-ATP8, dim(dat1): 18041 x 8626 
dat1 = readRDS(file.path(data.dir.github, sprintf("ct_mtx/%s.rds", grp)))

library(Matrix)

# Convert the sparse matrix to a dense matrix (regular matrix)
dat1 <- as.matrix(dat1)

# Convert the dense matrix to a data frame
dat1 <- as.data.frame(dat1)

dat1 = dat1[1:5000, 1:150]

dat1_norm <- scran_normalization(dat1)

# Get column index for each experimental condition
cond_A_column_index <- c(1:50) #c(1:8000) Condition A column indices: from column 1 to column 50
cond_B_column_index <- c(51:100) #c(8001:8600) Condition B column indices: from column 51 to column 100
cond_C_column_index <- c(101:150) #c(8601:8626) Condition C column indices: from column 101 to column 150

# Create conditions param
dat1_conditions <- list(cond_A = cond_A_column_index, 
                                        cond_B = cond_B_column_index, 
                                        cond_C = cond_C_column_index)

# Create SPARSim simulation parameter through the estimation from an existing count matrix
dat1_param <- SPARSim_estimate_parameter_from_data(raw_data = dat1, 
                                                          norm_data = dat1_norm, 
                                                          conditions = dat1_conditions)

# Run SPARSim simulation using the just created simulation parameter
dat1_result <- SPARSim_simulation(dataset_parameter = dat1_param)



dat2 = readRDS(file.path(data.dir.github, sprintf("ct_mtx/%s.rds", "PFC_L4")))








