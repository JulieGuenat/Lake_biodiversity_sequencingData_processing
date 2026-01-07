###########################################################################
################ Merging Euka 01 03 datasets ##############################
###############    Julie Guenat     #######################################
##########################################################################

# 1. Set directory ####
setwd("C:/Users/jguenat/Documents/PhD_Analyses/Cleaning_data/Output_metabar_data")

# 2. Load packages ####
library(tidyverse)

# 3. Load required functions ####

## 3.1. Merging euka01 and euka03 ####
#This function allows for merging the datasets from the same PCR plate (e.g. P1, P2, P3..)
# with different primers, here Euka01 and Euka03. 
#first, we determine the missing samples in one and the other datasets. 
#second, we add the missing samples in the datasets and add "0"
#third, we reorder the columns to ensure the merging is smooth
#fourth, we Combine unique taxonomic paths.
#fifth, we create an empty matrix that is filled with presence or absence of the different
#taxa according to dataset1 and dataset2. 

# Generalized function for aligning and merging datasets
align_and_merge_datasets <- function(dataset1, dataset2, samp, sample_start_col = 19) {
  # Extract sample columns
  sample_cols1 <- colnames(dataset1)[sample_start_col:ncol(dataset1)]
  sample_cols2 <- colnames(dataset2)[sample_start_col:ncol(dataset2)]
  
  # Identify missing samples in both datasets
  missing_in_dataset1 <- setdiff(samp, sample_cols1)
  missing_in_dataset2 <- setdiff(samp, sample_cols2)
  
  # Add missing columns with zeros to both datasets
  for (sample in missing_in_dataset1) {
    dataset1[[sample]] <- 0
  }
  for (sample in missing_in_dataset2) {
    dataset2[[sample]] <- 0
  }
  
  # Ensure columns in both datasets are ordered according to samp
  dataset1 <- dataset1[, c(colnames(dataset1)[1:(sample_start_col - 1)], samp)]
  dataset2 <- dataset2[, c(colnames(dataset2)[1:(sample_start_col - 1)], samp)]
  
  # Combine unique taxonomic paths
  taxa_X <- unique(c(dataset1$taxonomic_path, dataset2$taxonomic_path))
  
  # Create an empty matrix for presence/absence data
  combined_matrix <- matrix(0, nrow = length(taxa_X), ncol = length(samp))
  colnames(combined_matrix) <- samp
  rownames(combined_matrix) <- taxa_X
  
  # Populate the matrix with presence/absence data from both datasets
  for (j in seq_along(samp)) {
    sample <- samp[j]
    
    # Update matrix with dataset1
    taxa_present1 <- dataset1$taxonomic_path[dataset1[[sample]] > 0]
    combined_matrix[rownames(combined_matrix) %in% taxa_present1, j] <- 1
    
    # Update matrix with dataset2
    taxa_present2 <- dataset2$taxonomic_path[dataset2[[sample]] > 0]
    combined_matrix[rownames(combined_matrix) %in% taxa_present2, j] <- 1
  }
  
  # Return the combined matrix and taxonomic paths
  return(list(matrix = combined_matrix, taxa = taxa_X))
}

## 3.2 Extracting the taxonomical path from Obitools ####
#first we need to defined the ranks we would like to extract.
ranks <- c('clade', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species')

# Function to process each taxonomic path created by Fred ! 
process_taxonomic_path <- function(s) {
  # Split the string on "|" to get taxonomic elements
  taxonomic_elements <- strsplit(s, split = "|", fixed = TRUE)[[1]]
  
  # Further split on "@" to separate ranks and labels
  element_splits <- lapply(taxonomic_elements, function(t) strsplit(t, split = "@", fixed = TRUE)[[1]])
  
  # Filter out elements that match the desired ranks
  filtered_elements <- Filter(function(v) v[3] %in% ranks, element_splits)
  
  # If no matching ranks are found, return an empty dataframe with taxid and taxo_path
  if (length(filtered_elements) == 0) {
    taxid <- element_splits[[length(element_splits)]][1]
    return(data.frame(taxid = taxid, taxo_path = s, matrix(NA, nrow = 1, ncol = length(ranks), dimnames = list(NULL, ranks))))
  }
  
  # Create a named vector to store ranks
  rank_vector <- setNames(rep(NA, length(ranks)), ranks)
  
  # Populate rank_vector with actual values from filtered_elements
  for (element in filtered_elements) {
    rank_vector[element[3]] <- element[2]
  }
  
  # Extract the taxid (ID before "@")
  taxid <- element_splits[[length(element_splits)]][1]
  
  # Return a dataframe with one row, including taxo_path
  result <- data.frame(t(c(taxid, rank_vector)), stringsAsFactors = FALSE)
  colnames(result) <- c("taxid", ranks)
  
  # Add the original taxo_path as a column
  result$taxo_path <- s
  
  return(result)
}

####### P1 ######################################
## P1 DNA Euka01 ####
P1.DNA.euka01.reads <- read.table("P1_euka01_DNA/files/P1_euka01_DNA_reads_mean_final.csv", header = T, sep = ",")
P1.DNA.euka01.motus <- read.table("P1_euka01_DNA/files/P1_euka01_DNA_motus_final.csv", header = T, sep = ",")

rownames(P1.DNA.euka01.reads) <- P1.DNA.euka01.reads$X
P1.DNA.euka01.reads = subset(P1.DNA.euka01.reads, select = -c(1))
P1.DNA.euka01.reads <- as.data.frame(t(P1.DNA.euka01.reads))
P1.DNA.euka01.reads$id <- row.names(P1.DNA.euka01.reads)
P1.DNA.euka01.reads$id[1] <- gsub(".*(MN00722.*)", "\\1", P1.DNA.euka01.reads$id[1])

P1.DNA.euka01.motus <- P1.DNA.euka01.motus %>%
  mutate(id = X %>% 
           gsub(":", ".", .) %>% 
           gsub("-", ".", .) %>% 
           gsub("\\[", ".", .) %>% 
           gsub("\\]", ".", .))
P1.DNA.euka01.motus$id[1] <- gsub(".*(MN00722.*)", "\\1", P1.DNA.euka01.motus$id[1])
P1_euka01_DNA_reads_motus = left_join(P1.DNA.euka01.motus, P1.DNA.euka01.reads, by="id")

rm(P1.DNA.euka01.reads, P1.DNA.euka01.motus)


## P1 DNA euka03 ####
P1.DNA.euka03.reads <- read.table("P1_euka03_DNA/files/P1_euka03_DNA_reads_mean_final.csv", header = T, sep = ",")
P1.DNA.euka03.motus <- read.table("P1_euka03_DNA/files/P1_euka03_DNA_motus_final.csv", header = T, sep = ",")

rownames(P1.DNA.euka03.reads) <- P1.DNA.euka03.reads$X
P1.DNA.euka03.reads = subset(P1.DNA.euka03.reads, select = -c(1))
P1.DNA.euka03.reads <- as.data.frame(t(P1.DNA.euka03.reads))
P1.DNA.euka03.reads$id <- row.names(P1.DNA.euka03.reads)
P1.DNA.euka03.reads$id[1] <- gsub(".*(MN00722.*)", "\\1", P1.DNA.euka03.reads$id[1])

P1.DNA.euka03.motus <- P1.DNA.euka03.motus %>%
  mutate(id = X %>% 
           gsub(":", ".", .) %>% 
           gsub("-", ".", .) %>% 
           gsub("\\[", ".", .) %>% 
           gsub("\\]", ".", .))
P1.DNA.euka03.motus$id[1] <- gsub(".*(MN00722.*)", "\\1", P1.DNA.euka03.motus$id[1])
P1_euka03_DNA_reads_motus = left_join(P1.DNA.euka03.motus, P1.DNA.euka03.reads, by="id")

rm(P1.DNA.euka03.reads, P1.DNA.euka03.motus)

## Merging Euka01 and euka03 ####
#First, we need to define the list of the samples presen in the pcr plate
samp<-as.character(1:58)

#then we run the command to merge both datasets
P1_merge <- align_and_merge_datasets(P1_euka01_DNA_reads_motus, P1_euka03_DNA_reads_motus, samp)

# Extract the combined presence/absence matrix and taxonomic paths
P1_euka <- P1_merge$matrix
taxa_X <- P1_merge$taxa

# Verify all samples are included
all_samples_included <- setequal(colnames(P1_euka), samp)
print(all_samples_included)  # Should print TRUE

# Adding Taxonomic Ranks
P1 <- as.data.frame(P1_euka)
P1$taxo_path <- rownames(P1)

# Parse taxonomic paths using the defined function
dff <- do.call(rbind, lapply(P1$taxo_path, process_taxonomic_path))
P1_euka_final <- merge(P1, dff, by = "taxo_path", all.x = TRUE)

# Removing columns corresponding to negative controls in the field == column for which
# colsum=0
# In this plate ctrl should be: "13", "24", "34", "57"
P1_euka_final<-P1_euka_final %>%
  mutate(across(2:59, as.numeric)) %>%
  select(c(1, where(~ is.numeric(.) && sum(.) > 0), 60:68))

# Verification of merging
final_samples <- colnames(P1_euka_final)[2:(length(samp) + 1)]  # Adjust for column indexing
missing_samples <- setdiff(samp, final_samples)
print(missing_samples)# Output missing samples (should be neg ctrl)

rm(taxa_X, samp, missing_samples, final_samples, all_samples_included, P1_merge, 
   P1_euka01_DNA_reads_motus, P1_euka03_DNA_reads_motus, P1_euka, P1, dff)


####### P2 ######################################
## P2 DNA Euka01 ####
P2.DNA.euka01.reads <- read.table("P2_euka01_DNA/files/P2_euka01_DNA_reads_mean_final.csv", header = T, sep = ",")
P2.DNA.euka01.motus <- read.table("P2_euka01_DNA/files/P2_euka01_DNA_motus_final.csv", header = T, sep = ",")

rownames(P2.DNA.euka01.reads) <- P2.DNA.euka01.reads$X
P2.DNA.euka01.reads = subset(P2.DNA.euka01.reads, select = -c(1))
P2.DNA.euka01.reads <- as.data.frame(t(P2.DNA.euka01.reads))
P2.DNA.euka01.reads$id <- row.names(P2.DNA.euka01.reads)
P2.DNA.euka01.reads$id[1] <- gsub(".*(MN00722.*)", "\\1", P2.DNA.euka01.reads$id[1])

P2.DNA.euka01.motus <- P2.DNA.euka01.motus %>%
  mutate(id = X %>% 
           gsub(":", ".", .) %>% 
           gsub("-", ".", .) %>% 
           gsub("\\[", ".", .) %>% 
           gsub("\\]", ".", .))
P2.DNA.euka01.motus$id[1] <- gsub(".*(MN00722.*)", "\\1", P2.DNA.euka01.motus$id[1])
P2_euka01_DNA_reads_motus = left_join(P2.DNA.euka01.motus, P2.DNA.euka01.reads, by="id")

rm(P2.DNA.euka01.reads, P2.DNA.euka01.motus)


## P2 DNA euka03 ####
P2.DNA.euka03.reads <- read.table("P2_euka03_DNA/files/P2_euka03_DNA_reads_mean_final.csv", header = T, sep = ",")
P2.DNA.euka03.motus <- read.table("P2_euka03_DNA/files/P2_euka03_DNA_motus_final.csv", header = T, sep = ",")

rownames(P2.DNA.euka03.reads) <- P2.DNA.euka03.reads$X
P2.DNA.euka03.reads = subset(P2.DNA.euka03.reads, select = -c(1))
P2.DNA.euka03.reads <- as.data.frame(t(P2.DNA.euka03.reads))
P2.DNA.euka03.reads$id <- row.names(P2.DNA.euka03.reads)
P2.DNA.euka03.reads$id[1] <- gsub(".*(MN00722.*)", "\\1", P2.DNA.euka03.reads$id[1])

P2.DNA.euka03.motus <- P2.DNA.euka03.motus %>%
  mutate(id = X %>% 
           gsub(":", ".", .) %>% 
           gsub("-", ".", .) %>% 
           gsub("\\[", ".", .) %>% 
           gsub("\\]", ".", .))
P2.DNA.euka03.motus$id[1] <- gsub(".*(MN00722.*)", "\\1", P2.DNA.euka03.motus$id[1])
P2_euka03_DNA_reads_motus = left_join(P2.DNA.euka03.motus, P2.DNA.euka03.reads, by="id")

rm(P2.DNA.euka03.reads, P2.DNA.euka03.motus)

## Merging Euka01 and euka03 ####
#First, we need to define the list of the samples presen in the pcr plate
samp<-as.character(59:116)

#then we run the command to merge both datasets
P2_merge <- align_and_merge_datasets(P2_euka01_DNA_reads_motus, P2_euka03_DNA_reads_motus, samp)

# Extract the combined presence/absence matrix and taxonomic paths
P2_euka <- P2_merge$matrix
taxa_X <- P2_merge$taxa

# Verify all samples are included
all_samples_included <- setequal(colnames(P2_euka), samp)
print(all_samples_included)  # Should print TRUE

# Adding Taxonomic Ranks
P2 <- as.data.frame(P2_euka)
P2$taxo_path <- rownames(P2)

# Parse taxonomic paths using the defined function
dff <- do.call(rbind, lapply(P2$taxo_path, process_taxonomic_path))
P2_euka_final <- merge(P2, dff, by = "taxo_path", all.x = TRUE)

# Removing columns corresponding to negative controls in the field == column for which
# colsum=0
# In this plate ctrl should be: "68", "89", "95", "116"
P2_euka_final<-P2_euka_final %>%
  mutate(across(2:59, as.numeric)) %>%
  select(c(1, where(~ is.numeric(.) && sum(.) > 0), 60:68))

# Verification of merging
final_samples <- colnames(P2_euka_final)[2:(length(samp) + 1)]  # Adjust for column indexing
missing_samples <- setdiff(samp, final_samples)
print(missing_samples)# Output missing samples (should be neg ctrl)

rm(taxa_X, samp, missing_samples, final_samples, all_samples_included, P2_merge, 
   P2_euka01_DNA_reads_motus, P2_euka03_DNA_reads_motus, P2_euka, P2, dff)

####### P3 ######################################
## P3 DNA Euka01 ####
P3.DNA.euka01.reads <- read.table("P3_euka01_DNA/files/P3_euka01_DNA_reads_mean_final.csv", header = T, sep = ",")
P3.DNA.euka01.motus <- read.table("P3_euka01_DNA/files/P3_euka01_DNA_motus_final.csv", header = T, sep = ",")

rownames(P3.DNA.euka01.reads) <- P3.DNA.euka01.reads$X
P3.DNA.euka01.reads = subset(P3.DNA.euka01.reads, select = -c(1))
P3.DNA.euka01.reads <- as.data.frame(t(P3.DNA.euka01.reads))
P3.DNA.euka01.reads$id <- row.names(P3.DNA.euka01.reads)
P3.DNA.euka01.reads$id[1] <- gsub(".*(MN00722.*)", "\\1", P3.DNA.euka01.reads$id[1])

P3.DNA.euka01.motus <- P3.DNA.euka01.motus %>%
  mutate(id = X %>% 
           gsub(":", ".", .) %>% 
           gsub("-", ".", .) %>% 
           gsub("\\[", ".", .) %>% 
           gsub("\\]", ".", .))
P3.DNA.euka01.motus$id[1] <- gsub(".*(MN00722.*)", "\\1", P3.DNA.euka01.motus$id[1])
P3_euka01_DNA_reads_motus = left_join(P3.DNA.euka01.motus, P3.DNA.euka01.reads, by="id")

rm(P3.DNA.euka01.reads, P3.DNA.euka01.motus)


## P3 DNA euka03 ####
P3.DNA.euka03.reads <- read.table("P3_euka03_DNA/files/P3_euka03_DNA_reads_mean_final.csv", header = T, sep = ",")
P3.DNA.euka03.motus <- read.table("P3_euka03_DNA/files/P3_euka03_DNA_motus_final.csv", header = T, sep = ",")

rownames(P3.DNA.euka03.reads) <- P3.DNA.euka03.reads$X
P3.DNA.euka03.reads = subset(P3.DNA.euka03.reads, select = -c(1))
P3.DNA.euka03.reads <- as.data.frame(t(P3.DNA.euka03.reads))
P3.DNA.euka03.reads$id <- row.names(P3.DNA.euka03.reads)
P3.DNA.euka03.reads$id[1] <- gsub(".*(MN00722.*)", "\\1", P3.DNA.euka03.reads$id[1])

P3.DNA.euka03.motus <- P3.DNA.euka03.motus %>%
  mutate(id = X %>% 
           gsub(":", ".", .) %>% 
           gsub("-", ".", .) %>% 
           gsub("\\[", ".", .) %>% 
           gsub("\\]", ".", .))
P3.DNA.euka03.motus$id[1] <- gsub(".*(MN00722.*)", "\\1", P3.DNA.euka03.motus$id[1])
P3_euka03_DNA_reads_motus = left_join(P3.DNA.euka03.motus, P3.DNA.euka03.reads, by="id")

rm(P3.DNA.euka03.reads, P3.DNA.euka03.motus)

## Merging Euka01 and euka03 ####
#First, we need to define the list of the samples presen in the pcr plate
samp<-as.character(117:174)

#then we run the command to merge both datasets
P3_merge <- align_and_merge_datasets(P3_euka01_DNA_reads_motus, P3_euka03_DNA_reads_motus, samp)

# Extract the combined presence/absence matrix and taxonomic paths
P3_euka <- P3_merge$matrix
taxa_X <- P3_merge$taxa

# Verify all samples are included
all_samples_included <- setequal(colnames(P3_euka), samp)
print(all_samples_included)  # Should print TRUE

# Adding Taxonomic Ranks
P3 <- as.data.frame(P3_euka)
P3$taxo_path <- rownames(P3)

# Parse taxonomic paths using the defined function
dff <- do.call(rbind, lapply(P3$taxo_path, process_taxonomic_path))
P3_euka_final <- merge(P3, dff, by = "taxo_path", all.x = TRUE)

# Removing columns corresponding to negative controls in the field == column for which
# colsum=0
# In this plate ctrl should be: "135", "140", "160", "172"
P3_euka_final<-P3_euka_final %>%
  mutate(across(2:59, as.numeric)) %>%
  select(c(1, where(~ is.numeric(.) && sum(.) > 0), 60:68))

# Verification of merging
final_samples <- colnames(P3_euka_final)[2:(length(samp) + 1)]  # Adjust for column indexing
missing_samples <- setdiff(samp, final_samples)
print(missing_samples)# Output missing samples (should be neg ctrl)
#In the case of this PCR plates some samples are not working at all
#"122" "139" "169"; "135" "140" "160" "172"

rm(taxa_X, samp, missing_samples, final_samples, all_samples_included, P3_merge, 
   P3_euka01_DNA_reads_motus, P3_euka03_DNA_reads_motus, P3_euka, P3, dff)

####### P4 ######################################
## P4 DNA Euka01 ####
P4.DNA.euka01.reads <- read.table("P4_euka01_DNA/files/P4_euka01_DNA_reads_mean_final.csv", header = T, sep = ",")
P4.DNA.euka01.motus <- read.table("P4_euka01_DNA/files/P4_euka01_DNA_motus_final.csv", header = T, sep = ",")

rownames(P4.DNA.euka01.reads) <- P4.DNA.euka01.reads$X
P4.DNA.euka01.reads = subset(P4.DNA.euka01.reads, select = -c(1))
P4.DNA.euka01.reads <- as.data.frame(t(P4.DNA.euka01.reads))
P4.DNA.euka01.reads$id <- row.names(P4.DNA.euka01.reads)
P4.DNA.euka01.reads$id[1] <- gsub(".*(MN00722.*)", "\\1", P4.DNA.euka01.reads$id[1])

P4.DNA.euka01.motus <- P4.DNA.euka01.motus %>%
  mutate(id = X %>% 
           gsub(":", ".", .) %>% 
           gsub("-", ".", .) %>% 
           gsub("\\[", ".", .) %>% 
           gsub("\\]", ".", .))
P4.DNA.euka01.motus$id[1] <- gsub(".*(MN00722.*)", "\\1", P4.DNA.euka01.motus$id[1])
P4_euka01_DNA_reads_motus = left_join(P4.DNA.euka01.motus, P4.DNA.euka01.reads, by="id")

rm(P4.DNA.euka01.reads, P4.DNA.euka01.motus)


## P4 DNA euka03 ####
P4.DNA.euka03.reads <- read.table("P4_euka03_DNA/files/P4_euka03_DNA_reads_mean_final.csv", header = T, sep = ",")
P4.DNA.euka03.motus <- read.table("P4_euka03_DNA/files/P4_euka03_DNA_motus_final.csv", header = T, sep = ",")

rownames(P4.DNA.euka03.reads) <- P4.DNA.euka03.reads$X
P4.DNA.euka03.reads = subset(P4.DNA.euka03.reads, select = -c(1))
P4.DNA.euka03.reads <- as.data.frame(t(P4.DNA.euka03.reads))
P4.DNA.euka03.reads$id <- row.names(P4.DNA.euka03.reads)
P4.DNA.euka03.reads$id[1] <- gsub(".*(MN00722.*)", "\\1", P4.DNA.euka03.reads$id[1])

P4.DNA.euka03.motus <- P4.DNA.euka03.motus %>%
  mutate(id = X %>% 
           gsub(":", ".", .) %>% 
           gsub("-", ".", .) %>% 
           gsub("\\[", ".", .) %>% 
           gsub("\\]", ".", .))
P4.DNA.euka03.motus$id[1] <- gsub(".*(MN00722.*)", "\\1", P4.DNA.euka03.motus$id[1])
P4_euka03_DNA_reads_motus = left_join(P4.DNA.euka03.motus, P4.DNA.euka03.reads, by="id")

rm(P4.DNA.euka03.reads, P4.DNA.euka03.motus)


## Merging Euka01 and euka03 ####
#First, we need to define the list of the samples presen in the pcr plate
samp<-as.character(175:232)

#then we run the command to merge both datasets
P4_merge <- align_and_merge_datasets(P4_euka01_DNA_reads_motus, P4_euka03_DNA_reads_motus, samp)

# Extract the combined presence/absence matrix and taxonomic paths
P4_euka <- P4_merge$matrix
taxa_X <- P4_merge$taxa

# Verify all samples are included
all_samples_included <- setequal(colnames(P4_euka), samp)
print(all_samples_included)  # Should print TRUE

# Adding Taxonomic Ranks
P4 <- as.data.frame(P4_euka)
P4$taxo_path <- rownames(P4)

# Parse taxonomic paths using the defined function
dff <- do.call(rbind, lapply(P4$taxo_path, process_taxonomic_path))
P4_euka_final <- merge(P4, dff, by = "taxo_path", all.x = TRUE)

# Removing columns corresponding to negative controls in the field == column for which
# colsum=0
# In this plate ctrl should be: "187", "214"
P4_euka_final<-P4_euka_final %>%
  mutate(across(2:59, as.numeric)) %>%
  select(c(1, where(~ is.numeric(.) && sum(.) > 0), 60:68))

# Verification of merging
final_samples <- colnames(P4_euka_final)[2:(length(samp) + 1)]  # Adjust for column indexing
missing_samples <- setdiff(samp, final_samples)
print(missing_samples)# Output missing samples (should be neg ctrl)
#In the case of this PCR plates some samples are not working at all
#"182" "187" "214"

rm(taxa_X, samp, missing_samples, final_samples, all_samples_included, P4_merge, 
   P4_euka01_DNA_reads_motus, P4_euka03_DNA_reads_motus, P4_euka, P4, dff)

####### P5 ######################################
## P5 DNA Euka01 ####
P5.DNA.euka01.reads <- read.table("P5_euka01_DNA/files/P5_euka01_DNA_reads_mean_final.csv", header = T, sep = ",")
P5.DNA.euka01.motus <- read.table("P5_euka01_DNA/files/P5_euka01_DNA_motus_final.csv", header = T, sep = ",")

rownames(P5.DNA.euka01.reads) <- P5.DNA.euka01.reads$X
P5.DNA.euka01.reads = subset(P5.DNA.euka01.reads, select = -c(1))
P5.DNA.euka01.reads <- as.data.frame(t(P5.DNA.euka01.reads))
P5.DNA.euka01.reads$id <- row.names(P5.DNA.euka01.reads)
P5.DNA.euka01.reads$id[1] <- gsub(".*(MN00722.*)", "\\1", P5.DNA.euka01.reads$id[1])

P5.DNA.euka01.motus <- P5.DNA.euka01.motus %>%
  mutate(id = X %>% 
           gsub(":", ".", .) %>% 
           gsub("-", ".", .) %>% 
           gsub("\\[", ".", .) %>% 
           gsub("\\]", ".", .))
P5.DNA.euka01.motus$id[1] <- gsub(".*(MN00722.*)", "\\1", P5.DNA.euka01.motus$id[1])
P5_euka01_DNA_reads_motus = left_join(P5.DNA.euka01.motus, P5.DNA.euka01.reads, by="id")

rm(P5.DNA.euka01.reads, P5.DNA.euka01.motus)


## P5 DNA euka03 ####
P5.DNA.euka03.reads <- read.table("P5_euka03_DNA/files/P5_euka03_DNA_reads_mean_final.csv", header = T, sep = ",")
P5.DNA.euka03.motus <- read.table("P5_euka03_DNA/files/P5_euka03_DNA_motus_final.csv", header = T, sep = ",")

rownames(P5.DNA.euka03.reads) <- P5.DNA.euka03.reads$X
P5.DNA.euka03.reads = subset(P5.DNA.euka03.reads, select = -c(1))
P5.DNA.euka03.reads <- as.data.frame(t(P5.DNA.euka03.reads))
P5.DNA.euka03.reads$id <- row.names(P5.DNA.euka03.reads)
P5.DNA.euka03.reads$id[1] <- gsub(".*(MN00722.*)", "\\1", P5.DNA.euka03.reads$id[1])

P5.DNA.euka03.motus <- P5.DNA.euka03.motus %>%
  mutate(id = X %>% 
           gsub(":", ".", .) %>% 
           gsub("-", ".", .) %>% 
           gsub("\\[", ".", .) %>% 
           gsub("\\]", ".", .))
P5.DNA.euka03.motus$id[1] <- gsub(".*(MN00722.*)", "\\1", P5.DNA.euka03.motus$id[1])
P5_euka03_DNA_reads_motus = left_join(P5.DNA.euka03.motus, P5.DNA.euka03.reads, by="id")

rm(P5.DNA.euka03.reads, P5.DNA.euka03.motus)

## Merging Euka01 and euka03 ####
#First, we need to define the list of the samples presen in the pcr plate
samp<-as.character(233:290)

#then we run the command to merge both datasets
P5_merge <- align_and_merge_datasets(P5_euka01_DNA_reads_motus, P5_euka03_DNA_reads_motus, samp)

# Extract the combined presence/absence matrix and taxonomic paths
P5_euka <- P5_merge$matrix
taxa_X <- P5_merge$taxa

# Verify all samples are included
all_samples_included <- setequal(colnames(P5_euka), samp)
print(all_samples_included)  # Should print TRUE

# Adding Taxonomic Ranks
P5 <- as.data.frame(P5_euka)
P5$taxo_path <- rownames(P5)

# Parse taxonomic paths using the defined function
dff <- do.call(rbind, lapply(P5$taxo_path, process_taxonomic_path))
P5_euka_final <- merge(P5, dff, by = "taxo_path", all.x = TRUE)

# Removing columns corresponding to negative controls in the field == column for which
# colsum=0
# In this plate ctrl should be: "233", "252", "271", "290"
P5_euka_final<-P5_euka_final %>%
  mutate(across(2:59, as.numeric)) %>%
  select(c(1, where(~ is.numeric(.) && sum(.) > 0), 60:68))

# Verification of merging
final_samples <- colnames(P5_euka_final)[2:(length(samp) + 1)]  # Adjust for column indexing
missing_samples <- setdiff(samp, final_samples)
print(missing_samples)# Output missing samples (should be neg ctrl)


rm(taxa_X, samp, missing_samples, final_samples, all_samples_included, P5_merge, 
   P5_euka01_DNA_reads_motus, P5_euka03_DNA_reads_motus, P5_euka, P5, dff)

###### MERGING PRES/ABS P1:P5 ################################
# Pres abs combination ####
# List of datasets
datasets <- list(P1_euka_final, P2_euka_final, P3_euka_final, P4_euka_final, P5_euka_final)

# Function to merge two datasets, keeping all rows and handling taxonomic columns to avoid .x and .y suffixes
merge_with_zeros_and_combine <- function(x, y) {
  # Perform the merge, keeping all rows and columns
  merged <- merge(x, y, by = "taxo_path", all = TRUE)
  
  # Combine the taxonomic columns, prioritizing non-NA values
  for (rank in c("taxid", "clade", "kingdom", "phylum", "class", "order", "family", "genus", "species")) {
    merged[[rank]] <- ifelse(is.na(merged[[paste0(rank, ".x")]]), 
                             merged[[paste0(rank, ".y")]], 
                             merged[[paste0(rank, ".x")]])
  }
  
  # Drop the duplicated columns (.x and .y)
  merged <- merged[, !grepl("\\.x$|\\.y$", names(merged))]
  
  # Identify sample columns (non-taxonomic)
  taxo_columns <- c("taxo_path", "taxid", "clade", "kingdom", "phylum", "class", "order", "family", "genus", "species")
  sample_columns <- setdiff(names(merged), taxo_columns)
  
  # Replace NA in sample columns with 0
  merged[sample_columns][is.na(merged[sample_columns])] <- 0
  
  return(merged)
}

# Use Reduce to iteratively merge all datasets
merged_df <- Reduce(merge_with_zeros_and_combine, datasets)

###### MERGING FIELD REPLICATES ##############################
# Merge REPLICATES ####

# Load the dataset
pa <- merged_df
samples_pla<-read.csv("../sample_names.csv", sep= ";", header=T)

# Define a function to simplify the sample names for merging
simplify_name <- function(sample_name) {
  # Split the name by '_' and remove the middle part
  parts <- unlist(strsplit(sample_name, "_"))
  paste(parts[1], parts[length(parts)], sep="_")
}

# Apply the simplify_name function to create simplified names
simplified_names <- sapply(samples_pla$Samples_names, simplify_name)

# Create a mapping of sample IDs to simplified names
sample_mapping <- data.frame(Sample_ID = samples_pla$X, Simplified_Name = simplified_names)

# Initialize a new data frame to store results with simplified names
pa2 <- pa

# Loop through the simplified names and merge the replicates
for (simplified_name in unique(sample_mapping$Simplified_Name)) {
  
  # Get the corresponding sample IDs
  sample_ids <- sample_mapping$Sample_ID[sample_mapping$Simplified_Name == simplified_name]
  
  # Ensure these IDs exist in the dataset columns
  replicate_cols <- as.character(sample_ids[sample_ids %in% colnames(pa)])
  
  if (length(replicate_cols) > 1) {  # Proceed only if more than one replicate exists
    
    # Combine the replicate columns
    pa2[[simplified_name]] <- apply(pa[, replicate_cols], 1, function(x) {
      if (sum(x > 0) >= 1) 1 else 0
    })
    
  } else if (length(replicate_cols) == 1) {
    # If only one column exists, just rename it in the new dataset
    pa2[[simplified_name]] <- pa[[replicate_cols]]
  }
}

pa_final <- pa2[,c(1, 271:464)]

write.csv(pa_final, "../final_datasets/pres_abs_mergedrep_euka_DNA.csv", row.names = F)
 rm(datasets, merged_df, P1_euka_final, P2_euka_final, P3_euka_final, P4_euka_final,
    P5_euka_final, pa, pa2, sample_mapping, samples_pla)

###### SPCEIES RICHNESS DATASET ##############################
# Species richness computation ####
#data loading 
pa <- read.csv("../final_datasets/pres_abs_mergedrep_euka_DNA.csv", header = T, sep = ",")

pa_rich<- pa[,c(11:195)]

sp_rich <- as.data.frame(colSums(pa_rich))
colnames(sp_rich) <- "Sp_rich"
sp_rich$sample_id <- row.names(sp_rich) 

# Merge simplified samples names with sp. richness ####
samp<-read.csv("../final_datasets/sample_mapping_names.csv", header = T, sep = ",")

sampID_sprich <- merge(samp, sp_rich, by.x = "Simplified_Name", by.y = "sample_id")%>%
  group_by(Simplified_Name) %>%
  slice(1) 

# merge metadat with sp_richness ####
samples_pla<-read.csv("../sample_names.csv", sep= ";", header=T)

sp_rich_pla <- merge(sampID_sprich, samples_pla[, c("X", "Samples_names")], 
                     by.x = "Sample_ID", by.y = "X", 
                     all.x = TRUE)

#metadata
metadata<- read.csv("../Samples_info_Planaqua.csv", header= T, sep=";", dec= ",")

Sp_rich_metadata<-merge(metadata, sp_rich_pla, by.x = "Samples_ID", by.y = "Samples_names")

write.csv(Sp_rich_metadata, "../final_datasets/Species_richness_metadata_euka_DNA.csv", row.names = F)
