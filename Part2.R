#Part 2 ----
# This script is to extract the counts for solarnum pennelli


# Here we load the libraries ----
library(tidyverse)
library(tximport)
library(datapasta)
library(ape)
library(rtracklayer)

# We need the paths ----
path2 <- file.path("/home/adejoro/Desktop/introgression_lines/SRA_PROJECTS/Spenn/SraRunTable", study_design$Run, "abundance.tsv")
path2

all(file.exists(path2))

# We proceed to read in the annotation file for solarnum pennelli

pennelli_annotation_file <- read.gff("spenn_v2.0_gene_models_annot.gff")

pennelli_annotation_file

# Here, work is done on a different script: rough_work.R to extract only the target_id and gene_name from the annotation file

pennelli_annotation_selected

#Importing reads and annotation file into iximport

sol_pennelli_genes <- tximport(path2,
                               type = "kallisto",
                               tx2gene = pennelli_annotation_selected,
                               txOut = FALSE,
                               countsFromAbundance = "lengthScaledTPM",
                               ignoreTxVersion = FALSE)
sol_pennelli_counts_matrix <- sol_pennelli_genes$counts


#Let's generate DGEList for sol_pennelli from the counts
sol_pennelli_DGEList <- DGEList(sol_pennelli_genes$counts)


#Let's generate a cpm object from this DGEList
sol_pennelli_cpm <- cpm(sol_pennelli_DGEList)
sol_pennelli_cpm

#We can now use this cpm object to create a filter condition by which we'll subset the DGEList we created the cpm from
sol_pennelli_filter_condition <- rowSums(sol_pennelli_cpm > 1)>=300

#We can now use our filter condition to subset our DGEList rows
sol_pennelli_DGEList_subset <- sol_pennelli_DGEList[sol_pennelli_filter_condition,]

#We now have the DGEList subset, but now we must proceed with normalizing this DGEList subset
sol_pennelli_DGEList_subset_normalized <- calcNormFactors(sol_pennelli_DGEList_subset, method = "TMM") 

#With that, we have a normlaized subset of the original DGEList data from which to generate new cpms
sol_pennelli_cpms_norm.filtered.logT <- cpm(sol_pennelli_DGEList_subset_normalized, log = TRUE)

#Now that we have this cpm, we can jump to carrying out PCA
sol_pennelli_PCA <- prcomp(t(sol_pennelli_cpms_norm.filtered.logT), scale. = F, retx = T)
ls(sol_pennelli_PCA)

screeplot(sol_pennelli_PCA)

sol_pennelli_PCA_variance <- sol_pennelli_PCA$sdev^2

sol_pennelli_PCA_variance_perPC <- round(sol_pennelli_PCA_variance / sum(sol_pennelli_PCA_variance) * 100, 1)

sol_pennelli_PCA_of_samples <- as_tibble(sol_pennelli_PCA$x)

ggplot(sol_pennelli_PCA_of_samples) +
  aes(x = PC1, y = PC2, color = group) +
  #geom_text(aes(label = myLabel)) +
  geom_point(size = 4) +
  stat_ellipse() +
  xlab(paste0("PC1 (", sol_pennelli_PCA_variance_perPC[1], "%", ")")) +
  ylab(paste0("PC2 (", sol_pennelli_PCA_variance_perPC[2], "%", ")")) +
  labs(title = "PCA Plot",
       caption = paste0("Produced on", Sys.time())) +
  coord_fixed() +
  theme_bw()

