library(tidyverse)
library(stringr)
library(rtracklayer)
persicum_annotation_file

persicum_annotation_file.filter.tmp <- filter(persicum_annotation_file, type== "mRNA")


persicum_annotation_filter.ids <- persicum_annotation_file.filter.tmp           %>% 
                      mutate(attributes = strsplit(attributes, ";"))            %>%
                      unnest(attributes)                                        %>%
                      separate(attributes, into = c("key", "value"), sep = "=") %>%
                      filter(key %in% c("ID", "Parent"))                        %>%
                      spread(key, value)
  
  
persicum_annotation_filter.ids  

persicum_annotation_pure.ids <- persicum_annotation_filter.ids %>%
                                mutate(target_id = str_replace(ID, "mRNA:", ""),
                                       gene_name = str_replace(Parent, "gene:", ""))

persicum_annotation_pure.ids    

persicum_annotation_selected <- persicum_annotation_pure.ids %>% select(target_id, gene_name)

persicum_annotation_selected

pennelli_annotation_file

### demarcation for pennelli

penelli_annotation_file_mRNA <- filter(pennelli_annotation_file, type == "mRNA")



pennelli_annotation_file.ids <- penelli_annotation_file_mRNA %>%
                                mutate(attributes = strsplit(attributes, ";")) %>%
                                unnest(attributes) %>%
                                separate(attributes, into = c("key", "value"), sep = "=") %>%
                                filter(key %in% c("ID", "Parent")) %>%
                                group_by(ID, Parent) %>%
                                mutate(unique_id = row_number()) %>%
                                spread(key, value)

pennelli_annotation_file.ids <- penelli_annotation_file_mRNA %>%
  mutate(attributes = strsplit(attributes, ";"))%>%
  unnest(attributes) %>%
  separate(attributes, into = c("key", "value"), sep = "=") %>%
  filter(key %in% c("ID", "Parent")) #%>%


spread_pennelli_annotation_file.ids <- pennelli_annotation_file.ids %>% mutate(unique_id = row_number()) %>% spread(key, value)

spread_pennelli_annotation_file.ids.na_removed <- spread_pennelli_annotation_file.ids %>% 
                                                  mutate(ID = ifelse(!is.na(ID), ID, Parent),
                                                         Parent = ifelse(!is.na(Parent), Parent, ID))

pennelli_annotation_selected <- spread_pennelli_annotation_file.ids.na_removed %>%
                                mutate(target_id = ID, gene_name = Parent)
                                #rename(target_id = ID, gene_name = Parent) %>%
                                #select(target_id, gene_name)

##### Trying rtracklayer

browseVignettes("rtracklayer")

trial2_pennelli_annotation_file <- readGFF("spenn_v2.0_gene_models_annot.gff")

trial3_pennelli_annotation_file <- as_tibble(trial2_pennelli_annotation_file)

fill_parents <- trial3_pennelli_annotation_file$Parent

trial4_mRNA <- filter(trial3_pennelli_annotation_file, type == "mRNA")

trial4_parents <- trial4_mRNA %>% mutate(Parent = map_chr(Parent, 1))

trial4_parents

trial4_parents.renamed <- trial4_parents %>% mutate(target_id = ID, gene_name = Parent)
trial4_parents.renamed
pennelli_annotation_selected <- trial4_parents.renamed %>% select(target_id, gene_name)

pennelli_annotation_selected

fill_parents

trial2_pennelli_annotation_file.df <- as.data.frame(trial2_pennelli_annotation_file)

trial2_pennelli_annotation_file.df.mRNA <- trial2_pennelli_annotation_file.df %>% filter(type == "mRNA")
trial2_pennelli_annotation_file.df.mRNA

trial2_pennelli_annotation_file.df.mRNA.renamed <- trial2_pennelli_annotation_file.df.mRNA %>% mutate(target_id = ID, gene_name = Parent) 
trial2_pennelli_annotation_file.df.mRNA.renamed

trial2_pennelli_annotation_file.df.mRNA.renamed.picked <- trial2_pennelli_annotation_file.df.mRNA.renamed %>% select(target_id, gene_name)

trial2_pennelli_annotation_file.df.mRNA.renamed.picked

write_csv(trial2_pennelli_annotation_file.df.mRNA.renamed.picked, "~/Downloads/pennelli_annotation.csv")
head(trial2_pennelli_annotation_file.df.mRNA.renamed.picked)

check <- read_csv("~/Downloads/pennelli_annotation.csv")
check

study_design

study_design_meta_data <- study_design[, c("Run", "population")]
study_design_meta_data

study_design_meta_data.renamed <- study_design_meta_data %>% dplyr::rename(sample_id = Run)

write_csv(study_design_meta_data.renamed, "~/Downloads/sample_id.csv")

sol_persicum_counts_matrix

sol_persicum_counts_df <- as_tibble(sol_persicum_counts_matrix, rownames="gene_id")

colnames(sol_persicum_counts_df) <- c("gene_id", study_design$Run)
sol_persicum_counts_df

write_csv(sol_persicum_counts_df, "~/Downloads/sample_id_sol_persicum_counts.csv")


sol_pennelli_counts_df <- as_tibble(sol_pennelli_counts_matrix, rownames="gene_id")
colnames(sol_pennelli_counts_df) <- c("gene_id", study_design$Run)

write_csv(sol_pennelli_counts_df, "~/Downloads/sample_id_sol_pennelli_counts.csv")
