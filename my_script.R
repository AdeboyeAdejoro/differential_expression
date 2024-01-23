library(tidyverse)

targets

sra_accession <- targets$sra_accession
myCounts

myCounts.df <- as_tibble(counts_kal, rownames="Gene ID")

colnames(myCounts.df) <- c("Gene ID", sra_accession)

myCounts.df

write_csv(myCounts.df, "~/Desktop/course_dataset/kallisto_output/myCounts.df.csv")


targets

write_csv(targets, "~/Desktop/course_dataset/kallisto_output/experiment-info.csv")

targets_2 <- targets[, 2-3]

colnames(targets_2) = c("sample_id", "group")

targets_2

write_csv(targets_2, "~/Desktop/course_dataset/kallisto_output/conditions.csv")

Txi_gene
my_raw_counts <- Txi_gene$counts

my_raw_counts <- as_tibble(my_raw_counts, rownames="Gene ID")

 colnames(my_raw_counts) <- c("Gene ID", sra_accession)
 
 my_raw_counts
 
 write_csv(my_raw_counts, "~/Desktop/course_dataset/kallisto_output/raw_counts.csv")
