
library(CalCEN)

parameters <- CalCEN::load_parameters()


cmd <- paste0(
  parameters$software_paths$rsem_path, "/rsem-prepare-reference \\
  --num-threads 10 \\
  --bowtie2 --bowtie2-path ", parameters$software_paths$bowtie2_path, " \\
  ", parameters$source_data$genome$transcript_fasta_path, " \\
  ", parameters$data_paths$reference_genome_path)
cat(cmd, "\n", sep = "")

system(cmd)
