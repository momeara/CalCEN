
library(CalCEN)
parameters <- source("parameters.R")

parameters <- CalCEN::load_parameters()
genome_params <- parameters$source_data$genome

base_dir <- paste0(parameters$data_paths$scratch_dir, "/reference_genome")

if (!dir.exists(paths = base_dir)) {
    cat("Creating ", base_dir, "\n", sep = "")
    dir.create(base_dir, recursive = TRUE)
}

cmd <- paste0(
parameters$software_paths$rsem_path, "/rsem-prepare-reference \\
  --num-threads 10 \\
  --bowtie2 --bowtie2-path ", parameters$software_paths$bowtie2_path, " \\
  ", genome_params$transcript_fasta_path, " \\
  ", base_dir, "/", genome_params$species_short_name)

cat(cmd, "\n", sep = "")
system(cmd)
