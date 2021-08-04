# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:

source("parameters.R")

base_dir <- paste0(scratch_dir, "/reference_genome")

system(paste0("
mkdir -p ", base_dir, "
cd ", base_dir, "

curl -O http://www.candidagenome.org/download/sequence/C_albicans_SC5314/Assembly22/current/C_albicans_SC5314_version_A22-s07-m01-r57_default_coding.fasta.gz

rsem-prepare-reference \\
  --num-threads 10 \\
  --bowtie2 --bowtie2-path ", bowtie2_path, " \\
  C_albicans_SC5314_version_A22-s07-m01-r57_default_coding.fasta  \\
  ", base_dir, "/SC5314_reference
"))
