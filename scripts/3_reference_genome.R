# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:


base_dir="/nfs/ex7/work/momeara/ca_coexp/reference_genome"

cmd(paste0("
mkdir -p ", base_dir, "
cd ", base_dir, "

curl -O http://www.candidagenome.org/download/sequence/C_albicans_SC5314/Assembly22/current/C_albicans_SC5314_version_A22-s07-m01-r57_default_coding.fasta.gz

rsem-prepare-reference \
  --num-threads 10 \
  --bowtie2 --bowtie2-path /mnt/nfs/home/momeara/opt/bowtie2-2.3.4.1-linux-x86_64 \
  C_albicans_SC5314_version_A22-s07-m01-r57_default_coding.fasta  \
  /nfs/ex7/work/momeara/ca_coexp/reference_genome/SC5314_reference
")
