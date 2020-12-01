
# Candida albicans
ncbi_taxon <- "5476"


base_dir <- "/home/maom/opt/ca_coexp"

# where raw and processed sequence data is stored
scratch_dir <- "/scratch/maom_root/maom99/maom/CalCEN"

# cluster scratch used for computing expression profiles
expression_scratch_dir <- "/scratch/maom_root/maom99/maom"

bowtie2_path <- "/home/maom/opt/bowtie2-2.3.4.1-linux-x86_64"
fastq_dump_program <- "/home/maom/opt/bin/fastq-dump" 
rsem_calculate_expression_program <- "/home/maom/opt/RSEM/rsem-calculate-expression"

candida_genome_assembly <- "C_albicans_SC5314_version_A22-s07-m01-r57"
reference_genome_path <- "/scratch/maom_root/maom99/maom/CalCEN/reference_genome/SC5314_reference"

crossmap_path <- "/home/maom/opt/crossmap"


cluster_type <- "SLURM"
slurm_account <- "maom99" 
slurm_mail_user <- "maom@umich.edu"
slurm_partition <- "standard"

