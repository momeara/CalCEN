

install_prerequisites:
	# sratoolkit
	# bowtie2
	# RSEM


make ca_coexp:
	Rscript scripts/1_define_runs.R
	Rscript scripts/2_download_runs.R
	Rscript scripts/3_reference_genome.R
	Rscript scripts/4_submit_estimate_expression.R
