# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:

library(plyr)
library(dplyr)
library(stringr)

source("parameters.R")




estimate_expression <- function(
	run_accession,
	sra_fname,
	is_paired,
	results_dir,
	logs_dir,
	work_dir,
	n_threads=1,
	reference_genome_path="/nfs/ex7/work/momeara/ca_coexp/reference_genome/SC5314_reference",
	fastq_dump_program="~/opt/sratoolkit.2.9.0-ubuntu64/bin/fastq-dump",
	rsem_calculate_expression_program="~/opt/bin/rsem-calculate-expression",
	bowtie2_path="/mnt/nfs/home/momeara/opt/bowtie2-2.3.4.1-linux-x86_64"
){
	log_fname <-  paste0(work_dir, "/", run_accession, ".log")
	cat("Writing logs to ", log_fname, "\n", sep="")

	cat(
		"# Begin estimating the gene expression levels:\n",
		"# run_accession: ", run_accession, "\n",
		"# sra_fname: ", sra_fname, "\n",
		"# read_layout: ", ifelse(is_paired, "paired-end", "single-end"), "\n",
		"# reference_genome: ", reference_genome_path, "\n",
		"# results_dir: ", results_dir, "\n",
		"# logs_dir: ", logs_dir, "\n",
		"# work_dir: ", work_dir, "\n",
		"# n_threads: ", n_threads, "\n\n",
		sep="", file=log_fname, append=TRUE)

	run_cmd <- function(info, cmd_str, log_fname){
		cat("# ", info, "\n", sep="", file=log_fname, append=TRUE)
		cmd_str <- paste0("cd ", work_dir, " && ", cmd_str)
		cat(cmd_str, "\n", sep="", file=log_fname, append=TRUE)
		system(paste0(cmd_str, " >> ", log_fname, " 2>&1"))
		cat("\n\n", file=log_fname, append=TRUE)
	}

	timing <- system.time({
		run_cmd(
			info="Copying SRA file to working directory",
			cmd_str=paste0("cp ", sra_fname, " ", work_dir),
			log_fname=log_fname)

		if(!is_paired){
			run_cmd(
				info="Convert .sra to .fastq assuming single-end reads",
				cmd_str=paste0(fastq_dump_program, " --gzip --skip-technical  --readids --read-filter pass --dumpbase --clip ", run_accession, ".sra"),
				log_fname=log_fname)
		} else {
			run_cmd(
				info="Convert .sra to .fastq assuming paired-end reads",
				cmd_str=paste0(fastq_dump_program, " --gzip --skip-technical  --readids --read-filter pass --dumpbase --split-3 --clip ", run_accession, ".sra"),
				log_fname=log_fname)
		}


		if(!is_paired){
			run_cmd(
				info="Estimate expression levels assuming single-end reads",
				cmd_str=paste0("~/opt/bin/rsem-calculate-expression -p ", n_threads, " --no-bam-output --estimate-rspd --bowtie2 --bowtie2-path ", bowtie2_path, " --append-names ", run_accession, "_pass.fastq.gz ", reference_genome_path, " ", run_accession),
				log_fname=log_fname)
		} else {
			run_cmd(
				info="Estimate expression levels assuming paired-end reads",
				cmd_str=paste0("~/opt/bin/rsem-calculate-expression -p ", n_threads, " --paired-end --no-bam-output --estimate-rspd --bowtie2 --bowtie2-path ", bowtie2_path, " --append-names ", run_accession, "_pass_1.fastq.gz ", run_accession, "_pass_2.fastq.gz ", reference_genome_path, " ", run_accession),
				log_fname=log_fname)
		}


	  run_cmd(
	  	info="Copy results",
	  	cmd_str=paste0(
	  		"cp ", work_dir, "/", run_accession, ".genes.results ",
	  		results_dir, "/"),
			log_fname=log_fname)

	}) # timing

	cat("# Runtime: ", timing[3], "\n", sep="", file=log_fname, append=TRUE)


	run_cmd(
		info="Copy log file",
		cmd_str=paste0("cp ", log_fname, " ", logs_dir),
		log_fname="")

	run_cmd(
	  info="Removing working files",
		cmd_str=paste0("rm ", work_dir, "/", run_accession, ".*"),
		log_fname="")
}
