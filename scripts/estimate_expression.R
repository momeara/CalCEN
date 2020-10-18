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
	reference_genome_path,
	fastq_dump_program,
	rsem_calculate_expression_program,
	bowtie2_path,
	n_threads=1
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

	run_cmd <- function(info, cmd_str, log_fname=NULL){
		cmd_str <- paste0("cd ", work_dir, " && ", cmd_str)
		if(!is.null(log_fname)){
				cat("\n", sep = "", file = log_fname, append = TRUE)
				cat("# ", info, "\n", sep = "", file = log_fname, append=TRUE)
				cat(cmd_str, "\n", sep = "", file = log_fname, append=TRUE)
				system(paste0(cmd_str, " >> ", log_fname, " 2>&1"))
				cat("\n\n", file = log_fname, append=TRUE)
		} else {
				cat("# ", info, "\n", sep="")
				cat(cmd_str, "\n", sep="")
				system(cmd_str)
		}
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
				# there is a typo in the fastq-dump command line arguments for version 2.10.8
				# https://github.com/ncbi/sra-tools/issues/381
				good_sra_version <- system(
						command = paste0(fastq_dump_program, " --version"),
						intern = TRUE)[2] %>%
						stringr::str_extract("[0-9.]+$") %>%
						compareVersion("2.10.8")
			if(good_sra_version){
				split_flag <- "--split-3"
			} else {
				split_flag <- "--split-e"
			}
			run_cmd(
				info="Convert .sra to .fastq assuming paired-end reads",
				cmd_str=paste0(fastq_dump_program, " --gzip --skip-technical  --readids --read-filter pass --dumpbase ", split_flag, " --clip ", run_accession, ".sra"),
				log_fname=log_fname)
		}


		if(!is_paired){
			run_cmd(
				info="Estimate expression levels assuming single-end reads",
				cmd_str=paste0(rsem_calculate_expression_program, " -p ", n_threads, " --no-bam-output --estimate-rspd --bowtie2 --bowtie2-path ", bowtie2_path, " --append-names ", run_accession, "_pass.fastq.gz ", reference_genome_path, " ", run_accession),
				log_fname=log_fname)
		} else {
			run_cmd(
				info="Estimate expression levels assuming paired-end reads",
				cmd_str=paste0(rsem_calculate_expression_program, " -p ", n_threads, " --paired-end --no-bam-output --estimate-rspd --bowtie2 --bowtie2-path ", bowtie2_path, " --append-names ", run_accession, "_pass_1.fastq.gz ", run_accession, "_pass_2.fastq.gz ", reference_genome_path, " ", run_accession),
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
		log_fname=NULL)

	run_cmd(
	  info="Removing working files",
		cmd_str=paste0("rm ", work_dir, "/", run_accession, ".*"),
		log_fname=NULL)
}
