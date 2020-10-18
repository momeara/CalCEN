#!/usr/bin/env Rscript
# -*- tab-width:2;indent-tabs-mode:f;show-trailing-whitespace:t;rm-trailing-spaces:t -*-




library(optparse)
library(plyr)
library(dplyr)
library(readr)

cat("Current working directory: ", getwd(), "\n", sep = "")

source("parameters.R")
source("scripts/estimate_expression.R")

option_list <- list(
  optparse::make_option(
		c("--runs_fname"), action="store", dest="runs_fname",
		help="TSV file of all RNA-seq runs [Default \"%default\"]"),
  optparse::make_option(
		c("--run_id"), action="store", dest="run_id",
		help="index into the runs table to select which run to analyze in this job [Default \"%default\"]"),
  optparse::make_option(
		c("--results_dir"), action="store", dest="results_dir",
		help="Where to copy the expression results to once they have been estimated [Default \"%default\"]"),
  optparse::make_option(
		c("--logs_dir"), action="store", dest="logs_dir",
		help="Where to write logs [Default \"%default\"]"),
  optparse::make_option(
		c("--work_dir"), action="store", dest="work_dir",
		help="Where to store temporary data while estimating expression [Default \"%default\"]"))


opt <- optparse::OptionParser(option_list=option_list) %>%
		optparse::parse_args()

ca_run <- readr::read_tsv(opt$runs_fname) %>%
	dplyr::slice(as.numeric(opt$run_id))

timing <- system.time({

	estimate_expression(
		run_accession = ca_run$run_accession[1],
		sra_fname = ca_run$sra_fname[1],
		is_paired = ca_run$is_paired[1],
		results_dir = opt$results_dir[1],
		logs_dir = opt$logs_dir[1],
		work_dir = opt$work_dir[1],
		reference_genome_path = reference_genome_path,
		fastq_dump_program = fastq_dump_program,
		rsem_calculate_expression_program= rsem_calculate_expression_program,
		bowtie2_path = bowtie2_path)

})

cat("Run time: ", timing[3], "\n", sep="")
