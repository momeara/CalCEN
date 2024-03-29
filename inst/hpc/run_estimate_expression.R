#!/usr/bin/env Rscript
# -*- tab-width:2;indent-tabs-mode:f;show-trailing-whitespace:t;rm-trailing-spaces:t -*-




library(optparse)
library(assertthat)
library(plyr)
library(dplyr)
library(readr)
library(CalCEN)

cat("Current working directory: ", getwd(), "\n", sep = "")

parameters <- CalCEN::load_parameters()


option_list <- list(
  optparse::make_option(
		c("--runs_fname"), action = "store", dest = "runs_fname",
		help = "TSV file of all RNA-seq runs [Default \"%default\"]"),
  optparse::make_option(
		c("--run_id"), action = "store", dest = "run_id",
		help = "index into the runs table to select which run to analyze in this job [Default \"%default\"]"),
  optparse::make_option(
		c("--results_dir"), action = "store", dest = "results_dir",
		help = "Where to copy the expression results to once they have been estimated [Default \"%default\"]"),
  optparse::make_option(
		c("--logs_dir"), action = "store", dest = "logs_dir",
		help = "Where to write logs [Default \"%default\"]"),
  optparse::make_option(
		c("--work_dir"), action = "store", dest = "work_dir",
		help = "Where to store temporary data while estimating expression [Default \"%default\"]"))

opt <- optparse::OptionParser(option_list = option_list) %>%
		optparse::parse_args()

# validate input
assertthat::assert_that(
		file.exists(opt$runs_fname),
		msg = paste0("Runs file '", opt$runs_fname, "' does not exist."))

ca_run <- readr::read_tsv(opt$runs_fname) %>%
	dplyr::slice(as.numeric(opt$run_id))

assertthat::assert_that(
	"run_accession" %in% names(ca_run),
	"sra_fname" %in% names(ca_run),
	"is_paired" %in% names(ca_run),
	msg = paste0(
			"Runs file '", opt$runs_fname, "' must have columns\n  ['run_accession', 'sra_fname', 'is_paired']\n",
			"Instead it has columns\n  ['", paste0(names(ca_run), collapse = "', '"), "']\n"))

assertthat::assert_that(
		file.exists(
				paste0(parameters$data_paths$reference_genome_path, ".transcripts.fa")),
		msg = paste0(
				"reference genome path ",
				"'", paste0(parameters$data_paths$reference_genome_path, ".transcripts.fa"), "'",
				" does not exist. Make sure to run "))

assertthat::assert_that(
		file.exists(parameters$software_paths$fastq_dump_program),
		msg = paste0(
				"fastq dump program ",
				"'", parameters$software_paths$fastq_dump_program, "'",
				" does not exist."))

assertthat::assert_that(
		file.exists(parameters$software_paths$rsem_calculate_expression_program),
		msg = paste0(
				"RSEM program ",
				"'", parameters$software_paths$rsem_calculate_expression_program, "'",
				" does not exist."))

assertthat::assert_that(
		file.exists(parameters$software_paths$bowtie2_path),
		msg = paste0(
				"Bowtie2 program ",
				"'", parameters$software_paths$bowtie2_path, "'",
				" does not exist."))


timing <- system.time({

	CalCEN::estimate_expression(
		run_accession = ca_run$run_accession[1],
		sra_fname = ca_run$sra_fname[1],
		is_paired = ca_run$is_paired[1],
		results_dir = opt$results_dir[1],
		logs_dir = opt$logs_dir[1],
		work_dir = opt$work_dir[1],
		reference_genome_path = parameters$data_paths$reference_genome_path,
		fastq_dump_program = parameters$software_paths$fastq_dump_program,
		rsem_calculate_expression_program =
				parameters$software_paths$rsem_calculate_expression_program,
		bowtie2_path = parameters$software_paths$bowtie2_path)

})

cat("Run time: ", timing[3], "\n", sep = "")
