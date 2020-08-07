# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:

library(plyr)
library(dplyr)
library(stringr)
library(readr)

source("scripts/parameters.R")


ca_runs <- readr::read_tsv("product/ca_runs_180608.tsv")

get_done_runs <- function(path){
	cat("Getting done runs from '", path, "' ...\n", sep="")

	done_run_results <- list.files(
		path=path,
		pattern="*.genes.results") %>%
		stringr::str_extract("^[^.]+") %>%
		tibble::data_frame(run_accession = .)

	done_run_logs <- list.files(
		path=paste0(path, "/logs"),
		pattern="*.log") %>%
		stringr::str_extract("^[^.]+") %>%
		tibble::data_frame(run_accession = .)

	done_runs <- done_run_results %>%
		dplyr::inner_join(done_run_logs, by="run_accession")
}

done_runs <- rbind(
	get_done_runs("intermediate_data/estimated_expression"),
	get_done_runs("intermediate_data/estimated_expression_180611"))

todo_runs <- ca_runs %>%
	dplyr::anti_join(
		done_runs, by="run_accession")

# this is read by run_estimate_expression.R and
# indexed into by the SGE_TASK_ID from the run_estimate_expression-wrapper.sh
todo_runs %>%
	readr::write_tsv("intermediate_data/todo_runs_180611.tsv")



#### SUBMIT THEM IN PARALLEL ####

n_jobs <- nrow(todo_runs)
tag <- "180611"

cat("# Submitting ", n_jobs, " estimate expression jobs in parallel\n", sep="")
cmd_str <- paste0("qsub -t 1-", n_jobs, " scripts/run_estimate_expression-wrapper.sh ", getwd(), " ", tag)
cat(cmd_str, "\n")
system(cmd_str)
cat("Monitor progress: qstat\n")
cat("Check results when done: intermediate_data/estimated_expression_180611/logs\n")

#### SUBMIT THEM IN SERIAL ####

dir.create(paste0("/nfs/home/momeara/work/collaborations/ca_coexp/intermediate_data/estimated_expression_", tag))
dir.create(paste0("/nfs/home/momeara/work/collaborations/ca_coexp/intermediate_data/estimated_expression_", tag, "/logs"))


todo_runs %>% plyr::a_ply(1, function(ca_run){
	cat("Estimating expression for ", ca_run$run_accession[1], "\n", sep="")

	work_dir <- paste0("/scratch/momeara/ca_coexp/estimated_expression_180611/estimated_expression_", ca_run$run_accession)
	dir.create(work_dir)

	estimate_expression(
			run_accession=ca_run$run_accession[1],
			sra_fname=ca_run$sra_fname[1],
			is_paired=ca_run$is_paired[1],
			results_dir=paste0("/nfs/home/momeara/work/collaborations/ca_coexp/intermediate_data/estimated_expression_", tag),
			logs_dir=paste0("/nfs/home/momeara/work/collaborations/ca_coexp/intermediate_data/estimated_expression_", tag, "/logs",
			work_dir=work_dir)
})

z <- list.files(paste0("/nfs/home/momeara/work/collaborations/ca_coexp/intermediate_data/estimated_expression_", tag, "/logs")) %>%
	stringr::str_extract("^[^.]+") %>%
	data.frame(run_accession = .) %>%
	dplyr::mutate(
		temp_results = paste0("/nfs/home/momeara/work/collaborations/ca_coexp/intermediate_data/estimated_expression_", tag, "/", run_accession, ".genes.results"),
		perm_results = paste0("/nfs/home/momeara/work/collaborations/ca_coexp/intermediate_data/estimated_expression/", run_accession, ".genes.results"),
		temp_log = paste0("/nfs/home/momeara/work/collaborations/ca_coexp/intermediate_data/estimated_expression_", tag, "/logs/", run_accession, ".log"),
		perm_log = paste0("/nfs/home/momeara/work/collaborations/ca_coexp/intermediate_data/estimated_expression/logs/", run_accession, ".log")) %>%
	plyr::a_ply(1, function(result){
		file.copy(result$temp_results[1], result$perm_results[1])
		file.copy(result$temp_log[1], result$perm_log[1])
	})



##### Clean failed runs ####

done_run_logs <- list.files(
	path="intermediate_data/estimated_expression/logs",
	pattern="*.log") %>%
	plyr::adply(1, function(log_fname){
		lines <- readr::read_lines(
			file=paste0("intermediate_data/estimated_expression/logs/", log_fname))

		browser()
		fail <- lines %>%
			stringr::str_detect("An error occurred during processing.") %>%
			any

		if(fail){
			cat("Found failed run: '", log_fname, "'\n", sep="")
		}
	})



