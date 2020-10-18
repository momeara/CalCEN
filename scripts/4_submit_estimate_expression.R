# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:

library(plyr)
library(dplyr)
library(stringr)
library(readr)

source("parameters.R")

source("scripts/estimate_expression.R")

#ca_runs <- readr::read_tsv("product/ca_runs_180608.tsv")
ca_runs <- readr::read_tsv("product/ca_runs_200928.tsv")

tag <- "20201007"


get_done_runs <- function(path, check_for_logs=TRUE){
	cat("Getting done runs from '", path, "' ...\n", sep="")

	done_run_results <- list.files(
		path=path,
		pattern="*.genes.results") %>%
		stringr::str_extract("^[^.]+") %>%
		tibble::tibble(run_accession = .)

	if(check_for_logs){
		done_run_logs <- list.files(
			path=paste0(path, "/logs"),
			pattern="*.log") %>%
			stringr::str_extract("^[^.]+") %>%
			tibble::tibble(run_accession = .)

		done_runs <- done_run_results %>%
				dplyr::inner_join(done_run_logs, by="run_accession")
	} else{
		done_runs <- done_run_results
	}
}

# some old runs
done_runs <- rbind(
	get_done_runs("intermediate_data/estimated_expression"),
	get_done_runs("intermediate_data/estimated_expression_180611"))

# somehow the logs got lost for these so far
done_runs <- get_done_runs(
		paste0("intermediate_data/estimated_expression_", tag),
		check_for_logs = FALSE)

todo_runs <- ca_runs %>%
	dplyr::anti_join(
		done_runs, by="run_accession")

# this is read by run_estimate_expression.R and
# indexed into by the arrayed from the run_estimate_expression_<cluster_type>_wrapper.sh
todo_runs %>%
		readr::write_tsv(file = paste0("intermediate_data/todo_runs_", tag, ".tsv"))
n_jobs <- nrow(todo_runs)

cat("# Submitting ", n_jobs, " estimate expression jobs in parallel\n", sep="")

if(cluster_type == "SGE"){
		cmd_str <- paste0(
				"qsub -t 1-", n_jobs, " ",
				"scripts/run_estimate_expression_SGE_wrapper.sh ",
				getwd(), " ",
				tag)
		info_message <- "Monitor progress with 'qstat'"

		cat(cmd_str, "\n")
		system(cmd_str)
		cat(info_message, "\n", sep = "")
		cat("Check results when done: intermediate_data/estimated_expression_", tag, "/logs\n", sep = "")

} else if(cluster_type == "SLURM"){

		job_dir <- paste0(scratch_dir, "/estimated_expression_", tag)
		if(!dir.exists(job_dir)){
				cat("Creating job directory: ", job_dir, "\n", sep = "")
				dir.create(job_dir, recursive = TRUE)
		}

		cmd_str <- paste0(
				"sbatch ",
				"--account=", slurm_account, " ",
				"--mail-user=maom@umich.edu ",
				"--mail-type=BEGIN,END,FAIL ",
				"--array=1-", n_jobs, " ",
				"--output='", job_dir, "/%j.log' ",
				"--export=",
				  "TAG='", tag, "',",
				  "BASE_DIR='", base_dir, "',",
				  "JOB_DIR='", job_dir, "' ",
				"scripts/run_estimate_expression_SLURM_wrapper.sh")
		info_message <- "Monitor progress with 'squeue'"

		cat(cmd_str, "\n")
		system(cmd_str)
		cat(info_message, "\n", sep = "")
		cat("Check results when done: intermediate_data/estimated_expression_", tag, "/logs\n", sep = "")

} else if(cluster_type == "LOCAL"){

    dir.create(paste0(getwd(), "/intermediate_data/estimated_expression_", tag))
    dir.create(paste0(getwd(), "/intermediate_data/estimated_expression_", tag, "/logs"))

    todo_runs %>% plyr::a_ply(1, function(ca_run){
    	cat("Estimating expression for ", ca_run$run_accession[1], "\n", sep="")

    	work_dir <- paste0(scratch_dir, "/estimated_expression_", tag, "/estimated_expression_", ca_run$run_accession)
    	dir.create(work_dir, recursive=TRUE)

    	estimate_expression(
    			run_accession = ca_run$run_accession[1],
    			sra_fname = ca_run$sra_fname[1],
    			is_paired = ca_run$is_paired[1],
    			results_dir = paste0(getwd(), "/intermediate_data/estimated_expression_", tag),
    			logs_dir = paste0(getwd(), "/intermediate_data/estimated_expression_", tag, "/logs"),
    			work_dir = work_dir,
    			reference_genome_path = reference_genome_path,
    			fastq_dump_program = fastq_dump_program,
    			rsem_calculate_expression_program= rsem_calculate_expression_program,
    			bowtie2_path = bowtie2_path)
    })

    z <- list.files(paste0(getwd(), "/intermediate_data/estimated_expression_", tag, "/logs")) %>%
    	stringr::str_extract("^[^.]+") %>%
    	data.frame(run_accession = .) %>%
    	dplyr::mutate(
    		temp_results = paste0(getwd(), "/intermediate_data/estimated_expression_", tag, "/", run_accession, ".genes.results"),
    		perm_results = paste0(getwd(), "/intermediate_data/estimated_expression/", run_accession, ".genes.results"),
    		temp_log = paste0(getwd(), "/intermediate_data/estimated_expression_", tag, "/logs/", run_accession, ".log"),
    		perm_log = paste0(getwd(), "/intermediate_data/estimated_expression/logs/", run_accession, ".log")) %>%
    	plyr::a_ply(1, function(result){
    		file.copy(result$temp_results[1], result$perm_results[1])
    		file.copy(result$temp_log[1], result$perm_log[1])
    	})
}

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



