#' Estimate expression data from SRA run and reference genome
#'
#'
#' Example
#'
#' tag <- "20201209"
#' estimate_expression(
#'    run_accession = "SRR513613",
#'    sra_fname = "/scratch/maom_root/maom99/maom/CryptoCEN/sra/SRR513613/SRR513613.sra",
#'    is_paired = TRUE,
#'    results_dir = "intermediate_data/estimated_expression_20201209",
#'    results_dir = "intermediate_data/estimated_expression_20201209/logs",
#'    work_dir = "/scratch/maom_root/maom99/maom/CryptoCEN/estimated_expression_20201209",
#'    reference_genome_path = "/scratch/maom_root/maom99/maom/CryptoCEN/reference_geneome",
#'    fastq_dump_program = "
#' @param run_accession identifier for exprsesion run
#' @param sra_fname path to .sra file for the run
#' @param is_paired is it paired end read expression data?
#' @param results_dir path where the results should be stored
#' @param logs_dir path where the logs should be stored
#' @param work_dir path where intermediate data should be stored
#' @param reference_genome_path path to reference genome
#' @param fastq_dump_program path to fastq_dump
#' @param rsem_calculate_expression_program path to rsem-calculate-expression
#'
#' 
#' @export
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
    n_threads = 1) {
    log_fname <-  paste0(work_dir, "/", run_accession, ".log")
    cat("Writing logs to ", log_fname, "\n", sep = "")

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
        sep = "", file = log_fname, append = TRUE)

    run_cmd <- function(info, cmd_str, log_fname=NULL) {
        cmd_str <- paste0("cd ", work_dir, " && ", cmd_str)
        if (!is.null(log_fname)) {
                cat("\n", sep = "", file = log_fname, append = TRUE)
                cat("# ", info, "\n", sep = "", file = log_fname, append=TRUE)
                cat(cmd_str, "\n", sep = "", file = log_fname, append=TRUE)
                system(paste0(cmd_str, " >> ", log_fname, " 2>&1"))
                cat("\n\n", file = log_fname, append = TRUE)
        } else {
                cat("# ", info, "\n", sep = "")
                cat(cmd_str, "\n", sep = "")
                system(cmd_str)
        }
    }

    timing <- system.time({
        run_cmd(
            info = "Copying SRA file to working directory",
            cmd_str = paste0("cp ", sra_fname, " ", work_dir),
            log_fname = log_fname)

        if (!is_paired) {
            run_cmd(
                info = "Convert .sra to .fastq assuming single-end reads",
                cmd_str = paste0(fastq_dump_program, " --gzip --skip-technical  --readids --read-filter pass --dumpbase --clip ", run_accession, ".sra"),
                log_fname = log_fname)
        } else {
                # there is a typo in the fastq-dump command line arguments for version 2.10.8
                # https://github.com/ncbi/sra-tools/issues/381
                good_sra_version <- system(
                        command = paste0(fastq_dump_program, " --version"),
                        intern = TRUE)[2] %>%
                        stringr::str_extract("[0-9.]+$") %>%
                        compareVersion("2.10.8")
            if (good_sra_version) {
                split_flag <- "--split-3"
            } else {
                split_flag <- "--split-e"
            }
            run_cmd(
                info = "Convert .sra to .fastq assuming paired-end reads",
                cmd_str = paste0(fastq_dump_program, " --gzip --skip-technical  --readids --read-filter pass --dumpbase ", split_flag, " --clip ", run_accession, ".sra"),
                log_fname = log_fname)
        }


        if (!is_paired) {
            run_cmd(
                info = "Estimate expression levels assuming single-end reads",
                cmd_str = paste0(rsem_calculate_expression_program, " -p ", n_threads, " --no-bam-output --estimate-rspd --bowtie2 --bowtie2-path ", bowtie2_path, " --append-names ", run_accession, "_pass.fastq.gz ", reference_genome_path, " ", run_accession),
                log_fname = log_fname)
        } else {
            run_cmd(
                info = "Estimate expression levels assuming paired-end reads",
                cmd_str = paste0(rsem_calculate_expression_program, " -p ", n_threads, " --paired-end --no-bam-output --estimate-rspd --bowtie2 --bowtie2-path ", bowtie2_path, " --append-names ", run_accession, "_pass_1.fastq.gz ", run_accession, "_pass_2.fastq.gz ", reference_genome_path, " ", run_accession),
                log_fname = log_fname)
        }


      run_cmd(
          info = "Copy results",
          cmd_str = paste0(
              "cp ", work_dir, "/", run_accession, ".genes.results ",
              results_dir, "/"),
            log_fname = log_fname)

    }) # timing

    cat("# Runtime: ", timing[3], "\n", sep = "", file = log_fname, append = TRUE)


    run_cmd(
        info = "Copy log file",
        cmd_str=paste0("cp ", log_fname, " ", logs_dir),
        log_fname = NULL)

    run_cmd(
      info = "Removing working files",
        cmd_str = paste0("rm ", work_dir, "/", run_accession, ".*"),
        log_fname = NULL)
}

#' Submit estimate gene expression job to SLURM cluster
#'
#'   Reads gene expression estimate jobs from
#' 
#'      intermediate_data/todo_runs_<tag>.tsv
#'          # tab separated table with column [run_accession]
#'
#'   outputs results to
#'
#'       intermediate_data/estimate_expression_<tag>
#'
#'   This submission script executes
#'
#'       cd <base_dir>
#'       sbatch \
#'         --array 1-<n_runs> \
#'         --export=TAG='<ta>',BASE_DIR='<base_dir>',JOB_DIR='<job_dir>' \
#'         ~/opt/CalCEN/inst/run_estimate_expression_SLURM_wrapper.sh . <tag>
#'
#'   A demon of running from the shell
#' 
#'      export BASE_DIR=/home/maom/opt/CalCEN
#'      export JOB_DIR=/scratch/maom_root/maom99/maom/CalCEN/estimate_expression_20201007
#'      export SLURM_ARRAY_TASK_ID=1
#'      ~/opt/CalCEN/inst/run_estimate_expression_SLURM_wrapper.sh
#' 
#'   this will set:
#' 
#'      TAG: 20201007
#'      BASE_DIR: /home/maom/opt/CalCEN/vignettes/<network>
#'      SLURM_ARRAY_TASK_ID: 1
#'      Rscript: /home/maom/opt/bin/Rscript
#'      RUNS_FNAME: /home/maom/opt/CalCEN/intermediate_data/todo_runs_20201007.tsv
#'      LOGS_DIR: /home/maom/opt/CalCEN/intermediate_data/estimated_expression_20201007/logs
#'      TASK_NAME: estimated_expression_1
#'      JOB_DIR: /scratch/maom_root/maom99/maom/CalCEN/estimated_expression_20201007
#'      TASK_DIR: /scratch/maom_root/maom99/maom/CalCEN/estimated_expression_20201007/estimated_expression_1
#' 
#'  and then call
#' 
#'      Rscript ~/opt/CalCEN/scripts/run_estimate_expression.R \
#'          --runs_fname ${BASE_DIR}intermediate_data/todo_runs_20201007.tsv \
#'          --run_id 1 \
#'          --results_dir ${BASE_DIR}/intermediate_data/estimated_expression_20201007 \
#'          --logs_dir ${BASE_DIR}/intermediate_data/estimated_expression_20201007/logs \
#'          --work_dir ${TASK_DIR}/estimated_expression_1
#' 
#' @param tag label for scratch and result directories
#' @param base_dir directory for project
#' @param scratch_dir scratch directory for the project
#' @param slurm_account slurm account
#' @param slurm_mail_user email for slurm job updates at BEGIN,END
#' @param slurm_partition partition on slurm cluster where to run the job
#' 
#' @export
submit_estimate_expression_slurm <- function(
    tag,
    base_dir,
    scratch_dir,
    slurm_account,
    slurm_mail_user,
    slurm_partition) {

    if (!dir.exists(job_dir)) {
            cat("Creating job directory: ", job_dir, "\n", sep = "")
            dir.create(job_dir, recursive = TRUE)
    }

    todo_runs
    cmd_str <- paste0(
            "sbatch ",
            "--account=", slurm_account, " ",
            "--mail-user=", slurm_mail_user, " ",
            "--mail-type=BEGIN,END,FAIL ",
            "--array=1-", n_jobs, " ",
            "--output='", job_dir, "/%j.log' ",
            "--time=03:00:00 ",
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
}        



#' Validate estimated expression runs
#'
#' After running estimate_expression(..., results_dir = <results_dir>, ...)
#' call validate_estimated_runs(results_dir = <results_dir>)
#' to get the runs that have been completed
#' 
#' @param results_dir directory to look for completed estimated expression runs
#' @param check_for_logs check for logs
#' @return data.frame with column [run_accession] of runs that have been completed
#' @export
validate_estimated_runs <- function(
    path,
    check_for_logs = TRUE) {
    cat("Getting estimated expression runs from '", path, "' ...\n", sep = "")
    done_run_results <- list.files(
        path = path,
        pattern = "*.genes.results") %>%
        stringr::str_extract("^[^.]+") %>%
        tibble::tibble(run_accession = .)
    if (check_for_logs) {
        done_run_logs <- list.files(
            path = paste0(path, "/logs"),
            pattern = "*.log") %>%
            stringr::str_extract("^[^.]+") %>%
            tibble::tibble(run_accession = .)
        done_runs <- done_run_results %>%
                dplyr::inner_join(done_run_logs, by="run_accession")
    } else{
        done_runs <- done_run_results
    }
}
