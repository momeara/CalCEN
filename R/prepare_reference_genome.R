
#' Use RSEM to prepare reference genome
#'
#' use rsem-prepare-reference to generate reference genome
#' for mapping RNAseq transcripts
#' 
#' # Usage example:
#' in parameters.yaml
#'
#'   source_data:
#'     genome:
#'        transcript_fasta_path: "raw_data/FungiDB-49_CneoformansH99_AnnotatedTranscripts.fasta"
#'   data_paths:
#'     reference_genome_path: "/scratch/maom_root/maom99/maom/CryptoCEN/reference_genome/CneoformansH99"
#'   software_paths:
#'     bowtie2_path: "/home/maom/opt/bowtie2-2.3.4.1-linux-x86_64"
#'     fastq_dump_program: "/home/maom/opt/bin/fastq-dump"
#'     rsem_path: "/home/maom/opt/RSEM/"
#'
#' In R:
#'
#'   library(CalCEN)
#'   parameters <- CalCEN::load_parameters()
#'   CalCEN::prepare_reference_genome()
#' 
#' @param parameters as read in by load_parameters()
#' @param verbose (default: TRUE)
#' @export
prepare_reference_genome <- function(
  parameters,
  verbose = TRUE) {

    assertthat::assert_that(
      file.exists(paste0(
        parameters$software_paths$rsem_path,
        "/rsem-prepare-reference")),
      msg = paste0(
        "RSEM program ",
        "'", paste0(
          parameters$software_paths$rsem_path,
          "/rsem-prepare-reference"),
        " does not exist."))
    
    assertthat::assert_that(
      file.exists(parameters$software_paths$bowtie2_path),
      msg = paste0(
          "Bowtie2 program ",
          "'", parameters$software_paths$bowtie2_path, "'",
          " does not exist."))

    assertthat::assert_that(
      file.exists(parameters$source_data$genome$transcript_fasta_path),
      msg = paste0(
          "Genome transcripts ",
          "'", parameters$source_data$genome$transcript_fasta_path, "'",
          " does not exist."))

    base_dir <- parameters$data_paths$reference_genome_path %>%
        dirname()
    if (!dir.exists(paths = base_dir)) {
        cat("Creating output path '", base_dir, "'\n", sep = "")
        dir.create(base_dir, recursive = TRUE)
    } else {
        cat(
            "WARNING: reference genome base dir ",
            "'", base_dir, "' already exists.\n", sep = "")
    }
    assertthat::assert_that(
      dir.exists(base_dir),
      msg = paste0(
          "reference genome base dir ",
          "'", base_dir, "'",
          " does not exist."))

    cmd_str <- paste0(
        parameters$software_paths$rsem_path, "/rsem-prepare-reference \\\n",
        "  --num-threads 10 \\\n",
        "  --bowtie2 ",
        "--bowtie2-path ", parameters$software_paths$bowtie2_path, " \\\n",
        "  ", parameters$source_data$genome$transcript_fasta_path, " \\\n",
        "  ", parameters$data_paths$reference_genome_path)

    if (verbose) {
        cat(cmd_str, "\n", sep = "")
    }

    system(cmd_str)
}
