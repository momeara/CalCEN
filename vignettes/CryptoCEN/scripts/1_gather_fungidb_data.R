

library(tidyverse)
library(seqinr)
library(CalCEN)

parameters <- CalCEN::load_parameters()

# H99 transcripts
cmd <- paste0(
    "cd raw_data && ",
    "wget ", parameters$source_data$genome$transcript_fasta_url)
cat(cmd, "\n", sep = "")
system(cmd)

h99_transcripts <- seqinr::read.fasta(
    file = parameters$source_data$genome$transcript_fasta_url,
    seqtype = "AA")
    
h99_transcript_annotations <- h99_transcripts %>%
    purrr::map_chr(~seqinr::getAnnot(.)) %>%
    data.frame(annotation = .) %>%
    tidyr::separate(
        col = annotation,
        into = c("cnag_id", "gene", "organism", "gene_product", "transcript_product", "location", "length", "sequence_SO", "SO", "is_pseudo"),
        sep = " [|] ") %>%
    dplyr::mutate(
        dplyr::across(
            .cols = everything(),
            ~stringr::str_replace(., "^[a-zA-Z_]+=", ""))) %>%
    dplyr::mutate(
        variant = cnag_id %>% stringr::str_extract("t[0-9]+_[0-9]+$"),
        is_pseudo = ifelse(is_pseudo == "true", TRUE, FALSE)) %>%
    dplyr::select(-gene)

save(
    h99_transcript_annotations,
    file = "intermediate_data/h99_transcript_annotations.Rdata")

h99_transcript_annotations %>%
    readr::write_tsv(file = "product/h99_transcript_annotations_20210724.tsv")
