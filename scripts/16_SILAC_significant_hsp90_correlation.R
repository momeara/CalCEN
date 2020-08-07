# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:

library(plyr)
library(dplyr)

load("intermediate_data/ca_silac_hsp90_fold_change.Rdata")


threshold <- 1.5
silac_signif_hits <- ca_silac_hsp90_fold_change %>%
    dplyr::mutate(
        pharm_diff = (
            ((wildtype == 0) & (pharmacological > 0)) |
            ((wildtype > 0) & (pharmacological/wildtype >= threshold)) |
            ((pharmacological == 0) & (wildtype > 0)) |
            ((pharmacological > 0) & (wildtype/pharmacological >= threshold))),
        gene_diff = (
            ((wildtype == 0) & (genetic > 0)) |
            ((wildtype > 0) & (genetic/wildtype >= threshold)) |
            ((genetic == 0) & (wildtype > 0)) |
            ((genetic > 0) & (wildtype/genetic >= threshold)))) %>%
    dplyr::filter(pharm_diff > 0 | gene_diff > 0) %>%
    dplyr::group_by(gene_label) %>%
    dplyr::summarize(
        pharm_n_diff = sum(pharm_diff),
        gene_n_diff = sum(gene_diff),
        pharm_diffs = paste0("(", fraction_index[pharm_diff], ",", pharmacological[pharm_diff], ",", pharmacological[pharm_diff]/wildtype[pharm_diff], ")", collapse="", sep=""),
        gene_diffs = paste0("(", fraction_index[gene_diff], ",", genetic[gene_diff], ",", genetic[gene_diff]/wildtype[gene_diff], ")", collapse="", sep=""))

silac_signif_hits %>% dplyr::filter(pharm_n_diff > 1) %>% nrow
silac_signif_hits %>% dplyr::filter(gene_n_diff > 1) %>% nrow
silac_signif_hits %>% dplyr::filter(gene_n_diff > 1, pharm_n_diff > 1) %>% nrow

silac_signif_hits %>%
    dplyr::filter(pharm_n_diff > 0 | gene_n_diff > 0) %>%
    readr::write_tsv("product/silac_signif_diff_190407.tsv")
