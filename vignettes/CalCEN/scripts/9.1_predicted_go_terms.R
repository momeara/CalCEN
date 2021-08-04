# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:

library(plyr)
library(dplyr)
library(tidyr)
library(tibble)
library(EGAD)
library(GO.db)

source("scripts/evaluate_network_gba.R")

load("intermediate_data/ca_go_annotations.Rdata")
load("product/CalCEN_network_full_20201007.Rdata")

go_terms <- GO.db::GO_dbconn() %>%
    dplyr::tbl("go_term") %>%
    dplyr::collect(n = Inf) %>%
    dplyr::select(-`_id`)

CalCEN_full_go_pred <- EGAD::run_GBA(CalCEN_network_full, ca_go_annotations)

CalCEN_go_pred <- CalCEN_full_go_pred[[2]] %>%
    tibble::as_tibble(rownames = "feature_name") %>%
    tidyr::pivot_longer(
        cols = -c("feature_name"),
        names_to = "go_id",
        values_to = "score") %>%
    dplyr::left_join(
        go_terms,
        by = "go_id") %>%
    dplyr::filter(!is.na(term)) %>%
    dplyr::group_by(go_id) %>%
    dplyr::mutate(score = score / sum(score)) %>%
    dplyr::mutate(score = score * 100000) %>%
    dplyr::ungroup()
    

CalCEN_go_pred %>%
    readr::write_tsv("product/CalCEN_full_predicted_go_terms_20210121.tsv")

CalCEN_go_pred_gene_top5 <- CalCEN_go_pred %>%
    dplyr::group_by(feature_name) %>%
    dplyr::arrange(dplyr::desc(score)) %>%
    dplyr::slice(1:5) %>%
    dplyr::ungroup()

CalCEN_go_pred_gene_top5 %>%
    readr::write_tsv("product/CalCEN_full_predicted_go_terms_gene_top5_20210121.tsv")

