# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:

library(plyr)
library(reshape2)
library(dplyr)
library(magrittr)
library(stringr)
library(readr)
library(reshape2)
library(EGAD)
library(future)
Sys.setenv(DEBUGME = "batchtools")
library(batchtools)
library(future.batchtools)

future::plan(
    list(
        tweak(
            future.batchtools::batchtools_slurm,
                resources = list(
                    account = "maom99",
                    ntasks=1,
                    ncpus=1L,
                    memory="10GB"),
            template = "inst/batchtools.greatlakes.tmpl")))


load("intermediate_data/ca_runs_final.Rdata")
load("intermediate_data/chromosome_features.Rdata")
load("intermediate_data/estimated_expression.Rdata")
load("intermediate_data/ca_genes.Rdata")

load("intermediate_data/ca_go_annotations.Rdata")

exprs <- reshape2::acast(
    data = estimated_expression %>%
            dplyr::semi_join(ca_runs_final, by = c("run_accession")),
    formula = gene_id ~ run_accession,
    value.var = "FPKM") %>%
    magrittr::extract(ca_genes, )

#############################
# spearman rank correlation #
#############################

studies <- unique(ca_runs_final$study_accession)
n_studies <- ca_runs_final %>%
    dplyr::distinct(study_accession) %>%
    nrow

########
# test #
########
orders <- 3:3
networks_per_order <- 3

ablation_scores <- orders %>%
  	purrr::map_dfr(function(order) {
				cat("Building network sets of order ", order, ":\n", sep = "")
  			utils::combn(x = studies, m = order, simplify = FALSE) %>%
  					sample(size = networks_per_order, replace = FALSE) %>%
  					purrr::map_dfr(function(study_set) {
								run_set <- ca_runs_final %>%
										dplyr::filter(study_accession %in% study_set)
  							cat("  set: [", paste0(study_set, collapse = ", "), "]",
										" nruns: ", nrow(run_set), sep = "")
								network <- EGAD::build_coexp_network(
  									exprs = exprs[, colnames(exprs) %in% run_set$run_accession],
  									gene.list = ca_genes)
  							gba <- EGAD::run_GBA(network, ca_go_annotations)
  							cat(
  									"  GBA: mean: ", gba[[1]][, 1] %>% mean(na.rm = TRUE),
  									" std: ", gba[[1]][, 1] %>% sd(na.rm = TRUE), "\n", sep = "")
  							tibble::tibble(
  									order = order,
										study_set = study_set %>% paste0(collapse = "|"),
  									gba_mean = gba[[1]][, 1] %>% mean(na.rm = TRUE),
  									gba_std = gba[[1]][, 1] %>% sd(na.rm = TRUE))
  					})
  	})



##############
# production #
##############
# use the futures package to calculate the GBA enrichment in parallel

summarize_network_gba <- function(ablation_exprs, study_set, ca_go_annotations) {
		network <- EGAD::build_coexp_network(
				exprs = ablation_exprs,
				gene.list = ca_genes)
		gba <- EGAD::run_GBA(network, ca_go_annotations)
		cat(
				"  GBA: mean: ", gba[[1]][, 1] %>% mean(na.rm = TRUE),
				" std: ", gba[[1]][, 1] %>% sd(na.rm = TRUE), "\n", sep = "")
		tibble::tibble(
				order = length(study_set),
				study_set = study_set %>% paste0(collapse = "|"),
				gba_mean = gba[[1]][, 1] %>% mean(na.rm = TRUE),
				gba_std = gba[[1]][, 1] %>% sd(na.rm = TRUE))
}

summarize_network_gba <- function(ablation_exprs, study_set, ca_go_annotations) {
		tibble::tibble(
				order = 1,
				study_set = "abc",
				gba_mean = .45,
				gba_std = .1)
}


orders <- 1:n_studies
networks_per_order <- n_studies
gba_futures <- orders %>%
  	purrr::map(function(order) {
				cat("Building network sets of order ", order, ":\n", sep = "")
  			utils::combn(x = studies, m = order, simplify = FALSE) %>%
  					sample(size = networks_per_order, replace = FALSE) %>%
  					purrr::map(function(study_set) {
								run_set <- ca_runs_final %>%
										dplyr::filter(study_accession %in% study_set)
  							cat("  set: [", paste0(study_set, collapse = ", "), "]",
										" nruns: ", nrow(run_set), "\n", sep = "")
								future::future(
										summarize_network_gba(
												ablation_exprs = exprs[, colnames(exprs) %in% run_set$run_accession],
												study_set = study_set,
												ca_go_annotations = ca_go_annotations),
										packages = "EGAD")
  					})
  	})

ablation_scores <- gba_futures %>%
		unlist() %>%
		purrr::map_dfr(function(gba_future) {
				future::value(gba_future) %>%
						dplyr::select(
								order, study_set, gba_mean, gba_std)
		})

ablation_scores %>%
		readr::write_tsv("product/go_pred_summary_experiment_ablation_order_20201121.tsv")

plot <- ggplot2::ggplot() +
		ggplot2::theme_bw() +
		ggplot2::geom_smooth(
				data = ablation_scores,
				mapping = ggplot2::aes(
						x = order,
						y = gba_mean)) +
		ggplot2::geom_jitter(
				data = ablation_scores,
				mapping = ggplot2::aes(
						x = order,
						y = gba_mean),
				height = 0,
				width = .1,
				color = "grey20") +
		ggplot2::ggtitle(
				label = "Co-Expression GBA by study subset size") +
		ggplot2::scale_x_continuous(
				"Number of Studies",
				breaks = c(1, 5, 10, 15, 18)) +
		ggplot2::scale_y_continuous(
				"Go Prediction AUROC",
				limits = c(.625, .8))

ggplot2::ggsave(
		"product/figures/go_pred_coexp_study_ablation_20201121.pdf",
		useDingbats=FALSE,
		width = 5,
		height = 4)

ggplot2::ggsave(
		"product/figures/go_pred_coexp_study_ablation_20201121.png",
		width = 5,
		height = 4)
