# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:

library(plyr)
library(dplyr)
library(magrittr)
library(stringr)
library(readr)
library(ggplot2)
library(scales)

ca_runs <- dplyr::bind_rows(
		readr::read_tsv("product/ca_runs_180608.tsv") %>% dplyr::mutate(year_collected = 2018),
		readr::read_tsv("product/ca_runs_200928.tsv") %>% dplyr::mutate(year_collected = 2020))

load("intermediate_data/estimated_expression.Rdata")
load("intermediate_data/estimated_expression_meta.Rdata")

run_progress <- ca_runs %>%
		dplyr::count(year_collected, study_accession, name="sra_runs_per_study") %>%
		dplyr::filter(sra_runs_per_study >= 20) %>%
		dplyr::left_join(
				ca_runs %>% dplyr::distinct(run_accession, study_accession),
				by = c("study_accession")) %>%
		dplyr::left_join(
				estimated_expression %>%
				dplyr::distinct(run_accession) %>%
				dplyr::mutate(got_expression = TRUE),
				by = c("run_accession")) %>%
		dplyr::left_join(
				estimated_expression_meta %>%
				dplyr::distinct(run_accession) %>%
				dplyr::mutate(got_meta = TRUE),
				by = c("run_accession"))

# which runs do we have but for which we're missing the meta info
run_progress %>%
		dplyr::filter(got_expression, is.na(got_meta))




n_zero_genes <- estimated_expression %>%
		dplyr::group_by(study_accession, run_accession) %>%
		dplyr::summarize(
				n_nonzero_expression = sum(FPKM != 0),
				.groups = "drop")






ca_runs_final <- ca_runs %>%
		dplyr::semi_join(
				ca_runs %>%
				dplyr::count(study_accession) %>%
				dplyr::filter(n >= 20),
				by = "study_accession") %>%
		dplyr::semi_join(
				n_zero_genes %>%
				dplyr::filter(n_nonzero_expression >= 3113), # half the genes
				by = c("study_accession", "run_accession"))

ca_runs_final %>%
		save(file="intermediate_data/ca_runs_final.Rdata")
ca_runs_final %>%
		save(file="product/ca_runs_20201024.tsv")



###########################################################################
# get the minimum fraction of reads aligned exactly 1 time for each study #
###########################################################################
n_reads_aligned <- ca_runs %>%
		dplyr::mutate(layout = ifelse(is_paired, "Paired", "Single")) %>%
		dplyr::select(study_accession, run_accession, layout) %>%
		dplyr::left_join(
				estimated_expression_meta,
				by=c("study_accession", "run_accession")) %>%
		dplyr::filter(n_reads %>% is.na %>% magrittr::not()) %>%
		dplyr::mutate(frac_exact_align = n_reads_aligned_1_time/n_reads) %>%
		dplyr::arrange(frac_exact_align) %>%
		dplyr::select(
				study_accession,
				run_accession,
				layout,
				n_reads_aligned_1_time,
				frac_exact_align)

p <- ggplot2::ggplot() +
	ggplot2::theme_bw() +
	ggplot2::geom_dotplot(
		data = n_reads_aligned,
		mapping = ggplot2::aes(x=study_accession, y=frac_exact_align),
		binaxis = 'y', stackdir = 'center',
		stackratio = .4,
		dotsize = .3) +
	ggplot2::facet_wrap(~layout, nrow=2, scales="free_x") +
	ggplot2::scale_x_discrete(name="SRA Study Accession") +
	ggplot2::scale_y_continuous(
		name="Percent reads mapped eactly once to SC5314 Assembly 22",
		label=scales::percent)

ggplot2::ggsave(
	filename="product/figures/percent_mapped_once_by_study_20201024.pdf",
	height=5, width=5,
	useDingbats=FALSE)

ggplot2::ggsave(
	filename="product/figures/percent_mapped_once_by_study_20201024.png",
	height=5, width=5)


###################################################
# N zero expression genes vs percent reads mapped #
###################################################


n_zero_genes <- estimated_expression %>%
		dplyr::group_by(study_accession, run_accession) %>%
		dplyr::summarize(
				n_nonzero_expression = sum(FPKM != 0),
				.groups = "drop")

data <- n_zero_genes %>%
		dplyr::left_join(
				n_reads_aligned,
				by = c("study_accession", "run_accession"))


ggplot2::ggplot() +
		ggplot2::theme_bw() +
		ggplot2::geom_point(
				data = data,
				mapping = ggplot2::aes(
						x = n_reads_aligned_1_time + 1,
						y = n_nonzero_expression + 1),
				color = "grey50") +
		ggplot2::geom_point(
				data = data %>% dplyr::filter(n_nonzero_expression >= 3113),
				mapping = ggplot2::aes(
						x = n_reads_aligned_1_time + 1,
						y = n_nonzero_expression + 1)) +
		ggplot2::ggtitle(
				label = "Expression depth vs coverage by RNA-seq run") +
		ggplot2::scale_x_log10(
				"Number reads that map exactly once to a candida Albicans gene") +
		ggplot2::scale_y_log10(
				"Number of candida Alibcans genes that have non-zero expression")

ggplot2::ggsave(
		filename="product/figures/expression_depth_vs_coverage_by_study_20201024.pdf",
		height=8, width=8,
		useDingbats=FALSE)

ggplot2::ggsave(
		filename="product/figures/expression_depth_vs_coverage_by_study_20201024.png",
		height=8, width=8)

data %>%
		readr::write_tsv("product/figures/expression-depth_vs_coverage_by_study_20201024.tsv")
