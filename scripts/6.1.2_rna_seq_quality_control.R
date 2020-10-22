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




# get the minimum fraction of reads aligned exactly 1 time for each study
z <- ca_runs %>%
	dplyr::mutate(layout = ifelse(is_paired, "Paired", "Single")) %>%
	dplyr::select(study_accession, run_accession, layout) %>%
	dplyr::left_join(estimated_expression_meta, by="run_accession") %>%
	dplyr::filter(n_reads %>% is.na %>% magrittr::not()) %>%
	dplyr::mutate(frac_exact_align = n_reads_aligned_1_time/n_reads) %>%
	dplyr::arrange(frac_exact_align) %>%
	dplyr::select(study_accession, run_accession, layout, frac_exact_align)

p <- ggplot2::ggplot() +
	ggplot2::theme_bw() +
	ggplot2::geom_dotplot(
		data = z,
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
	filename="product/figures/percent_mapped_once_by_study.pdf",
	height=5, width=5,
	useDingbats=FALSE)

ggplot2::ggsave(
	filename="product/figures/percent_mapped_once_by_study.png",
	height=5, width=5)


