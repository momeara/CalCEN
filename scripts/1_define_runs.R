# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:


library(SRAdb)
library(plyr)
library(dplyr)
library(stringr)
library(magrittr)
source("scripts/parameters.R")

sra_fname <- paste0(base_dir, "/sra_meta/SRAmetadb.sqlite")
if(!file.exists(sra_fname)){
	SRAdb::getSRAdbFile(destfile = sra_fname)
}

sra_con <- DBI::dbConnect(RSQLite::SQLite(), sra_fname)

candida_albicans_ncbi_taxon <- "5476"
ca_runs <- sra_con %>%
	dplyr::tbl("sra") %>%
	dplyr::filter(
		taxon_id == candida_albicans_ncbi_taxon,
		library_strategy == "RNA-Seq") %>%
	dplyr::collect()

# study_accession, ref, strain background, perturbations, runs per sample
ca_runs <- ca_runs %>%
#	dplyr::semi_join(
#		ca_runs %>%
#			dplyr::group_by(study_accession) %>%
#			dplyr::distinct(sample_accession) %>%
#			dplyr::summarize(n_samples = n()) %>%
#			dplyr::filter(n_samples >= 20),
#		by="study_accession") %>%
#	dplyr::filter(
#		study_accession != "SRP058281") %>%
	dplyr::filter(
		study_accession != "ERP014953",    # failed to parse SRA files
		run_accession != "SRR604746",      # only 742 (0.14%) aligned exactly 1 time C. albicans SC5314 transcriptome
		run_accession != "SRR772101") %>%  # only 191 reads total
	dplyr::mutate(
		is_paired = library_layout %>% stringr::str_detect("^PAIRED"),
		background_strain = sample_attribute %>%
			stringr::str_match("strain: ([^ ]+)") %>%
			unlist() %>%
			magrittr::extract(,2),
		background_strain = dplyr::case_when(
			study_accession == "SRP064163" ~ "SC5314",
			study_accession == "ERP013259" ~ "SN95",
			study_accession == "SRP080770" ~ "SC5314",
			study_accession == "SRP108757" ~ "SC5314",
			study_accession == "ERP003535" ~ "CAI4",
			study_accession == "ERP006610" ~ "CAI4",
			TRUE ~ background_strain)) %>%
	dplyr::mutate(
		sra_fname = paste0(base_dir, "/sra_80601/", run_accession, ".sra")) %>%
	dplyr::select(
		submission_accession,
		submission_lab,
		updated_date,
		sradb_updated,

		study_accession,
		study_alias,
		study_name,
		study_title,
		study_type,
		study_abstract,
		center_project_name,
		study_description,
		study_attribute,
		description,
		design_description,

		sample_accession,
		sample_alias,
		sample_name,
		sample_attribute,
		taxon_id,
		background_strain,

		platform,
		platform_parameters,
		instrument_model,
		library_name,
		library_strategy,
		library_source,
		library_selection,
		library_layout,
		library_construction_protocol,
		is_paired,
		read_spec,

		experiment_accession,
		experiment_alias,
		experiment_name,
		experiment_title,
		experiment_attribute,

		run_accession,
		run_alias,
		run_center,
		spots,
		bases,
		sra_fname)

ca_runs %>%
	readr::write_tsv("product/ca_runs_180608.tsv")


# this is the raw data for studies table in the manuscript
ca_runs %>%
	dplyr::count(study_accession, updated_date, is_paired, background_strain) %>%
	dplyr::arrange(updated_date %>% as.Date(format="%Y-%m-%d")) %>%
	dplyr::rename(n_runs=n) %>%
	readr::write_tsv("product/ca_studies_180608.tsv")

