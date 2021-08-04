# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:


library(SRAdb)
library(plyr)
library(dplyr)
library(stringr)
library(magrittr)
source("parameters.R")

ncbi_taxon <- "5476"

# this is the Sequence Read Archive metadata database
# if it doesn't exist, download it
sra_fname <- paste0(scratch_dir, "/sra_meta/SRAmetadb.sqlite")
if(!file.exists(sra_fname)){
	SRAdb::getSRAdbFile(
			destdir = paste0(scratch_dir, "/sra_meta"),
			destfile = "SRAmetadb.sqlite.gz")
}

sra_con <- DBI::dbConnect(RSQLite::SQLite(), sra_fname)

# identify all RNA-Seq datasets in the relevant species
runs <- sra_con %>%
	dplyr::tbl("sra") %>%
	dplyr::filter(
		taxon_id == ncbi_taxon,
		library_strategy == "RNA-Seq") %>%
	dplyr::collect(n = Inf)

# identify all RNA-seq datasets that have no taxon_id but
# have the species name in the abstract. This may capture
# for example co-culture experiments 
runs2 <- sra_con %>%
	dplyr::tbl("sra") %>%
	dplyr::filter(
		is.na(taxon_id)) %>%
	dplyr::count(library_strategy) %>%
	dplyr::collect(n = Inf)

runs_no_species <- sra_con %>%
	dplyr::tbl("sra") %>%
	dplyr::filter(
#		is.na(taxon_id),
		library_strategy == "RNA-Seq") %>%
	dplyr::collect(n = Inf)

studies2 <- runs_no_species %>%
	dplyr::filter(
		study_abstract %>% stringr::str_detect("[aA]lbicans")) %>%
		dplyr::group_by(study_accession) %>%
		dplyr::distinct(sample_accession) %>%
		dplyr::summarize(n_samples = n())

studies3 <- runs_no_species %>%
	dplyr::filter(
		study_abstract %>% stringr::str_detect("[aA]lbicans")) %>%
	dplyr::group_by(study_accession) %>%
	dplyr::distinct(sample_accession) %>%
	dplyr::summarize(n_samples = n()) %>%
	dplyr::anti_join(studies2, by="study_accession")

studies4 <- runs_no_species %>%
	dplyr::filter(
		study_title %>% stringr::str_detect("[aA]lbicans")) %>%
	dplyr::group_by(study_accession) %>%
	dplyr::distinct(sample_accession) %>%
	dplyr::summarize(n_samples = n()) %>%
	dplyr::anti_join(studies2, by="study_accession") %>%
	dplyr::anti_join(studies3, by="study_accession")

studies5 <- runs_no_species %>%
	dplyr::filter(
		experiment_title %>% stringr::str_detect("[aA]lbicans")) %>%
	dplyr::group_by(study_accession) %>%
	dplyr::distinct(sample_accession) %>%
	dplyr::summarize(n_samples = n()) %>%
	dplyr::anti_join(studies2, by="study_accession") %>%
	dplyr::anti_join(studies3, by="study_accession") %>%
	dplyr::anti_join(studies4, by="study_accession")


studies_a <- runs_no_species <- sra_con %>%
	dplyr::tbl("sra") %>%
	dplyr::filter(
		!is.na(taxon_id),
		taxon_id != ncbi_taxon,
		library_strategy == "RNA-Seq") %>%
	dplyr::collect(n = Inf) %>%
	dplyr::filter(
		experiment_title %>% stringr::str_detect("[aA]lbicans") |
		study_title %>% stringr::str_detect("[aA]lbicans") |
    study_abstract %>% stringr::str_detect("[aA]lbicans")) %>%
	dplyr::group_by(study_accession) %>%
	dplyr::distinct(sample_accession) %>%
	dplyr::summarize(n_samples = n())








######################################################################
# filter down to the list of studies/runs to include in the analysis #
######################################################################

studies_to_exclude <- data.frame(
		study_accession = c(
				"ERP014953"))              # failed to parse SRA files

runs_to_exclude <- data.frame(
		run_accession = c(
				"SRR604746",      # only 742 (0.14%) aligned exactly 1 time C. albicans SC5314 transcriptome
				"SRR772101"))     # only 191 reads total

runs <- runs %>%
		dplyr::anti_join(studies_to_exclude, by = "study_accession") %>%
		dplyr::anti_join(runs_to_exclude, by = "run_accession")


studies <- runs %>%
		dplyr::group_by(study_accession) %>%
		dplyr::distinct(sample_accession) %>%
		dplyr::summarize(n_samples = n())

CalCEN_datasets <- readr::read_tsv("raw_data/CalCEN_datasets.tsv")

sra_con %>%
		dplyr::copy_to(CalCEN_datasets)

ca_runs <- sra_con %>%
		dplyr::tbl("sra") %>%
		dplyr::semi_join(
				sra_con %>% dplyr::tbl("CalCEN_datasets"),
				by = "study_accession") %>%
		dplyr::collect(n = Inf)


# study_accession, ref, strain background, perturbations, runs per sample
ca_runs <- ca_runs %>%
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
		sra_fname = paste0(scratch_dir, "/sra/sra/", run_accession, ".sra")) %>%
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

ca_runs <- ca_runs %>%
		dplyr::filter(
				run_accession != "SRR1251735",
				run_accession != "SRR1251736")

# duplicate runs
ca_runs %>%
		dplyr::semi_join(
				ca_runs %>%
						dplyr::count(run_accession) %>%
						dplyr::filter(n > 1))


#ca_runs %>%
#	readr::write_tsv("product/ca_runs_180608.tsv")

# this is the raw data for studies table in the manuscript
ca_runs %>%
	dplyr::count(study_accession, updated_date, is_paired, background_strain) %>%
	dplyr::arrange(updated_date %>% as.Date(format="%Y-%m-%d")) %>%
	dplyr::rename(n_runs=n) %>%
	readr::write_tsv("product/ca_studies_180608.tsv")


ca_runs %>%
	readr::write_tsv("product/ca_runs_200928.tsv")



