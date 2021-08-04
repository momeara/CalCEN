# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:

library(plyr)
library(dplyr)
library(stringr)
library(readr)
library(tidyr)
library(magrittr)
library(readxl)
library(googledrive)
googledrive::drive_auth(email='mattjomeara@gmail.com')
load("intermediate_data/chromosome_features.Rdata")

source("scripts/apms_prey_id_map.R")


gdrive_path <- "~/Collaborations/Candida Functional Genomics/HSP90 Physical Interactors/Datasets/Sac Physical Interactions SGD 04-08-2019"
raw_path <- "raw_data/SGD_physical_190408"
dir.create(raw_path, recursive=TRUE, showWarnings=FALSE)


googledrive::drive_download(
	file=paste0(gdrive_path, "/AHA1_interactions.txt"),
  path=paste0(raw_path, "/AHA1_interactions.txt"))

googledrive::drive_download(
	file=paste0(gdrive_path, "/CDC37_interactions.txt"),
  path=paste0(raw_path, "/CDC37_interactions.txt"))

googledrive::drive_download(
	file=paste0(gdrive_path, "/CNS1_interactions.txt"),
  path=paste0(raw_path, "/CNS1_interactions.txt"))

googledrive::drive_download(
	file=paste0(gdrive_path, "/CPR6_interactions.txt"),
  path=paste0(raw_path, "/CPR6_interactions.txt"))

googledrive::drive_download(
	file=paste0(gdrive_path, "/CPR7_interactions.txt"),
  path=paste0(raw_path, "/CPR7_interactions.txt"))

googledrive::drive_download(
	file=paste0(gdrive_path, "/HCH1_interactions.txt"),
  path=paste0(raw_path, "/HCH1_interactions.txt"))

googledrive::drive_download(
	file=paste0(gdrive_path, "/HSC82_interactions.txt"),
  path=paste0(raw_path, "/HSC82_interactions.txt"))

googledrive::drive_download(
	file=paste0(gdrive_path, "/HSP82_interactions.txt"),
  path=paste0(raw_path, "/HSP82_interactions.txt"))

googledrive::drive_download(
	file=paste0(gdrive_path, "/SBA1_interactions.txt"),
  path=paste0(raw_path, "/SBA1_interactions.txt"))

googledrive::drive_download(
	file=paste0(gdrive_path, "/SGT1_interactions.txt"),
  path=paste0(raw_path, "/SGT1_interactions.txt"))

googledrive::drive_download(
	file=paste0(gdrive_path, "/STI1_interactions.txt"),
  path=paste0(raw_path, "/STI1_interactions.txt"))


parse_SGD_interactions_file <- function(fname){
	readr::read_tsv(
		file=fname,
		skip=9,
		col_names=c(
				"gene_symbol_1",          #Interactor
				"feature_name_1",     #Interactor Systematic Name
				"gene_symbol_2",          #Interactor
				"featur_name_2",      #Interactor Systematic Name
				"interaction_type",   #Type
				"experimental_system", #Assay
				"annotation",         #Annotation
				"action",             #Action
				"modification",       #Modification
				"phenotype",          #Phenotype
				"source",             #Source
				"reference",          #Reference
				"note"),              #Note
			col_types=readr::cols(
				gene_symbol_1 = readr::col_character(),
				feature_name_1 = readr::col_character(),
				gene_symbol_2 = readr::col_character(),
				featur_name_2 = readr::col_character(),
				interaction_type = readr::col_character(),
				experimental_system = readr::col_character(),
				annotation = readr::col_character(),
				action = readr::col_character(),
				modification = readr::col_character(),
				phenotype = readr::col_character(),
				source = readr::col_character(),
				reference = readr::col_character(),
				note = readr::col_character())) %>%
		dplyr::mutate(
			experimental_system_abbreviation = dplyr::case_when(
				experimental_system == "Dosage Lethality" ~ "DL",
				experimental_system == "Dosage Growth Defect" ~ "DGD",
				experimental_system == "Dosage Rescue" ~ "DR",
				experimental_system == "Synthetic Rescue" ~ "SR",
				experimental_system == "Phenotypic Suppression" ~ "PS",
				experimental_system == "Phenotypic Enhancement" ~ "PE",
				experimental_system == "Synthetic Lethality" ~ "SL",
				experimental_system == "Synthetic Growth Defect" ~ "SGD",
				experimental_system == "Positive Genetic" ~ "PG",
				experimental_system == "Negative Genetic" ~ "NG",
				experimental_system == "Synthetic Haploinsufficiency" ~ "HIP",
				experimental_system == "Affinity Capture-Luminescence" ~ "APL",
				experimental_system == "Proximity Label-MS" ~ "PL-MS",
				experimental_system == "Far Western" ~ "FW",
				experimental_system == "FRET" ~ "FRET",
				experimental_system == "Protein-RNA" ~ "P-RNA",
				experimental_system == "Co-localization" ~ "CoL",
				experimental_system == "Protein-peptide" ~ "Pro-Pep",
				experimental_system == "Co-crystal Structure" ~ "CoX",
				experimental_system == "Co-fractionation" ~ "CoF",
				experimental_system == "Co-purification" ~ "CoP",
				experimental_system == "Biochemical Activity" ~ "BA",
				experimental_system == "PCA" ~ "PCA",
				experimental_system == "Reconstituted Complex" ~ "ReC",
				experimental_system == "Two-hybrid" ~ "2H",
				experimental_system == "Affinity Capture-Western" ~ "AP-W",
				experimental_system == "Affinity Capture-RNA" ~ "AP-RNA",
				experimental_system == "Affinity Capture-MS" ~ "AP-MS",
				TRUE ~ experimental_system),
			annotation_abbreviation = dplyr::case_when(
				annotation == "high-throughput" ~ "ht",
				annotation == "manually curated" ~ "mc"))
}


sac_sge_interactions <- dplyr::bind_rows(
	parse_SGD_interactions_file(paste0(raw_path, "/AHA1_interactions.txt")),
	parse_SGD_interactions_file(paste0(raw_path, "/CDC37_interactions.txt")),
	parse_SGD_interactions_file(paste0(raw_path, "/CNS1_interactions.txt")),
	parse_SGD_interactions_file(paste0(raw_path, "/CPR6_interactions.txt")),
	parse_SGD_interactions_file(paste0(raw_path, "/CPR7_interactions.txt")),
	parse_SGD_interactions_file(paste0(raw_path, "/HCH1_interactions.txt")),
	parse_SGD_interactions_file(paste0(raw_path, "/HSC82_interactions.txt")),
	parse_SGD_interactions_file(paste0(raw_path, "/HSP82_interactions.txt")),
	parse_SGD_interactions_file(paste0(raw_path, "/SBA1_interactions.txt")),
	parse_SGD_interactions_file(paste0(raw_path, "/SGT1_interactions.txt")),
	parse_SGD_interactions_file(paste0(raw_path, "/STI1_interactions.txt"))) %>%
	dplyr::left_join(
		chromosome_features %>%
			dplyr::select(
				ca_feature_name_1 = feature_name, gene_symbol_1 = sac_ortholog),
		by="gene_symbol_1") %>%
	dplyr::left_join(
		chromosome_features %>%
			dplyr::select(
				ca_feature_name_2 = feature_name, gene_symbol_2 = sac_ortholog),
		by="gene_symbol_2") %>%
	dplyr::mutate(
		ca_feature_name_1 = ifelse(gene_symbol_1 == "HSP82", "C7_02030W_A", ca_feature_name_1),
		ca_feature_name_2 = ifelse(gene_symbol_2 == "HSP82", "C7_02030W_A", ca_feature_name_2))

save(sac_sge_interactions, file="intermediate_data/sac_sge_interactions.Rdata")
