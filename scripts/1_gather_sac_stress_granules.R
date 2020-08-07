# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:

library(plyr)
library(dplyr)
library(stringr)
library(readr)
library(tidyr)
library(magrittr)
library(googledrive)
googledrive::drive_auth("~/.R/googledrive_auth.rds")

library(googlesheets)
googlesheets::gs_auth("~/.R/googledrive_auth.rds")

gdrive_path <- "~/Collaborations/Candida Functional Genomics/HSP90 Physical Interactors/Case Studies/Stress Granules"
raw_path <- "raw_data/stress_granules"
dir.create(raw_path, recursive=TRUE, showWarnings=FALSE)


# retrieved from SGD 18/08/22
googledrive::drive_download(
	file=paste0(gdrive_path, "/cytoplasmic_stress_granule_annotations.txt"),
	path=paste0(raw_path, "/cytoplasmic_stress_granule_annotations.txt"))

# retrieved from SGD 18/08/22
googledrive::drive_download(
	file=paste0(gdrive_path, "/Pbody_annotations.txt"),
	path=paste0(raw_path, "/Pbody_annotations.txt"))

sac_stress_granule_annotations <- readr::read_tsv(
	file=paste0(raw_path, "/cytoplasmic_stress_granule_annotations.txt"),
	skip=7,
	col_types=readr::cols(
		Gene = readr::col_character(),
		`Gene Systematic Name` = readr::col_character(),
		`Gene Ontology Term` = readr::col_character(),
		`Gene Ontology Term ID` = readr::col_character(),
		Qualifier = readr::col_character(),
		Aspect = readr::col_character(),
		Evidence = readr::col_character(),
		Method = readr::col_character(),
		Source = readr::col_character(),
		`Assigned On` = readr::col_date(format = ""),
		`Annotation Extension` = readr::col_character(),
		Reference = readr::col_character())) %>%
	dplyr::rename(
		sac_ortholog=Gene,
		sac_feature_name=`Gene Systematic Name`,
		go_term=`Gene Ontology Term`,
		go_id=`Gene Ontology Term ID`,
		qualifier=Qualifier,
		aspect=Aspect,
		evidence=Evidence,
		method=Method,
		source=Source,
		assigned_on=`Assigned On`,
		annotation_extention=`Annotation Extension`,
		reference=Reference) %>%
	dplyr::left_join(
		chromosome_features %>%
			dplyr::select(
				feature_name,
				gene_name,
				sac_ortholog,
				feature_type,
				feature_status,
				description),
		by=c("sac_ortholog"))
save(sac_stress_granule_annotations, file="intermediate_data/sac_stress_granule_annotations.Rdata")

sac_pbody_annotations <- readr::read_tsv(
	file=paste0(raw_path, "/Pbody_annotations.txt"),
	skip=7,
	col_types=readr::cols(
		Gene = readr::col_character(),
		`Gene Systematic Name` = readr::col_character(),
		`Gene Ontology Term` = readr::col_character(),
		`Gene Ontology Term ID` = readr::col_character(),
		Qualifier = readr::col_character(),
		Aspect = readr::col_character(),
		Evidence = readr::col_character(),
		Method = readr::col_character(),
		Source = readr::col_character(),
		`Assigned On` = readr::col_date(format = ""),
		`Annotation Extension` = readr::col_character(),
		Reference = readr::col_character())) %>%
	dplyr::rename(
		sac_ortholog=Gene,
		sac_feature_name=`Gene Systematic Name`,
		go_term=`Gene Ontology Term`,
		go_id=`Gene Ontology Term ID`,
		qualifier=Qualifier,
		aspect=Aspect,
		evidence=Evidence,
		method=Method,
		source=Source,
		assigned_on=`Assigned On`,
		annotation_extention=`Annotation Extension`,
		reference=Reference) %>%
	dplyr::left_join(
		chromosome_features %>%
			dplyr::select(
				feature_name,
				gene_name,
				sac_ortholog,
				feature_type,
				feature_status,
				description),
		by=c("sac_ortholog"))
save(sac_pbody_annotations, file="intermediate_data/sac_pbody_annotations.Rdata")



