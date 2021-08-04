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

gdrive_path <- "~/Collaborations/Candida Functional Genomics/HSP90 Physical Interactors/Datasets/AP-MS datasets"
raw_path <- "raw_data/APMS"
dir.create(raw_path, recursive=TRUE, showWarnings=FALSE)


### Download AP-MS datasets

# HSP90 and Co-chaperones as Bait
googledrive::drive_download(
	file=paste0(gdrive_path, "/list_tmp_TAP.xlsx"),
	path=paste0(raw_path, "/list_tmp_TAP.xlsx"))

# HSP90 as Bait under drug/no-drug conditions
googledrive::drive_download(
	file=paste0(gdrive_path, "/list_tmp_GFP.xlsx"),
	path=paste0(raw_path, "/list_tmp_GFP.xlsx"))

# HSP90
googledrive::drive_download(
	file=paste0(gdrive_path, "/bait-prey_cytoscape_TAP.txt"),
	path=paste0(raw_path, "/bait-prey_cytoscape_TAP.txt"))


# Process TAP Tag dataset
apms_tap <- readxl::read_xlsx(
	path=paste0(raw_path, "/list_tmp_TAP.xlsx")) %>%
	dplyr::mutate(
		bait_gene = dplyr::case_when(
				Bait == "E36A_HSP90" ~ "HSP90",
				Bait == "Sgt1" ~ "SGT1",
				TRUE ~ Bait),
		prey_gene = PreyGene %>% stringr::str_replace("\"","")) %>%
	dplyr::mutate(
		hit = `include as hit` == 'yes') %>%
	dplyr::left_join(
		chromosome_features %>%
			dplyr::filter(feature_name %>% stringr::str_detect("A$")) %>%
			dplyr::select(
				bait_feature_name = feature_name,
				bait_gene = gene_name,
				bait_feature_status=feature_status,
				bait_feature_type=feature_type),
			by=c("bait_gene")) %>%
	dplyr::mutate(
		bait_feature_name = dplyr::case_when(
			Bait == "CPR7" ~ "C3_00950C_A",
			TRUE ~ bait_feature_name))

apms_tap <- apms_tap %>% prey_id_map()
save(apms_tap, file="intermediate_data/apms_tap.Rdata")


###########################
# Process GFP tag dataset
apms_gfp <- readxl::read_xlsx(
	path=paste0(raw_path, "/list_tmp_GFP.xlsx")) %>%
	dplyr::mutate(
		bait_gene = "HSP90",
		bait_feature_name = "C7_02030W_A",
		condition = Bait %>% stringr::str_replace("Candida_HSP90_", ""),
		prey_gene = PreyGene,
		hit = include == "yes")

apms_gfp <- apms_gfp %>%
		prey_id_map()
save(apms_gfp, file="intermediate_data/apms_gfp.Rdata")


###################3
# Process TAP Tag datset for Cyctoscape

apms_tap_cytoscape <- readr::read_tsv(
	paste0(raw_path, "/bait-prey_cytoscape_TAP.txt"),
	col_types=readr::cols(
		Bait = readr::col_character(),
		Prey = readr::col_character(),
		AvgSpec = readr::col_double(),
		BFDR = readr::col_double())) %>%
	dplyr::mutate(
		bait_gene = dplyr::case_when(
				Bait == "Sgt1" ~ "SGT1",
				TRUE ~ Bait),
		prey_gene = Prey) %>%
	dplyr::arrange(dplyr::desc(AvgSpec)) %>%
	dplyr::distinct(bait_gene, prey_gene, .keep_all=TRUE) %>%
	dplyr::left_join(
		chromosome_features %>%
			dplyr::filter(feature_name %>% stringr::str_detect("A$")) %>%
			dplyr::select(
				bait_feature_name = feature_name,
				bait_gene = gene_name,
				bait_feature_status=feature_status,
				bait_feature_type=feature_type),
			by=c("bait_gene")) %>%
	dplyr::mutate(
		bait_feature_name = dplyr::case_when(
			Bait == "CPR7" ~ "C3_00950C_A",
			TRUE ~ bait_feature_name))

apms_tap_cytoscape <- apms_tap_cytoscape %>%
	prey_id_map()

save(apms_tap_cytoscape, file="intermediate_data/apms_tap_cytoscape.Rdata")

