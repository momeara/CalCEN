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
googledrive::drive_auth("~/.R/googledrive_auth.rds")


source("scripts/uniprot_to_cgd_id.R")

gdrive_path <- "~/Collaborations/Candida Functional Genomics/HSP90 Physical Interactors/Datasets/SILAC dataset/Ca-Sc-Sp"
raw_path <- "raw_data/silac/Ca-Sc_sp"
dir.create(raw_path, recursive=TRUE, showWarnings=FALSE)

####  load raw silac data from googledrive and store it in the raw_data directory ####

googledrive::drive_download(
	file=paste0(gdrive_path, "/raw_data/Yeast_IEF_IEX_CA_msb_prot_count_uniqpeps_FDR0010_no_singletons.txt"),
	path=paste0(raw_path, "/Yeast_IEF_IEX_CA_msb_prot_count_uniqpeps_FDR0010_no_singletons.txt"))

ca_silac_profile <- readr::read_tsv("raw_data/silac/Ca-Sc-Sp/Yeast_IEF_IEX_CA_msb_prot_count_uniqpeps_FDR0010_no_singletons.txt") %>%
	dplyr::rename(uniprot_accn = ProteinID)

ca_silac_profile <- ca_silac_profile %>% map_ca_uniprot_accn_to_feature_name()

##### complexes predicted by EPIC using  SILAC co-fractation, seeded with Sac complexes and YeastNet ####

googledrive::drive_download(
	file=paste0(gdrive_path, "/intermediate_data/Out.rf.comb.clust.xlsx"),
	path=paste0(raw_path, "/Out.rf.comb.clust.xlsx"))

ca_silac_clusters <- readxl::read_xlsx(
	path=paste0(raw_path, "/Out.rf.comb.clust.xlsx"),
	sheet="With Uniprot IDS") %>%
	dplyr::transmute(
		complex_index = `Complex ID` %>%
			stringr::str_extract("[0-9]+$") %>%
			as.numeric(),
		uniprot_accn = Members %>% stringr::str_split(", ")) %>%
	tidyr::unnest(uniprot_accn) %>%
	map_ca_uniprot_accn_to_feature_name() %>%
	dplyr::mutate(
		feature_name = dplyr::case_when(
			uniprot_accn == "Q9Y7F0"     ~ "C3_06180C_A", # TSA1
			uniprot_accn == "A0A1D8PPN6" ~ "C6_01700W_A", # RPL32
			uniprot_accn == "A0A1D8PFG4" ~ "C1_12390C_A", # RPL27A
			uniprot_accn == "A0A1D8PH21" ~ "C2_03960W_A", # RPL39
			uniprot_accn == "A0A1D8PCW6" ~ "C1_03030W_A", # RSP16A
			uniprot_accn == "A0A1D8PES3" ~ "C1_10100C_A", # NTF2
			uniprot_accn == "A0A1D8PLJ3" ~ "C4_02320C_A", # SOD1
			uniprot_accn == "A0A1D8PNC7" ~ "C5_02080C_A", # HSP12
			uniprot_accn == "Q59NX9"     ~ "C5_03640W_A", # DPH5
			TRUE ~ feature_name))
save(ca_silac_clusters, file="intermediate_data/ca_silac_clusters.Rdata")
ca_silac_clusters %>% readr::write_tsv("product/ca_silac_clusters_180712.tsv")
