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

source("scripts/uniprot_to_cgd_id.R")
load("intermediate_data/cgd_id_2_uniprot_accn.Rdata")
load("intermediate_data/chromosome_features.Rdata")


gdrive_path <- "~/Collaborations/Candida Functional Genomics/HSP90 Physical Interactors/Datasets/SILAC dataset/HSP90 Depletion"
raw_path <-  "raw_data/silac/HSP90_depletion"

dir.create(raw_path, recursive=TRUE, showWarnings=FALSE)

####  load raw silac data from googledrive and store it in the raw_data directory ####

googledrive::drive_download(
	file=paste0(gdrive_path, "/Calbicans_new_fasta_corrected_labels_fc_pval.cdt"),
	path=paste0(raw_path, "/Calbicans_new_fasta_corrected_labels_fc_pval.cdt")


###########################################
# Download the raw intensities and compute fold changes from there
###########################################

googledrive::drive_download(
	file=paste0(gdrive_path, "/Calbicans_Fractions_Intensities _and_FC.xlsx"),
	path=paste0(raw_path, "/Calbicans_Fractions_Intensities _and_FC.xlsx"))

ca_silac_hsp90_intensities <- readxl::read_xlsx(
	path=paste0(raw_path, "/Calbicans_Fractions_Intensities _and_FC.xlsx"),
	sheet="Intensities") %>%
	dplyr::rename(
		uniprot_accn=Protein.IDs,
		description=Desription) %>%
	tidyr::separate(
		col=description,
		into=c("gene_name", "uniprot_entry", "uniprot_accn2", "description"),
		sep="[|][|]") %>%
	dplyr::select(-uniprot_accn2) %>%
	map_ca_uniprot_accn_to_feature_name() %>%
	dplyr::mutate(
		gene_name = ifelse(gene_name == "---", NA, gene_name)) %>%
	tidyr::gather(
		key="fraction_id",
		value="intensity",
		-uniprot_accn,
		-uniprot_entry,
		-gene_name,
		-description,
		-feature_name) %>%
	tidyr::separate(
		col="fraction_id",
		into=c("unit", "label", "fraction_index"),
		sep="[.]") %>%
	dplyr::select(-unit) %>%
	dplyr::mutate(
		fraction_index =  fraction_index %>% as.numeric(),
		condition = dplyr::case_when(
			label == "L" ~ "pharmacological",
			label == "M" ~ "wildtype",
			label == "H" ~ "genetic")) %>%
	dplyr::left_join(
		chromosome_features %>%
			dplyr::select(
				feature_name,
				ca_gene_name=gene_name,
				sac_ortholog),
			by="feature_name") %>%
	dplyr::mutate(
		gene_label =
			ifelse( !is.na(ca_gene_name), ca_gene_name,
			ifelse( !is.na(sac_ortholog), sac_ortholog,
			ifelse( !is.na(feature_name), feature_name,
			uniprot_accn)))) %>%
	dplyr::mutate(
		gene_label = ifelse(
			gene_label %in% c("GLY1", "GRE2", "PRB1", "UGA2"),
			feature_name, gene_label)) %>%
	dplyr::mutate(
		feature_name = ifelse(uniprot_accn == "A0A1D8PPN6", "C6_01700W_A", feature_name),
		gene_label = ifelse(uniprot_accn == "A0A1D8PPN6", "RPL32", gene_label)) %>%
	dplyr::mutate(
		feature_name = ifelse(uniprot_accn == "A0A1D8PRS1", "CR_00690C_A", feature_name),
		gene_label = ifelse(uniprot_accn == "A0A1D8PRS1", "GTB1", gene_label)) %>%
	dplyr::mutate(
		feature_name = ifelse(uniprot_accn == "A0A1D8PS38", "CR_01970C_A", feature_name),
		gene_label = ifelse(uniprot_accn == "A0A1D8PS38", "VMA4", gene_label)) %>%
	dplyr::mutate(
		feature_name = ifelse(uniprot_accn == "Q59NX9", "C5_03640W_A", feature_name),
		gene_label = ifelse(uniprot_accn == "Q59NX9", "DPH5", gene_label)) %>%
	dplyr::mutate(
		feature_name = ifelse(uniprot_accn == "Q59VN4", "C1_04240C_A", feature_name),
		gene_label = ifelse(uniprot_accn == "Q59VN4", "HHF1", gene_label)) %>%
	dplyr::mutate(
		feature_name = ifelse(uniprot_accn == "Q59XW4", "C2_10350C_A", feature_name),
		gene_label = ifelse(uniprot_accn == "Q59XW4", "ACS1", gene_label)) %>%
	dplyr::mutate(
		feature_name = ifelse(uniprot_accn == "Q5AB15", "C1_06760C_A", feature_name),
		gene_label = ifelse(uniprot_accn == "Q5AB15", "IPI3", gene_label)) %>%
	dplyr::mutate(
		feature_name = ifelse(uniprot_accn == "Q5AHF9", "C2_03870W_A", feature_name),
		gene_label = ifelse(uniprot_accn == "Q5AHF9", "GNA1", gene_label)) %>%
	dplyr::mutate(
		feature_name = ifelse(uniprot_accn == "Q9Y7F0", "C3_06180C_A", feature_name),
		gene_label = ifelse(uniprot_accn == "Q9Y7F0", "TSA1", gene_label)) %>%
	dplyr::mutate(
		feature_name = ifelse(uniprot_accn == "A0A1D8PNC7", "C5_02080C_A", feature_name),
		gene_label = ifelse(uniprot_accn == "A0A1D8PNC7", "HSP12", gene_label)) %>%
	dplyr::mutate(
		feature_name = ifelse(uniprot_accn == "A0A1D8PCG7", "C1_01370C_A", feature_name),
		gene_label = ifelse(uniprot_accn == "A0A1D8PCG7", "RPS21B", gene_label)) %>%
	dplyr::mutate(
		feature_name = ifelse(uniprot_accn == "A0A1D8PEN1", "C1_09530W_A", feature_name),
		gene_label = ifelse(uniprot_accn == "A0A1D8PEN1", "RAM2", gene_label)) %>%
	dplyr::mutate(
		feature_name = ifelse(uniprot_accn == "A0A1D8PFG4", "C1_12390C_A", feature_name),
		gene_label = ifelse(uniprot_accn == "A0A1D8PFG4", "RPL27A", gene_label)) %>%
	dplyr::mutate(
		feature_name = ifelse(uniprot_accn == "A0A1D8PFK5", "C1_13110C_A", feature_name),
		gene_label = ifelse(uniprot_accn == "A0A1D8PFK5", "CHS3", gene_label)) %>%
	dplyr::mutate(
		feature_name = ifelse(uniprot_accn == "A0A1D8PLJ3", "C4_02320C_A", feature_name),
		gene_label = ifelse(uniprot_accn == "A0A1D8PLJ3", "SOD1", gene_label)) %>%
	dplyr::mutate(
		feature_name = ifelse(uniprot_accn == "A0A1D8PMD0", "C4_05640C_A", feature_name),
		gene_label = ifelse(uniprot_accn == "A0A1D8PMD0", "HIS6", gene_label)) %>%
	dplyr::mutate(
		feature_name = ifelse(uniprot_accn == "A0A1D8PME6", "C4_05760W_A", feature_name),
		gene_label = ifelse(uniprot_accn == "A0A1D8PME6", "CGT1", gene_label)) %>%
	dplyr::mutate(
		feature_name = ifelse(uniprot_accn == "A0A1D8PPD4", "C6_00620W_A", feature_name),
		gene_label = ifelse(uniprot_accn == "A0A1D8PPD4", "FCA1", gene_label)) %>%
	dplyr::mutate(
		gene_label = ifelse(feature_name == "C7_03860W_A", "PBR1", gene_label)) %>%
	dplyr::mutate(
		gene_label = ifelse(feature_name == "C4_02340W_A", "PBI2", gene_label)) %>%
	dplyr::select(-gene_name, -ca_gene_name, -sac_ortholog)
save(ca_silac_hsp90_intensities, file="intermediate_data/ca_silac_hsp90_intensities.Rdata")


ca_silac_hsp90_fold_change <- ca_silac_hsp90_intensities %>%
	dplyr::select(-label) %>%
	tidyr::spread(condition, intensity) %>%
	dplyr::mutate(
		pharmacological_log_fc = log(pharmacological) - log(wildtype),
		pharmacological_log_fc = ifelse(is.finite(pharmacological_log_fc), pharmacological_log_fc, 0),
		genetic_log_fc = log(genetic) - log(wildtype),
		genetic_log_fc = ifelse(is.finite(genetic_log_fc), genetic_log_fc, 0))
save(ca_silac_hsp90_fold_change, file="intermediate_data/ca_silac_hsp90_fold_change.Rdata")

###############33
#### import silac data matrix and do some preliminary cleanup ###
#silac_matrix <- readr::read_tsv(
#	file=paste0(raw_path, "/Calbicans_new_fasta_corrected_labels_fc_pval.cdt",
#	col_types=readr::cols(
#	  .default = readr::col_double(),
#  	GID = readr::col_character(),
#  	GID_1 = readr::col_character(),
#		NAME = readr::col_character())) %>%
#	dplyr::rename(
#		protein_index = GID,
#		uniprot_accn = GID_1,
#		name = NAME) %>%
#
#	#get rid of blank columns
#	dplyr::select(-GWEIGHT, -Blank) %>%
#
#	#EWEIGHT is always 1
#	dplyr::filter(protein_index != "EWEIGHT") %>%
#
#	# replace "GENE<protein_index>X" with <protein_index>
#	dplyr::mutate(
#		protein_index = protein_index %>%
#			stringr::str_extract("[0-9]+") %>%
#			as.numeric) %>%
#
#	# separate out the annotations for the protein in to separate columns
#	tidyr::separate(
#		col=name,
#		into=c("name", "localization", "description", "go_annotations"),
#		sep="[|][|]") %>%
#	dplyr::mutate(
#			name = ifelse(name == "---", NA, name),
#		localization = ifelse(localization == "---", NA, localization),
#		description = ifelse(description == "---", NA, description),
#		go_annotations = ifelse(go_annotations == "---", NA, go_annotations))
#
#silac_matrix <- silac_matrix %>% map_ca_uniprot_accn_to_feature_name()
#
## tidy silac data
#silac <- silac_matrix %>%
#	tidyr::gather(
#		key="fraction",
#		value="fold_change",
#		-protein_index,
#		-uniprot_accn,
#		-feature_name,
#		-name,
#		-localization,
#		-description,
#		-go_annotations) %>%
#
#	# cleanup the fraction names into [perturbation, run_index, depletion, label]
#	tidyr::separate(
#		col="fraction",
#		into=c("value", "pertubation", "fraction_index"),
#		sep="[.]") %>%
#	dplyr::select(-value) %>%  # this is just "FC" but we're already calling the value column in spread "fold_change"
#	dplyr::mutate(fraction_index = as.numeric(fraction_index)) %>%
#	tidyr::separate(
#		col="pertubation",
#		into=c("depletion", "label"),
#		sep="_") %>%
#	dplyr::mutate(
#		depletion = dplyr::case_when(
#			depletion == "GeneticallyDepleted" ~ "genetic",
#			depletion == "PharmacologicallyDepleted" ~ "pharmacological",
#			TRUE ~ NA_character_)) %>%
#	dplyr::select(
#		protein_index,
#		uniprot_accn,
#		feature_name,
#		name,
#		description,
#		localization,
#		go_annotations,
#		depletion,
#		label,
#		fraction_index,
#		fold_change)
#
#save(silac, file="intermediate_data/silac.Rdata")
#
