# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:

library(plyr)
library(dplyr)
library(magrittr)
library(stringr)
library(readr)

load("intermediate_data/CalCEN.Rdata")

source("scripts/gene_report.R")

top_CalCEN <- CalCEN %>%
		dplyr::filter(feature_name_1 != feature_name_2) %>%
		dplyr::filter(score > .95) %>%
		dplyr::group_by(feature_name_1) %>%
		dplyr::arrange(desc(score)) %>%
		dplyr::filter(dplyr::row_number() <= 50) %>%
		dplyr::ungroup() %>%
		gene_gene_report()

CalCEN_sig_neighbors <- CalCEN %>%
		dplyr::filter(feature_name_1 != feature_name_2) %>%
		dplyr::filter(score > .98) %>%
		dplyr::group_by(feature_name_1) %>%
		dplyr::arrange(dplyr::desc(score)) %>%
		dplyr::slice(1:50) %>%
		dplyr::select(-score) %>%
		dplyr::ungroup() %>%
		gene_gene_report() %>%
		dplyr::select(
				-blastp_EValue,
				-genetic_interaction,
				-physical_interaction,
				-sac_genetic_interaction,
				-sac_physical_interaction) %>%
		dplyr::rename(CalCEN = coexp_score)

CalCEN_sig_neighbors %>%
		readr::write_tsv("product/CalCEN_sig_neighbors_20210121.tsv")


PGA52 <- CalCEN %>%
	dplyr::filter(feature_name_1 == "C2_00100C_A") %>%
	dplyr::filter(feature_name_2 != feature_name_1) %>%
	dplyr::arrange(desc(score)) %>%
	dplyr::select(-score) %>%
	gene_gene_report()

C7_00310C_A <- CalCEN %>%
	dplyr::filter(feature_name_1 == "C7_00310C_A") %>%
	dplyr::filter(feature_name_2 != feature_name_1) %>%
	dplyr::arrange(desc(score)) %>%
	dplyr::select(-score) %>%
	gene_gene_report() %>%
	dplyr::select(
		-description_1,
		-gene_name_1, -sac_ortholog_1, -blastp_EValue, -physical_interaction, -genetic_interaction, -sac_genetic_interaction, -sac_physical_interaction
	) %>%
	dplyr::mutate(description_2 = description_2 %>% stringr::str_sub(1,50)) %>%
	head(50)



FOX2 <- CalCEN %>%
	dplyr::filter(feature_name_1 == "C3_00810C_A") %>%
	dplyr::filter(feature_name_2 != feature_name_1) %>%
	dplyr::arrange(desc(score)) %>%
	dplyr::select(-score) %>%
	gene_gene_report() %>%
	dplyr::select(
		-description_1, -gene_name_1, -sac_ortholog_1,
#		-blastp_EValue, -physical_interaction, -genetic_interaction, -sac_genetic_interaction, -sac_physical_interaction
	) %>%
	dplyr::mutate(description_2 = description_2 %>% stringr::str_sub(1,80)) %>%
	head(50)


#
seeds <- readr::read_tsv("raw_data/deORFinizing_seeds_20201006.tsv")

# all present
missing_seeds <- seeds %>% dplyr::anti_join(
		chromosome_features,
		by = c("Seed" = "feature_name"))



# and in the top 30 associations to a seed
# and itself not a seed
seed_to_neighbor <- CalCEN %>%
		dplyr::filter(feature_name_1 != feature_name_2) %>%
		dplyr::semi_join(
				seeds,
				by = c("feature_name_1" = "Seed")) %>%
		dplyr::group_by(feature_name_1) %>%
		dplyr::arrange(desc(score)) %>%
		dplyr::slice(1:50) %>%
		dplyr::ungroup() %>%
		dplyr::transmute(
				seed = feature_name_1,
				source_type = "seed",
				target_type = "neighbor",
				feature_name_1,
				feature_name_2)

neighbor_to_neighbor <- seed_to_neighbor %>%
		dplyr::group_by(feature_name_1) %>%
		dplyr::do({
				data <- .
				tidyr::expand_grid(
						feature_name_1 = data$feature_name_2,
						feature_name_2 = data$feature_name_2) %>%
						dplyr::mutate(
								seed = data$feature_name_1[1])}) %>%
		dplyr::ungroup() %>%
		dplyr::filter(feature_name_1 < feature_name_2) %>%
		dplyr::transmute(
				seed,
				source_type = "neighbor",
				target_type = "neighbor",
				feature_name_1,
				feature_name_2)


interactions <- dplyr::bind_rows(
		seed_to_neighbor,
		neighbor_to_neighbor)

interactions <- interactions %>%
		dplyr::left_join(
				interactions %>%
						dplyr::distinct(feature_name_1, feature_name_2) %>%
						gene_gene_report() %>%
						dplyr::distinct(feature_name_1, feature_name_2, .keep_all=T),
				by=c("feature_name_1", "feature_name_2"))

interactions <- interactions %>%
		dplyr::filter(!(
				coexp_score < .9 &
				is.na(blastp_EValue) &
				is.na(genetic_interaction) &
				is.na(physical_interaction) &
				is.na(sac_genetic_interaction) &
				is.na(sac_physical_interaction)))

# the sac intractions are all na to remove them
interactions <- interactions %>%
		dplyr::select(
				-sac_genetic_interaction,
				-sac_physical_interaction)

interactions <- interactions %>%
		dplyr::mutate(
				source_label = ifelse(!is.na(gene_name_1), gene_name_1, feature_name_1),
				target_label = ifelse(!is.na(gene_name_2), gene_name_2, feature_name_2))

interactions %>%
		readr::write_tsv("product/deORFinizing_seed_to_neighborhood_20201007.tsv")


### seeds for kyla
CalCEN_for_kyla <- CalCEN %>%
	dplyr::semi_join(
		readr::read_tsv("raw_data/seeds_selvig_20201013.tsv"),
		by = c("feature_name_1" = "feature_name")) %>%
	dplyr::filter(feature_name_1 != feature_name_2) %>%
	dplyr::arrange(desc(score)) %>%
	dplyr::group_by(feature_name_1) %>%
	dplyr::slice(1:50) %>%
	dplyr::ungroup() %>%
	gene_gene_report()

CalCEN_for_kyla %>%
		readr::write_tsv("product/seeds_for_selvid_coexp_20201013.tsv")

##########
# ERG11

report <- CalCEN %>%
	dplyr::filter(feature_name_1 == "C5_00660C_A") %>%
	dplyr::filter(feature_name_2 != feature_name_1) %>%
	dplyr::arrange(desc(score)) %>%
	dplyr::slice(1:50) %>%
	dplyr::select(-score) %>%
	gene_gene_report()

report %>%
	readr::write_tsv("product/ERG11_gene_gene_report_20201015.tsv")

report <- CalCEN %>%
	dplyr::filter(feature_name_1 == "C3_06660C_A") %>%
	dplyr::filter(feature_name_2 != feature_name_1) %>%
	dplyr::arrange(desc(score)) %>%
	dplyr::slice(1:50) %>%
	dplyr::select(-score) %>%
	gene_gene_report()

report %>%
	readr::write_tsv("product/C3_06660C_A_gene_gene_report_20201027.tsv")


####################
#
grace_null_goterm <- CalCEN %>%
		dplyr::semi_join(
				readr::read_csv("raw_data/GRACE_NULLGOTERM.csv"),
				by = c("feature_name_1" = "Gene ID")) %>%
		dplyr::filter(feature_name_1 != feature_name_2) %>%
		dplyr::arrange(desc(score)) %>%
		dplyr::group_by(feature_name_1) %>%
		dplyr::slice(1:50) %>%
		dplyr::ungroup() %>%
		dplyr::select(-score) %>%
		gene_gene_report()

grace_null_goterm %>%
		readr::write_tsv("product/grace_null_goterm_gene_gene_report_20201101.tsv")


# hsp90
report <- CalCEN %>%
	dplyr::filter(feature_name_1 == "C7_02030W_A") %>%
	dplyr::filter(feature_name_2 != feature_name_1) %>%
	dplyr::arrange(desc(score)) %>%
	dplyr::slice(1:50) %>%
	dplyr::select(-score) %>%
	gene_gene_report()

report %>%
	readr::write_tsv("product/C7_02030W_A_HSP90_gene_gene_report_20201027.tsv")

## cell wall
grace_null_goterm <- CalCEN %>%
		dplyr::semi_join(
				readr::read_csv("raw_data/cellwall_20201105.csv"),
				by = c("feature_name_1" = "Gene ID")) %>%
		dplyr::filter(feature_name_1 != feature_name_2) %>%
		dplyr::arrange(desc(score)) %>%
		dplyr::group_by(feature_name_1) %>%
		dplyr::slice(1:50) %>%
		dplyr::ungroup() %>%
		dplyr::select(-score) %>%
		gene_gene_report()

grace_null_goterm %>%
		readr::write_tsv("product/cellwall_gene_gene_report_20201105.tsv")

# CR_06140W_A
report <- CalCEN %>%
	dplyr::filter(feature_name_1 == "CR_06140W_A") %>%
	dplyr::filter(feature_name_2 != feature_name_1) %>%
	dplyr::arrange(desc(score)) %>%
	dplyr::slice(1:50) %>%
	dplyr::select(-score) %>%
	gene_gene_report()

report %>%
	readr::write_tsv("product/CR_06140W_A_gene_gene_report_20201112.tsv")


# C1_14340C_A
report <- CalCEN %>%
	dplyr::filter(feature_name_1 == "C1_14340C_A") %>%
	dplyr::filter(feature_name_2 != feature_name_1) %>%
	dplyr::arrange(desc(score)) %>%
	dplyr::slice(1:50) %>%
	dplyr::select(-score) %>%
	gene_gene_report()

report %>%
	readr::write_tsv("product/C1_14340C_A_gene_gene_report_20201112.tsv")


# DnaJ from SAC:
sac_dnaj <- tibble::tibble(sac_ortholog = c(
		"YDJ1", "XDJ1", "APJ1", "SIS1", "DJP1",
		"ZUO1", "SWA2", "JJJ1", "JJJ2", "JJJ3",
		"CAJ1", "CWC23", "MDJ1", "MDJ2", "PAM18",
		"JAC1", "JID1", "SCJ1", "HLJ1", "JEM1",
		"SEC63", "ERJ5"))
candida_dnaj <- 	sac_dnaj %>%
		dplyr::left_join(
				chromosome_features,
				by = "sac_ortholog")

candida_dnaj %>%
	readr::write_tsv("product/dnaj_sac_orthologs_20201114.tsv")


report <- CalCEN %>%
		dplyr::semi_join(
		    candida_dnaj %>% dplyr::filter(!is.na(feature_name)),
				by = c("feature_name_1" = "feature_name")) %>%
		dplyr::filter(feature_name_1 != feature_name_2) %>%
		dplyr::arrange(desc(score)) %>%
		dplyr::group_by(feature_name_1) %>%
		dplyr::slice(1:50) %>%
		dplyr::ungroup() %>%
		dplyr::select(-score) %>%
		gene_gene_report()
report %>%
	readr::write_tsv("product/dnaj_gene_gene_report_20201114.tsv")


# C2_01370C_A
report <- CalCEN %>%
	dplyr::filter(feature_name_1 == "C2_01370C_A") %>%
	dplyr::filter(feature_name_2 != feature_name_1) %>%
	dplyr::arrange(desc(score)) %>%
	dplyr::slice(1:50) %>%
	dplyr::select(-score) %>%
	gene_gene_report()

report %>%
	readr::write_tsv("product/C2_01370C_A_gene_gene_report_20201122.tsv")


# C1_01800W_A
report <- CalCEN %>%
	dplyr::filter(feature_name_1 == "C1_01800W_A") %>%
	dplyr::filter(feature_name_2 != feature_name_1) %>%
	dplyr::arrange(desc(score)) %>%
	dplyr::slice(1:50) %>%
	dplyr::select(-score) %>%
	gene_gene_report()

report %>%
	readr::write_tsv("product/C1_01800W_A_gene_gene_report_20201122.tsv")


feature_name <- "C6_03200W_A"
report <- CalCEN %>%
	dplyr::filter(feature_name_1 == feature_name) %>%
	dplyr::filter(feature_name_2 != feature_name_1) %>%
	dplyr::arrange(desc(score)) %>%
	dplyr::slice(1:50) %>%
	dplyr::select(-score) %>%
	gene_gene_report()

report %>%
	readr::write_tsv(paste0("product/", feature_name, "_gene_gene_report_20201204.tsv"))

# use 51 because we're including the self-self which has correlation 1 at the top
neighborhood <- CalCEN %>%
	dplyr::filter(feature_name_1 == feature_name) %>%
	dplyr::arrange(desc(score)) %>%
	dplyr::slice(1:51)

neighborhood_report <- CalCEN %>%
		dplyr::semi_join(neighborhood, by = c("feature_name_1" = "feature_name_2")) %>%
		dplyr::semi_join(neighborhood, by = c("feature_name_2" = "feature_name_2")) %>%
		dplyr::filter(feature_name_1 < feature_name_2) %>%
		dplyr::filter(score > .95) %>%
		dplyr::select(-score) %>%
		gene_gene_report()

neighborhood_report %>%
	readr::write_tsv(paste0("product/", feature_name, "_neighborhood_report_20210109.tsv"))


# C1_01070C_A
report <- CalCEN %>%
	dplyr::filter(feature_name_1 == "C1_01070C_A") %>%
	dplyr::filter(feature_name_2 != feature_name_1) %>%
	dplyr::arrange(desc(score)) %>%
	dplyr::slice(1:50) %>%
	dplyr::select(-score) %>%
	gene_gene_report()

report %>%
	readr::write_tsv("product/C1_01070C_A_gene_gene_report_20201206.tsv")

# C2_04370W_A
report <- CalCEN %>%
	dplyr::filter(feature_name_1 == "C2_04370W_A") %>%
	dplyr::filter(feature_name_2 != feature_name_1) %>%
	dplyr::arrange(desc(score)) %>%
	dplyr::slice(1:50) %>%
	dplyr::select(-score) %>%
	gene_gene_report()

report %>%
	readr::write_tsv("product/C2_04370W_A_gene_gene_report_20201206.tsv")

# C1_09670C_A
report <- CalCEN %>%
	dplyr::filter(feature_name_1 == "C1_09670C_A") %>%
	dplyr::filter(feature_name_2 != feature_name_1) %>%
	dplyr::arrange(desc(score)) %>%
	dplyr::slice(1:50) %>%
	dplyr::select(-score) %>%
	gene_gene_report()

report %>%
	readr::write_tsv("product/C1_09670C_A_gene_gene_report_20201206.tsv")


feature_name <- "C1_05840W_A"
report <- CalCEN %>%
	dplyr::filter(feature_name_1 == feature_name) %>%
	dplyr::filter(feature_name_2 != feature_name_1) %>%
	dplyr::arrange(desc(score)) %>%
	dplyr::slice(1:50) %>%
	dplyr::select(-score) %>%
	gene_gene_report()

report %>%
	readr::write_tsv(paste0("product/", feature_name, "_gene_gene_report.tsv"))

feature_name <- "CR_04080C_A"
report <- CalCEN %>%
	dplyr::filter(feature_name_1 == feature_name) %>%
	dplyr::filter(feature_name_2 != feature_name_1) %>%
	dplyr::arrange(desc(score)) %>%
	dplyr::slice(1:50) %>%
	dplyr::select(-score) %>%
	gene_gene_report()

report %>%
	readr::write_tsv(paste0("product/", feature_name, "_gene_gene_report.tsv"))


feature_name <- "C4_06590W_A"
report <- CalCEN %>%
	dplyr::filter(feature_name_1 == feature_name) %>%
	dplyr::filter(feature_name_2 != feature_name_1) %>%
	dplyr::arrange(desc(score)) %>%
	dplyr::slice(1:50) %>%
	dplyr::select(-score) %>%
	gene_gene_report()

report %>%
	readr::write_tsv(paste0("product/", feature_name, "_gene_gene_report.tsv"))


feature_name <- "C1_08550C_A"
report <- CalCEN %>%
	dplyr::filter(feature_name_1 == feature_name) %>%
	dplyr::filter(feature_name_2 != feature_name_1) %>%
	dplyr::arrange(desc(score)) %>%
	dplyr::slice(1:50) %>%
	dplyr::select(-score) %>%
	gene_gene_report()

report %>%
	readr::write_tsv(paste0("product/", feature_name, "_gene_gene_report.tsv"))

##
# compare GCF1 and C6_03200W_A
feature_name_a <- "C6_03200W_A"
feature_name_b <- "C1_08550C_A" # GCF1
overlap_report <- dplyr::inner_join(
		CalCEN %>%
		    dplyr::filter(feature_name_1 == feature_name_a) %>%
		    dplyr::filter(feature_name_2 != feature_name_a) %>%
    		dplyr::filter(score > .98) %>%
		    dplyr::rename(score_to_a = score),
		CalCEN %>%
		    dplyr::filter(feature_name_1 == feature_name_b) %>%
		    dplyr::filter(feature_name_2 != feature_name_b) %>%
    		dplyr::filter(score > .98) %>%
		    dplyr::rename(score_to_b = score),
		by = c("feature_name_2")) %>%
		dplyr::mutate(feature_name_1 = feature_name_1.x) %>%
		gene_gene_report()

# one hit C4_02880C_A (SDO1 oroholog)

# compare C6_03200W_A and SDO1 at .99
feature_name_a <- "C6_03200W_A"
feature_name_b <- "C4_02880C_A"
overlap_report_C6_SDO1 <- dplyr::inner_join(
		CalCEN %>%
		    dplyr::filter(feature_name_1 == feature_name_a) %>%
		    dplyr::filter(feature_name_2 != feature_name_a) %>%
    		dplyr::filter(score > .99) %>%
		    dplyr::rename(score_to_a = score),
		CalCEN %>%
		    dplyr::filter(feature_name_1 == feature_name_b) %>%
		    dplyr::filter(feature_name_2 != feature_name_b) %>%
    		dplyr::filter(score > .99) %>%
		    dplyr::rename(score_to_b = score),
		by = c("feature_name_2")) %>%
		dplyr::mutate(feature_name_1 = feature_name_1.x) %>%
		gene_gene_report()


# compare GCF1 and SDO1 at .99
feature_name_a <- "C1_08550C_A"
feature_name_b <- "C4_02880C_A"
overlap_report_GCF1_SDO1 <- dplyr::inner_join(
		CalCEN %>%
		    dplyr::filter(feature_name_1 == feature_name_a) %>%
		    dplyr::filter(feature_name_2 != feature_name_a) %>%
    		dplyr::filter(score > .99) %>%
		    dplyr::rename(score_to_a = score),
		CalCEN %>%
		    dplyr::filter(feature_name_1 == feature_name_b) %>%
		    dplyr::filter(feature_name_2 != feature_name_b) %>%
    		dplyr::filter(score > .99) %>%
		    dplyr::rename(score_to_b = score),
		by = c("feature_name_2")) %>%
		dplyr::mutate(feature_name_1 = feature_name_1.x) %>%
		gene_gene_report()

feature_name <- "C6_03200W_A"


feature_name <- "C1_01060W_A"
report <- CalCEN %>%
	dplyr::filter(feature_name_1 == feature_name) %>%
	dplyr::filter(feature_name_2 != feature_name_1) %>%
	dplyr::arrange(desc(score)) %>%
	dplyr::slice(1:50) %>%
	dplyr::select(-score) %>%
	gene_gene_report()

report %>%
	readr::write_tsv(paste0("product/", feature_name, "_gene_gene_report.tsv"))


feature_name <- "C2_09940W_A"
report <- CalCEN %>%
	dplyr::filter(feature_name_1 == feature_name) %>%
	dplyr::filter(feature_name_2 != feature_name_1) %>%
	dplyr::arrange(desc(score)) %>%
	dplyr::slice(1:50) %>%
	dplyr::select(-score) %>%
	gene_gene_report()

report %>%
	readr::write_tsv(paste0("product/", feature_name, "_gene_gene_report.tsv"))


feature_name_a <- "C1_01060W_A"
feature_name_b <- "C2_09940W_A"
score_threshold <- .95
overlap_report_MBF1_GCN4 <- dplyr::inner_join(
		CalCEN %>%
		    dplyr::filter(feature_name_1 == feature_name_a) %>%
		    dplyr::filter(feature_name_2 != feature_name_a) %>%
    		dplyr::filter(score > score_threshold) %>%
		    dplyr::rename(score_to_a = score),
		CalCEN %>%
		    dplyr::filter(feature_name_1 == feature_name_b) %>%
		    dplyr::filter(feature_name_2 != feature_name_b) %>%
    		dplyr::filter(score > score_threshold) %>%
		    dplyr::rename(score_to_b = score),
		by = c("feature_name_2")) %>%
		dplyr::mutate(feature_name_1 = feature_name_1.x) %>%
		gene_gene_report()


feature_name <- "C1_09420W_A"
report <- CalCEN %>%
	dplyr::filter(feature_name_1 == feature_name) %>%
	dplyr::filter(feature_name_2 != feature_name_1) %>%
	dplyr::arrange(desc(score)) %>%
	dplyr::slice(1:50) %>%
	dplyr::select(-score) %>%
	gene_gene_report()

report %>%
	readr::write_tsv(paste0("product/", feature_name, "_gene_gene_report.tsv"))

feature_name <- "C7_00260C_A"
report <- CalCEN %>%
	dplyr::filter(feature_name_1 == feature_name) %>%
	dplyr::filter(feature_name_2 != feature_name_1) %>%
	dplyr::arrange(desc(score)) %>%
	dplyr::slice(1:50) %>%
	dplyr::select(-score) %>%
	gene_gene_report()

report %>%
	readr::write_tsv(paste0("product/", feature_name, "_gene_gene_report.tsv"))


feature_name <- "C5_03330C_A"
report <- CalCEN %>%
	dplyr::filter(feature_name_1 == feature_name) %>%
	dplyr::filter(feature_name_2 != feature_name_1) %>%
	dplyr::arrange(desc(score)) %>%
	dplyr::slice(1:50) %>%
	dplyr::select(-score) %>%
	gene_gene_report()

report %>%
	readr::write_tsv(paste0("product/", feature_name, "_gene_gene_report.tsv"))
