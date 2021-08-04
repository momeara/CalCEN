# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:

library(plyr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(alluvial)


load("intermediate_data/silac_clusters.Rdata")

# #v1
# wt_order <- tibble::tibble(
#    wildtype=rev(c(
#				"YJR142W",
#				"YFR016C",
#				"UBP15",
#				"TSA1",
#				"TRM82",
#				"SEC17",
#				"RPN5",
#				"CKA2",
#				"PAA1",
#				"YPR089W",
#				"RBP1",
#				"C7_00250C_A",
#				"C2_05550W_A",
#				"C1_10500W_A",
#				"AGM1",
#				"HTS1",
#				"MIS11",
#				"RPN11",
#				"RPL20B",
#				"RAD27",
#				"PRM15",
#				"NNR1",
#				"NIT3",
#				"MUD2",
#				"MTR10",
#				"MSC6",
#				"MCT1",
#				"HAM1",
#				"ENG1",
#				"ECM4",
#				"C5_01190W_A",
#				"C3_07330W_A",
#				"ACT1",
#				"ACO1",
#				"ABP140",
#				"CHC1",
#				"GCY1"))) %>%
#     dplyr::mutate(
#         wt_sort_order=1:n())
#
# gt_order <- tibble::tibble(
#     genetic=rev(c(
#				"SEC65",
#				"RPS21B",
#             "RPP2B",
#				"PTP3",
#				"POL30",
#				"NMT1",
#				"MRPL35",
#				"EBP7",
#				"DRS1",
#				"C2_00400C_A",
#				"BNA5",
#				"SPS20",
#				"LPD1",
#				"DAP2",
#				"BMT6",
#				"YPR089W",
#				"YFR016C",
#				"UBP15",
#				"TRM82",
#				"SEC17",
#				"RPN5",
#				"RPN11",
#				"RPL20B",
#				"RAD27",
#				"PRM15",
#				"NNR1",
#				"NIT3",
#				"MUD2",
#				"MTR10",
#				"MSC6",
#				"MCT1",
#				"HAM1",
#				"ECM4",
#				"C5_01190W_A",
#				"C3_07330W_A",
#				"AGM1",
#				"ACT1",
#				"ACO1",
#				"ABP140",
#				"CHC1",
#				"GCY1"
#     ))) %>%
#     dplyr::mutate(
#         gt_sort_order=1:n())
#
#

#v2
gt_order <- tibble::tibble(
		genetic=rev(c(
			"PTP3",
			"SEC17",
			"SEC65",
			"DAP2",
			"MUD2",
			"UBP15",
			"RPP2B",
			"POL30",
			"DRS1",
			"NMT1",
			"NNR1",
			"C2_00400C_A",
			"BNA5",
			"SPS20",
			"MRPL35",
			"AGM1",
			"EBP7",
			"BMT6",
			"TRM82",
			"RPN11",
			"C5_01190W_A",
			"RAD27",
			"LPD1",
			"RPL20B",
			"PRM15",
			"RPS21B",
			"NIT3",
			"MTR10",
			"RPN5",
			"HAM1",
			"ECM4",
			"C3_07330W_A",
			"ACO1",
			"ABP140",
			"YPR089W",
			"YFR016C",
			"ACT1",
			"MCT1",
			"CHC1",
			"MSC6",
			"GCY1"
		))) %>%
		dplyr::mutate(
				gt_sort_order=1:n())
wt_order <- tibble::tibble(
	 wildtype=rev(c(
			"SEC17",
			"YJR142W",
			"MUD2",
			"UBP15",
			"CKA2",
			"PAA1",
			"C2_05550W_A",
			"C1_10500W_A",
			"NNR1",
			"RBP1",
			"C7_00250C_A",
			"AGM1",
			"HTS1",
				"TRM82",
			"RPN11",
			"C5_01190W_A",
			"RAD27",
			"TSA1",
			"MIS11",
			"RPL20B",
			"PRM15",
			"ENG1",
			"NIT3",
			"MTR10",
			"RPN5",
			"HAM1",
			"ECM4",
			"C3_07330W_A",
			"ACO1",
			"ABP140",
			"YPR089W",
			"YFR016C",
			"ACT1",
			"MCT1",
			"CHC1",
			"MSC6",
			"GCY1"
	 ))) %>%
		dplyr::mutate(
				wt_sort_order=1:n())
ph_order <- tibble::data_frame(
	 pharmacological=rev(c(
			"SEC17",
			"MAS2",
			"DAP2",
			"MUD2",
			"UBP15",
			"CDC73",
			"PAA1",
			"C2_05550W_A",
			"NNR1",
				"RBP1",
			"SPS20",
			"AGM1",
			"OSH2",
			"BMT6",
			"TRM82",
			"RPN11",
			"RAD27",
			"LPD1",
			"RPL20B",
			"PRM15",
			"NIT3",
			"MTR10",
			"RPN5",
			"HAM1",
			"ECM4",
			"C3_07330W_A",
			"ACO1",
			"ABP140",
			"YPR089W",
			"YFR016C",
			"CHC1",
			"GCY1"
	 ))) %>%
		dplyr::mutate(
				ph_sort_order=1:n())
z <- silac_clusters %>%
		dplyr::distinct(gene_label, condition, exemplar, .keep_all=TRUE) %>%
		dplyr::select(gene_label, condition, exemplar) %>%
		tidyr::spread(condition, exemplar) %>%
		dplyr::count(genetic, wildtype, pharmacological) %>%
		dplyr::mutate(
				cannonical_order = 1:n())
z <- z %>%
		dplyr::left_join(wt_order, by="wildtype") %>%
		dplyr::arrange(wt_sort_order) %>%
		dplyr::mutate(wt_order = 1:n()) %>%
		dplyr::left_join(gt_order, by="genetic") %>%
		dplyr::arrange(gt_sort_order) %>%
		dplyr::mutate(gt_order = 1:n()) %>%
		dplyr::left_join(ph_order, by="pharmacological") %>%
		dplyr::arrange(ph_sort_order) %>%
		dplyr::mutate(ph_order = 1:n()) %>%
		dplyr::arrange(cannonical_order)
source("scripts/alluvial.R")
pdf(
		file="product/figures/silac_clusters_alluvial_gt_wt_ph_180917.pdf",
		height=25,
		useDingbats=FALSE)
p <- alluvial(
		z[,c("genetic", "wildtype", "pharmacological")],
		col=ifelse(z$wildtype == "RPL20B" | z$genetic == "RPL20B", "orange", "grey"),
		freq=z$n,
		cw=.2,
		ordering=list(gt_order=z$gt_order, wt_order=z$wt_order, ph_order=z$ph_order))
dev.off()



pdf(
		file="product/figures/silac_clusters_alluvial.pdf",
		height=25,
		useDingbats=FALSE)

alluvial::alluvial(
		z[1:3],
		col=ifelse(z$wildtype == "RPL20B", "orange", "grey"),
		freq=z$n,
		cw=.2)

dev.off()


pdf(
		file="product/figures/silac_clusters_alluvial_wt_ph.pdf",
		height=25,
		useDingbats=FALSE)
alluvial::alluvial(
		z[,c(2,3)],
		col=ifelse(z$wildtype == "RPL20B", "orange", "grey"),
		freq=z$n,
		cw=.2)
dev.off()



pdf(
		file="product/figures/test_alluvial.pdf")
tit <- as.data.frame(Titanic, stringsAsFactors = FALSE)
tit %>% group_by(Class, Survived) %>%
		summarise(n = sum(Freq)) -> tit2d
alluvial::alluvial(tit2d[,1:2], freq=tit2d$n)
dev.off()



class_sort <- c("Crew", "2nd", "3rd", "1st")
survived_sort <- c("No", "Yes")

tit <- as.data.frame(Titanic, stringsAsFactors = FALSE)
tit2d <- tit %>%
		dplyr::group_by(Class, Survived) %>%
		dplyr::summarise(n = sum(Freq)) %>%
		dplyr::ungroup()

tit2d <- tit2d %>%
	dplyr::mutate(canonical_order=1:n()) %>%
	# add custom Class sort
	dplyr::right_join(
			tibble::data_frame(
					Class=class_sort) %>%
					dplyr::mutate(class_sort_order=1:n()),
			by="Class") %>%
	dplyr::arrange(class_sort_order, canonical_order) %>%
	dplyr::mutate(class_order=1:n()) %>%
	# add custom Survived sort
	dplyr::right_join(
			tibble::data_frame(
					Survived=survived_sort) %>%
					dplyr::mutate(survived_sort_order=1:n()),
			by="Survived") %>%
	dplyr::arrange(survived_sort_order, canonical_order) %>%
	dplyr::mutate(survived_order=1:n()) %>%
	dplyr::arrange(canonical_order)


source("scripts/alluvial.R")
pdf(
		file="product/figures/test_alluvial.pdf")
alluvial(
		tit2d[,1:2],
		freq=tit2d$n,
		ordering=list(rev(tit2d$class_order), rev(tit2d$survived_order)))
dev.off()
