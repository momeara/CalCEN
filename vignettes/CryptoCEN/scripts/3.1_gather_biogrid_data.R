# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:

library(plyr)
library(dplyr)
library(stringr)
library(readr)
library(CalCEN)

load("intermediate_data/h9_transcript_annotations.Rdata")

system("\\
  mkdir -p raw_data/biogrid && \\
  mkdir cd raw_data/biogrid && \\
  wget https://downloads.thebiogrid.org/Download/BioGRID/Release-Archive/BIOGRID-3.4.161/BIOGRID-ORGANISM-3.4.161.tab2.zip && \\
  unzip BIOGRID-3.4.161/BIOGRID-ORGANISM-3.4.161.tab2.zip && \\
  ls | grep -v -e 'Candida_albicans' -e 'Saccharomyces_cerevisiae' | xargs rm")

ca_biogrid <- read_biogrid_tab2(
	fname="raw_data/biogrid/BIOGRID-ORGANISM-Candida_albicans_SC5314-3.4.161.tab2.txt",
	taxon=237561)

sac_biogrid <- read_biogrid_tab2(
	fname="raw_data/biogrid/BIOGRID-ORGANISM-Saccharomyces_cerevisiae_S288c-3.4.161.tab2.txt",
	taxon=559292) %>%
	dplyr::left_join(
		chromosome_features %>%
			dplyr::select(
				ca_feature_name_1 = feature_name, gene_symbol_1 = sac_ortholog),
		by="gene_symbol_1") %>%
	dplyr::left_join(
		chromosome_features %>%
			dplyr::select(
				ca_feature_name_2 = feature_name, gene_symbol_2 = sac_ortholog),
		by="gene_symbol_2")

# Sac HSP82 and HSC82 are both functionally compensating HSP90 paralog
sac_biogrid <- sac_biogrid %>%
	dplyr::mutate(
		ca_feature_name_1 = ifelse(gene_symbol_1 == "HSP82", "C7_02030W_A", ca_feature_name_1),
		ca_feature_name_2 = ifelse(gene_symbol_2 == "HSP82", "C7_02030W_A", ca_feature_name_2))



save(ca_biogrid, file="intermediate_data/ca_biogrid.Rdata")
save(sac_biogrid, file="intermediate_data/sac_biogrid.Rdata")
