# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:

library(plyr)
library(dplyr)
library(magrittr)
library(stringr)
library(readr)
library(seqinr)

load("intermediate_data/chromosome_features.Rdata")


ca_protein_sequences <- seqinr::read.fasta(
	file = "raw_data/C_albicans_SC5314_A22_current_default_protein.fasta",
	seqtype = "AA")

ca_genes <- ca_protein_sequences %>% seqinr::getName()

save(ca_genes, file="intermediate_data/ca_genes.Rdata")
