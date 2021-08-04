# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:

library(plyr)
library(dplyr)
library(stringr)
library(readr)


### CHROMOSOMAL FEATURES ###
system("cd raw_data && wget https://downloads.yeastgenome.org/curation/chromosomal_feature/SGD_features.tab")

sac_chromosome_features <- readr::read_tsv(
	file="raw_data/SGD_features.tab",
	col_names=c(
		"primary_sgd_id",				# Primary SGDID (mandatory)
		"feature_type",					# Feature type (mandatory)
		"feature_qualifier",		# Feature qualifier (optional)
		"feature_name",					# Feature name (optional)
		"standard_gene_name",		# Standard gene name (optional)
		"alias",								# Alias (optional, multiples separated by |)
		"parent_feature_name",	# Parent feature name (optional)
		"secondary_sgd_id",			# Secondary SGDID (optional, multiples separated by |)
		"chromosome",						# Chromosome (optional)
		"start_coordinate",			# Start_coordinate (optional)
		"stop_coordinate",			# Stop_coordinate (optional)
		"strand",								# Strand (optional)
		"genetic_position",			# Genetic position (optional)
		"coodinate_version",		# Coordinate version (optional)
		"sequence_version",			# Sequence version (optional)
		"description"),					# Description (optional)
	col_types=readr::cols(
		primary_sgd_id = readr::col_character(),
		feature_type = readr::col_character(),
		feature_qualifier = readr::col_character(),
		feature_name = readr::col_character(),
		standard_gene_name = readr::col_character(),
		alias = readr::col_character(),
		parent_feature_name = readr::col_character(),
		secondary_sgd_id = readr::col_character(),
		chromosome = readr::col_double(),
		start_coordinate = readr::col_double(),
		stop_coordinate = readr::col_double(),
		strand = readr::col_character(),
		genetic_position = readr::col_double(),
		coodinate_version = readr::col_date(format = ""),
		sequence_version = readr::col_date(format = ""),
		description = readr::col_character()))

save(sac_chromosome_features, file="intermediate_data/sac_chromosome_features.Rdata")
