# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:

library(plyr)
library(dplyr)
library(readr)
library(httr)
library(glue)


load("intermediate_data/cgd_id_2_uniprot_accn.Rdata")
load("intermediate_data/chromosome_features.Rdata")


uniprot_to_cgd_id <- function(uniprot_accn, verbose=TRUE) {
	if(verbose){
		"Looking up {count} uniprot accessions for associated CGD IDs with uniprot...\n\n" %>%
			glue::glue(count=length(uniprot_accn)) %>%
			cat()
	}
	httr::GET(
		url="http://www.uniprot.org/uniprot/",
		httr::user_agent("httr mattjomeara@gmail.com"),
		query=list(
			query=uniprot_accn %>% paste0(collapse=" or "),
			format='tab',
			compress='yes',
			columns="id,database(CGD)")) %>%
		httr::content() %>%
		rawConnection() %>%
		gzcon() %>%
		readr::read_tsv(
			col_types=readr::cols(
				Entry = readr::col_character(),
				`Cross-reference (CGD)` = readr::col_character())) %>%
		dplyr::rename(
			uniprot_accn = Entry,
			cgd_id=`Cross-reference (CGD)`) %>%
		dplyr::filter(!is.na(cgd_id)) %>%
		dplyr::mutate(cgd_id = cgd_id %>% stringr::str_extract("^CAL[0-9]+"))
}

map_ca_uniprot_accn_to_feature_name <- function(data, verbose=TRUE){
	if(!("uniprot_accn" %in% names(data))){
		"ERROR: uniprot_accn must be a column in the input data.frame\n" %>% cat
		stop()
	}

	if(verbose){
		"# n distinct uniprot accession identifers : {count} \n\n" %>%
			glue::glue(
				count=data %>% dplyr::distinct(uniprot_accn) %>% nrow) %>%
			cat()
	}

	cgd_id_map <- cgd_id_2_uniprot_accn %>%
		dplyr::semi_join(data, by="uniprot_accn")
	if(verbose){
		"# n proteins in CGD: {count}\n\n" %>%
			glue::glue(count=cgd_id_map %>% distinct(uniprot_accn) %>% nrow) %>%
			cat()
	}

	uniprot_id_map <- data %>%
		dplyr::anti_join(cgd_id_2_uniprot_accn, by="uniprot_accn") %>%
		magrittr::extract2("uniprot_accn") %>%
		uniprot_to_cgd_id()
	if(verbose){
	 	"# n proteins not in CGD but in uniprot: {count}\n\n" %>%
			glue::glue(count=	uniprot_id_map %>% distinct(uniprot_accn) %>% nrow) %>%
			cat()
	}

	id_map <- rbind(cgd_id_map, uniprot_id_map) %>%
		dplyr::inner_join(
			chromosome_features %>% dplyr::select(feature_name, cgd_id=primary_cgd_id),
			by=c("cgd_id")) %>%
		dplyr::select(uniprot_accn, feature_name)
	if(verbose){
	 	"# n proteins mapped in total: {count}\n\n" %>%
			glue::glue(count=	id_map %>% distinct(uniprot_accn) %>% nrow) %>%
			cat()
	}

	id_map <- id_map %>%
		dplyr::semi_join(
			id_map %>%
				dplyr::count(uniprot_accn) %>%
				dplyr::filter(n==1),
			by="uniprot_accn")
	if(verbose){
	 	"# n proteins that map to exactly 1 gene: {count}\n\n" %>%
			glue::glue(count=id_map %>%  nrow) %>%
			cat()
	}

	data %>% dplyr::left_join(id_map, by="uniprot_accn")
}


