# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:

library(plyr)
library(dplyr)
library(stringr)
library(readr)


### CHROMOSOMAL FEATURES ###
system("cd raw_data && wget http://www.candidagenome.org/download/chromosomal_feature_files/C_albicans_SC5314/C_albicans_SC5314_A22_current_chromosomal_feature.tab")

chromosome_features <- readr::read_tsv(
	file="raw_data/C_albicans_SC5314_A22_current_chromosomal_feature.tab",
	skip=8,
	col_names=c(
		"feature_name",
		"gene_name",
		"aliases",
		"feature_type",
		"chromosome",
		"start_coordinate",
		"stop_coordinate",
		"strand",
		"primary_cgd_id",
		"secondary_cgd_id",
		"description",
		"creation_date",
		"sequence_coordinate_date",
		"null_1",
		"null_2",
		"gene_name_date",
		"gene_name_is_standard",
		"sac_ortholog"),
	col_types=readr::cols(
	  feature_name = readr::col_character(),
	  gene_name = readr::col_character(),
	  aliases = readr::col_character(),
	  feature_type = readr::col_character(),
	  chromosome = readr::col_character(),
	  start_coordinate = readr::col_double(),
	  stop_coordinate = readr::col_double(),
	  strand = readr::col_character(),
	  primary_cgd_id = readr::col_character(),
	  secondary_cgd_id = readr::col_character(),
	  description = readr::col_character(),
	  creation_date = readr::col_date(format = ""),
	  sequence_coordinate_date = readr::col_date(format = ""),
	  null_1 = readr::col_logical(),
	  null_2 = readr::col_logical(),
	  gene_name_date = readr::col_date(format = ""),
	  gene_name_is_standard = readr::col_character(),
	  sac_ortholog = readr::col_character())) %>%
	dplyr::select(-null_1, -null_2) %>%
	dplyr::mutate(
		feature_class = dplyr::case_when(
			feature_type %>% stringr::str_detect("^ORF") ~ "ORF",
			feature_type %>% stringr::str_detect("RNA") ~ "RNA",
			feature_type %>% stringr::str_detect("repeat") ~ "repeat",
			feature_type %>% stringr::str_detect("gene") ~ "gene_like",
			feature_type %>% stringr::str_detect("blocked_reading_frame") ~ "gene_like",
			feature_type %>% stringr::str_detect("centromere") ~ "centromere",
			TRUE ~ NA_character_),
		feature_status = dplyr::case_when(
			feature_type %>% stringr::str_detect("Uncharacterized") ~ "Uncharacterized",
			feature_type %>% stringr::str_detect("Verified") ~ "Verified",
			feature_type %>% stringr::str_detect("Dubious") ~ "Dubious",
			TRUE ~ NA_character_),
		is_alternatively_spliced = feature_type %>%
			stringr::str_detect("Alternatively spliced"),
		is_transposable_element = feature_type %>%
			stringr::str_detect("transposable element gene"),
		rna_type = dplyr::case_when(
			feature_class == "RNA" ~ feature_type %>% str_detect("^[a-z]+"),
			TRUE ~ NA),
		is_blocked_reading_frame = feature_type %>%
			stringr::str_detect("blocked_reading_frame"),
		gene_name_is_standard = gene_name_is_standard == "Y")
save(chromosome_features, file="intermediate_data/chromosome_features.Rdata")


### MAP ENTREZ GENE ID ###
system("cd raw_data && wget http://www.candidagenome.org/download/External_id_mappings/CGDID_2_GeneID.tab.gz")
cgd_id_2_gene_id <- readr::read_tsv(
	file="raw_data/CGDID_2_GeneID.tab.gz",
	col_types=readr::cols(
		CGDID = readr::col_character(),
		`Entrez GeneID` = readr::col_double())) %>%
	dplyr::rename(
		cgd_id = CGDID,
		gene_id = `Entrez GeneID`)
save(cgd_id_2_gene_id, file="intermediate_data/cgd_id_2_gene_id.Rdata")


### MAP UNIPROT ENTRY ###
system("cd raw_data && wget http://www.candidagenome.org/download/External_id_mappings/gp2protein.cgd.gz")

cgd_id_2_uniprot_accn <- readr::read_tsv(
	file="raw_data/gp2protein.cgd.gz",
	col_names=c("cgd_id", "uniprot_accn"),
	col_types=readr::cols(
		cgd_id=readr::col_character(),
		uniprot_accn=readr::col_character())) %>%
	dplyr::transmute(
		primary_cgd_id = cgd_id %>% stringr::str_extract("[^:]+$"),
		uniprot_accn = uniprot_accn %>% stringr::str_extract("[^:]+$")) %>%
	dplyr::left_join(
		chromosome_features %>%
			dplyr::select(
				feature_name,
				gene_name,
				primary_cgd_id),
			by="primary_cgd_id") %>%
	dplyr::select(
		feature_name,
		gene_name,
		primary_cgd_id,
		uniprot_accn)

save(cgd_id_2_uniprot_accn, file="intermediate_data/cgd_id_2_uniprot_accn.Rdata")



chromosome_feature_to_alias <- chromosome_features %>%
	dplyr::select(feature_name, aliases) %>%
	plyr::adply(1, function(feature){
		tibble::tibble(
			alias=feature$aliases[1] %>%
				stringr::str_split("[|]") %>%
				unlist())
	}) %>%
	dplyr::select(-aliases)
save(chromosome_feature_to_alias, file="intermediate_data/chromosome_feature_to_alias.Rdata")


### MAP GO ANNOTATIONS ###
system("cd raw_data && wget http://www.candidagenome.org/download/go/gene_association.cgd.gz")
cgd_go <- readr::read_tsv(
	file="raw_data/gene_association.cgd.gz",
	skip=19,
	col_names=c(
		"db",
		"db_object_id",
		"db_object_symbol",
		"qualifier",
		"go_id",
		"db_reference",
		"evidence",
		"with_or_from",
		"aspect",
		"db_object_name",
		"db_object_synonym",
		"db_object_type",
		"taxon",
		"date",
		"assigned_by"),
	col_types=readr::cols(
		db = readr::col_character(),
		db_object_id = readr::col_character(),
		db_object_symbol = readr::col_character(),
		qualifier = readr::col_character(),
		go_id = readr::col_character(),
		db_reference = readr::col_character(),
		evidence = readr::col_character(),
		with_or_from = readr::col_character(),
		aspect = readr::col_character(),
		db_object_name = readr::col_logical(),
		db_object_synonym = readr::col_character(),
		db_object_type = readr::col_character(),
		taxon = readr::col_character(),
		date = readr::col_double(),
		assigned_by = readr::col_character()))
save(cgd_go, file="intermediate_data/cgd_go.Rdata")


#### Sac Orthologs #####
system("cd raw_data && wget http://www.candidagenome.org/download/homology/orthologs/C_albicans_SC5314_S_cerevisiae_by_CGOB/C_albicans_SC5314_S_cerevisiae_orthologs.txt")
ca_sac_orthologs <- readr::read_tsv(
	file="raw_data/C_albicans_SC5314_S_cerevisiae_orthologs.txt",
	skip=8,
	col_names=c(
		"feature_name",
		"gene_name",
		"cgd_id",
		"sac_feature_name",
		"sac_gene_name",
		"ygd_id"),
	col_types=c(
		feature_name = readr::col_character(),
		gene_name = readr::col_character(),
		cgd_id = readr::col_character(),
		sac_feature_name = readr::col_character(),
		sac_gene_name = readr::col_character(),
		ygd_id = readr::col_character()))
save(ca_sac_orthologs, file="intermediate_data/ca_sac_orthologs.Rdata")


### Candida albicans SC5314 sequence assemblies in the Generic Feature Format (GFF) ###
system("cd raw_data && wget http://www.candidagenome.org/download/gff/C_albicans_SC5314/Assembly22/C_albicans_SC5314_version_A22-s07-m01-r123_features.gff")
system("cd raw_data ln -s C_albicans_SC5314_version_A22-s07-m01-r123_features.gff C_albicans_SC5314_A22_current_features.gff")
system("cd raw_data && wget http://www.candidagenome.org/download/gff/C_albicans_SC5314/Assembly22/C_albicans_SC5314_version_A22-s07-m01-r123_features.gtf")
system("cd raw_data ln -s C_albicans_SC5314_version_A22-s07-m01-r123_features.gtf C_albicans_SC5314_A22_current_features.gtf")

system("cd raw_data && wget http://www.candidagenome.org/download/gff/C_albicans_SC5314/Assembly22/C_albicans_SC5314_version_A22-s07-m01-r123_intergenic.gff")
system("cd raw_data ln -s C_albicans_SC5314_version_A22-s07-m01-r123_intergenic.gff C_albicans_SC5314_A22_current_intergenic.gff")


# this is used for cross map analysis among other things
