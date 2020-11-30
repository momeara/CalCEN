# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:

library(plyr)
library(dplyr)
library(purrr)
library(GO.db)

load("intermediate_data/cgd_go.Rdata")


# get all the go terms
go_terms <- GO.db::GO_dbconn() %>%
	dplyr::tbl("go_term") %>%
	dplyr::collect(n=Inf)

# get a map of each go term to its parents
go_parents <- GO.db::GO_dbconn() %>%
	dplyr::tbl("go_cc_parents") %>%
	dplyr::collect(n=Inf) %>%
	dplyr::left_join(go_terms %>% dplyr::select(`_id`, go_id), by="_id") %>%
	dplyr::left_join(go_terms %>% dplyr::select(`_parent_id`=`_id`, parent_go_id=go_id), by="_parent_id") %>%
	dplyr::select(-`_id`, -`_parent_id`) %>%
	dplyr::filter(parent_go_id != 'all')

ca_go <- cgd_go %>%
	dplyr::filter(taxon=="taxon:5476") %>%
	dplyr::inner_join(
		chromosome_features %>%
			dplyr::filter(feature_class == "ORF") %>%
			dplyr::select(feature_name, primary_cgd_id),
		by=c("db_object_id"="primary_cgd_id"))

ca_go_propagated <- ca_go %>%
	dplyr::inner_join(go_parents, by="go_id") %>%
	dplyr::anti_join(cgd_go, by=c("db_object_id"="db_object_id", "parent_go_id"="go_id")) %>%
	dplyr::mutate(go_id = parent_go_id) %>%
	dplyr::select(-parent_go_id, -relationship_type) %>%
	rbind(ca_go) %>%
	dplyr::left_join(go_terms %>% dplyr::select(-`_id`), by="go_id") %>%
	dplyr::filter(!is.na(ontology))

save(ca_go_propagated, file="intermediate_data/ca_go_propagated.Rdata")

ca_go_propagated_filtered <- ca_go_propagated %>%
	dplyr::semi_join(
		ca_go_propagated %>%
		dplyr::count(go_id) %>%
		dplyr::filter(20 <= n, n <= 1000),
	by=c("go_id"))

save(
	ca_go_propagated_filtered,
	file="intermediate_data/ca_go_propagated_filtered.Rdata")

# warning by filtering out NOT qualifiers some terms have less than 20 annotations
# default EGAD::run_GBA re-filters for terms with min=20 max=1000, so some terms
# will be excluded when the analysis is run.
# FIXME: move this filter before ca_go_propagated_filtered is computed
ca_go_positive <-  ca_go_propagated_filtered %>%
	dplyr::filter(is.na(qualifier) || qualifier != "NOT") %>%
	dplyr::distinct(feature_name, go_id, .keep_all=TRUE)


spread_annotations <- function(annotations){
	annotations %>%
		dplyr::transmute(feature_name, go_id, value=1) %>%
		tidyr::spread(go_id, value, fill=0) %>%
		as.data.frame() %>%
		tibble::column_to_rownames("feature_name")
}


ca_go_annotations <- ca_go_positive %>% spread_annotations
save(ca_go_annotations, file="intermediate_data/ca_go_annotations.Rdata")


ca_go_annotations_by_subontology <- list(
	all = ca_go_positive,
	MF = ca_go_positive %>% dplyr::filter(ontology == "MF"),
	BP = ca_go_positive %>% dplyr::filter(ontology == "BP"),
	CC = ca_go_positive %>%	dplyr::filter(ontology == "CC")) %>%
	plyr::llply(spread_annotations)
save(ca_go_annotations_by_subontology, file="intermediate_data/ca_go_annotations_by_subontology.Rdata")


ca_go_annotations_by_evidence <- list(
	  all = ca_go_positive,

		experimental = ca_go_positive %>%
			dplyr::filter(evidence %in% c("EXP", "IDA", "IPI", "IMP", "IGI", "IEP")),
		EXP = ca_go_positive %>% dplyr::filter(evidence == "EXP"),
		IDA = ca_go_positive %>% dplyr::filter(evidence == "IDA"),
		IPI = ca_go_positive %>% dplyr::filter(evidence == "IPI"),
		IMP = ca_go_positive %>% dplyr::filter(evidence == "IMP"),
		IGI = ca_go_positive %>% dplyr::filter(evidence == "IGI"),
		IEP = ca_go_positive %>% dplyr::filter(evidence == "IEP"),

		high_throughput = ca_go_positive %>%
			dplyr::filter(evidence %in% c("HTP", "HDA", "HMP", "HGI", "HEP")),
		HTP = ca_go_positive %>% dplyr::filter(evidence == "HTP"),
		HDA = ca_go_positive %>% dplyr::filter(evidence == "HDA"),
		HMP = ca_go_positive %>% dplyr::filter(evidence == "HMP"),
		HGI = ca_go_positive %>% dplyr::filter(evidence == "HGI"),
		HEP = ca_go_positive %>% dplyr::filter(evidence == "HEP"),

		computational_analysis = ca_go_positive %>%
			dplyr::filter(evidence %in% c("ISS", "ISO", "ISA", "ISM", "IGC", "IBA", "IBD", "IKR", "IRD", "RCA")),
		ISS = ca_go_positive %>% dplyr::filter(evidence == "ISS"),
		ISO = ca_go_positive %>% dplyr::filter(evidence == "ISO"),
		ISA = ca_go_positive %>% dplyr::filter(evidence == "ISA"),
		ISM = ca_go_positive %>% dplyr::filter(evidence == "ISM"),
		IGC = ca_go_positive %>% dplyr::filter(evidence == "IGC"),
		IBA = ca_go_positive %>% dplyr::filter(evidence == "IBA"),
		IBD = ca_go_positive %>% dplyr::filter(evidence == "IBD"),
		IKR = ca_go_positive %>% dplyr::filter(evidence == "IKR"),
		IRD = ca_go_positive %>% dplyr::filter(evidence == "IRD"),
		RCA = ca_go_positive %>% dplyr::filter(evidence == "RCA"),

		author_statement = ca_go_positive %>% dplyr::filter(evidence %in% c("TAS", "NAS")),
		TAS = ca_go_positive %>% dplyr::filter(evidence == "TAS"),

		curator_statement = ca_go_positive %>% dplyr::filter(evidence %in% c("IC", "ND")),
		ID = ca_go_positive %>% dplyr::filter(evidence == "ID"),

		electronic_annotation = ca_go_positive %>% dplyr::filter(evidence %in% c("IEA")),
		IEA = ca_go_positive %>% dplyr::filter(evidence == "IEA")) %>%
	purrr::keep(~nrow(.x) > 10) %>%
	plyr::llply(spread_annotations)
save(ca_go_annotations_by_evidence, file="intermediate_data/ca_go_annotations_by_evidence.Rdata")



