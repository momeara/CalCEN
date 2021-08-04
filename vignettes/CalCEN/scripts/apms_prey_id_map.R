# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:

library(plyr)
library(dplyr)
library(stringr)



prey_id_map <- function(apms_set) {
	map <- apms_set %>%
		dplyr::distinct(prey_gene) %>%
		dplyr::left_join(
			chromosome_features %>%
				dplyr::filter(feature_name %>% stringr::str_detect("A$")) %>%
				dplyr::select(
					prey_feature_name = feature_name,
					prey_gene = gene_name,
					prey_feature_status = feature_status,
					prey_feature_type = feature_type),
				by=c("prey_gene")) %>%
		dplyr::mutate(
			prey_feature_name = dplyr::case_when(
				prey_gene %>% stringr::str_detect("C[0-9R]_[0-9]+[CW]_A") ~ prey_gene,
				prey_gene == "orf19.12401" ~  "C5_04000W_A",
				prey_gene == "orf19.3377" ~ "C1_07880C_A",
				prey_gene == "ALD99" ~ "C2_02970C_A",
				prey_gene == "BUL2" ~ "C1_00540C_A",
				prey_gene == "CAP001" ~ "CR_04190W_A",
				prey_gene == "DOR14" ~ "C1_13860C_A",
				prey_gene == "LPF2" ~ NA_character_,
				prey_gene == "LPF36" ~ "C7_02850W_A",
				prey_gene == "LYS211" ~ "C2_04460W_A",
				prey_gene == "PHO131" ~ "C1_07230W_A",
				prey_gene == "PKH1" ~ "C1_12410C_A",
				prey_gene == "RPB140" ~ "C1_01590C_A",
				prey_gene == "RPL4" ~ "C1_14110C_A",
				prey_gene == "TDH1" ~ "C3_06870W_A",
				prey_gene == "TEF3" ~ "C5_01580C_A",
				prey_gene == "TIF1" ~ "C1_01350C_A",
				prey_gene == "TKT1" ~ "C1_08320W_A",
				prey_gene == "VIG9" ~ "C3_07950C_A",
				prey_gene == "FUN11" ~ "C1_08100W_A",
				prey_gene == "SDS23" ~ "C1_08370W_A",
				prey_gene == "PRSS1" ~ "NA_character_",
				prey_gene == "orf19.4517" ~ "C2_04370W_A",
				prey_gene == "orf19.6396" ~ "CR_08300C_A",
				prey_gene == "orf19.12956" ~ "C6_02430W_A",
				prey_gene == "orf19.9458" ~ "C2_07340W_A", # NOC4
				prey_gene == "orf19.10395" ~ "C4_06570C_A", # PDC11
				prey_gene == "RPL1" ~ "C6_02240C_A", # RPL10A
				prey_gene == "orf19.8442" ~ "C2_04010C_A", # HSP21
				prey_gene == "orf19.5328" ~ "C2_10550C_A", # GCN1
				prey_gene == "orf19.10642" ~ "C4_06790W_A",
				prey_gene == "orf19.5949" ~ "C3_04830C_A", # FAS2
				prey_gene == "orf19.6701" ~ "C7_03660C_A",
				prey_gene == "orf19.7522" ~ "CR_00130C_A",
				prey_gene == "RPS26" ~ "C2_01610C_A", # RPS26A
				prey_gene == "orf19.11759" ~ "C5_02660C_A",
				prey_gene == "orf19.7459" ~ "C3_06700C_A",
				prey_gene == "orf19.6463" ~ "C7_02460C_A",
				prey_gene == "orf19.9952" ~ "CR_03120W_A",
				prey_gene == "orf19.8595" ~ "C5_00180W_A", # AFL1
				prey_gene == "orf19.12091" ~ "C4_01730C_A",
				prey_gene == "orf19.1956" ~ "C5_01140C_A",
				prey_gene == "orf19.5801" ~ "C2_03010C_A", # RNR21
				prey_gene == "orf19.12932" ~ "C2_06170C_A", # ECM17
				prey_gene == "orf19.6457" ~ "C7_02530C_A",
				prey_gene == "orf19.5356" ~ "C2_10740C_A",
				prey_gene == "orf19.11525" ~ "C5_05440C_A",
				prey_gene == "orf19.8065" ~ "NA_character_",
				prey_gene == "orf19.9943" ~ "CR_03200C_A", # POP1
				prey_gene == "orf19.5429" ~ "C3_00380C_A",
				prey_gene == "orf19.7913" ~ "C3_02900W_A",
				prey_gene == "orf19.3192" ~ "C5_01820W_A", # STI1
				prey_gene == "LPF20" ~ "C1_06200W_A",
				prey_gene == "orf19.12723" ~ "C1_12070C_A",
				TRUE ~ prey_feature_name))

	map <- rbind(
		map %>%
			dplyr::filter(!is.na(prey_feature_name)),
		map %>%
			dplyr::filter(is.na(prey_feature_name)) %>%
			dplyr::select(prey_gene) %>%
			dplyr::left_join(
				chromosome_features %>%
					dplyr::select(
						prey_gene = sac_ortholog,
						prey_feature_name = feature_name,
						prey_feature_status = feature_status,
						prey_feature_type = feature_type),
				by="prey_gene")) %>%
		dplyr::filter(!is.na(prey_feature_name))

	apms_set <- apms_set %>%
		dplyr::left_join(map, by="prey_gene") %>%
		dplyr::mutate(
			prey_gene = dplyr::case_when(
				prey_gene == "orf19.9458" ~ "NOC4",
				prey_gene == "orf19.10395" ~ "PDC11",
				prey_gene == "RPL1" ~ "RPL10A",
				prey_gene == "orf19.8442" ~ "HSP21",
				prey_gene == "orf19.5328" ~ "GCN1",
				prey_gene == "orf19.5949" ~ "FAS2",
				prey_gene == "RPS26" ~ "RPS26A",
				prey_gene == "orf19.8595" ~ "AFL1",
				prey_gene == "orf19.5801" ~ "RNR21",
				prey_gene == "orf19.12932" ~ "ECM17",
				prey_gene == "orf19.9943" ~ "POP1",
				prey_gene == "orf19.3192" ~ "STI1",
				prey_gene == "orf19.5429" ~ "DIA2",
				prey_gene == "orf19.10642" ~ "GPN3",
				prey_gene == "orf19.7459" ~ "YBR238C",
				prey_gene == "orf19.6457" ~ "YBL086C",
				prey_gene == "orf19.6463" ~ "NPA3",
				prey_gene == "orf19.12956" ~ "ITC1",
				prey_gene == "orf19.6701" ~ "YHR020W",
				prey_gene == "orf19.9952" ~ "MPM1",
				prey_gene == "orf19.6396" ~ "NTE1",
				prey_gene == "orf19.12091" ~ "PBY1",
				prey_gene == "orf19.1956" ~ "MSC6",
				prey_gene == "orf19.5356" ~ "EMW1",
				prey_gene == "orf19.12401" ~ "RFC3",
				prey_gene == "orf19.3377" ~ "GSH1",
				TRUE ~ prey_gene))
}
