#' Embed a network in a symmmetric network over given set of nodes
#'
#' This function filters the input network to those in genes.list
#' and assigns the values to the cells corresponding to row/columns in genes.list
#'
#' @param network represented as a dense matrix
#' @param genes.list representing the nodes in the target network
#' @param missing values to assign to edges that are not present in the input network
#' @return network matrix of type missing
#'
#' @keywords network
#'
#' @examples
#' network <- matrix(1:25, ncol=5)
#' rownames(network) <- c("a", "b", "c", "d", "e")
#' colnames(network) <- c("a", "b", "c", "d", "e")
#' genes.list <- c("a", "b",  "d", "e", "f", "g")
#' embed_network(network, genes.list)
#'
#' @export
#' 
embed_network <- function(network, genes.list, missing=0){
	n <- length(genes.list)
	target <- matrix(missing, ncol=n, nrow=n)
        dimnames(target) <- list(genes.list, genes.list)
	target_row_match <- na.omit(match(rownames(network), genes.list))
	target_col_match <- na.omit(match(colnames(network), genes.list))
	network_row_match <- na.omit(match(genes.list, rownames(network)))
	network_col_match <- na.omit(match(genes.list, colnames(network)))
	target[target_row_match,target_col_match] <- network[network_row_match, network_col_match]
	target
}
