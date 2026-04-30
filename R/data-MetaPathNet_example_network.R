#' Example tryptophan-related MetaPathNet network
#'
#' A precomputed host-microbe example network constructed from human and
#' Escherichia coli KEGG pathways related to tryptophan metabolism and
#' immune/signaling context. This dataset is intended for runnable examples,
#' tests, and package demonstrations without requiring live KEGG queries.
#'
#' @format A character matrix with 1769 rows and 3 columns:
#' \describe{
#'   \item{source}{Source node identifier.}
#'   \item{target}{Target node identifier.}
#'   \item{interaction_type}{Interaction type between source and target.}
#' }
#'
#' @return A character matrix containing source nodes, target nodes, and
#'   interaction types.
#'
#' @source Constructed using \code{MPN_keggNetwork()} and
#' \code{MPN_mergeNetworks()} from KEGG pathways \code{hsa00380},
#' \code{hsa04060}, \code{hsa04630}, \code{hsa04064}, \code{hsa04660},
#' \code{hsa04659}, \code{eco00380}, and \code{eco00400}.
#'
#' @usage data(MetaPathNet_example_network)
#'
#' @keywords datasets
"MetaPathNet_example_network"
