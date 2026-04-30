#' Community Detection in KEGG-Derived Networks
#'
#' Applies graph-based community detection (Louvain or Leiden) to a
#' KEGG-derived network and returns an annotated edge list with cluster
#' membership for each source and target node. This supports module-based
#' interpretation and visualisation of large metabolic / signalling networks.
#'
#' @param network_table A matrix or data.frame with at least 3 columns,
#'   representing a directed KEGG-based network (typically from
#'   \code{MPN_keggNetwork()} or \code{MPN_crossSpeciesNetwork()}). The first three
#'   columns are interpreted as \code{source}, \code{target}, and
#'   \code{interaction_type}.
#' @param method Character. Community detection algorithm:
#'   \itemize{
#'     \item \code{"louvain"} — Louvain modularity optimisation.
#'     \item \code{"leiden"} — Leiden algorithm (improved modularity-based clustering).
#'   }
#' @param resolution Numeric (default: \code{0.1}). Resolution parameter
#'   for modularity-based clustering. Smaller values give fewer, larger
#'   communities; larger values give more, smaller communities.
#'
#' @details
#' Community detection is performed on the full graph (all nodes together),
#' using an undirected collapsed version of the input network.
#'
#' In the returned edge list, \code{source_cluster} and \code{target_cluster}
#' are endpoint annotations: they indicate the community membership of the
#' source node and target node of each edge, respectively.
#'
#' This edge-level annotation format is used to keep the output compatible
#' with downstream visualisation and filtering functions (e.g.
#' \code{MPN_viewClusterCy()} and \code{MPN_viewNetworkR()}).
#'
#' In MetaPathNet, these cluster assignments can also be used to define
#' modules for downstream enrichment or ego-network extraction.
#'
#' @return
#' A matrix with the same rows as \code{network_table} and the original
#' first three columns:
#' \describe{
#'   \item{source}{Source node identifier.}
#'   \item{target}{Target node identifier.}
#'   \item{interaction_type}{Interaction or edge annotation.}
#' }
#' plus two additional columns:
#' \describe{
#'   \item{source_cluster}{Community label of the source node.}
#'   \item{target_cluster}{Community label of the target node.}
#' }
#'
#' @examples
#' ## 1) Load the precomputed example network; see MPN_keggNetwork() for its construction
#' data(MetaPathNet_example_network)
#'
#' net_trp <- MetaPathNet_example_network
#'
#' ## 2) Run community detection on the example network
#' net_trp_clusters <- MPN_clusterNetwork(
#'   network_table = net_trp,
#'   method        = "leiden",
#'   resolution    = 0.01
#' )
#'
#' ## 3) (Optional) Visualise one cluster in R
#' ##    Here, source_cluster and target_cluster are node-cluster labels
#' ##    attached to each edge endpoint.
#' subnet_cluster1 <- net_trp_clusters[
#'   net_trp_clusters[, "source_cluster"] == 1 |
#'     net_trp_clusters[, "target_cluster"] == 1,
#' ]
#'
#' p_cluster1 <- MPN_viewNetworkR(
#'   network_table = subnet_cluster1,
#'   category      = TRUE,
#'   name          = TRUE,
#'   layout        = "fr",
#'   node_size     = 3,
#'   network_title = "Cluster 1 subnetwork in example network"
#' )
#'
#' print(p_cluster1)
#'
#' \dontrun{
#' ## Requires a running Cytoscape session for network layout and rendering
#' ## 4) Recommended: full cluster-aware visualisation in Cytoscape
#' MPN_viewClusterCy(
#'   cluster_matrix    = net_trp_clusters,
#'   label             = "all",   # convert KEGG IDs to biological names and label all nodes
#'   category          = TRUE,    # shape by node type (Compound / KO / Contextual)
#'   network_title     = "Example network clusters",
#'   collection_title  = "MetaPathNet_Examples"
#' )
#' }
#'
#' @seealso
#' \code{\link{MPN_keggNetwork}}, \code{\link{MPN_crossSpeciesNetwork}},
#' \code{\link{MPN_viewClusterCy}}, \code{\link{MPN_viewNetworkR}}
#'
#' @references
#' Blondel, V. D., Guillaume, J.-L., Lambiotte, R., and Lefebvre, E. (2008).
#' Fast unfolding of communities in large networks.
#' \emph{Journal of Statistical Mechanics: Theory and Experiment}, 2008(10), P10008.
#' \doi{10.1088/1742-5468/2008/10/P10008}
#'
#' Traag, V. A., Waltman, L., and van Eck, N. J. (2019).
#' From Louvain to Leiden: guaranteeing well-connected communities.
#' \emph{Scientific Reports}, 9, 5233.
#' \doi{10.1038/s41598-019-41695-z}
#'
#' @importFrom igraph graph_from_data_frame as.undirected
#'   cluster_louvain cluster_leiden membership
#' @export
MPN_clusterNetwork <- function(
    network_table,
    method     = c("louvain", "leiden"),
    resolution = 0.1
) {
  #==============================================================
  # Step 1 — Normalize and validate input
  #==============================================================
  if (!is.matrix(network_table)) {
    network_table <- as.matrix(network_table)
  }

  network_table <- unique(network_table)

  if (ncol(network_table) < 3L) {
    stop("network_table must be a matrix or data.frame with at least 3 columns: source, target, interaction_type.")
  }

  net_df <- as.data.frame(network_table, stringsAsFactors = FALSE)
  colnames(net_df)[seq_len(3)] <- c("source", "target", "interaction_type")

  method <- match.arg(method)

  if (!is.numeric(resolution) || length(resolution) != 1L || is.na(resolution) ||
      resolution <= 0) {
    stop("resolution must be a single positive numeric value.")
  }

  #==============================================================
  # Step 2 — Build undirected graph for community detection
  #==============================================================
  g_dir <- igraph::graph_from_data_frame(net_df[, seq_len(2)], directed = TRUE)
  g_u   <- igraph::as.undirected(g_dir, mode = "collapse")

  if (igraph::vcount(g_u) == 0) {
    stop("The network graph is empty after conversion. No community detection can be performed.")
  }

  #==============================================================
  # Step 3 — Run chosen clustering method
  #==============================================================
  if (method == "louvain") {
    cl <- igraph::cluster_louvain(g_u, resolution = resolution)
  } else { # "leiden"
    cl <- igraph::cluster_leiden(g_u, resolution_parameter = resolution)
  }

  membership_vec <- igraph::membership(cl)  # named vector: node -> cluster ID

  #==============================================================
  # Step 4 — Annotate network with source/target cluster IDs
  #==============================================================
  cluster_map <- setNames(as.integer(membership_vec), names(membership_vec))

  net_df$source_cluster <- cluster_map[net_df$source]
  net_df$target_cluster <- cluster_map[net_df$target]

  out <- as.matrix(net_df)
  out
}
