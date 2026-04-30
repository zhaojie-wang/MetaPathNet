#' Extract an Ego Subnetwork from a Metabolic / Signalling Network
#'
#' Builds an ego (neighbourhood) subnetwork around one or more seed nodes in a
#' directed network, up to a user-defined distance (\code{order}). In graph
#' terminology, \emph{ego} refers to the focal node(s) used as the starting
#' point(s), and the function returns the local network around them. This is
#' useful to inspect which nodes connect to a selected compound, KO, or other
#' node(s) within a chosen radius in a pathway- or organism-specific network.
#'
#' The default output is a 3-column MetaPathNet-style edge list
#' (\code{source}, \code{target}, \code{interaction_type}); alternatively, the
#' function can return only the unique node vector of the ego neighbourhood.
#'
#' @param network_table A three-column matrix or data.frame representing a
#'   directed network. The first three columns must correspond to
#'   \code{source}, \code{target}, and \code{interaction_type} (e.g. output of
#'   \code{MPN_keggNetwork()}).
#' @param nodes Character vector of one or more node IDs (e.g.
#'   \code{"cpd:C00078"}). Values must match the node names used in the network
#'   (case sensitive); leading/trailing whitespace is trimmed. Compound IDs
#'   supplied as \code{"Cxxxxx"} are internally harmonised to
#'   \code{"cpd:Cxxxxx"}.
#' @param order Integer (default \code{2}). Maximum graph distance (in number
#'   of edges) from the ego nodes to include in the neighbourhood. Must be a
#'   single positive integer \eqn{\ge 1}.
#' @param mode Character. Directionality of the search, passed to
#'   \code{igraph::ego()}:
#'   \itemize{
#'     \item \code{"out"} (default): include nodes reachable by following edge
#'       direction outward from the ego node(s).
#'     \item \code{"in"}: include nodes that can reach the ego node(s)
#'       (against edge direction).
#'     \item \code{"all"}: ignore edge direction for neighbourhood expansion
#'       (treat the graph as undirected for this search).
#'   }
#' @param mindist Integer (default \code{0}). Minimum graph distance from the
#'   ego nodes to include (also passed to \code{igraph::ego()}). Must satisfy
#'   \eqn{0 \le mindist \le order}.
#' @param output Character. Type of result:
#'   \itemize{
#'     \item \code{"matrix"} (default): 3-column edge list
#'       (\code{source}, \code{target}, \code{interaction_type}) for the
#'       induced subgraph on the ego neighbourhood.
#'     \item \code{"vector"}: character vector of unique node IDs present in
#'       the ego subnetwork.
#'   }
#'
#' @return
#' If \code{output = "matrix"}, a character matrix with columns
#' \code{source}, \code{target}, \code{interaction_type} representing the
#' induced subnetwork on all nodes within \code{[mindist, order]} of the ego
#' nodes (according to \code{mode}). If no edges are found, an empty
#' 3-column matrix with these columns is returned.
#'
#' If \code{output = "vector"}, a character vector of unique node IDs in the
#' ego subnetwork is returned.
#'
#' @details
#' Internally, the function:
#' \enumerate{
#'   \item checks and normalizes \code{network_table} and input \code{nodes};
#'   \item builds a directed \code{igraph} object from the first three columns;
#'   \item maps the requested ego nodes to existing vertex names (with a message
#'         for nodes not found);
#'   \item uses \code{igraph::ego()} to collect all vertices within
#'         \code{[mindist, order]} of the mapped egos;
#'   \item returns either the induced edge list (\code{"matrix"}) or the node
#'         set (\code{"vector"}).
#' }
#'
#' For large networks and high values of \code{order}, ego subnetworks can
#' become large and slower to generate. It is often practical to start with
#' \code{order = 1} or \code{2} and increase only if needed.
#'
#' @examples
#' ## Example — Ego subnetwork around tryptophan-related metabolites in a mixed hsa/eco network
#'
#' ## 1) Build a mixed host–microbe network for tryptophan-related metabolism
#' ## Load the precomputed example network; see MPN_keggNetwork() for its construction
#' data(MetaPathNet_example_network)
#'
#' ## Combined host–microbe network
#' net_trp_mixed <- MetaPathNet_example_network
#'
#' ## 2) Extract an ego subnetwork around key metabolites used as seed nodes
#' ## (tryptophan, serotonin, and indole are central compounds in the example storyline)
#' ego_trp_compounds <- MPN_egoNetwork(
#'   network_table = net_trp_mixed,
#'   nodes         = c(
#'     "cpd:C00078",  # L-tryptophan
#'     "cpd:C00780",  # serotonin
#'     "cpd:C00463"   # indole
#'   ),
#'   order   = 2,
#'   mode    = "out",
#'   mindist = 0,
#'   output  = "matrix"
#' )
#'
#' @references
#' Ego-network (ego graph) concept is standard in social/network analysis; see
#' Wasserman S, Faust K (1994). \emph{Social Network Analysis: Methods and Applications}.
#' Cambridge University Press.
#'
#' The neighbourhood extraction in this function is implemented with
#' \code{igraph::ego()} from the \pkg{igraph} package.
#'
#' @seealso
#' \code{\link[igraph]{ego}},
#' \code{\link[igraph]{induced_subgraph}},
#' \code{\link{MPN_enrichPathway}},
#' \code{\link{MPN_viewNetworkR}}
#'
#' @export
MPN_egoNetwork <- function(
    network_table,
    nodes,
    order   = 2,
    mode    = c("out", "all", "in"),
    mindist = 0,
    output  = c("matrix", "vector")
) {
  #==============================================================
  # Step 0 — Basic input checks
  #==============================================================
  # Accept both matrix and data.frame input, then standardise to matrix
  if (is.data.frame(network_table)) {
    network_table <- as.matrix(network_table)
  }

  # Ensure network_table is a proper matrix with at least 3 columns
  network_table <- unique(network_table)
  if (!is.matrix(network_table) || ncol(network_table) < 3) {
    stop("network_table must be a matrix with at least 3 columns")
  }

  # Work with character data
  storage.mode(network_table) <- "character"

  # Check nodes
  if (missing(nodes) || length(nodes) == 0) {
    stop("You must supply at least one node in 'nodes'.")
  }

  nodes <- trimws(as.character(nodes))
  nodes <- nodes[!is.na(nodes) & nodes != ""]
  nodes <- unique(nodes)

  if (length(nodes) == 0) {
    stop("No valid node IDs provided after removing NA/empty values.")
  }

  # Harmonise compound-like IDs so "Cxxxxx" and "cpd:Cxxxxx" are treated equally
  is_compound <- grepl("^(cpd:)?[Cc][0-9]{5}$", nodes)
  if (any(is_compound)) {
    core_ids <- toupper(gsub("^(cpd:)?", "", nodes[is_compound]))
    nodes[is_compound] <- paste0("cpd:", core_ids)
  }

  # Match mode and output arguments
  mode   <- match.arg(mode)
  output <- match.arg(output)

  #==============================================================
  # Step 1 — Validate order and mindist
  #==============================================================
  # order: single, non-NA, positive integer (>= 1)
  if (length(order) != 1 || is.na(order) || !is.numeric(order) ||
      order < 1 || order != as.integer(order)) {
    stop("'order' must be a single positive integer (>= 1).")
  }
  order <- as.integer(order)

  # mindist: single, non-NA, integer >= 0 and <= order
  if (length(mindist) != 1 || is.na(mindist) || !is.numeric(mindist) ||
      mindist < 0 || mindist != as.integer(mindist)) {
    stop("'mindist' must be a single integer >= 0.")
  }
  mindist <- as.integer(mindist)

  if (mindist > order) {
    stop("'mindist' cannot be greater than 'order'.")
  }

  #==============================================================
  # Step 2 — Build igraph object and map input nodes
  #==============================================================
  network_df <- as.data.frame(network_table[, seq_len(3), drop = FALSE], stringsAsFactors = FALSE)
  network_df <- unique(network_df)
  colnames(network_df) <- c("source", "target", "interaction_type")

  # Directed metabolic/signalling network
  g <- igraph::graph_from_data_frame(network_df, directed = TRUE)
  graph_names <- igraph::V(g)$name

  # Report nodes that are not present in the network
  missing_nodes <- setdiff(nodes, graph_names)
  if (length(missing_nodes) > 0) {
    message("These nodes are not present in the network: ",
            paste(missing_nodes, collapse = ", "))
  }

  # Restrict to nodes that actually exist in the graph
  mapped_nodes <- intersect(nodes, graph_names)
  if (length(mapped_nodes) == 0) {
    stop("None of the input nodes are present in the network.")
  }

  #==============================================================
  # Step 3 — Compute ego neighborhoods
  #==============================================================
  ego_nodes_list <- igraph::ego(
    graph   = g,
    order   = order,
    nodes   = mapped_nodes,
    mode    = mode,
    mindist = mindist
  )

  # Collect all nodes appearing in any ego neighborhood
  all_ego_nodes <- unique(unlist(lapply(ego_nodes_list, igraph::as_ids)))

  if (length(all_ego_nodes) == 0) {
    stop("No nodes found in ego network for the given parameters.")
  }

  #==============================================================
  # Step 4 — Return as node vector or induced edge subnetwork
  #==============================================================
  if (output == "vector") {
    # Just the set of nodes in the ego subnetwork
    return(all_ego_nodes)
  } else {
    # Induced subgraph on the ego nodes, returned as a MetaPathNet-style matrix
    sub_g <- igraph::induced_subgraph(g, vids = all_ego_nodes)
    df_out <- igraph::as_data_frame(sub_g, what = "edges")

    # If no edges are found, return an empty 3-column matrix
    if (nrow(df_out) == 0) {
      out <- matrix(character(0), ncol = 3)
      colnames(out) <- c("source", "target", "interaction_type")
      return(out)
    }

    # Ensure interaction_type column exists
    if (!"interaction_type" %in% colnames(df_out)) {
      df_out$interaction_type <- NA_character_
    }

    # Keep only the standard 3 columns and rename them
    df_out <- df_out[, c("from", "to", "interaction_type"), drop = FALSE]
    colnames(df_out) <- c("source", "target", "interaction_type")

    df_out <- unique(df_out)

    out <- as.matrix(df_out)
    rownames(out) <- NULL

    return(out)
  }
}
