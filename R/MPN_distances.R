#' Calculate Shortest-Path Distances Between Node Sets in a KEGG-Based Network
#'
#' Compute shortest-path distances between selected source and target nodes
#' in a KEGG-based network produced by \code{MPN_keggNetwork()} (or related
#' functions). Source nodes are often KEGG KOs (upstream biological functions)
#' and target nodes are often KEGG compounds (downstream metabolites), but this
#' is a recommended usage pattern rather than a strict requirement.
#'
#' @param network_table A 3-column matrix or data.frame of edges
#'   (\code{source}, \code{target}, \code{interaction_type}), typically produced
#'   by \code{MPN_keggNetwork()}. Entries should use KEGG-style node identifiers.
#' @param mode Character, passed to \code{igraph::distances()}.
#'   Default \code{"out"}. One of:
#'   \describe{
#'     \item{\code{"out"}}{Distances follow edge direction.}
#'     \item{\code{"in"}}{Distances are computed against edge direction (reverse direction).}
#'     \item{\code{"all"}}{Distances ignore edge direction (treat the graph as undirected).}
#'   }
#' @param source_nodes Character vector of source nodes. Use \code{"all"}
#'   to include all non-compound nodes if \code{target_nodes = "all"}.
#'   Otherwise, \code{"all"} includes all nodes present in the network.
#' @param target_nodes Character vector of target nodes. Use \code{"all"}
#'   to include all compound-like nodes (\code{cpd:}, \code{dr:},
#'   \code{gl:}) if \code{source_nodes = "all"}. Otherwise, \code{"all"}
#'   includes all nodes present in the network.
#' @param name Logical (default \code{FALSE}). If \code{TRUE}, row and column
#'   labels are converted to KEGG names via \code{KEGGREST::keggGet()}.
#'
#' @return A numeric matrix of shortest-path distances. Rows correspond to
#'   source nodes and columns to target nodes.
#'
#' @details
#' The function builds a directed \pkg{igraph} object from the first two columns
#' of \code{network_table} and computes pairwise shortest-path distances using
#' \code{igraph::distances()}. For unweighted networks such as MetaPathNet edge
#' lists, distance lengths are computed with a breadth-first search strategy.
#'
#' Compound-like node IDs supplied as \code{"Cxxxxx"} or \code{"cpd:Cxxxxx"}
#' are treated equivalently and internally harmonised to \code{"cpd:Cxxxxx"}.
#'
#' If both \code{source_nodes = "all"} and \code{target_nodes = "all"}, the
#' function automatically assigns non-compound nodes (e.g. KOs/genes) as sources
#' and compound-like nodes (\code{cpd:}, \code{dr:}, \code{gl:}) as targets.
#'
#' If any node is present in both the mapped source and target sets, the function
#' stops and reports the overlapping node(s), because distance from a node to
#' itself is 0 and is not allowed in this workflow.
#'
#' Distances equal to \code{Inf} indicate unreachable source-target pairs under
#' the selected \code{mode}, i.e. no valid path exists between those nodes in
#' the current network direction setting.
#'
#' If \code{name = TRUE}, \code{KEGGREST::keggGet()} is used to retrieve KEGG
#' names. This may be slow for large node sets.
#'
#' @examples
#' ## 1) Load the precomputed example network; see MPN_keggNetwork() for its construction
#' data(MetaPathNet_example_network)
#'
#' ## 2) Define node sets of interest in this network
#' ##    - source_nodes: KOs representing enzymes/genes of interest in Trp/kynurenine biology
#' ##    - target_nodes: KEGG compounds representing metabolites of interest in the same context
#' ##    Note: MPN_distances can accept any node types, but source/target must be disjoint:
#' ##          if a node is in both sets, the function stops because self-distance is 0.
#' source_nodes <- c("K00486", "K00453", "K00463", "K09093", "K01667")
#' target_nodes <- c("cpd:C00078", "cpd:C00328", "cpd:C00780", "cpd:C00463")
#'
#' ## Keep only nodes that are present in the network (avoid "not mapped" errors)
#' map_src <- MPN_findMappedNodes(source_nodes, MetaPathNet_example_network)
#' map_tgt <- MPN_findMappedNodes(target_nodes, MetaPathNet_example_network)
#'
#' ## 3) Compute the shortest-path distance matrix
#' dist_trp <- MPN_distances(
#'   network_table = MetaPathNet_example_network,
#'   source_nodes  = map_src$mapped_nodes,
#'   target_nodes  = map_tgt$mapped_nodes,
#'   name          = FALSE
#' )
#'
#' dist_trp
#'
#' @seealso \code{\link{MPN_keggNetwork}}, \code{\link[igraph]{distances}}
#' @importFrom igraph graph_from_data_frame distances
#' @importFrom KEGGREST keggGet
#' @export
MPN_distances <- function(
    network_table,
    mode = c("out", "in", "all"),
    source_nodes,
    target_nodes,
    name = FALSE
) {
  #==============================================================
  # Step 1 — Basic input checks
  #==============================================================
  # Check mode argument
  mode <- match.arg(mode)

  # Check source_nodes / target_nodes arguments
  if (!is.character(source_nodes) || !length(source_nodes)) {
    stop("source_nodes must be a non-empty character vector.", call. = FALSE)
  }
  if (!is.character(target_nodes) || !length(target_nodes)) {
    stop("target_nodes must be a non-empty character vector.", call. = FALSE)
  }

  # Check name argument
  if (!is.logical(name) || length(name) != 1L || is.na(name)) {
    stop("name must be a single logical value (TRUE/FALSE).", call. = FALSE)
  }

  # Accept both matrix and data.frame input, then standardise to matrix
  if (is.data.frame(network_table)) {
    network_table <- as.matrix(network_table)
  }

  # Ensure network_table is a 3-column matrix (edge list)
  network_table <- unique(network_table)
  if (!is.matrix(network_table) || ncol(network_table) < 3) {
    stop("network_table must be a matrix with at least 3 columns")
  }

  # Work with character data
  storage.mode(network_table) <- "character"

  # Check that the network has nodes
  all_nodes <- unique(as.vector(network_table[, seq_len(2)]))
  if (length(all_nodes) == 0) {
    stop("No nodes found in the network_table.")
  }

  #==============================================================
  # Step 2 — Build graph and compute distance matrix
  #==============================================================
  g <- igraph::graph_from_data_frame(network_table[, seq_len(2)], directed = TRUE)
  distance_matrix <- igraph::distances(g, mode = mode)

  #==============================================================
  # Step 3 — Map requested source and target nodes
  #==============================================================
  # Helper: normalise compound-like IDs so "Cxxxxx" and "cpd:Cxxxxx" are treated equally
  normalize_node_ids <- function(x) {
    x <- unique(gsub(" ", "", as.character(x)))
    is_compound <- grepl("^(cpd:)?[Cc][0-9]{5}$", x)

    if (any(is_compound)) {
      core_ids <- toupper(gsub("^(cpd:)?", "", x[is_compound]))
      x[is_compound] <- paste0("cpd:", core_ids)
    }

    x
  }

  graph_nodes <- rownames(distance_matrix)

  # Identify compound-like vs non-compound nodes in the graph (for the special all/all case)
  compound_like_nodes <- graph_nodes[grepl("^(cpd:|dr:|gl:)", graph_nodes)]
  non_compound_nodes  <- setdiff(graph_nodes, compound_like_nodes)

  # Special case: if both are "all", use KO/gene-like nodes as sources and compound-like nodes as targets
  if (identical(source_nodes[1], "all") && identical(target_nodes[1], "all")) {
    mapped_sources   <- non_compound_nodes
    unmapped_sources <- NULL
    mapped_targets   <- compound_like_nodes
    unmapped_targets <- NULL

    if (length(mapped_sources) == 0) {
      stop("No non-compound nodes found in the network_table for source_nodes = 'all'.")
    }
    if (length(mapped_targets) == 0) {
      stop("No compound-like nodes found in the network_table for target_nodes = 'all'.")
    }

  } else {

    # Sources: default = all nodes in the network
    if (identical(source_nodes[1], "all")) {
      mapped_sources   <- graph_nodes
      unmapped_sources <- NULL
    } else {
      source_nodes     <- normalize_node_ids(source_nodes)
      mapped_sources   <- intersect(source_nodes, graph_nodes)
      unmapped_sources <- setdiff(source_nodes, graph_nodes)

      if (length(unmapped_sources) > 0) {
        message("Some source_nodes were not mapped onto the network: ",
                paste(unmapped_sources, collapse = ", "))
      }
      if (length(mapped_sources) == 0) {
        stop("None of the source_nodes was mapped onto the network")
      }
    }

    # Targets: default = all nodes in the network
    if (identical(target_nodes[1], "all")) {
      mapped_targets   <- graph_nodes
      unmapped_targets <- NULL
    } else {
      target_nodes     <- normalize_node_ids(target_nodes)
      mapped_targets   <- intersect(target_nodes, graph_nodes)
      unmapped_targets <- setdiff(target_nodes, graph_nodes)

      if (length(unmapped_targets) > 0) {
        message("Some target_nodes were not mapped onto the network: ",
                paste(unmapped_targets, collapse = ", "))
      }
      if (length(mapped_targets) == 0) {
        stop("None of the target_nodes was mapped onto the network")
      }
    }
  }

  # Stop if any node is present in both source and target sets (distance to itself = 0)
  overlap_nodes <- intersect(mapped_sources, mapped_targets)
  if (length(overlap_nodes) > 0) {
    stop(
      "Some nodes are present in both source_nodes and target_nodes. ",
      "Distance from a node to itself is 0, which is not allowed here. ",
      "Overlapping node(s): ", paste(overlap_nodes, collapse = ", "),
      call. = FALSE
    )
  }

  # Extract source → target submatrix
  submatrixGM <- distance_matrix[mapped_sources, mapped_targets, drop = FALSE]

  #==============================================================
  # Step 4 — Optional KEGG label conversion
  #==============================================================

  # Column labels: convert compound IDs to KEGG names
  if (isTRUE(name)) {
    target_ids <- colnames(submatrixGM)
    target_names <- vapply(target_ids, function(id) {
      Sys.sleep(0.05)
      res <- tryCatch(KEGGREST::keggGet(id), error = function(e) list(NULL))
      if (!is.null(res[[1]]) && !is.null(res[[1]]$NAME)) {
        gsub(";", "", res[[1]]$NAME[1])
      } else {
        id
      }
    }, character(1))
    colnames(submatrixGM) <- make.unique(target_names)
    message("Node name conversion uses the KEGG API and may take time for large lists.")
  }

  # Row labels: convert KO/gene IDs to KEGG names
  if (isTRUE(name)) {
    source_ids <- rownames(submatrixGM)
    source_names <- vapply(source_ids, function(id) {
      Sys.sleep(0.05)
      res <- tryCatch(KEGGREST::keggGet(id), error = function(e) list(NULL))

      if (!is.null(res[[1]]) && !is.null(res[[1]]$NAME)) {
        res[[1]]$NAME[1]
      } else {
        id
      }
    }, character(1))
    rownames(submatrixGM) <- make.unique(source_names)
  }

  return(submatrixGM)
}
