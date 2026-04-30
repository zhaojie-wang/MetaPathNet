#' Visualise a KEGG-Based Network Matrix in R
#'
#' Plots a KEGG-derived network matrix (e.g. from \code{MPN_keggNetwork()} or
#' \code{MPN_shortestPaths()}) directly in R using \pkg{ggraph}. Nodes can be
#' coloured by type (compound, KO, contextual) and labelled with KEGG IDs or
#' KEGG-derived names.
#'
#' @param network_table A network matrix with at least three columns: source, target,
#'   interaction type (e.g. output of \code{MPN_keggNetwork()},
#'   \code{MPN_crossSpeciesNetwork()}, or \code{MPN_shortestPaths()} with
#'   \code{output = "network_matrix"}).
#' @param category Logical. If \code{TRUE}, colour nodes by type
#'   (\code{"compound"}, \code{"KO"}, \code{"contextual"}); if \code{FALSE},
#'   all nodes share the same colour.
#' @param name Logical (default \code{FALSE}). If \code{TRUE}, node labels are
#'   retrieved from KEGG (compound names and KO symbols/names) via \pkg{KEGGREST}.
#'   If \code{FALSE}, node labels remain the raw KEGG IDs.
#' @param layout Character. Layout algorithm for \pkg{ggraph}
#'   (e.g. \code{"randomly"}, \code{"fr"}, \code{"kk"}).
#' @param node_size Numeric. Base size for nodes (and label text) in the plot.
#' @param node_compound Colour for compound nodes (default \code{"#9DC7DD"}).
#' @param node_gene Colour for KO nodes (default \code{"#9ED17B"}).
#' @param network_title Character. Optional plot title.
#' @details
#' For large networks (> 2000 edges), only the largest connected
#' component is displayed. If more than ~300 nodes are present, a message
#' suggests using \code{MPN_viewNetworkCy()} for Cytoscape-based visualisation.
#' When \code{name = TRUE}, KEGG REST queries are issued for each node, which
#' can be slow for large graphs.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' ## Load the precomputed example network; see MPN_keggNetwork() for its construction
#' data(MetaPathNet_example_network)
#'
#' ## Define source and target node sets, then keep only nodes mapped in the network
#' source_nodes <- c("K00453", "K00463", "K01667")
#' target_nodes <- c("cpd:C00078", "cpd:C00328", "cpd:C02700")
#'
#' map_src <- MPN_findMappedNodes(source_nodes, MetaPathNet_example_network)
#' map_tgt <- MPN_findMappedNodes(target_nodes, MetaPathNet_example_network)
#'
#' ## Extract a shortest-path subnetwork for visualisation
#' paths_trp <- MPN_shortestPaths(
#'   network_table      = MetaPathNet_example_network,
#'   source_nodes       = map_src$mapped_nodes,
#'   target_nodes       = map_tgt$mapped_nodes,
#'   mode               = "out",
#'   output             = "network_matrix",
#'   name               = FALSE,
#'   distance_threshold = 6
#' )
#'
#' ## Visualise the shortest-path subnetwork in R
#' p_trp <- MPN_viewNetworkR(
#'   network_table = paths_trp,
#'   category      = TRUE,
#'   name          = FALSE,
#'   layout        = "fr",
#'   node_size     = 5,
#'   network_title = "Tryptophan shortest-path subnetwork"
#' )
#' print(p_trp)
#'
#' @seealso
#' \code{\link{MPN_viewNetworkCy}}, \code{\link{MPN_crossSpeciesNetwork}},
#' \code{\link{MPN_keggNetwork}}, \code{\link{MPN_shortestPaths}}
#'
#' @importFrom tidygraph tbl_graph as_tbl_graph
#' @importFrom ggraph ggraph geom_edge_link geom_node_point geom_node_text
#' @importFrom dplyr case_when
#' @importFrom igraph components induced_subgraph
#' @importFrom ggplot2 aes ggtitle theme_void theme element_text margin scale_color_manual
#' @importFrom grid unit
#' @importFrom KEGGREST keggGet
#' @export
MPN_viewNetworkR <- function(
    network_table,
    category      = TRUE,
    name          = FALSE,
    layout        = "randomly",
    node_size     = 2,
    node_compound = "#9DC7DD",
    node_gene     = "#9ED17B",
    network_title = "MetaPathNet Network"
) {

  #==============================================================
  # STEP 1 — Build edge and node tables
  #==============================================================
  # Coerce input to a data frame for safe manipulation
  edges <- as.data.frame(network_table, stringsAsFactors = FALSE)
  # Standardise first three columns as edge from / to / interaction
  colnames(edges)[seq_len(3)] <- c("from", "to", "interaction_type")

  # Collect all unique node identifiers present in the edges
  node_ids <- unique(c(edges$from, edges$to))
  # Initialise node table with a single name column
  nodes <- data.frame(name = node_ids, stringsAsFactors = FALSE)

  # Classify nodes by identifier pattern: compound / KO / contextual
  nodes$type <- dplyr::case_when(
    grepl("^cpd:|^C\\d{5}", nodes$name, ignore.case = FALSE) ~ "compound",
    grepl("^K\\d+",        nodes$name, ignore.case = TRUE)  ~ "KO",
    TRUE                                                    ~ "contextual"
  )

  #==============================================================
  # STEP 2 — Retrieve biological labels (optional)
  #==============================================================
  # If name = FALSE, directly use raw node IDs as labels
  nodes$label <- nodes$name

  if (isTRUE(name)) {
    message(
      "Retrieving biological node annotations from the KEGG API; ",
      "this may take several minutes depending on network size."
    )

    # For each node, fetch a human-readable label from KEGG when possible
    nodes$label <- vapply(
      seq_len(nrow(nodes)),
      function(i) {
        id   <- nodes$name[i]
        type <- nodes$type[i]

        res <- tryCatch(KEGGREST::keggGet(id), error = function(e) list(NULL))
        if (length(res) == 0 || is.null(res[[1]])) return(id)

        if (type == "compound") {
          ## Use first compound name (before any ';') when available
          if (!is.null(res[[1]]$NAME)) {
            return(strsplit(res[[1]]$NAME[1], ";")[[1]][1])
          } else {
            return(id)
          }

        } else if (type == "KO") {
          ## Prefer KO symbol, then KO name, fall back to ID
          if (!is.null(res[[1]]$SYMBOL)) {
            return(strsplit(res[[1]]$SYMBOL, ",")[[1]][1])
          } else if (!is.null(res[[1]]$NAME)) {
            return(res[[1]]$NAME[1])
          } else {
            return(id)
          }

        } else {
          ## Contextual nodes keep their original ID
          return(id)
        }
      },
      character(1)
    )
  }

  #==============================================================
  # STEP 3 — Create graph object
  #==============================================================
  # Build a tbl_graph object compatible with ggraph
  g <- tbl_graph(nodes = nodes, edges = edges, directed = TRUE)

  #==============================================================
  # STEP 4 — Handle large networks (performance safeguards)
  #==============================================================
  # For very dense edge sets, keep only the largest weakly connected component
  if (nrow(edges) > 2000) {
    comp    <- igraph::components(g, mode = "weak")$membership
    biggest <- which.max(tabulate(comp))
    g       <- igraph::induced_subgraph(g, which(comp == biggest))
    g       <- as_tbl_graph(g)
  }

  # For many nodes, suggest Cytoscape-based visualisation
  if (nrow(nodes) > 300) {
    message(
      "For clearer visualization of large networks, ",
      "consider using MPN_viewNetworkCy() with Cytoscape."
    )
  }

  #==============================================================
  # STEP 5 — Customise visualisation
  #==============================================================
  # Draw edges with a neutral grey colour
  p <- suppressMessages(
    suppressWarnings(
      ggraph(g, layout = layout) +
        geom_edge_link(
          edge_alpha  = 0.7,
          edge_colour = "#839192",
          linetype    = "solid"
        )
    )
  )

  # Draw nodes, optionally coloured by category
  if (category) {
    p <- p +
      geom_node_point(aes(color = type), size = node_size) +
      scale_color_manual(
        values = c(
          compound   = node_compound,
          KO         = node_gene,
          contextual = "#839192"
        ),
        breaks = c("compound", "KO", "contextual"),
        labels = c("compound", "KEGG Ortholog (KO)", "contextual"),
        name   = "Node category"
      )
  } else {
    p <- p + geom_node_point(color = "#bfc9ca", size = node_size)
  }

  # Add node labels (KEGG IDs or KEGG names)
  p <- suppressMessages(
    suppressWarnings(
      p + geom_node_text(
        aes(label = label),
        repel        = TRUE,
        max.overlaps = Inf,
        vjust        = 1.5,
        size         = node_size + 1
      )
    )
  )

  #==============================================================
  # STEP 6 — Final styling and output
  #==============================================================
  # Apply minimal theme and centred title, then return ggplot object
  p <- p +
    ggtitle(network_title) +
    theme_void() +
    theme(
      plot.title  = element_text(
        hjust  = 0.5,
        face   = "bold",
        size   = 14,
        margin = margin(b = 10)
      ),
      plot.margin = unit(c(2, 2, 2, 2), "cm")
    )

  return(p)
}
