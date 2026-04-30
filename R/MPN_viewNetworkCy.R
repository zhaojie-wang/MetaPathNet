#' Visualize a KEGG-Derived Network in Cytoscape via RCy3
#'
#' Imports a MetaPathNet-style network table (typically from \code{MPN_keggNetwork()}
#' or \code{MPN_shortestPaths()}) into Cytoscape and applies basic visual styling.
#' Nodes can be coloured by type (compound / KO / other). Node labels can be kept
#' as KEGG IDs (default) or converted to KEGG-derived names.
#'
#' @param network_table A 3-column matrix or data.frame representing the network
#'   (columns: \code{source}, \code{target}, \code{interaction_type}).
#' @param name Logical (default \code{FALSE}). If \code{TRUE}, node labels are
#'   converted from KEGG IDs to KEGG-derived names/symbols via \pkg{KEGGREST}.
#'   If \code{FALSE}, node labels remain KEGG IDs and no KEGG REST calls are made.
#' @param category Logical. If \code{TRUE}, nodes are coloured by type
#'   (\code{"Compound"}, \code{"KO"}, \code{"Other"}). If \code{FALSE}, all nodes
#'   use a default colour.
#' @param style_interaction Logical (default \code{FALSE}). If \code{TRUE}, basic
#'   edge styling is applied based on \code{interaction_type}:
#'   \code{"k_compound:reversible"} edges are drawn with parallel line style, and
#'   \code{"k_compound:irreversible"} edges receive a target arrow. Other
#'   interaction types are displayed with default Cytoscape edge styling.
#' @param node_compound Colour for compound nodes (default \code{"#9DC7DD"}).
#' @param node_gene Colour for KO/gene nodes (default \code{"#9ED17B"}).
#' @param network_title Character. Title of the network in Cytoscape.
#' @param collection_title Character. Cytoscape collection name.
#'
#' @details
#' Cytoscape (v\eqn{>=}3.9.0) must be running and reachable via \pkg{RCy3}
#' before calling this function.
#'
#' Reverse duplicated edges are collapsed for Cytoscape display when they
#' share the same metabolic interaction type
#'
#' When \code{name = TRUE}, KEGG IDs are resolved once per node. This can be slow
#' for large networks and depends on KEGG API response time.
#'
#' @return
#' Invisibly returns \code{NULL}. The function acts by side-effect in Cytoscape,
#' creating a styled network ready for manual or scripted exploration.
#'
#' @examples
#' ## 1) Load the precomputed example network; see MPN_keggNetwork() for its construction
#' data(MetaPathNet_example_network)
#'
#' net_trp_mixed <- MetaPathNet_example_network
#'
#' ## 2) Extract a shortest-path subnetwork (paths_trp_mixed)
#' source_nodes <- c("K00486", "K00453", "K00463",
#'                   "K01667", "K01696", "K01610", "K01695")
#' target_nodes <- c("cpd:C00078", "cpd:C00328", "cpd:C00780",
#'                   "cpd:C00463", "cpd:C00108", "cpd:C00458", "cpd:C00336")
#'
#' map_src <- MPN_findMappedNodes(source_nodes, net_trp_mixed)
#' map_tgt <- MPN_findMappedNodes(target_nodes, net_trp_mixed)
#'
#' paths_trp_mixed <- MPN_shortestPaths(
#'   network_table      = net_trp_mixed,
#'   source_nodes       = map_src$mapped_nodes,
#'   target_nodes       = map_tgt$mapped_nodes,
#'   mode               = "out",
#'   output             = "network_matrix",
#'   name               = FALSE,
#'   distance_threshold = 6
#' )
#'
#' \dontrun{
#' ## Requires a running Cytoscape session for network layout and rendering
#' ## Example 1 — Shortest-path subnetwork with KEGG name conversion + interaction styling
#' ## (k_compound:reversible shown as parallel lines; k_compound:irreversible with arrows)
#' MPN_viewNetworkCy(
#'   network_table     = paths_trp_mixed,
#'   name              = TRUE,
#'   category          = TRUE,
#'   style_interaction = TRUE,
#'   node_compound     = "#9DC7DD",
#'   node_gene         = "#9ED17B",
#'   network_title     = "Trp-immune shortest paths",
#'   collection_title  = "MetaPathNet_Examples"
#' )
#'
#' ## Example 2 — Full Trp network without name conversion (faster)
#' MPN_viewNetworkCy(
#'   network_table     = net_trp_mixed,
#'   name              = FALSE,
#'   category          = TRUE,
#'   style_interaction = FALSE,
#'   node_compound     = "#9DC7DD",
#'   node_gene         = "#9ED17B",
#'   network_title     = "Trp-immune network",
#'   collection_title  = "MetaPathNet_Examples"
#' )
#' }
#'
#' @seealso \code{\link{MPN_keggNetwork}}, \code{\link{MPN_viewNetworkR}},
#'   \code{\link{MPN_shortestPaths}}
#' @importFrom RCy3 cytoscapePing createNetworkFromDataFrames setNodeColorMapping
#' @importFrom RCy3 setNodeLabelMapping setNodePropertyBypass getTableColumns
#' @importFrom RCy3 setVisualPropertyDefault setEdgeLineStyleMapping
#' @importFrom RCy3 setEdgeTargetArrowShapeMapping
#' @importFrom KEGGREST keggGet
#' @export
MPN_viewNetworkCy <- function(
    network_table,
    name = FALSE,
    category = TRUE,
    style_interaction = FALSE,
    node_compound = "#9DC7DD",
    node_gene = "#9ED17B",
    network_title = "",
    collection_title = ""
) {
  #==============================================================
  # Step 1 — Check dependencies and Cytoscape connection
  #==============================================================
  cy_status <- tryCatch({
    suppressMessages(RCy3::cytoscapePing())
    TRUE
  }, error = function(e) FALSE)

  if (!cy_status) {
    message("Cytoscape is not running. Please start Cytoscape first.")
    return(invisible(NULL))
  } else {
    message("Cytoscape is connected")
  }

  #==============================================================
  # Step 2 — Normalize input and extract node set
  #==============================================================
  df_net <- as.data.frame(network_table, stringsAsFactors = FALSE)

  if (ncol(df_net) < 3) {
    stop("network_table must contain at least 3 columns: source, target, interaction_type.")
  }
  df_net <- df_net[, seq_len(3), drop = FALSE]
  colnames(df_net) <- c("source", "target", "interaction_type")

  all_nodes <- unique(c(df_net$source, df_net$target))

  #==============================================================
  # Step 3 — KEGG name conversion (optional)
  #==============================================================
  node_name <- all_nodes
  if (isTRUE(name)) {
    message("Converting KEGG IDs to names, may take time depending on network size.")
    for (i in seq_along(all_nodes)) {
      info <- tryCatch({
        KEGGREST::keggGet(all_nodes[i])
      }, error = function(e) list(list(NAME = all_nodes[i])))

      nm <- if (!is.null(info[[1]]$NAME)) info[[1]]$NAME[1] else all_nodes[i]

      if (grepl("^K", all_nodes[i])) {
        ## Remove bracketed KEGG annotation for EC numbers
        nm <- gsub("\\[.*?\\]", "", nm)
        nm <- trimws(nm)
      } else if (grepl("^cpd:", all_nodes[i])) {
        nm <- gsub(";", "", nm)
        nm <- trimws(nm)
      }

      node_name[i] <- nm
    }
  }

  #==============================================================
  # Step 4 — Define node category and color
  #==============================================================
  node_type  <- rep("Other", length(all_nodes))
  node_color <- rep("#95a5a6", length(all_nodes))

  if (isTRUE(category)) {
    node_type <- ifelse(
      grepl("^cpd:", all_nodes), "Compound",
      ifelse(grepl("^K", all_nodes), "KO", "Other")
    )
    node_color <- ifelse(
      node_type == "Compound", node_compound,
      ifelse(node_type == "KO", node_gene, "#95a5a6")
    )
  }

  #==============================================================
  # Step 5 — Prepare Cytoscape node and edge tables
  #==============================================================
  nodes_df <- data.frame(
    id    = all_nodes,
    label = if (isTRUE(name)) node_name else all_nodes,
    type  = node_type,
    color = node_color,
    stringsAsFactors = FALSE
  )

  # Collapse duplicated reverse edges for Cytoscape display:
  # if both A->B and B->A exist under the same metabolic interaction type,
  # keep only one edge. Different interaction types are kept separate.
  df_edges_vis <- df_net[, c("source", "target", "interaction_type")]

  collapse_types <- c(
    "k_compound:reversible",
    "custom:reversible",
    "k_compound:irreversible",
    "custom:irreversible"
  )
  collapse_idx <- which(df_edges_vis$interaction_type %in% collapse_types)

  if (length(collapse_idx) > 0L) {
    collapse_edges <- df_edges_vis[collapse_idx, , drop = FALSE]
    other_edges    <- df_edges_vis[-collapse_idx, , drop = FALSE]
    collapse_key <- paste(
      pmin(collapse_edges$source, collapse_edges$target),
      pmax(collapse_edges$source, collapse_edges$target),
      collapse_edges$interaction_type,
      sep = "__"
    )
    collapse_edges <- collapse_edges[!duplicated(collapse_key), , drop = FALSE]
    df_edges_vis <- rbind(other_edges, collapse_edges)
  }

  edges_df <- data.frame(
    source           = df_edges_vis$source,
    target           = df_edges_vis$target,
    interaction_type = df_edges_vis$interaction_type,
    stringsAsFactors = FALSE
  )

  #==============================================================
  # Step 6 — Create network in Cytoscape
  #==============================================================
  RCy3::createNetworkFromDataFrames(
    nodes      = nodes_df,
    edges      = edges_df,
    title      = network_title,
    collection = collection_title
  )

  #==============================================================
  # Step 7 — Apply Cytoscape visual styling
  #==============================================================
  if (isTRUE(category)) {
    RCy3::setNodeColorMapping(
      table.column        = "type",
      table.column.values = c("Compound", "KO", "Other"),
      colors              = c(node_compound, node_gene, "#95a5a6"),
      mapping.type        = "d"
    )
  } else {
    RCy3::setNodeColorMapping(
      table.column        = "type",
      table.column.values = unique(node_type),
      colors              = "#95a5a6",
      mapping.type        = "d"
    )
  }

  # Node label mapping
  RCy3::setNodeLabelMapping(table.column = "label")

  # Improve readability (increase node label font size)
  cy_nodes <- RCy3::getTableColumns("node")$name
  RCy3::setNodePropertyBypass(
    node.names       = cy_nodes,
    new.values       = 20,
    visual.property  = "NODE_LABEL_FONT_SIZE"
  )

  #==============================================================
  # Step 8 — Optional: edge styling by interaction_type
  #==============================================================
  if (isTRUE(style_interaction)) {

    ## Set default edge line type
    RCy3::setVisualPropertyDefault(
      list(visualProperty = "EDGE_LINE_TYPE", value = "SOLID")
    )

    ## Set default target arrow shape
    RCy3::setVisualPropertyDefault(
      list(visualProperty = "EDGE_TARGET_ARROW_SHAPE", value = "NONE")
    )

    ## Line style mapping: reversible vs everything else
    ## (default SOLID already covers "all others")
    RCy3::setEdgeLineStyleMapping(
      table.column        = "interaction_type",
      table.column.values = c("k_compound:reversible"),
      line.styles         = c("PARALLEL_LINES"),
      default.line.style  = "SOLID"
    )

    ## Arrow mapping: irreversible gets an ARROW target
    tryCatch(
      RCy3::setEdgeTargetArrowShapeMapping(
        table.column        = "interaction_type",
        table.column.values = c("k_compound:irreversible"),
        shapes              = c("ARROW"),
        default.shape       = "NONE"
      ),
      error = function(e) {
        message("Arrow mapping not available in this RCy3 version. Arrow mapping skipped.")
      }
    )
  }

  invisible(NULL)
}

