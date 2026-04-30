#' Visualise Clustered MetaPathNet Networks in Cytoscape
#'
#' Sends a cluster-annotated MetaPathNet network (typically the output of
#' \code{MPN_clusterNetwork()}) to Cytoscape for interactive exploration.
#' Nodes are coloured by community (cluster), optionally shaped by category
#' (compound / KO / contextual), and can be labelled with KEGG-derived names.
#'
#' Cytoscape (version \strong{>= 3.9.0}) must be running and reachable via
#' \pkg{RCy3} before calling this function.
#'
#' @param cluster_matrix A matrix or data.frame with at least five columns:
#'   \code{source}, \code{target}, \code{interaction_type},
#'   \code{source_cluster}, and \code{target_cluster}. Typically the direct
#'   output of \code{MPN_clusterNetwork()}.
#' @param label Character; controls KEGG-based node labelling. One of:
#'   \itemize{
#'     \item \code{"none"} – no labels are shown on the Cytoscape network.
#'     \item \code{"all"} – convert KEGG IDs to biological names and label all nodes.
#'     \item \code{"non_background"} – convert and label only nodes in non-background
#'           clusters. The background cluster is defined as the largest
#'           community and is left unlabelled to reduce KEGG REST calls and visual clutter.
#'   }
#' @param category Logical. If \code{TRUE}, node shapes are mapped to
#'   basic categories: \code{"Compound"}, \code{"KO"}, \code{"Contextual"}.
#'   If \code{FALSE}, all nodes use the same shape (ellipse).
#' @param network_title Character. Cytoscape network title. If empty or
#'   \code{NULL}, defaults to \code{"cluster_network"}.
#' @param collection_title Character. Cytoscape collection title. If empty or
#'   \code{NULL}, defaults to \code{"cluster_collection"}.
#'
#' @details
#' The function reconstructs node-level cluster membership from the
#' \code{source_cluster} and \code{target_cluster} columns, identifies the
#' largest community as a background cluster (for \code{label = "non_background"}),
#' optionally converts KEGG IDs to names via \pkg{KEGGREST}, and then creates a
#' Cytoscape network via \pkg{RCy3} with:
#' \itemize{
#'   \item colour mapping by cluster,
#'   \item shape mapping by node type (if \code{category = TRUE}),
#'   \item node labels from the chosen \code{label} mode.
#' }
#'
#' KEGG name conversion can be slow on large networks, especially with
#' \code{label = "all"}. For big MetaPathNet networks, \code{label = "non_background"}
#' is often a good compromise.
#'
#' @return
#' Invisibly returns \code{NULL}. The function acts via side effects in
#' Cytoscape, creating and styling a new network.
#'
#' @seealso
#' \code{\link{MPN_clusterNetwork}} for community detection, and
#' \code{\link{MPN_viewNetworkCy}} for non-clustered KEGG network visualisation.
#'
#' @examples
#' ## 1) Load the precomputed example network; see MPN_keggNetwork() for its construction
#' data(MetaPathNet_example_network)
#'
#' net_trp <- MetaPathNet_example_network
#'
#' ## 2) Detect communities on the full network graph
#' net_trp_clusters <- MPN_clusterNetwork(
#'   network_table = net_trp,
#'   method        = "leiden",
#'   resolution    = 0.01
#' )
#'
#' \dontrun{
#' ## Requires a running Cytoscape session for network layout and rendering
#' ## 3) Recommended: full cluster-aware visualisation in Cytoscape
#' MPN_viewClusterCy(
#'   cluster_matrix    = net_trp_clusters,
#'   label             = "all",   # convert KEGG IDs to biological names and label all nodes
#'   category          = TRUE,    # shape by node type (Compound / KO / Contextual)
#'   network_title     = "Example network clusters",
#'   collection_title  = "MetaPathNet_Examples"
#' )
#'
#' ## 4) Faster option for larger networks: label only non-background clusters
#' MPN_viewClusterCy(
#'   cluster_matrix    = net_trp_clusters,
#'   label             = "non_background",
#'   category          = TRUE,
#'   network_title     = "Example network clusters (foreground labels)",
#'   collection_title  = "MetaPathNet_Examples"
#' )
#' }
#'
#' @importFrom RCy3 cytoscapePing createNetworkFromDataFrames setNodeColorMapping
#' @importFrom RCy3 setNodeShapeMapping setNodeLabelMapping setNodePropertyBypass getTableColumns
#' @importFrom KEGGREST keggGet
#' @export
MPN_viewClusterCy <- function(
    cluster_matrix,
    label            = c("none", "all", "non_background"),
    category         = TRUE,
    network_title    = "",
    collection_title = ""
) {
  #==============================================================
  # Step 1 — Normalize input and basic checks
  #==============================================================
  cm <- as.data.frame(cluster_matrix, stringsAsFactors = FALSE)

  if (ncol(cm) < 5L) {
    stop("cluster_matrix must have at least 5 columns: ",
         "source, target, interaction_type, source_cluster, target_cluster.")
  }

  colnames(cm)[seq_len(5)] <- c(
    "source",
    "target",
    "interaction_type",
    "source_cluster",
    "target_cluster"
  )

  label_mode <- match.arg(label)

  cy_status <- tryCatch({
    suppressMessages(RCy3::cytoscapePing())
    TRUE
  }, error = function(e) FALSE)

  if (!cy_status) {
    message("Cytoscape is not running or not reachable. Please start Cytoscape first.")
    return(invisible(NULL))
  } else {
    message("Cytoscape is connected.")
  }

  #==============================================================
  # Step 2 — Build node-level table with cluster membership
  #==============================================================
  all_nodes <- unique(c(cm$source, cm$target))

  # Recover one cluster label per node from the edge-level cluster columns
  source_cluster_df <- unique(cm[, c("source", "source_cluster"), drop = FALSE])
  target_cluster_df <- unique(cm[, c("target", "target_cluster"), drop = FALSE])

  node_cluster <- source_cluster_df$source_cluster[
    match(all_nodes, source_cluster_df$source)
  ]

  fill_idx <- is.na(node_cluster)
  node_cluster[fill_idx] <- target_cluster_df$target_cluster[
    match(all_nodes[fill_idx], target_cluster_df$target)
  ]

  # Cluster IDs arrive as character values from the edge-list matrix; convert back to integers
  node_cluster <- as.integer(node_cluster)
  names(node_cluster) <- all_nodes

  # Identify largest cluster = considered "background"
  valid_clusters <- node_cluster[!is.na(node_cluster)]
  bg_cluster <- if (length(valid_clusters)) as.integer(names(sort(table(valid_clusters), TRUE))[1]) else NA_integer_

  # Basic node type classification
  node_type <- ifelse(grepl("^cpd:", all_nodes), "Compound",
                      ifelse(grepl("^K", all_nodes), "KO", "Contextual"))

  #==============================================================
  # Step 3 — KEGG name conversion according to 'label' mode
  #==============================================================
  node_name <- all_nodes  # default: keep KEGG IDs / no labels

  if (label_mode != "none") {

    # Select nodes to convert depending on mode
    ids_to_convert <- if (label_mode == "all") {
      all_nodes
    } else {
      # Convert only non-background clusters
      names(node_cluster)[!is.na(node_cluster) & node_cluster != bg_cluster]
    }

    ids_to_convert <- unique(ids_to_convert)

    if (length(ids_to_convert) > 0L) {
      message("Converting ", length(ids_to_convert),
              " nodes to biological names. This step may take some time.")

      # name_map stores KEGG → biological name
      name_map <- setNames(all_nodes, all_nodes)

      for (id in ids_to_convert) {
        info <- tryCatch(KEGGREST::keggGet(id), error = function(e) NULL)
        new_label <- id  # default fallback

        # Compound: use NAME[1] if available
        if (!is.null(info) && grepl("^cpd:", id) && !is.null(info[[1]]$NAME)) {
          new_label <- trimws(gsub(";", "", info[[1]]$NAME[1]))
        }

        # KO: prefer SYMBOL, else NAME[1]
        else if (!is.null(info) && grepl("^K", id)) {
          entry <- info[[1]]
          if (!is.null(entry$SYMBOL)) {
            new_label <- trimws(strsplit(entry$SYMBOL, ",")[[1]][1])
          } else if (!is.null(entry$NAME)) {
            new_label <- trimws(gsub("\\[.*?\\]", "", entry$NAME[1]))
          }
        }

        name_map[id] <- new_label
      }

      node_name <- name_map[all_nodes]
    }
  }

  #==============================================================
  # Step 4 — Compose Cytoscape node and edge tables
  #==============================================================
  nodes_df <- data.frame(
    id      = all_nodes,
    label   = if (identical(label_mode, "none")) "" else node_name,
    type    = node_type,
    cluster = node_cluster,
    stringsAsFactors = FALSE
  )

  # Collapse duplicated reverse edges for Cytoscape display:
  # if both A->B and B->A exist under the same metabolic interaction type,
  # keep only one edge. Different interaction types are kept separate.
  cm_edges_vis <- cm[, c("source", "target", "interaction_type")]

  collapse_types <- c(
    "k_compound:reversible",
    "custom:reversible",
    "k_compound:irreversible",
    "custom:irreversible"
  )

  collapse_idx <- which(cm_edges_vis$interaction_type %in% collapse_types)

  if (length(collapse_idx) > 0L) {
    collapse_edges <- cm_edges_vis[collapse_idx, , drop = FALSE]
    other_edges    <- cm_edges_vis[-collapse_idx, , drop = FALSE]
    collapse_key <- paste(
      pmin(collapse_edges$source, collapse_edges$target),
      pmax(collapse_edges$source, collapse_edges$target),
      collapse_edges$interaction_type,
      sep = "__"
    )
    collapse_edges <- collapse_edges[!duplicated(collapse_key), , drop = FALSE]
    cm_edges_vis <- rbind(other_edges, collapse_edges)
  }

  edges_df <- data.frame(
    source      = cm_edges_vis$source,
    target      = cm_edges_vis$target,
    interaction = cm_edges_vis$interaction_type,
    stringsAsFactors = FALSE
  )

  #==============================================================
  # Step 5 — Create network in Cytoscape
  #==============================================================
  # Fallback to generic names if user leaves them empty
  if (is.null(network_title) || identical(network_title, "")) {
    network_title <- "cluster_network"
  }
  if (is.null(collection_title) || identical(collection_title, "")) {
    collection_title <- "cluster_collection"
  }

  message("Creating Cytoscape network. This may take some time for large networks.")

  suppressMessages(
    RCy3::createNetworkFromDataFrames(
      nodes      = nodes_df,
      edges      = edges_df,
      title      = network_title,
      collection = collection_title
    )
  )

  #==============================================================
  # Step 6 — Visual mapping: color by cluster, shape by type
  #==============================================================
  cluster_vals <- sort(unique(nodes_df$cluster))
  cluster_vals <- cluster_vals[!is.na(cluster_vals)]

  ## Color mapping by cluster
  if (length(cluster_vals) > 0L) {
    cluster_colors <- grDevices::rainbow(length(cluster_vals))

    message("Customizing node colors according to cluster membership.")
    suppressMessages(
      RCy3::setNodeColorMapping(
        table.column        = "cluster",
        table.column.values = cluster_vals,
        colors              = cluster_colors,
        mapping.type        = "d"
      )
    )
  }

  ## Shape mapping by category (node_type)
  if (isTRUE(category)) {
    message("Customizing node shapes according to node category (Compound / KO / Contextual).")
    suppressMessages(
      RCy3::setNodeShapeMapping(
        table.column        = "type",
        table.column.values = c("Compound", "KO", "Contextual"),
        shapes              = c("ELLIPSE", "DIAMOND", "HEXAGON")
      )
    )
  } else {
    message("Applying uniform node shape (ELLIPSE) for all node types.")
    suppressMessages(
      RCy3::setNodeShapeMapping(
        table.column        = "type",
        table.column.values = unique(nodes_df$type),
        shapes              = rep("ELLIPSE", length(unique(nodes_df$type)))
      )
    )
  }

  ## Node labels
  message("Labelling nodes according to the 'label' column.")
  suppressMessages(
    RCy3::setNodeLabelMapping(table.column = "label")
  )

  ## Increase label font size for readability
  suppressMessages({
    cy_nodes <- RCy3::getTableColumns("node")$name
    RCy3::setNodePropertyBypass(
      node.names      = cy_nodes,
      new.values      = 20,
      visual.property = "NODE_LABEL_FONT_SIZE"
    )
  })

  invisible(NULL)
}
