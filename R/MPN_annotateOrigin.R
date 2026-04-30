#' Annotate KO Origin (Host vs Bacteria) and Optionally Export to Cytoscape
#'
#' Annotates nodes in a MetaPathNet-style network with inferred biological
#' origin, focusing on KEGG Orthology (KO) nodes. KO origin is classified as
#' human-only (\code{"hsa"}), bacteria-only (\code{"bacteria"}),
#' shared (\code{"hsa&bacteria"}), or \code{"other"} based on KEGG annotations.
#' Compound nodes are assigned the separate class \code{"compound"}.
#'
#' The primary output is a node-level annotation table (ID, type, display name,
#' origin). Cytoscape export is optional and can be enabled to visualise the
#' annotated network with origin-based node colours.
#'
#' KO origin is inferred from KEGG \code{GENES} annotations via
#' \code{KEGGREST::keggGet()}:
#' \itemize{
#'   \item \code{"hsa"} if the KO has human genes only.
#'   \item \code{"bacteria"} if the KO has bacterial genes only.
#'   \item \code{"hsa&bacteria"} if both human and bacterial genes are present.
#'   \item \code{"other"} if none of the above can be established.
#' }
#'
#' Bacterial organism codes are derived from \code{KEGGREST::keggList("organism")}
#' using the \code{"Prokaryotes;Bacteria"} phylogeny prefix. The
#' \code{bacteria_codes} argument can be used to restrict the bacterial origin
#' class to a user-defined organism panel (e.g. cohort-specific microbes).
#' If \code{bacteria_codes = NULL} (default), all KEGG bacterial organism codes
#' are used.
#'
#' @param network_table A MetaPathNet-style network with at least two columns:
#'   \code{source} and \code{target}. If a third column is present, it is
#'   interpreted as \code{interaction_type}.
#' @param bacteria_codes Optional character vector of KEGG organism codes.
#'   Restricts the \code{"bacteria"} and \code{"hsa&bacteria"} classes to this
#'   set. If \code{NULL} (default), all KEGG bacterial organism codes are used.
#' @param name Logical. If \code{TRUE}, node display names are converted from
#'   KEGG entries (KO: \code{SYMBOL} preferred, then \code{NAME}; non-KO:
#'   cleaned KEGG \code{NAME}). If \code{FALSE}, KEGG IDs are kept as labels.
#' @param export_cytoscape Logical (default \code{FALSE}). If \code{TRUE}, the
#'   annotated network is exported to Cytoscape and styled by origin class.
#'   If \code{FALSE}, only the annotation table is returned.
#' @param compound_color Hex colour for compound nodes.
#' @param hsa_color Hex colour for KO nodes classified as \code{"hsa"}.
#' @param micro_color Hex colour for KO nodes classified as \code{"bacteria"}.
#' @param hsa_micro_color Hex colour for KO nodes classified as \code{"hsa&bacteria"}.
#' @param other_color Hex colour for KO nodes classified as \code{"other"}.
#' @param network_title Character. Cytoscape network title (used only when
#'   \code{export_cytoscape = TRUE}).
#' @param collection_title Character. Cytoscape collection title (used only when
#'   \code{export_cytoscape = TRUE}).
#'
#' @details
#' Cytoscape (v\eqn{>=}3.9.0) must be running
#'
#' Reverse duplicated edges are collapsed for Cytoscape display when they
#' share the same metabolic interaction type
#'
#' @return
#' A node-level data frame with columns:
#' \itemize{
#'   \item \code{id}: KEGG node identifier.
#'   \item \code{type}: coarse node type (\code{"KO"}, \code{"Compound"}, \code{"Other"}).
#'   \item \code{name}: display label (KEGG ID or converted KEGG name/symbol).
#'   \item \code{origin}: inferred origin class (\code{"compound"}, \code{"hsa"},
#'         \code{"bacteria"}, \code{"hsa&bacteria"}, \code{"other"}).
#' }
#'
#' If \code{export_cytoscape = TRUE}, the function also creates and styles a
#' Cytoscape network as a side effect, and returns the same annotation table
#' invisibly.
#'
#' @examples
#' ## 1) Build a mixed human-E. coli Trp-focused network (used below)
#' ## Load the precomputed example network; see MPN_keggNetwork() for its construction
#' data(MetaPathNet_example_network)
#' net_trp_mixed <- MetaPathNet_example_network
#'
#' ## 2) Build a shortest-path subnetwork (paths_trp_mixed) for origin annotation
#' source_nodes <- c("K00486", "K00453", "K00463", "K01667", "K01696", "K01610", "K01695")
#' target_nodes <- c("cpd:C00078", "cpd:C00328", "cpd:C00780", "cpd:C00463",
#'                   "cpd:C00108", "cpd:C00458", "cpd:C00336")
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
#' ## 3) Return annotation table only (main use)
#' origin_paths_trp_mixed <- MPN_annotateOrigin(
#'   network_table   = paths_trp_mixed,
#'   bacteria_codes  = NULL,   # NULL = use all KEGG bacterial organism codes
#'   name            = FALSE,   # convert IDs to KEGG names/symbols in the output table
#'   export_cytoscape = FALSE
#' )
#'
#' @seealso
#' \code{\link{MPN_annotateKoClass}} for KO functional class annotation, and
#' \code{\link{MPN_viewNetworkCy}} for generic Cytoscape rendering.
#'
#' @importFrom KEGGREST keggGet keggList
#' @importFrom RCy3 cytoscapePing createNetworkFromDataFrames
#'   setNodeColorMapping setNodeLabelMapping
#' @export
MPN_annotateOrigin <- function(
    network_table,
    bacteria_codes = NULL,
    name  = FALSE,
    export_cytoscape = FALSE,
    compound_color   = "#9DC7DD",
    hsa_color        = "#C45745",
    micro_color      = "#9ED17B",
    hsa_micro_color  = "#B57F60",
    other_color      = "#D1C2C2",
    network_title    = "",
    collection_title = ""
) {
  #==============================================================
  # STEP 1 - Basic input checks and matrix normalisation
  #==============================================================

  # Ensure the input network matrix is provided
  if (missing(network_table)) {
    stop("Argument 'network_table' (MetaPathNet-style 3-column network) is required.")
  }

  # Validate 'name' argument
  if (!is.logical(name) || length(name) != 1 || is.na(name)) {
    stop("Argument 'name' must be TRUE or FALSE.")
  }

  # Convert matrix to data frame and check structure
  net_df <- as.data.frame(network_table, stringsAsFactors = FALSE)
  if (ncol(net_df) < 2) {
    stop("Input 'network_table' must have at least 2 columns (source, target).")
  }

  # Standardise column names and add missing 'interaction_type' if needed
  colnames(net_df)[seq_len(2)] <- c("source", "target")
  if (ncol(net_df) >= 3) {
    colnames(net_df)[3] <- "interaction_type"
  } else {
    net_df$interaction_type <- NA_character_
  }

  #==============================================================
  # STEP 2 - Classify node types (Compound, KO, Other)
  #==============================================================

  # Extract all unique node IDs from the network
  node_ids <- unique(c(net_df$source, net_df$target))

  # Assign coarse node types based on identifier patterns
  node_type <- character(length(node_ids))
  node_type[grepl("^(cpd:|gl:|dr:)", node_ids)] <- "Compound"
  node_type[grepl("^K\\d{5}$", node_ids)]       <- "KO"
  node_type[node_type == ""]                    <- "Other"

  #==============================================================
  # STEP 3 - If Cytoscape export is requested, check RCy3 connection
  #==============================================================

  if (isTRUE(export_cytoscape)) {
    ok <- tryCatch(
      {
        RCy3::cytoscapePing()
        TRUE
      },
      error = function(e) FALSE
    )

    if (!ok) {
      stop("Cytoscape is not running or RCy3 cannot connect.")
    }
  }

  #==============================================================
  # STEP 4 - Convert node IDs to KEGG names (optional)
  #==============================================================

  name_map_nonko <- NULL
  if (isTRUE(name)) {
    non_ko_idx <- node_type != "KO"
    non_ko_ids <- node_ids[non_ko_idx]

    if (length(non_ko_ids) > 0) {
      name_map_nonko <- setNames(rep(NA_character_, length(non_ko_ids)), non_ko_ids)

      for (i in seq_along(non_ko_ids)) {
        id_i <- non_ko_ids[i]

        ## Pause every 10 queries to reduce KEGG API pressure
        if (i %% 10 == 0) {
          Sys.sleep(1)
        }

        # Query KEGG REST for compound/drug ID
        info_i <- tryCatch(KEGGREST::keggGet(id_i), error = function(e) NULL)
        if (is.null(info_i) || length(info_i) == 0 ||
            is.null(info_i[[1]]$NAME) || length(info_i[[1]]$NAME) == 0) {
          next
        }

        # Use first name field; strip trailing semicolon if present
        nm <- info_i[[1]]$NAME[1]
        if (grepl("^cpd:", id_i)) nm <- gsub(";", "", nm)
        name_map_nonko[id_i] <- trimws(nm)
      }
    }
  }

  #==============================================================
  # STEP 5 - Parse KO nodes: labels, associated organisms, origin class
  #==============================================================

  kos <- unique(node_ids[node_type == "KO"])

  # Prepare a data frame to store KO-level metadata
  ko_df <- data.frame(
    KO             = kos,
    organism_codes = I(vector("list", length(kos))),
    clean_code     = NA_character_,
    stringsAsFactors = FALSE
  )

  # Prepare optional KO name mapping (SYMBOL or NAME)
  ko_name_map <- if (isTRUE(name) && length(kos) > 0)
    setNames(rep(NA_character_, length(kos)), kos) else NULL

  if (length(kos) > 0) {
    for (i in seq_along(kos)) {
      ko_id <- kos[i]

      ## Pause every 10 queries to reduce KEGG API pressure
      if (i %% 10 == 0) {
        Sys.sleep(1)
      }

      # Query KEGG for KO entry
      res <- tryCatch(KEGGREST::keggGet(ko_id), error = function(e) NULL)
      if (is.null(res) || length(res) == 0) next
      entry <- res[[1]]

      # 5.1 - KO label (SYMBOL preferred, NAME fallback)
      if (isTRUE(name)) {
        used_label <- NA_character_

        if (!is.null(entry$SYMBOL) && length(entry$SYMBOL) > 0) {
          used_label <- trimws(strsplit(entry$SYMBOL[1], ",", fixed = TRUE)[[1]][1])
        }

        if (is.na(used_label) || nchar(used_label) == 0) {
          if (!is.null(entry$NAME) && length(entry$NAME) > 0) {
            used_label <- trimws(gsub("\\[.*?\\]", "", entry$NAME[1]))
          }
        }

        ko_name_map[ko_id] <- used_label
      }

      # 5.2 - Extract organism codes from KO -> GENES
      genes_field <- entry$GENES
      if (is.null(genes_field) || length(genes_field) == 0) next

      genes_vec <- if (is.list(genes_field)) {
        unlist(genes_field, use.names = TRUE)
      } else {
        genes_field
      }

      org_codes <- names(genes_vec)
      if (is.null(org_codes) || !any(nzchar(org_codes))) {
        org_codes <- sub("^([^:]+):.*$", "\\1", genes_vec)
      }

      org_codes <- unique(tolower(org_codes[nzchar(org_codes)]))
      ko_df$organism_codes[[i]] <- org_codes
      ko_df$clean_code[i] <- paste(org_codes, collapse = ",")
    }
  }

  #==============================================================
  # STEP 6 - Classify KO origin using KEGG organism table
  #==============================================================

  # Identify bacterial organisms from KEGG taxonomy
  info <- tryCatch(
    KEGGREST::keggList("organism"),
    error = function(e) NULL
  )

  if (is.null(info)) {
    stop("Failed to retrieve KEGG organism information.")
  }

  info <- as.data.frame(info, stringsAsFactors = FALSE)

  info <- info[, seq_len(4), drop = FALSE]
  colnames(info) <- c("T_number", "organism", "species", "phylogeny")

  parts <- strsplit(info$phylogeny, ";", fixed = TRUE)
  info$phylo2 <- vapply(
    parts,
    function(x) {
      x <- trimws(x)
      if (length(x) >= 2) {
        paste(x[seq_len(2)], collapse = ";")
      } else {
        NA_character_
      }
    },
    FUN.VALUE = character(1)
  )

  organism_bacteria <- tolower(info$organism[info$phylo2 == "Prokaryotes;Bacteria"])
  # Optional restriction to user-provided bacterial codes
  if (!is.null(bacteria_codes)) {
    organism_bacteria <- intersect(tolower(bacteria_codes), organism_bacteria)
  }

  # Assign origin class for each KO: hsa, bacteria, hsa&bacteria, other
  ko_df$origin <- vapply(
    ko_df$clean_code,
    function(cc) {
      if (is.na(cc) || nchar(cc) == 0) return("other")
      codes <- strsplit(cc, ",", fixed = TRUE)[[1]]
      has_hsa   <- "hsa" %in% codes
      has_micro <- any(codes %in% organism_bacteria)

      if (has_hsa && has_micro) "hsa&bacteria"
      else if (has_hsa) "hsa"
      else if (has_micro) "bacteria"
      else "other"
    },
    FUN.VALUE = character(1)
  )
  ko_origin_map <- setNames(ko_df$origin, ko_df$KO)

  #==============================================================
  # STEP 7 - Build node table and assign final labels
  #==============================================================

  nodes_df <- data.frame(
    id     = node_ids,
    type   = node_type,
    origin = "other",
    stringsAsFactors = FALSE
  )

  # Assign origin class for compound and KO nodes
  nodes_df$origin[nodes_df$type == "Compound"] <- "compound"
  is_ko <- nodes_df$type == "KO"
  if (any(is_ko)) {
    nodes_df$origin[is_ko] <- ko_origin_map[nodes_df$id[is_ko]]
  }

  # Assign node labels based on 'name' mode
  nodes_df$label <- nodes_df$id
  if (isTRUE(name)) {
    if (!is.null(name_map_nonko)) {
      idx <- nodes_df$type != "KO"
      nodes_df$label[idx] <- ifelse(
        is.na(name_map_nonko[nodes_df$id[idx]]),
        nodes_df$id[idx],
        name_map_nonko[nodes_df$id[idx]]
      )
    }

    if (!is.null(ko_name_map)) {
      nodes_df$label[is_ko] <- ifelse(
        is.na(ko_name_map[nodes_df$id[is_ko]]),
        nodes_df$id[is_ko],
        ko_name_map[nodes_df$id[is_ko]]
      )
    }
  }

  #==============================================================
  # STEP 8 - Return annotation table, optionally export to Cytoscape
  #==============================================================

  out_df <- nodes_df[, c("id", "type", "label", "origin"), drop = FALSE]
  colnames(out_df)[colnames(out_df) == "label"] <- "name"

  if (!isTRUE(export_cytoscape)) {
    return(out_df)
  }

  # Build edge table from input matrix
  edges_df <- data.frame(
    source           = net_df$source,
    target           = net_df$target,
    interaction_type = net_df$interaction_type,
    stringsAsFactors = FALSE
  )

  # Collapse duplicated reverse edges for Cytoscape display:
  # if both A->B and B->A exist under the same metabolic interaction type,
  # keep only one edge. Different interaction types are kept separate.
  collapse_types <- c(
    "k_compound:reversible",
    "custom:reversible",
    "k_compound:irreversible",
    "custom:irreversible"
  )
  collapse_idx <- which(edges_df$interaction_type %in% collapse_types)

  if (length(collapse_idx) > 0L) {
    collapse_edges <- edges_df[collapse_idx, , drop = FALSE]
    other_edges    <- edges_df[-collapse_idx, , drop = FALSE]
    collapse_key <- paste(
      pmin(collapse_edges$source, collapse_edges$target),
      pmax(collapse_edges$source, collapse_edges$target),
      collapse_edges$interaction_type,
      sep = "__"
    )

    collapse_edges <- collapse_edges[!duplicated(collapse_key), , drop = FALSE]
    edges_df <- rbind(other_edges, collapse_edges)
  }

  # Create Cytoscape network via RCy3
  RCy3::createNetworkFromDataFrames(
    nodes      = nodes_df,
    edges      = edges_df,
    title      = network_title,
    collection = collection_title
  )

  # Define color mapping by origin class
  origin_colors <- c(
    compound        = compound_color,
    hsa             = hsa_color,
    bacteria        = micro_color,
    "hsa&bacteria"  = hsa_micro_color,
    other           = other_color
  )

  # Apply visual styling to Cytoscape network
  RCy3::setNodeColorMapping(
    table.column        = "origin",
    table.column.values = names(origin_colors),
    colors              = unname(origin_colors),
    mapping.type        = "d"
  )

  RCy3::setNodeLabelMapping(table.column = "label")

  invisible(out_df)
}
