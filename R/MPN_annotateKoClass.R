#' Annotate KO functional class (metabolic vs signalling)
#'
#' Annotates KEGG Orthology (KO) nodes in a MetaPathNet-style network by
#' pathway context: \code{"metabolic"}, \code{"signaling"}, or
#' \code{"metabolic&signaling"}. Classification is derived from KEGG BRITE
#' pathway hierarchy and KO-to-pathway links. The function returns a KO-level
#' annotation table (compound nodes excluded) and can optionally export a
#' Cytoscape network with KO class colouring.
#'
#' Cytoscape (version \strong{>= 3.9.0}) must be open and reachable via
#' \pkg{RCy3} only when \code{export_cytoscape = TRUE}.
#'
#' @param network_table A MetaPathNet-style network with at least two columns:
#'   \code{source} and \code{target}. If a third column is present, it is
#'   interpreted as \code{interaction_type}; otherwise it is set to
#'   \code{NA_character_}.
#' @param name Logical. If \code{TRUE}, convert node IDs to KEGG-derived labels
#'   (KO: \code{SYMBOL} preferred, then \code{NAME}; non-KO: cleaned
#'   \code{NAME}). If \code{FALSE} (default), keep original node IDs.
#' @param export_cytoscape Logical. If \code{FALSE} (default), return only the
#'   annotation data frame. If \code{TRUE}, also create and style a Cytoscape
#'   network using the annotated node table.
#' @param compound_color Hex colour for compound nodes (used only when
#'   \code{export_cytoscape = TRUE}).
#' @param metabolic_color Hex colour for KO nodes classified as
#'   \code{"metabolic"} (used only when \code{export_cytoscape = TRUE}).
#' @param signaling_color Hex colour for KO nodes classified as
#'   \code{"signaling"} (used only when \code{export_cytoscape = TRUE}).
#' @param metabolic_signaling_color Hex colour for KO nodes classified as
#'   \code{"metabolic&signaling"} (used only when
#'   \code{export_cytoscape = TRUE}).
#' @param other_color Hex colour for KO nodes whose class cannot be determined
#'   and for non-compound, non-KO nodes (used only when
#'   \code{export_cytoscape = TRUE}).
#' @param network_title Character. Cytoscape network title (used only when
#'   \code{export_cytoscape = TRUE}).
#' @param collection_title Character. Cytoscape collection title (used only when
#'   \code{export_cytoscape = TRUE}).
#'
#' @details
#' Reverse duplicated edges are collapsed for Cytoscape display when they
#' share the same metabolic interaction type
#' KO functional class is derived from the KEGG BRITE pathway hierarchy
#' (\code{br:br08901}). Pathways under the \code{"Metabolism"} branch are
#' classified as metabolic, whereas all other linked pathways are treated as
#' signaling for this annotation workflow.
#'
#' @return
#' A data frame with columns \code{id}, \code{type}, \code{name}, and
#' \code{ko_class}. Compound nodes are excluded from the returned table.
#'
#' If \code{export_cytoscape = TRUE}, the same annotation is also applied to a
#' Cytoscape network as a side effect, with node colours mapped to
#' \code{ko_class}; the function invisibly returns the annotation data frame.
#'
#' @examples
#' ## 1) Prepare a mixed human–E. coli tryptophan-focused network
#' ## Human: tryptophan metabolism + selected immune/inflammatory signaling pathways
#' ## Load the precomputed example network; see MPN_keggNetwork() for its construction
#' data(MetaPathNet_example_network)
#'
#' ## Combined host–microbe network
#' net_trp_mixed <- MetaPathNet_example_network
#'
#' ## 2) Build a shortest-path subnetwork (used here as the input network_table)
#' source_nodes <- c("K00486", "K00453", "K00463", "K01667", "K01696", "K01695")
#' target_nodes <- c("cpd:C00078", "cpd:C00328", "cpd:C00780", "cpd:C00463", "cpd:C00108")
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
#' ## 3) KO class annotation table (data frame output)
#' ko_class_df <- MPN_annotateKoClass(
#'   network_table = paths_trp_mixed,
#'   name          = FALSE
#' )
#'
#' @seealso
#' \code{\link{MPN_annotateOrigin}} for host vs bacteria origin, and
#' \code{\link{MPN_viewNetworkCy}} for generic Cytoscape visualisation of
#' MetaPathNet networks.
#'
#' @importFrom KEGGREST keggGet keggList keggLink
#' @importFrom RCy3 cytoscapePing createNetworkFromDataFrames
#'   setNodeColorMapping setNodeLabelMapping
#' @export
MPN_annotateKoClass <- function(
    network_table,
    name  = FALSE,
    export_cytoscape = FALSE,
    compound_color             = "#9DC7DD",
    metabolic_color            = "#9ED17B",
    signaling_color            = "#C45745",
    metabolic_signaling_color  = "#B57F60",
    other_color                = "#D1C2C2",
    network_title    = "",
    collection_title = ""
) {
  #==============================================================
  # STEP 1 - Validate input network and detect node types
  #==============================================================

  # Input must be an edge table with at least two columns (source, target)
  if (missing(network_table)) {
    stop("Argument 'network_table' (MetaPathNet-style network) is required.")
  }

  if (!is.logical(name) || length(name) != 1 || is.na(name)) {
    stop("Argument 'name' must be TRUE or FALSE.")
  }

  # Coerce to data.frame and standardise column names
  net_df <- as.data.frame(network_table, stringsAsFactors = FALSE)
  if (ncol(net_df) < 2) {
    stop("Input 'network_table' must have at least 2 columns (source, target).")
  }

  colnames(net_df)[seq_len(2)] <- c("source", "target")

  # Keep interaction_type if present, otherwise create it
  if (ncol(net_df) >= 3) {
    colnames(net_df)[3] <- "interaction_type"
  } else {
    net_df$interaction_type <- NA_character_
  }

  # Build the unique node set from edge endpoints
  node_ids <- unique(c(net_df$source, net_df$target))

  # Classify node "type" using KEGG ID patterns
  node_type <- character(length(node_ids))
  node_type[grepl("^(cpd:|gl:|dr:)", node_ids)] <- "Compound"
  node_type[grepl("^K\\d{5}$", node_ids)]       <- "KO"
  node_type[node_type == ""]                    <- "Other"

  #==============================================================
  # STEP 2 - Confirm Cytoscape is reachable (RCy3) if export is requested
  #==============================================================
  # cytoscapePing() must succeed before creating the network
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
  # STEP 3 - Optional KEGG label conversion (KO: SYMBOL => NAME; others: NAME)
  #==============================================================
  # Build a label map only when KEGG lookups are requested
  name_map <- NULL
  if (isTRUE(name)) {
    message("Converting KEGG IDs to labels")

    # Preallocate a named character vector for later lookup
    name_map <- setNames(rep(NA_character_, length(node_ids)), node_ids)

    # Query KEGG once per node and store the best display label
    for (i in seq_along(node_ids)) {
      id_i <- node_ids[i]

      ## Pause every 10 queries to reduce KEGG API pressure
      if (i %% 10 == 0) {
        Sys.sleep(1)
      }

      # Retrieve KEGG entry; skip cleanly on any API error
      info_i <- tryCatch(
        KEGGREST::keggGet(id_i),
        error = function(e) NULL
      )
      if (is.null(info_i) || length(info_i) == 0) next
      entry <- info_i[[1]]

      # KO nodes: prefer SYMBOL (first token), fall back to NAME
      if (grepl("^K\\d{5}$", id_i)) {
        sym <- entry$SYMBOL

        # SYMBOL can contain multiple symbols in one string; keep the first one
        if (!is.null(sym) && length(sym) > 0) {
          sym1 <- trimws(strsplit(sym[1], ",", fixed = TRUE)[[1]][1])
          if (nchar(sym1) > 0) {
            name_map[id_i] <- sym1
            next
          }
        }

        # NAME fallback for KO (strip bracket annotations)
        nm <- entry$NAME
        if (!is.null(nm) && length(nm) > 0) {
          nm1 <- trimws(gsub("\\[.*?\\]", "", nm[1]))
          if (nchar(nm1) > 0) {
            name_map[id_i] <- nm1
          }
        }

      } else {
        # Non-KO nodes: use NAME (clean punctuation for compounds)
        nm <- entry$NAME
        if (is.null(nm) || length(nm) == 0) next

        nm1 <- nm[1]
        if (grepl("^cpd:", id_i)) nm1 <- gsub(";", "", nm1)
        nm1 <- trimws(nm1)

        if (nchar(nm1) > 0) {
          name_map[id_i] <- nm1
        }
      }
    }
  }

  #==============================================================
  # STEP 4 - Build metabolic vs signalling pathway ID sets (BRITE)
  #==============================================================
  # Retrieve the BRITE pathway map hierarchy and split it into lines
  br <- tryCatch(
    KEGGREST::keggGet("br:br08901")[[1]],
    error = function(e) stop("Failed to retrieve KEGG BRITE hierarchy (br08901).")
  )

  # Some KEGG REST responses return a single string; enforce one line per element
  if (length(br) == 1L) {
    lines <- strsplit(br, "\n", fixed = TRUE)[[1]]
  } else {
    lines <- br
  }

  # Normalise line endings and leading spaces for stable parsing
  lines <- sub("\r$", "", lines)
  trim  <- sub("^\\s+", "", lines)

  # Track the closest preceding A-level category for each line
  current_A <- character(length(trim))
  lastA <- NA_character_
  for (i in seq_along(trim)) {
    if (grepl("^A", trim[i])) lastA <- trim[i]
    current_A[i] <- lastA
  }

  # Keep only pathway ID lines and extract their 5-digit identifiers
  isC     <- grepl("^C\\s+[0-9]{5}\\b", trim)
  C_lines <- trim[isC]
  A_for_C <- current_A[isC]

  ids   <- sub("^C\\s+([0-9]{5}).*$", "\\1", C_lines)
  A_cat <- sub("^A", "", A_for_C)

  # Classification rule: Metabolism = metabolic, everything else = signalling
  metabolic_paths <- unique(ids[A_cat == "Metabolism"])
  signaling_paths <- unique(ids[A_cat != "Metabolism"])

  #==============================================================
  # STEP 5 - Map KO => pathways (keggLink) and assign KO class
  #==============================================================
  # Only KO nodes are classified; other nodes remain as "other"
  kos <- unique(node_ids[node_type == "KO"])

  # Store the pathway IDs (comma-separated) plus the final KO class
  ko_path_df <- data.frame(
    KO         = kos,
    Pathway_ID = NA_character_,
    KO_class   = NA_character_,
    stringsAsFactors = FALSE
  )

  # Query KEGG once per KO and classify based on pathway set membership
  if (length(kos) > 0) {
    for (i in seq_along(kos)) {
      ko_id <- kos[i]

      ## Pause every 10 queries to reduce KEGG API pressure
      if (i %% 10 == 0) {
        Sys.sleep(1)
      }

      # Link KO to pathways; skip cleanly on any API error
      link_res <- tryCatch(
        KEGGREST::keggLink("pathway", ko_id),
        error = function(e) NULL
      )
      if (is.null(link_res) || length(link_res) == 0L) next

      # Extract linked pathway IDs and remove generic KO-pathway entries
      ids_full <- sub("^path:", "", as.vector(link_res))
      ids_full <- ids_full[!grepl("^ko", ids_full)]

      # Reduce to 5-digit identifiers and keep unique non-empty values
      ids_num <- unique(sub("^[a-zA-Z]+", "", ids_full))
      ids_num <- ids_num[ids_num != ""]
      if (length(ids_num) == 0L) next

      ko_path_df$Pathway_ID[i] <- paste(ids_num, collapse = ",")

      # Classify KO by whether it maps to metabolic, signalling, or both
      in_met <- ids_num %in% metabolic_paths
      in_sig <- ids_num %in% signaling_paths

      if (any(in_met) && !any(in_sig)) {
        ko_path_df$KO_class[i] <- "metabolic"
      } else if (!any(in_met) && any(in_sig)) {
        ko_path_df$KO_class[i] <- "signaling"
      } else if (any(in_met) && any(in_sig)) {
        ko_path_df$KO_class[i] <- "metabolic&signaling"
      } else {
        # Leave NA when none of the pathways match the BRITE-derived sets
        message("KO ", ko_id, ": pathways not recognised in BRITE sets - leaving NA.")
      }
    }
  }

  # Build a lookup from KO ID to KO_class
  ko_class_map <- setNames(ko_path_df$KO_class, ko_path_df$KO)

  #==============================================================
  # STEP 6 - Build annotation table, optionally export to Cytoscape
  #==============================================================
  # Node table must include an "id" column; KO class is stored as a node attribute
  nodes_df <- data.frame(
    id    = node_ids,
    type  = node_type,
    stringsAsFactors = FALSE
  )

  # Default: everything is "other" unless assigned below
  nodes_df$ko_class <- "other"

  # Compounds use a separate visual category
  nodes_df$ko_class[nodes_df$type == "Compound"] <- "compound"

  # KO nodes receive the computed class; fall back to "other" when NA
  is_ko <- nodes_df$type == "KO"
  if (any(is_ko)) {
    ko_ids <- nodes_df$id[is_ko]
    mapped <- ko_class_map[ko_ids]
    mapped[is.na(mapped)] <- "other"
    nodes_df$ko_class[is_ko] <- unname(mapped)
  }

  # Label column controls the Cytoscape node text display and dataframe names
  nodes_df$label <- nodes_df$id
  if (isTRUE(name) && !is.null(name_map)) {
    repl <- name_map[nodes_df$id]
    repl[is.na(repl)] <- nodes_df$id[is.na(repl)]
    nodes_df$label <- repl
  }

  # Dataframe output is the primary return (compound nodes excluded)
  out_df <- nodes_df[, c("id", "type", "label", "ko_class"), drop = FALSE]
  out_df <- out_df[out_df$type != "Compound", , drop = FALSE]
  colnames(out_df) <- c("id", "type", "name", "ko_class")

  if (!isTRUE(export_cytoscape)) {
    return(out_df)
  }

  # Edge table uses standard source/target plus optional interaction_type
  edges_df <- data.frame(
    source           = net_df$source,
    target           = net_df$target,
    interaction_type = net_df$interaction_type,
    stringsAsFactors = FALSE
  )

  # Collapse duplicated reverse edges for Cytoscape display:
  # if both A->B and B->A exist under the same metabolic interaction type,
  # keep only one edge. Different interaction types are kept separate.
  edges_vis <- edges_df[, c("source", "target", "interaction_type")]

  collapse_types <- c(
    "k_compound:reversible",
    "custom:reversible",
    "k_compound:irreversible",
    "custom:irreversible"
  )
  collapse_idx <- which(edges_vis$interaction_type %in% collapse_types)

  if (length(collapse_idx) > 0L) {
    collapse_edges <- edges_vis[collapse_idx, , drop = FALSE]
    other_edges    <- edges_vis[-collapse_idx, , drop = FALSE]
    collapse_key <- paste(
      pmin(collapse_edges$source, collapse_edges$target),
      pmax(collapse_edges$source, collapse_edges$target),
      collapse_edges$interaction_type,
      sep = "__"
    )
    collapse_edges <- collapse_edges[!duplicated(collapse_key), , drop = FALSE]
    edges_vis <- rbind(other_edges, collapse_edges)
  }

  # Create the Cytoscape network in the requested collection
  RCy3::createNetworkFromDataFrames(
    nodes      = nodes_df,
    edges      = edges_vis,
    title      = network_title,
    collection = collection_title
  )

  # Apply colour mapping to the KO class attribute
  class_colors <- c(
    compound              = compound_color,
    metabolic             = metabolic_color,
    signaling             = signaling_color,
    "metabolic&signaling" = metabolic_signaling_color,
    other                 = other_color
  )

  RCy3::setNodeColorMapping(
    table.column        = "ko_class",
    table.column.values = names(class_colors),
    colors              = unname(class_colors),
    mapping.type        = "d"
  )

  # Apply node label mapping last
  RCy3::setNodeLabelMapping(table.column = "label")

  invisible(out_df)
}
