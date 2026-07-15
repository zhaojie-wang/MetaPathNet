#' Annotate KO Origin (Host vs Bacteria) and Optionally Export to Cytoscape
#'
#' Annotates nodes in a MetaPathNet-style network with inferred biological
#' origin, focusing on KEGG Orthology (KO) nodes. KO origin is classified as
#' human-only (\code{"hsa"}), bacteria-only (\code{"bacteria"}),
#' shared (\code{"hsa&bacteria"}), or \code{"other"} based on KEGG annotations.
#' Compound nodes are assigned the separate class \code{"compound"}.
#'
#' The primary output is a node-level annotation table containing the node ID,
#' type, display name, and origin. Cytoscape export is optional and can be
#' enabled to visualise the annotated network with origin-based node colours.
#'
#' KO origin is inferred from KEGG \code{GENES} annotations retrieved with
#' \code{KEGGREST::keggGet()}:
#' \itemize{
#'   \item \code{"hsa"} if \code{hsa} is present and no selected bacterial
#'     organism code is detected.
#'   \item \code{"bacteria"} if at least one selected bacterial organism code
#'     is detected and \code{hsa} is absent.
#'   \item \code{"hsa&bacteria"} if both \code{hsa} and at least one selected
#'     bacterial organism code are detected.
#'   \item \code{"other"} if none of the above can be established.
#' }
#'
#' The \code{bacteria_codes} argument can be used to restrict the bacterial
#' origin class to a user-defined organism panel, such as cohort-specific
#' microbes. If \code{bacteria_codes = NULL}, all KEGG bacterial organism
#' codes are used.
#'
#' @param network_table A MetaPathNet-style network with at least two columns:
#'   \code{source} and \code{target}. If a third column is present, it is
#'   interpreted as \code{interaction_type}.
#' @param bacteria_codes Optional character vector of KEGG bacterial organism
#'   codes. Restricts the \code{"bacteria"} and \code{"hsa&bacteria"} classes
#'   to this set. If \code{NULL}, all KEGG bacterial organism codes are used.
#' @param name Logical. If \code{TRUE}, node display names are converted from
#'   KEGG entries. For KO nodes, \code{SYMBOL} is preferred and \code{NAME}
#'   is used as a fallback. For non-KO nodes, the KEGG \code{NAME} field is
#'   used. If \code{FALSE}, KEGG IDs are retained as labels.
#' @param export_cytoscape Logical. If \code{TRUE}, the annotated network is
#'   exported to Cytoscape and styled by origin class. If \code{FALSE}, only
#'   the annotation table is returned.
#' @param compound_color Hex colour for compound nodes.
#' @param hsa_color Hex colour for KO nodes classified as \code{"hsa"}.
#' @param micro_color Hex colour for KO nodes classified as
#'   \code{"bacteria"}.
#' @param hsa_micro_color Hex colour for KO nodes classified as
#'   \code{"hsa&bacteria"}.
#' @param other_color Hex colour for KO nodes classified as \code{"other"}.
#' @param network_title Character string specifying the Cytoscape network
#'   title. Used only when \code{export_cytoscape = TRUE}.
#' @param collection_title Character string specifying the Cytoscape collection
#'   title. Used only when \code{export_cytoscape = TRUE}.
#'
#' @details
#' If the KEGG genome or bacterial genome-group information cannot be
#' retrieved, the function issues a warning instead of stopping. In this case,
#' KO nodes linked exclusively to \code{hsa} are classified as \code{"hsa"},
#' while other KO nodes are classified as \code{"other"}.
#'
#' Cytoscape version 3.9.0 or later must be running when
#' \code{export_cytoscape = TRUE}.
#'
#' Reverse duplicated edges are collapsed for Cytoscape display when they
#' share the same supported metabolic interaction type.
#'
#' @return A node-level data frame with columns:
#' \itemize{
#'   \item \code{id}: KEGG node identifier.
#'   \item \code{type}: coarse node type (\code{"KO"}, \code{"Compound"},
#'     or \code{"Other"}).
#'   \item \code{name}: display label containing the KEGG ID or converted
#'     KEGG name or symbol.
#'   \item \code{origin}: inferred origin class (\code{"compound"},
#'     \code{"hsa"}, \code{"bacteria"}, \code{"hsa&bacteria"}, or
#'     \code{"other"}).
#' }
#'
#' If \code{export_cytoscape = TRUE}, the function also creates and styles a
#' Cytoscape network as a side effect and returns the annotation table
#' invisibly.
#'
#' @examples
#' \donttest{
#' ## 1) Load the precomputed mixed human-E. coli example network
#' data(MetaPathNet_example_network)
#' net_trp_mixed <- MetaPathNet_example_network
#'
#' ## 2) Build a shortest-path subnetwork for origin annotation
#' source_nodes <- c(
#'   "K00486", "K00453", "K00463", "K01667",
#'   "K01696", "K01610", "K01695"
#' )
#'
#' target_nodes <- c(
#'   "cpd:C00078", "cpd:C00328", "cpd:C00780", "cpd:C00463",
#'   "cpd:C00108", "cpd:C00458", "cpd:C00336"
#' )
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
#' ## 3) Return the annotation table while retaining KEGG IDs as labels
#' origin_paths_trp_mixed <- MPN_annotateOrigin(
#'   network_table    = paths_trp_mixed,
#'   bacteria_codes   = NULL,
#'   name             = FALSE,
#'   export_cytoscape = FALSE
#' )
#'
#' head(origin_paths_trp_mixed)
#' }
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
    name = FALSE,
    export_cytoscape = FALSE,
    compound_color = "#9DC7DD",
    hsa_color = "#C45745",
    micro_color = "#9ED17B",
    hsa_micro_color = "#B57F60",
    other_color = "#D1C2C2",
    network_title = "",
    collection_title = ""
) {
  #==============================================================
  # STEP 1 - Basic input checks and matrix normalisation
  #==============================================================

  # Ensure the input network matrix is provided
  if (missing(network_table)) {
    stop(
      "Argument 'network_table' ",
      "(MetaPathNet-style 3-column network) is required."
    )
  }

  # Validate 'export_cytoscape' argument
  if (
    !is.logical(export_cytoscape) ||
    length(export_cytoscape) != 1L ||
    is.na(export_cytoscape)
  ) {
    stop("Argument 'export_cytoscape' must be TRUE or FALSE.")
  }

  # Validate 'name' argument
  if (!is.logical(name) || length(name) != 1L || is.na(name)) {
    stop("Argument 'name' must be TRUE or FALSE.")
  }

  # Validate and normalise optional bacterial organism codes
  if (!is.null(bacteria_codes)) {
    if (!is.character(bacteria_codes)) {
      stop(
        "Argument 'bacteria_codes' must be NULL or a character vector ",
        "of KEGG organism codes."
      )
    }

    bacteria_codes <- unique(tolower(trimws(bacteria_codes)))
    bacteria_codes <- bacteria_codes[
      !is.na(bacteria_codes) & nzchar(bacteria_codes)
    ]

    if (length(bacteria_codes) == 0L) {
      stop(
        "Argument 'bacteria_codes' contains no usable organism codes ",
        "after trimming."
      )
    }
  }

  # Convert matrix to data frame and check structure
  net_df <- as.data.frame(network_table, stringsAsFactors = FALSE)

  if (ncol(net_df) < 2L) {
    stop(
      "Input 'network_table' must have at least 2 columns ",
      "(source, target)."
    )
  }

  # Standardise column names and add missing 'interaction_type' if needed
  colnames(net_df)[seq_len(2L)] <- c("source", "target")

  if (ncol(net_df) >= 3L) {
    colnames(net_df)[3L] <- "interaction_type"
  } else {
    net_df$interaction_type <- NA_character_
  }

  # Normalise and validate source/target node IDs
  net_df$source <- trimws(as.character(net_df$source))
  net_df$target <- trimws(as.character(net_df$target))
  net_df$interaction_type <- as.character(net_df$interaction_type)

  invalid_edges <- is.na(net_df$source) |
    !nzchar(net_df$source) |
    is.na(net_df$target) |
    !nzchar(net_df$target)

  if (any(invalid_edges)) {
    stop(
      "Columns 'source' and 'target' must contain non-missing, ",
      "non-empty node identifiers."
    )
  }

  #==============================================================
  # STEP 2 - Classify node types (Compound, KO, Other)
  #==============================================================

  # Extract all unique node IDs from the network
  node_ids <- unique(c(net_df$source, net_df$target))

  # Assign coarse node types based on identifier patterns
  node_type <- character(length(node_ids))
  node_type[grepl("^(cpd:|gl:|dr:)", node_ids)] <- "Compound"
  node_type[grepl("^K\\d{5}$", node_ids)] <- "KO"
  node_type[node_type == ""] <- "Other"

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

    if (length(non_ko_ids) > 0L) {
      name_map_nonko <- setNames(
        rep(NA_character_, length(non_ko_ids)),
        non_ko_ids
      )

      for (i in seq_along(non_ko_ids)) {
        id_i <- non_ko_ids[i]

        ## Pause every 10 queries to reduce KEGG API pressure
        if (i %% 10L == 0L) {
          Sys.sleep(1)
        }

        # Query KEGG REST for compound/drug ID
        info_i <- tryCatch(
          KEGGREST::keggGet(id_i),
          error = function(e) NULL
        )

        if (
          is.null(info_i) ||
          length(info_i) == 0L ||
          is.null(info_i[[1]]$NAME) ||
          length(info_i[[1]]$NAME) == 0L
        ) {
          next
        }

        # Use first name field; strip trailing semicolon if present
        nm <- as.character(info_i[[1]]$NAME[1L])
        nm <- trimws(sub(";.*$", "", nm))

        if (nzchar(nm)) {
          name_map_nonko[id_i] <- nm
        }
      }
    }
  }

  #==============================================================
  # STEP 5 - Parse KO nodes: labels, associated organisms, origin class
  #==============================================================

  kos <- unique(node_ids[node_type == "KO"])

  # Prepare a data frame to store KO-level metadata
  ko_df <- data.frame(
    KO = kos,
    organism_codes = I(vector("list", length(kos))),
    clean_code = NA_character_,
    stringsAsFactors = FALSE
  )

  # Prepare optional KO name mapping (SYMBOL or NAME)
  ko_name_map <- if (isTRUE(name) && length(kos) > 0L) {
    setNames(rep(NA_character_, length(kos)), kos)
  } else {
    NULL
  }

  if (length(kos) > 0L) {
    for (i in seq_along(kos)) {
      ko_id <- kos[i]

      ## Pause every 10 queries to reduce KEGG API pressure
      if (i %% 10L == 0L) {
        Sys.sleep(1)
      }

      # Query KEGG for KO entry
      res <- tryCatch(
        KEGGREST::keggGet(ko_id),
        error = function(e) NULL
      )

      if (is.null(res) || length(res) == 0L) {
        next
      }

      entry <- res[[1]]

      # 5.1 - KO label (SYMBOL preferred, NAME fallback)
      if (isTRUE(name)) {
        used_label <- NA_character_

        if (!is.null(entry$SYMBOL) && length(entry$SYMBOL) > 0L) {
          used_label <- trimws(
            strsplit(
              entry$SYMBOL[1L],
              ",",
              fixed = TRUE
            )[[1L]][1L]
          )
        }

        if (is.na(used_label) || nchar(used_label) == 0L) {
          if (!is.null(entry$NAME) && length(entry$NAME) > 0L) {
            used_label <- trimws(
              gsub("\\[.*?\\]", "", entry$NAME[1L])
            )
          }
        }

        ko_name_map[ko_id] <- used_label
      }

      # 5.2 - Extract organism codes from KO -> GENES
      genes_field <- entry$GENES

      if (is.null(genes_field) || length(genes_field) == 0L) {
        next
      }

      genes_vec <- if (is.list(genes_field)) {
        unlist(genes_field, use.names = TRUE)
      } else {
        genes_field
      }

      ## Prefer top-level GENES names when KEGG returns a list
      if (is.list(genes_field) && !is.null(names(genes_field))) {
        org_codes <- names(genes_field)
      } else {
        org_codes <- names(genes_vec)
      }

      if (
        is.null(org_codes) ||
        !any(!is.na(org_codes) & nzchar(org_codes))
      ) {
        org_codes <- sub(
          "^([^:]+):.*$",
          "\\1",
          as.character(genes_vec)
        )
      }

      # Clean organism codes extracted from KEGG GENES annotations
      org_codes <- tolower(trimws(org_codes))
      org_codes <- sub(":.*$", "", org_codes)
      org_codes <- sub("\\..*$", "", org_codes)

      org_codes <- unique(
        org_codes[
          !is.na(org_codes) &
            nzchar(org_codes) &
            grepl("^[a-z][a-z0-9_]*$", org_codes)
        ]
      )

      ko_df$organism_codes[[i]] <- org_codes
      ko_df$clean_code[i] <- paste(org_codes, collapse = ",")
    }
  }

  #==============================================================
  # STEP 6 - Classify KO origin using KEGG genome groups
  #==============================================================

  # Retrieve the complete KEGG genome table and bacterial genome group
  genome <- tryCatch(
    KEGGREST::keggList("genome"),
    error = function(e) NULL
  )

  bacteria <- tryCatch(
    KEGGREST::keggList("genome", "bacteria"),
    error = function(e) NULL
  )

  taxonomy_available <- TRUE
  organism_bacteria <- character(0)

  if (
    is.null(genome) ||
    length(genome) == 0L ||
    is.null(names(genome)) ||
    is.null(bacteria) ||
    length(bacteria) == 0L
  ) {
    taxonomy_available <- FALSE

    warning(
      "KEGG bacterial genome classification could not be retrieved. ",
      "KO nodes linked exclusively to hsa will be classified as \"hsa\"; ",
      "other KO nodes will be classified as \"other\"."
    )
  } else {
    genome_value <- trimws(as.character(genome))
    bacteria_value <- trimws(as.character(bacteria))

    valid_genome <- !is.na(names(genome)) &
      nzchar(names(genome)) &
      !is.na(genome_value) &
      nzchar(genome_value) &
      grepl(";", genome_value, fixed = TRUE)

    valid_bacteria <- !is.na(bacteria_value) &
      nzchar(bacteria_value) &
      grepl(";", bacteria_value, fixed = TRUE)

    if (!any(valid_genome) || !any(valid_bacteria)) {
      taxonomy_available <- FALSE

      warning(
        "KEGG genome groups returned no parseable bacterial entries. ",
        "KO nodes linked exclusively to hsa will be classified as \"hsa\"; ",
        "other KO nodes will be classified as \"other\"."
      )
    } else {
      genome <- genome[valid_genome]
      genome_value <- genome_value[valid_genome]

      # Reconstruct the four-column KEGG organism table
      organism_df <- data.frame(
        taxon_id = names(genome),
        organism = trimws(sub(";.*$", "", genome_value)),
        species = trimws(sub("^[^;]+;\\s*", "", genome_value)),
        phylogeny = NA_character_,
        stringsAsFactors = FALSE
      )

      organism_df <- organism_df[
        nzchar(organism_df$taxon_id) &
          nzchar(organism_df$organism) &
          nzchar(organism_df$species),
        ,
        drop = FALSE
      ]

      bacteria_kegg_codes <- unique(tolower(trimws(
        sub(
          ";.*$",
          "",
          bacteria_value[valid_bacteria]
        )
      )))

      organism_df$phylogeny[
        tolower(organism_df$organism) %in% bacteria_kegg_codes
      ] <- "Prokaryotes;Bacteria"

      organism_bacteria <- tolower(
        organism_df$organism[
          !is.na(organism_df$phylogeny) &
            organism_df$phylogeny == "Prokaryotes;Bacteria"
        ]
      )

      organism_bacteria <- unique(organism_bacteria)

      if (length(organism_bacteria) == 0L) {
        taxonomy_available <- FALSE

        warning(
          "No bacterial organism codes could be matched to the KEGG ",
          "genome table. KO nodes linked exclusively to hsa will be ",
          "classified as \"hsa\"; other KO nodes will be classified ",
          "as \"other\"."
        )
      }
    }
  }

  # Optional restriction to user-provided bacterial codes
  if (!is.null(bacteria_codes) && taxonomy_available) {
    unrecognised_codes <- setdiff(
      bacteria_codes,
      organism_bacteria
    )

    if (length(unrecognised_codes) > 0L) {
      message(
        "The following bacteria_codes were not recognised and were excluded:\n  ",
        paste(unrecognised_codes, collapse = ", ")
      )
    }

    organism_bacteria <- intersect(
      bacteria_codes,
      organism_bacteria
    )

    if (length(organism_bacteria) == 0L) {
      warning(
        "None of the supplied bacteria_codes matched KEGG bacterial ",
        "organism codes."
      )
    }
  }

  # Assign origin class for each KO: hsa, bacteria, hsa&bacteria, other
  ko_df$origin <- vapply(
    ko_df$clean_code,
    function(cc) {
      if (is.na(cc) || nchar(cc) == 0L) {
        return("other")
      }

      codes <- strsplit(cc, ",", fixed = TRUE)[[1L]]

      if (!taxonomy_available) {
        if (length(codes) == 1L && codes == "hsa") {
          return("hsa")
        }

        return("other")
      }

      has_hsa <- "hsa" %in% codes
      has_micro <- any(codes %in% organism_bacteria)

      if (has_hsa && has_micro) {
        "hsa&bacteria"
      } else if (has_hsa) {
        "hsa"
      } else if (has_micro) {
        "bacteria"
      } else {
        "other"
      }
    },
    FUN.VALUE = character(1)
  )

  ko_origin_map <- setNames(ko_df$origin, ko_df$KO)

  #==============================================================
  # STEP 7 - Build node table and assign final labels
  #==============================================================

  nodes_df <- data.frame(
    id = node_ids,
    type = node_type,
    origin = "other",
    stringsAsFactors = FALSE
  )

  # Assign origin class for compound and KO nodes
  nodes_df$origin[nodes_df$type == "Compound"] <- "compound"
  is_ko <- nodes_df$type == "KO"

  if (any(is_ko)) {
    nodes_df$origin[is_ko] <- ko_origin_map[
      nodes_df$id[is_ko]
    ]
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

  out_df <- nodes_df[
    ,
    c("id", "type", "label", "origin"),
    drop = FALSE
  ]

  colnames(out_df)[colnames(out_df) == "label"] <- "name"

  if (!isTRUE(export_cytoscape)) {
    return(out_df)
  }

  # Build edge table from input matrix
  edges_df <- data.frame(
    source = net_df$source,
    target = net_df$target,
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

  collapse_idx <- which(
    edges_df$interaction_type %in% collapse_types
  )

  if (length(collapse_idx) > 0L) {
    collapse_edges <- edges_df[
      collapse_idx,
      ,
      drop = FALSE
    ]

    other_edges <- edges_df[
      -collapse_idx,
      ,
      drop = FALSE
    ]

    collapse_key <- paste(
      pmin(collapse_edges$source, collapse_edges$target),
      pmax(collapse_edges$source, collapse_edges$target),
      collapse_edges$interaction_type,
      sep = "__"
    )

    collapse_edges <- collapse_edges[
      !duplicated(collapse_key),
      ,
      drop = FALSE
    ]

    edges_df <- rbind(other_edges, collapse_edges)
  }

  # Create Cytoscape network via RCy3
  RCy3::createNetworkFromDataFrames(
    nodes = nodes_df,
    edges = edges_df,
    title = network_title,
    collection = collection_title
  )

  # Define color mapping by origin class
  origin_colors <- c(
    compound = compound_color,
    hsa = hsa_color,
    bacteria = micro_color,
    "hsa&bacteria" = hsa_micro_color,
    other = other_color
  )

  # Apply visual styling to Cytoscape network
  RCy3::setNodeColorMapping(
    table.column = "origin",
    table.column.values = names(origin_colors),
    colors = unname(origin_colors),
    mapping.type = "d"
  )

  RCy3::setNodeLabelMapping(
    table.column = "label"
  )

  invisible(out_df)
}
