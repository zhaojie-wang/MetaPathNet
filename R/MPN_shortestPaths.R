#' Identify Shortest Paths Between Node Sets in a KEGG-Derived Network
#'
#' Computes shortest paths between specified source and target node sets
#' in a KEGG-based network. Source nodes are often KEGG KOs/genes
#' (upstream biological functions) and target nodes are often KEGG compounds
#' (downstream metabolites), but this is a recommended usage pattern rather
#' than a strict requirement.
#'
#' Valid source-target pairs are first filtered by graph distance, and paths
#' can be returned either as path lists or as a condensed edge list.
#'
#' @param network_table A 3-column matrix or data.frame of edges with columns
#'   \code{source}, \code{target}, \code{interaction_type}, typically from
#'   \code{MPN_keggNetwork()} or \code{MPN_crossSpeciesNetwork()}.
#' @param source_nodes Character vector of KEGG node IDs to use as path sources.
#' @param target_nodes Character vector of KEGG node IDs to use as path targets.
#' @param mode Character; direction mode passed to \code{igraph} shortest-path
#'   and distance functions. One of:
#'   \describe{
#'     \item{\code{"out"}}{Shortest-path search follows edge direction.}
#'     \item{\code{"in"}}{Shortest-path search is computed against edge direction (reverse direction).}
#'     \item{\code{"all"}}{Ignore edge direction during shortest-path search
#'     (graph treated as undirected). For \code{output = "network_matrix"},
#'     interaction types are still recovered from the original network edge
#'     annotation.}
#'   }
#' @param output Character; output format:
#'   \itemize{
#'     \item \code{"path_list"} - list of shortest paths per source-target pair.
#'     \item \code{"network_matrix"} - 3-column edge list (source, target, interaction_type).
#'   }
#' @param name Logical; if \code{TRUE} and \code{output = "path_list"}, node
#'   IDs are converted to KEGG names or symbols via \pkg{KEGGREST}.
#' @param distance_threshold Numeric; only source-target pairs with
#'   shortest-path distance \eqn{\le} this value are processed (default \code{Inf}).
#' @param betweenness Logical (default \code{FALSE}). If \code{TRUE}, when more than
#'   one shortest path of the same length exists for a given source-target pair,
#'   a single path is selected by ranking candidate shortest paths by node betweenness.
#' @param reactant_filter Logical (default \code{FALSE}). If \code{TRUE},
#'   shortest paths are filtered to keep only paths whose consecutive
#'   compound-compound transitions are supported by KEGG reaction data.
#'   This step issues repeated KEGG API queries and may be slow for large
#'   path sets.
#'
#' @details
#' The function:
#' \itemize{
#'   \item Accepts \code{matrix} or \code{data.frame} input and internally
#'         standardizes it to a character edge table.
#'   \item Builds a directed \pkg{igraph} graph from \code{network_table} and
#'         computes pairwise shortest-path distances (breadth-first search on
#'         this unweighted graph).
#'   \item Stops if any node is present in both \code{source_nodes} and
#'         \code{target_nodes}, because the distance from a node to itself is 0
#'         and is not allowed in this workflow.
#'   \item Keeps only mapped source-target pairs with finite distance within
#'         \code{distance_threshold}; unreachable pairs (\code{Inf}) and pairs
#'         longer than the threshold are filtered out and not processed further.
#'   \item If \code{mode = "all"}, shortest-path search ignores edge direction.
#'         However, for \code{output = "network_matrix"}, the returned
#'         \code{interaction_type} values are still recovered from the original
#'         directed edge annotation in \code{network_table}, including reverse
#'         matching when needed.
#'   \item If \code{reactant_filter = TRUE}, candidate shortest paths are
#'         further filtered by checking whether consecutive compound-compound
#'         transitions are supported by KEGG reaction data.
#' }
#'
#' @return
#' Depending on \code{output}:
#' \itemize{
#'   \item \code{"path_list"} - named list; each element is one or more shortest
#'         paths (character vectors of node IDs or names).
#'   \item \code{"network_matrix"} - 3-column matrix of unique edges traversed
#'         by the identified shortest paths. Edges are returned in path order,
#'         and \code{interaction_type} is assigned from the original network
#'         edge annotation.
#' }
#'
#' @note
#' Name/symbol conversion relies on the KEGG REST API and can be slow for large
#' node sets. The internal 6-step cutoff is empirical and intended to prevent
#' excessive enumeration of very long shortest paths. When \code{betweenness = TRUE},
#' enumerating all shortest paths for tie-breaking may increase runtime and memory use.
#'
#' @examples
#' ## Load the precomputed example network; see MPN_keggNetwork() for its construction
#' data(MetaPathNet_example_network)
#' ## Use the precomputed object as the example network_table
#' net_trp_mixed <- MetaPathNet_example_network
#'
#' ## 2) Define source and target nodes of interest in this network
#' ## Source/target can be any KEGG node types, but they must be mapped in the network
#' ## and should not overlap (same node cannot be both source and target).
#' source_nodes <- c("K00486", "K00453", "K00463",
#'                   "K01667", "K01696", "K01610", "K01695")
#' target_nodes <- c("cpd:C00078", "cpd:C00328", "cpd:C00780",
#'                   "cpd:C00463", "cpd:C00108", "cpd:C00458", "cpd:C00336")
#'
#' map_src <- MPN_findMappedNodes(source_nodes, net_trp_mixed)
#' map_tgt <- MPN_findMappedNodes(target_nodes, net_trp_mixed)
#'
#' ## 3) Optional: inspect shortest-path distances between mapped node sets
#' ## (useful for checking connectivity before extracting shortest paths)
#' dist_trp <- MPN_distances(
#'   network_table = net_trp_mixed,
#'   source_nodes  = map_src$mapped_nodes,
#'   target_nodes  = map_tgt$mapped_nodes,
#'   name          = FALSE
#' )
#'
#' ## 4) Extract shortest paths as annotated path list
#' ## Pairs that are unreachable (Inf distance) or above the distance threshold
#' ## are filtered internally.
#' paths_trp_mixed_list <- MPN_shortestPaths(
#'   network_table      = net_trp_mixed,
#'   source_nodes       = map_src$mapped_nodes,
#'   target_nodes       = map_tgt$mapped_nodes,
#'   mode               = "out",
#'   output             = "path_list",
#'   name               = FALSE,
#'   distance_threshold = 6
#' )
#'
#' ## 5) Extract a subnetwork containing all edges used in the selected shortest paths
#' paths_trp_mixed <- MPN_shortestPaths(
#'   network_table      = net_trp_mixed,
#'   source_nodes       = map_src$mapped_nodes,
#'   target_nodes       = map_tgt$mapped_nodes,
#'   mode               = "out",
#'   output             = "network_matrix",
#'   name               = FALSE
#' )
#'
#' @export
MPN_shortestPaths <- function(
    network_table,
    source_nodes,
    target_nodes,
    mode   = c("out", "in", "all"),
    output = c("path_list", "network_matrix"),
    name   = FALSE,
    distance_threshold = Inf,
    betweenness = FALSE,
    reactant_filter = FALSE
) {
  ## ==============================================================
  ## Step 1 - Input checks
  ## ==============================================================
  mode   <- match.arg(mode)
  output <- match.arg(output)

  # Accept both matrix and data.frame input, then standardise to matrix
  if (is.data.frame(network_table)) {
    network_table <- as.matrix(network_table)
  }

  if (!is.matrix(network_table) || ncol(network_table) < 3) {
    stop("network_table must be a matrix with at least 3 columns: source, target, interaction_type.")
  }

  # Work with character data
  storage.mode(network_table) <- "character"

  network_df <- as.data.frame(network_table, stringsAsFactors = FALSE)
  colnames(network_df)[seq_len(3)] <- c("source", "target", "interaction_type")

  if (!is.character(source_nodes) || !length(source_nodes)) {
    stop("source_nodes must be a non-empty character vector.")
  }
  if (!is.character(target_nodes) || !length(target_nodes)) {
    stop("target_nodes must be a non-empty character vector.")
  }

  ## ==============================================================
  ## Step 2 - Build network graph and map nodes
  ## ==============================================================
  g <- igraph::graph_from_data_frame(network_df[, seq_len(2)], directed = TRUE)
  graph_nodes <- igraph::V(g)$name

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

  source_nodes <- normalize_node_ids(source_nodes)
  target_nodes <- normalize_node_ids(target_nodes)

  mapped_sources <- intersect(source_nodes, graph_nodes)
  mapped_targets <- intersect(target_nodes, graph_nodes)

  if (!length(mapped_sources)) stop("No source_nodes mapped onto the network.")
  if (!length(mapped_targets)) stop("No target_nodes mapped onto the network.")

  unmapped_sources <- setdiff(source_nodes, graph_nodes)
  unmapped_targets <- setdiff(target_nodes, graph_nodes)

  if (length(unmapped_sources))
    message("Some source_nodes were not mapped: ", paste(unmapped_sources, collapse = ", "))
  if (length(unmapped_targets))
    message("Some target_nodes were not mapped: ", paste(unmapped_targets, collapse = ", "))

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

  ## ==============================================================
  ## Step 3 - Compute distances and extract valid pairs
  ## ==============================================================
  dist_mat <- igraph::distances(g, mode = mode)
  ms <- intersect(mapped_sources, rownames(dist_mat))
  mt <- intersect(mapped_targets, colnames(dist_mat))

  subdist <- dist_mat[ms, mt, drop = FALSE]
  ok <- is.finite(subdist) & (subdist <= distance_threshold)
  if (!any(ok)) {
    message("No valid source-target pairs found within distance_threshold = ", distance_threshold)
    if (output == "path_list") return(list())
    empty_out <- matrix(character(0), nrow = 0, ncol = 3)
    colnames(empty_out) <- c("source", "target", "interaction_type")
    return(empty_out)
  }

  valid_idx <- which(ok, arr.ind = TRUE)
  valid_pairs <- data.frame(
    source   = rownames(subdist)[valid_idx[, 1]],
    target   = colnames(subdist)[valid_idx[, 2]],
    distance = subdist[valid_idx],
    stringsAsFactors = FALSE
  )

  message("Found ", nrow(valid_pairs),
          " valid source-target pairs (distance <= ", distance_threshold,
          "). Proceeding to shortest-path identification")

  ## ==============================================================
  ## Step 4 - Split pairs by 'magic number' 6 if threshold > 6
  ## ==============================================================
  ENUM_CUTOFF <- 6L     # Empirical cutoff for enumerating all shortest-path ties

  if (isTRUE(betweenness)) {
    ## Betweenness selection requires enumerating shortest-path ties.
    ## Here we always compute all shortest paths for all valid pairs.
    pairs_short <- valid_pairs
    pairs_long  <- valid_pairs[0, ]
  } else {
    use_split <- is.finite(distance_threshold) && distance_threshold > ENUM_CUTOFF
    if (use_split) {
      pairs_short <- subset(valid_pairs, distance <= ENUM_CUTOFF)
      pairs_long  <- subset(valid_pairs, distance > ENUM_CUTOFF)
    } else {
      pairs_short <- valid_pairs
      pairs_long  <- valid_pairs[0, ]
    }
  }

  ## ==============================================================
  ## Step 5 - Helper: collect paths for a set of pairs with a solver
  ## ==============================================================
  collect_paths <- function(pairs_df, solver = c("all", "first")) {
    solver <- match.arg(solver)
    pl <- list()
    if (!nrow(pairs_df)) return(pl)

    by_src <- split(pairs_df, pairs_df$source)

    for (src in names(by_src)) {
      tgt_vec <- by_src[[src]]$target

      if (solver == "all") {
        res <- igraph::all_shortest_paths(g, from = src, to = tgt_vec, mode = mode)$res
        if (!length(res)) next
        for (vseq in res) {
          ids <- igraph::as_ids(vseq)
          if (!length(ids)) next
          tgt <- ids[length(ids)]
          key <- paste0(src, "_", tgt)
          pl[[key]] <- c(pl[[key]], list(ids))
        }
      } else { # "first"
        sp <- igraph::shortest_paths(g, from = src, to = tgt_vec, mode = mode)$vpath
        if (!length(sp)) next
        for (i in seq_along(sp)) {
          ids <- igraph::as_ids(sp[[i]])
          if (!length(ids)) next
          tgt <- tgt_vec[i]
          key <- paste0(src, "_", tgt)
          pl[[key]] <- c(pl[[key]], list(ids))
        }
      }
    }
    pl
  }

  ## ==============================================================
  ## Step 6 - Compute paths for both subsets
  ## ==============================================================
  path_list <- list()
  if (nrow(pairs_short)) {
    pl_short <- collect_paths(pairs_short, solver = "all")
    if (length(pl_short)) path_list[names(pl_short)] <- pl_short
  }
  if (nrow(pairs_long)) {
    pl_long <- collect_paths(pairs_long, solver = "first")
    if (length(pl_long)) {
      for (k in names(pl_long)) {
        path_list[[k]] <- c(path_list[[k]], pl_long[[k]])
      }
    }
  }

  if (!length(path_list)) {
    message("No valid paths found between mapped sources and targets.")
    if (output == "path_list") return(list())
    empty_out <- matrix(character(0), nrow = 0, ncol = 3)
    colnames(empty_out) <- c("source", "target", "interaction_type")
    return(empty_out)
  }

  ## ==============================================================
  ## Step 7 - Optional: filter paths by KEGG-supported reactant pairs
  ## ==============================================================
  if (isTRUE(reactant_filter)) {
    message(
      "Reactant-pair filtering processing. ",
      "This step may be slow for large path sets."
    )

    ## Flatten path_list into a data frame
    path_rows <- lapply(names(path_list), function(nm) {
      lapply(seq_along(path_list[[nm]]), function(i) {
        data.frame(
          pair_id = nm,
          path_id = i,
          path    = paste(path_list[[nm]][[i]], collapse = ","),
          stringsAsFactors = FALSE
        )
      })
    })
    path_df <- do.call(rbind, unlist(path_rows, recursive = FALSE))
    rownames(path_df) <- NULL

    ## Extract consecutive compound-compound pairs from each path
    path_df$cpd_pairs <- NA_character_

    for (i in seq_len(nrow(path_df))) {
      nodes_i <- trimws(unlist(strsplit(path_df$path[i], ",")))
      cpd_i   <- nodes_i[grepl("^cpd:C\\d{5}$", nodes_i)]

      if (length(cpd_i) >= 2) {
        pair_vec <- vapply(seq_len(length(cpd_i) - 1), function(j) {
          paste(cpd_i[j], cpd_i[j + 1], sep = "_")
        }, character(1))
        path_df$cpd_pairs[i] <- paste(pair_vec, collapse = ",")
      }
    }

    ## Collect all unique compound pairs involved in the candidate paths
    pair_vec_all <- unique(unlist(strsplit(na.omit(path_df$cpd_pairs), ",")))
    pair_vec_all <- trimws(pair_vec_all)
    pair_vec_all <- pair_vec_all[pair_vec_all != ""]

    ## Validate each compound pair against KEGG reactions
    pair_valid_df <- data.frame(
      cpd_pair = pair_vec_all,
      valid    = "no",
      stringsAsFactors = FALSE
    )

    for (i in seq_len(nrow(pair_valid_df))) {
      pair_i <- pair_valid_df$cpd_pair[i]
      comps  <- unlist(strsplit(pair_i, "_"))
      c1 <- comps[1]
      c2 <- comps[2]

      rxn1 <- tryCatch(KEGGREST::keggLink("reaction", c1), error = function(e) NULL)
      rxn2 <- tryCatch(KEGGREST::keggLink("reaction", c2), error = function(e) NULL)

      if (is.null(rxn1) || is.null(rxn2) || length(rxn1) == 0 || length(rxn2) == 0) {
        next
      }

      rxn1_ids <- unique(as.character(rxn1))
      rxn2_ids <- unique(as.character(rxn2))
      common_rxn <- intersect(rxn1_ids, rxn2_ids)

      if (length(common_rxn) == 0) {
        next
      }

      pair_is_valid <- FALSE

      for (rx in common_rxn) {
        info <- tryCatch(KEGGREST::keggGet(rx), error = function(e) NULL)
        if (is.null(info) || length(info) == 0 || is.null(info[[1]]$EQUATION)) {
          next
        }

        eq <- info[[1]]$EQUATION[1]

        if (grepl("<=>", eq, fixed = TRUE)) {
          parts <- strsplit(eq, "<=>", fixed = TRUE)[[1]]
        } else if (grepl("=>", eq, fixed = TRUE)) {
          parts <- strsplit(eq, "=>", fixed = TRUE)[[1]]
        } else if (grepl("<=", eq, fixed = TRUE)) {
          parts <- strsplit(eq, "<=", fixed = TRUE)[[1]]
        } else {
          next
        }

        if (length(parts) != 2) {
          next
        }

        left_ids  <- unique(regmatches(parts[1], gregexpr("C\\d{5}", parts[1]))[[1]])
        right_ids <- unique(regmatches(parts[2], gregexpr("C\\d{5}", parts[2]))[[1]])

        left_ids  <- paste0("cpd:", left_ids)
        right_ids <- paste0("cpd:", right_ids)

        if ((c1 %in% left_ids && c2 %in% right_ids) ||
            (c2 %in% left_ids && c1 %in% right_ids)) {
          pair_is_valid <- TRUE
          break
        }
      }

      if (pair_is_valid) {
        pair_valid_df$valid[i] <- "yes"
      }
    }

    ## Mark each path as valid only if all compound pairs are valid
    ## Paths without any compound-compound transition are kept
    path_df$match <- NA_character_

    for (i in seq_len(nrow(path_df))) {
      if (is.na(path_df$cpd_pairs[i]) || path_df$cpd_pairs[i] == "") {
        path_df$match[i] <- "yes"
        next
      }

      vals <- trimws(unlist(strsplit(path_df$cpd_pairs[i], ",")))
      valid_flags <- pair_valid_df$valid[match(vals, pair_valid_df$cpd_pair)]

      if (length(valid_flags) > 0 && all(!is.na(valid_flags) & valid_flags == "yes")) {
        path_df$match[i] <- "yes"
      } else {
        path_df$match[i] <- "no"
      }
    }

    ## Rebuild filtered path_list
    path_df_valid <- path_df[path_df$match == "yes", , drop = FALSE]
    path_list_valid <- list()

    for (i in seq_len(nrow(path_df_valid))) {
      nm <- path_df_valid$pair_id[i]
      path_nodes <- trimws(unlist(strsplit(path_df_valid$path[i], ",")))

      if (is.null(path_list_valid[[nm]])) {
        path_list_valid[[nm]] <- list(path_nodes)
      } else {
        path_list_valid[[nm]] <- c(path_list_valid[[nm]], list(path_nodes))
      }
    }

    path_list <- path_list_valid

    if (!length(path_list)) {
      message("No valid paths remained after reactant-pair filtering.")
      if (output == "path_list") return(list())
      empty_out <- matrix(character(0), nrow = 0, ncol = 3)
      colnames(empty_out) <- c("source", "target", "interaction_type")
      return(empty_out)
    }
  }

  ## ==============================================================
  ## Step 8 - Optional: select a betweenness-ranked shortest path per pair
  ## ==============================================================
  if (isTRUE(betweenness)) {
    message("Betweenness selection processing")

    ## In 'all' mode the graph is treated as undirected for centrality scoring
    directed_flag <- !(identical(mode, "all"))

    ## Compute node betweenness on the full network graph
    bw_vec <- igraph::betweenness(
      graph      = g,
      directed   = directed_flag,
      weights    = NULL,
      normalized = TRUE
    )

    ## Build a simple lookup: node_id -> betweenness score
    bw_map <- setNames(as.numeric(bw_vec), names(bw_vec))

    ## Score a candidate path by summing betweenness across its nodes
    score_path <- function(p) {
      if (length(p) == 0) return(-Inf)
      sum(bw_map[p], na.rm = TRUE)
    }

    ## For each source-target pair, keep only the top-ranked shortest path
    path_list_bw <- list()
    for (k in names(path_list)) {
      cand <- path_list[[k]]
      if (!length(cand)) next

      ## If there is only one shortest path, keep it as-is
      if (length(cand) == 1) {
        path_list_bw[[k]] <- cand
      } else {
        scores <- vapply(cand, score_path, numeric(1))
        best_i <- which(scores == max(scores, na.rm = TRUE))[1]
        path_list_bw[[k]] <- list(cand[[best_i]])
      }
    }

    ## Replace the original path list with the filtered betweenness-ranked paths
    path_list <- path_list_bw
  }

  ## ==============================================================
  ## Step 9 - Output type handling
  ## ==============================================================
  if (output == "path_list") {
    if (!name) return(path_list)

    # ---- Name conversion ----
    message("Converting node identifiers to names or symbols. ",
            "This step may take some time depending on network size")

    all_nodes_in_paths <- unique(unlist(unlist(path_list, recursive = FALSE)))
    annot <- setNames(all_nodes_in_paths, all_nodes_in_paths)

    annot <- vapply(all_nodes_in_paths, function(nid) {
      Sys.sleep(0.1)
      info <- tryCatch(KEGGREST::keggGet(nid), error = function(e) NULL)
      if (is.null(info) || is.null(info[[1]])) return(nid)
      e <- info[[1]]
      if (grepl("^cpd:|^dr:|^gl:", nid)) {
        if (!is.null(e$NAME)) {
          sub(";.*$", "", e$NAME[1])
        } else {
          nid
        }
      } else {
        if (!is.null(e$SYMBOL)) {
          strsplit(e$SYMBOL, ",")[[1]][1]
        } else if (!is.null(e$NAME)) {
          e$NAME[1]
        } else {
          nid
        }
      }
    }, character(1), USE.NAMES = TRUE)

    # Prevent duplicate display names after KEGG conversion
    annot_unique <- annot
    names(annot_unique) <- names(annot)
    dup_vals <- duplicated(annot_unique)
    if (any(dup_vals)) {
      annot_unique[] <- make.unique(annot_unique)
    }

    path_list_named <- lapply(path_list, function(pl) lapply(pl, function(p) unname(annot_unique[p])))

    message("Node name conversion completed successfully.")
    return(path_list_named)
  }

  ## ==============================================================
  ## Step 10 - network_matrix output (MetaPathNet-style edge list)
  ## ==============================================================
  edges <- list()
  for (k in names(path_list)) {
    for (pth in path_list[[k]]) {
      if (length(pth) > 1) {
        edges[[length(edges) + 1]] <- data.frame(
          source = head(pth, -1),
          target = tail(pth, -1),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  if (!length(edges)) {
    message("No edges could be formed from computed paths.")
    empty_out <- matrix(character(0), nrow = 0, ncol = 3)
    colnames(empty_out) <- c("source", "target", "interaction_type")
    return(empty_out)
  }

  edge_df <- unique(do.call(rbind, edges))

  key_net <- paste(network_df$source, network_df$target, sep = "_")
  lookup  <- setNames(network_df$interaction_type, key_net)

  keys_forward <- paste(edge_df$source, edge_df$target, sep = "_")
  edge_df$interaction_type <- lookup[keys_forward]

  if (identical(mode, "all")) {
    na_idx <- is.na(edge_df$interaction_type)

    if (any(na_idx)) {
      keys_reverse <- paste(edge_df$target[na_idx], edge_df$source[na_idx], sep = "_")
      edge_df$interaction_type[na_idx] <- lookup[keys_reverse]
    }
  }

  out <- as.matrix(edge_df[, c("source", "target", "interaction_type")])
  colnames(out) <- c("source", "target", "interaction_type")
  rownames(out) <- NULL
  out
}
