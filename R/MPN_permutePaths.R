#' Permutation Test of Shortest-Path Structure Between Node Sets
#'
#' Performs a permutation-based test on the mean shortest-path length
#' between a set of source nodes and target nodes within a KEGG-derived
#' network. The observed mean shortest-path distance is compared to a null
#' distribution obtained by repeatedly sampling random source and target sets
#' of the same size and node-type composition from the same network.
#'
#' @param network_table A 3-column matrix or data.frame representing
#'   directed network edges. Columns must correspond to
#'   \code{source}, \code{target}, and \code{interaction_type}.
#' @param source_nodes Character vector of KEGG node IDs used as
#'   \emph{sources} in the permutation test.
#' @param target_nodes Character vector of KEGG node IDs used as
#'   \emph{targets} in the permutation test.
#' @param n_perm Integer. Number of permutations used to build the null
#'   distribution of mean shortest-path lengths (default: \code{10000}).
#' @param mode Character. Direction parameter passed to
#'   \code{igraph::distances()}. One of:
#'   \describe{
#'     \item{\code{"out"}}{Distances follow edge direction.}
#'     \item{\code{"in"}}{Distances are computed against edge direction (reverse direction).}
#'     \item{\code{"all"}}{Distances ignore edge direction (graph treated as undirected).}
#'   }
#' @param plot Logical (default: \code{FALSE}). If \code{TRUE}, a histogram
#'   and density curve of permuted mean shortest-path lengths is plotted,
#'   with a vertical line for the observed mean and an annotated p-value.
#' @param return_perm Logical (default: \code{FALSE}). If \code{FALSE},
#'   only summary statistics are returned. If \code{TRUE}, the raw
#'   permutation means are returned in addition to the summary.
#'
#' @return
#' If \code{return_perm = FALSE}, a data.frame with one row:
#' \describe{
#'   \item{observed_sp}{Observed mean shortest-path length
#'         (finite distances only; unreachable pairs are excluded).}
#'   \item{mean_perm_sp}{Mean of permuted mean shortest-path lengths
#'         under the null (finite values only).}
#'   \item{p_value}{Empirical one-sided p-value
#'         \eqn{P(\text{random mean} < \text{observed mean})}.}
#'   \item{n_perm}{Number of permutations used.}
#' }
#'
#' If \code{return_perm = TRUE}, a list with:
#' \describe{
#'   \item{stats}{The summary data.frame described above.}
#'   \item{random_sp}{Numeric vector of permuted mean shortest-path
#'         lengths (finite values only).}
#' }
#'
#' @details
#' Internally, the function:
#' \itemize{
#'   \item identifies supported node types in \code{network_table},
#'         currently KEGG Orthology nodes (\code{"K"} prefix) and
#'         KEGG compound nodes (\code{"cpd:"} prefix);
#'   \item maps the user-supplied source and target nodes onto the network;
#'   \item removes unreachable pairs (distance = \code{Inf}) when calculating
#'         observed and permuted mean distances;
#'   \item derives an empirical one-sided p-value from the proportion of
#'         permuted means that are smaller than the observed mean.
#' }
#'
#' Source and target node sets are user-defined and may contain supported KEGG
#' node types according to the biological question being addressed. A smaller
#' mean shortest-path length indicates tighter connectivity. Therefore, the
#' reported p-value tests whether random source-target sets tend to be
#' \emph{more tightly connected} than the observed set
#' (\eqn{P(\text{random mean} < \text{observed mean})}).
#'
#' @section Plotting:
#' If \code{plot = TRUE}, a \pkg{ggplot2}-based histogram + density overlay
#' of the permutation distribution is printed. For custom visualisation,
#' set \code{return_perm = TRUE} and use the \code{random_sp} vector
#' together with \code{stats} in your own plotting code.
#'
#' @examples
#' ## 1) Prepare a mixed human–E. coli network focused on tryptophan-related pathways
#' ## Human: tryptophan metabolism + selected signaling pathways
#' ## Load the precomputed example network; see MPN_keggNetwork() for its construction
#' data(MetaPathNet_example_network)
#'
#' ## Combine host and microbe networks into one edge list
#' net_trp_mixed <- MetaPathNet_example_network
#'
#' ## 2) Define source and target node sets for the Trp/kynurenine/serotonin context
#' source_nodes <- c(
#'   "K00486",  # kynurenine 3-monooxygenase
#'   "K00453",  # tryptophan 2,3-dioxygenase
#'   "K00463",  # indoleamine 2,3-dioxygenase
#'   "K01667",  # kynureninase
#'   "K01696",  # tryptophan synthase beta chain
#'   "K01695"   # tryptophan synthase alpha chain
#' )
#'
#' target_nodes <- c(
#'   "cpd:C00078",  # L-tryptophan
#'   "cpd:C00328",  # L-kynurenine
#'   "cpd:C00780",  # serotonin
#'   "cpd:C00463",  # indole
#'   "cpd:C00108"   # anthranilate
#' )
#'
#' ## 3) Permutation test on mean shortest-path length between the two node sets
#' perm_trp <- MPN_permutePaths(
#'   network_table = net_trp_mixed, # mixed hsa/eco Trp + immune-related network
#'   source_nodes  = source_nodes,
#'   target_nodes  = target_nodes,
#'   n_perm        = 100,           # increase (e.g. 10000) for more stable inference
#'   mode          = "out",
#'   plot          = FALSE,         # show permutation distribution
#'   return_perm   = TRUE           # keep individual permuted means
#' )
#'
#' ## 4) View summary statistics
#' perm_trp$stats
#'
#' @seealso
#' \code{\link{MPN_distances}} for shortest-path distance matrices,
#' \code{\link{MPN_shortestPaths}} for explicit shortest paths,
#' \code{\link{MPN_keggNetwork}}, \code{\link{MPN_crossSpeciesNetwork}} for
#' KEGG-based network construction.
#'
#' @importFrom igraph distances graph_from_data_frame V
#' @importFrom ggplot2 ggplot aes geom_histogram after_stat geom_density
#'   geom_vline annotate labs theme_minimal element_text theme
#'   coord_cartesian margin
#' @export
MPN_permutePaths <- function(
    network_table,
    source_nodes,
    target_nodes,
    n_perm      = 10000,
    mode        = c("out", "in", "all"),
    plot        = FALSE,
    return_perm = FALSE
) {
  #==============================================================
  # Step 1 — Input checks and graph construction
  #==============================================================
  mode <- match.arg(mode)

  # Accept matrix/data.frame input and remove duplicated edges
  if (!is.matrix(network_table)) {
    network_table <- as.matrix(network_table)
  }
  network_table <- unique(network_table)

  # Basic structure checks
  if (ncol(network_table) < 3L) {
    stop("network_table must be a matrix with at least 3 columns: source, target, interaction_type.")
  }

  if (!is.character(source_nodes) || !length(source_nodes)) {
    stop("source_nodes must be a non-empty character vector.")
  }
  if (!is.character(target_nodes) || !length(target_nodes)) {
    stop("target_nodes must be a non-empty character vector.")
  }

  # Build directed graph from source/target columns
  g <- igraph::graph_from_data_frame(
    as.data.frame(network_table[, seq_len(2)], stringsAsFactors = FALSE),
    directed = TRUE
  )
  graph_nodes <- igraph::V(g)$name

  #==============================================================
  # Step 2 — Define KO and compound node pools in the network
  #==============================================================
  all_nodes <- unique(graph_nodes)
  ko_nodes  <- all_nodes[grep("^K", all_nodes)]
  cpd_nodes <- all_nodes[grep("^cpd:", all_nodes)]

  # Require at least one supported node type in the network
  if (!length(ko_nodes) && !length(cpd_nodes)) {
    stop("No supported node types detected in the network_table. Expected KO nodes ('Kxxxxx') and/or compound nodes ('cpd:Cxxxxx').")
  }

  #==============================================================
  # Step 3 — Map user-specified source/target nodes and define type-matched pools
  #==============================================================
  source_nodes <- unique(source_nodes)
  target_nodes <- unique(target_nodes)

  # Harmonise KEGG compound IDs
  normalize_compound_ids <- function(x) {
    is_compound <- grepl("^(cpd:)?[Cc][0-9]{5}$", x)
    if (any(is_compound)) {
      # Remove optional "cpd:" prefix, standardise case, then rebuild as "cpd:Cxxxxx"
      core_ids <- toupper(gsub("^(cpd:)?", "", x[is_compound]))
      x[is_compound] <- paste0("cpd:", core_ids)
    }
    x
  }

  # Standardise source and target compound IDs before mapping to the network
  source_nodes <- normalize_compound_ids(source_nodes)
  target_nodes <- normalize_compound_ids(target_nodes)

  # Global mapping check across all provided source + target IDs
  input_nodes_all <- unique(c(source_nodes, target_nodes))
  unmapped_nodes_all <- setdiff(input_nodes_all, graph_nodes)
  if (length(unmapped_nodes_all)) {
    show_ids <- head(unmapped_nodes_all, 10L)
    message(
      length(unmapped_nodes_all),
      " input nodes were not mapped onto the network. Showing first 10: ",
      paste(show_ids, collapse = ", ")
    )
  }

  # Keep only source and target nodes that are actually present in the graph
  mapped_sources <- intersect(source_nodes, graph_nodes)
  mapped_targets <- intersect(target_nodes, graph_nodes)

  # Stop if either side has no valid mapped nodes left after filtering
  if (!length(mapped_sources)) {
    stop("None of the source_nodes mapped onto nodes in the network.")
  }
  if (!length(mapped_targets)) {
    stop("None of the target_nodes mapped onto nodes in the network.")
  }

  # Classify mapped nodes into supported MetaPathNet node types
  get_node_type <- function(x) {
    out <- rep(NA_character_, length(x))
    out[grepl("^K\\d{5}$", x)] <- "KO"
    out[grepl("^cpd:C\\d{5}$", x)] <- "compound"
    out
  }

  # Infer source-side and target-side node-type composition from mapped inputs
  source_types <- get_node_type(mapped_sources)
  target_types <- get_node_type(mapped_targets)

  # Stop if mapped source/target nodes include unsupported node types
  if (any(is.na(source_types))) {
    bad <- mapped_sources[is.na(source_types)]
    stop(
      "Mapped source_nodes contain unsupported node types. Supported types are KO ('Kxxxxx') and compound ('cpd:Cxxxxx'). ",
      "First unsupported source_nodes: ", paste(head(bad, 10), collapse = ", ")
    )
  }

  if (any(is.na(target_types))) {
    bad <- mapped_targets[is.na(target_types)]
    stop(
      "Mapped target_nodes contain unsupported node types. Supported types are KO ('Kxxxxx') and compound ('cpd:Cxxxxx'). ",
      "First unsupported target_nodes: ", paste(head(bad, 10), collapse = ", ")
    )
  }

  # Count how many KO and compound nodes are present on each side
  source_type_table <- table(source_types)
  target_type_table <- table(target_types)

  n_source_KO       <- if ("KO" %in% names(source_type_table)) unname(source_type_table["KO"]) else 0L
  n_source_compound <- if ("compound" %in% names(source_type_table)) unname(source_type_table["compound"]) else 0L
  n_target_KO       <- if ("KO" %in% names(target_type_table)) unname(target_type_table["KO"]) else 0L
  n_target_compound <- if ("compound" %in% names(target_type_table)) unname(target_type_table["compound"]) else 0L

  # Check that the network contains enough nodes of each required type
  if (n_source_KO > length(ko_nodes) || n_target_KO > length(ko_nodes)) {
    stop("The network does not contain enough KO nodes to support permutation sampling with the requested source/target composition.")
  }
  if (n_source_compound > length(cpd_nodes) || n_target_compound > length(cpd_nodes)) {
    stop("The network does not contain enough compound nodes to support permutation sampling with the requested source/target composition.")
  }

  # Sample a random node set with the same KO/compound composition
  sample_type_matched_nodes <- function(n_ko, n_cpd, ko_pool, cpd_pool) {
    out <- character(0)

    if (n_ko > 0L) {
      out <- c(out, sample(ko_pool, n_ko))
    }
    if (n_cpd > 0L) {
      out <- c(out, sample(cpd_pool, n_cpd))
    }

    out
  }

  #==============================================================
  # Step 4 — Compute observed shortest-path statistic
  #==============================================================
  dist_obs <- igraph::distances(
    g,
    v    = mapped_sources,
    to   = mapped_targets,
    mode = mode
  )

  # Use only finite distances (reachable pairs) for the observed mean
  finite_obs <- dist_obs[is.finite(dist_obs)]
  if (!length(finite_obs)) {
    stop("No finite shortest paths between the specified source and target nodes.")
  }

  observed_sp <- mean(finite_obs)

  #==============================================================
  # Step 5 — Permutation loop (null distribution)
  #==============================================================
  n_perm <- as.integer(n_perm)
  if (is.na(n_perm) || n_perm < 1L) {
    stop("n_perm must be a positive integer.")
  }

  random_sp <- numeric(n_perm)

  # Re-sample source and target sets with the same node-type composition and compute mean distances
  for (b in seq_len(n_perm)) {
    rand_sources <- sample_type_matched_nodes(
      n_ko   = n_source_KO,
      n_cpd  = n_source_compound,
      ko_pool  = ko_nodes,
      cpd_pool = cpd_nodes
    )

    rand_targets <- sample_type_matched_nodes(
      n_ko   = n_target_KO,
      n_cpd  = n_target_compound,
      ko_pool  = ko_nodes,
      cpd_pool = cpd_nodes
    )

    dist_rand <- igraph::distances(
      g,
      v    = rand_sources,
      to   = rand_targets,
      mode = mode
    )

    finite_rand <- dist_rand[is.finite(dist_rand)]
    if (!length(finite_rand)) {
      random_sp[b] <- NA_real_
    } else {
      random_sp[b] <- mean(finite_rand)
    }
  }

  #==============================================================
  # Step 6 — Summarise permutation results and empirical p-value
  #==============================================================
  valid_rand <- random_sp[is.finite(random_sp)]
  if (!length(valid_rand)) {
    warning("All permuted shortest-path means are NA. Returning NA statistics.")
    mean_perm_sp <- NA_real_
    p_value      <- NA_real_
  } else {
    mean_perm_sp <- mean(valid_rand)
    p_value <- (sum(valid_rand < observed_sp) + 1) / (length(valid_rand) + 1)
  }

  stats_df <- data.frame(
    observed_sp  = observed_sp,
    mean_perm_sp = mean_perm_sp,
    p_value      = p_value,
    n_perm       = n_perm,
    stringsAsFactors = FALSE
  )

  #==============================================================
  # Step 7 — Optional plot: permutation distribution
  #==============================================================
  if (plot && length(valid_rand)) {
    df_plot <- data.frame(random_sp = valid_rand)

    p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = random_sp)) +
      ggplot2::geom_histogram(
        ggplot2::aes(y = ggplot2::after_stat(density)),
        bins  = 40,
        fill  = "#dbe4f0",
        color = "white",
        alpha = 0.8
      ) +
      ggplot2::geom_density(
        linewidth = 1,
        color     = "#34495e"
      ) +
      ggplot2::geom_vline(
        xintercept = observed_sp,
        color      = "#d35400",
        linewidth  = 1
      ) +
      ggplot2::annotate(
        "text",
        x     = observed_sp,
        y     = Inf,
        vjust = 1.2,
        hjust = 1.1,
        label = sprintf("Observed = %.3f\np = %.4f", observed_sp, p_value),
        size  = 4,
        color = "#d35400"
      )  +
      ggplot2::coord_cartesian(clip = "off") +
      ggplot2::labs(
        title = "Permutation distribution of mean shortest path length",
        x     = "Mean shortest path length (permutation)",
        y     = "Density"
      ) +
      ggplot2::theme_minimal(base_size = 13) +
      ggplot2::theme(
        plot.title       = ggplot2::element_text(hjust = 0.5, face = "bold"),
        panel.grid.minor = ggplot2::element_blank(),
        plot.margin      = ggplot2::margin(t = 10, r = 10, b = 10, l = 35)
      )

    print(p)
  }

  #==============================================================
  # Step 8 — Return statistics, optionally with permutation values
  #==============================================================
  if (return_perm) {
    return(list(
      stats     = stats_df,
      random_sp = valid_rand
    ))
  }

  stats_df
}
