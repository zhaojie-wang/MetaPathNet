#' Pathway Enrichment Analysis (Over-representation Analysis) on KEGG-Based Node Sets
#'
#' Performs pathway enrichment analysis (over-representation analysis, ORA)
#' using the hypergeometric test with Benjamini–Hochberg correction on nodes
#' from KEGG-derived networks. Supports KO-level, compound-level, or integrated
#' (KO + compound) analyses.
#'
#' @param test_set Character vector of node IDs (e.g. KEGG compounds, KO IDs),
#'   or a network matrix/data.frame, in which case unique nodes from the first
#'   two columns are used.
#' @param background_set Character vector of node IDs or a network matrix/data.frame
#'   defining the background universe for enrichment analysis.
#' @param pathway_set Either a character vector of KEGG reference pathway IDs
#'   or a named list where each element contains the node members of one pathway.
#' @param entity Character string specifying which entity type is tested:
#'   \code{"ko"} (KOs only), \code{"compound"} (compounds only), or
#'   \code{"integrated"} (both). Default \code{"ko"}.
#' @param pvalueCutoff Numeric. Significance cutoff for adjusted p-values
#'   (default \code{0.05}, after Benjamini–Hochberg correction).
#' @param add_hits Logical. If \code{TRUE}, adds a \code{Test_Hits} column
#'   listing the test nodes found in each pathway (default \code{FALSE}).
#'
#' @details
#' If \code{test_set} or \code{background_set} is a matrix/data.frame, unique
#' nodes are taken from the first two columns. The test set is intersected with
#' the background (\code{test_set <= background_set}) before enrichment.
#'
#' Enrichment is computed with the hypergeometric distribution (\code{phyper})
#' followed by Benjamini–Hochberg FDR adjustment.
#'
#' @return A data frame with columns:
#' \describe{
#'   \item{Pathway_ID}{KEGG pathway ID (e.g. \code{"map00380"}).}
#'   \item{Description}{Pathway description from KEGG.}
#'   \item{bg_overlap}{Number of background nodes in the pathway.}
#'   \item{test_overlap}{Number of test nodes in the pathway.}
#'   \item{Test_Hits}{(Optional) Semicolon-separated test nodes in the pathway.}
#'   \item{p_value}{Raw hypergeometric p-value.}
#'   \item{p_adj}{Adjusted p-value (Benjamini–Hochberg).}
#'   \item{Enrichment_Ratio}{\code{(test_overlap / |test_set|) / (bg_overlap / |background_set|)}.}
#'   \item{significance}{\code{"significant"} or \code{"NS"} based on \code{pvalueCutoff}.}
#' }
#'
#' @note
#' Pathways with zero test hits are excluded automatically. Very small test
#' sets or test sets almost identical to the background may yield uninformative
#' results (p-values close to 1).
#'
#' @examples
#' ## Load the precomputed example network; see MPN_keggNetwork() for its construction
#' data(MetaPathNet_example_network)
#'
#' ## Define a small test set from the example network
#' test_nodes <- c("K00463", "K00453", "cpd:C00078", "cpd:C02700")
#'
#' ## Run integrated pathway over-representation analysis
#' enrich_res <- MPN_enrichPathway(
#'   test_set       = test_nodes,
#'   background_set = MetaPathNet_example_network,
#'   pathway_set    = c("map00380", "map00400"),
#'   entity         = "integrated",
#'   pvalueCutoff   = 1,
#'   add_hits       = TRUE
#' )
#'
#' ## Preview enrichment results
#' head(enrich_res)
#'
#' @seealso \code{\link{MPN_egoNetwork}}, \code{\link{MPN_getPathIDs}},
#'   \code{\link[KEGGREST]{keggLink}}, \code{\link[KEGGREST]{keggCompounds}}
#' @importFrom KEGGREST keggLink keggCompounds
#' @importFrom stats phyper p.adjust
#' @importFrom utils read.delim
#' @importFrom curl curl_fetch_memory new_handle
#' @export
MPN_enrichPathway <- function(test_set, background_set, pathway_set,
                             entity = c("ko", "compound", "integrated"),
                             pvalueCutoff = 0.05, add_hits = FALSE) {

  #==============================================================
  # Step 1 — Match arguments and normalise input sets
  #==============================================================
  entity <- match.arg(entity)

  # Accept edge tables (matrix / data.frame) for bg_set and test_set
  if (is.matrix(background_set)) {
    bg_vals <- unique(as.vector(background_set[, seq_len(2)]))

    if (entity == "compound") {
      background_set <- bg_vals[grepl("^cpd", bg_vals)]
    } else if (entity == "ko") {
      background_set <- bg_vals[grepl("^K",  bg_vals)]
    } else {
      background_set <- bg_vals
    }

  } else if (is.data.frame(background_set)) {
    bg_vals <- unique(unlist(background_set[, seq_len(2)]))

    if (entity == "compound") {
      background_set <- bg_vals[grepl("^cpd", bg_vals)]
    } else if (entity == "ko") {
      background_set <- bg_vals[grepl("^K", bg_vals)]
    } else {
      background_set <- bg_vals
    }
  }

  if (is.matrix(test_set) || is.data.frame(test_set)) {
    test_set <- unique(unlist(test_set[, seq_len(2)]))
  }

  # Strip prefixes consistently
  test_set       <- sub("^(cpd:|ko:)", "", test_set)
  background_set <- sub("^(cpd:|ko:)", "", background_set)

  # Enforce ORA assumption: test_set must be a subset of bg_set
  if (!all(test_set %in% background_set)) {
    message("Adjusting: intersecting test_set with background_set to satisfy ORA assumptions.")
    test_set <- intersect(test_set, background_set)
  }

  # Basic size / power checks
  if (setequal(sort(test_set), sort(background_set))) {
    message("Test set is identical to background set. Enrichment analysis will not be informative (all p-values ~1).")
  }
  if (length(test_set) < 5) {
    message("Test set is very small (<5 nodes). Results may lack statistical power.")
  } else if (length(test_set) / length(background_set) < 0.05) {
    message("Test set is much smaller than background set (<5%). Results may lack statistical power. Consider a larger test set.")
  }

  #==============================================================
  # Step 2 — Build pathway membership list (KOs + compounds)
  #==============================================================
  if (is.character(pathway_set) && !is.list(pathway_set)) {
    reference_pathways <- paste0("map", sub("^[a-zA-Z]+", "", pathway_set))
    pl <- vector("list", length(reference_pathways))
    names(pl) <- pathway_set

    for (i in seq_along(reference_pathways)) {
      kos <- tryCatch(
        sub("^ko:", "", KEGGREST::keggLink("ko", reference_pathways[i])),
        error = function(e) character(0)
      )
      comps <- tryCatch(
        KEGGREST::keggCompounds(reference_pathways[i]),
        error = function(e) character(0)
      )
      pl[[i]] <- unique(c(kos, comps))
    }
    pathway_list <- pl

  } else if (is.list(pathway_set)) {
    pathway_list <- pathway_set

  } else {
    stop("pathway_set must be either a character vector of pathway IDs or a named list.")
  }

  #==============================================================
  # Step 3 — Filter pathway members and inputs by entity type
  #==============================================================
  if (entity == "ko") {
    pathway_list <- lapply(pathway_list, function(x) x[grepl("^K", x)])
    test_set       <- test_set[grepl("^K", test_set)]
    background_set <- background_set[grepl("^K", background_set)]

  } else if (entity == "compound") {
    pathway_list <- lapply(pathway_list, function(x) x[grepl("^C", x)])
    test_set       <- test_set[grepl("^C", test_set)]
    background_set <- background_set[grepl("^C", background_set)]
  }

  #==============================================================
  # Step 4 — Hypergeometric ORA per pathway
  #==============================================================
  N <- length(background_set)   # background size
  n <- length(test_set)         # test-set size

  results_list <- lapply(names(pathway_list), function(pathway) {
    pathway_members <- intersect(pathway_list[[pathway]], background_set)
    bg_overlap      <- length(pathway_members)
    test_in_pathway <- intersect(pathway_list[[pathway]], test_set)
    test_overlap    <- length(test_in_pathway)

    if (test_overlap == 0) return(NULL)

    p_value <- phyper(test_overlap - 1, bg_overlap, N - bg_overlap, n, lower.tail = FALSE)

    data.frame(
      Pathway_ID    = pathway,
      bg_overlap    = bg_overlap,
      test_overlap  = test_overlap,
      p_value       = p_value,
      stringsAsFactors = FALSE
    )
  })

  results_list <- Filter(Negate(is.null), results_list)

  # Empty result with full structure
  if (length(results_list) == 0) {
    return(data.frame(
      Pathway_ID       = character(0),
      bg_overlap       = integer(0),
      test_overlap     = integer(0),
      p_value          = numeric(0),
      p_adj            = numeric(0),
      Enrichment_Ratio = numeric(0),
      significance     = character(0),
      Description      = character(0),
      stringsAsFactors = FALSE
    ))
  }

  ora_df <- do.call(rbind, results_list)

  #==============================================================
  # Step 5 — Multiple testing, enrichment metrics, significance
  #==============================================================
  ora_df$p_adj <- p.adjust(ora_df$p_value, method = "BH")
  ora_df$Enrichment_Ratio <- (ora_df$test_overlap / n) / (ora_df$bg_overlap / N)
  ora_df$significance <- ifelse(ora_df$p_adj < pvalueCutoff, "significant", "NS")
  ora_df <- ora_df[order(ora_df$p_adj), ]
  rownames(ora_df) <- NULL

  #==============================================================
  # Step 6 — Pathway descriptions, optional hits, final layout
  #==============================================================
  # KEGG pathway descriptions
  response <- tryCatch(
    rawToChar(
      curl::curl_fetch_memory(
        "https://rest.kegg.jp/list/pathway/",
        handle = curl::new_handle()
      )$content
    ),
    error = function(e) NULL
  )
  if (is.null(response)) {
    stop("Failed to retrieve KEGG pathway descriptions from /list/pathway.")
  }
  pathway_df <- read.delim(
    text = response,
    header = FALSE,
    sep = "\t",
    quote = "",
    stringsAsFactors = FALSE
  )
  colnames(pathway_df) <- c("Pathway_ID", "Description")
  pathway_df$Digits <- sub(".*(\\d{5})$", "\\1", pathway_df$Pathway_ID)

  ora_digits <- sub(".*(\\d{5})$", "\\1", ora_df$Pathway_ID)
  ora_df$Description <- pathway_df$Description[match(ora_digits, pathway_df$Digits)]

  # Optional: add test hits per pathway
  if (add_hits) {
    ora_df$Test_Hits <- vapply(ora_df$Pathway_ID, function(pid) {
      hits <- intersect(pathway_list[[pid]], test_set)
      if (length(hits) == 0) NA_character_ else paste(hits, collapse = ";")
    }, character(1))
  }

  # Column ordering
  col_order <- c(
    "Pathway_ID",
    "Description",
    "bg_overlap",
    "test_overlap",
    if ("Test_Hits" %in% colnames(ora_df)) "Test_Hits",
    "p_value",
    "p_adj",
    "Enrichment_Ratio",
    "significance"
  )

  ora_df <- ora_df[, col_order]

  return(ora_df)
}
