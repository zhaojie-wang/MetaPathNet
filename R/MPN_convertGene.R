#' Convert Ensembl or Entrez Gene IDs to KEGG Gene IDs and KO Terms
#'
#' Converts Ensembl or Entrez gene identifiers to KEGG gene IDs and,
#' optionally, to KEGG Orthology (KO) IDs using \pkg{mygene} and
#' \pkg{KEGGREST}. The input gene type must be specified.
#'
#' This function supports integration of commonly used gene identifier systems
#' (Ensembl/Entrez) into KO-based and KEGG network workflows.
#'
#' The returned data frame keeps all input IDs, including those that cannot
#' be mapped; unmapped entries are labelled as \code{"unmapped"} in the
#' corresponding columns. This is useful when building KO- or pathway-centric
#' views of host genes, for example in tryptophan–kynurenine metabolism or
#' other immune–metabolic axes.
#'
#' @param genes A character vector of gene IDs. All IDs must be of the same type,
#'   specified by \code{gene_type}.
#' @param gene_type Type of gene ID supplied. One of \code{"ensembl"} or
#'   \code{"entrez"}.
#' @param orthology Logical; if \code{TRUE} (default), KEGG gene IDs are further
#'   mapped to KO IDs.
#'
#' @return A \code{data.frame} with columns:
#' \itemize{
#'   \item \code{Input_ID}: original gene ID (Ensembl or Entrez).
#'   \item \code{Entrez_ID}: mapped Entrez ID, or \code{"unmapped"}.
#'   \item \code{KEGG_ID}: mapped KEGG gene ID, or \code{"unmapped"}.
#'   \item \code{KO_ID}: mapped KEGG Orthology ID, or \code{"unmapped"}
#'         (present only when \code{orthology = TRUE}).
#' }
#'
#' @details
#' \itemize{
#'   \item Only one gene type can be used per call; mixed Ensembl/Entrez input
#'         must be split beforehand.
#'   \item When \code{gene_type = "ensembl"}, Ensembl IDs are first converted to
#'         Entrez IDs using \pkg{mygene}, then mapped to KEGG gene IDs via
#'         \pkg{KEGGREST}.
#'   \item When \code{gene_type = "entrez"}, input IDs are checked for a numeric
#'         pattern; non-numeric entries are retained and marked as
#'         \code{"unmapped"}.
#'   \item KEGG gene IDs are optionally mapped to KO terms using
#'         \code{\link[KEGGREST]{keggLink}} when \code{orthology = TRUE}.
#'   \item If any step fails (no Entrez, no KEGG gene, or no KO), the
#'         corresponding field is set to \code{"unmapped"}; no rows are dropped.
#' }
#'
#' @examples
#' ## Example 1 — Ensembl IDs for human tryptophan–kynurenine genes
#' trp_kynurenine_ens <- c("ENSG00000151790", "ENSG00000131203")
#'
#' MPN_convertGene(
#'   genes      = trp_kynurenine_ens,
#'   gene_type  = "ensembl",
#'   orthology  = TRUE
#' )
#'
#' ## Example 2 — Entrez IDs linked to downstream MetaPathNet network use
#' ## Human enzymes TDO2 (6999) and IDO1 (3620), plus tnaA (948221; a typical microbial enzyme),
#' ## are included here to illustrate how converted KO terms can be checked in a human network.
#' trp_kynurenine_genes <- c("6999", "3620", "948221")
#'
#' trp_kyn_conversion <- MPN_convertGene(
#'   genes      = trp_kynurenine_genes,
#'   gene_type  = "entrez",
#'   orthology  = TRUE
#' )
#'
#' @seealso
#' \code{\link{MPN_keggFinder}}
#'
#' @importFrom KEGGREST keggConv keggLink
#' @importFrom mygene query
#' @export
MPN_convertGene <- function(genes,
                           gene_type  = c("ensembl", "entrez"),
                           orthology = TRUE) {

  if (!is.logical(orthology) || length(orthology) != 1L || is.na(orthology)) {
    stop("orthology must be a single logical value (TRUE/FALSE).")
  }
  gene_type <- match.arg(gene_type)

  # Basic cleaning and de-duplication of input IDs
  genes <- gsub(" ", "", as.character(genes))
  genes[is.na(genes)] <- ""
  genes <- genes[genes != ""]
  genes <- unique(genes)

  if (length(genes) == 0) {
    stop("No valid gene IDs provided after removing NA/empty values.")
  }

  #==============================================================
  # Step 1 — Map input IDs to Entrez IDs
  #==============================================================

  if (gene_type == "ensembl") {
    # Ensembl IDs should start with "ENS"
    if (!all(grepl("^ENS", genes))) {
      message("Some input Ensembl IDs may be invalid (do not start with 'ENS').")
    }

    # Convert Ensembl -> Entrez using mygene
    entrez_ids <- vapply(genes, function(id) {
      res <- tryCatch(
        mygene::query(id, fields = "entrezgene"),
        error = function(e) NULL
      )
      if (!is.null(res) && !is.null(res$hits$entrezgene)) {
        as.character(res$hits$entrezgene[1])
      } else {
        "unmapped"
      }
    }, character(1), USE.NAMES = FALSE)

  } else {  # gene_type == "entrez"

    # We keep ALL original Entrez-like input IDs and mark invalid ones as "unmapped"
    numeric_idx <- grepl("^\\d+$", genes)

    if (any(!numeric_idx)) {
      message("Some input IDs are not valid numeric Entrez IDs and were marked as 'unmapped'.")
    }

    entrez_ids <- ifelse(numeric_idx, genes, "unmapped")
  }

  #==============================================================
  # Step 2 — Map Entrez IDs to KEGG gene IDs
  #==============================================================

  kegg_gene_ids <- vapply(entrez_ids, function(eid) {
    if (eid == "unmapped" || is.na(eid)) {
      return("unmapped")
    }
    ncbi_gene <- paste0("ncbi-geneid:", eid)

    res <- tryCatch(
      KEGGREST::keggConv("genes", ncbi_gene),
      error = function(e) NULL
    )

    if (is.null(res) || length(res) == 0) {
      "unmapped"
    } else {
      as.character(res[[1]])
    }
  }, character(1), USE.NAMES = FALSE)

  #==============================================================
  # Step 3 — Map KEGG gene IDs to KO IDs (optional)
  #==============================================================

  ko_ids <- if (orthology) {
    vapply(kegg_gene_ids, function(kg) {
      if (kg == "unmapped" || is.na(kg)) return("unmapped")

      res <- tryCatch(
        KEGGREST::keggLink("ko", kg),
        error = function(e) NULL
      )

      if (is.null(res) || length(res) == 0) {
        "unmapped"
      } else {
        sub("ko:", "", res[[1]])
      }
    }, character(1), USE.NAMES = FALSE)
  } else {
    NULL
  }

  #==============================================================
  # Step 4 — Assemble final result
  #==============================================================

  res_df <- data.frame(
    Input_ID  = genes,
    Entrez_ID = entrez_ids,
    KEGG_ID   = kegg_gene_ids,
    stringsAsFactors = FALSE
  )

  if (orthology) {
    res_df$KO_ID <- ko_ids
  }

  #==============================================================
  # Step 5 — Mapping summary (transparency)
  #==============================================================
  final_col <- if (orthology) "KO_ID" else "KEGG_ID"

  message(
    "Mapping summary | Input: ", nrow(res_df),
    " | Mapped: ", sum(res_df[[final_col]] != "unmapped"),
    " | Unmapped: ", sum(res_df[[final_col]] == "unmapped")
  )

  return(res_df)
}
