## metabolic_matrix ####
#' Build metabolic edges from parsed KEGG pathways
#' @noRd
metabolic_matrix <- function(path_names, list_parsed_paths, organism_code,
                             expand_genes = FALSE) {
  #==============================================================
  # STEP 1 — Build reaction–metabolite edges (network_from_path)
  #==============================================================
  metabolic_table_list <- lapply(path_names, function(path) {
    lines <- list_parsed_paths[[path]]
    if (is.null(lines)) return(NULL)

    ## Select reactions
    global_lines_reactions <- intersect(grep("rn", lines), grep("Name", lines))
    Reactions <- lines[global_lines_reactions]
    Reactions <- substr(Reactions, 11, nchar(Reactions))

    ## Select substrates
    global_lines_substrates <- grep("Substrate Name", lines)
    Substrates <- lines[global_lines_substrates]
    Substrates <- substr(Substrates, 21, nchar(Substrates))

    ## Select products
    global_lines_products <- grep("Product Name", lines)
    Products <- lines[global_lines_products]
    Products <- substr(Products, 19, nchar(Products))

    ## Select type
    global_lines_type <- grep("Type", lines)
    Type <- lines[global_lines_type]
    Type <- substr(Type, 11, nchar(Type))

    ## Build network edges (inline add_reaction_edge + new_edge)
    all_edges <- mapply(
      FUN = function(Reaction, Substrate, Product, Type_one) {
        reactions  <- unlist(strsplit(Reaction, " "))
        substrates <- unlist(strsplit(Substrate, ";"))
        substrates <- unlist(strsplit(substrates, " "))
        products   <- unlist(strsplit(Product, ";"))
        products   <- unlist(strsplit(products, " "))

        # Build edges for each reaction (inline new_edge)
        edges_list <- lapply(reactions, function(reaction) {
          substrate_lines <- cbind(substrates, rep(reaction, length(substrates)))
          product_lines   <- cbind(rep(reaction, length(products)), products)
          edges_reaction  <- rbind(substrate_lines, product_lines)

          if (Type_one == "reversible") {
            reverse <- cbind(edges_reaction[, 2], edges_reaction[, 1])
            edges_reaction <- rbind(edges_reaction, reverse)
          }
          edges_reaction
        })

        unique(do.call(rbind, edges_list))
      },
      Reaction  = Reactions,
      Substrate = Substrates,
      Product   = Products,
      Type_one  = Type,
      SIMPLIFY  = FALSE
    )

    unique(do.call(rbind, all_edges))
  })

  metabolic_table_list <- Filter(Negate(is.null), metabolic_table_list)
  if (!length(metabolic_table_list)) {
    return(NULL)
  }

  metabolic_table <- unique(do.call(rbind, metabolic_table_list))
  metabolic_table <- matrix(metabolic_table, ncol = 2)
  colnames(metabolic_table) <- c("node1", "node2")
  rownames(metabolic_table) <- NULL

  #==============================================================
  # STEP 2 — Link reactions to genes (enzyme → reaction, KO)
  #          (equivalent to enzymeTable/reactionTable/koTable logic)
  #==============================================================
  file_enzyme   <- paste0("https://rest.kegg.jp/link/enzyme/", organism_code)
  response_enzyme <- tryCatch(RCurl::getURL(file_enzyme), error = function(e) NULL)
  enzymeTable   <- convertTable(response_enzyme)

  file_reaction   <- "https://rest.kegg.jp/link/reaction/enzyme"
  response_reaction <- tryCatch(RCurl::getURL(file_reaction), error = function(e) NULL)
  reactionTable   <- convertTable(response_reaction)

  file_ko       <- paste0("https://rest.kegg.jp/link/ko/", organism_code)
  response_ko <- tryCatch(RCurl::getURL(file_ko), error = function(e) NULL)
  koTable       <- convertTable(response_ko)
  if (is.matrix(koTable)) {
    koTable[, 2] <- substr(koTable[, 2], 4, 9)
  }
  metabolic_table_RG <- NULL

  if (is.matrix(enzymeTable) && is.matrix(reactionTable) && is.matrix(koTable)) {
    message("Linking reactions to genes")

    genes   <- enzymeTable[, 1]
    enzymes <- enzymeTable[, 2]
    ko      <- vector(mode = "character", length = length(genes))

    for (i in seq_along(genes)) {
      index <- which(koTable[, 1] == genes[i])
      if (length(index) > 0) {
        ko[i] <- koTable[index[1], 2]
      } else {
        ko[i] <- genes[i]
      }
    }

    ## Link reactions to genes (inline link_reaction_gene)
    reactions_genes <- mapply(
      FUN = function(enzyme, gene, ko_one) {
        index <- which(reactionTable[, 1] == enzyme)
        if (length(index) > 0) {
          reactions <- reactionTable[index, 2]
          line <- cbind(
            reactions,
            rep(enzyme, length(reactions)),
            rep(gene,   length(reactions)),
            rep(ko_one, length(reactions))
          )
          matrix(line, ncol = 4, nrow = length(reactions))
        } else {
          NULL
        }
      },
      enzyme = enzymes,
      gene   = genes,
      ko_one = ko,
      SIMPLIFY = FALSE
    )

    reactions_genes <- Filter(Negate(is.null), reactions_genes)
    if (!length(reactions_genes)) {
      return(NULL)
    }

    reactionM   <- unique(do.call(rbind, reactions_genes))
    reactionko  <- unique(na.omit(reactionM[, c(1, 4)]))
    reactiongene <- unique(na.omit(reactionM[, c(1, 3)]))

    #============================================================
    # STEP 3 — Replace reactions by KO / gene in edges
    #          (inline bind_reaction_gene + build_RG_edge)
    #============================================================
    all_edges_RG <- lapply(
      X = split(metabolic_table, row(metabolic_table)),
      FUN = function(edge_row) {
        edge <- as.character(edge_row)
        index_rn <- grep("rn", edge)
        if (length(index_rn) != 1) {
          return(NULL)
        }

        indexM <- if (index_rn == 1) 2 else 1
        reaction   <- edge[index_rn]
        metabolite <- edge[indexM]

        # bind_reaction_gene logic
        if (isTRUE(expand_genes)) {
          mat <- reactiongene
        } else {
          mat <- reactionko
        }

        index_reaction <- which(mat[, 1] == reaction)
        if (length(index_reaction) > 0) {
          reaction_genes_all <- paste(reaction, mat[index_reaction, 2], sep = "_")
        } else {
          reaction_genes_all <- reaction
        }

        # build_RG_edge logic
        if (index_rn == 1) {
          new_edges <- cbind(
            reaction_genes_all,
            rep(metabolite, length(reaction_genes_all))
          )
        } else {
          new_edges <- cbind(
            rep(metabolite, length(reaction_genes_all)),
            reaction_genes_all
          )
        }
        new_edges
      }
    )

    all_edges_RG <- Filter(Negate(is.null), all_edges_RG)
    if (!length(all_edges_RG)) {
      return(NULL)
    }

    metabolic_table_RG <- unique(do.call(rbind, all_edges_RG))
    metabolic_table_RG <- matrix(metabolic_table_RG, ncol = 2)
    colnames(metabolic_table_RG) <- c("node1", "node2")
    rownames(metabolic_table_RG) <- NULL
  }

  return(metabolic_table_RG)
}


## signaling_matrix ####
#' Build signaling edges from a KEGG graph object
#' @noRd
signaling_matrix <- function(global_network_all) {

  signaling_table <- NULL

  #==============================================================
  # STEP 1 — Create raw network edges from all nodes
  #==============================================================

  nodes <- graph::nodes(global_network_all)
  edges <- KEGGgraph::edges(global_network_all)

  ## For a given node, it retrieves all connected nodes from the edges list
  network_edges <- lapply(nodes, function(node) {
    linked_nodes <- edges[[node]]

    if (length(linked_nodes) >= 1) {
      # Node with edges
      new_edge <- cbind(rep(node, length(linked_nodes)), edges[[node]])
    } else {
      new_edge <- NULL
    }
    return(new_edge)
  })

  network_edges <- unique(do.call(rbind, network_edges))

  #==============================================================
  # STEP 2 — Filter invalid or unwanted nodes
  #==============================================================

  if (length(network_edges) == 0) {
    to_print <- "Impossible to build a signaling network"
    warning(to_print)
  } else {
    network_edges <- matrix(network_edges, ncol = 2)

    ## Remove the edges that contain paths
    all_lines_paths <- c(
      grep("path", network_edges[, 1]),
      grep("path", network_edges[, 2])
    )

    if (length(all_lines_paths) == nrow(network_edges)) {
      # All edges contain paths.
      to_print <- "Impossible to build a signaling network"
      warning(to_print)
    } else {

      if (length(all_lines_paths) >= 1) {
        network_edges <- network_edges[-(all_lines_paths), ]
        network_edges <- matrix(network_edges, ncol = 2)
      }

      ## Remove non-biological nodes from KEGG
      nonbio_nodes <- paste(
        "fibrates", "non-steroids", "agents",
        "Antiinflammatory", "BR:br08303",
        "Thiazolidinediones", sep = "|"
      )
      nonbio_edges <- c(
        grep(nonbio_nodes, network_edges[, 1]),
        grep(nonbio_nodes, network_edges[, 2])
      )
      nonbio_edges <- unique(nonbio_edges)

      if (length(nonbio_edges) >= 1) {
        network_edges <- network_edges[-c(nonbio_edges), ]
      }

      signaling_table <- unique(network_edges)
      signaling_table <- matrix(signaling_table, ncol = 2)
      colnames(signaling_table) <- c("node1", "node2")
    }
  }

  return(signaling_table)
}


## map_interactions ####
#' Map interaction types onto network edges
#' @noRd
map_interactions <- function(MetaPathNet_table, interaction_type) {
  #---------------------------------------------------------------
  # Use a unique, regex-safe delimiter that won’t appear in KEGG IDs
  #---------------------------------------------------------------
  delimiter <- "##"

  # Build encoded edge representations
  edge_net     <- paste(MetaPathNet_table[, 1], MetaPathNet_table[, 2], sep = delimiter)
  edge_net_rev <- paste(MetaPathNet_table[, 2], MetaPathNet_table[, 1], sep = delimiter)

  # Detect compound / glycan / drug edges
  cpd_ind <- grep("cpd:|gl:|dr:", edge_net)

  #---------------------------------------------------------------
  # 1. Mixed network (both gene and compound edges)
  #---------------------------------------------------------------
  if (length(cpd_ind) > 0 && (length(cpd_ind) < length(edge_net))) {
    edge_gene <- edge_net[-cpd_ind]

    # Annotate signaling interactions for gene–gene edges
    signaling_interactions <- vapply(edge_gene, function(st) {
      collapsed_edges <- paste(interaction_type[, 1], interaction_type[, 2], sep = delimiter)
      ind <- which(collapsed_edges == st)

      if (length(ind) == 0) {
        "unknown"
      } else {
        subtypes <- unique(gsub(";", "/", interaction_type[ind, 3]))
        subtypes[1]
      }
    }, character(1))

    # Handle compound edges (reversible / irreversible)
    edge_compound     <- edge_net[grep("cpd:|gl:|dr:", edge_net)]
    edge_compound_rev <- edge_net_rev[grep("cpd:|gl:|dr:", edge_net)]
    reversible_edges  <- intersect(edge_compound, edge_compound_rev)

    if (length(reversible_edges) > 0) {
      irreversible_edges <- setdiff(edge_compound, reversible_edges)
      rev_edges <- do.call(rbind, strsplit(reversible_edges, delimiter, fixed = TRUE))
      edges_direc <- cbind(rev_edges, "compound:reversible")

      if (length(irreversible_edges) > 0) {
        irev_edges <- do.call(rbind, strsplit(irreversible_edges, delimiter, fixed = TRUE))
        irev_edges <- cbind(irev_edges, "compound:irreversible")
        edges_direc <- rbind(edges_direc, irev_edges)
      }
    } else {
      irev_edges <- do.call(rbind, strsplit(edge_compound, delimiter, fixed = TRUE))
      edges_direc <- cbind(irev_edges, "compound:irreversible")
    }

    colnames(edges_direc) <- c("source", "target", "interaction_type")
    metabo_net <- edges_direc

    # Rebuild gene–gene table and merge both parts
    gene_edges_split <- do.call(rbind, strsplit(edge_gene, delimiter, fixed = TRUE))
    gene_net <- cbind(gene_edges_split, signaling_interactions)
    colnames(gene_net) <- c("source", "target", "interaction_type")

    merged_net <- rbind(metabo_net, gene_net)

    #---------------------------------------------------------------
    # 2. Only compound edges
    #---------------------------------------------------------------
  } else if (length(cpd_ind) == length(edge_net)) {
    edge_compound     <- edge_net
    edge_compound_rev <- edge_net_rev
    reversible_edges  <- intersect(edge_compound, edge_compound_rev)

    if (length(reversible_edges) > 0) {
      irreversible_edges <- setdiff(edge_compound, reversible_edges)
      rev_edges <- do.call(rbind, strsplit(reversible_edges, delimiter, fixed = TRUE))
      edges_direc <- cbind(rev_edges, "compound:reversible")

      if (length(irreversible_edges) > 0) {
        irev_edges <- do.call(rbind, strsplit(irreversible_edges, delimiter, fixed = TRUE))
        irev_edges <- cbind(irev_edges, "compound:irreversible")
        edges_direc <- rbind(edges_direc, irev_edges)
      }
    } else {
      irev_edges <- do.call(rbind, strsplit(edge_compound, delimiter, fixed = TRUE))
      edges_direc <- cbind(irev_edges, "compound:irreversible")
    }

    colnames(edges_direc) <- c("source", "target", "interaction_type")
    merged_net <- edges_direc

    #---------------------------------------------------------------
    # 3. Only gene–gene edges
    #---------------------------------------------------------------
  } else {
    edge_gene <- edge_net

    signaling_interactions <- vapply(edge_gene, function(st) {
      collapsed_edges <- paste(interaction_type[, 1], interaction_type[, 2], sep = delimiter)
      ind <- which(collapsed_edges == st)

      if (length(ind) == 0) {
        "unknown"
      } else {
        subtypes <- unique(gsub(";", "/", interaction_type[ind, 3]))
        subtypes[1]
      }
    }, character(1))

    gene_edges_split <- do.call(rbind, strsplit(edge_gene, delimiter, fixed = TRUE))
    gene_net <- cbind(gene_edges_split, signaling_interactions)
    colnames(gene_net) <- c("source", "target", "interaction_type")
    merged_net <- gene_net
  }

  #---------------------------------------------------------------
  # Final cleanup and conversion to matrix
  #---------------------------------------------------------------
  merged_net <- as.matrix(merged_net)
  merged_net <- apply(merged_net, 2, as.character)
  rownames(merged_net) <- NULL
  merged_net[, "interaction_type"] <- paste0("k_", merged_net[, "interaction_type"])

  return(merged_net)
}


## get_metabonet ####
#' Retrieve and parse a KEGG metabolic pathway
#' @noRd
get_metabonet <- function(path, all_paths, organism_code) {
  message(path)

  if (!(path %in% all_paths)) {
    to_print <- paste(path, "-incorrect path ID:path removed", sep = "")
    message(to_print)
    return(NULL)
  }

  if (substr(path, nchar(organism_code) + 1, nchar(path)) == "01100") {
    ## Remove metabolic overview map
    return(NULL)
  }

  file <- paste("https://rest.kegg.jp/get/", path, "/kgml", sep = "")

  pathway <- tryCatch(
    RCurl::getURL(file),
    error = function(e) NULL
  )

  if (is.null(pathway)) {
    to_print <- paste(path, "-path ID without XML:path removed", sep = "")
    message(to_print)
    return(NULL)
  }

  parsed_kgml <- tryCatch(
    KEGGgraph::parseKGML(pathway),
    error = function(e) NULL
  )

  if (is.null(parsed_kgml)) {
    to_print <- paste(path, "-path ID without XML:path removed", sep = "")
    message(to_print)
    return(NULL)
  }

  reactions <- tryCatch(
    KEGGgraph::getReactions(parsed_kgml),
    error = function(e) NULL
  )

  if (is.null(reactions)) {
    to_print <- paste(path, "-path ID without XML:path removed", sep = "")
    message(to_print)
    return(NULL)
  }

  parsed_path <- capture.output(reactions, file = NULL)
  return(parsed_path)
}


## get_signalnet ####
#' Retrieve and parse a KEGG signaling pathway
#' @noRd
get_signalnet <- function(path, all_paths) {
  message(path)

  if (!(path %in% all_paths)) {
    to_print <- paste(path, "-incorrect path ID:path removed", sep = "")
    message(to_print)
    return("bad_path")
  }

  file <- paste("https://rest.kegg.jp/get/", path, "/kgml", sep = "")

  pathway <- tryCatch(
    RCurl::getURL(file),
    error = function(e) NULL
  )

  if (is.null(pathway)) {
    to_print <- paste(path, "-path ID without XML:path removed", sep = "")
    message(to_print)
    return("bad_path")
  }

  path_parsed <- tryCatch(
    KEGGgraph::parseKGML(pathway),
    error = function(e) NULL
  )

  if (is.null(path_parsed)) {
    to_print <- paste(path, "-path ID without XML:path removed", sep = "")
    message(to_print)
    return("bad_path")
  }

  path_network <- tryCatch(
    KEGGgraph::KEGGpathway2Graph(
      path_parsed,
      genesOnly = FALSE,
      expandGenes = TRUE
    ),
    error = function(e) NULL
  )

  if (is.null(path_network)) {
    to_print <- paste(path, "-path ID without XML:path removed", sep = "")
    message(to_print)
    return("bad_path")
  }

  urlcheck <- tryCatch(
    KEGGgraph::edges(path_network),
    error = function(e) NULL
  )

  if (is.null(urlcheck)) {
    to_print <- paste(path, "-path ID without XML:path removed", sep = "")
    message(to_print)
    return("bad_path")
  }

  return(path_network)
}


## convertTable ####
#' Convert tab-delimited KEGG responses to a matrix
#' @noRd
convertTable <- function(res) {
  if (is.null(res) || length(res) == 0 || nchar(res) == 0) {
    message("no result")
    result <- NULL
  } else {
    rows <- strsplit(res, "\n")
    rows.len <- length(rows[[1]])
    result <- matrix(
      unlist(lapply(rows, strsplit, "\t")),
      nrow = rows.len,
      byrow = TRUE
    )
  }
  return(result)
}


## network_features ####
#' Print basic network features
#' @noRd
network_features <- function(network_table) {
  all_nodes <- unique(as.vector(network_table))
  nodes_number <- length(all_nodes)
  message("Network features:")
  to_print <- paste("Number of nodes:", nodes_number, sep = "")
  message(to_print)
  to_print <- paste("Number of edges:", nrow(network_table), sep = "")
  message(to_print)
}

