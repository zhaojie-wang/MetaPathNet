## Create a fake MetaPathNet-style network (no KEGG access needed)
fake_net <- matrix(
  c("cpd:C00078", "K00453", "k_compound:irreversible",
    "K00453", "cpd:C00328", "k_compound:irreversible",
    "cpd:C00328", "K00486", "k_compound:reversible",
    "K00486", "cpd:C00328", "k_compound:reversible",
    "K00486", "cpd:C00780", "k_compound:irreversible",
    "dr:D00001", "K00453", "k_compound:irreversible"),
  ncol = 3, byrow = TRUE
)
colnames(fake_net) <- c("source", "target", "interaction_type")

## -----------------------------------------------------------
## MPN_findMappedNodes
## -----------------------------------------------------------
test_that("MPN_findMappedNodes finds mapped and unmapped nodes", {
  res <- MPN_findMappedNodes(
    nodes = c("cpd:C00078", "K00453", "FAKE_NODE"),
    network_table = fake_net
  )

  expect_true("cpd:C00078" %in% res$mapped_nodes)
  expect_true("K00453" %in% res$mapped_nodes)
  expect_true("FAKE_NODE" %in% res$unmapped_nodes)
})

test_that("MPN_findMappedNodes normalises Cxxxxx to cpd:Cxxxxx", {
  res <- MPN_findMappedNodes(
    nodes = c("C00078"),
    network_table = fake_net
  )

  expect_true("cpd:C00078" %in% res$mapped_nodes)
})

## -----------------------------------------------------------
## MPN_mergeNetworks
## -----------------------------------------------------------
test_that("MPN_mergeNetworks merges two networks", {
  net_a <- fake_net[1:3, , drop = FALSE]
  net_b <- fake_net[4:6, , drop = FALSE]

  merged <- MPN_mergeNetworks(net_a, net_b)

  expect_true(is.matrix(merged))
  expect_equal(ncol(merged), 3)
  expect_true(nrow(merged) >= max(nrow(net_a), nrow(net_b)))
})

test_that("MPN_mergeNetworks removes duplicate edges", {
  merged <- MPN_mergeNetworks(fake_net, fake_net)

  expect_equal(nrow(merged), nrow(unique(fake_net)))
})

test_that("MPN_mergeNetworks rejects single input", {
  expect_error(MPN_mergeNetworks(fake_net))
})

## -----------------------------------------------------------
## MPN_removeNode
## -----------------------------------------------------------
test_that("MPN_removeNode removes a node and its edges", {
  result <- MPN_removeNode(
    nodes_to_remove = "K00453",
    network_table   = fake_net
  )

  all_nodes <- unique(as.vector(result[, 1:2]))
  expect_false("K00453" %in% all_nodes)
  expect_true(nrow(result) < nrow(fake_net))
})

test_that("MPN_removeNode normalises compound IDs", {
  result <- MPN_removeNode(
    nodes_to_remove = "C00078",
    network_table   = fake_net
  )

  all_nodes <- unique(as.vector(result[, 1:2]))
  expect_false("cpd:C00078" %in% all_nodes)
})

## -----------------------------------------------------------
## MPN_removeDrugs
## -----------------------------------------------------------
test_that("MPN_removeDrugs removes dr: nodes", {
  result <- MPN_removeDrugs(fake_net)

  all_nodes <- unique(as.vector(result[, 1:2]))
  has_drug <- any(grepl("^dr:", all_nodes))
  expect_false(has_drug)
  expect_true(nrow(result) < nrow(fake_net))
})

test_that("MPN_removeDrugs returns unchanged network when no drugs present", {
  no_drug_net <- fake_net[!grepl("dr:", fake_net[, 1]) & !grepl("dr:", fake_net[, 2]), , drop = FALSE]
  result <- MPN_removeDrugs(no_drug_net)

  expect_equal(nrow(result), nrow(no_drug_net))
})

## -----------------------------------------------------------
## MPN_replaceNode
## -----------------------------------------------------------
test_that("MPN_replaceNode replaces a node ID", {
  result <- MPN_replaceNode(
    nodes_to_replace = "K00453",
    replacement_node = "K99999",
    network_table    = fake_net
  )

  all_nodes <- unique(as.vector(result[, 1:2]))
  expect_false("K00453" %in% all_nodes)
  expect_true("K99999" %in% all_nodes)
})

test_that("MPN_replaceNode normalises compound IDs", {
  result <- MPN_replaceNode(
    nodes_to_replace = "C00078",
    replacement_node = "C99999",
    network_table    = fake_net
  )

  all_nodes <- unique(as.vector(result[, 1:2]))
  expect_false("cpd:C00078" %in% all_nodes)
  expect_true("cpd:C99999" %in% all_nodes)
})
