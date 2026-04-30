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
## MPN_netSimilarity
## -----------------------------------------------------------
test_that("MPN_netSimilarity computes node-level Jaccard", {
  net_a <- fake_net[1:3, , drop = FALSE]
  net_b <- fake_net[3:6, , drop = FALSE]

  sim <- MPN_netSimilarity(net_a, net_b, type = "node")

  expect_true(is.matrix(sim))
  expect_equal(nrow(sim), 2)
  expect_equal(ncol(sim), 2)
  expect_equal(unname(diag(sim)), c(1, 1))
  expect_true(sim[1, 2] >= 0 && sim[1, 2] <= 1)
})

test_that("MPN_netSimilarity computes edge-level Jaccard", {
  sim <- MPN_netSimilarity(fake_net, fake_net, type = "edge")

  expect_true(is.matrix(sim))
  expect_equal(nrow(sim), 2)
  expect_equal(ncol(sim), 2)
  expect_equal(sim[1, 2], 1)
})

test_that("MPN_netSimilarity rejects single input", {
  expect_error(MPN_netSimilarity(fake_net))
})
