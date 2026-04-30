## Create a fake MetaPathNet-style network-independent reaction table
## (no KEGG access needed)

## -----------------------------------------------------------
## MPN_customReaction
## -----------------------------------------------------------
test_that("MPN_customReaction builds edges with KO layer", {
  rt <- data.frame(
    substrates = "cpd:C00078",
    products   = "cpd:C00328",
    ko         = "K00453",
    direction  = "irreversible",
    stringsAsFactors = FALSE
  )

  result <- MPN_customReaction(
    reaction_table = rt,
    substrate_col  = "substrates",
    product_col    = "products",
    ko_col         = "ko",
    direction_col  = "direction"
  )

  expect_true(is.matrix(result))
  expect_equal(ncol(result), 3)
  expect_true(nrow(result) >= 2)
  expect_true(any(result[, "interaction_type"] == "custom:irreversible"))
})

test_that("MPN_customReaction builds reversible edges in both directions", {
  rt <- data.frame(
    substrates = "cpd:C00078",
    products   = "cpd:C00328",
    ko         = "K00453",
    direction  = "reversible",
    stringsAsFactors = FALSE
  )

  result <- MPN_customReaction(
    reaction_table = rt,
    substrate_col  = "substrates",
    product_col    = "products",
    ko_col         = "ko",
    direction_col  = "direction"
  )

  forward <- any(result[, "source"] == "cpd:C00078" & result[, "target"] == "K00453")
  reverse <- any(result[, "source"] == "K00453" & result[, "target"] == "cpd:C00078")
  expect_true(forward)
  expect_true(reverse)
})

test_that("MPN_customReaction works without KO column", {
  rt <- data.frame(
    substrates = "cpd:C00078",
    products   = "cpd:C00328",
    stringsAsFactors = FALSE
  )

  result <- MPN_customReaction(
    reaction_table = rt,
    substrate_col  = "substrates",
    product_col    = "products"
  )

  expect_true(is.matrix(result))
  expect_equal(ncol(result), 3)
  expect_true(any(result[, "source"] == "cpd:C00078" & result[, "target"] == "cpd:C00328"))
})
