test_that("can parse example get_rev", {
    seqs <- get_rev("ATTCCGNNGCC")
    expect_equal(seqs, "CCGNNGCCTTA")
})
