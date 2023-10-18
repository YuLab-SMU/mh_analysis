test_that("can parse example get_comple", {
    seqs <- get_comple("ATTCCGNNGCC")
    expect_equal(seqs, "TAAGGCNNCGG")
})
