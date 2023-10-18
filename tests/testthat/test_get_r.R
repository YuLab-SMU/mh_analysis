test_that("can parse example get_r", {
    reads <- get_r(c(68,74))
    expect_equal(reads, c("read1", "read1"))
})
