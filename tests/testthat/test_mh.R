test_that("can parse example mh", {
    data(BSgenomehpv16)
    SurVirus_dir <- system.file(file.path("extdata", "survirus_result"), package = "mhAnalysis")
    bam_dir <- file.path(SurVirus_dir, "readsx") 
    result_t1 <- read.table(file.path(SurVirus_dir, "results.t1.txt"), sep = " ", header = FALSE, fill = TRUE)
    results <- read.table(file.path(SurVirus_dir, "results.txt"), sep = " ", header = FALSE, fill = TRUE)
    rownames(results) <- paste("ID", results[, 1], sep = "=")
    results <- results[result_t1[, 1], ]
    result_rel <- add_strand(results, bam_dir, result_t1)
    expect_true("strand_host" %in% colnames(result_rel))
    result_rel <- get_real_loc(result_rel, BSgenome_host=BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38, 
        BSgenome_virus=BSgenomehpv16)   
    result_rel <- get_seq(result_rel, len = 10, BSgenome_host=BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38, 
        BSgenome_virus=BSgenomehpv16)
    expect_true("seq_host" %in% colnames(result_rel))
    result_rel <- add_mh_flank(result_rel, len = 5)
    expect_true("mh" %in% colnames(result_rel))
})
