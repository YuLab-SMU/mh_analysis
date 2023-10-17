#' Get microhomology(mh) length
#'
#' @param seq_df a dataframe of three rows, each containing the host sequence,
#' integration sequence, and virus sequence
#' @param unite If TRUE, calculate the lengths of each microsome separately, 
#' otherwise calculate the sum of microsome lengths
#'
#' @importFrom stringr str_count
#'
#' @return a numeric of mh length
#' @export
#'
#' @examples
#' seq_df <- data.frame(host = c("A", "T", "G", "G", "C", "T", "A", "A"),
#'                      inte = c("A", "T", "G", "G", "A", "T", "A", "C"),
#'                      virus= c("C", "A", "T", "T", "A", "T", "A", "C"))
#' seq_df <- t(seq_df) |> as.data.frame()
#' get_mh(seq_df)
get_mh <- function(seq_df, unite = TRUE) {
    seqs <- c(1,2,3,4,5)
    names(seqs) <- c("A", "T", "G", "C", "N")
    for (i in seq_len(ncol(seq_df))) {
        seq_df[, i] <- seqs[seq_df[, i]]
    }
    aa1 <- seq_df[1, ] - seq_df[2, ]
    aa2 <- seq_df[3, ] - seq_df[2, ]
    aa <- abs(aa1) + abs(aa2)
    aa <- as.numeric(aa)
    # 统计两个以上连续0的数量
    bb <- aa[2:length(aa)] + aa[1:(length(aa)-1)]
    bb <- c(1, as.numeric(bb))
    # 将含有0的都改成非零, 最大是16, 全改成9，为一位，不影响后续分析，保持位数不变
    bb[bb >9 ] <- 9

    # 将0的左右两边都是非零的当作一坨
    # 将连续的多个0替换成一个0，就可以得到坨数
    cc <- paste(bb, collapse = "")
    if (unite) {
        cc <- strsplit(cc, split = "[1-9]+") |> unlist() |> nchar()
        mh <- cc[cc > 0] + 1
        mh <- paste(mh, collapse = "_")
        
    } else {
        cc <- gsub("[0]+", "0", cc)
        mh <- str_count(cc, "0") + sum(bb == 0)
    }   
    mh
}


#' add mh value to SurVirus result (result of get_seq())
#'
#' @param result_rel SurVirus result (result of get_seq())
#' @param unite If TRUE, calculate the lengths of each microsome separately, 
#' otherwise calculate the sum of microsome lengths
#' @return dataframe
#' @export
#'
#' @examples
#' data(BSgenomehpv16)
#' SurVirus_dir <- system.file(file.path("extdata", "survirus_result"), 
#'     package = "mhAnalysis")
#' bam_dir <- file.path(SurVirus_dir, "readsx")
#'
#' result_t1 <- read.table(file.path(SurVirus_dir, "results.t1.txt"), sep = " ", header = FALSE, fill = TRUE)
#' results <- read.table(file.path(SurVirus_dir, "results.txt"), sep = " ", header = FALSE, fill = TRUE)
#' rownames(results) <- paste("ID", results[, 1], sep = "=")
#' results <- results[result_t1[, 1], ]
#' result_rel <- add_strand(results, bam_dir, result_t1)
#' result_rel <- get_real_loc(result_rel, BSgenome_host=BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38, 
#'     BSgenome_virus=BSgenomehpv16)
#' result_rel <- get_seq(result_rel, len = 10, BSgenome_host=BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38, 
#'     BSgenome_virus=BSgenomehpv16)
#' result_rel <- add_mh(result_rel)
add_mh <- function(result_rel, unite = TRUE) {
    mh <- rep(0, nrow(result_rel))
    for (i in seq_len(length(mh))) {
        seq_host <- result_rel[i, "seq_host40"]
        seq_hpv <- result_rel[i, "seq_hpv40"]
        seq_integration <- paste0(result_rel[i, "seq_host"], result_rel[i, "seq_hpv"])
        seq_df <- strsplit(c(seq_host, seq_integration, seq_hpv), "")
        seq_df <- do.call(rbind, seq_df) |> as.data.frame()
        
        mh[i] <- get_mh(seq_df, unite = unite)
    }
    result_rel$mh <- mh
    return(result_rel)
}



#' get mh value in flank region of SurVirus result (result of get_seq())
#'
#' @param result_rel SurVirus result (result of get_seq())
#' @param len flanking region size
#' @param unite If TRUE, calculate the lengths of each microsome separately, 
#' otherwise calculate the sum of microsome lengths
#'
#' @return dataframe
#' @export
#'
#' @examples
#' data(BSgenomehpv16)
#' SurVirus_dir <- system.file(file.path("extdata", "survirus_result"), 
#'     package = "mhAnalysis")
#' bam_dir <- file.path(SurVirus_dir, "readsx")
#'
#' result_t1 <- read.table(file.path(SurVirus_dir, "results.t1.txt"), sep = " ", 
#'     header = FALSE, fill = TRUE)
#' results <- read.table(file.path(SurVirus_dir, "results.txt"), sep = " ", 
#'     header = FALSE, fill = TRUE)
#' rownames(results) <- paste("ID", results[, 1], sep = "=")
#' results <- results[result_t1[, 1], ]
#' result_rel <- add_strand(results, bam_dir, result_t1)
#' result_rel <- get_real_loc(result_rel, BSgenome_host=BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38, 
#'     BSgenome_virus=BSgenomehpv16)
#' result_rel <- get_seq(result_rel, len = 10, BSgenome_host=BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38, 
#'     BSgenome_virus=BSgenomehpv16)
#' result_rel <- add_mh_flank(result_rel, len = 5)
add_mh_flank <- function(result_rel, len, unite = TRUE) {
    result_rel$mh <- rep(0, nrow(result_rel))
    for (i in seq_len(nrow(result_rel))) {
        n1 <- n2 <-  len
        seq_host <- result_rel[i, "seq_host40"]
        seq_hpv <- result_rel[i, "seq_hpv40"]
        seq_integration <- paste0(result_rel[i, "seq_host"], result_rel[i, "seq_hpv"])
        seq_df <- strsplit(c(seq_host, seq_integration, seq_hpv), "")
        seq_df <- do.call(rbind, seq_df) |> as.data.frame()
        seqs <- c(1,2,3,4,5)
        names(seqs) <- c("A", "T", "G", "C", "N")
        for (j in seq_len(ncol(seq_df))) {
            seq_df[, j] <- seqs[seq_df[, j]]
        }
        aa1 <- seq_df[1, ] - seq_df[2, ]
        aa2 <- seq_df[3, ] - seq_df[2, ]
        aa <- abs(aa1) + abs(aa2)
        aa <- as.numeric(aa)
        # 统计两个以上连续0的数量
        bb <- aa[2:length(aa)] + aa[1:(length(aa)-1)]
        bb <- c(1, as.numeric(bb))

        # flank size = 5

        cc <- aa[(length(aa) / 2 - n1 + 1): (length(aa)/2 + n2)]
        # 从中心flank往两侧延申，

        while(cc[1] == 0 || cc[length(cc)] == 0) {
            if (cc[1] == 0) n1 <- n1 + 1
            if (cc[length(cc)] == 0) n2 <- n2 + 1
            if (length(aa) / 2 - n1 < 0 || length(aa)/2 + n2 > length(aa)) {
                break
            }
            n11 <- max(1, (length(aa) / 2 - n1 + 1))
            n22 <- min(ncol(seq_df), (length(aa)/2 + n2))
            cc <- aa[(length(aa) / 2 - n1): (length(aa)/2 + n2)]
        }
        n11 <- max(1, (length(aa) / 2 - n1 + 1))
        n22 <- min(ncol(seq_df), (length(aa)/2 + n2))
        seq_df <- seq_df[, n11:n22]
        result_rel$mh[i] <- get_mh(seq_df, unite = unite)
    }
    return(result_rel)
}
