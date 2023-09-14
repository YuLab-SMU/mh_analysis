#' get_r
#'
#' @param y flag of bam/sam file, a vector.
#'
#' @return a vector of character, shows the reads are "read1" or "read2"
#' @export
#'
#' @examples
#' get_r(c(68,74))
get_r <- function(y) {
    result <- rep("read2", length(y))
    for (i in 1:length(y)) {
        x <- y[i]
        aa <- as.numeric(intToBits(x))
        cc <- rep(0, length(aa))
        names(cc) <- 2^(0:(length(aa)-1))
        r <- names(cc)[aa == 1]
        if ("64" %in% r) result[i] <- "read1"
    }
    result
}


#' Get complementary sequences
#'
#' @param x a sequence contains A T G C N
#'
#' @return a complementary sequence contains A T G C N
#' @export
#'
#' @examples
#' get_comple("ATTCCGNNGCC")
get_comple <- function(x) {
    seqs <- c("A", "T", "G", "C", "N")
    names(seqs) <- c("T", "A", "C", "G", "N")
    y <- strsplit(x, "") |> unlist()
    y2 <- seqs[y]
    return(paste(y2, collapse = ""))

}


#' Get reverse sequences
#'
#' @param x a sequence contains A T G C N
#' @return a reverse sequence
#' @export
#'
#' @examples
#' get_rev("ATTCCGNNGCC")
get_rev <- function(x) {
    y <- strsplit(x, "") |> unlist()
    return(paste(rev(y), collapse = ""))
}




#' Add strand information
#'
#' @param results dataframe of results.txt in SurVirus software
#' @param bam_dir The directory of readsx folder (contains bam files) of
#' @param result_t1 result_t1.txt of SurVirus software
#' @importFrom Rsamtools ScanBamParam
#' @importFrom Rsamtools scanBamFlag
#' @importFrom GenomicAlignments readGAlignments
#' @importFrom stats na.omit
#' @importFrom magrittr `%>%`
#' @return datafrmae contaions strand information
#' @export
#'
#' @examples
#' \dontrun{
#' bidui_dir <- "SurVirus_result"
#' bam_dir <- "SurVirus_result/readsx"
#'
#' result_t1 <- fread(file.path(bidui_dir, "results.t1.txt"), sep = " ", header = FALSE, fill = TRUE)
#' results <- fread(file.path(bidui_dir, "results.txt"), sep = " ", header = FALSE, fill = TRUE)
#' class(results) <- class(result_t1) <- "data.frame"
#' rownames(results) <- paste("ID", results[, 1], sep = "=")
#' results <- results[result_t1[, 1], ]
#' result_rel <- add_strand(results, result_t1, result)
#' result_rel <- get_real_loc(result_rel)
#' }
add_strand <- function(results, bam_dir, result_t1) {
    `.` <- NULL
    result_rel <- data.frame(loc_host = result_t1[, 2],
                             loc_hpv = result_t1[, 3],
                             range_host = results[, 2],
                             range_hpv = results[, 3],
                             strand_host = rep(NA, nrow(results)),
                             strand_hpv = rep(NA, nrow(results)),
                             seqid = rep(NA, nrow(results)))
    rownames(result_rel) <- results[, 1]
    param <- ScanBamParam(flag=scanBamFlag(isDuplicate=FALSE,
                                       isSecondaryAlignment=FALSE),
                      what=c("qual", "flag"))
    for (i in seq_len(nrow(results))) {
        filename <- file.path(bam_dir, paste0(results[i, 1], ".bam"))
        bamFile <- readGAlignments(filename, param=param, use.names = TRUE) |> as.data.frame()
        bamFile$ID <- gsub("\\..*", "", rownames(bamFile))

        aa <- bamFile
        aa$reads <- get_r(aa$flag)
        # aa$strand <- get_zf(aa$flag)
        tt <- table(aa$seqnames) |> as.data.frame()
        tt <- tt[order(tt$Freq, decreasing = TRUE), ]
        chrs <- tt[1:2, 1] |> as.character()
        host_chr <- chrs[grep("chr", chrs)]
        hpv_chr <- chrs[grep("NC", chrs)]

        aa1 <- aa[aa[, "seqnames"] == host_chr, ]
        aa2 <- aa[aa[, "seqnames"] == hpv_chr, ]

        aa1$id <- paste(aa1$ID, aa1$reads, sep = "_")
        aa2$id <- paste(aa2$ID, aa2$reads, sep = "_")
        ids <- intersect(aa1$id, aa2$id)

        if (length(ids) == 0) next

        aa1 <- aa1[aa1$id %in% ids, ]
        aa2 <- aa2[aa2$id %in% ids, ]
        aa1 <- aa1[order(aa1$id), ]
        aa2 <- aa2[order(aa2$id), ]


        t_loc <- as.data.frame(table(aa1$start))
        t_loc <- t_loc[order(t_loc[, 2], decreasing = TRUE), ]
        aa1 <- aa1[aa1$start == t_loc[1,1], ]
        aa2 <- aa2[aa2$id %in% aa1$id, ]
        t_loc <- as.data.frame(table(aa2$start))
        t_loc <- t_loc[order(t_loc[, 2], decreasing = TRUE), ]
        aa2 <- aa2[aa2$start == t_loc[1,1], ]
        ids <- intersect(aa1$id, aa2$id)
        if (length(ids) == 0) next
        aa1 <- aa1[aa1$id %in% ids, ]
        aa2 <- aa2[aa2$id %in% ids, ]
        aa1 <- aa1[order(aa1$id), ]
        aa2 <- aa2[order(aa2$id), ]
        result_rel[i, "seqid"] <- aa1$id[1]
        aa1$strand <- as.character(aa1$strand)
        aa2$strand <- as.character(aa2$strand)
        result_rel[i, "strand_host"] <- aa1$strand[1]
        result_rel[i, "strand_hpv"] <- aa2$strand[1]
    }
    result_rel <- na.omit(result_rel)
    rownames(result_rel) <- paste("ID", rownames(result_rel), sep = "=")
    result_rel[, 1] <- gsub("[+-]", "", result_rel[, 1])
    result_rel[, 2] <- gsub("[+-]", "", result_rel[, 2])
    result_rel$start_host <- strsplit(result_rel$range_host, ":") %>% do.call(rbind, .) %>% .[, 3] %>% as.integer()
    result_rel$end_host <- strsplit(result_rel$range_host, ":") %>% do.call(rbind, .) %>% .[, 4] %>% as.integer()
    result_rel$start_hpv <- strsplit(result_rel$range_hpv, ":") %>% do.call(rbind, .) %>% .[, 3] %>% as.integer()
    result_rel$end_hpv <- strsplit(result_rel$range_hpv, ":") %>% do.call(rbind, .) %>% .[, 4] %>% as.integer()
    result_rel$chr <- strsplit(result_rel$loc_host, ":") %>% do.call(rbind, .) %>% .[, 1]
    result_rel$loc_host <- strsplit(result_rel$loc_host, ":") %>% do.call(rbind, .) %>% .[, 2] %>% as.integer()
    result_rel$loc_hpv <- strsplit(result_rel$loc_hpv, ":") %>% do.call(rbind, .) %>% .[, 2] %>% as.integer()
    return(result_rel)
}



#' Get the true insertion site of the virus
#'
#' @param result_rel dataframe of results.txt in SurVirus software,
#' which contains strand information (result of add_strand())
#'
#' @param BSgenome_host BSgenome object of host
#' @param BSgenome_virus BSgenome object of virus
#' @return data.frame
#' @export
#'
#' @examples
#' \dontrun{
#' bidui_dir <- "SurVirus_result"
#' bam_dir <- "SurVirus_result/readsx"
#'
#' result_t1 <- fread(file.path(bidui_dir, "results.t1.txt"), sep = " ", header = FALSE, fill = TRUE)
#' results <- fread(file.path(bidui_dir, "results.txt"), sep = " ", header = FALSE, fill = TRUE)
#' class(results) <- class(result_t1) <- "data.frame"
#' rownames(results) <- paste("ID", results[, 1], sep = "=")
#' results <- results[result_t1[, 1], ]
#' result_rel <- add_strand(results)
#' result_rel <- get_real_loc(result_rel)
#' }
get_real_loc <- function(result_rel, BSgenome_host,BSgenome_virus) {
    hostloc_start <- which(result_rel$start_host == result_rel$loc_host)
    hostloc_end <- which(result_rel$end_host == result_rel$loc_host)
    result_rel[hostloc_start, "end_host"] <- result_rel[hostloc_start, "start_host"] + 19
    result_rel[hostloc_end, "start_host"] <- result_rel[hostloc_end, "end_host"] - 19

    hpvloc_start <- which(result_rel$start_hpv == result_rel$loc_hpv)
    hpvloc_end <- which(result_rel$end_hpv == result_rel$loc_hpv)
    result_rel[hpvloc_start, "end_hpv"] <- result_rel[hpvloc_start, "start_hpv"] + 19
    result_rel[hpvloc_end, "start_hpv"] <- result_rel[hpvloc_end, "end_hpv"] - 19

    result_rel$seq_host <- result_rel$seq_hpv <- rep("1", nrow(result_rel))



    for (i in 1:nrow(result_rel)) {
        
        seq_host <- getSeq(BSgenome_host, result_rel$chr[i], start = result_rel$start_host[i], end = result_rel$end_host[i]) |> as.character()
        seq_hpv <- getSeq(BSgenome_virus, "NC_001526.2", start = result_rel$start_hpv[i], end = result_rel$end_hpv[i]) |> as.character()
        if (result_rel$strand_host[i] == "-") seq_host <- get_comple(seq_host)
        if (result_rel$strand_hpv[i] == "-") seq_hpv <- get_comple(seq_hpv)
        result_rel$seq_host[i] <- seq_host
        result_rel$seq_hpv[i] <- seq_hpv
    }


    for (i in hostloc_start) {
        result_rel$seq_host[i] <- get_rev(result_rel$seq_host[i])
    }


    # Loop search, host from back to front, hpv from front to back, find all sequences that can match, and pick the longest one
    ins_mh <- rep(0, nrow(result_rel))
    get_ins_mh_len <- function(seq_host, seq_hpv) {
        for (j in 20:1) {
            host <- substr(seq_host, (21-j),20)
            hpv <- substr(seq_hpv, 1,j)
            if (identical(host, hpv)) {
                return(j)
            }
        }
        return(0)
    }
    for (i in seq_len(nrow(result_rel))) {
        ins_mh[i] <- get_ins_mh_len(result_rel[i, "seq_host"], result_rel[i, "seq_hpv"])
    }
    result_rel$ins_mh <- ins_mh

  
    result_rel[hpvloc_start, "start_hpv"] <- result_rel[hpvloc_start, "start_hpv"] + result_rel[hpvloc_start, "ins_mh"]
    result_rel[hpvloc_start, "end_hpv"] <- result_rel[hpvloc_start, "start_hpv"] + 19
    result_rel[hpvloc_end, "end_hpv"] <- result_rel[hpvloc_end, "end_hpv"] - result_rel[hpvloc_end, "ins_mh"]
    result_rel[hpvloc_end, "start_hpv"] <- result_rel[hpvloc_end, "end_hpv"] - 19

    result_rel[hpvloc_start, "loc_hpv"] <- result_rel[hpvloc_start, "start_hpv"]
    result_rel[hpvloc_end, "loc_hpv"] <- result_rel[hpvloc_end, "end_hpv"]
    result_rel$seq_hpv <- result_rel$seq_host <- result_rel$range_host <- result_rel$range_hpv <-NULL
    return(result_rel)
}



#' Get sequence of SurVirus result, and update the location(start and end) of
#' host sequence and virus sequence
#'
#' @param result_rel result of get_real_loc()
#' @param len flank region size
#' @param BSgenome_host BSgenome object of host
#' @param BSgenome_virus BSgenome object of virus
#' @importFrom BSgenome getSeq
#'
#' @return dataframe
#' @export
#'
#' @examples
#' \dontrun{
#' bidui_dir <- "SurVirus_result"
#' bam_dir <- "SurVirus_result/readsx"
#'
#' result_t1 <- fread(file.path(bidui_dir, "results.t1.txt"), sep = " ", header = FALSE, fill = TRUE)
#' results <- fread(file.path(bidui_dir, "results.txt"), sep = " ", header = FALSE, fill = TRUE)
#' class(results) <- class(result_t1) <- "data.frame"
#' rownames(results) <- paste("ID", results[, 1], sep = "=")
#' results <- results[result_t1[, 1], ]
#' result_rel <- add_strand(results)
#' result_rel <- get_real_loc(result_rel)
#' result_rel <- get_seq(result_rel)
#' }
get_seq <- function(result_rel, len, BSgenome_host, BSgenome_virus) {
    # 根据len重新修改start和end
    hostloc_start <- which(result_rel$start_host == result_rel$loc_host) # host -
    result_rel[hostloc_start, "end_host"] <- result_rel[hostloc_start, "start_host"] + (len-1)
    result_rel[hostloc_start, "start_host40"] <- result_rel[hostloc_start, "start_host"] - len
    result_rel[hostloc_start, "end_host40"] <- result_rel[hostloc_start, "start_host"] + (len-1)

    hostloc_end <- which(result_rel$end_host == result_rel$loc_host) # host +
    result_rel[hostloc_end, "start_host"] <- result_rel[hostloc_end, "end_host"] - (len-1)
    result_rel[hostloc_end, "start_host40"] <- result_rel[hostloc_end, "end_host"] - (len-1)
    result_rel[hostloc_end, "end_host40"] <- result_rel[hostloc_end, "end_host"] + len

    hpvloc_start <- which(result_rel$start_hpv == result_rel$loc_hpv) # hpv -
    result_rel[hpvloc_start, "end_hpv"] <- result_rel[hpvloc_start, "start_hpv"] + (len-1)
    result_rel[hpvloc_start, "start_hpv40"] <- result_rel[hpvloc_start, "start_hpv"] - len
    result_rel[hpvloc_start, "end_hpv40"] <- result_rel[hpvloc_start, "start_hpv"] + (len-1)


    hpvloc_end <- which(result_rel$end_hpv == result_rel$loc_hpv) # hpv +
    result_rel[hpvloc_end, "start_hpv"] <- result_rel[hpvloc_end, "end_hpv"] - (len-1)
    result_rel[hpvloc_end, "start_hpv40"] <- result_rel[hpvloc_end, "end_hpv"] - (len-1)
    result_rel[hpvloc_end, "end_hpv40"] <- result_rel[hpvloc_end, "end_hpv"] + len
    # len不能太离谱而导致跳出基因组范围外。
    for (i in 1:nrow(result_rel)) {
        result_rel$start_hpv[i] <- max(1, result_rel$start_hpv[i])
        result_rel$start_hpv40[i] <- max(1, result_rel$start_hpv40[i])
        result_rel$end_hpv[i] <- min(7905, result_rel$end_hpv[i])
        result_rel$end_hpv40[i] <- min(7905, result_rel$end_hpv40[i])
    }


    result_rel$seq_host <- result_rel$seq_host40 <- result_rel$seq_hpv <- result_rel$seq_hpv40 <- rep("1", nrow(result_rel))
    for (i in 1:nrow(result_rel)) {
        # 根据截断点，找到宿主和病毒序列
        seq_host <- getSeq(BSgenome_host, result_rel$chr[i], start = result_rel$start_host[i], end = result_rel$end_host[i]) |> as.character()
        seq_hpv <- getSeq(BSgenome_virus, "NC_001526.2", start = result_rel$start_hpv[i], end = result_rel$end_hpv[i]) |> as.character()
        if (result_rel$strand_host[i] == "-") seq_host <- get_comple(seq_host)
        if (result_rel$strand_hpv[i] == "-") seq_hpv <- get_comple(seq_hpv)

        seq_host40 <- getSeq(BSgenome_host, result_rel$chr[i], start = result_rel$start_host40[i], end = result_rel$end_host40[i]) |> as.character()
        seq_hpv40 <- getSeq(BSgenome_virus, "NC_001526.2", start = result_rel$start_hpv40[i], end = result_rel$end_hpv40[i]) |> as.character()
        if (result_rel$strand_host[i] == "-") seq_host40 <- get_comple(seq_host40)
        if (result_rel$strand_hpv[i] == "-") seq_hpv40 <- get_comple(seq_hpv40)
        result_rel$seq_host[i] <- seq_host
        result_rel$seq_host40[i] <- seq_host40
        result_rel$seq_hpv[i] <- seq_hpv
        result_rel$seq_hpv40[i] <- seq_hpv40
    }

    for (i in hostloc_start) { # host -
        result_rel$seq_host[i] <- get_rev(result_rel$seq_host[i])
        result_rel$seq_host40[i] <- get_rev(result_rel$seq_host40[i])
    }


    for (i in hpvloc_end) { # hpv +
        result_rel$seq_hpv[i] <- get_rev(result_rel$seq_hpv[i])
        result_rel$seq_hpv40[i] <- get_rev(result_rel$seq_hpv40[i])
    }
    return(result_rel)

}



#' Get microhomology(mh) length
#'
#' @param seq_df a dataframe of three rows, each containing the host sequence,
#' integration sequence, and virus sequence
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
#' @return dataframe
#' @export
#'
#' @examples
#' \dontrun{
#' bidui_dir <- "SurVirus_result"
#' bam_dir <- "SurVirus_result/readsx"
#'
#' result_t1 <- fread(file.path(bidui_dir, "results.t1.txt"), sep = " ", header = FALSE, fill = TRUE)
#' results <- fread(file.path(bidui_dir, "results.txt"), sep = " ", header = FALSE, fill = TRUE)
#' class(results) <- class(result_t1) <- "data.frame"
#' rownames(results) <- paste("ID", results[, 1], sep = "=")
#' results <- results[result_t1[, 1], ]
#' result_rel <- add_strand(results)
#' result_rel <- get_real_loc(result_rel)
#' result_rel <- get_seq(result_rel)
#' result_rel <- add_mh(result_rel)
#' }
add_mh <- function(result_rel) {
    mh <- rep(0, nrow(result_rel))
    for (i in seq_len(length(mh))) {
        seq_host <- result_rel[i, "seq_host40"]
        seq_hpv <- result_rel[i, "seq_hpv40"]
        seq_integration <- paste0(result_rel[i, "seq_host"], result_rel[i, "seq_hpv"])
        # seq_df <- strsplit(c(seq_host, seq_integration, seq_hpv), "") %>% do.call(rbind, .) %>% as.data.frame()
        seq_df <- strsplit(c(seq_host, seq_integration, seq_hpv), "")
        seq_df <- do.call(rbind, seq_df) |> as.data.frame()
        
        mh[i] <- get_mh(seq_df)
    }
    result_rel$mh <- mh
    return(result_rel)
}



#' get mh value in flank region of SurVirus result (result of get_seq())
#'
#' @param result_rel SurVirus result (result of get_seq())
#'
#' @param len flanking region size
#'
#' @return dataframe
#' @export
#'
#' @examples
#' \dontrun{
#' bidui_dir <- "SurVirus_result"
#' bam_dir <- "SurVirus_result/readsx"
#'
#' result_t1 <- fread(file.path(bidui_dir, "results.t1.txt"), sep = " ", header = FALSE, fill = TRUE)
#' results <- fread(file.path(bidui_dir, "results.txt"), sep = " ", header = FALSE, fill = TRUE)
#' class(results) <- class(result_t1) <- "data.frame"
#' rownames(results) <- paste("ID", results[, 1], sep = "=")
#' results <- results[result_t1[, 1], ]
#' result_rel <- add_strand(results)
#' result_rel <- get_real_loc(result_rel)
#' result_rel <- get_seq(result_rel)
#' result_rel <- add_mh_flank(result_rel, len = 5)
#' }
add_mh_flank <- function(result_rel, len) {
    result_rel$mh <- rep(0, nrow(result_rel))
    for (i in seq_len(nrow(result_rel))) {
        n1 <- n2 <-  len
        seq_host <- result_rel[i, "seq_host40"]
        seq_hpv <- result_rel[i, "seq_hpv40"]
        seq_integration <- paste0(result_rel[i, "seq_host"], result_rel[i, "seq_hpv"])
        # seq_df <- strsplit(c(seq_host, seq_integration, seq_hpv), "") %>% do.call(rbind, .) %>% as.data.frame()
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
        result_rel$mh[i] <- get_mh(seq_df)
    }
    return(result_rel)
}
