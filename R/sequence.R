#' Use the flag value to detect which of the two ends the sequence belongs to
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
        
        seq_host <- getSeq(BSgenome_host, result_rel$chr[i], start = result_rel$start_host[i], 
            end = result_rel$end_host[i]) |> as.character()
        seq_hpv <- getSeq(BSgenome_virus, "NC_001526.2", start = result_rel$start_hpv[i], 
            end = result_rel$end_hpv[i]) |> as.character()
        if (result_rel$strand_host[i] == "-") seq_host <- get_comple(seq_host)
        if (result_rel$strand_hpv[i] == "-") seq_hpv <- get_comple(seq_hpv)
        result_rel$seq_host[i] <- seq_host
        result_rel$seq_hpv[i] <- seq_hpv
    }


    for (i in hostloc_start) {
        result_rel$seq_host[i] <- get_rev(result_rel$seq_host[i])
    }


    # Loop search, host from back to front, hpv from front to back, 
    # find all sequences that can match, and pick the longest one
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

  
    result_rel[hpvloc_start, "start_hpv"] <- result_rel[hpvloc_start, "start_hpv"] + 
        result_rel[hpvloc_start, "ins_mh"]
    result_rel[hpvloc_start, "end_hpv"] <- result_rel[hpvloc_start, "start_hpv"] + 19
    result_rel[hpvloc_end, "end_hpv"] <- result_rel[hpvloc_end, "end_hpv"] - 
        result_rel[hpvloc_end, "ins_mh"]
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
#' data(BSgenomehpv16)
#' SurVirus_dir <- system.file(file.path("extdata", "survirus_result"), package = "mhAnalysis")
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
get_seq <- function(result_rel, len = 10, BSgenome_host, BSgenome_virus) {
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
        seq_host <- getSeq(BSgenome_host, result_rel$chr[i], start = result_rel$start_host[i], 
            end = result_rel$end_host[i]) |> as.character()
        seq_hpv <- getSeq(BSgenome_virus, "NC_001526.2", start = result_rel$start_hpv[i], 
            end = result_rel$end_hpv[i]) |> as.character()
        if (result_rel$strand_host[i] == "-") seq_host <- get_comple(seq_host)
        if (result_rel$strand_hpv[i] == "-") seq_hpv <- get_comple(seq_hpv)

        seq_host40 <- getSeq(BSgenome_host, result_rel$chr[i], start = result_rel$start_host40[i], 
            end = result_rel$end_host40[i]) |> as.character()
        seq_hpv40 <- getSeq(BSgenome_virus, "NC_001526.2", start = result_rel$start_hpv40[i], 
            end = result_rel$end_hpv40[i]) |> as.character()
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
