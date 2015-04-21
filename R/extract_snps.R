extract_snps <- function(snps, indir, chunkmap, chunkmap_cols = 1:3,
                         keep = NULL, idfile = NULL, pattern = "\\.gz$",
                         ncore = 1L, chr_chunk = ".*chr([^_]+)_(\\d+)")
{
    ## =================================================================
    ## Validate arguments.
    ## =================================================================

    ## Check for existence of mandatory arguments.
    force(indir)
    force(snps)
    force(chunkmap)

    ## Check that input directory exists.
    if (!is.directory(indir))
        stop("`indir' must be an existing directory with chunk files.")

    ## Check that all chunk map files exist.
    not_there <- ! file.exists(chunkmap)
    if (any(not_there))
        stop("The following files in `chunkmap' don't exist:\n",
             paste("\t", chunkmap[not_there], collapse = "\n"))
    rm(not_there)

    ## Check that at least 1 snp was selected for extraction.
    if (length(snps) == 0L)
        stop("You must specify at least 1 snp in `snps'.")

    ## Check that `keep' and `idfile' are consistent.
    if (length(keep)) {
        if (is.null(idfile))
            stop("`idfile' missing.")

        if (!file.exists(idfile))
            stop("`idfile' doesn't exist.")

        ids <- scan(idfile, character(), quiet = TRUE)

        not_there <- ! keep %in% ids
        if (any(not_there))
            stop("Some ids from `keep' are not among the ids in ",
                 "`idfile': ", paste(keep[not_there], collapse = ", "))
    }

    ## =================================================================
    ## Local helper functions.
    ## =================================================================

    ## Number of columns before imputed probability triples in impute
    ## output files.
    leading_columns <- 5L

    ## Build functions for extracting chromosome and chunk number from
    ## input chunk file names.
    chr   <- function(xs) submatch(chr_chunk, xs)[, 1L]
    chunk <- function(xs) submatch(chr_chunk, xs)[, 2L]

    ## Format snps and their corresponding chromosome and chunk in a
    ## table-like string.
    format_snps <- function(snp, chr, chunk)
    {
        fmt <- "%11s %3s %5s"
        x <- paste(c(sprintf(fmt, "snp", "chr", "chunk"),
                     sprintf(fmt, snp, chr, chunk)), collapse = "\n")
        x
    }

    ## =================================================================
    ## Get list of imputed chunk files.
    ## =================================================================
    chunk_files <- list.files(indir, pattern = pattern,
                              full.names = TRUE)
    files <- data.frame(
        chr   = chr(chunk_files),
        chunk = chunk(chunk_files),
        file  = chunk_files,
        stringsAsFactors = FALSE)

    ## =================================================================
    ## Read chunk map files.
    ## =================================================================
    read_3_columns <- function(f)
    {
        ## Only read snp, chr, and chunk columns.
        d <- data.table::fread(f, data.table = FALSE,
                               select = chunkmap_cols,
                               colClasses = "character")

        cols <- c("snp", "chr", "chunk")
        names(d) <- cols[order(chunkmap_cols)]

        d[d$snp %in% snps, cols, drop = FALSE]
    }

    pr("Reading `chunkmap' files ...")
    snp2chunk <- do.call(
        rbind, parallel::mclapply(chunkmap, read_3_columns,
                                  mc.cores = ncore))

    not_there <- ! snps %in% snp2chunk$snp
    if (any(not_there))
        stop("Some `snps' were not found in `chunkmap' files: ",
             paste(snps[not_there], collapse = ", "))
    rm(not_there)

    ## =================================================================
    ## Map snps to chunk files.
    ## =================================================================

    pr("Mapping snps to chunk files ...")
    d <- merge(snp2chunk, files, by = c("chr", "chunk"))

    not_there <- ! snp2chunk$snp %in% d$snp
    if (any(not_there))
        stop("Some `snps' have chr-chunk values in `chunkmap' files ",
             "that don't match any chunk file in `indir':\n",
             format_snps(snp2chunk$snp[not_there],
                         snp2chunk$chr[not_there],
                         snp2chunk$chunk[not_there]))
    rm(not_there)

    d <- d[order(d$chr, d$chunk), , drop = FALSE]

    ## =================================================================
    ## Extract snps from chunk files.
    ## =================================================================

    ## Check that all chunk files have the same number of columns.
    ncols <- vapply(
        vapply(d$file, gzhead, character(1L)),
        function(x)
        {
            length(strsplit(x, " ", fixed = TRUE)[[1L]])
        },
        integer(1L))
    stopifnot(all_neighbors(`==`, ncols))

    ## Check that the number of ids in `idfile' corresponds to the
    ## number of probability triples in the imputed chunk files.
    if (length(keep)) {
        nids <- length(ids)
        ntriples <- (ncols[1L] - leading_columns) %/% 3L
        if (nids != ntriples)
            stop("The number of individuals implied by the number of ",
                 "ids in `idfile' (", nids, ") differs from the ",
                 "number of probability triples in the imputed chunk ",
                 "files (", ntriples, ").")
        rm(nids, ntriples)
    }

    ## Group snps by chunk file.
    by_chunk <- split(d$snp, d$file)

    pr("Number of chunk files to search through: ", length(by_chunk))

    extract_from_chunk <- function(k)
    {
        snps <- by_chunk[[k]]
        chunk_file <- names(by_chunk)[k]
        con <- gzfile(chunk_file, "r")
        on.exit(close(con))
        d <- read.table(con, colClasses = "character",
                        stringsAsFactors = FALSE)
        pr1(".")
        d[d[[2]] %in% snps, , drop = FALSE]
    }

    ds <- parallel::mclapply(
        seq_along(by_chunk), extract_from_chunk, mc.preschedule = FALSE,
        mc.cores = ncore, mc.silent = FALSE)
    pr()

    ## Issue a warning unless all snps were found.
    missing_snps <- unlist(Map(setdiff, by_chunk, lapply(ds, `[[`, 2L)),
                           use.names = FALSE)
    if (length(missing_snps)) {
        not_there <- snp2chunk$snp %in% missing_snps
        warning("Some snps could not be found in the corresponding ",
                "chunk files:\n",
                format_snps(snp2chunk$snp[not_there],
                            snp2chunk$chr[not_there],
                            snp2chunk$chunk[not_there]))
    }

    do.call(rbind, ds)
}
