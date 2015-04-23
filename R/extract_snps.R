extract_snps <- function(snps, indir, chunkmap, chunkmap_cols = 1:3,
                         keep = NULL, idfile = NULL, pattern = "\\.gz$",
                         ncore = 1L, chr_chunk = ".*chr([^_]+)_(\\d+)")
{
    t0 <- Sys.time()                    # record starting time

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

    ## Check that `pattern' contain exactly two parenthesized
    ## subexpressions.
    if (nsubexp(chr_chunk) != 2L)
        stop("`chr_chunk' must contain exactly 2 parenthesized ",
             "subexpressions, one for chromosome and one for chunk.")

    ## Check that `chunkmap' is a valid data frame or file name.
    if (nuniq(chunkmap_cols) != 3L)
        stop("`chunkmap_cols' must give the positions of the columns ",
             "corresponding to snp, chr, and chunk (in this order).")

    if (is.data.frame(chunkmap)) {
        if (ncol(chunkmap) < 3L)
            stop("`chunkmap' data frame must have at least 3 colums.")

        if (!all(chunkmap_cols %in% seq_len(ncol(chunkmap))))
            stop("`chunkmap_cols' must be a valid 3-element index for ",
                 "`chunkmap'.")
    }
    else if (is.character(chunkmap)) {
        chunkmap_files <- chunkmap
        not_there <- ! file.exists(chunkmap_files)
        if (any(not_there))
            stop("The following `chunkmap' files don't exist:\n",
                 paste("\t", chunkmap_files[not_there], collapse = "\n"))
        rm(not_there)
    }
    else {
        stop("`chunkmap' must be a data frame or character vector of ",
             "file names.")
    }

    ## Check that at least 1 snp was selected for extraction.
    if (length(snps) == 0L)
        stop("You must specify at least 1 snp in `snps'.")

    ## Check that `keep' and `idfile' are consistent.
    extract_individuals <- FALSE
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

        extract_individuals <- TRUE
    }

    ## =================================================================
    ## Local helper functions.
    ## =================================================================

    ## Build functions for extracting chromosome and chunk number from
    ## input chunk file names.
    chr   <- function(xs) submatch(chr_chunk, xs)[, 1L]
    chunk <- function(xs) submatch(chr_chunk, xs)[, 2L]

    ## Format snps and their corresponding chromosome and chunk in a
    ## table-like string.
    format_snps <- function(snp, chr, chunk)
    {
        ## Choose optimal width for snp column.
        fmt <- paste0("%", max(nchar(snp)), "s %3s %5s")

        ## Truncate after `max_lines' of output.
        max_lines <- 10L
        idx <- seq_len(min(length(snp), max_lines))
        too_long <- length(snp) > max_lines

        paste(c(sprintf(fmt, "snp", "chr", "chunk"),
                sprintf(fmt, snp[idx], chr[idx], chunk[idx]),
                if (too_long) "[... truncated]"),
              collapse = "\n")
    }

    ## Returns the columns in the chunk files corresponding to the
    ## individual ids that can be found in rows `xs' of `idfile'.
    pos2col <- function(xs)
    {
        step <- 3L         # probabilities come in triples
        offset <- 5L       # number of leading non-probability columns

        g <- function(x)
        {
            seq(step * (x-1L) + 1L, step * x)
        }

        offset + unlist(lapply(xs, g), use.names = FALSE)
    }

    take_3_columns <- function(d)
    {
        if (is.character(d)) {
            ## Only read snp, chr, and chunk columns.  `fread' doesn't
            ## care about the order of columns in `select'.  Whether
            ## it's 1:3 or 3:1, either way it'll come out as 1:3.  Use
            ## an explicit `sort' to dispell false expectations.  The
            ## snp, chr, and chunk columns will be in the same
            ## (relative) order as they were in the original table
            ## from which they were read.
            d <- data.table::fread(d, data.table = FALSE,
                                   select = sort(chunkmap_cols),
                                   colClasses = "character")
        }
        else if (is.data.frame(d)) {
            ## Remove superfluous columns.  Use `sort(chunkmap_cols)'
            ## for subsetting so as not to modify the order of the
            ## remaining columns.
            if (ncol(d) > 3L)
                d <- d[sort(chunkmap_cols)]
        }
        else {
            stop("`d' must be either a data frame or a character vector.")
        }

        ## Assign column names according to `chunkmap_cols'.  For this
        ## to work it was crucial that we did not modify the order of
        ## the snp, chr, and chunk columns in the above if-statement.
        cols <- c("snp", "chr", "chunk")
        names(d) <- cols[order(chunkmap_cols)]

        d[d$snp %in% snps, cols, drop = FALSE]
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
    if (is.data.frame(chunkmap)) {
        chunkmap <- take_3_columns(chunkmap)
    }
    else {
        pr("Reading `chunkmap' file",
           if (length(chunkmap_files) > 1L) "s", " ...")
        chunkmap <- do.call(
            rbind, parallel::mclapply(chunkmap_files, take_3_columns,
                                      mc.cores = ncore,
                                      mc.allow.recursive = FALSE))
    }

    not_there <- ! snps %in% chunkmap$snp
    snps_not_in_chunkmap <- character(0L)
    if (any(not_there)) {
        snps_not_in_chunkmap <- snps[not_there]
        warning("Some `snps' were not in `chunkmap': ",
                paste(snps_not_in_chunkmap, collapse = ", "))
    }
    rm(not_there)

    ## =================================================================
    ## Map snps to chunk files.
    ## =================================================================

    pr("Mapping snps to chunk files ...")
    d <- merge(chunkmap, files, by = c("chr", "chunk"))

    not_there <- ! chunkmap$snp %in% d$snp
    snps_in_nonexistent_chunks <- character(0L)
    if (any(not_there)) {
        snps_in_nonexistent_chunks <- chunkmap$snp[not_there]
        warning("Some `snps' have chr-chunk values in `chunkmap' files ",
                "that don't match any chunk file in `indir':\n",
                format_snps(chunkmap$snp[not_there],
                            chunkmap$chr[not_there],
                            chunkmap$chunk[not_there]))
    }
    rm(not_there)

    d <- d[order(d$chr, d$chunk), , drop = FALSE]

    ## =================================================================
    ## Determine columns to be extracted from chunk files.
    ## =================================================================

    ## Number of columns before imputed probability triples in impute
    ## chunk files.
    leading_columns <- 5L

    if (extract_individuals) {
        selected_columns <- c(seq_len(leading_columns),
                              pos2col(match(keep, ids)))
    }

    ## =================================================================
    ## Extract snps from chunk files.
    ## =================================================================

    ## Check that all chunk files have the same number of columns.
    ncols <- vapply(unique(d$file),
                    function(f)
                    {
                        con <- gzfile(f, "rb")
                        on.exit(close(con))
                        nfields(con, sep = " ")
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

    ## Find path to perl script that does the heavy lifting.
    perl_script <- file.path(find.package("genFun"),
                             "perl", "extract_from_chunk.pl")

    ## Create temporary file names before running in parallel.
    snp_files <- tempfile(rep_len("extract_snps_", length(by_chunk)),
                          tmpdir = "/tmp")

    cmd <- paste(perl_script, "-s", snp_files, "-c", names(by_chunk))

    ## Read the first few columns as type "character" and the
    ## probability triples as type "double".
    colCl <- rep_len("double", ncols[1L])
    colCl[seq_len(leading_columns)] <- "character"

    extract_from_chunk <- function(k)
    {
        cat(by_chunk[[k]], file = snp_files[k], sep = "\n")

        ## Create a closed connection.  `read.table' will open it,
        ## read from it, and close it when done reading.
        con <- pipe(cmd[k])

        tryCatch({
            d <- read.table(con, colClasses = colCl)
        }, finally = {
            file.remove(snp_files[k])
        })

        pr1(".")

        if (extract_individuals)
            d[d[[2]] %in% snps, selected_columns, drop = FALSE]
        else
            d[d[[2]] %in% snps, , drop = FALSE]
    }

    pr("Using ", ncore, " core", if (ncore > 1) "s", " to search ",
       length(snps), " snp", if (length(snps) > 1) "s", " in ",
       length(by_chunk), " chunk", if (length(by_chunk) > 1) "s")

    pr1("Chunks: ")
    ds <- parallel::mclapply(
        seq_along(by_chunk), extract_from_chunk, mc.preschedule = FALSE,
        mc.cores = ncore, mc.silent = FALSE, mc.allow.recursive = FALSE)
    pr(" [done]")

    ## Label chunks as `chr<chromosome>_<chunk>'.
    names(ds) <- paste0("chr", chr(names(by_chunk)),
                        "_", chunk(names(by_chunk)))

    ## Warn if some snps could not be found in "their" chunks.
    found_snps <- unique(unlist(lapply(ds, `[[`, 2L), use.names = FALSE))
    snps_not_found_in_chunks <- setdiff(d$snp, found_snps)
    if (length(snps_not_found_in_chunks)) {
        not_there <- d$snp %in% snps_not_found_in_chunks
        warning("Some snps could not be found in the corresponding ",
                "chunk files:\n",
                format_snps(d$snp[not_there],
                            d$chr[not_there],
                            d$chunk[not_there]))
        rm(not_there)
    }

    ## Attach data frame with snps that could not be found to return
    ## value.  Mention the reason why a snp was not found in a
    ## `comment' column.
    if (length(snps_not_in_chunkmap) |
        length(snps_in_nonexistent_chunks) |
        length(snps_not_found_in_chunks)) {

        cols <- c("snp", "chr", "chunk")
        d1 <- chunkmap[chunkmap$snp %in% snps_in_nonexistent_chunks,
                        cols, drop = FALSE]
        d1$comment <- rep_len("nonexistent chunk file", nrow(d1))

        d2 <- chunkmap[chunkmap$snp %in% snps_not_found_in_chunks,
                        cols, drop = FALSE]
        d2$comment <- rep_len("not found in chunk file", nrow(d2))

        d3 <- data.frame(snp = snps_not_in_chunkmap,
                         stringsAsFactors = FALSE)
        d3$chr     <- rep_len(NA_character_,     nrow(d3))
        d3$chunk   <- rep_len(NA_character_,     nrow(d3))
        d3$comment <- rep_len("not in chunkmap", nrow(d3))

        d <- rbind(d1, d2, d3)
        rownames(d) <- NULL
        attr(ds, "missing_snps") <- d
    }

    print(Sys.time() - t0)

    ds
}
