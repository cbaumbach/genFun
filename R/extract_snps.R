extract_snps <- function(snps, indir, outdir, chunkmap, idfile,
                         keep = NULL, chunkmap_cols = 1:3,
                         pattern = "\\.gz$", ncore = 1L,
                         chr_chunk = ".*chr([^_]+)_(\\d+)")
{
    t0 <- Sys.time()                    # record starting time

    ## =================================================================
    ## Validate arguments.
    ## =================================================================

    ## Check for existence of mandatory arguments.
    force(indir)
    force(snps)
    force(outdir)
    force(chunkmap)
    force(idfile)

    ## Check that input directory exists.
    if (!is.directory(indir))
        stop("`indir' must be an existing directory with chunk files.")

    if (is.factor(snps))
        snps <- as.character(snps)

    if (missing(outdir))
        stop("You must specify an directory for output via `outdir'.")
    else
        outdir <- unslash(outdir)

    ## Check if file with order of individuals already exists.
    id_outfile <- file.path(outdir, "order_of_individuals.txt")

    if (file.exists(outdir) && is.directory(outdir)) {
        if (file.exists(id_outfile))
            stop("File already exists: ", id_outfile)
    }
    else mkdir(outdir)

    ## Check that `pattern' contains exactly two parenthesized
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

        ## Convert factor snp ids to character.
        if (is.factor(chunkmap[[chunkmap_cols[1L]]]))
            chunkmap[[chunkmap_cols[1L]]] <- as.character(chunkmap[[chunkmap_cols[1L]]])
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
        ## Choose optimal column width.
        fmt <- paste0("%", max(nchar(snp),   nchar("snp")),   "s ",
                      "%", max(nchar(chr),   nchar("chr")),   "s ",
                      "%", max(nchar(chunk), nchar("chunk")), "s")

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

    ## Create chunk -> outfile map.
    chunk_files <- unique(files$file)
    chunk2outfile <- file.path(outdir, paste0("chr", chr(chunk_files),
                                              "_", chunk(chunk_files),
                                              ".txt.gz"))
    names(chunk2outfile) <- chunk_files
    rm(chunk_files)

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
                paste(snps_not_in_chunkmap, collapse = ", "),
                immediate. = TRUE)
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
                            chunkmap$chunk[not_there]),
                immediate. = TRUE)
    }
    rm(not_there)

    d <- d[order(d$chr, d$chunk), , drop = FALSE]

    ## Stop if `outdir' contains files that would be overwritten.
    future_outfiles <- chunk2outfile[unique(d$file)]
    existing_outfiles <- future_outfiles[file.exists(future_outfiles)]
    if (length(existing_outfiles) > 0L)
        stop("Some chunks seem to have already been extracted into ",
             "`outdir':\n\t", paste0(existing_outfiles, collapse = "\n\t"))
    rm(future_outfiles, existing_outfiles)

    ## =================================================================
    ## Determine columns to be extracted from chunk files.
    ## =================================================================

    ## Number of columns before imputed probability triples in impute
    ## chunk files.
    leading_columns <- 5L

    if (extract_individuals) {
        selected_columns <- c(seq_len(leading_columns),
                              pos2col(which(ids %in% keep)))

        ## Write selected columns to temporary file to be read by perl
        ## script.
        column_file <- tempfile("selected_columns_", tmpdir = "/tmp")
        cat(selected_columns, sep = "\n", file = column_file)
    }


    ## Write order of individuals to file in output directory.
    if (extract_individuals) {
        cat(ids[ids %in% keep], sep = "\n", file = id_outfile)
    }
    else {
        cat(ids, sep = "\n", file = id_outfile)
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
    perl_script <- system.file(file.path("perl", "extract_from_chunk.pl"),
                               package = "genFun", mustWork = TRUE)

    ## Create temporary file names before running in parallel.
    snp_files <- tempfile(
        paste0("extract_snps_from_chunk_", seq_along(by_chunk),
               "_of_", length(by_chunk), "_"),
        tmpdir = "/tmp")

    cmd <- paste(perl_script,
                 "-s", snp_files,
                 "-c", names(by_chunk),
                 "-o", chunk2outfile[names(by_chunk)])

    if (extract_individuals)
        cmd <- paste(cmd, "-i", column_file)

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
            snps_not_found <- scan(con, character(), quiet = TRUE)
        }, finally = {
            file.remove(snp_files[k])
        })

        pr1(".")

        snps_not_found
    }

    pr("Using ", ncore, " core", if (ncore > 1) "s", " to search ",
       length(snps), " snp", if (length(snps) > 1) "s", " in ",
       length(by_chunk), " chunk", if (length(by_chunk) > 1) "s")

    pr1("Chunks: ")
    snps_not_found_in_chunks <- unlist(parallel::mclapply(
        seq_along(by_chunk), extract_from_chunk, mc.preschedule = FALSE,
        mc.cores = ncore, mc.silent = FALSE, mc.allow.recursive = FALSE))
    pr(" [done]")

    if (extract_individuals)
        file.remove(column_file)

    ## Warn if some snps could not be found in "their" chunks.
    if (length(snps_not_found_in_chunks)) {
        not_there <- d$snp %in% snps_not_found_in_chunks
        warning("Some snps could not be found in the corresponding ",
                "chunk files:\n",
                format_snps(d$snp[not_there],
                            d$chr[not_there],
                            d$chunk[not_there]),
                immediate. = TRUE)
        rm(not_there)
    }

    ## Attach data frame with snps that could not be found to return
    ## value.  Mention the reason why a snp was not found in a
    ## `comment' column.
    missing_snps <- NULL
    if (length(snps_not_in_chunkmap) ||
        length(snps_in_nonexistent_chunks) ||
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

        missing_snps <- rbind(d1, d2, d3)
        rownames(missing_snps) <- NULL
    }

    print(Sys.time() - t0)

    missing_snps
}
