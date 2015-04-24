snptest <- function(indir, sample_file, outdir, pheno, covs,
                    exclusion_file = NULL, add_args = NULL, ncore = 1L,
                    pattern = "\\.gz$",
                    chr_chunk = ".*chr([^_]+)_(\\d+)",
                    executable = "snptest")
{
    ## =================================================================
    ## Local helper functions.
    ## =================================================================

    ## Build functions for extracting chromosome and chunk number from
    ## input chunk file names.
    extract <- function(k)
    {
        function(xs) maybe(as.integer, submatch(chr_chunk, xs)[, k])
    }
    chr   <- extract(1L)
    chunk <- extract(2L)

    ## Return TRUE if `logfile' ends in "finito".
    finito <- function(logfile)
    {
        x <- readLines(logfile, warn = FALSE)
        if (length(x) == 0L) FALSE
        else x[length(x)] == "finito"
    }

    ## Return TRUE if the chunk with snptest output file `outfile' and
    ## snptest log file `logfile' was successfully snptest'ed,
    ## otherwise return FALSE.
    chunk_finished <- function(outfile, logfile)
    {
        file.exists(outfile) && finito(logfile)
    }

    ## Return TRUE for every successfully tested chunk in `d', and
    ## FALSE otherwise.
    find_finished_chunks <- function(d)
    {
        unlist(parallel::mclapply(
            seq_len(nrow(d)),
            function(k) chunk_finished(d$output[k], d$log[k]),
            mc.cores = ncore))
    }

    ## =================================================================
    ## Validate arguments.
    ## =================================================================

    ## Check that all mandatory arguments were supplied.
    force(indir)
    force(sample_file)
    force(outdir)
    force(pheno)

    ## Force user to explicitly specify all covariates via `covs'.
    if (missing(covs) ||                           # not specified
        length(covs) == 0L ||                      # 0-length
        all(grepl("^\\s*$", covs, perl = TRUE)) || # all empty
        ## attempt to override `covs'
        (!is.null(add_args) && any(grepl("-cov_all", add_args))))
        stop("You must specify the names of all covariates as used in ",
             "the `sample_file' using the `covs' argument to ",
             "`snptest'.  Something like \"-cov_all\" won't work :)")

    indir  <- unslash(indir)
    outdir <- unslash(outdir)

    stopifnot(file.exists(sample_file),
              is.directory(indir))

    if (!is.null(exclusion_file))
        stopifnot(file.exists(exclusion_file))

    ## Check that snptest executable is known.
    if (system(paste("which", executable),
               ignore.stdout = TRUE,
               ignore.stderr = TRUE) != 0L)
        stop("Cannot find snptest executable: ", executable)

    ## =================================================================
    ## Create table of all file names involved.
    ## =================================================================
    logdir <- file.path(outdir, "log")
    infiles <- list.files(indir, pattern = pattern, full.names = TRUE)
    files <- data.frame(
        chr   = chr(infiles),
        chunk = chunk(infiles),
        input = infiles,
        stringsAsFactors = FALSE)
    files$stub   <- paste0("chr", files$chr, "_", files$chunk)
    files$output <- file.path(outdir, paste0(files$stub, ".txt"))
    files$log    <- file.path(logdir, paste0(files$stub, ".log"))
    files$stub   <- NULL

    ## Sort by chromosome and chunk.
    files <- files[order(files$chr, files$chunk), ]

    ## Create directories for the various kinds of output files.
    mkdir(c(outdir, logdir))

    ## =================================================================
    ## Run snptest on as yet unfinished chunks.
    ## =================================================================

    pr("Checking for successfully snptest'ed chunks ...")
    files$done <- find_finished_chunks(files)

    if (!all(files$done)) {

        pr("Number of chunks already successfully snptest'ed: ", sum( files$done))
        pr("Number of chunks not yet successfully snptest'ed: ", sum(!files$done))

        d <- files[!files$done, , drop = FALSE]

        ## Assemble snptest commands, one per chunk.
        covsum <- paste(covs, collapse = " + ")
        analysis_name <- shQuote(paste(pheno, "~ 1 +", covsum))

        cmd <- paste(executable,
                     "-analysis_name", analysis_name,
                     "-data", d$input, sample_file,
                     "-log", d$log,
                     "-o", d$output,
                     "-pheno", pheno,
                     if (!is.null(exclusion_file))
                         paste("-exclude_samples", exclusion_file),
                     "-cov_names", paste(covs, collapse = " "),
                     if (!missing(add_args))
                         paste(add_args, collapse = " "))

        ncore <- min(ncore, length(cmd))
        pr("Running snptest on ", ncore, " cores ...")
        parallel::mclapply(cmd, system, mc.preschedule = FALSE,
                           mc.silent = TRUE, mc.cores = ncore)

        pr("Updating list of successfully snptest'ed chunks ...")
        files$done <- find_finished_chunks(files)

        ## Warn about failed chunks.
        for (f in files$input[!files$done])
            warning("snptest failed for ", f, immediate. = TRUE)

        pr("Percentage of chunks successfully snptest'ed: ",
           sprintf("%.1f%%", 100 * sum(files$done) / nrow(files)))
    }
    else {
        pr("All chunks have already been successfully snptest'ed.")
    }

    files
}

summarize_snptest <- function(filename, chr = NULL, old2new = NULL,
                              select = !is.null(old2new), hook = NULL)
{
    ## =================================================================
    ## Local functions.
    ## =================================================================

    ## Create column names.
    variable <- function(summary_measure)
    {
        ## `test_type' and `genetic_model' are defined and looked up
        ## in the environment in which `variable' was defined (see
        ## below).
        paste0(test_type, "_", genetic_model, "_", summary_measure)
    }
    if (select && is.null(old2new))
        stop("You must specify `old2new' to select columns.")

    ## =================================================================
    ## Read snptest table.
    ## =================================================================

    ## Let `fread' detect column classes.  This is safe, since `fread'
    ## will never coerce a column from a higher to a lower type.  In
    ## other words, if we get values of a type lower than the `true'
    ## type of the column, it's because all values in the column can
    ## be represented in the lower type without data loss.  Since the
    ## values in the lower type represent the same data, we can treat
    ## them as being equivalent to the higher type for all practical
    ## purposes, e.g., arithmetic, character operations, rbind, etc.
    ## Letting `fread' select the lowest possible type ensures that
    ## the memory footprint of the resulting table will be minimal.
    d <- data.table::fread(filename, data.table = FALSE)

    ## =================================================================
    ## Determine type of snptest analysis performed.
    ## =================================================================
    test_types <- c("frequentist", "bayesian")
    genetic_models <- c("add", "dom", "rec", "gen", "het")

    x <- paste0("^", test_types)
    test_type <- test_types[find_first_match(x, names(d))]

    if (is.na(test_type))
        stop("Unable to determine \"test type\" from: ",
             paste(names(d), collapse = ", "))

    x <- paste0("^", test_type, "_", genetic_models)
    genetic_model <- genetic_models[find_first_match(x, names(d))]

    if (is.na(genetic_model))
        stop("Unable to determine \"genetic model\" from: ",
             paste(names(d), collapse = ", "))

    ## =================================================================
    ## Add additional columns.
    ## =================================================================
    if (!is.null(chr))
        d$chromosome <- chr

    ## Allele frequency of coding allele, i.e., `alleleB' in snptest.
    d$freq_alleleB <- with(d, (all_AB + (2 * all_BB))
                           / (2 * (all_AA + all_AB + all_BB)))

    d$callrate <- with(d, 1 - all_NULL
                       / (all_AA + all_AB + all_BB + all_NULL))

    d$imputed <- as.integer(d$alternate_ids == "---")

    ## =================================================================
    ## Rename, select, and reorder.
    ## =================================================================
    hwe <- "cohort_1_hwe"
    cols <- c("rsid", "chromosome", "position", "alleleA", "alleleB",
              "average_maximum_posterior_call", "info", "all_total",
              "missing_data_proportion", if (hwe %in% names(d)) hwe,
              variable(c("beta_1", "se_1", "pvalue")),
              "freq_alleleB", "imputed", "callrate")

    d[, cols, drop = FALSE]
    if (!is.null(old2new)) {
        ## If `select' is TRUE, there might be unnamed elements in
        ## `old2new' whose only purpose is to select and reorder
        ## columns.  In that case, we want to ignore any warnings due
        ## to non
        if (select && (is.null(names(old2new))
                       || "" %in% names(old2new))) {
            names(d) <- rename(d, old2new, warn = FALSE)
        }
        else {
            ## Warnings that are not due to the above circumstances
            ## should not be muffled.
            names(d) <- rename(d, old2new)
        }

        ## Select and reorder a subset of colums.
        if (select)
            d <- d[, old2new, drop = FALSE]
    }

    ## =================================================================
    ## Feed data frame to user-supplied function before returning.
    ## =================================================================
    if (!is.null(hook))
        return(hook(d))

    d
}

combine_snptest <- function(indir, outdir, pattern, ncore = 1L,
                              old2new, reorder = FALSE, gzip = FALSE)
{
    ## =================================================================
    ## Validate user-supplied arguments.
    ## =================================================================

    ## Remove trailing slashes from directories.
    indir  <- unslash(indir)
    outdir <- unslash(outdir)

    ## Check that directories exist.
    stopifnot(is.directory(indir))
    mkdir(outdir)
    stopifnot(is.directory(outdir))

    if (reorder && missing(old2new))
        stop("You must specify `old2new' to reorder columns.")

    if (gzip && !grepl("\\.gz$", pattern))
        stop("You want to use gzip: `pattern' must end in \".gz\"")

    ## =================================================================
    ## Find output file prefix and suffix from `pattern'.
    ## =================================================================
    placeholder <- "CHROMOSOME"
    stopifnot(grepl(placeholder, pattern, fixed = TRUE))

    x <- unlist(strsplit(pattern, placeholder, fixed = TRUE))
    stopifnot(length(x) == 2L)
    prefix <- x[1L]
    suffix <- x[2L]

    ## =================================================================
    ## Put snptest files into data frame to simplify manipulation.
    ## =================================================================
    infiles <- list.files(indir, "\\.txt", full.names = TRUE)
    files <- data.frame(
        chr   = chr(infiles),
        chunk = chunk(infiles),
        input = infiles,
        stringsAsFactors = FALSE)
    files <- files[order(files$chr, files$chunk), ]

    ## =================================================================
    ## Process snptest tables by chromosome.
    ## =================================================================
    for (chr in sort(unique(files$chr))) {

        pr("Processing snptest output files for chromosome ", chr)

        d <- do.call(rbind, parallel::mclapply(
            files$input[files$chr == chr],
            summarize_snptest, chr, mc.preschedule = FALSE,
            mc.cores = ncore, mc.silent = FALSE))

        ## Rename and reorder columns.
        if (!missing(old2new)) {
            names(d) <- rename(d, old2new)
            if (reorder)
                d <- d[, old2new, drop = FALSE]
        }

        ## Write to disk.
        outfile <- file.path(outdir, paste0(prefix, chr, suffix))

        if (gzip) con <- gzfile(outfile)
        else      con <- file(outfile)

        write.table(d, con, quote = FALSE, sep = "\t",
                    row.names = FALSE, col.names = TRUE)
    }
}
