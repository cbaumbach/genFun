snptest <- function(indir, sample_file, outdir, pheno, covs = NULL,
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

    ## Force the user to be explicit about which covariates he wants
    ## to include in the model.  Intercept attempts to specify all
    ## covariates via "-cov_all".
    if (!is.null(add_args) && any(grepl("-cov_all", add_args)))
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

        ## Assemble shell commands, one per input chunk file.
        cmd <- paste(
            executable,

            "-analysis_name",
            if (is.null(covs))
                shQuote(paste(pheno, "~ 1 + snp"))
            else
                shQuote(paste(pheno, "~ 1 +",
                              paste(covs, collapse = " + "),
                              "+ snp")),

            "-data", d$input, sample_file,
            "-log", d$log,
            "-o", d$output,
            "-pheno", pheno,

            if (!is.null(exclusion_file))
                paste("-exclude_samples", exclusion_file),

            if (!is.null(covs))
                paste("-cov_names", paste(covs, collapse = " ")),

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
    d <- data.table::fread(filename, data.table = FALSE,
                           verbose = FALSE, showProgress = FALSE)

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

combine_snptest <- function(indir, outdir, ncore = 1L, old2new = NULL,
                            select = !is.null(old2new),
                            pattern = "^chr.*\\.txt$",
                            template = "chr<CHROMOSOME>.txt",
                            chr_chunk = ".*chr([^_]+)_(\\d+)",
                            gzip = FALSE, hook = NULL,
                            overwrite = FALSE)
{
    ## =================================================================
    ## Validate user-supplied arguments.
    ## =================================================================

    ## Remove trailing slashes from directories.
    indir  <- unslash(indir)
    outdir <- unslash(outdir)

    if (nsubexp(chr_chunk) != 2L)
        stop("`chr_chunk' must contain exactly 2 parenthesized ",
             "subexpressions.")

    if (select && is.null(old2new))
        stop("You must specify `old2new' to select columns.")

    if (gzip && !grepl("\\.gz$", template))
        stop("You want to use gzip: `template' must end in \".gz\"")

    stopifnot(is.directory(indir))
    mkdir(outdir)

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

    ## =================================================================
    ## Find output file prefix and suffix from `template'.
    ## =================================================================
    placeholder <- "<CHROMOSOME>"
    stopifnot(grepl(placeholder, template, fixed = TRUE))

    x <- unlist(strsplit(template, placeholder, fixed = TRUE))
    stopifnot(length(x) == 2L)
    prefix <- x[1L]
    suffix <- x[2L]

    ## =================================================================
    ## Put snptest files into data frame to simplify manipulation.
    ## =================================================================
    infiles <- list.files(indir, pattern, full.names = TRUE)
    files <- data.frame(
        chr   = chr(infiles),
        chunk = chunk(infiles),
        input = infiles,
        stringsAsFactors = FALSE)
    files <- files[order(files$chr, files$chunk), ]

    ## =================================================================
    ## Process snptest output.
    ## =================================================================

    pr1("Chromosomes: ")
    for (chrom in sort(unique(files$chr))) {

        pr1(chrom, " ")

        outfile <- file.path(outdir, paste0(prefix, chrom, suffix))

        if (!overwrite && file.exists(outfile)) {
            pr()
            pr("File already exists: ", outfile)
            pr("Skipping ...")
            next
        }

        d <- do.call(rbind, parallel::mclapply(
            files$input[files$chr == chrom],
            summarize_snptest,
            chr     = chrom,
            old2new = old2new,
            select  = select,
            hook    = hook,
            mc.cores           = ncore,
            mc.silent          = FALSE,
            mc.preschedule     = FALSE,
            mc.allow.recursive = FALSE))

        if (gzip) con <- gzfile(outfile)
        else      con <- file(outfile)

        write.table(d, con, quote = FALSE, sep = "\t",
                    row.names = FALSE, col.names = TRUE)
    }
    pr()
}
