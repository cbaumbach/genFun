snptest <- function(indir, sample_file, exclusion_file, outdir, pheno,
                    covs, add_args = NULL, ncore = 1L,
                    executable = "snptest")
{
    ## =================================================================
    ## Local helper functions.
    ## =================================================================

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
    force(exclusion_file)
    force(outdir)
    force(pheno)
    force(covs)

    ## Remove trailing slash in directory names.
    indir  <- unslash(indir)
    outdir <- unslash(outdir)

    ## Check that all input files and directories exist.
    stopifnot(all(vapply(c(indir, sample_file, exclusion_file),
                         file.exists, logical(1L))),
              is.directory(indir))

    ## Check that snptest executable is known.
    if (system(paste("which", executable),
               ignore.stdout = TRUE,
               ignore.stderr = TRUE) != 0L)
        stop("Cannot find snptest executable: ", executable)

    ## =================================================================
    ## Create table of all file names involved.
    ## =================================================================
    logdir <- file.path(outdir, "log")
    infiles <- list.files(indir, pattern = "\\.gz$", full.names = TRUE)
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
                     "-data", d$input, sample_file,
                     "-log", d$log,
                     "-o", d$output,
                     "-pheno", pheno,
                     "-exclude_samples", exclusion_file,
                     "-cov_names", paste(covs, collapse = " "),
                     "-analysis_name", analysis_name,
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
