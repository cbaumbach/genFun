context("Extracting snps from gz chunk files.")

test_that("extract_snps works", {

    ## ===============================================================
    ## Setup
    ## ===============================================================

    snps <- paste0("rs", c(1:4, 10:14))
    keep <- c(1:10, 20:30, 40:50)

    ## Files.
    indir <- "data"
    infiles <- list.files(path = indir, pattern = "^chr.*\\.txt$", full.names = TRUE)
    infiles_gz <- paste0(infiles, ".gz")
    outdir <- tempdir()
    chunkmap_file <- file.path(indir, "list_of_snps.txt")
    idfile   <- file.path(indir, "list_of_individuals.txt")
    outfiles <- file.path(outdir, c(basename(infiles_gz), "order_of_individuals.txt"))

    ## Create compressed chunk files.
    cmds <- paste("gzip -c", infiles, ">", paste0(infiles, ".gz"))
    stopifnot(all(vapply(cmds, system, integer(1L)) == 0L))


    ## ===============================================================
    ## Construct expected output tables from input files.
    ## ===============================================================

    stopifnot(all(file.exists(infiles)))

    din1 <- utils::read.table(gzfile(infiles[1L]), stringsAsFactors = FALSE)
    din2 <- utils::read.table(gzfile(infiles[2L]), stringsAsFactors = FALSE)

    ## After a leading set of 5 columns, individual-specific genotype
    ## probabilities come in triples.
    columns <- c(1:5, 5L + rep((keep - 1L) * 3L, each = 3L) + 1:3)
    ddin1 <- din1[, columns, drop = FALSE]
    ddin2 <- din2[, columns, drop = FALSE]


    ## ---------------------------------------------------------------
    ## Test 1: reading chunkmap from file
    ## ---------------------------------------------------------------
    missing_snps <- extract_snps(
        snps   = snps,
        indir  = indir,
        outdir = outdir,
        keep   = keep,
        idfile = idfile,
        ncore  = 2L,
        chunkmap = chunkmap_file,
        chunkmap_cols = c(1L, 2L, 4L),
        quiet = TRUE)

    expect_that(is.null(missing_snps), is_true())
    stopifnot(all(file.exists(outfiles)))

    dout1 <- utils::read.table(gzfile(outfiles[1L]), stringsAsFactors = FALSE)
    dout2 <- utils::read.table(gzfile(outfiles[2L]), stringsAsFactors = FALSE)

    file.remove(outfiles)

    expect_that(dout1, is_equivalent_to(ddin1))
    expect_that(dout2, is_equivalent_to(ddin2))

    rm(missing_snps, dout1, dout2)

    ## ---------------------------------------------------------------
    ## Test 2: reading chunkmap table with snp as factor column
    ## ---------------------------------------------------------------

    chunkmap <- utils::read.table(chunkmap_file, header = TRUE)

    stopifnot(is.factor(chunkmap$rsid))

    missing_snps <- extract_snps(
        snps   = snps,
        indir  = indir,
        outdir = outdir,
        keep   = keep,
        idfile = idfile,
        ncore  = 2L,
        chunkmap = chunkmap,            # snp column of class factor
        chunkmap_cols = c(1L, 2L, 4L),
        quiet = TRUE)

    expect_that(is.null(missing_snps), is_true())
    stopifnot(all(file.exists(outfiles)))

    dout1 <- utils::read.table(gzfile(outfiles[1L]), stringsAsFactors = FALSE)
    dout2 <- utils::read.table(gzfile(outfiles[2L]), stringsAsFactors = FALSE)

    file.remove(outfiles)

    expect_that(dout1, is_equivalent_to(ddin1))
    expect_that(dout2, is_equivalent_to(ddin2))

    rm(missing_snps, dout1, dout2)

    ## ---------------------------------------------------------------
    ## Test 3: specifying snps as factor
    ## ---------------------------------------------------------------

    chunkmap <- utils::read.table(chunkmap_file, stringsAsFactor = FALSE, header = TRUE)

    missing_snps <- extract_snps(
        snps   = as.factor(snps),       # snps of class factor
        indir  = indir,
        outdir = outdir,
        keep   = keep,
        idfile = idfile,
        ncore  = 2L,
        chunkmap = chunkmap,
        chunkmap_cols = c(1L, 2L, 4L),
        quiet = TRUE)

    expect_that(is.null(missing_snps), is_true())
    stopifnot(all(file.exists(outfiles)))

    dout1 <- utils::read.table(gzfile(outfiles[1L]), stringsAsFactors = FALSE)
    dout2 <- utils::read.table(gzfile(outfiles[2L]), stringsAsFactors = FALSE)

    file.remove(outfiles)

    expect_that(dout1, is_equivalent_to(ddin1))
    expect_that(dout2, is_equivalent_to(ddin2))

    rm(missing_snps, dout1, dout2)

    ## ===============================================================
    ## Clean-up
    ## ===============================================================

    file.remove(infiles_gz)
})
