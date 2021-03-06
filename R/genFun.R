## genFun.R

find_genes <- function(d, genes, chr1 = "chr", pos = "pos", out = "genes",
                       chr2 = chr1, start = "start", end = "end",
                       id = "id", quiet = FALSE)
{
    ## Check whether all required column are present.
    names1 <- names(d)
    names2 <- names(genes)

    if (! chr1 %in% names1)
        stop("column missing in d: ", chr1)
    if (! pos %in% names1)
        stop("column missing in d: ", pos)
    if (! chr2 %in% names2)
        stop("column missing in genes: ", chr2)
    if (! start %in% names2)
        stop("column missing in genes: ", start)
    if (! end %in% names2)
        stop("column missing in genes: ", end)
    if (! id %in% names2)
        stop("column missing in genes: ", id)

    ## Drop rows with missings from genes.
    any_missings <- is.na(genes[[chr2]])  |
                    is.na(genes[[start]]) |
                    is.na(genes[[end]])   |
                    is.na(genes[[id]])
    genes <- genes[!any_missings, , drop = FALSE]

    d[[out]] <- NA
    for (k in sort(unique(d[[chr1]]))) {
        if (!quiet)
            pr("Chromosome ", k)
        idx <- d[[chr1]] == k
        idx[is.na(idx)] <- FALSE
        gidx <- genes[[chr2]] == k
        if (!any(gidx)) next
        d[[out]][idx] <- match_intervals(
            d[[pos]][idx],
            genes[[start]][gidx],
            genes[[end]][gidx],
            genes[[id]][gidx],
            quiet = quiet)
    }
    d[[out]][d[[out]] == ""] <- NA
    d
}

chr2int <- function(x, prefix = NULL)
{
    crop_to_range <- function(x)
    {
        x[x < 1L | 25 < x] <- NA
        x
    }
    if (typeof(x) == "integer")
        return(crop_to_range(x))
    if (!is.null(prefix))
        x <- sub(prefix, "", x)
    x <- toupper(x)
    x[x == "X"]  <- "23"
    x[x == "Y"]  <- "24"
    x[x == "MT"] <- "25"
    x[!grepl("^[1-9][0-9]?$", x)] <- NA
    crop_to_range(as.integer(x))
}

int2chr <- function(x)
{
    if (typeof(x) != "integer")
        stop("argument must be of type 'integer'")
    x[x < 1L | 25 < x] <- NA
    x[x == 23L] <- "X"
    x[x == 24L] <- "Y"
    x[x == 25L] <- "MT"
    as.character(x)
}

chr <- function(xs, pattern = ".*chr([^_]+)")
{
    if (nsubexp(pattern) != 1L)
        stop("`pattern' must contain exactly one parenthesized ",
             "subexpression: ", single_quote(pattern))

    maybe(as.integer, submatch(pattern, xs, drop = TRUE))
}

single_quote <- function(x) {
    paste0("'", x, "'")
}

chunk <- function(xs, pattern = ".*chr[^_]+_(\\d+)")
{
    if (nsubexp(pattern) != 1L)
        stop("`pattern' must contain exactly one parenthesized ",
             "subexpression: ", single_quote(pattern))

    maybe(as.integer, submatch(pattern, xs, drop = TRUE))
}

unslash <- function(dirs)
{
    sub("/+$", "", dirs)
}
