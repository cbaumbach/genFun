## genFun.R

find_genes <- function(d, genes, chr1 = "chr", pos = "pos", out = id,
                       chr2 = chr1, start = "start", end = "end",
                       id = "id")
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
    genes <- subset(genes, !any_missings)

    d[[out]] <- NA
    for (k in sort(unique(d[[chr1]]))) {
        idx <- d[[chr1]] == k
        idx[is.na(idx)] <- FALSE
        gidx <- genes[[chr2]] == k
        if (!any(gidx)) next
        d[[out]][idx] <- match_intervals(
            d[[pos]][idx],
            genes[[start]][gidx],
            genes[[end]][gidx],
            genes[[id]][gidx])
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
