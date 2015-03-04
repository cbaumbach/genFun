## genFun.R

find_genes <- function(d, genes, chr1 = "chr", pos = "pos", out = id,
                       chr2 = chr1, start = "start", end = "end",
                       id = "id")
{
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
