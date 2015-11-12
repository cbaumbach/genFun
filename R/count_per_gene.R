count_per_gene <- function(pos, chr, gene_start, gene_end, gene_chr)
{
    if (!same_length(pos, chr))
        stop("pos and chr must be of the same length")
    if (!same_length(gene_start, gene_end, gene_chr))
        stop("gene_start, gene_end, and gene_chr must be of the same length")
    if (length(pos) == 0L)
        return(rep_len(0L, length(gene_start)))
    if (length(gene_start) == 0L)
        return(integer())
    if (typeof(pos) != "integer")
        pos <- as.integer(pos)
    if (typeof(gene_start) != "integer")
        gene_start <- as.integer(gene_start)
    if (typeof(gene_end) != "integer")
        gene_end <- as.integer(gene_end)
    count <- function(idx) {
        current_chr <- gene_chr[idx][1L]
        positions <- pos[chr == current_chr]
        start <- gene_start[idx]
        end <- gene_end[idx]
        ord <- order_intervals(start, end)
        old_order <- order(ord)
        count_points_per_interval(positions,
            start[ord], end[ord])[old_order]
    }
    counts <- ave2(seq_along(gene_chr), gene_chr, count)
    counts[is.na(gene_start) | is.na(gene_end)] <- NA
    counts
}
