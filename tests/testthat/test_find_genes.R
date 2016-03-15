context("Finding genes containing given SNPs")

test_that("everything works fine in the vanilla use case", {

    d <- utils::read.table(text = "
    snp chr pos
    rs1 1 1
    rs2 1 5
    rs3 2 3
    rs4 3 4
    ", header = TRUE, stringsAsFactors = FALSE)

    genes <- utils::read.table(text = "
    id chr start end
    g1 1 1 3
    g2 1 1 6
    g3 2 1 4
    g4 4 1 5
    ", header = TRUE, stringsAsFactors = FALSE)

    dout <- utils::read.table(text = "
    snp chr pos genes
    rs1 1 1 g1,g2
    rs2 1 5 g2
    rs3 2 3 g3
    rs4 3 4 NA
    ", header = TRUE, stringsAsFactors = FALSE)

    expect_that(find_genes(d, genes, quiet = TRUE), equals(dout))
})

test_that("missing values in snp data frame are propagated", {

    d <- utils::read.table(text = "
    snp chr pos
    rs1 NA 1
    rs2 1 5
    rs3 2 NA
    rs4 3 4
    ", header = TRUE, stringsAsFactors = FALSE)

    genes <- utils::read.table(text = "
    id chr start end
    g1 1 1 3
    g2 1 1 6
    g3 2 1 4
    g4 4 1 5
    ", header = TRUE, stringsAsFactors = FALSE)

    dout <- utils::read.table(text = "
    snp chr pos genes
    rs1 NA 1 NA
    rs2 1 5 g2
    rs3 2 NA NA
    rs4 3 4 NA
    ", header = TRUE, stringsAsFactors = FALSE)

    expect_that(find_genes(d, genes, quiet = TRUE), equals(dout))
})

test_that("NAs in genes data frame are ignored", {

    d <- utils::read.table(text = "
    snp chr pos
    rs1 1 1
    rs2 1 5
    rs3 2 3
    rs4 3 4
    ", header = TRUE, stringsAsFactors = FALSE)

    genes <- utils::read.table(text = "
    id chr start end
    g1 NA 1 3
    g2 1 NA 6
    g3 2 1 4
    g4 4 1 NA
    ", header = TRUE, stringsAsFactors = FALSE)

    dout <- utils::read.table(text = "
    snp chr pos genes
    rs1 1 1 NA
    rs2 1 5 NA
    rs3 2 3 g3
    rs4 3 4 NA
    ", header = TRUE, stringsAsFactors = FALSE)

    expect_that(find_genes(d, genes, quiet = TRUE), equals(dout))
})

test_that("missing columns throw an error", {

    d <- utils::read.table(text = "
    snp chr pos
    rs1 1 1
    rs2 1 5
    rs3 2 3
    rs4 3 4
    ", header = TRUE, stringsAsFactors = FALSE)

    genes <- utils::read.table(text = "
    id chr start end
    g1 1 1 3
    g2 1 1 6
    g3 2 1 4
    g4 4 1 5
    ", header = TRUE, stringsAsFactors = FALSE)

    d_chr1 <- utils::read.table(text = "
    snp x pos
    rs1 1 1
    rs2 1 5
    rs3 2 3
    rs4 3 4
    ", header = TRUE, stringsAsFactors = FALSE)

    d_pos <- utils::read.table(text = "
    snp chr x
    rs1 1 1
    rs2 1 5
    rs3 2 3
    rs4 3 4
    ", header = TRUE, stringsAsFactors = FALSE)

    genes_chr <- utils::read.table(text = "
    id x start end
    g1 1 1 3
    g2 1 1 6
    g3 2 1 4
    g4 4 1 5
    ", header = TRUE, stringsAsFactors = FALSE)

    genes_start <- utils::read.table(text = "
    id chr x end
    g1 1 1 3
    g2 1 1 6
    g3 2 1 4
    g4 4 1 5
    ", header = TRUE, stringsAsFactors = FALSE)

    genes_end <- utils::read.table(text = "
    id chr start x
    g1 1 1 3
    g2 1 1 6
    g3 2 1 4
    g4 4 1 5
    ", header = TRUE, stringsAsFactors = FALSE)

    genes_id <- utils::read.table(text = "
    x chr start end
    g1 1 1 3
    g2 1 1 6
    g3 2 1 4
    g4 4 1 5
    ", header = TRUE, stringsAsFactors = FALSE)

    expect_that(find_genes(d_chr1, genes, quiet = TRUE),  throws_error("column missing in d: chr"))
    expect_that(find_genes(d_pos,  genes, quiet = TRUE),  throws_error("column missing in d: pos"))
    expect_that(find_genes(d, genes_chr, quiet = TRUE),   throws_error("column missing in genes: chr"))
    expect_that(find_genes(d, genes_start, quiet = TRUE), throws_error("column missing in genes: start"))
    expect_that(find_genes(d, genes_end, quiet = TRUE),   throws_error("column missing in genes: end"))
    expect_that(find_genes(d, genes_id, quiet = TRUE),    throws_error("column missing in genes: id"))
})
