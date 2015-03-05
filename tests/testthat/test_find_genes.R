context("Finding genes containing given SNPs")

test_that("everything works fine in the vanilla use case", {

    d <- read.table(textConnection("\
snp chr pos
rs1 1 1
rs2 1 5
rs3 2 3
rs4 3 4
", "r"), header = TRUE, stringsAsFactors = FALSE)

    genes <- read.table(textConnection("\
id chr start end
g1 1 1 3
g2 1 1 6
g3 2 1 4
g4 4 1 5
", "r"), header = TRUE, stringsAsFactors = FALSE)

    dout <- read.table(textConnection("\
snp chr pos genes
rs1 1 1 g1,g2
rs2 1 5 g2
rs3 2 3 g3
rs4 3 4 NA
", "r"), header = TRUE, stringsAsFactors = FALSE)

    expect_that(find_genes(d, genes), equals(dout))
})

test_that("missing values in snp data frame are propagated", {

    d <- read.table(textConnection("\
snp chr pos
rs1 NA 1
rs2 1 5
rs3 2 NA
rs4 3 4
", "r"), header = TRUE, stringsAsFactors = FALSE)

    genes <- read.table(textConnection("\
id chr start end
g1 1 1 3
g2 1 1 6
g3 2 1 4
g4 4 1 5
", "r"), header = TRUE, stringsAsFactors = FALSE)

    dout <- read.table(textConnection("\
snp chr pos genes
rs1 NA 1 NA
rs2 1 5 g2
rs3 2 NA NA
rs4 3 4 NA
", "r"), header = TRUE, stringsAsFactors = FALSE)

    expect_that(find_genes(d, genes), equals(dout))
})

test_that("NAs in genes data frame are ignored", {

    d <- read.table(textConnection("\
snp chr pos
rs1 1 1
rs2 1 5
rs3 2 3
rs4 3 4
", "r"), header = TRUE, stringsAsFactors = FALSE)

    genes <- read.table(textConnection("\
id chr start end
g1 NA 1 3
g2 1 NA 6
g3 2 1 4
g4 4 1 NA
", "r"), header = TRUE, stringsAsFactors = FALSE)

    dout <- read.table(textConnection("\
snp chr pos genes
rs1 1 1 NA
rs2 1 5 NA
rs3 2 3 g3
rs4 3 4 NA
", "r"), header = TRUE, stringsAsFactors = FALSE)

    expect_that(find_genes(d, genes), equals(dout))
})

test_that("missing columns throw an error", {

    d <- read.table(textConnection("\
snp chr pos
rs1 1 1
rs2 1 5
rs3 2 3
rs4 3 4
", "r"), header = TRUE, stringsAsFactors = FALSE)

    genes <- read.table(textConnection("\
id chr start end
g1 1 1 3
g2 1 1 6
g3 2 1 4
g4 4 1 5
", "r"), header = TRUE, stringsAsFactors = FALSE)

    d_chr1 <- read.table(textConnection("\
snp x pos
rs1 1 1
rs2 1 5
rs3 2 3
rs4 3 4
", "r"), header = TRUE, stringsAsFactors = FALSE)

    d_pos <- read.table(textConnection("\
snp chr x
rs1 1 1
rs2 1 5
rs3 2 3
rs4 3 4
", "r"), header = TRUE, stringsAsFactors = FALSE)

    genes_chr <- read.table(textConnection("\
id x start end
g1 1 1 3
g2 1 1 6
g3 2 1 4
g4 4 1 5
", "r"), header = TRUE, stringsAsFactors = FALSE)

    genes_start <- read.table(textConnection("\
id chr x end
g1 1 1 3
g2 1 1 6
g3 2 1 4
g4 4 1 5
", "r"), header = TRUE, stringsAsFactors = FALSE)

    genes_end <- read.table(textConnection("\
id chr start x
g1 1 1 3
g2 1 1 6
g3 2 1 4
g4 4 1 5
", "r"), header = TRUE, stringsAsFactors = FALSE)

    genes_id <- read.table(textConnection("\
x chr start end
g1 1 1 3
g2 1 1 6
g3 2 1 4
g4 4 1 5
", "r"), header = TRUE, stringsAsFactors = FALSE)

    expect_that(find_genes(d_chr1, genes),  throws_error("column missing in d: chr"))
    expect_that(find_genes(d_pos,  genes),  throws_error("column missing in d: pos"))
    expect_that(find_genes(d, genes_chr),   throws_error("column missing in genes: chr"))
    expect_that(find_genes(d, genes_start), throws_error("column missing in genes: start"))
    expect_that(find_genes(d, genes_end),   throws_error("column missing in genes: end"))
    expect_that(find_genes(d, genes_id),    throws_error("column missing in genes: id"))
})
