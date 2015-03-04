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
snp chr pos id
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
snp chr pos id
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
snp chr pos id
rs1 1 1 NA
rs2 1 5 NA
rs3 2 3 g3
rs4 3 4 NA
", "r"), header = TRUE, stringsAsFactors = FALSE)

    expect_that(find_genes(d, genes), equals(dout))
})

