context("Chromosome to integer conversion")

test_that("valid chromosomes map to correct integers", {
    x <- c("1", "22", "X", "Y", "MT")
    chrx <- paste0("chr", x)

    expect_that(chr2int(x),                    equals(c(1, 22:25)))
    expect_that(chr2int(chrx, prefix = "chr"), equals(c(1, 22:25)))
})

test_that("invalid chromosomes map to NA", {
    x_char <- c("-1", "99", "A", "B", NA)
    x_int  <- c(-1L, 99L, NA_integer_)
    x_num  <- c(-1, 99, NA_real_)
    chrx   <- paste0("chr", x_char)

    expect_that(chr2int(x_char),               equals(rep(NA_integer_, 5L)))
    expect_that(chr2int(chrx, prefix = "chr"), equals(rep(NA_integer_, 5L)))
    expect_that(chr2int(chrx, prefix = "foo"), equals(rep(NA_integer_, 5L)))
    expect_that(chr2int(x_int),                equals(rep(NA_integer_, 3L)))
    expect_that(chr2int(x_num),                equals(rep(NA_integer_, 3L)))
})

test_that("conversion is case insensitive", {
    x <- c("1", "X", "Y", "MT")

    test_that(chr2int(toupper(x)),                equals(c(1L,23:25)))
    test_that(chr2int(tolower(x)),                equals(c(1L,23:25)))
    test_that(chr2int(paste0("chr", toupper(x))), equals(c(1L,23:25)))
    test_that(chr2int(paste0("chr", tolower(x))), equals(c(1L,23:25)))
})
