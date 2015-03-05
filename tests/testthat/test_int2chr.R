context("Integer to chromosome conversion")

test_that("valid integers map to the correct chromosomes", {
    expect_that(int2chr(c(1L,22:25)), equals(c("1","22","X","Y","MT")))
})

test_that("invalid integers map to NA", {
    expect_that(int2chr(c(-1L,99L,NA_integer_)), equals(rep(NA_character_, 3L)))
})

test_that("non-integer input maps to NA", {
    expect_that(int2chr(-1),       throws_error("argument must be of type 'integer'"))
    expect_that(int2chr(99),       throws_error("argument must be of type 'integer'"))
    expect_that(int2chr(NA_real_), throws_error("argument must be of type 'integer'"))
})
