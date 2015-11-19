context("count_per_gene")

# Notation used in test cases:
#     points:  chr|pos
#     genes:   chr[start,end]

f <- count_per_gene

test_that("gene_start, gene_end, and gene_chr must be of the same length", {
    msg <- "gene_start, gene_end, and gene_chr must be of the same length"
    expect_error(f(0, 0, 0, 0:1, 0:1), msg)
})

test_that("pos and chr must be of the same length", {
    msg <- "pos and chr must be of the same length"
    expect_error(f(0, 0:1, 0, 0, 0), msg)
})

test_that("0|0 -> 0[0,0] -> 1", expect_equal(f(0, 0, 0, 0, 0), 1))
test_that("0|0 -> 0[1,1] -> 0", expect_equal(f(0, 0, 1, 1, 0), 0))
test_that("0|0 -> 0[NA,0] -> NA", expect_equal(f(0, 0, NA, 0, 0), NA_integer_))
test_that("0|0 -> 0[0,NA] -> NA", expect_equal(f(0, 0, 0, NA, 0), NA_integer_))
test_that("0|0 -> 0[NA,NA] -> NA", expect_equal(f(0, 0, NA, NA, 0), NA_integer_))
test_that("0|0 -> 1[0,0] -> 0", expect_equal(f(0, 0, 0, 0, 1), 0))
test_that("0|0 -> 0['0',0] -> 1", expect_equal(f(0, 0, "0", 0, 0), 1))
test_that("0|0 -> 0[0,'0'] -> 1", expect_equal(f(0, 0, 0, "0", 0), 1))
test_that("0|0 -> 0['0','0'] -> 1", expect_equal(f(0, 0, "0", "0", 0), 1))
test_that("0|'0' -> 0[0,0] -> 1", expect_equal(f("0", 0, 0, 0, 0), 1))
test_that("0|0 -> 'foo'[0,0] -> 0", expect_equal(f(0, 0, 0, 0, "foo"), 0))

test_that("0|0 -> 0[0,0] 0[1,1] 0[0,0] -> 1 0 1", {
    expect_equal(f(0, 0, c(0,1,0), c(0,1,0), c(0,0,0)), c(1,0,1))
})

test_that("0|0 1|0 -> 1[0,0] 0[0,0] 1[1,1] 0[1,1] 1[0,0] 0[0,0] -> 1 1 0 0 1 1", {
    expect_equal(f(pos = c(0,0),
                   chr = c(0,1),
                   gene_start = c(0,0,1,1,0,0),
                   gene_end   = c(0,0,1,1,0,0),
                   gene_chr   = c(1,0,1,0,1,0)),
        c(1,1,0,0,1,1))
})

test_that("NULL|NULL -> 0[0,0] -> 0", {
    expect_equal(f(NULL, NULL, 0, 0, 0), 0)
})

test_that("0|0 -> NULL[NULL,NULL] -> empty", {
    expect_equal(f(0, 0, NULL, NULL, NULL), integer())
})

test_that("NULL|NULL -> NULL[NULL,NULL] -> empty", {
    expect_equal(f(NULL, NULL, NULL, NULL, NULL), integer())
})

test_that("empty|empty -> 0[0,0] -> 0", {
    expect_equal(f(integer(), integer(), 0, 0, 0), 0)
})

test_that("0|0 -> empty[empty,empty] -> empty", {
    expect_equal(f(0, 0, integer(), integer(), integer()), integer())
})

test_that("empty|empty -> empty[empty,empty] -> empty", {
    expect_equal(f(integer(), integer(), integer(), integer(), integer()), integer())
})
