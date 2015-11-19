from <- function(...) {
    f(rbind(do.call("c", list(...))))
}
to <- rbind
expect <- function(from, to) {
    expect_equal(from, to)
}

# ====================================================================
context("convert_to_dosages: additive model")

f <- function(m) convert_to_dosages(m, "additive")

test_that("[100] -> [0]", expect(from(1,0,0), to(0)))
test_that("[010] -> [1]", expect(from(0,1,0), to(1)))
test_that("[001] -> [2]", expect(from(0,0,1), to(2)))

test_that("[.1 .2 .7] -> [1.6]", expect(from(.1,.2,.7), to(1.6)))

test_that("[100 100] -> [00]", {
    expect_equal(f(rbind(c(1,0,0, 1,0,0))), rbind(c(0,0)))
})

test_that("[100][100] -> [0][0]", {
    expect_equal(f(rbind(c(1,0,0),
                         c(1,0,0))),
        rbind(0, 0))
})

test_that("[100 010 001][001 010 100] -> [012][210]", {
    expect_equal(f(rbind(c(1,0,0, 0,1,0, 0,0,1),
                         c(0,0,1, 0,1,0, 1,0,0))),
        rbind(c(0,1,2),
              c(2,1,0)))
})

# ====================================================================
context("convert_to_dosages: recessive model")

f <- function(m) convert_to_dosages(m, "recessive")

test_that("[100] -> [0]", expect(from(1,0,0), to(0)))
test_that("[010] -> [0]", expect(from(0,1,0), to(0)))
test_that("[001] -> [1]", expect(from(0,0,1), to(1)))

test_that("[.1 .2 .7] -> [.7]", expect(from(.1,.2,.7), to(.7)))

# ====================================================================
context("convert_to_dosages: dominant model")

f <- function(m) convert_to_dosages(m, "dominant")

test_that("[100] -> [0]", expect(from(1,0,0), to(0)))
test_that("[010] -> [1]", expect(from(0,1,0), to(1)))
test_that("[001] -> [1]", expect(from(0,0,1), to(1)))

test_that("[.1 .2 .7] -> [.9]", expect(from(.1,.2,.7), to(.9)))

# ====================================================================
context("convert_to_dosages in terms of first allele: additive model:")

f <- function(m) convert_to_dosages(m, "additive", "A")

test_that("[100] -> [2]", expect(from(1,0,0), to(2)))
test_that("[010] -> [1]", expect(from(0,1,0), to(1)))
test_that("[001] -> [0]", expect(from(0,0,1), to(0)))

test_that("[.1 .2 .7] -> [.4]", expect(from(.1,.2,.7), to(.4)))

# ====================================================================
context("convert_to_dosages in terms of first allele: recessive model:")

f <- function(m) convert_to_dosages(m, "recessive", "A")

test_that("[100] -> [1]", expect(from(1,0,0), to(1)))
test_that("[010] -> [0]", expect(from(0,1,0), to(0)))
test_that("[001] -> [0]", expect(from(0,0,1), to(0)))

test_that("[.1 .2 .7] -> [.1]", expect(from(.1,.2,.7), to(.1)))

# ====================================================================
context("convert_to_dosages in terms of first allele: dominant model:")

f <- function(m) convert_to_dosages(m, "dominant", "A")

test_that("[100] -> [1]", expect(from(1,0,0), to(1)))
test_that("[010] -> [1]", expect(from(0,1,0), to(1)))
test_that("[001] -> [0]", expect(from(0,0,1), to(0)))

test_that("[.1 .2 .7] -> [.3]", expect(from(.1,.2,.7), to(.3)))

# ====================================================================
context("convert_to_dosages: treatment of NAs")

f <- function(m) convert_to_dosages(m, "additive")

test_that("[NA NA NA] -> [NA]", expect(from(NA,NA,NA), to(as.numeric(NA))))
test_that("[NA**] -> [NA]", expect(from(NA,0,0), to(as.numeric(NA))))
test_that("[NA** 001] -> [NA 2]", expect(from(NA,0,0, 0,0,1), to(c(NA, 2))))

test_that("[NA**][001] -> [NA][2]", {
    expect_equal(f(rbind(c(NA,0,0), c(0,0,1))), rbind(NA, 2))
})
