convert_to_dosages <- function(aMatrix,
    model = c("additive", "recessive", "dominant"),
    in_terms_of_allele = c("B", "A"),
    transpose = FALSE)
{
    weights <- switch(EXPR = model[1L],
        additive  = c(0L, 1L, 2L),
        recessive = c(0L, 0L, 1L),
        dominant  = c(0L, 1L, 1L))
    if (in_terms_of_allele[1L] == "A")
        weights <- rev(weights)
    row_to_dosages <- function(aRow) {
        y <- aRow * weights
        X <- TRUE; O <- FALSE
        y[c(X,O,O)] + y[c(O,X,O)] + y[c(O,O,X)]
    }
    dosages <- apply(aMatrix, 1L, row_to_dosages)
    if (is.vector(dosages))
        dosages <- t(as.matrix(dosages))
    if (transpose)
        dosages
    else
        t(dosages)
}
