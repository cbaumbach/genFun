read.obo <- function(filename, colClasses = NULL)
{
    ## Find path to perl script for parsing obo file.
    package_dir <- find.package("genFun")
    perl_script <- file.path(package_dir, "perl", "obo2table.pl")

    cmd <- paste(perl_script, filename)

    if (is.null(colClasses))
        colClasses <- "character"

    con <- pipe(cmd, "r")
    on.exit(close(con))

    read.delim(con, colClasses = colClasses)
}

dot_attributes <- function(...)
{
    x <- list(...)
    n <- names(x)
    x <- as.character(x)
    names(x) <- n
    class(x) <- c("dot_attributes", class(x))
    x
}

print.dot_attributes <- function(x, ...)
{
    pr1("[")
    y <- sapply(names(x), function(tag)
        sprintf("%s = %s", tag, double_quote(x[tag])))
    pr1(paste(y, collapse = ", "))
    pr1("]")
}
