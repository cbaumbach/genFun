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
