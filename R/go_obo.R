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

read.gene_ontology <- function(filename, ontology = c(
                                         "biological_process",
                                         "cellular_component",
                                         "molecular_function"))
{
    d <- read.obo(filename, colClasses = "character")

    ## Whenever a relationship `R' holds between two nodes `a' and
    ## `b', i.e., whenever we have `aRb', where `R' is one of `is_a',
    ## `part_of', `regulates', `negatively_regulates', or
    ## `positively_regulates', we will interpret this to mean `a is a
    ## child of b', or equivalently, `b is a parent of a'.
    relationships <- c("is_a", "part_of", "regulates",
                       "negatively_regulates",
                       "positively_regulates")
    d$parents <- apply(d[, relationships], 1L,
                       function(x) paste(x[!is.na(x)], collapse = ","))

    ## Create a separate data frame for each ontology.
    ds <- lapply(ontology, function(x) d[d$namespace == x, ])
    names(ds) <- ontology

    ## Convert data frames into trees.
    trees <- lapply(ds, make_tree)

    ## Make trees instances of class "gene_ontology".
    lapply(trees,
           function(tr)
           {
               class(tr) <- c("gene_ontology", class(tr))
               tr
           })
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

go_node <- function(id, data, attrib)
{
    name <- data$name[data$id == id]
    n <- wrap_lines(name, 15, max_lines = 2L, sep = "\\n")

    print(dot_attributes(label = n))
}

go_edge <- function(from, to, data, attrib)
{
    x <- data[data$id == to, ]

    d <- list()
    for (i in c("is_a", "part_of", "regulates",
                "positively_regulates",
                "negatively_regulates")) {
        d[[i]] <- strsplit(x[, i], ",", fixed = TRUE)[[1L]]
    }

    ## IS_A
   if (from %in% d[["is_a"]]) {

        print(dot_attributes(
            style = "solid",
            color = "black",
            dir   = "back"))
    }
    ## PART_OF
    else if (from %in% d[["part_of"]]) {

        print(dot_attributes(
            style = "dashed",
            color = "black",
            dir   = "back"))
    }
    ## REGULATES
    else if (from %in% d[["regulates"]]) {

        print(dot_attributes(
            style = "dotted",
            color = "gray",
            dir   = "back"))
    }
    ## NEGATIVELY_REGULATES
    else if (from %in% d[["negatively_regulates"]]) {

        print(dot_attributes(
            style = "dotted",
            color = "red",
            dir   = "back"))
    }
    ## POSITIVELY_REGULATES
    else if (from %in% d[["positively_regulates"]]) {

        print(dot_attributes(
            style = "dotted",
            color = "green",
            dir   = "back"))
    }
}

print.gene_ontology <- function(x, ...)
{
    NextMethod("print", x, nodef = go_node, edgef = go_edge)
}
