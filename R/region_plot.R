arrange_regions <- function(label, start, end = start)
{
    ## Check whether there is an active graphics device.
    if (dev.cur() == 1L)
        stop("Need active graphics device other than null device.")

    ## Sort regions.
    d <- data.frame(label = label,
                    start = start,
                    end   = end,
                    stringsAsFactors = FALSE)
    d <- d[order(d$start, d$end), , drop = FALSE]

    d$level <- NA

    ## For every level used so far we save the rightmost position.
    pos_on_level <- NULL

    for (i in seq_len(nrow(d))) {

        ## Compute effective start and end positions taking into
        ## account the width of the label.
        start <- d$start[i]
        end   <- d$end[i]
        mid   <- (start + end) / 2
        len   <- strwidth(d$label[i])
        real_start <- min(mid - len/2, start)
        real_end   <- max(mid + len/2, end)

        ## Find the first free level.
        if (is.null(pos_on_level)) {
            pos_on_level <- real_end
            level <- 1L
        }
        else {
            free <- pos_on_level < real_start
            if (!any(free)) {
                ## In order to keep the code as simple as possible we
                ## want pos_on_level to be sorted by level, with level
                ## 1 coming first.  In return we have to append new
                ## levels to the end of pos_on_level.  Inefficient,
                ## but since the number of labels is not going to be
                ## very high we won't notice.
                pos_on_level <- c(pos_on_level, real_end)
                level <- length(pos_on_level)
            }
            else {
                level <- which.max(free)
                pos_on_level[level] <- real_end
            }
        }
        d$level[i] <- level
    }

    d
}

## +--------------------------------------+
## |       --  --  --  --  --  --         |
## |    +---|---|---|---|---|---|----+    |  upper part
## | -- |   .           . .          |    |
## | -- |   ..          ...     .    |    | (fixed size)
## | -- |  ...   .      ....    . .  |    |
## |----+----------------------------+----| - - - - - -
## |                                      |
## |  LABEL1    LABEL3         LABEL6     |
## |   ----  -------------     ------     |   lower part
## |     LABEL2       LABEL4     LABEL7   |
## |      ----     -----------    -----   | (variable size)
## |                     LABEL5           |
## |                   -----------        |
## +--------------------------------------+

region_plot <- function(file, x, y, start, end, label, width,
                        ymax = NULL, transform = NULL, col = NULL,
                        main = NULL)
{
    ## ===============================================================
    ## Argument checking.
    ## ===============================================================

    ## Check for existence of mandatory arguments.
    force(file)
    force(x)
    force(y)
    force(start)
    force(end)
    force(label)
    force(width)

    if (!all_neighbors(`==`, vapply(list(start, end, label), length,
                                    integer(1L))))
        stop("`start', `end', and `label' must have same length.")

    plot_regions <- length(label) > 0L

    ## ===============================================================
    ## Graphical parameters.
    ## ===============================================================
    label_cex <- .7                # size of region annotation labels
    point_cex <- .3                # size of points in plotting region
    xaxis_cex <- label_cex         # size of x-axis labels
    yaxis_cex <- label_cex         # size of y-axis labels

    ## ===============================================================
    ## Determine dimensions of upper part.
    ## ===============================================================
    if (!is.null(transform))
        y <- transform(y)

    if (is.null(ymax))
        ymax <- max(y[is.finite(y)])

    y[is.infinite(y)] <- ymax

    ## Increase x-range by 1 percent on either side.
    xlim <- range(x) + .01 * diff(range(x)) * c(-1, 1)
    ylim <- c(0, ymax)

    ## ===============================================================
    ## Determine number of lines of region annotation.
    ## ===============================================================
    if (plot_regions) {
        tryCatch({
            pdf("/dev/null", width = width)
            plot(xlim, 0:1, type = "n")
            par(cex = label_cex)
            region <- arrange_regions(label = label, start = start, end = end)
        }, finally = dev.off())
    }
    nlevel <- if (plot_regions) max(region$level) else 0L

    ## ===============================================================
    ## Compute device dimensions.
    ## ===============================================================
    label_ht <- .2         # height of label annotation line in inches
    upper_ht <- 3          # height of upper device region in inches
    device_ht <- upper_ht + label_ht * nlevel # total height
    upper_frac <- upper_ht / device_ht # relative height of upper part

    ## ===============================================================
    ## Create graph.
    ## ===============================================================
    pdf(file, width, device_ht)
    on.exit(dev.off())
    if (plot_regions) {
        par(
            mai = c(0, .5, 1, .5),      # no bottom margin
            xaxs = "i"                  # inner x-axis
        )
        layout(matrix(1:2), heights = c(upper_frac, 1 - upper_frac))
    }
    else {
        par(
            mai = c(.1, .5, 1, .5),     # small no bottom margin
            xaxs = "i"                  # inner x-axis
        )
    }

    ## xy-plot in upper region.
    the_col <- rep_len(if (is.null(col)) "black" else col, length(x))
    plot(x, y, xlim = xlim, ylim = ylim, cex = point_cex,
        axes = FALSE, col = the_col)
    axis(2, las = 1, hadj = 1, cex.axis = yaxis_cex)
    xat <- pretty(xlim, n = min(length(x), 10L))
    if (length(xat) >= 3L)
        xat <- xat[-c(1, length(xat))]  # drop first and last value
    axis(3, at = xat, cex.axis = xaxis_cex,
         labels = format(xat, scientific = FALSE, big.mark = ","))
    box()
    if (!is.null(main))
        title(main = main, line = 3)

    if (plot_regions) {
        ## Plot horizontal segments representing regions.
        par(
            mai = c(0, .5, 0, .5),      # no bottom/top margin
            xaxs = "i",                 # inner x-axis
            yaxs = "i"                  # inner y-axis
        )
        plot(xlim[1L], 0, xlim = xlim, ylim = c(nlevel+.5, -.5),
             ann = FALSE, axes = FALSE, type = "n")
        segments(x0 = region$start, y0 = region$level, x1 = region$end,
                 xpd = NA, lwd = 2)

        ## Add labels.
        text(x      = (region$start + region$end) / 2,
             y      = region$level,
             labels = region$label,
             pos    = 3,
             cex    = label_cex,
             offset = .15,
             xpd    = TRUE)
    }
}
