mplot=function (pval, bp, chr, groups = NULL, cutoff = NULL, xlab = "Chromosome (base-pair position)", 
    ylab = expression(paste(-log[10](italic(p)))), transform = TRUE, 
    cex = 0.5, ...) 
{
    lu <- function(...) length(unique(...))
    if (is.character(chr)) {
        chr <- as.integer(swap(chr, c("x", "X", "y", "Y", "M"), 
            c(23, 23, 24, 24, 25)))
    }
    dat <- data.frame(P = pval, BP = bp, CHR = chr)
    if (!is.null(groups)) {
        dat$GROUP <- factor_(groups)
    }
    if (any(!complete.cases(dat))) {
        warning("There are NAs in your data; these points are removed for the plot")
        dat <- dat[complete.cases(dat), ]
    }
    dat <- dat[order(dat$CHR), ]
    reptimes <- counts(dat$CHR)
    maxes <- with(dat, tapply(BP, CHR, max))
    maxes.cumsum <- cumsum(as.numeric(maxes))
    tmp <- c(0, maxes.cumsum)
    diffs <- rep(0, length(maxes.cumsum))
    for (i in 1:(length(tmp) - 1)) {
        diffs[i] <- (tmp[i + 1] - tmp[i])/2 + tmp[i]
    }
    nCHR <- lu(dat$CHR)
    adj <- rep(c(0, maxes.cumsum[-nCHR]), times = reptimes)
    dat$relBP <- dat$BP + adj
    cols <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", 
        "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A")
    cols <- cols[c(1, 2, 7, 8, 9, 10, 3, 4, 5, 6, 11, 12)]
    if (is.null(groups)) {
        cols <- c("grey50", "grey15")
    }
    else {
        cols <- c(cols[seq(1, length.out = lu(dat$GROUP), by = 2)], 
            cols[seq(2, length.out = lu(dat$GROUP), by = 2)])
    }
    dat$COL_GROUP <- paste(dat$CHR%%2 == 0, dat$GROUP, sep = "_")
    dat$COL <- swap(dat$COL_GROUP, unique(dat$COL_GROUP), cols)
    axis_tick_labels <- 1:nCHR
    axis_tick_labels[axis_tick_labels > 12 & axis_tick_labels%%2 == 
        1] <- ""
    if (is.null(groups)) {
        myKey <- NULL
    }
    else {
        myKey <- list(text = list(as.character(unique(dat$GROUP))), 
            points = list(col = cols[1:lu(dat$GROUP)], pch = 20), 
            points = list(col = cols[(lu(dat$GROUP) + 1):(2 * 
                lu(dat$GROUP))], pch = 20))
    }
    if (transform) {
        dat$P <- -log10(dat$P)
    }
    print(with(dat, xyplot(P ~ relBP, pch = 20, col = COL, cex = cex, 
        scales = list(x = list(relation = "same", tck = c(1, 
            0), at = diffs, labels = axis_tick_labels)), panel = function(x, 
            y, ...) {
            panel.xyplot(x, y, ...)
            panel.abline(h = -log10(0.05), col = "red", 
                lty = "dashed")
            for (i in maxes.cumsum) {
                panel.segments(i, 0, i, 1000, col = "grey85", 
                  lty = "dashed")
            }
        }, key = myKey, xlim = c(0, max(dat$relBP)), ylim = c(0, 
            100), origin = 0, xlab = xlab, ylab = ylab, 
        ...)))
}