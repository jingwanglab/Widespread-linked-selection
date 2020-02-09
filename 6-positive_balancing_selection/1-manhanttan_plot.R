manhattan <- function (x, chr = "CHR", bp = "BP", p = "P", snp = "SNP", col = c("gray10",
    "gray60"), chrlabs = NULL, highlight = NULL, logp = FALSE,ylabel=NULL,
    ...)
{
    CHR = POS = P = index = NULL
    if (!(chr %in% names(x)))
        stop(paste("Column", chr, "not found!"))
    if (!(bp %in% names(x)))
        stop(paste("Column", bp, "not found!"))
    if (!(p %in% names(x)))
        stop(paste("Column", p, "not found!"))
    if (!(snp %in% names(x)))
        warning(paste("No SNP column found. OK unless you're trying to highlight."))
    if (!is.numeric(x[[chr]]))
        stop(paste(chr, "column should be numeric. Do you have 'X', 'Y', 'MT', etc? If so change to numbers and try again."))
    if (!is.numeric(x[[bp]]))
        stop(paste(bp, "column should be numeric."))
    if (!is.numeric(x[[p]]))
        stop(paste(p, "column should be numeric."))
    d = data.frame(CHR = x[[chr]], BP = x[[bp]], P = x[[p]])
    if (!is.null(x[[snp]]))
        d = transform(d, SNP = x[[snp]])
    d <- subset(d, (is.numeric(CHR) & is.numeric(BP) & is.numeric(P)))
    d <- d[order(d$CHR, d$BP), ]
    if (logp) {
        d$logp <- -log10(d$P)
    }
    else {
        d$logp <- d$P
    }

    ##background colours
    bg_col = c("grey90","grey95")
    d$pos = NA
    d$index = NA
    ind = 0
    for (i in unique(d$CHR)) {
        ind = ind + 1
        d[d$CHR == i, ]$index = ind
    }
    nchr = length(unique(d$CHR))
    if (nchr == 1) {
        options(scipen = 999)
        d$pos = d$BP/1e+05
        ticks = floor(length(d$pos))/2 + 1
        xlabel = paste("Chromosome", unique(d$CHR), "position(Mb)")
        labs = ticks
    }
    else {
        lastbase = 0
        ticks = NULL
        for (i in unique(d$index)) {
            if (i == 1) {
                d[d$index == i, ]$pos = d[d$index == i, ]$BP
            }
            else {
                lastbase = lastbase + tail(subset(d, index ==
                  i - 1)$BP, 1)
                d[d$index == i, ]$pos = d[d$index == i, ]$BP +
                  lastbase
            }
            ticks = c(ticks, (min(d[d$CHR == i, ]$pos) + max(d[d$CHR ==
                i, ]$pos))/2 + 1)
        }
        xlabel = "Chromosome"
        labs <- unique(d$CHR)
    }
    print(ceiling(max(d$logp)))
    xmax = ceiling(max(d$pos) * 1.03)
    xmin = floor(max(d$pos) * -0.03)
    def_args <- list(yaxt="n", xaxt = "n", bty = "n", xaxs = "i", yaxs = "i",cex.lab=2,cex.axis=1.5,
        las = 1, pch = 20, xlim = c(xmin, xmax), ylim = c(0,
            ceiling(max(d$logp))), xlab = xlabel, ylab = ylabel)
    dotargs <- list(...)
    do.call("plot", c(NA, dotargs, def_args[!names(def_args) %in%
        names(dotargs)]))
    if (!is.null(chrlabs)) {
        if (is.character(chrlabs)) {
            if (length(chrlabs) == length(labs)) {
                labs <- chrlabs
            }
            else {
                warning("You're trying to specify chromosome labels but the number of labels != number of chromosomes.")
            }
        }
        else {
            warning("If you're trying to specify chromosome labels, chrlabs must be a character vector")
        }
    }
    
    col = rep(col, max(d$CHR))
    bg_col = rep(bg_col, max(d$CHR))
    if (nchr == 1) {
        with(d, points(pos, logp, pch = 20, col = col[1], ...))
    }
    else {
        icol = 1
        for (i in unique(d$index)) {
            with(d[d$index == unique(d$index)[i], ], rect(xleft=min(pos),xright=max(pos),ybottom=0,ytop=ceiling(max(d$logp)),col = bg_col[icol], border=NA,...))
            icol = icol + 1
        }
        icol = 1
        for (i in unique(d$index)) {
            with(d[d$index == unique(d$index)[i], ], points(pos,logp, col = col[icol], pch = 20, ...))
            icol = icol + 1
        }
    }
    if (nchr == 1) {
        axis(1, ...)
    }
    else {
        ## axis(2, lwd = 2, las = 1, ...)
        axis(2, lwd = 2, las = 1,cex.axis=1.5)
        axis(1, at = ticks, labels = labs, lwd = 2,cex.axis=1.5, ...)
        ##axis(1, at = ticks, labels = NA, lwd = 2, ...)

    }

        if (!is.null(highlight)) {
        if (any(!(highlight %in% d$SNP)))
            warning("You're trying to highlight SNPs that don't exist in your results.")
        d.highlight = d[which(d$SNP %in% highlight), ]
        with(d.highlight, points(pos, logp, col = "black", pch = 20,
            ...))
    }
}

