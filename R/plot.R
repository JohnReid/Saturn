#' Plot a motif scan
#'
plot.motif.scan <- function(motif.scan, range.to.plot) {
    # Subset hits that overlap the range to plot
    overlaps <- Filter(length, lapply(motif.scan, function(ms) subsetByOverlaps(ms, range.to.plot)))
    # Merge hits adding motif name as meta data
    motifs.merged <- do.call(c, lapply(names(overlaps),
                                       function(n) {
                                           mcols(overlaps[[n]])$motif <- n
                                           overlaps[[n]]
                                       }))
    if (length(motifs.merged)) {
      # Plot motifs
      autoplot(motifs.merged, aes(fill = logBF), facets = motif ~ .)
    } else {
      return(NULL)
    }
}


#' Plot ranges
#'
plot.ranges <- function(range.to.plot, cell, motif.scans) {
  heights <- NULL
  plot.args <- list()
  #
  # Load DNase
  suppressWarnings(
    dnase.wiggle <- rtracklayer::import(dnase.signal.file(cell),
                                        selection = range.to.plot))
  plot.args$DNase <- autoplot(dnase.wiggle, geom = "step", aes(y = score))
  heights <- append(heights, 2)
  #
  # Load DGF
  dgf <- load.wellington(cell.valid)
  dgf.range <- subsetByOverlaps(dgf, range.to.plot)
  if (length(dgf.range)) {
    heights <- append(heights, 1)
    plot.args$Wellington <- autoplot(dgf.range, aes(fill = score))
  }
  #
  # Motifs
  motif.plots <- lapply(motif.scans, function(ms) plot.motif.scan(ms, range.to.plot))
  names(motif.plots) <- names(motif.scans)
  motif.plots <- Filter(length, motif.plots)
  plot.args <- append(plot.args, motif.plots)
  plot.args$heights <- c(heights, rep(7, length(motif.plots)))
  # Plot
  do.call(tracks, plot.args) +
    theme_few() +
    theme(strip.text.y = element_text(angle=0))
}
