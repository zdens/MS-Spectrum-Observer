library(MALDIquantForeign)
library(data.table)
library(ggsci)

mirror <- function(x, direction = 'left') {
  if (direction == 'right') {
    x <- rev(x)
  }
  .ps <- x[[1]]
  mirrored <- 2 * .ps - x[-1]
  if (direction == 'left') {
    mirrored <- rev(mirrored)
  }
  return(mirrored)
}
remove_negative_snr <- function(peaks,
                                spectra_aligned,
                                zeros_num = 100,
                                halfWindowSize = 2e-4,
                                snr = 1) {
  neg_snr_peaks_idx <- which(sapply(peaks, function(p) {
    length(which(p@snr < 0))
  }) > 0)
  if (length(neg_snr_peaks_idx) > 0) {
    peaks[neg_snr_peaks_idx] <-
      lapply(neg_snr_peaks_idx, function(neg_idx) {
        .spectra <- spectra_aligned[[neg_idx]]
        .mass <- c(
          mirror(.spectra@mass[1:zeros_num]),
          .spectra@mass,
          mirror(.spectra@mass[(length(.spectra@mass) - zeros_num + 1):length(.spectra@mass)], direction = 'right')
        )
        .intensity <- c(rep(0, zeros_num - 1),
                        .spectra@intensity,
                        rep(0, zeros_num - 1))
        .spectra <-
          MALDIquant::createMassSpectrum(
            mass = .mass,
            intensity = .intensity,
            metaData = .spectra@metaData
          )
        detectPeaks(
          .spectra,
          SNR = snr,
          halfWindowSize = halfWindowSize,
          method = "SuperSmoother"
        )
      })
  }
  return(peaks)
}

scientific_10 <- function(x) {
  #browser()
  cat("X: ", x, "\n")
  y <- gsub("e-0", "e-", format(x, scientific = TRUE, trim = TRUE))
  y <- gsub("e\\+0", "e", y)
  text_labs <- gsub("e", " %*% 10^\"", paste0(y, "\""))
  zeros <- stri_startswith_fixed(text_labs, pattern = "0 %")
  na_s <- stri_startswith_fixed(text_labs, pattern = "NA")
  text_labs[which(zeros)] <- "0"
  text_labs[which(na_s)] <- "NA"
  str2expression(text_labs)
}

plot_empty <- function(){
  text <- "Selected spectra file not found"
  p <- ggplot() +
    theme_void() +
    geom_text(aes(0, 0, label = text))
  return(p)
}
