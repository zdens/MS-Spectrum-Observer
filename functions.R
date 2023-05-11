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
                                halfWindowSize = 2e-4) {
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
          SNR = 1,
          halfWindowSize = halfWindowSize,
          method = "SuperSmoother"
        )
      })
  }
  return(peaks)
}
parse_spectrum <- function(spectrum, protocol, absTicThreshold = 1000, ticThreshold = 0.01, verbose = FALSE) {
  rt <- sapply(spectrum,function(.x)metaData(.x)$retentionTime)
  num <- sapply(spectrum,function(.x)metaData(.x)$number)
  tic <- sapply(spectrum, totalIonCurrent)
  mdL <- lapply(1:length(rt),function(.x){
    .idx <- which(protocol$Acq.start <= rt[[.x]] & protocol$Acq.end >= rt[[.x]])
    if (length(.idx) > 0) {
      list(number = num[.x],
           Polarity = as.character(protocol$Polarity[.idx]),
           Detector = as.character(protocol$Detector[.idx]),
           Range = as.character(protocol$Range[.idx]),
           protocol = paste(protocol$Polarity[.idx],
                            protocol$Detector[.idx],
                            protocol$Range[.idx],sep = '.'))
    } else {
      list(number = num[.x],
           Polarity = NULL,
           Detector = NULL,
           Range = NULL,
           protocol = "Mode.Change")
    }
  })
  pr <- sapply(mdL,function(.x).x$protocol)
  change_mode_scans_idx <- which(pr == "Mode.Change")
  if (verbose) {
    cat('Number of scans obtained while mode was changed', 
        length(change_mode_scans_idx), 
        "\n")
  }
  spectra <- calibrateIntensity(spectrum[-change_mode_scans_idx], method = "TIC")
  pr <- pr[-change_mode_scans_idx]
  return(list(spectra = spectra,
              mdL = mdL[-change_mode_scans_idx],
              change.mode.idx = change_mode_scans_idx,
              pr = pr,
              tic = tic[-change_mode_scans_idx],
              rt = rt[-change_mode_scans_idx],
              modes = as.factor(unique(pr))))
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
get_feature_matrix <- function(spectra,
                               hw_size,
                               tol_align,
                               tol_bin_peaks,
                               tic,
                               labels,
                               SNR = 1) {
  spectra_aligned <- alignSpectra(
    spectra,
    reference = spectra[[which.max(tic)]],
    halfWindowSize = hw_size,
    tolerance = tol_align
  )
  peaks <-
    detectPeaks(
      spectra_aligned,
      SNR = SNR,
      halfWindowSize = hw_size,
      method = "SuperSmoother"
    )
  peaks <-
    remove_negative_snr(peaks, spectra_aligned, halfWindowSize = hw_size)
  idx_gt_10 <- which(sapply(peaks, length) > 10)
  peaks2 <- binPeaks(peaks[idx_gt_10],
                     method = "strict",
                     tolerance = tol_bin_peaks)
  fpeaks <-
    filterPeaks(
      peaks2,
      minNumber = 10,
      labels = labels[idx_gt_10],
      mergeWhitelists = FALSE
    )
  featureMatrix <-
    intensityMatrix(fpeaks, spectra_aligned[idx_gt_10])
  return(featureMatrix)
}

plot_df <- function(df, title.info, sub_title, scan.df, scan.df.l, scan.df.r) {
  ggplot(df) +
    geom_segment(aes(
      x = mz,
      y = 0,
      xend = mz,
      yend = intensity,
      colour = scan
    )) +
    geom_line(data = scan.df, aes(x = mz, y = intensity, color = scan)) +
    geom_line(data = scan.df.l, aes(x = mz, y = intensity, color = scan)) +
    geom_line(data = scan.df.r, aes(x = mz, y = intensity, color = scan)) +
    scale_y_continuous() +
    labs(
      title = title.info,
      x = "M/Z",
      y = "Интенсивность",
      subtitle = sub_title
    ) +
    theme_bw(base_size = 12, base_family = "Helvetica") +
    theme(
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12)
    ) +
    scale_fill_npg()
}
plot_empty <- function(){
  text <- "Selected spectra file not found"
  p <- ggplot() +
    theme_void() +
    geom_text(aes(0, 0, label = text))
  return(p)
}
#' Experiment steps protocol
#'
#' Read experiment steps parameters from file
#'
#' @param fpath Path to the spectra file
#'
#' @return Data.frame with experiment steps parameters
#' @export
get_protocol <- function(fpath){
  dir_path <- dirname(fpath)
  protocol_name <- "protocol.txt"
  protocol <- read.delim(file.path(dir_path, protocol_name))
  protocol$rt <- protocol$timestamp * 60
  protocol$Detector <- gsub("LTQ", "LowRes", protocol$Detector)
  protocol$Detector <- gsub("FT", "HighRes", protocol$Detector)
  protocol$Acq.start <- .parse_time_period(protocol$Acq..time..min, 1) * 60
  protocol$Acq.end <- c(.parse_time_period(protocol$Acq..time..min[1:(nrow(protocol) - 1)], 2) * 60,
                        protocol$rt[[nrow(protocol)]] + 120)
  protocol$Acq.start[[nrow(protocol)]] <- protocol$rt[[nrow(protocol)]]
  return(protocol)
}

.parse_time_period <- function(period, ret_idx){
  return(as.double(sapply(period,
                          function(x) stri_split_fixed(x, pattern = "-")[[1]][[ret_idx]])))
}
