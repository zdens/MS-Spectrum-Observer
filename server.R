library(stringi)
library(plotly)
reactiveConsole(TRUE)

source('db_connect.R')
source('functions.R')

data_dir_path <- 'archive'

shinyServer(function(input, output, session) {
  
  # to relay the height/width of the plot's container, we'll query this 
  # session's client data http://shiny.rstudio.com/articles/client-data.html
  cdata <- session$clientData
  spec_info_react <- reactiveVal(NULL)
  
  files_from_db <- reactive({
    spec_files <- get_files_from_db(data_dir_path)
    spec_files$label <- paste0("File ",
                               spec_files$filename)
    sorted_labels <- sort(spec_files$filename, index.return = TRUE)
    spec_files <- spec_files[sorted_labels$ix, ]
    spec_files$label <- as.factor(spec_files$filename)
    return(list(spec_files = spec_files,
                labels = c("Select spectra file",
                           sorted_labels$x)))
  })
  
  readSpectra <- reactive({
    #browser()
    spec_files <- files_from_db()
    loaded.file <- input$input.select.file
    if (is.null(loaded.file) || loaded.file == "Select spectra file") {
      return(list(files_inf = spec_files,
                  spectra = NULL))
    }
    idx <- which(spec_files$spec_files$label == loaded.file)
    if (length(idx) == 0) {
      return(list(files_inf = spec_files,
                  spectra = NULL))
    }
    filename <- spec_files$spec_files$filename[[idx]]
    fpath <- file.path(data_dir_path, filename)
    if (file.exists(fpath)) {
      spec_data <- readRDS(fpath)
      spectra_info <- list(modes = spec_data$modes,
                           spectra = spec_data$spectra,
                           pr = spec_data$pr)
      spectra_info$files_inf <- spec_files
      spec_info_react(spectra_info)
      return(spectra_info)
    } else {
      return(list(files_inf = spec_files,
                  spectra = NULL))
    }
    
  })
  
  observe({
    spectra_info <- spec_info_react()
    if (is.null(spectra_info$spectra)) {
      return()
    }
    updateSelectInput(session, "input.select.mode",
                      choices = spectra_info$modes,
                      selected = spectra_info$modes[1])
  })
  observe({
    # browser()
    files_list <- readSpectra()$files_inf
    if (is.null(files_list)) {
      updateSelectInput(session, "input.select.file",
                        choices = "Failed to get files from DB")
    } else {
      selected <- input$input.select.file
      updateSelectInput(session, "input.select.file",
                        choices = files_list$labels,
                        selected = selected)
    }
  })
  
  spectraList <- reactive({
    if (input$input.select.mode == "Select spectra file") {
      return(NULL)
    }
    spectra_info <- spec_info_react()
    if (is.null(spectra_info$spectra)) {
      return(NULL)
    }
    spec_mode_idx <- which(spectra_info$pr == input$input.select.mode)
    spectra_list <- spectra_info$spectra[spec_mode_idx]
    return(list(spectra_list = spectra_list,
                pr = spectra_info$pr[spec_mode_idx]))
  })
  scanLines <- reactive({
    spectra_list <- spectraList()
    if (is.null(spectra_list)) {
      return(NULL)
    }
    scan.num.l <- input$input.scan.num - input$input.scan.num.l
    scan.num.r <- input$input.scan.num + input$input.scan.num.r
    scan.df <- data.frame(intensity = spectra_list$spectra_list[[input$input.scan.num]]@intensity,
                          mz = spectra_list$spectra_list[[input$input.scan.num]]@mass,
                          scan = "current_scan")
    scan.df.l <- data.frame(intensity = spectra_list$spectra_list[[scan.num.l]]@intensity,
                            mz = spectra_list$spectra_list[[scan.num.l]]@mass,
                            scan = "left_scan")
    scan.df.r <- data.frame(intensity = spectra_list$spectra_list[[scan.num.r]]@intensity,
                            mz = spectra_list$spectra_list[[scan.num.r]]@mass,
                            scan = "right_scan")
    .md <- spectra_list$spectra_list[[input$input.scan.num]]@metaData
    .md.l <- spectra_list$spectra_list[[scan.num.l]]@metaData
    .md.r <- spectra_list$spectra_list[[scan.num.r]]@metaData
    title.info <-
      paste0("Raw scan# ",
             paste0(c(.md.l$num, .md$num, .md.r$num), collapse = ", "),
             ", retention time: ",
             paste0(c(.md.l$retentionTime, .md$retentionTime, .md.r$retentionTime), collapse = ", "),
             ", points count# ",
             length(spectra_list$spectra_list[[input$input.scan.num]]@intensity > 0))
    list(center = scan.df,
         left = scan.df.l,
         right = scan.df.r,
         title = title.info)
  })
  peaksDf <- reactive({
    spectra_list <- spectraList()
    if (is.null(spectra_list)) {
      return(NULL)
    }
    scan.num.l <- input$input.scan.num - input$input.scan.num.l
    scan.num.r <- input$input.scan.num + input$input.scan.num.r
    peaks <- detectPeaks(spectra_list$spectra_list, halfWindowSize = input$input.halfWindowSize,
                         method = "SuperSmoother", SNR = input$input.snr)
    p.df <- data.frame(intensity = peaks[[input$input.scan.num]]@intensity,
                       mz = peaks[[input$input.scan.num]]@mass,
                       scan = "current_scan")
    p.df <- rbind(p.df,
                  data.frame(intensity = peaks[[scan.num.l]]@intensity,
                             mz = peaks[[scan.num.l]]@mass,
                             scan = "left_scan"),
                  data.frame(intensity = peaks[[scan.num.r]]@intensity,
                             mz = peaks[[scan.num.r]]@mass,
                             scan = "right_scan"))
    p.df$scan <- factor(p.df$scan)
    .md <- peaks[[input$input.scan.num]]@metaData
    .md.l <- peaks[[scan.num.l]]@metaData
    .md.r <- peaks[[scan.num.r]]@metaData
    title.info <-
      paste0("Peaks scan# ",
             paste0(c(.md.l$num, .md$num, .md.r$num), collapse = ", "),
             ", retention time: ",
             paste0(c(.md.l$retentionTime, .md$retentionTime, .md.r$retentionTime), collapse = ", "),
             ", peaks count# ",
             length(peaks[[input$input.scan.num]]@intensity > 0))
    return(list(peaks = peaks,
                peaksDf = p.df,
                title = title.info))
  })
  spectraAligned <- reactive({
    spectra_list <- spectraList()
    if (is.null(spectra_list)) {
      return(NULL)
    }
    scan.num.l <- input$input.scan.num - input$input.scan.num.l
    scan.num.r <- input$input.scan.num + input$input.scan.num.r
    spec.tic <- sapply(spectra_list$spectra_list, totalIonCurrent)
    spectra_aligned <- alignSpectra(
      spectra_list$spectra_list,
      reference = spectra_list$spectra_list[[which.max(spec.tic)]],
      halfWindowSize = input$input.halfWindowSize,
      tolerance = input$input.tol.align,
      SNR = input$input.snr
    )
    p.df <- data.frame(intensity = spectra_aligned[[input$input.scan.num]]@intensity,
                       mz = spectra_aligned[[input$input.scan.num]]@mass,
                       scan = "current_scan")
    p.df <- rbind(p.df,
                  data.frame(intensity = spectra_aligned[[scan.num.l]]@intensity,
                             mz = spectra_aligned[[scan.num.l]]@mass,
                             scan = "left_scan"),
                  data.frame(intensity = spectra_aligned[[scan.num.r]]@intensity,
                             mz = spectra_aligned[[scan.num.r]]@mass,
                             scan = "right_scan"))
    .md <- spectra_aligned[[input$input.scan.num]]@metaData
    .md.l <- spectra_aligned[[scan.num.l]]@metaData
    .md.r <- spectra_aligned[[scan.num.r]]@metaData
    p.df$scan <- factor(p.df$scan)
    title.info <-
      paste0("Aligned spectra scan# ",
             paste0(c(.md.l$num, .md$num, .md.r$num), collapse = ", "),
             ", retention time: ",
             paste0(c(.md.l$retentionTime, .md$retentionTime, .md.r$retentionTime), collapse = ", "),
             ", points count# ",
             length(spectra_aligned[[input$input.scan.num]]@intensity > 0))
    return(list(aligned_spectra = spectra_aligned,
                aligned_spectra_df = p.df,
                title = title.info,
                pr = spectra_list$pr))
  })
  detectedAlignedPeaks <- reactive({
    scan.num.l <- input$input.scan.num - input$input.scan.num.l
    scan.num.r <- input$input.scan.num + input$input.scan.num.r
    spec_aligned_list <- spectraAligned()
    if (is.null(spec_aligned_list)) {
      return(NULL)
    }
    peaks <- detectPeaks(
      spec_aligned_list$aligned_spectra,
      SNR = input$input.snr,
      halfWindowSize = input$input.halfWindowSize,
      method = "SuperSmoother"
    )
    p.df <- data.frame(intensity = peaks[[input$input.scan.num]]@intensity,
                       mz = peaks[[input$input.scan.num]]@mass,
                       scan = "current_scan")
    p.df <- rbind(p.df,
                  data.frame(intensity = peaks[[scan.num.l]]@intensity,
                             mz = peaks[[scan.num.l]]@mass,
                             scan = "left_scan"),
                  data.frame(intensity = peaks[[scan.num.r]]@intensity,
                             mz = peaks[[scan.num.r]]@mass,
                             scan = "right_scan"))
    .md <- peaks[[input$input.scan.num]]@metaData
    .md.l <- peaks[[scan.num.l]]@metaData
    .md.r <- peaks[[scan.num.r]]@metaData
    p.df$scan <- factor(p.df$scan)
    title.info <-
      paste0("Detected peaks in aligned scan# ",
             paste0(c(.md.l$num, .md$num, .md.r$num), collapse = ", "),
             ", retention time: ",
             paste0(c(.md.l$retentionTime, .md$retentionTime, .md.r$retentionTime), collapse = ", "),
             ", peaks count# ",
             length(peaks[[input$input.scan.num]]@intensity > 0))
    return(list(detected_peaks = peaks,
                detected_peaks_df = p.df,
                title = title.info,
                aligned_spectra = spec_aligned_list$aligned_spectra,
                pr = spec_aligned_list$pr))
  })
  binPeaksPlot <- reactive({
    scan.num.l <- input$input.scan.num - input$input.scan.num.l
    scan.num.r <- input$input.scan.num + input$input.scan.num.r
    detected_peaks <- detectedAlignedPeaks()
    if (is.null(detected_peaks)) {
      return(NULL)
    }
    peaks <-
      remove_negative_snr(detected_peaks$detected_peaks,
                          detected_peaks$aligned_spectra, 
                          halfWindowSize = input$input.halfWindowSize,
                          snr = input$input.snr)
    idx_gt_10 <- which(sapply(peaks, length) > 10)
    peaks2 <- binPeaks(peaks[idx_gt_10],
                       method = "strict",
                       tolerance = input$input.tol.bin.peaks)
    labels <- detected_peaks$pr
    fpeaks <-
      filterPeaks(
        peaks2,
        minNumber = 10,
        labels = labels[idx_gt_10],
        mergeWhitelists = FALSE
      )
    p.df <- data.frame(intensity = fpeaks[[input$input.scan.num]]@intensity,
                       mz = fpeaks[[input$input.scan.num]]@mass,
                       scan = "current_scan")
    p.df <- rbind(p.df,
                  data.frame(intensity = fpeaks[[scan.num.l]]@intensity,
                             mz = fpeaks[[scan.num.l]]@mass,
                             scan = "left_scan"),
                  data.frame(intensity = fpeaks[[scan.num.r]]@intensity,
                             mz = fpeaks[[scan.num.r]]@mass,
                             scan = "right_scan"))
    .md <- fpeaks[[input$input.scan.num]]@metaData
    .md.l <- fpeaks[[scan.num.l]]@metaData
    .md.r <- fpeaks[[scan.num.r]]@metaData
    p.df$scan <- factor(p.df$scan)
    title.info <-
      paste0("Detected filtered peaks scan# ",
             paste0(c(.md.l$num, .md$num, .md.r$num), collapse = ", "),
             ", retention time: ",
             paste0(c(.md.l$retentionTime, .md$retentionTime, .md.r$retentionTime), collapse = ", "),
             ", peaks count# ",
             length(fpeaks[[input$input.scan.num]]@intensity > 0))
    mdL = lapply(detected_peaks$aligned_spectra[idx_gt_10], function(x)metaData(x))
    return(list(bin_f_peaks = fpeaks,
                title = title.info,
                bin_f_peaks_df = p.df,
                spec_aligned = detected_peaks$aligned_spectra[idx_gt_10],
                mdL = mdL))
  })
  featureMatrixPlot <- reactive({
    scan.num.l <- input$input.scan.num - input$input.scan.num.l
    scan.num.r <- input$input.scan.num + input$input.scan.num.r
    fpeaks <- binPeaksPlot()
    if (is.null(fpeaks)) {
      return(NULL)
    }
    featureMatrix <-
      intensityMatrix(fpeaks$bin_f_peaks, fpeaks$spec_aligned)
    mz.list <- attr(featureMatrix, "mass")
    .df <- data.frame(intensity = featureMatrix[input$input.scan.num, ],
                      mz = mz.list,
                      scan = "current_scan")
    .df <- rbind(.df,
                 data.frame(intensity = featureMatrix[scan.num.l, ],
                            mz = mz.list,
                            scan = "left_scan"),
                 data.frame(intensity = featureMatrix[scan.num.r, ],
                            mz = mz.list,
                            scan = "right_scan"),
                 data.frame(intensity = apply(featureMatrix[-input$input.scan.num, ], 2, max),
                            mz = mz.list,
                            scan = "others"))
    .df$scan <- factor(.df$scan)
    .md <- fpeaks$mdL[[input$input.scan.num]]
    ncols <- length(which(featureMatrix[input$input.scan.num,] > 0))
    title.info <-
      paste0("Intensities matrix, scan# ",
             .md$num,
             ", retention time: ",
             .md$retentionTime,
             ", features count# ",
             ncols)
    return(list(fm = featureMatrix,
                fm_df = .df,
                title = title.info))
  })
  
  output$pid.raw <- renderPlotly({
    #browser()
    spectra_list <- spectraList()
    if (is.null(spectra_list)) {
      return(plot_empty())
    }
    scan_lines <- scanLines()
    scan.num.l <- input$input.scan.num - input$input.scan.num.l
    scan.num.r <- input$input.scan.num + input$input.scan.num.r
    scan.df <- rbind(scan_lines$center,
                     scan_lines$left,
                     scan_lines$right)
    scan.df$scan <- factor(scan.df$scan)
    get_plotly <- function(segment_df, lines_df, title_text, showlegend = FALSE) {
      ggp <- plot_ly(segment_df, legendgroup = ~scan,
                     showlegend = FALSE) %>%
        add_segments(x = ~mz,
                     y = 0,
                     xend = ~mz,
                     yend = ~intensity,
                     color = ~scan,
                     showlegend = showlegend) %>%
        add_lines(data = lines_df, x = ~mz, y = ~intensity, color = ~scan,
                  showlegend = FALSE) %>%
        add_annotations(text = title_text, x = 0.85,
                        y = 0.9,
                        yref = "paper",
                        xref = "paper",
                        yanchor = "bottom",
                        showarrow = FALSE,
                        font = list(size = 13))
      return(ggp)
    }
    ggp1 <- get_plotly(scan.df, scan.df, scan_lines$title)
    
    detected_peaks <- detectedAlignedPeaks()
    if (is.null(detected_peaks)) {
      return(plot_empty())
    }
    ggp4 <- get_plotly(detected_peaks$detected_peaks_df, scan.df, detected_peaks$title)
    
    data_to_plot <- binPeaksPlot()
    if (is.null(data_to_plot)) {
      return(plot_empty())
    }
    ggp5 <- get_plotly(data_to_plot$bin_f_peaks_df, scan.df, data_to_plot$title)
    
    data_to_plot <- featureMatrixPlot()
    if (is.null(data_to_plot)) {
      return(plot_empty())
    }
    sub_title <- paste0(
      "SNR: ",
      input$input.snr,
      ", ",
      "HalfWinSize: ",
      input$input.halfWindowSize,
      ", ",
      "TolAlign: ",
      input$input.tol.align,
      ", ",
      "TolBinPeaks: ",
      input$input.tol.bin.peaks
    )
    
    ggp6 <- get_plotly(data_to_plot$fm_df, scan.df, data_to_plot$title, TRUE)
    plotly_title <- paste0("Spectra file: ",
                           input$input.select.file,
                           "<br><sup>",
                           sub_title,
                           "</sup>")
    
    subplot(list(ggp1, ggp4, ggp5, ggp6), nrows = 4, shareX = TRUE,
            heights = rep(1/4, length.out = 4), margin = 0.01) %>% 
      layout(legend=list(title=list(text='<b> Scans </b>')),
             title = list(text = plotly_title, y = 0.99, yanchor = "top", x = 0.02))
  })
})
