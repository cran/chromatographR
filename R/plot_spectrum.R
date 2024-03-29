#' Plot spectrum from peak table
#' 
#' Plots the trace and/or spectrum for a given peak in peak.table object, or
#' plots the spectrum a particular retention time for a given chromatogram.
#' 
#' Can be used to confirm the identity of a peak or check that a particular
#' column in the peak table represents a single compound. Retention times can
#' also be selected by clicking on the plotted trace if what == 'click'.
#' 
#' @importFrom scales rescale
#' @importFrom graphics identify title text
#' @importFrom utils head tail
#' @param loc The name of the peak or retention time for which you wish to
#' extract spectral data.
#' @param peak_table The peak table (output from \code{\link{get_peaktable}}
#' function).
#' @param chrom_list A list of chromatograms in matrix form (timepoints x
#' wavelengths).
#' @param chr Numerical index of chromatogram you wish to plot, or "max" to
#' automatically plot the chromatogram with the largest signal.
#' @param lambda The wavelength you wish to plot the trace at if plot_trace ==
#' TRUE and/or the wavelength to be used for the determination of signal
#' abundance.
#' @param plot_spectrum Logical. If TRUE, plots the spectrum of the chosen
#' peak. Defaults to TRUE.
#' @param plot_trace Logical. If TRUE, plots the trace of the chosen peak at
#' lambda. Defaults to TRUE.
#' @param spectrum_labels Logical. If TRUE, plots labels on maxima in spectral
#' plot. Defaults to TRUE.
#' @param scale_spectrum Logical. If TRUE, scales spectrum to unit height.
#' Defaults to FALSE.
#' @param export_spectrum Logical. If TRUE, exports spectrum to console.
#' Defaults to FALSE.
#' @param verbose Logical. If TRUE, prints verbose output to console. Defaults
#' to TRUE.
#' @param what What to look for. Either "peak" to extract spectral information
#' for a certain peak, "rt" to scan by retention time, or "click" to manually
#' select retention time by clicking on the chromatogram. Defaults to "peak"
#' mode.
#' @param ... Additional arguments.
#' @return If \code{export_spectrum} is TRUE, returns the spectrum as a \code{
#' data.frame} with wavelengths as rows and a single column encoding the
#' absorbance (or normalized absorbance, if \code{scale_spectrum} is TRUE)
#' at each wavelength. Otherwise, there is no return value.
#' @section Side effects:
#' If \code{plot_trace} is TRUE, plots the chromatographic trace of the specified
#' chromatogram (\code{chr}), at the specified wavelength (\code{lambda}) with a
#' dotted red line to indicate the retention time given by \code{loc}. The
#' trace is a single column from the chromatographic matrix.
#'
#' If \code{plot_spectrum} is TRUE, plots the spectrum for the specified chromatogram
#' at the specified retention time. The spectrum is a single row from the chromatographic
#' matrix.
#' @author Ethan Bass
#' @examplesIf interactive()
# data(Sa)
# pks <- get_peaks(Sa,lambda="220.00000")
# pk_tab <- get_peaktable(pks)
# oldpar <- par(no.readonly = TRUE)
# par(mfrow=c(2,1))
# plot_spectrum(loc = "V10", peak_table = pk_tab, what="peak")
# par(oldpar)
#' @export plot_spectrum
#' @md
plot_spectrum <- function(loc, peak_table, chrom_list,
                          chr = 'max', lambda = 'max',
                          plot_spectrum = TRUE, plot_trace = TRUE,
                          spectrum_labels=TRUE, scale_spectrum=FALSE,
                          export_spectrum = FALSE, verbose=TRUE, 
                          what=c("peak", "rt", "click"), ...){
  what <- match.arg(what, c("peak", "rt", "click"))
  if (what == "click")
    plot_trace <- TRUE
  if (missing(chrom_list) & missing(peak_table))
    stop("Must provide either a peak_table or a chrom_list.")
  if (missing(chrom_list)){
    chrom_list <- get_chrom_list(peak_table)
  } else get_chrom_list(peak_table, chrom_list)
  if (!(class(chrom_list) %in% c("list", "chrom_list", "matrix")))
    stop("Chrom_list does not appear to be valid. Check chrom_list argument")
  if (is.matrix(chrom_list)){
    chrom_list <- list(chrom_list)
    chr <- 1
  }
  new.ts <- round(as.numeric(rownames(chrom_list[[1]])),2)
  new.lambdas <- as.numeric(colnames(chrom_list[[1]]))
  sig <- max(nchar(gsub(".*\\.","",rownames(chrom_list[[1]]))))
  if (what == "rt" | what == "click"){
    if (chr == "max")
      stop("Chromatogram must be specified for scan function.")
    if (is.null(chrom_list))
      stop("List of chromatograms must be provided for scan function.")
  }
    if (what == "peak"){
      if(missing(peak_table)){
        stop("Peak table must be provided to locate peak.")}
      if (!(loc %in% colnames(peak_table$tab))){
        stop(paste0("No match found for peak \'", loc, "\' in peak table."))}
      RT <- round(peak_table$pk_meta['rt', loc], sig)
    } else if (what == "rt"){
      RT <- round(as.numeric(loc), sig)
    }
    if (what != "click"){
      if (!is.numeric(RT))
        stop("Retention time not found!")
      if (RT > tail(new.ts, 1) | RT < head(new.ts, 1))
        stop("The supplied retention time falls outside the bounds of the chromatogram.")
      idx <- which.min(abs(RT - new.ts))
      if (chr == 'max'){
        chr <- which.max(peak_table$tab[,loc])
      }
    }
  if (is.character(chr) & !(chr %in% names(chrom_list))){
    stop("Chromatogram not found. check `chr` argument!")
  }
  if (lambda == 'max'){
    if (what == "click")
      stop("Wavelength (`lambda`) must be specified for interactive scanning.")
    y <- unlist(chrom_list[[chr]][idx, , drop=TRUE])
    lambda <- names(which.max(y))
    lambda.idx <- which(new.lambdas == lambda)
  } else{
    lambda <- as.numeric(lambda)
    lambda.idx <- which(new.lambdas == lambda) 
    if (length(lambda.idx) == 0)
      stop("The specified wavelength (`lambda`) could not be found!")
  }
  if (plot_trace){
    y_trace <- chrom_list[[chr]][,lambda.idx]
    matplot(x = new.ts, y = y_trace,type='l',
            ylab='', xlab='')
    if (what == "click"){
      message("Click trace to select timepoint")
      idx <- identify(new.ts, y_trace, n = 1, plot = FALSE)
      RT <- new.ts[idx]
    }
    abline(v = RT,col='red', lty=3)
    title(bquote(paste("\n\n Chr ", .(chr),  " ;   RT: ", .(RT), " ;  ", lambda, ": ", .(lambda), " nm",
                       #" abs: ", .(round(y_trace[idx], 2))
                       )))
    
  }
  if (verbose){
    message(paste0("chrome no. ", chr, " (`", names(chrom_list)[chr], "`) \n",
                   "RT: ", RT, "; \n",
                   "lambda = ", lambda, " nm; \n",
                   "abs = ", round(y_trace[idx], 2)))
    ### report closest match ###
    if (what != "peak" & !is.null(peak_table$tab)){
      pk <- names(which.min(abs(peak_table$pk_meta["rt",] - RT)))
      message(paste("nearest peak:", pk))
    }
  }
  y <- unlist(chrom_list[[chr]][idx, , drop=TRUE])
  if (scale_spectrum){
    y <- rescale(y)
  }
  if (plot_spectrum){
    plot_spec(y = y, spectrum_labels = spectrum_labels, ...)
  }
  if (export_spectrum){
    y <- data.frame(y)
    colnames(y) <- names(chrom_list)[chr]
    y
  }
}


#' Elementwise all equal function
#' @author Brian Diggs
#' @references https://stackoverflow.com/questions/9508518/why-are-these-numbers-not-equal
#' @noRd
elementwise.all.equal <- Vectorize(function(x, y, ...) {isTRUE(all.equal(x, y, ...))})

#' Plot all spectra for chosen peak.
#' 
#' Plot multiple for a given peak in peak table. Wrapper for
#' \code{\link{plot_spectrum}}.
#' 
#' @param peak The name of a peak to plot (in character
#' format)
#' @param peak_table The peak table (output from \code{\link{get_peaktable}}
#' function)
#' @param chrom_list A list of profile matrices, each of the same dimensions
#' (timepoints x components).
#' @param chrs Vector of chromatograms to plot.
#' @param plot_spectrum Logical. If TRUE, plots the spectrum of the chosen
#' peak.
#' @param export_spectrum Logical. If TRUE, exports spectrum to console.
#' Defaults to FALSE.
#' @param scale_spectrum Logical. If TRUE, scales spectrum to unit height.
#' @param overlapping Logical. If TRUE, plot spectra in single plot.
#' @param verbose Logical. If TRUE, prints verbose output to console.
#' @param \dots Additional arguments to plot_spectrum.
#' @return If \code{export_spectrum} is TRUE, returns the spectra as a \code{
#' data.frame} with wavelengths as rows and one column for each sample in the
#' \code{chrom_list} encoding the absorbance (or normalized absorbance, if
#' \code{scale_spectrum} is TRUE) at each wavelength. Otherwise, there is no
#' return value.
#' @section Side effects:
#' If \code{plot_spectrum} is TRUE, plots the spectra for the specified chromatogram
#' (\code{chr}) of the given \code{peak}. The spectrum is a single row
#' from the chromatographic matrix.
#' @author Ethan Bass
#' @seealso \code{\link{plot_spectrum}}
#' @examplesIf interactive()
#' data(Sa_warp)
#' pks <- get_peaks(Sa_warp, lambda="220")
#' pk_tab <- get_peaktable(pks)
#' plot_all_spectra(peak="V13", peak_table = pk_tab, overlapping=TRUE)
#' @export plot_all_spectra
#'
plot_all_spectra <- function(peak, peak_table, chrom_list, chrs="all", 
                             plot_spectrum = TRUE, export_spectrum=TRUE,
                             scale_spectrum=TRUE, overlapping=TRUE,
                             verbose=FALSE, ...){
  if (missing(chrom_list)){
    chrom_list <- get_chrom_list(peak_table)
  } else get_chrom_list(peak_table, chrom_list)
  if (!(inherits(chrom_list, "list") | inherits(chrom_list, "chrom_list")))
    stop("chrom_list is not a list")
  new.lambdas <- as.numeric(colnames(chrom_list[[1]]))
  if ("all" %in% chrs) chrs <- seq_along(chrom_list)
  sp <- sapply(chrs, function(chr){
    plot_spectrum(loc = peak, peak_table = peak_table, chrom_list=chrom_list,
                  chr=chr, plot_spectrum=FALSE, plot_trace=FALSE, 
                  export_spectrum = TRUE, scale_spectrum=scale_spectrum,
                  verbose=verbose, what="peak")
  })
  sp<-as.data.frame(do.call(cbind, sp))
  colnames(sp) <- names(chrom_list)[chrs]
  rownames(sp) <- colnames(chrom_list[[1]])
  if (plot_spectrum){
    if(overlapping){
      matplot(new.lambdas, sp, type='l', xlab='wavelength', ylab='intensity',las=2)
    } else{
      apply(sp, 2,function(spp){
        plot(new.lambdas, spp, type='l', xlab='', ylab='',las=2)
      })
  }
  }
  if(export_spectrum){
    sp
  }
}

#' Scan spectrum
#' 
#' Convenience function to call plot_spectrum with what = "click".
#' 
#' @importFrom scales rescale
#' @importFrom graphics identify title text abline
#' @param peak_table The peak table (output from \code{\link{get_peaktable}}
#' function).
#' @param chrom_list A list of chromatograms in matrix form (timepoints x
#' wavelengths).
#' @param chr Numerical index of chromatogram you wish to plot.
#' @param lambda The wavelength to plot the trace at.
#' @param spectrum_labels Logical. If TRUE, plots labels on maxima in spectral
#' plot. Defaults to TRUE.
#' @param scale_spectrum Logical. If TRUE, scales spectrum to unit height.
#' Defaults to FALSE.
#' @param export_spectrum Logical. If TRUE, exports spectrum to console.
#' Defaults to FALSE.
#' @param ... Additional arguments.
#' @return If \code{export_spectrum} is TRUE, returns the spectrum as a \code{
#' data.frame} with wavelengths as rows and a single column encoding the
#' absorbance (or normalized absorbance, if \code{scale_spectrum} is TRUE)
#' at each wavelength. Otherwise, there is no return value.
#' @section Side effects:
#' Plots a chromatographic trace from the specified chromatogram (\code{chr}),
#' at the specified wavelength (\code{lambda}) with a dotted red line to indicate
#' the user-selected retention time. The trace is a single column from the
#' chromatographic matrix.
#' 
#' If \code{plot_spectrum} is TRUE, plots the spectrum for the specified
#' chromatogram at the user-specified retention time. The spectrum is a single
#" row from the chromatographic matrix.
#' @md
#' @author Ethan Bass
#' @examplesIf interactive()
#' data(Sa_pr)
#' scan_chrom(Sa_pr, lambda="210", chr=2, export_spectrum=TRUE)
#' @export scan_chrom
scan_chrom <- function(chrom_list, lambda, chr, peak_table=NULL, 
                       scale_spectrum = FALSE, spectrum_labels = TRUE,
                       export_spectrum = FALSE, ...){
  if (!(inherits(chrom_list, "list") | inherits(chrom_list, "chrom_list")))
    stop("chrom_list is not a list")
  if (missing(lambda))
    stop("Please specify wavelength (`lambda`) to proceed with scanning.")
  if (missing(chr)){
    chr <- as.numeric(readline(
      prompt="Which chromatogram do you wish to plot? \n"))
  }
  plot_spectrum(chrom_list = chrom_list, peak_table=peak_table,
                            chr = chr, lambda = lambda, what="click",
                            scale_spectrum = scale_spectrum,
                            spectrum_labels = spectrum_labels, 
                            export_spectrum = export_spectrum, ...)
}

#' Retrieve chrom_list from peak_table object
#' @param peak_table Peak table object
#' @author Ethan Bass
#' @noRd
get_chrom_list <- function(x, chrom_list, verbose = FALSE){
  if (inherits(x, "peak_table")){
    if (missing(chrom_list)){
      chrom_list <- try(get(x$args["chrom_list"]))
      if (inherits(chrom_list, "try-error")) stop("Chromatograms not found!")
    }
    if (length(chrom_list) != nrow(x$tab)){
      stop("Dimensions of chrom_list and peak_table do not match.")
    } else{
      if (verbose & any(names(chrom_list) != rownames(x$tab)))
        warning("Names of chromatograms do not match peak_table")
    }
  } else if (inherits(x, "peak_list")){
    if (missing(chrom_list)){
      chrom_list <- try(get(attr(x, "chrom_list")))
      if (inherits(chrom_list, "try-error")) stop("Chromatograms not found!")
    }
    if (length(chrom_list) != length(x)){
      stop("Dimensions of chrom_list and peak_list do not match.")
    } else{
      if (verbose & any(names(chrom_list) != names(x)))
        warning("Names of chromatograms do not match peak_list")
    }
  }
  chrom_list
}

#' Plot spectrum
#' @param y Numeric vector containing spectral data.
#' @param spectrum_labels Logical. Whether to label peaks in spectrum.
#' @author Ethan Bass
#' @noRd
plot_spec <- function(y, spectrum_labels = TRUE, ...){
  matplot(x = as.numeric(names(y)), y = as.numeric(y), type='l',
          ylab = 'Intensity', xlab = 'Wavelength (nm)',
          ylim=c(0,max(y, na.rm = TRUE)*1.2), ...)
  if (spectrum_labels){
    suppressWarnings(pks <- find_peaks(y, slope_thresh = .00001, bounds = FALSE))
    if (length(pks)>0){
      pks <- data.frame(round(as.numeric(names(y)[pks]),0), y[pks],stringsAsFactors = FALSE)
      text(pks[,1], pks[,2], pks[,1], pos=3, offset=.3, cex = .8)
    }
  }
}