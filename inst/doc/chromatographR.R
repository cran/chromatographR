## ----setup, include=FALSE-----------------------------------------------------
library("knitr", quietly=TRUE)
library(chromatographR)
library(parallel)
opts_chunk$set(prompt = TRUE, highlight = FALSE, comment=NA, 
               fig.width=6)
suppressMessages(require(chromatographR, quiet=TRUE))

## ----eval=F-------------------------------------------------------------------
#  # single folder
#  load_chroms(paths = path, format.in="csv")
#  
#  # multiple folders
#  path = 'foo'
#  folders <- list.files(path = path, pattern = "EXPORT3D")
#  dat <- load_chroms(folders)

## ---- data--------------------------------------------------------------------
data(Sa)

## ---- choosing parameters-----------------------------------------------------
i=2 # chromatogram number in list of data
tpoints <- as.numeric(rownames(Sa[[i]]))
lambda='200.00000'

matplot(x=tpoints, y=Sa[[i]][,lambda],
       type='l', ylab='Abs (mAU)', xlab='Time (min)')
matplot(x=tpoints, y = ptw::baseline.corr(Sa[[i]][,lambda],p=.001,lambda=1e5),
       type='l', add = T, col='blue', lty = 3)

## ---- choosing dimensions-----------------------------------------------------
new.ts <- seq(10,18.66,by=.01) # choose time-points
new.lambdas <- seq(200, 318, by = 2) # choose wavelengths

## ---- preprocessing, eval=T---------------------------------------------------
dat.pr <- preprocess(Sa, dim1=new.ts, dim2=new.lambdas,
           parallel=F, p=.001, lambda=1e5)

## ---- warp models-------------------------------------------------------------
warping.models <- correct_rt(dat.pr, what = "models", lambdas=c("210"), scale=TRUE)
warp <- correct_rt(chrom_list=dat.pr, models=warping.models, what="corrected.values")

## ---- alignment_plot, fig.height=6--------------------------------------------
par(mfrow=c(2,1))
lambdas=c('210','260')
plot.new()
ts <- as.numeric(rownames(warp[[i]]))
plot.window(xlim=c(head(ts,1), tail(ts,1)),ylim=c(0,1000))
for (i in 1:length(warp)){
  matplot(ts, warp[[i]][,lambdas],type='l',add=T)
}
legend("topright", legend="ptw", bty = "n")

plot.new()
ts <- as.numeric(rownames(dat.pr[[i]]))
plot.window(xlim=c(head(ts,1),tail(ts,1)),ylim=c(0,1000))
for (i in 1:length(dat.pr)){
  matplot(ts, dat.pr[[i]][,lambdas],type='l',add=T)
}
legend("topright", legend="raw", bty = "n")

## ---- eval=FALSE--------------------------------------------------------------
#  install.packages("VPdtw", repos="https://ethanbass.github.io/drat")
#  warp <- correct_rt(chrom_list=dat.pr, alg="vpdtw", lambda="210", what="corrected.values")

## ---- get_peaks, message=F, warning=F-----------------------------------------
pks_gauss <- get_peaks(warp, lambdas = c('210','260'), sd.max=40, fit="gaussian")
pks_egh <- get_peaks(warp, lambdas = c('210', '260'), sd.max=40, fit="egh")

## ---- plot_peaks, warning=F, fig.height=6-------------------------------------
par(mfrow=c(2,1))
plot(pks_gauss, index=1, lambda='210')
plot(pks_egh, index=1, lambda='210')

## ---- get_peaktable-----------------------------------------------------------
pk_tab <- get_peaktable(pks_egh, response = "area")
head(pk_tab$tab[,1:6])

## -----------------------------------------------------------------------------
path <- system.file("extdata", "Sa_metadata.csv", package = "chromatographR")
meta <- read.csv(path)
pk_tab <- attach_metadata(peak_table = pk_tab, metadata = meta, column="vial")
pk_tab <- normalize_data(peak_table = pk_tab, column="mass")

## -----------------------------------------------------------------------------
mirror_plot(pk_tab, lambdas = c("210","260"), var = "trt", legend_size=2)

## ---- plot spectra,fig.height=6-----------------------------------------------
oldpar <- par(no.readonly = TRUE)
par(mfrow=c(2,1))
peak="V7"
plot_spectrum(peak, peak_table = pk_tab, chrom_list=warp, verbose=F)
par(oldpar)

## ---- plot_all----------------------------------------------------------------
peak="V13"
plot_all_spectra(peak, peak_table=pk_tab, export=F, overlapping=T)

## ----eval=T-------------------------------------------------------------------
plot(pk_tab, loc = "V13", box_plot = TRUE, vars = "trt", verbose = FALSE)

## -----------------------------------------------------------------------------
sessionInfo()

