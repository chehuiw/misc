###Use for stacking abberant profiles of multiple samples in circos plot
###Huiwen Che 
###2018/03/16

.libPaths(c(.libPaths(), "/home/hche0/R/x86_64-pc-linux-gnu-library/3.2"))

list.of.packages <- c( "circlize")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(circlize)


ReadSingleSample <- function(infh, winsize) {
    winsize = as.numeric(as.character(winsize))
    fh = read.delim(infh, head = TRUE, sep=",")
    fh = fh[,c("CHR", "BIN.START", "Z")]
    fh[,c("BIN.START", "Z")] = apply(fh[,c("BIN.START", "Z")], 2, function(x) as.numeric(as.character(x)))
    fh$WIN.END = fh$BIN.START - winsize
    fh = fh[,c("CHR", "BIN.START", "WIN.END", "Z")]
    fh$Z[fh$Z<3 & fh$Z>-3] = 0
    fh$Z[fh$Z>=3] = 1
    fh$Z[fh$Z<=-3] = -1
    fh$CHR = paste0("chr", fh$CHR)
    colnames(fh) = c("chr", "start", "end", "Z")
    return(fh)
}


FetchMultiSamplesNegZ <- function(samplelist, winsize) {
    sl = read.delim(samplelist, header = FALSE)
    dflist = list()
    for (i in (1:nrow(sl))) {
	zfh = ReadSingleSample(as.character(sl[i,]), winsize)
	zfhneg = subset(zfh, zfh$Z == -1)
	dflist[[i]] = zfhneg
    }
    return(dflist)
}


FetchMultiSamplesPosZ <- function(samplelist, winsize) {
    sl = read.delim(samplelist, header = FALSE)
    dflist = list()
    for (i in (1:nrow(sl))) {
	zfh = ReadSingleSample(as.character(sl[i,]), winsize)
	zfhpos = subset(zfh, zfh$Z == 1)
	dflist[[i]] = zfhpos
    }
    return(dflist)
}


circosplot <- function(workdir, prefix, samplelist, winsize) {
    fo = paste0(workdir, "/", prefix, ".circos.pdf")
    mfhneg = FetchMultiSamplesNegZ(samplelist, winsize)
    mfhpos = FetchMultiSamplesPosZ(samplelist, winsize)
    pdf(fo, width = 8, height = 8)
    circos.par("start.degree" = 90)
    circos.initializeWithIdeogram(chromosome.index = paste0("chr", 1:22), plotType = c("ideogram", "axis", "labels"),track.height = 0.03, ideogram.height = 0.03)
    circos.genomicTrack(mfhneg, stack = TRUE, bg.border = "#e5e5e5", bg.col = "#f9e5e5", track.height = 0.32,
        panel.fun = function(region, value, ...) {
	        circos.genomicLines(region, value, col = "#cc0000", lwd = 1.8, ...)
	})
    circos.genomicTrack(mfhpos, stack = TRUE, bg.border = "#e5e5e5", bg.col = "#e5e5f9", track.height = 0.32,
        panel.fun = function(region, value, ...) {
                circos.genomicLines(region, value, col = "#0000cc", lwd = 1.8, ...)
        })
    dev.off()
}
