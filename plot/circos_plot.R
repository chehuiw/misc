.libPaths(c(.libPaths(), "/home/hche0/R/x86_64-pc-linux-gnu-library/3.2"))

list.of.packages <- c( "circlize")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(circlize)

circoplot <- paste0(outputdir, sample.name, ".circos.pdf")
     pdf(circoplot, width=8, height=8)
     circos.par("start.degree" = 90)
     circos.initializeWithIdeogram(plotType = c("ideogram", "axis", "labels"),track.height = 0.05, ideogram.height = 0.05)
     circos.genomicTrackPlotRegion(my.sam.gipseq.profile.gain, panel.fun = function(region, value, ...) {
         circos.genomicLines(region, value, type = "l", lwd=0.1, area=TRUE, col="blue", ...)
     })
     circos.genomicTrackPlotRegion(my.sam.gipseq.profile.loss, panel.fun = function(region, value, ...) {
         circos.genomicLines(region, value, type = "l", lwd=0.1, area=TRUE, col="red", baseline="top",...)
     })
#     circos.genomicTrackPlotRegion(my.sam.gipseq.profile.list, panel.fun = function(region, value, ...) {
#         i=getI(...)
#         circos.genomicLines(region, value, type = "l", area=TRUE, col=i, ...)
#     })

     circos.clear()
