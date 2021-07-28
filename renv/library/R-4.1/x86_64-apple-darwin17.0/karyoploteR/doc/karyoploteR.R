## ---- include=FALSE-----------------------------------------------------------
library(knitr)
library(karyoploteR)
library(regioneR)
library(magrittr)
opts_chunk$set(concordance=FALSE)
set.seed(21666641)

## -----------------------------------------------------------------------------
  gains <- toGRanges(data.frame(chr=c("chr1", "chr5", "chr17", "chr22"), start=c(1, 1000000, 4000000, 1),
                      end=c(5000000, 3200000, 80000000, 1200000)))
  losses <- toGRanges(data.frame(chr=c("chr3", "chr9", "chr17"), start=c(80000000, 20000000, 1),
                       end=c(170000000, 30000000, 25000000)))
  kp <- plotKaryotype(genome="hg19")
  kpPlotRegions(kp, gains, col="#FFAACC")
  kpPlotRegions(kp, losses, col="#CCFFAA")


## -----------------------------------------------------------------------------
  kp <- plotKaryotype(genome="mm10", plot.type=1, main="The mm10 genome")

## -----------------------------------------------------------------------------
  kp <- plotKaryotype(genome="hg19", plot.type=2, chromosomes=c("chr1", "chr2", "chr3"))

## -----------------------------------------------------------------------------
  rand.data <- createRandomRegions(genome="hg19", nregions=1000, length.mean=1, length.sd=0,
                      mask=NA, non.overlapping=TRUE) 
  rand.data <- toDataframe(sort(rand.data))
  rand.data <- cbind(rand.data, y=runif(n=1000, min=-1, max=1))
  
  #Select some data points as "special ones"
  sel.data <- rand.data[c(7, 30, 38, 52),] 
  head(rand.data)

## -----------------------------------------------------------------------------
  kp <- plotKaryotype(genome="hg19", plot.type=2, chromosomes=c("chr1", "chr2", "chr3"))
  
  kpDataBackground(kp, data.panel = 1, r0=0, r1=0.45)
  kpAxis(kp, ymin=-1, ymax=1, r0=0.05, r1=0.4, col="gray50", cex=0.5)
  kpPoints(kp, chr=rand.data$chr, x=rand.data$start, y=rand.data$y,
           ymin=-1, ymax=1, r0=0.05, r1=0.4, col="black", pch=".", cex=2)
  kpPoints(kp, chr=sel.data$chr, x=sel.data$start, y=sel.data$y,
           ymin=-1, ymax=1, r0=0.05, r1=0.4, col="red")
  kpText(kp, chr=sel.data$chr, x=sel.data$start, y=sel.data$y,
         ymin=-1, ymax=1, r0=0.05, r1=0.4, labels=c("A", "B", "C", "D"), col="red",
         pos=4, cex=0.8)
  
  
  #Upper part: data.panel=1
  kpDataBackground(kp, data.panel = 1, r0=0.5, r1=1)
  kpAxis(kp, ymin=-1, ymax=1, r0=0.5, r1=1, col="gray50", cex=0.5, numticks = 5)
  kpAbline(kp, h=c(-0.5, 0, 0.5), col="gray50", ymin=-1, ymax=1, r0=0.5, r1=1)
  kpLines(kp, chr=rand.data$chr, x=rand.data$start, y=rand.data$y,
          col="#AA88FF", ymin=-1, ymax=1, r0=0.5, r1=1)
  #Use kpSegments to add small tic to the line
  kpSegments(kp, chr=rand.data$chr, x0=rand.data$start, x1=rand.data$start,
             y0=rand.data$y-0.1, y1=rand.data$y+0.1,
             col="#8866DD", ymin=-1, ymax=1, r0=0.5, r1=1)
  #Plot the same line but inverting the data by pssing a r0>r1
  kpLines(kp, chr=rand.data$chr, x=rand.data$start, y=rand.data$y,
          col="#FF88AA", ymin=-1, ymax=1, r0=1, r1=0.5)
  
  
  #Lower part: data.panel=2
  kpDataBackground(kp, r0=0, r1=0.29, color = "#EEFFEE", data.panel = 2)
  kpAxis(kp, col="#AADDAA", ymin=-1, ymax=1, r0=0, r1=0.29, data.panel = 2,
         numticks = 2, cex=0.5, tick.len = 0)
  kpAbline(kp, h=0, col="#AADDAA", ymin=-1, ymax=1, r0=0, r1=0.29, data.panel = 2)
  kpBars(kp, chr=rand.data$chr, x0=rand.data$start, x1=rand.data$end, y1 = rand.data$y,
          col="#AADDAA", ymin=-1, ymax=1, r0=0, r1=0.29, data.panel = 2, border="#AADDAA" )
  
  kpDataBackground(kp, r0=0.34, r1=0.63, color = "#EEEEFF", data.panel = 2)
  kpAxis(kp, col="#AAAADD", ymin=-1, ymax=1, r0=0.34, r1=0.63, data.panel = 2, 
         numticks = 2, cex=0.5, tick.len = 0)
  kpAbline(kp, h=0, col="#AAAADD", ymin=-1, ymax=1, r0=0.34, r1=0.63, data.panel = 2)
  kpSegments(kp, chr=rand.data$chr, x0=rand.data$start, x1=rand.data$end, 
             y0=rand.data$y-0.2, y1=rand.data$y, 
             col="#AAAADD", ymin=-1, ymax=1, r0=0.34, r1=0.63, data.panel = 2, lwd=2)
  
  kpDataBackground(kp, r0=0.68, r1=0.97, color = "#FFEEEE", data.panel = 2)
  kpAxis(kp, col="#DDAAAA", ymin=-1, ymax=1, r0=0.68, r1=0.97, data.panel = 2,
         numticks = 2, cex=0.5, tick.len = 0)
  kpPoints(kp, chr=rand.data$chr, x=rand.data$start, y=rand.data$y,
          col="#DDAAAA", ymin=-1, ymax=1, r0=0.68, r1=0.97, data.panel = 2, pch=".", cex=3)
 

## -----------------------------------------------------------------------------
  n <- 10 #number of sets of random regions
  #Create random regions
  all.regions <- lapply(seq_len(n), function(i) {
    filterChromosomes(createRandomRegions(nregions=2000, length.mean=100000))
  })
  #join all regions into a single GRanges for coverage plotting
  joined.regions <- do.call(c, all.regions)
  
  #Create the plot
  kp <- plotKaryotype(plot.type=2)
  #Plot the masked regions
  kpPlotRegions(kp, data=getGenomeAndMask("hg19")$mask, r0=0,
               r1=1, col="lightgray", border="lightgray")
  #Plot the random regions (using invisible to ignore the return value)
  invisible(lapply(seq_len(n), function(i) {
    kpPlotRegions(kp, all.regions[[i]],
                  r0=(1/n)*(i-1), r1=(1/n)*i)
  }))
  #Plot the coverage of random regions
  kpPlotCoverage(kp, data=joined.regions, data.panel=2, col="#AAAAFF")

## -----------------------------------------------------------------------------
  kp <- plotKaryotype()

## -----------------------------------------------------------------------------
  kp <- plotKaryotype(chromosomes=c("autosomal"))

## -----------------------------------------------------------------------------
  kp <- plotKaryotype(genome="mm10", chromosomes=c("chr10", "chr13", "chr2"))

## -----------------------------------------------------------------------------
  kp <- plotKaryotype(chromosomes = c("chr1", "chr2"), plot.type = 1)
  kpDataBackground(kp, data.panel = 1)

## -----------------------------------------------------------------------------
 kp <- plotKaryotype(chromosomes = c("chr1", "chr2"), plot.type = 2)
 kpDataBackground(kp, data.panel = 1)
 kpDataBackground(kp, data.panel = 2)

## -----------------------------------------------------------------------------
 kp <- plotKaryotype(chromosomes = c("chr1", "chr2"), plot.type = 3)
 kpDataBackground(kp, data.panel = 1)
 kpDataBackground(kp, data.panel = 2)

## -----------------------------------------------------------------------------
 kp <- plotKaryotype(chromosomes = c("chr1", "chr2"), plot.type = 4)
 kpDataBackground(kp, data.panel = 1)

## -----------------------------------------------------------------------------
 kp <- plotKaryotype(chromosomes = c("chr1", "chr2"), plot.type = 5)
 kpDataBackground(kp, data.panel = 1)

## ---- eval=FALSE--------------------------------------------------------------
#   # NOT RUN
#   pt <- permTest(A=my.regions, B=repeats, randomize.function=randomizeRegions,
#   evaluate.function=overlapRegions)

## -----------------------------------------------------------------------------
  kp <- plotKaryotype(chromosomes=c("chr1", "chr2"), plot.type=2)
 
  #data.panel=1
  kpDataBackground(kp)
  #Default axis
    kpAxis(kp) 
  #Axis on the right side of the data.panel
    kpAxis(kp, side = 2) 
  
  #data.panel=2
  kpDataBackground(kp, r1=0.47, data.panel=2)
  #Changing the limits and having more ticks, with a smaller font size
  kpAxis(kp, r1=0.47, ymin=-5000, ymax = 5000, numticks = 5, cex=0.5, data.panel=2)  
  #and a different scale on the right
  kpAxis(kp, r1=0.47, ymin=-2, ymax = 2, numticks = 3, cex=0.5, data.panel=2, side=2)
  
  kpDataBackground(kp, r0=0.53, data.panel=2)
  #Changing the colors and labels and tick positions
  kpAxis(kp, r0=0.53, tick.pos = c(0.3, 0.6, 1), labels = c("A", "B", "C"), col="#66AADD",
         cex=0.5, data.panel=2)
  

## -----------------------------------------------------------------------------
plotDefaultPlotParams(plot.type=2)

## -----------------------------------------------------------------------------
plotDefaultPlotParams(plot.type=3)

## -----------------------------------------------------------------------------
plot.params <- getDefaultPlotParams(plot.type=2)
plot.params$ideogramheight <- 5
plot.params$data2height <- 50
plot.params$leftmargin <- 0.05
plot.params$bottommargin <- 20
plot.params$topmargin <- 20

plotDefaultPlotParams(plot.type=2, plot.params=plot.params)

## -----------------------------------------------------------------------------
plot.params <- getDefaultPlotParams(plot.type=3)
plot.params$ideogramheight <- 50
plot.params$data1height <- 50
plot.params$data1inmargin <- 1
plot.params$data2height <- 400

plotDefaultPlotParams(plot.type = 3, plot.params=plot.params)


## -----------------------------------------------------------------------------
  data.points <- data.frame(chr="chr1", pos=(1:240)*1e6, value=rnorm(240, 0.5, 0.1))
  kp <- plotKaryotype(plot.type = 4, chromosomes = "chr1")
  kpDataBackground(kp, data.panel = 1)
  kpAddBaseNumbers(kp)
  kpPoints(kp, chr = data.points$chr, x=data.points$pos, y=data.points$value, col=rainbow(240))

## -----------------------------------------------------------------------------
  detail.region <- toGRanges(data.frame("chr1", 180e6, 240e6))
  data.points <- data.frame(chr="chr1", pos=(1:240)*1e6, value=rnorm(240, 0.5, 0.1))
  kp <- plotKaryotype(plot.type = 4, zoom=detail.region)
  kpDataBackground(kp, data.panel = 1)
  kpAddBaseNumbers(kp)
  kpPoints(kp, chr = data.points$chr, x=data.points$pos, y=data.points$value, col=rainbow(240))

## -----------------------------------------------------------------------------
  kp <- plotKaryotype(plot.type = 4, chromosomes = "chr1")
  kpDataBackground(kp, data.panel = 1)
  kpAddBaseNumbers(kp)
  kpPoints(kp, chr = "chr1", x=180e6, y=0.5, pch=16, cex=8, col="#FAC67A", clipping=TRUE)
  kpArrows(kp, chr = "chr1", x0 = 50e6, x1=170e6, y0=0.2, y1=1.1)
  kpArrows(kp, chr = "chr1", x0 = 160e6, x1=240e6, y0=0.9, y1=0.5)

## -----------------------------------------------------------------------------
  detail.region <- toGRanges(data.frame("chr1", 135e6, 180e6))
  kp <- plotKaryotype(plot.type = 4, zoom = detail.region)
  kpDataBackground(kp, data.panel = 1)
  kpAddBaseNumbers(kp)
  kpPoints(kp, chr = "chr1", x=180e6, y=0.5, pch=16, cex=8, col="#FAC67A")
  kpArrows(kp, chr = "chr1", x0 = 50e6, x1=170e6, y0=0.2, y1=1.1)
  kpArrows(kp, chr = "chr1", x0 = 160e6, x1=240e6, y0=0.9, y1=0.5)

## -----------------------------------------------------------------------------
  detail.region <- toGRanges(data.frame("chr1", 135e6, 180e6))
  kp <- plotKaryotype(plot.type = 4, zoom = detail.region)
  kpDataBackground(kp, data.panel = 1)
  kpAddBaseNumbers(kp)
  kpPoints(kp, chr = "chr1", x=180e6, y=0.5, pch=16, cex=8, col="#FAC67A")
  kpArrows(kp, chr = "chr1", x0 = 50e6, x1=170e6, y0=0.2, y1=1.1, clipping=FALSE)
  kpArrows(kp, chr = "chr1", x0 = 160e6, x1=240e6, y0=0.9, y1=0.5, clipping=FALSE)

## -----------------------------------------------------------------------------
  kp <- plotKaryotype(chromosomes=c("chr1"))
  kpDataBackground(kp)
  kpAxis(kp)
  kpPoints(kp, chr="chr1", x=30000000, y=0.2)  
  kpRect(kp, chr="chr1", x0=100000000, x1=120000000, y0=0.2, y1=0.4)

## -----------------------------------------------------------------------------
  dd <- toGRanges(data.frame(chr="chr1", start=30000000, end=30000000, y=0.2))
  dd2 <- toGRanges(data.frame(chr="chr1", start=100000000, end=120000000, y0=0.2, y1=0.4))
  
  kp <- plotKaryotype(chromosomes=c("chr1")) 
  kpDataBackground(kp)
  kpAxis(kp)
  kpPoints(kp, data=dd)
  kpPoints(kp, data=dd, y=0.8)  
  kpRect(kp, data=dd2)
  

## -----------------------------------------------------------------------------
  pp <- getDefaultPlotParams(plot.type = 1)
  pp$data1height=600
  
  tr.i <- 1/11
  tr.o <- 1/10
  
  kp <- plotKaryotype(chromosomes=c("chr1"), plot.params = pp) 
  
  dd <- toGRanges(data.frame(chr="chr1", start=end(kp$genome[1])/50*(0:49), end=end(kp$genome[1])/50*(1:50)))
  mcols(dd) <- data.frame(y=((sin(start(dd)) + rnorm(n=50, mean=0, sd=0.1))/5)+0.5)
  
  tn <- 0
  kpDataBackground(kp, r0=tr.o*tn, r1=tr.o*tn+tr.i)
  kpPoints(kp, dd, r0=tr.o*tn, r1=tr.o*tn+tr.i, pch=".", cex=2)
  kpRect(kp, chr="chr1", x0=5000000, x1=45000000, y0=0.2, y1=0.8, r0=tr.o*tn, r1=tr.o*tn+tr.i, col="#EEEEEE", border="#666666")
  kpText(kp, chr="chr1", x=25000000, y=0.5, col="red", r0=tr.o*tn, r1=tr.o*tn+tr.i, labels="kpPoints", cex=0.7)
  
  tn <- 1
  kpDataBackground(kp, r0=tr.o*tn, r1=tr.o*tn+tr.i)
  kpLines(kp, dd, r0=tr.o*tn, r1=tr.o*tn+tr.i, pch=".", cex=2)
  kpRect(kp, chr="chr1", x0=5000000, x1=45000000, y0=0.2, y1=0.8, r0=tr.o*tn, r1=tr.o*tn+tr.i, col="#EEEEEE", border="#666666")
  kpText(kp, chr="chr1", x=25000000, y=0.5, col="red", r0=tr.o*tn, r1=tr.o*tn+tr.i, labels="kpLines", cex=0.7)
  
  tn <- 2
  kpDataBackground(kp, r0=tr.o*tn, r1=tr.o*tn+tr.i)
  kpBars(kp, dd, y1=dd$y, r0=tr.o*tn, r1=tr.o*tn+tr.i, col="#AAFFAA", border="#66DD66")
  kpRect(kp, chr="chr1", x0=5000000, x1=45000000, y0=0.2, y1=0.8, r0=tr.o*tn, r1=tr.o*tn+tr.i, col="#EEEEEE", border="#666666")
  kpText(kp, chr="chr1", x=25000000, y=0.5, col="red", r0=tr.o*tn, r1=tr.o*tn+tr.i, labels="kpBars", cex=0.7)
  
  tn <- 3
  kpDataBackground(kp, r0=tr.o*tn, r1=tr.o*tn+tr.i)
  kpRect(kp, dd, y0=dd$y-0.3, y1=dd$y, r0=tr.o*tn, r1=tr.o*tn+tr.i, col="#AAAAFF", border="#6666DD")
  kpRect(kp, chr="chr1", x0=5000000, x1=45000000, y0=0.2, y1=0.8, r0=tr.o*tn, r1=tr.o*tn+tr.i, col="#EEEEEE", border="#666666")
  kpText(kp, chr="chr1", x=25000000, y=0.5, col="red", r0=tr.o*tn, r1=tr.o*tn+tr.i, labels="kpRect", cex=0.7)
  
  tn <- 4
  kpDataBackground(kp, r0=tr.o*tn, r1=tr.o*tn+tr.i)
  kpText(kp, dd, labels=as.character(1:50), cex=0.5, r0=tr.o*tn, r1=tr.o*tn+tr.i, col="#DDAADD")
  kpRect(kp, chr="chr1", x0=5000000, x1=45000000, y0=0.2, y1=0.8, r0=tr.o*tn, r1=tr.o*tn+tr.i, col="#EEEEEE", border="#666666")
  kpText(kp, chr="chr1", x=25000000, y=0.5, col="red", r0=tr.o*tn, r1=tr.o*tn+tr.i, labels="kpText", cex=0.7)
  
  tn <- 5
  kpDataBackground(kp, r0=tr.o*tn, r1=tr.o*tn+tr.i)
  kpSegments(kp, dd, y0=dd$y-0.3, y1=dd$y, r0=tr.o*tn, r1=tr.o*tn+tr.i)
  kpRect(kp, chr="chr1", x0=5000000, x1=45000000, y0=0.2, y1=0.8, r0=tr.o*tn, r1=tr.o*tn+tr.i, col="#EEEEEE", border="#666666")
  kpText(kp, chr="chr1", x=25000000, y=0.5, col="red", r0=tr.o*tn, r1=tr.o*tn+tr.i, labels="kpSegments", cex=0.7)

  tn <- 6
  kpDataBackground(kp, r0=tr.o*tn, r1=tr.o*tn+tr.i)
  kpArrows(kp, dd, y0=dd$y-0.3, y1=dd$y, r0=tr.o*tn, r1=tr.o*tn+tr.i, length=0.04)
  kpRect(kp, chr="chr1", x0=5000000, x1=45000000, y0=0.2, y1=0.8, r0=tr.o*tn, r1=tr.o*tn+tr.i, col="#EEEEEE", border="#666666")
  kpText(kp, chr="chr1", x=25000000, y=0.5, col="red", r0=tr.o*tn, r1=tr.o*tn+tr.i, labels="kpArrows", cex=0.7)

  tn <- 7
  kpDataBackground(kp, r0=tr.o*tn, r1=tr.o*tn+tr.i)
  kpHeatmap(kp, dd, r0=tr.o*tn+tr.i/4, r1=tr.o*tn+tr.i-tr.i/4, colors = c("green", "black", "red"))
  kpRect(kp, chr="chr1", x0=5000000, x1=45000000, y0=0.2, y1=0.8, r0=tr.o*tn, r1=tr.o*tn+tr.i, col="#EEEEEE", border="#666666")
  kpText(kp, chr="chr1", x=25000000, y=0.5, col="red", r0=tr.o*tn, r1=tr.o*tn+tr.i, labels="kpHeatmap", cex=0.7)

  tn <- 8
  kpDataBackground(kp, r0=tr.o*tn, r1=tr.o*tn+tr.i)
  kpPolygon(kp, dd, r0=tr.o*tn, r1=tr.o*tn+tr.i)
  kpRect(kp, chr="chr1", x0=5000000, x1=45000000, y0=0.2, y1=0.8, r0=tr.o*tn, r1=tr.o*tn+tr.i, col="#EEEEEE", border="#666666")
  kpText(kp, chr="chr1", x=25000000, y=0.5, col="red", r0=tr.o*tn, r1=tr.o*tn+tr.i, labels="kpPolygon", cex=0.7)
  
  tn <- 9
  kpDataBackground(kp, r0=tr.o*tn, r1=tr.o*tn+tr.i)
  kpAbline(kp, h=c(0.25, 0.5, 0.75), v=start(dd), r0=tr.o*tn, r1=tr.o*tn+tr.i)
  kpRect(kp, chr="chr1", x0=5000000, x1=45000000, y0=0.2, y1=0.8, r0=tr.o*tn, r1=tr.o*tn+tr.i, col="#EEEEEE", border="#666666")
  kpText(kp, chr="chr1", x=25000000, y=0.5, col="red", r0=tr.o*tn, r1=tr.o*tn+tr.i, labels="kpAbline", cex=0.7)
  
  

## ---- eval=FALSE--------------------------------------------------------------
#    kp <- plotKaryotype(genome="hg19")
#    kpDataBackground(kp)
#    kpAxis(kp)
#    kpPlotRegions(kp, gains, col="#FFAACC")
#    kpPlotRegions(kp, losses, col="#CCFFAA")
#  

## ---- eval=FALSE--------------------------------------------------------------
#    library(magrittr)
#    kp <- plotKaryotype(genome="hg19") %>% kpDataBackground() %>%  kpAxis() %>%
#      kpPlotRegions(gains, col="#FFAACC") %>%
#      kpPlotRegions(losses, col="#CCFFAA")

## ---- sessionInfo-------------------------------------------------------------
sessionInfo()

