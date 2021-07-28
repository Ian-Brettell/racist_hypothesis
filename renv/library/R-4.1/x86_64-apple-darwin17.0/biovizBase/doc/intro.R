### R code from vignette source 'intro.Rnw'

###################################################
### code chunk number 1: options
###################################################
options(width=72)


###################################################
### code chunk number 2: color-scheme-option
###################################################
library(biovizBase)
## library(scales)



###################################################
### code chunk number 3: color-pal-info
###################################################
head(blind.pal.info)


###################################################
### code chunk number 4: colorBlindPal
###################################################
## with no arguments, return blind.pal.info
head(colorBlindSafePal())
## 
mypalFun <- colorBlindSafePal("Set2")
## mypalFun(12, repeatable = FALSE) #only three
mypalFun(11, repeatable = TRUE)  #repeat


###################################################
### code chunk number 5: color-dichromat
###################################################
## for palette "Paried"
mypalFun <- colorBlindSafePal(21)   
par(mfrow = c(1, 3))
showColor(mypalFun(4))
library(dichromat)
showColor(dichromat(mypalFun(4), "deutan"))
showColor(dichromat(mypalFun(4), "protan"))


###################################################
### code chunk number 6: cytoband-color
###################################################
getOption("biovizBase")$cytobandColor
getBioColor("CYTOBAND")
## differece source from default or options.
opts <- getOption("biovizBase")
opts$DNABasesNColor[1] <- "red"
options(biovizBase = opts)
## get from option(default)
getBioColor("DNA_BASES_N")
## get default fixed color
getBioColor("DNA_BASES_N", source = "default")
seqs <- c("A", "C", "T", "G", "G", "G", "C")
## get colors for a sequence.
getBioColor("DNA_BASES_N")[seqs]


###################################################
### code chunk number 7: cytoband-color-legend
###################################################
cols <- getBioColor("CYTOBAND")
plotColorLegend(cols, title  = "cytoband")


###################################################
### code chunk number 8: strand-color
###################################################
par(mfrow = c(1, 3))
cols <- getBioColor("STRAND")
showColor(cols)
showColor(dichromat(cols, "deutan"))
showColor(dichromat(cols, "protan"))


###################################################
### code chunk number 9: base-color
###################################################
getBioColor("DNA_BASES_N")



###################################################
### code chunk number 10: ne-dichromat
###################################################
par(mfrow = c(1, 3))
cols <- getBioColor("DNA_BASES_N", "default")
showColor(cols, "name")
cols.deu <- dichromat(cols, "deutan")
names(cols.deu) <- names(cols)
cols.pro <- dichromat(cols, "protan")
names(cols.pro) <- names(cols)
showColor(cols.deu, "name")
showColor(cols.pro, "name")


###################################################
### code chunk number 11: sample-gr
###################################################
library(GenomicRanges)
set.seed(1)
N <- 500
gr <- GRanges(seqnames =
              sample(c("chr1", "chr2", "chr3", "chrX", "chrY"),
                     size = N, replace = TRUE),
              IRanges(
                      start = sample(1:300, size = N, replace = TRUE),
                      width = sample(70:75, size = N,replace = TRUE)),
              strand = sample(c("+", "-", "*"), size = N,
                replace = TRUE),
              value = rnorm(N, 10, 3), score = rnorm(N, 100, 30),
              group = sample(c("Normal", "Tumor"),
                size = N, replace = TRUE),
              pair = sample(letters, size = N,
                replace = TRUE))


###################################################
### code chunk number 12: addStepping-gr
###################################################
head(addStepping(gr))
head(addStepping(gr, group.name = "pair"))
gr.close <- GRanges(c("chr1", "chr1"), IRanges(c(10, 20), width = 9))
addStepping(gr.close)
addStepping(gr.close, extend.size = 5)


###################################################
### code chunk number 13: maxGap
###################################################
gr.temp <- GRanges("chr1", IRanges(start = c(100, 250),
                                 end = c(200, 300)))
maxGap(gaps(gr.temp, start = min(start(gr.temp))))
maxGap(gaps(gr.temp, start = min(start(gr.temp))), ratio = 0.5)


###################################################
### code chunk number 14: shrink-single
###################################################
gr1 <- GRanges("chr1", IRanges(start = c(100, 300, 600),
                               end = c(200, 400, 800)))
shrink.fun1 <- shrinkageFun(gaps(gr1), max.gap = maxGap(gaps(gr1), 0.15))
shrink.fun2 <- shrinkageFun(gaps(gr1), max.gap = 0)
head(shrink.fun1(gr1))
head(shrink.fun2(gr1))


###################################################
### code chunk number 15: shrinkageFun
###################################################
gr2 <- GRanges("chr1", IRanges(start = c(100, 350, 550),
                               end = c(220, 500, 900)))
gaps.gr <- intersect(gaps(gr1, start = min(start(gr1))),
                     gaps(gr2, start = min(start(gr2))))
shrink.fun <- shrinkageFun(gaps.gr, max.gap = maxGap(gaps.gr))
head(shrink.fun(gr1))
head(shrink.fun(gr2))


###################################################
### code chunk number 16: gc-content (eval = FALSE)
###################################################
## library(BSgenome.Hsapiens.UCSC.hg19)
## GCcontent(Hsapiens, GRanges("chr1", IRanges(1e6, 1e6 + 1000)))
## GCcontent(Hsapiens, GRanges("chr1", IRanges(1e6, 1e6 + 1000)), view.width = 300)


###################################################
### code chunk number 17: PileupAsGRanges (eval = FALSE)
###################################################
## library(Rsamtools)
## data(genesymbol)
## library(BSgenome.Hsapiens.UCSC.hg19)    
## bamfile <- system.file("extdata", "SRR027894subRBM17.bam", package="biovizBase")
## test <- pileupAsGRanges(bamfile, region = genesymbol["RBM17"])
## test.match <- pileupGRangesAsVariantTable(test, Hsapiens)
## head(test[,-7])
## head(test.match[,-5])


###################################################
### code chunk number 18: ucsc-genome (eval = FALSE)
###################################################
## library(rtracklayer)
## hg19IdeogramCyto <- getIdeogram("hg19", cytoband = TRUE)
## hg19Ideogram <- getIdeogram("hg19", cytoband = FALSE)
## unknowIdeogram <- getIdeogram()


###################################################
### code chunk number 19: ucsc-genome-db (eval = FALSE)
###################################################
## head(ucscGenomes()$db)


###################################################
### code chunk number 20: default-data
###################################################
data(hg19IdeogramCyto)
head(hg19IdeogramCyto)
data(hg19Ideogram)
head(hg19Ideogram)


###################################################
### code chunk number 21: check-ideogram
###################################################
isIdeogram(hg19IdeogramCyto)
isIdeogram(hg19Ideogram)
isSimpleIdeogram(hg19IdeogramCyto)
isSimpleIdeogram(hg19Ideogram)


###################################################
### code chunk number 22: genesymbol
###################################################
data(genesymbol)
head(genesymbol)
genesymbol["RBM17"]



###################################################
### code chunk number 23: R-session
###################################################
sessionInfo()


