## ---- include=FALSE-----------------------------------------------------------
library(knitr)
library(regioneR)
opts_chunk$set(concordance=FALSE)
set.seed(21666641)

## -----------------------------------------------------------------------------
  A <- createRandomRegions(nregions=50, length.mean=5000000, length.sd=3000000)
  B <- c(A[1:25], createRandomRegions(nregions=25, length.mean=500000, length.sd=30000))
  numOverlaps(A, B, count.once=TRUE)
  numOverlaps(randomizeRegions(A), B, count.once=TRUE)

## -----------------------------------------------------------------------------
  pt <- overlapPermTest(A=A, B=B, ntimes=50)
  pt

## ---- fig.height=6,fig.width=8------------------------------------------------
  plot(pt)

## ---- fig.height=6,fig.width=8------------------------------------------------
  lz <- localZScore(pt=pt, A=A, B=B)
  plot(lz)

## -----------------------------------------------------------------------------
  special <- toGRanges(system.file("extdata", "my.special.genes.txt", package="regioneR"))
  all.genes <- toGRanges(system.file("extdata", "all.genes.txt", package="regioneR"))
  altered <- toGRanges(system.file("extdata", "my.altered.regions.txt", package="regioneR"))

  length(special)
  length(all.genes)
  length(altered)

## -----------------------------------------------------------------------------
  numOverlaps(special, altered)

## -----------------------------------------------------------------------------
  random.RS <- resampleRegions(special, universe=all.genes)
  random.RS

## -----------------------------------------------------------------------------
  random.RS <- resampleRegions(special, universe=all.genes)
  numOverlaps(random.RS, altered)
  random.RS <- resampleRegions(special, universe=all.genes)
  numOverlaps(random.RS, altered)
  random.RS <- resampleRegions(special, universe=all.genes)
  numOverlaps(random.RS, altered)
  random.RS <- resampleRegions(special, universe=all.genes)
  numOverlaps(random.RS, altered)

## ---- eval=FALSE--------------------------------------------------------------
#    # NOT RUN
#    pt <- permTest(A=my.regions, B=repeats, randomize.function=randomizeRegions,
#    evaluate.function=numOverlaps)

## ---- eval=FALSE--------------------------------------------------------------
#    # NOT RUN
#    pt <- permTest(A=my.genes, randomize.function=resampleRegions, universe=all.genes,
#    evaluate.function=meanInRegions, x=methylation.levels.450K)

## -----------------------------------------------------------------------------
  pt <- permTest(A=special, ntimes=50, randomize.function=resampleRegions, universe=all.genes,
  evaluate.function=numOverlaps, B=altered, verbose=FALSE)

## -----------------------------------------------------------------------------
  pt
  summary(pt)

## ---- fig.height=6,fig.width=8------------------------------------------------
  plot(pt)

## ---- fig.height=6,fig.width=8------------------------------------------------
  regular <- toGRanges(system.file("extdata", "my.regular.genes.txt", package="regioneR"))

  length(regular)
  numOverlaps(regular, altered)
  pt.reg <- permTest(A=regular, ntimes=50, randomize.function=resampleRegions, universe=all.genes,
  evaluate.function=numOverlaps, B=altered, verbose=FALSE)
  pt.reg
  plot(pt.reg)

## ---- eval=FALSE--------------------------------------------------------------
#  #NOT RUN - See Figure 1
#  pt.5000 <- permTest(A=special, ntimes=5000, randomize.function=resampleRegions,
#  universe=all.genes, evaluate.function=numOverlaps, B=altered, verbose=TRUE)
#  plot(pt.5000)

## ---- eval=FALSE--------------------------------------------------------------
#    #NOT RUN - See Figure 2
#    pt.5000.reg <- permTest(A=regular, ntimes=5000, randomize.function=resampleRegions,
#    universe=all.genes, evaluate.function=numOverlaps, B=altered, verbose=TRUE)
#    plot(pt.5000.reg)

## -----------------------------------------------------------------------------
  A <- toGRanges(data.frame(chr=c("chr1", "chr1", "chr1"), start=c(20000, 50000, 100000), end=c(22000, 70000, 400000)))

## -----------------------------------------------------------------------------
  randomizeRegions(A, genome="hg19")

## -----------------------------------------------------------------------------
  randomizeRegions(A, genome="hg19", per.chromosome=TRUE)

## -----------------------------------------------------------------------------
circularRandomizeRegions(A, genome="hg19", mask=NA)

## -----------------------------------------------------------------------------
gcContent <- function(A, bsgenome, ...) {
A <- toGRanges(A)
reg.seqs <- getSeq(bsgenome, A)
base.frequency <- alphabetFrequency(reg.seqs)
gc.pct <- (sum(base.frequency[,"C"]) + sum(base.frequency[,"G"]))/sum(width(A))
return(gc.pct)
}

## -----------------------------------------------------------------------------
permuteRegionsMetadata <- function(A, ...) {
A <- toGRanges(A)
mcols(A)[,1] <- mcols(A)[sample(length(A)),1]
return(A)
}

## -----------------------------------------------------------------------------
A <- createRandomRegions(nregions=50, length.mean=5000000, length.sd=3000000)
B <- c(A[1:25], createRandomRegions(nregions=25, length.mean=500000, length.sd=30000))

#Without mc.set.seed=FALSE
set.seed(123)
overlapPermTest(A, B, ntimes=10, force.parallel=TRUE)
set.seed(123)
overlapPermTest(A, B, ntimes=10, force.parallel=TRUE)

#With mc.set.seed=FALSE
set.seed(123)
overlapPermTest(A, B, ntimes=10, mc.set.seed=FALSE, force.parallel=TRUE)
set.seed(123)
overlapPermTest(A, B, ntimes=10, mc.set.seed=FALSE, force.parallel=TRUE)


## ---- fig.height=6,fig.width=8------------------------------------------------
pt <- permTest(A=special, ntimes=50, randomize.function=resampleRegions, universe=all.genes,
         evaluate.function=numOverlaps, B=altered, verbose=FALSE)
lz <- localZScore(A=special, pt=pt, window=10*mean(width(special)), 
            step=mean(width(special))/2, B=altered)
plot(lz)

## ---- fig.height=6,fig.width=8------------------------------------------------

genome <- filterChromosomes(getGenome("hg19"), keep.chr="chr1")
B <- createRandomRegions(nregions=100, length.mean=10000, length.sd=20000, genome=genome,
                   non.overlapping=FALSE)
A <- B[sample(20)]

pt <- overlapPermTest(A=A, B=B, ntimes=100, genome=genome, non.overlapping=FALSE)
pt

lz <- localZScore(A=A, B=B, pt=pt, window=10*mean(width(A)), step=mean(width(A))/2 )
plot(lz)

## -----------------------------------------------------------------------------
A <- data.frame(chr=1, start=c(1,15,24), end=c(10,20,30),  x=c(1,2,3), y=c("a","b","c"))
B <- toGRanges(A)
B

## -----------------------------------------------------------------------------
toDataframe(B) 

## -----------------------------------------------------------------------------
human.genome <- getGenomeAndMask("hg19", mask=NA)$genome
human.canonical <- filterChromosomes(human.genome, organism="hg")
listChrTypes()
human.autosomal <- filterChromosomes(human.genome, organism="hg", chr.type="autosomal")
human.123 <- filterChromosomes(human.genome, keep.chr=c("chr1", "chr2", "chr3"))

## -----------------------------------------------------------------------------
overlaps.df <- overlapRegions(A, B)
overlaps.df

## -----------------------------------------------------------------------------
overlaps.df <- overlapRegions(A, B, type="BinA", get.pctA=TRUE)
overlaps.df

## -----------------------------------------------------------------------------
overlaps.df <- overlapRegions(A, B, min.bases=5)
overlaps.df

## -----------------------------------------------------------------------------
overlaps.bool <- overlapRegions(A, B, min.bases=5, only.boolean=TRUE)
overlaps.bool

## -----------------------------------------------------------------------------
overlaps.int <- overlapRegions(A, B, min.bases=5, only.count=TRUE)
overlaps.int

## ---- eval=FALSE--------------------------------------------------------------
#    #NOT RUN
#    set.seed(12345)
#    cpgHMM <- toGRanges("http://www.haowulab.org/software/makeCGI/model-based-cpg-islands-hg19.txt")
#    promoters <- toGRanges("http://gattaca.imppc.org/regioner/data/UCSC.promoters.hg19.bed")

## ---- eval=FALSE--------------------------------------------------------------
#    #NOT RUN
#    cpgHMM <- filterChromosomes(cpgHMM, organism="hg", chr.type="canonical")
#    promoters <- filterChromosomes(promoters, organism="hg", chr.type="canonical")

## ---- eval=FALSE--------------------------------------------------------------
#    #NOT RUN
#    cpgHMM_2K <- sample(cpgHMM, 2000)
#  
#    pt <- overlapPermTest(cpgHMM_2K, promoters, ntimes=1000, genome="hg19", count.once=TRUE)
#    pt
#      P-value: 0.000999000999000999
#      Z-score: 60.5237
#      Number of iterations: 1000
#      Alternative: greater
#      Evaluation of the original region set: 614
#      Evaluation function: numOverlaps
#      Randomization function: randomizeRegions
#    mean(pt$permuted)
#    79.087

## ---- eval=FALSE--------------------------------------------------------------
#    #NOT RUN
#    plot(pt)

## ---- eval=FALSE--------------------------------------------------------------
#    #NOT RUN
#    set.seed(12345)
#    download.file("http://hgdownload-test.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeAwgTfbsUniform/wgEncodeAwgTfbsSydhHepg2Rad21IggrabUniPk.narrowPeak.gz", "Rad21.gz")
#  
#    download.file("http://hgdownload-test.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeAwgTfbsUniform/wgEncodeAwgTfbsUwHepg2CtcfUniPk.narrowPeak.gz", "Ctcf.gz")
#  
#    HepG2_Rad21 <- toGRanges(gzfile("Rad21.gz"), header=FALSE)
#    HepG2_Ctcf <- toGRanges(gzfile("Ctcf.gz"), header=FALSE)
#  
#    promoters <- toGRanges("http://gattaca.imppc.org/regioner/data/UCSC.promoters.hg19.bed")

## ---- eval=FALSE--------------------------------------------------------------
#    #NOT RUN
#    promoters <- filterChromosomes(promoters, organism="hg19")
#    HepG2_Rad21_5K <- sample(HepG2_Rad21, 5000)

## ---- eval=FALSE--------------------------------------------------------------
#    #NOT RUN
#    pt_Rad21_5k_vs_Ctcf <- permTest(A=HepG2_Rad21_5K, B=HepG2_Ctcf, ntimes=1000,
#                                    randomize.function=circularRandomizeRegions,
#                                    evaluate.function=numOverlaps, count.once=TRUE,
#                                    genome="hg19", mc.set.seed=FALSE, mc.cores=4)
#  
#    pt_Rad21_5k_vs_Prom <- permTest(A=HepG2_Rad21_5K, B=promoters, ntimes=1000,
#                                    randomize.function=circularRandomizeRegions,
#                                    evaluate.function=numOverlaps, count.once=TRUE,
#                                    genome="hg19", mc.set.seed=FALSE, mc.cores=4)
#  
#    pt_Rad21_5k_vs_Ctcf
#  
#    pt_Rad21_5k_vs_Prom
#  
#    plot(pt_Rad21_5k_vs_Ctcf, main="Rad21_5K vs CTCF")
#    plot(pt_Rad21_5k_vs_Prom, main="Rad21_5K vs Promoters")

## ---- eval=FALSE--------------------------------------------------------------
#    #NOT RUN
#    lz_Rad21_vs_Ctcf_1 <- localZScore(A=HepG2_Rad21_5K, B=HepG2_Ctcf,
#                                      pt=pt_Rad21_5k_vs_Ctcf,
#                                      window=1000, step=50, count.once=TRUE)
#  
#    lz_Rad21_vs_Prom_1 <- localZScore(A=HepG2_Rad21_5K, B=promoters,
#                                      pt=pt_Rad21_5k_vs_Prom,
#                                      window=1000, step=50, count.once=TRUE)
#  
#    plot(lz_Rad21_vs_Ctcf_1, main="Rad21_5k vs CTCF (1Kbp)")
#    plot(lz_Rad21_vs_Prom_1, main="Rad21_5k vs promoters  (1Kbp)")
#  

## ---- eval=FALSE--------------------------------------------------------------
#    #NOT RUN
#    lz_Rad21_vs_Ctcf_2 <- localZScore(A=HepG2_Rad21_5K, B=HepG2_Ctcf,
#                                      pt=pt_Rad21_5k_vs_Ctcf,
#                                      window=10000, step=500, count.once=TRUE)
#  
#    lz_Rad21_vs_Prom_2 <- localZScore(A=HepG2_Rad21_5K, B=promoters,
#                                      pt=pt_Rad21_5k_vs_Prom,
#                                      window=10000, step=500, count.once=TRUE)
#  
#    plot(lz_Rad21_vs_Ctcf_2, main="Rad21_5k vs CTCF (10Kbp)")
#    plot(lz_Rad21_vs_Prom_2, main="Rad21_5k vs promoters  (10Kbp)")

## ---- sessionInfo-------------------------------------------------------------
sessionInfo()

