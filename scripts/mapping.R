
library(Rsubread)
# https://rockefelleruniversity.github.io/Bioconductor_Introduction/#Overview
# https://www.alzheimersworkbench.ucsd.edu/EndToEndAnalysis_RNASeq.html

#get fq.files and bam files
bamDir<-"bam"
fastqDir<-"../haas/NovaSeq/run96"
unmapped<-"unmapped"

list.files(path = fastqDir)

# create file names
fqFiles<-list.files(fastqDir, "fastq.gz",full.names = T)
fqFileNames<-gsub("\\_(L|[0-9]).+","", fqFiles)
bamFiles<-paste(fqFileNames, ".bam", sep="")
bamFiles<-gsub(fastqDir, bamDir, bamFiles)
toDo<-!file.exists(bamFiles) #only map files that have not yet been mapped
fqFiles<-fqFiles[toDo] #only map files that have not yet been mapped
bamFiles<-bamFiles[toDo] #only map files that have not yet been mapped, code only if in same directory!

#
# QC -----------------------------------------------------------------------------------
# https://www.alzheimersworkbench.ucsd.edu/EndToEndAnalysis_RNASeq.html
library(fastqcr)  # library(Rqc)  # alternative
# http://www.sthda.com/english/wiki/fastqcr-an-r-package-facilitating-quality-controls-of-sequencing-data-for-large-numbers-of-samples
# fastqc_install()  # tool needs to be installed at first use

fastq1 <- list.files(path = file.path(fastqDir), pattern = "*1.fastq.gz$", full.names = TRUE)
fastq2 <- list.files(path = file.path(fastqDir), pattern = "*2.fastq.gz$", full.names = TRUE)
all.equal(length(fastq1),length(fastq2))

fastqc_dir <- "rna-seq-cSBRT/fastqc_dir/"
fastqcr::fastqc(fq.dir = fastqDir,  # FASTQ files directory
                qc.dir = fastqc_dir,  # Results directory
                threads = 4  # default
                )

# Next, we can aggregate all the fastqc results and view them
qc <- qc_aggregate(fastqc_dir)
qc

# Let’s look at what different statistics are available for one file
qc <- qc_read(file.path(fastqc_dir, "file_fastqc.zip"))
names(qc)

# Let’s plot some of the statistics
par(mfrow=c(2,2))
qc_plot(qc, "Per base sequence quality")
qc_plot(qc, "Per sequence quality scores")
qc_plot(qc, "Per base sequence content")
qc_plot(qc, "sequence length distribution")

##
# preprocess raw sequencind data ------------------------------------------
#Using TrimGalore! and MultiQC to preprocess the raw sequencing data
#We’re ready to trim our reads.

rseqR::trim_fastq(fastq1, fastq2, 
                  illumina = TRUE,
                  minlength = 90, 
                  minqual = 30, 
                  trimN = TRUE,
                  retainUnpaired = FALSE,
                  dest.dir = "./TRIMMED_FASTQC",
                  trimgalore = "trim_galore")

#Next, we can run MultiQC on the trimmed fastqc files
run_multiqc(fastqc_dir, "./MULTIQC", multiqc = "multiqc")

#After viewing the MultiQC output, we can move to importing in the trimmed fastq files
trimmed_dir <- "./TRIMMED_FASTQC"
reads1 <- list.files(path = file.path(trimmed_dir), pattern = "*1.fq.gz$", full.names = TRUE)
reads2 <- list.files(path = file.path(trimmed_dir), pattern = "*2.fq.gz$", full.names = TRUE)
all.equal(length(reads1),length(reads2))

#
# align files -----------------------------------------------------------------------------------
# http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/

#runtime 98,8minutes
time1 <- Sys.time()
Rsubread::buildindex(basename = "my_index", reference = "reference_genome/GRCh38/hg38.fa.gz")
time2 <- Sys.time()
difftime <- time2-time1  # Time difference of 2.028082 hours

# align.stat <- Rsubread::align(
#      index="my_index",
#      readfile1=fqFiles[1], 
#      readfile2 = fqFiles[2], 
#      type="rna",
#      input_format = "gzFASTQ",
#      output_format = "BAM",
#      output_file=bamFiles[1],
#      sortReadsByCoordinates=T,
#      phredOffset = 33, # default: numeric value added to base-calling Phred scores to make quality scores 
#      nthreads=64 #use all cores (max 64...)
#    )

# align with subjunc (splice-aware!!)
time1 <- Sys.time()
for(i in seq(from=1,to=length(fqFiles), by=2)){
  j<-i+1
  Rsubread::subjunc(
    index="my_index",
    readfile1=fqFiles[i], 
    readfile2 = fqFiles[j], 
    # type="rna",   ## only in `align()` function
    input_format = "gzFASTQ",
    output_format = "BAM",
    output_file=bamFiles[i],
    sortReadsByCoordinates=T,  # otherwise afterwards with Rsamtools::sortBam() and then indexBAM() -> "SortedActB.bam.bai"
    phredOffset = 33, # default: numeric value added to base-calling Phred scores to make quality scores 
    nthreads=64 #use all cores (max 64...)
  ) #sort to create bigwig Files later
}
time2 <- Sys.time()
difftime <- time2-time1  # Time difference of 5.580616 hours

# We can get an overview of BAM file information using the quickBamFlagSummary() function. Count the number of unaligned reads. How many unique read IDs (QNAMEs) are there?
Rsamtools::quickBamFlagSummary("rna-seq-cSBRT/bam/SID16733.R1.fastq.gz.bam")
# We can now review our BAM file in IGV.

#####
# Read to gene counting with FeatureCounts ----------------------------------------------
# Now that we’ve aligned our reads, we can use FeatureCounts to count the reads to genes from the genomic alignment. We can use the inbuilt annotation.

# count features 
bamFiles <-list.files(path = "rna-seq-cSBRT/bam", pattern = ".bam$", full.names = TRUE)  # needs to be done since bam files do not include "^.R2.fastq.gz.bam" files

countsHg38Paired <- Rsubread::featureCounts(files=bamFiles,
                                          annot.inbuilt="hg38",
                                          annot.ext=NULL,
                                          isGTFAnnotationFile=F,
                                          chrAliases=NULL,
                                          isPairedEnd=T, #check if single or paired end reads are found
                                          countMultiMappingReads=T,
                                          requireBothEndsMapped=T,
                                          nthreads=6)

save(countsHg38Paired, file="rna-seq-cSBRT/rawcounts/rawCounts.rda")
saveRDS(countsHg38Paired, file="rna-seq-cSBRT/rawcounts/rawCounts.rds")
# lets take a look at the stats of our counts

countsHg38Paired$stat
samples$samplename <- colnames(countsHg38Paired$counts)
#samples

