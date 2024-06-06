###############################################################################################################################################

### NOTES: ####

### R codes to process raw DNA sequences from the Novaseq 6000 platform:PE250 (paired end 250bp) with 100K raw tags.
### These codes were customized and compiled based on the DADA2 and ErnakovichLab pipelines (Refer to 'Bioinformaics analysis' in the Method's section and Supplementary Information.)
### Output generates taxonomically assigned ASVs and phyloseq objects.

### Every section labelled 1-9, was split and run in separate chunk on the High Performance Computer servers.
### Section 10-13, could be executed together without loading any workspace if all required objects were saved as .rds files and then loaded with 'readRDS'

###################################################################################################################################################

### Load required packages

library(plyr)
library(data.table)
library(devtools)
library(janitor)
library(splitstackshape)
library(data.table)
library(tidyverse)
library(geosphere)
library(tidyr)
library(stringr)
library(zoo)
library(hms)
library(readxl)
library(readr)
library(ggplot2)
library(markdown)
library(reshape2)
library(reshape)
library(rmarkdown)
library(purrr)
library(Hmisc)
library(checkmate)
library(Formula)
library(agricolae)
library(latticeExtra)
library(htmlTable)
library(ggeasy)
library(forcats)
library(hrbrthemes)
library(viridis)
library(directlabels)
library(vegan)
library(dplyr)
library(lubridate)
library(colorspace)
library(gplots)
library(RColorBrewer)
library(ggrepel)
library(gtools)
library(ComplexHeatmap)
library(dada2); packageVersion("dada2") # the dada2 pipeline
library(ShortRead); packageVersion("ShortRead") # dada2 depends on this
library(cluster)
library(remote)
library(report)
library("phyloseq"); packageVersion("phyloseq")
library(colorRamp2)
library("Biostrings")
library(fantaxtic) # to assign NAs in taxa table as Unassigned and proper labels
library(microbiome)


##############################################################################################################

### 1. Set up pathways to input and output directories, sub-directories and file path for data processing ####

##############################################################################################################

### Add all raw compressed files to one directory
### Set path to directory with raw DNA reads
TMFZ.fq.r1.r2.fp <- "/path_to_directory_with_raw_reads/directory_with_raw_reads"

### List all files in directory and check path
list.files(TMFZ.fq.r1.r2.fp)

### CHANGE FILE PATH & PATTERN depending on file extension
### Match Forward and Reverse fastq filenames format:
fnFs <- sort(list.files(TMFZ.fq.r1.r2.fp, pattern="_1.fq.gz", full.names = TRUE)) 
fnRs <- sort(list.files(TMFZ.fq.r1.r2.fp, pattern="_2.fq.gz", full.names = TRUE))

### Check data: ensure reads are not mixed
fnFs # only contain files with forward reads with format_1.fq.gz 
fnRs # only contain files with reverse reads with format_2.fq.gz

### Extract sample names
### assuming filenames have format: SAMPLENAME_XX_YY.fq.gz. splits filenames using '_' and select the first 3 objects after the split
sample.names <- sapply(strsplit(as.character(basename(fnFs)),"_"), function(x){paste(x[[1]], x[[2]],x[[3]], sep="_")})
sample.names

### Create output directory - TMFZ
all.processed.fp <- "/path_to_output_directory/TMFZ"
dir.create(all.processed.fp)

### Set path to 01_preprocess which will be created in TMFZ later
preprocess.fp <- file.path(all.processed.fp, "01_preprocess") 

### Set path to filtN which will be created in TMFZ/01_preprocess/ later
filtN.fp <- file.path(preprocess.fp, "filtN")

### Set path to trimmed.cutadapt which will be created in TMFZ/01_preprocess/ later
trimmed.fp <- file.path(preprocess.fp, "trimmed.cutadapt") 

### Set path to 02_filter which will be created in TMFZ folder later
filter.fp <- file.path(all.processed.fp, "02_filter")  

### Set path to 03_tabletax which will be created in TMFZ folder later
table.fp <- file.path(all.processed.fp, "03_tabletax") 


### save R environment 
### Saving in home directory or directory from where the code is submitted to the server (change path to file as needed)
save.image(file = "TMFZ_Renv.RData")

###############################################################################################################################

### 2. Pre-filter to remove sequence reads with Ns ####

##############################################################################################################################

### DADA2 does not allow Ns. 
### Ambiguous bases will make it hard for cutadapt to find short primer sequences in the reads  - Section 3. 
### To solve this problem, remove sequences with ambiguous bases (Ns)

### Sequences with more than maxN Ns will be discarded. 
### And will be saved into the filtN directory 
### This creates the 01_preprocess and filtN directory and outputs N-trimmed R1 and R2 reads 
### DADA2 filtering is performed on the forward and reverse reads independently, and both reads must pass for the read pair to be output

### Load R environment saved at the end of section 1 if this section is executed in a separate script
## if not then there is no need to reload the R environment
load(file = "TMFZ_Renv.RData")

### DO NOT CHANGE
### Name the N-filtered files to be saved in filtN/ subdirectory
### Only provides the folder path- does not create the 01_preprocess and filtN subdirectories
fnFs.filtN <- file.path(preprocess.fp, "filtN", basename(fnFs))
fnRs.filtN <- file.path(preprocess.fp, "filtN", basename(fnRs))

### DO NOT CHANGE
### multithread = 6: using less processes than cores so that each has enough memory
pre_out <- filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = 6) 

### DO NOT CHANGE
### filt_out returns 2 columnns: reads_in and reads_our
pre_out

### check if fastq files differed before and after N trimming on Linux/Ubuntu terminal using:
### diff -arq /path_to_directory_with_raw_reads/directory_with_raw_reads /path_to_output_directory/TMFZ/01_preprocess/filtN
### both the forward and reverse files appear to differ before and after N trimming on terminal

### save R environment
save.image(file = "TMFZ_Renv.RData")

#####################################################################################################################################

### 3. Remove primers using Cutadapt ####

############################################################################################################################

### Load R environment saved at the end of section 2 if this section is executed in a separate script
## if not then there is no need to reload the R environment
load(file = "TMFZ_Renv.RData")

### Trim primer sequences using cutadapt
### PERFORM CUTADPT testing for every sample set######

### Download and Set up pathway to cutadapt
cutadapt <- "cutadapt" 
system2(cutadapt, args = "--version") # Check by running shell command from R

### CHANGE If using different primers
### Set up primer sequences to use with function to detect primers

### this is 515F - Forward primer
FWD <- "GTGCCAGCMGCCGCGGTAA" 

### this is 806R - Reverse primer
REV <- "GGACTACHVGGGTWTCTAAT" 

### DO NOT CHANGE : Function that creates a list of all orientations of the primers

allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer) # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna),
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString)) # Convert back to character vector
}

### DO NOT CHANGE
### To cater for reads in different orientations: Save primers in fwd and reverse compliment orientations as well
### variables stored as 'orients' are used with function to detect primers 
FWD.orients <- allOrients(FWD) 
REV.orients <- allOrients(REV)

### Get reverse compliments of forward and reverse primers
FWD2 <- FWD.orients[["RevComp"]]
REV2 <- REV.orients[["RevComp"]]

FWD.orients.2 <- allOrients(FWD2)
REV.orients.2 <- allOrients(REV2)

### sanity check
FWD.orients
REV.orients
FWD.orients.2
REV.orients.2

### DO NOT CHANGE
### Function to count the number of times the primers appear in the forward and reverse read, while considering all possible primer orientations. 

primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}

### Before running cutadapt, perform primer detection in all samples

rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN),
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN),
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN),
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN))

rbind(FWD.ForwardReads = sapply(FWD.orients.2, primerHits, fn = fnFs.filtN),
      FWD.ReverseReads = sapply(FWD.orients.2, primerHits, fn = fnRs.filtN),
      REV.ForwardReads = sapply(REV.orients.2, primerHits, fn = fnFs.filtN),
      REV.ReverseReads = sapply(REV.orients.2, primerHits, fn = fnRs.filtN))


### Save the reverse complements of the primers to variables
FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)

### Save the reverse complements of the primers to variables
FWD.RC.2 <- dada2:::rc(FWD2) #FWD primer
REV.RC.2 <- dada2:::rc(REV2) #Rev primer

### DO NOT CHANGE
### Create the cutadapt flags
### Using only R1.flags and R2.flags with cutadapt, resulted in some reads with untrimmed primers.
### So use all 4 flags to trim away primers from reads in mixed orientations.

### Trim FWD and the reverse-complement of REV off of R1 (forward reads)
### -g:Sequence of an adapter ligated to the 5' end which is fwd primer sequence found as exact match on 5' of R1
### -a:Sequence of an adapter ligated to the 3' end which is the rev complement of the reverse primer sequence found as exact match on 3' of R1
R1.flags <- paste("-g", FWD, "-a", REV.RC)

### Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
### -G :5' adapter to be removed from R2 which is reverse primer sequence found as exact match on 5' of R2
### -A: 3' adapter to be removed from R2 which is the same as the rev complement of the forward primer sequence found as exact match on 3' of R2
R2.flags <- paste("-G", REV, "-A", FWD.RC)

### Trim FWD2 (reverse complement of FWD) and the REV off of R1 (forward reads)
### -g:Sequence of an adapter ligated to the 5' end which is fwd primer sequence found as exact match on 5' of R1
### -a:Sequence of an adapter ligated to the 3' end which is the rev complement of the reverse primer sequence found as exact match on 3' of R1
R1.flags.2 <- paste("-g", FWD2, "-a", REV.RC.2)

### Trim REV2 (reverse complement of REV) and the FWD off of R2 (reverse reads)
### -G :5' adapter to be removed from R2 which is reverse primer sequence found as exact match on 5' of R2
### -A: 3' adapter to be removed from R2 which is the same as the rev complement of the forward primer sequence found as exact match on 3' of R2
R2.flags.2 <- paste("-G", REV2, "-A", FWD.RC.2)

### DO NOT CHANGE
### Create directory to hold the output from cutadapt. This creates the trimmed.cutadapt folder in 01_preprocess
if (!dir.exists(trimmed.fp)) dir.create(trimmed.fp)

### DO NOT CHANGE
### prepare files for cutadapt output
fnFs.cut <- file.path(trimmed.fp, basename(fnFs)) 
fnRs.cut <- file.path(trimmed.fp, basename(fnRs))

### Run Cutadapt:

for (i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags,R1.flags.2, R2.flags.2, "-n", 5, # -n 2 required to remove FWD and REV from reads
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             "--minimum-length", 150, #reads below 100 base pairs are removed
                             "-e", 0, # reduce maximum allowed error rate to 0 instead of default 0.1%
                             "--report", "minimal", # one line summary
                             "-j", 0,	
                             fnFs.filtN[i], fnRs.filtN[i])) # input files
}

### "-j", 0, # automatically detects the number of available cores
### "--report=minimal", # one line summary
### j : The detection takes into account computational resource restrictions that may be in place. 
### For example, if running Cutadapt as a batch job on a cluster system, the actual number of cores assigned to the job will be used.

### DO NOT CHANGE
### run sanity check 1

rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut),
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut),
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut),
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut))

### sanity check 2

rbind(FWD.ForwardReads = sapply(FWD.orients.2, primerHits, fn = fnFs.cut),
      FWD.ReverseReads = sapply(FWD.orients.2, primerHits, fn = fnRs.cut),
      REV.ForwardReads = sapply(REV.orients.2, primerHits, fn = fnFs.cut),
      REV.ReverseReads = sapply(REV.orients.2, primerHits, fn = fnRs.cut))


### All primer counts should be 0. check the output to make sure that there are no primers remaining in the samples.


### Save the R environment
save.image(file = "TMFZ_Renv.RData")

########################################################################################################################################

### 4. Inspect read quality profile: visualizing the quality profiles of the forward and reverse N-filtered & cutadapt trimmed reads ####

######################################################################################################################################

### Raw reads often have regions with low-quality reads. 
### Most Illumina sequencing data shows a trend of decreasing average quality towards the end of sequencing reads (Ben J. Callahan 2016). 
### To know which regions have low quality reads, plot it using plotQualityProfile() from DADA2 package.
### This function plots a visual summary of the distribution of quality scores for each sequence position.

### Plot description: The quality profile plot is a gray-scale heatmap of the frequency of each quality score at each base position.
### The mean quality score at each position is shown by the green line
### The median quality score is shown by the solid orange line 
### The quartiles of the quality score distribution by the dashed orange lines.
### The red line shows the scaled proportion of reads that extend to at least that position-(this is more useful for other sequencing technologies, 
### as Illumina reads are typically all the same lenghth, hence the flat red line).

### Load R environment saved at the end of section 3 if this section is executed in a separate script
## if not then there is no need to reload the R environment
load(file = "TMFZ_Renv.RData")

### For all further filtering after cutadapt use the processed/02_filter.
### Create file paths and subdirectories in 02_filter to separate forward and reverse from trimmed(cutadapt) object fnFs.cut 

### DO NOT CHANGE
### creates /path_to_output_directory/TMFZ/02_filter/
dir.create(filter.fp)

### creates path only 
subF.fp <- file.path(filter.fp, "preprocessed_F") 
subR.fp <- file.path(filter.fp, "preprocessed_R") 

### sanity check
subF.fp
subR.fp

### creates 02_filter/preprocessed_F
dir.create(subF.fp) 

### creates 02_filter/preprocessed_R
dir.create(subR.fp) 


### Now copy R1 (forward) and R2(reverse) reads from trimmed(cutadapt) object fnFs.cut into separate sub directories

### Creates path for files to be added to 02_filter/preprocessed_F
fnFs.Q <- file.path(subF.fp, basename(fnFs))

### Creates path for files to be added to 02_filter/preprocessed_R
fnRs.Q <- file.path(subR.fp, basename(fnRs))

### Create a copy of separate forward and reverse cutadapt trimmed reads in 02_filter as preprocessed_f and preprocessed_f
file.copy(from = fnFs.cut, to = fnFs.Q)
file.copy(from = fnRs.cut, to = fnRs.Q)

### List files 
fastqFs <- sort(list.files(subF.fp, pattern="_1.fq.gz",full.names = TRUE)) 
fastqRs <- sort(list.files(subR.fp, pattern="_2.fq.gz",full.names = TRUE)) 
if(length(fastqFs) != length(fastqRs)) stop("Forward and reverse files do not match.")

### sanity check
fastqFs 
fastqRs

### Chose random samples
rand_samples <- sample(size = 12, 1:length(fastqFs)) # grab 12 random samples to plot 

### check chosen samples
rand_samples

### Plot quality scores of forward and reverse reads
fwd.qual.plot <- plotQualityProfile(fastqFs[rand_samples]) + labs(x = "Sequence Position") 
rev.qual.plot <- plotQualityProfile(fastqRs[rand_samples]) + labs(x = "Sequence Position")

### Warning: 'package:stats' may not be available when loading is considered a bug! : IGNORE WARNING

### Save plots

pdf(file="/path_to_output_directory/TMFZ/fwd.qual.plot.pdf", width = 15, height = 8.5) 
fwd.qual.plot 
dev.off() 

pdf(file="/path_to_output_directory/TMFZ/rev.qual.plot.pdf", width = 15, height = 8.5) 
rev.qual.plot 
dev.off()

### Write plots to disk
saveRDS(fwd.qual.plot, "/path_to_output_directory/TMFZ/fwd.qual.plot.rds")
saveRDS(rev.qual.plot, "/path_to_output_directory/TMFZ/rev.qual.plot.rds")


### Save the R environment
save.image(file = "TMFZ_Renv.RData")

#########################################################################################################################

### 5. Trim & Filter reads ####

#########################################################################################################################

### Load R environment saved at the end of section 4 if this section is executed in a separate script
## if not then there is no need to reload the R environment
load(file = "TMFZ_Renv.RData")

### File parsing: create file path for filtered and trimmed reads separately.
### Create path for subdirectories : files go into "/path_to_output_directory/TMFZ/02_filter/"
filtpathF <- file.path(filter.fp, "filtered_F") 
filtpathR <- file.path(filter.fp, "filtered_R") 

dir.create(filtpathF)
dir.create(filtpathR)

### Creates path for files to be added to 02_filter/filtered_F
filtpathF.file <- file.path(filtpathF, basename(fnFs))

### creates path for files to be added to 02_filter/filtered_R
filtpathR.file <- file.path(filtpathR, basename(fnRs))

filt_out <- filterAndTrim(fastqFs, filtpathF.file, fastqRs, filtpathR.file,
                          maxEE=c(2,2), truncQ=2, maxN=0, truncLen=c(200,200),rm.phix=TRUE,
                          compress=FALSE, verbose=TRUE, multithread=24,matchIDs=TRUE)

### Although only sample names of the forward reads are shown both forward and reverse reads are filtered.
filt_out

### summary of samples in filt_out by percentage - shows an average of 96% filtered paired sequences

filt_out_summary <- filt_out %>% data.frame() %>% 
  mutate(Samples = rownames(.), percent_kept = 100*(reads.out/reads.in)) %>% 
  select(Samples, everything())

filt_out_summary

### Plot the quality of the filtered fastq files
### Plot 12 randomly selected samples

### list out filtered files
filtFs <- sort(list.files(filtpathF, pattern="_1.fq.gz",full.names = TRUE)) 
filtRs <- sort(list.files(filtpathR, pattern="_2.fq.gz",full.names = TRUE))

### sanity check
filtFs
filtRs  

### Grab 12 random samples to plot
rand_samples_filt <- sample(size = 12, 1:length(filtFs))  
fwd.fil.plot <- plotQualityProfile(filtFs[rand_samples_filt]) + labs(x = "Sequence Position")
rev.fil.plot <- plotQualityProfile(filtRs[rand_samples_filt]) + labs(x = "Sequence Position")

### Save plots

pdf(file="/path_to_output_directory/TMFZ/fwd.fil.plot.pdf", width = 15, height = 8.5) 
fwd.fil.plot 
dev.off() 

pdf(file="/path_to_output_directory/TMFZ/rev.fil.plot.pdf", width = 15, height = 8.5) 
rev.fil.plot
dev.off()

### write plots to disk
saveRDS(fwd.fil.plot, "/path_to_output_directory/TMFZ/fwd.fil.plot.rds")
saveRDS(rev.fil.plot, "/path_to_output_directory/TMFZ/rev.fil.plot.rds")


### Save the R environment
save.image(file = "TMFZ_Renv.RData")

#########################################################################################################################

### 6. Learn Error Rates to Infer sequence variants - test 5 different error models ####

########################################################################################################################

### Test 4 other error models apart from the traditional DADA2 error model to estimate errors in Novaseq reads with 
### “binned” quality scores

### INTERPRETTING ERROR PLOTS AND CHOOSING THE BEST ERROR MODEL:
### The error rates for each possible transition (A→C, A→G, …) are shown. 
### Points are the observed error rates for each consensus quality score. 
### The black line shows the estimated error rates after convergence of the machine-learning algorithm. 
### The red line shows the error rates expected under the nominal definition of the Q-score. 
### Choose model with:
### 1. good fit to the observed rates (plots that have points that mostly align with the black lines) and 
### 2. decreasing error rates with increased quality (the black line is continuously decreasing)

### DADA2 authors suggest implementing an enforced monotonicity option in the main package which aligns with Option2 
### https://github.com/benjjneb/dada2/issues/791, option 2 in github - https://github.com/benjjneb/dada2/issues/1307


### Load R environment saved at the end of section 5 if this section is executed in a separate script
## if not then there is no need to reload the R environment
load(file = "TMFZ_Renv.RData")

### Recall and cross-check file listing
filtFs
filtRs

### check if sample names are in order
sample.namesF <- basename(filtFs) # doesn't drop fastq.gz 
sample.namesF

sample.namesF <- gsub("_1.fq.gz", "", sample.namesF) 
sample.namesF

sample.namesR <- basename(filtRs) #doesn't drop fastq.gz 
sample.namesR

sample.namesR <- gsub("_2.fq.gz", "", sample.namesR)
sample.namesR

# Double check
if(!identical(sample.namesF, sample.namesR)) stop("Forward and reverse files do not match.") 

names(filtFs) <- sample.namesF
names(filtRs) <- sample.namesR

names(filtFs)
names(filtRs)

### Run error models

### Traditional way of learning error rates in DADA2 - learn error rates 1
### The DADA2 algorithm makes use of a parametric error model (err) 

### set seed to ensure that randomized steps are replicatable
set.seed(100) 

### nbases = 1e10 uses 10018061600 total bases in 50090308 reads from 384 samples will be used for learning the error rates in fwd reads and
###                    10016016400 total bases in 50080082 reads from 382 samples will be used for learning the error rates in reverse reads

errF <- learnErrors(filtFs, nbases = 1e10, multithread = 24, randomize = TRUE)
errR <- learnErrors(filtRs, nbases = 1e10, multithread = 24, randomize = TRUE)

### Plot errors

### In such large datasets, a warning message received is:
### Warning messages:
#   1: Transformation introduced infinite values in continuous y-axis
#   2: Transformation introduced infinite values in continuous y-axis
### From Ben at DADA2: 
#   This isn't an error, just a message from the plotting function to let you know that there were some zero values in the data plotted 
#   (which turn into infinities on the log-scale). That is completely expected, it results from the fact that not every combination of error type 
#   (e.g. A->C) and quality score (e.g. 33) is observed in your data which is normal.


pdf(file="/path_to_output_directoryTMFZ/errF.DADA2.pdf", width = 15, height = 8.5) 
plotErrors(errF, nominalQ=TRUE) 
dev.off() 

pdf(file="/path_to_output_directory/TMFZ/errR.DADA2.pdf", width = 15, height = 8.5)
plotErrors(errR, nominalQ=TRUE) 
dev.off()

saveRDS(errF, "/path_to_output_directory/TMFZ/errF.DADA2.rds")
saveRDS(errR, "/path_to_output_directory/TMFZ/errR.DADA2.rds")

### Four options for learning error rates with NovaSeq data####

### Option 1 from JacobRPrice alter loess arguments (weights and span and enforce monotonicity) benjjneb/dada2#1307####

loessErrfun_mod1 <- function(trans) {
  qq <- as.numeric(colnames(trans))
  est <- matrix(0, nrow=0, ncol=length(qq))
  for(nti in c("A","C","G","T")) {
    for(ntj in c("A","C","G","T")) {
      if(nti != ntj) {
        errs <- trans[paste0(nti,"2",ntj),]
        tot <- colSums(trans[paste0(nti,"2",c("A","C","G","T")),])
        rlogp <- log10((errs+1)/tot) # 1 psuedocount for each err, but if tot=0 will give NA
        rlogp[is.infinite(rlogp)] <- NA
        df <- data.frame(q=qq, errs=errs, tot=tot, rlogp=rlogp)
        
        # original
        # ###! mod.lo <- loess(rlogp ~ q, df, weights=errs) ###!
        # mod.lo <- loess(rlogp ~ q, df, weights=tot) ###!
        # #        mod.lo <- loess(rlogp ~ q, df)
        
        # Gulliem Salazar's solution https://github.com/benjjneb/dada2/issues/938
        mod.lo <- loess(rlogp ~ q, df, weights = log10(tot),span = 2)
        
        pred <- predict(mod.lo, qq)
        maxrli <- max(which(!is.na(pred)))
        minrli <- min(which(!is.na(pred)))
        pred[seq_along(pred)>maxrli] <- pred[[maxrli]]
        pred[seq_along(pred)<minrli] <- pred[[minrli]]
        est <- rbind(est, 10^pred)
      } # if(nti != ntj)
    } # for(ntj in c("A","C","G","T"))
  } # for(nti in c("A","C","G","T"))
  
  # HACKY
  MAX_ERROR_RATE <- 0.25
  MIN_ERROR_RATE <- 1e-7
  est[est>MAX_ERROR_RATE] <- MAX_ERROR_RATE
  est[est<MIN_ERROR_RATE] <- MIN_ERROR_RATE
  
  # enforce monotonicity
  # https://github.com/benjjneb/dada2/issues/791
  estorig <- est
  est <- est %>%
    data.frame() %>%
    mutate_all(funs(case_when(. < X40 ~ X40,
                              . >= X40 ~ .))) %>% as.matrix()
  rownames(est) <- rownames(estorig)
  colnames(est) <- colnames(estorig)
  
  # Expand the err matrix with the self-transition probs
  err <- rbind(1-colSums(est[1:3,]), est[1:3,],
               est[4,], 1-colSums(est[4:6,]), est[5:6,],
               est[7:8,], 1-colSums(est[7:9,]), est[9,],
               est[10:12,], 1-colSums(est[10:12,]))
  rownames(err) <- paste0(rep(c("A","C","G","T"), each=4), "2", c("A","C","G","T"))
  colnames(err) <- colnames(trans)
  # Return
  return(err)
}

### check what this looks like

errF_1 <- learnErrors(
  filtFs,
  multithread = TRUE,
  nbases = 1e10,
  errorEstimationFunction = loessErrfun_mod1,
  verbose = TRUE ) 

errR_1 <- learnErrors(
  filtRs,
  multithread = TRUE,
  nbases = 1e10,
  errorEstimationFunction = loessErrfun_mod1,
  verbose = TRUE )

### Error inferrence plots

pdf(file="/path_to_output_directoryTMFZ/errF_1.pdf", width = 15, height = 8.5)
plotErrors(errF_1, nominalQ=TRUE)
dev.off()

pdf(file="/path_to_output_directoryTMFZ/errR_1.pdf", width = 15, height = 8.5)
plotErrors(errR_1, nominalQ=TRUE)
dev.off()

saveRDS(errF_1, "/path_to_output_directory/TMFZ/errF_1.rds")
saveRDS(errR_1, "/path_to_output_directory/TMFZ/errR_1.rds")


### Option 2 enforce monotonicity only.Originally recommended by DADA2 developers in: benjjneb/dada2#791 ####

loessErrfun_mod2 <- function(trans) {
  qq <- as.numeric(colnames(trans))
  est <- matrix(0, nrow=0, ncol=length(qq))
  for(nti in c("A","C","G","T")) {
    for(ntj in c("A","C","G","T")) {
      if(nti != ntj) {
        errs <- trans[paste0(nti,"2",ntj),]
        tot <- colSums(trans[paste0(nti,"2",c("A","C","G","T")),])
        rlogp <- log10((errs+1)/tot)  # 1 psuedocount for each err, but if tot=0 will give NA
        rlogp[is.infinite(rlogp)] <- NA
        df <- data.frame(q=qq, errs=errs, tot=tot, rlogp=rlogp)
        
        # original
        # ###! mod.lo <- loess(rlogp ~ q, df, weights=errs) ###!
        mod.lo <- loess(rlogp ~ q, df, weights=tot) ###!
        # #        mod.lo <- loess(rlogp ~ q, df)
        
        # Gulliem Salazar's solution
        # https://github.com/benjjneb/dada2/issues/938
        # mod.lo <- loess(rlogp ~ q, df, weights = log10(tot),span = 2)
        
        pred <- predict(mod.lo, qq)
        maxrli <- max(which(!is.na(pred)))
        minrli <- min(which(!is.na(pred)))
        pred[seq_along(pred)>maxrli] <- pred[[maxrli]]
        pred[seq_along(pred)<minrli] <- pred[[minrli]]
        est <- rbind(est, 10^pred)
      } # if(nti != ntj)
    } # for(ntj in c("A","C","G","T"))
  } # for(nti in c("A","C","G","T"))
  
  # HACKY
  MAX_ERROR_RATE <- 0.25
  MIN_ERROR_RATE <- 1e-7
  est[est>MAX_ERROR_RATE] <- MAX_ERROR_RATE
  est[est<MIN_ERROR_RATE] <- MIN_ERROR_RATE
  
  # enforce monotonicity
  # https://github.com/benjjneb/dada2/issues/791
  estorig <- est
  est <- est %>%
    data.frame() %>%
    mutate_all(funs(case_when(. < X40 ~ X40,
                              . >= X40 ~ .))) %>% as.matrix()
  rownames(est) <- rownames(estorig)
  colnames(est) <- colnames(estorig)
  
  # Expand the err matrix with the self-transition probs
  err <- rbind(1-colSums(est[1:3,]), est[1:3,],
               est[4,], 1-colSums(est[4:6,]), est[5:6,],
               est[7:8,], 1-colSums(est[7:9,]), est[9,],
               est[10:12,], 1-colSums(est[10:12,]))
  rownames(err) <- paste0(rep(c("A","C","G","T"), each=4), "2", c("A","C","G","T"))
  colnames(err) <- colnames(trans)
  # Return
  return(err)
}


### check what this looks like

errF_2 <- learnErrors(
  filtFs,
  multithread = TRUE,
  nbases = 1e10,
  errorEstimationFunction = loessErrfun_mod2,
  verbose = TRUE
)

errR_2 <- learnErrors(
  filtRs,
  multithread = TRUE,
  nbases = 1e10,
  errorEstimationFunction = loessErrfun_mod2,
  verbose = TRUE
)

### Plot errors

pdf(file="/path_to_output_directory/TMFZ/errF_2.pdf", width = 15, height = 8.5)
plotErrors(errF_2, nominalQ=TRUE)
dev.off()

pdf(file="/path_to_output_directory/TMFZ/errR_2.pdf", width = 15, height = 8.5)
plotErrors(errR_2, nominalQ=TRUE)
dev.off()

saveRDS(errF_2, "/path_to_output_directory/TMFZ/errF_2.rds")
saveRDS(errR_2, "/path_to_output_directory/TMFZ/errR_2.rds")


### Option 3 alter loess function (weights only) and enforce monotonicity: From JacobRPrice benjjneb/dada2#1307 ####

loessErrfun_mod3 <- function(trans) {
  qq <- as.numeric(colnames(trans))
  est <- matrix(0, nrow=0, ncol=length(qq))
  for(nti in c("A","C","G","T")) {
    for(ntj in c("A","C","G","T")) {
      if(nti != ntj) {
        errs <- trans[paste0(nti,"2",ntj),]
        tot <- colSums(trans[paste0(nti,"2",c("A","C","G","T")),])
        rlogp <- log10((errs+1)/tot)  # 1 psuedocount for each err, but if tot=0 will give NA
        rlogp[is.infinite(rlogp)] <- NA
        df <- data.frame(q=qq, errs=errs, tot=tot, rlogp=rlogp)
        
        # original
        # ###! mod.lo <- loess(rlogp ~ q, df, weights=errs) ###!
        # mod.lo <- loess(rlogp ~ q, df, weights=tot) ###!
        # #        mod.lo <- loess(rlogp ~ q, df)
        
        # Gulliem Salazar's solution
        # https://github.com/benjjneb/dada2/issues/938
        # mod.lo <- loess(rlogp ~ q, df, weights = log10(tot),span = 2)
        
        
        # only change the weights
        mod.lo <- loess(rlogp ~ q, df, weights = log10(tot))
        
        pred <- predict(mod.lo, qq)
        maxrli <- max(which(!is.na(pred)))
        minrli <- min(which(!is.na(pred)))
        pred[seq_along(pred)>maxrli] <- pred[[maxrli]]
        pred[seq_along(pred)<minrli] <- pred[[minrli]]
        est <- rbind(est, 10^pred)
      } # if(nti != ntj)
    } # for(ntj in c("A","C","G","T"))
  } # for(nti in c("A","C","G","T"))
  
  # HACKY
  MAX_ERROR_RATE <- 0.25
  MIN_ERROR_RATE <- 1e-7
  est[est>MAX_ERROR_RATE] <- MAX_ERROR_RATE
  est[est<MIN_ERROR_RATE] <- MIN_ERROR_RATE
  
  # enforce monotonicity
  # https://github.com/benjjneb/dada2/issues/791
  estorig <- est
  est <- est %>%
    data.frame() %>%
    mutate_all(funs(case_when(. < X40 ~ X40,
                              . >= X40 ~ .))) %>% as.matrix()
  rownames(est) <- rownames(estorig)
  colnames(est) <- colnames(estorig)
  
  # Expand the err matrix with the self-transition probs
  err <- rbind(1-colSums(est[1:3,]), est[1:3,],
               est[4,], 1-colSums(est[4:6,]), est[5:6,],
               est[7:8,], 1-colSums(est[7:9,]), est[9,],
               est[10:12,], 1-colSums(est[10:12,]))
  rownames(err) <- paste0(rep(c("A","C","G","T"), each=4), "2", c("A","C","G","T"))
  colnames(err) <- colnames(trans)
  # Return
  return(err)
}

### check what this looks like

errF_3 <- learnErrors(
  filtFs,
  multithread = TRUE,
  nbases = 1e10,
  errorEstimationFunction = loessErrfun_mod3,
  verbose = TRUE
)

errR_3 <- learnErrors(
  filtRs,
  multithread = TRUE,
  nbases = 1e10,
  errorEstimationFunction = loessErrfun_mod3,
  verbose = TRUE
)

### Plot errors

pdf(file="/path_to_output_directory/TMFZ/errF_3.pdf", width = 15, height = 8.5)
plotErrors(errF_3, nominalQ=TRUE)
dev.off()

pdf(file="/path_to_output_directory/TMFZ/errR_3.pdf", width = 15, height = 8.5)
plotErrors(errR_3, nominalQ=TRUE)
dev.off()

saveRDS(errF_3, "/path_to_output_directory/TMFZ/errF_3.rds")
saveRDS(errR_3, "/path_to_output_directory/TMFZ/errR_3.rds")


### Option 4 Alter loess function arguments (weights and span and degree, also enforce monotonicity):From Jonalim’s comment in benjjneb/dada2#1307 ####

loessErrfun_mod4 <- function(trans) {
  qq <- as.numeric(colnames(trans))
  est <- matrix(0, nrow=0, ncol=length(qq))
  for(nti in c("A","C","G","T")) {
    for(ntj in c("A","C","G","T")) {
      if(nti != ntj) {
        errs <- trans[paste0(nti,"2",ntj),]
        tot <- colSums(trans[paste0(nti,"2",c("A","C","G","T")),])
        rlogp <- log10((errs+1)/tot)  # 1 psuedocount for each err, but if tot=0 will give NA
        rlogp[is.infinite(rlogp)] <- NA
        df <- data.frame(q=qq, errs=errs, tot=tot, rlogp=rlogp)
        
        # original
        # ###! mod.lo <- loess(rlogp ~ q, df, weights=errs) ###!
        # mod.lo <- loess(rlogp ~ q, df, weights=tot) ###!
        # #        mod.lo <- loess(rlogp ~ q, df)
        
        # jonalim's solution
        # https://github.com/benjjneb/dada2/issues/938
        mod.lo <- loess(rlogp ~ q, df, weights = log10(tot),degree = 1, span = 0.95)
        
        pred <- predict(mod.lo, qq)
        maxrli <- max(which(!is.na(pred)))
        minrli <- min(which(!is.na(pred)))
        pred[seq_along(pred)>maxrli] <- pred[[maxrli]]
        pred[seq_along(pred)<minrli] <- pred[[minrli]]
        est <- rbind(est, 10^pred)
      } # if(nti != ntj)
    } # for(ntj in c("A","C","G","T"))
  } # for(nti in c("A","C","G","T"))
  
  # HACKY
  MAX_ERROR_RATE <- 0.25
  MIN_ERROR_RATE <- 1e-7
  est[est>MAX_ERROR_RATE] <- MAX_ERROR_RATE
  est[est<MIN_ERROR_RATE] <- MIN_ERROR_RATE
  
  # enforce monotonicity
  # https://github.com/benjjneb/dada2/issues/791
  estorig <- est
  est <- est %>%
    data.frame() %>%
    mutate_all(funs(case_when(. < X40 ~ X40,
                              . >= X40 ~ .))) %>% as.matrix()
  rownames(est) <- rownames(estorig)
  colnames(est) <- colnames(estorig)
  
  # Expand the err matrix with the self-transition probs
  err <- rbind(1-colSums(est[1:3,]), est[1:3,],
               est[4,], 1-colSums(est[4:6,]), est[5:6,],
               est[7:8,], 1-colSums(est[7:9,]), est[9,],
               est[10:12,], 1-colSums(est[10:12,]))
  rownames(err) <- paste0(rep(c("A","C","G","T"), each=4), "2", c("A","C","G","T"))
  colnames(err) <- colnames(trans)
  # Return
  return(err)
}

### check what this looks like

errF_4 <- learnErrors(
  filtFs,
  multithread = TRUE,
  nbases = 1e10,
  errorEstimationFunction = loessErrfun_mod4,
  verbose = TRUE
)

errR_4 <- learnErrors(
  filtRs,
  multithread = TRUE,
  nbases = 1e10,
  errorEstimationFunction = loessErrfun_mod4,
  verbose = TRUE
)

### Plot errors

pdf(file="/path_to_output_directory/TMFZ/errF_4.pdf", width = 15, height = 8.5)
plotErrors(errF_4, nominalQ=TRUE)
dev.off()

pdf(file="/path_to_output_directory/TMFZ/errR_4.pdf", width = 15, height = 8.5)
plotErrors(errR_4, nominalQ=TRUE)
dev.off()

saveRDS(errF_4, "/path_to_output_directory/TMFZ/errF_4.rds")
saveRDS(errR_4, "/path_to_output_directory/TMFZ/errR_4.rds")


### Save the R environment
save.image(file = "TMFZ_Renv.RData")

################################################################################################################

### 7. Dereplication ####

##########################################################################################################################################################

### Dereplication combines identical sequencing reads into into “unique sequences” with a corresponding “abundance” equal to the number of reads 
### with that unique sequence. 

### Dereplication in the DADA2 pipeline has one crucial addition from other pipelines: 
#   DADA2 retains a summary of the quality information associated with each unique sequence. 
#   The consensus quality profile of a unique sequence is the average of the positional qualities from the dereplicated reads. These quality profiles 
#   informs the error model of the subsequent sample inference step, significantly increasing DADA2’s accuracy.

### Load R environment saved at the end of section 6 if this section is executed in a separate script
## if not then there is no need to reload the R environment
load(file = "TMFZ_Renv.RData")

### Dereplicate forward and reverse reads
derepF <- derepFastq(filtFs, verbose = TRUE) 
derepR <- derepFastq(filtRs, verbose = TRUE)

### sanity check on output 
### outputs the summary of dereplication for each file 
derepF 
derepR

### Name the derep-class objects by the sample names
names(derepF) <- sample.namesF 
names(derepR) <- sample.namesR


### Save the R environment
save.image(file = "TMFZ_Renv.RData")

########################################################################################################################

### 8. Infer sequence variants with Option 3 alter loess function (weights only) and enforce monotonicity ####

##########################################################################################################################################################

### Assigning sequences to ASVs (called “sequence inference”). 

### At this step, the core sample inference algorithm is applied to the dereplicated data.

### By default, the dada function processes each sample independently:
#   DADA2 resolves amplicon sequence variants exactly, and because exact DNA sequences are consistent labels, samples can be processed independently by
#   DADA2 and then combined, and this is the default behavior.

### Independent sample processing has two major advantages: 
#   Computation time is linear in the number of samples, and memory requirements are flat with the number of samples. 
#   However, pooling allows information to be shared across samples, which makes it easier to resolve rare variants that were 
#   present as singletons or doubletone in one sample but were present many times across samples. 


### https://benjjneb.github.io/dada2/pool.html :
#   De novo methods must pool samples before processing them, as without pooling the labels between samples are not consistent and cannot be compared,
#   i.e. OTU1 in sample 1 and OTU1 sample 2 won’t be the same.

### However, pooling information across samples can increase sensitivity to sequence variants that may be present at very low frequencies in 
#   multiple samples.

### The dada2 package offers two types of pooling:
#   dada(..., pool=TRUE) performs standard pooled processing, in which all samples are pooled together for sample inference.
#   dada(..., pool="pseudo") performs pseudo-pooling, in which samples are processed independently after sharing information between samples,
#   approximating pooled sample inference in linear time.

### The parameter pool = can be set to: pool = FALSE (default), pool = TRUE, or pool = psuedo 

### pool = FALSE: Sequence information is not shared between samples. Fast processing time, less sensitivity to rare taxa. 
### pool = pseudo: Sequence information is shared in a separate “prior” step. Intermediate processing time, intermediate sensitivity to rare taxa. 
### pool = TRUE: Sequence information from all samples is pooled together. Slow processing time, most sensitivity to rare taxa.

### We use pool = pseudo and error model for Option 3 (for selection of error model refer to the Supplementary Information)

### Load R environment saved at the end of section 7 if this section is executed in a separate script
## if not then there is no need to reload the R environment
load(file = "TMFZ_Renv.RData")

dadaF.3 <- dada(derepF, err = errF_3, multithread = 24, pool = "pseudo")
dadaR.3 <- dada(derepR, err = errR_3, multithread = 24, pool = "pseudo")

### Inspecting the returned dada-class object:only the first sample
dadaF.3[[1]]
dadaR.3[[1]]

saveRDS(dadaF.3, "/path_to_output_directory/TMFZ/03_tabletax/dadaF.3.rds")
saveRDS(dadaR.3, "/path_to_output_directory/TMFZ/03_tabletax/dadaR.3.rds")

### Merge paired-end reads
mergers.3 <- mergePairs(dadaF.3, derepF, dadaR.3, derepR, verbose = TRUE)

### Inspect the merger data.frame from the first sample 
head(mergers.3[[1]])

saveRDS(mergers.3, "/path_to_output_directory/TMFZ/03_tabletax/mergers.3.rds")

### Construct ASV table based on error model from Option 1
seqtab.3 <- makeSequenceTable(mergers.3)

saveRDS(seqtab.3, "/path_to_output_directory/TMFZ/03_tabletax/seqtab.3.rds")

dim.seqtab.3 <- dim(seqtab.3)
dim.seqtab.3

saveRDS(dim.seqtab.3, "/path_to_output_directory/TMFZ/03_tabletax/dim.seqtab.3.rds")

### Inspect distribution of sequence lengths
seqlen.3 <- table(nchar(getSequences(seqtab.3)))
seqlen.3

saveRDS(seqlen.3, "/path_to_output_directory/TMFZ/03_tabletax/seqlen.3.rds")

### Histogram of read distribution
hist.3 <- hist(nchar(getSequences(seqtab.3)), main="Distribution of sequence lengths using Option 3 alter loess function (weights only) and enforce monotonicity")

saveRDS(hist.3, "/path_to_output_directory/TMFZ/03_tabletax/hist.3.rds")

pdf(file="/hpctmp/chrisg25/R_out/TMFZ/hist.3.pdf", width = 15, height = 8.5)
hist.3
dev.off()


### Save the R environment specific for the traditional error model
save.image(file = "TMFZ_Renv.RData")

#################################################################################################

### 9. Remove Chimeras And Assign taxonomy for ASVs inferred with Option 3 error model ####

##################################################################################################

### Although dada2 has searched for indel errors and subsitutions, there may still be chimeric sequences in the dataset.
### Chimeric sequences that are derived from forward and reverse sequences from two different organisms becoming fused together during PCR and/or sequencing).

### To identify chimeras, search for rare sequence variants that can be reconstructed by combining left-hand and right-hand segments from two
#   more abundant “parent” sequences.

### it is not uncommon for a majority of sequence variants to be removed but the abundance of these variant should be minimal.

### Load R environment saved at the end of section 8.1 if this section is executed in a separate script
## if not then there is no need to reload the R environment
load(file = "TMFZ_Renv.RData")

### Read in RDS 
SeqTable.all.3  <- readRDS("/path_to_output_directory/TMFZ/03_tabletax/seqtab.3.rds")


### Remove chimeras
seqtab.nochim.3 <- removeBimeraDenovo(SeqTable.all.3, method="pooled", multithread=4, verbose = TRUE)

saveRDS(seqtab.nochim.3, "/path_to_output_directory/TMFZ/03_tabletax/seqtab.nochim.3.rds")

### Check chimera table & Print percentage of sequences that were not chimeric.

## read abundances of non-chimeric variants
total.reads.3 <- sum(seqtab.all.3)

total.nonchimeric.reads.3 <- sum(seqtab.nochim.3)

saveRDS(total.reads.3, "/path_to_output_directory/TMFZ/03_tabletax/total.reads.3.rds")

saveRDS(total.nonchimeric.reads.3, "/path_to_output_directory/TMFZ/03_tabletax/total.nonchimeric.reads.3.rds")

print( paste("Total reads: ", total.reads.3)) # 49432927 reads

print( paste("Total non-chimeric reads : ", total.nonchimeric.reads.3)) # 45947873 reads

nonchimera.3 <- 100*sum(seqtab.nochim.3)/sum(seqtab.all.3)

saveRDS(nonchimera.3, "/path_to_output_directory/TMFZ/03_tabletax/nonchimera.3.rds")

print( paste("% Abundance of non-chimeric reads: ", nonchimera.3)) # 92.93%

## gives total number of non-chimeric ASVs
dim.seqtab.nochime.3 <- dim(seqtab.nochim.3)

dim.seqtab.nochime.3 #   395 samples and 57589 ASVs

saveRDS(dim.seqtab.nochime.3, "/path_to_output_directory/TMFZ/03_tabletax/dim.seqtab.nochime.3.rds")


### TRACK READS THROUGH THE PIPELINE ####

getN <- function(x) sum(getUniques(x))

### filt_out returns 2 columnns: reads_in and reads_our

### track reads using  Traditional DADA2 error model
track.3 <- cbind(pre_out, filt_out, sapply(dadaF.3, getN), sapply(dadaR.3, getN), sapply(mergers.3, getN), rowSums(seqtab.nochim.3))
colnames(track.3) <- c("pre.DADA2.input", "pre.DADA2.filtered", "DADA2.input", "DADA2.filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track.3) <- sample.names

saveRDS(track.3, "/path_to_output_directory/TMFZ/03_tabletax/track.3.rds")

### Assign Taxonomy using naive Bayesian classifier method ####

### After removing chimeras, use a taxonomy database to train a classifer-algorithm to assign names to sequence variants

### Download DADA2 formated database from https://benjjneb.github.io/dada2/training.html: 
#   NOTE: As of Silva version 138, the official DADA2-formatted reference fastas are optimized for classification of Bacteria and Archaea, 
#   and are not suitable for classifying Eukaryotes.

### The assignTaxonomy function takes as input a set of sequences to be classified, and a training set of reference sequences with known taxonomy.Then 
#   outputs taxonomic assignments with at least minBoot bootstrap confidence (minboot = 50, default minBoot = 50) # minimal bootstrap confidence

### The dada2 package also implements a method to make species level assignments based on exact matching between ASVs and sequenced reference strains. 
#   Recent analysis suggests that exact matching (or 100% identitiy) is the only appropriate way to assign species to 16S gene fragments. 
#   Currently, species-assignment training FASTAS are available for the Silva and RDP 16S databases.

### Assign taxonomy with Silva db v138 using silva_nr99_v138.1_train_set.fa : with taxonomy till genus level (tutorial) and
#   add species assignment with silva_species_assignment_v128.fa.gz

### If reads do not seem to be appropriately assigned, for example lots of your bacterial 16S sequences are being assigned as Eukaryota NA NA NA NA NA,
#   reads may be in the opposite orientation as the reference database. 
#   Get dada2 to try the reverse-complement orientation with assignTaxonomy(..., tryRC=TRUE) 

### When using the outputBootstraps argument, the object returns both the taxa table[1] and the bootstrap values table[2]. 
#   outputBootstraps  returns a list with both the assigned taxonomy matrix, and a corresponding matrix of the bootstrap values
#   The bootstrap confidence is the number of times the subset of kmers (size =8, Naive Bayesian Classifier for Rapid Assignment of rRNA Sequences into the
#   New Bacterial Taxonomy) matched the same taxon (not necessarily the same reference sequence) as the full kmer representation of the query sequence.
#   Exact matching is not required, it is only a best match.
#   It is possible to get spurious assignments with high bootstrap confidence when the reference database does not include appropriate outgroups,
#   i.e. representatives of other taxa that would be the best match to the subsetted kmer profiles if they were present.

### Assign taxonomy at genus and species level for reads
taxa.nb.3 <- assignTaxonomy(seqtab.nochim.3, "/path_to_SILVA_database_directory/SILVA.FOR.DADA2/silva_nr99_v138.1_train_set.fa", tryRC = TRUE,  multithread=TRUE, outputBootstraps = TRUE, verbose = TRUE)
taxa.nb.3.1 <- addSpecies(taxa.nb.3[[1]], "/path_to_SILVA_database_directory/SILVA.FOR.DADA2/silva_species_assignment_v138.1.fa")
head(taxa.nb.3.1)

saveRDS(taxa.nb.3.1, "/path_to_output_directory/TMFZ/03_tabletax/taxa.nb.3.1.rds")
saveRDS(taxa.nb.3[[2]], "/path_to_output_directory/TMFZ/03_tabletax/taxa.nb.3.2.rds")

### Check Taxonomy
taxa.nb.3.print <- taxa.nb.3.1
rownames(taxa.nb.3.print) <- NULL
head(taxa.nb.3.print)

saveRDS(taxa.nb.3.print, "/path_to_output_directory/TMFZ/03_tabletax/taxa.nb.3.print.rds")

save.image(file = "TMFZ_Renv.RData")

######################################################################################################

### 10. Generate plots for tracked reads for all error models for comparison ####

########################################################################################################

### Table of tracked reads through DADA2 pipeline:
# Cutadapt filtered out and DADA2.input are the same.
# DenoisedF and DenoisedR are reads after dereplication and inference of sequence variants from dereplicated reads using trained error models.

track.3 <- readRDS("/path_to_output_directory/TMFZ/03_tabletax/track.3.rds")

### Get Column totals to get raw reads
track.total.3 <- data.frame(t(colSums(track.3)))


### Percentage of reads remaining at the end of each step of the pipeline (output of each step compared to raw reads -initial)
track.3.p1 <- as.data.frame(track.3) %>% data.frame() %>%
  dplyr::mutate (pre.DADA2.NFilt.pct = pre.DADA2.filtered/pre.DADA2.input*100,
                 Cutadapt.Filt.pct = DADA2.input/pre.DADA2.input*100,
                 DADA2.TrimNFilt.pct = DADA2.filtered/pre.DADA2.input*100,
                 Denoised.F.pct = denoisedF/pre.DADA2.input*100,
                 Denoised.R.pct = denoisedR/pre.DADA2.input*100,
                 merged.pct = merged/pre.DADA2.input*100,
                 nonchim.pct = nonchim/pre.DADA2.input*100)  %>%
                 dplyr::select(ends_with(".pct")
  )

 
### Percentage of reads lost at each step of the pipeline (output of each step subtracted from the input and compared to input)
### To determine at which step/s most reads are lost

track.3.p2 <- as.data.frame(track.3) %>% data.frame() %>%
  dplyr::mutate (pre.DADA2.NFilt.pct = (pre.DADA2.input-pre.DADA2.filtered)/pre.DADA2.input*100,
                 Cutadapt.Filt.pct = (pre.DADA2.filtered-DADA2.input)/pre.DADA2.filtered*100,
                 DADA2.TrimNFilt.pct = (DADA2.input-DADA2.filtered)/DADA2.input*100,
                 Denoised.F.pct = (DADA2.filtered-denoisedF)/DADA2.filtered*100,
                 Denoised.R.pct = (DADA2.filtered-denoisedR)/DADA2.filtered*100,
                 merged.pct = (((denoisedF-merged)/denoisedF*100)+((denoisedR-merged)/denoisedR*100))/2,
                 nonchim.pct = (merged-nonchim)/merged*100)  %>%
                 dplyr::select(ends_with(".pct")
  )


##Plot tracked reads

# summary stats of tracked reads averaged across samples to include in plots

track.3.p1.avg <- track.3.p1  %>% summarize_at(vars(ends_with(".pct")), list(avg = mean))
track.3.p2.avg <- track.3.p2  %>% summarize_at(vars(ends_with(".pct")), list(avg = mean))

### as boxplots

track.3.plot  <- as.data.frame(track.3) %>% data.frame() %>%
  mutate(Sample = rownames(.)) %>% gather(key = "Step", value = "Reads", -Sample) %>%
  mutate(Step = factor(Step, levels = c("pre.DADA2.input", "pre.DADA2.filtered", "DADA2.input", "DADA2.filtered", "denoisedF", "denoisedR", "merged", "nonchim"))) %>%
    ggplot(aes(x = Step, y = Reads)) +
    geom_boxplot (varwidth =TRUE ,fill = "cadetblue",alpha =0.3,outlier.colour="red") +
    stat_summary(fun.y = mean, geom = "point", group = 1, color = "blue", size = 3, alpha = 0.7,shape = 18) +
      geom_label(data = t(track.3.p1.avg[1:6]) %>% data.frame() %>% rename(Percent = 1) %>%
      mutate(Step = c("pre.DADA2.filtered", "DADA2.input", "DADA2.filtered", "denoisedF", "denoisedR", "merged"),
      Percent = paste(round(Percent, 2), "%")), aes(label = Percent), y = 1.1 * max(track.3[,3])) +
          geom_label(data = track.3.p1.avg[7] %>% data.frame() %>% rename(total = 1),
          aes(label = paste("Total\nRemaining:\n", round(track.3.p1.avg[1,7], 2), "%")),
          y = mean(track.3[,8]), x = 8.7) +
  geom_label(aes(label = paste("Total reads (%)")), y = 1.1 * max(track.3[,3]), x = 8.7) +
  geom_label(aes(label = paste("Loss in reads (%)")), y = 1.05 * max(track.3[,3]), x = 8.7) +
      geom_label(data = t(track.3.p2.avg[1:6]) %>% data.frame() %>% rename(Percent = 1) %>%
      mutate(Step = c("pre.DADA2.filtered", "DADA2.input", "DADA2.filtered", "denoisedF", "denoisedR", "merged"),
      Percent = paste(round(Percent, 2), "%")), aes(label = Percent), y = 1.05 * max(track.3[,3])) +
          geom_label(data = track.3.p2.avg[7] %>% data.frame() %>% rename(total = 1), 
          aes(label = paste("Lost\nReads:\n", round(track.3.p2.avg[1,7], 2), "%")),
          y = (mean(track.3[,8])), x = 9.3) +
  expand_limits(y = 1.1 * max(track.3[,3]), x = 9.7) +
  theme_classic()

track.3.plot

ggsave("track.3.plot.pdf",height = 8, width = 14)


#########################################################################################################

### 11. Create phyloseq object ####

#########################################################################################################

### Phyloseqs require ASV abundance table, Taxonomy table and Metadata 

### Read non-chimeric tables (abundance table)
seqtab.nochim.3 <- readRDS("/path_to_output_directory/TMFZ/03_tabletax/seqtab.nochim.3.rds")

## Taxonomy table using naive Bayesian classifier method and SILVA
### assignTaxonomy implements the RDP Naive Bayesian Classifier algorithm described in Wang et al. Applied and Environmental Microbiology 2007, with kmer size 8 and 100 bootstrap replicates.
taxa.nb.3.1 <- readRDS("/path_to_output_directory/TMFZ/03_tabletax/taxa.nb.3.1.rds")

### Water chemistry variables in Table 1 are provided to the 'sample_metadata' slot in the phyloseq
### Water chemistry data is loaded in R in excel format
### Add column with Sample_id so that each sample_id is a row with the respective water chemistry variables in separate columns
### Latitude and Latitude are to be separated out to 2 separate columns
### Add column '`Temp.adj (°C)`' with adjusted temperatures for each sample as reflected in the sample_id
metadata <- "/path_to_water_chemistry_metadata/water_chemistry_metadata.xlsx"

### Create Phyloseq
### Phyloseq with ASVs inferred from the Option 3 error model
ps.3 <- phyloseq(otu_table(seqtab.nochim.3, taxa_are_rows=FALSE), 
                 sample_data(metadata), 
                 tax_table(taxa.nb.3.1))

# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 57589 taxa and 395 samples ]
# sample_data() Sample Data:       [ 395 samples by 27 sample variables ]
# tax_table()   Taxonomy Table:    [ 57589 taxa by 7 taxonomic ranks ]
# refseq()      DNAStringSet:      [ 57589 reference sequences ]

### Use short names for ASVs instead of using entire DNA sequence. http://benjjneb.github.io/dada2/tutorial.html
### Store the DNA sequences of ASVs in the refseq slot of the phyloseq object, and then rename taxa to a short string. 
### That way, the short new taxa names will appear in tables and plots, and the DNA sequences corresponding to each ASV can be recovered as needed with 'refseq(ps)'.

dna.3 <- Biostrings::DNAStringSet(taxa_names(ps.3))
names(dna.3) <- taxa_names(ps.3)
ps.3 <- merge_phyloseq(ps.3, dna.3)
taxa_names(ps.3) <- paste0("ASV", seq(ntaxa(ps.3)))
ps.3

#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 57589 taxa and 395 samples ]
#sample_data() Sample Data:       [ 395 samples by 27 sample variables ]
#tax_table()   Taxonomy Table:    [ 57589 taxa by 7 taxonomic ranks ]
#refseq()      DNAStringSet:      [ 57589 reference sequences ]

saveRDS(ps.3, "/path_to_output_directory/TMFZ/03_tabletax/ps.3.rds")

#########################################################################################################

### 12. Pre-processing of phyloseq for all downstream analyses ####

#########################################################################################################

### Load phyloseq object
ps.3  <- readRDS("/path_to_output_directory/TMFZ/03_tabletax/ps.3.rds")

### 12.1. Check if there are ASVs with no counts and how many there are: ####
### If TRUE then these ASVs must be removed.
### No ASVs have a total of 0 in any of the datasets.
any(taxa_sums(ps.3) == 0)

### 12.2. Assign NAs in taxonomic ranks higher than domain using 'fantaxtic' with tax and rank #### 
### name_na_taxa assigns the name of the lowest known taxonomic rank for every NA value in each ASV. 
ps.edit.3 <- name_na_taxa(ps.3, na_label = "Unassigned <tax> (<rank>)")

### NA in the Kingdom column will lead to NAs across all 7 taxonomic ranks even with name_na_taxa
### Blanks that are present in tax_table after 'name_na_taxa' across all 7 taxonomic levels were assigned as 'Unassigned'. 
### This new new taxa table is then replaced in the  phyloseq
ASV.3.taxa.new <- (as.matrix(ps.edit.3@tax_table))
ASV.3.taxa.new[is.na(ASV.3.taxa.new)] = "Unassigned"

### Add new tax_table to phyloseq. 
### The new tax_table will not have any blanks
tax_table(ps.edit.3) <- ASV.3.taxa.new

### 12.3. Determine % of of 'Unassigned', 'Mitochondrial' and 'Chlororplast' ASV variants across all 7 taxonomic ranks ####
### Must be calculate before subsetting away and not after

### Subset 'Unassigned ASVs'
ASV.3.taxa.new.unassigned <- (as.data.frame(ASV.3.taxa.new)) %>%
  filter(if_any(everything(), ~ .x == "Unassigned"))

### % of 'Unassigned' ASV variants
Unassigned.percent.3 <- ((nrow(ASV.3.taxa.new.unassigned)/nrow(ASV.3.taxa.new))*100)

### Subset 'Mitochondria' & 'Chloroplast' ASVs
### Classified using : Order = Chlororplast, Family = Mitochondira
Unwanted.taxa.3 <- as.data.frame(ASV.3.taxa.new) %>% filter_all(any_vars(. %in% c('Chloroplast','Mitochondria')))

### % of 'Mitochondria' & 'Chloroplast' ASV variants
Mito.chloro.percent.3 <- ((nrow(Unwanted.taxa.3)/nrow(ASV.3.taxa.new))*100)

### Determine % of total read abundances of 'Unassigned', 'Chloroplast' and 'Mitochindiral' ASVs

### Melt phyloseq
ps.edit.melt.3 <- psmelt(ps.edit.3)

### Total read abundances of 'Unassigned', 'Chloroplast' and 'Mitochindiral' ASVs
ps.edit.unassigned.3 <- sum(ps.edit.melt.3[which(ps.edit.melt.3$Kingdom=='Unassigned'), 3])
ps.edit.chloroplast.3 <- sum(ps.edit.melt.3[which(ps.edit.melt.3$Order=='Chloroplast'), 3])
ps.edit.mitochondria.3 <- sum(ps.edit.melt.3[which(ps.edit.melt.3$Family=='Mitochondria'), 3])

### % of total read abundances of 'Unassigned', 'Chloroplast' and 'Mitochindiral' ASVs
ps.edit.unassigned.reads.percent.3 <- (ps.edit.unassigned.3/sum(ps.edit.melt.3$Abundance))*100 # 0.0049%
ps.edit.chloroplast.reads.percent.3 <- (ps.edit.chloroplast.3/sum(ps.edit.melt.3$Abundance))*100 # 0.23%
ps.edit.mitochondria.reads.percent.3 <- (ps.edit.mitochondria.3/sum(ps.edit.melt.3$Abundance))*100 # 0.67% 


### 12.4. Get 'clean' phyloseq = Phyloseq wihtout 'Unassighed', 'Mitochondiral' & 'Chlororplast' ASVs ####

### Subset phyloseq to exclude 'Unassigned'  across all 7 taxonomic levels, 
ps.no.unassign.3 <- subset_taxa(ps.edit.3, Kingdom !="Unassigned")

### Subset phyloseq to exclude mitochondria and chloroplast ASVs 
ps.clean.3 <- ps.no.unassign.3 %>%
  subset_taxa(Order != "Chloroplast" & Family  != "Mitochondria" | is.na(Family))

# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 56879 taxa and 395 samples ]
# sample_data() Sample Data:       [ 395 samples by 27 sample variables ]
# tax_table()   Taxonomy Table:    [ 56879 taxa by 7 taxonomic ranks ]
# refseq()      DNAStringSet:      [ 56879 reference sequences ]

### Sanity check : compare total number of reads before and after removal of unwanted taxa ('Unassighed', 'Mitochondiral' & 'Chlororplast' ASVs)
sum(sample_sums(ps.edit.3)) # same as first phyloseq - ps.3 : 45947873 reads 
sum(sample_sums(ps.clean.3)) # 45531317 reads

### Sanity check for ps.clean phyloseqs to ensure Unassigned, Mitochondria and Chloroplast ASVs are not present.
### resulting table should have 0 observations
check.nomito.nochloro.no.unassign.3 <- as.data.frame(ps.clean.3@tax_table) %>% filter_all(any_vars(. %in% c('chloroplast','Mitochondria','Unassigned')))

### Save clean phyloseqs
saveRDS(ps.clean.3, "/path_to_output_directory/TMFZ/03_tabletax/ps.clean.3.rds")


### 12.5. Add data for Site labelled AF to AP (collected from the same hotspring) and get new phyloseq ####

### Change name of Site "AF" to "AP" in otu table
ASV.wide.3 <- as.data.frame(ps.clean.3@otu_table) %>% rownames_to_column(var = "sample_id") %>% as_tibble()
ASV.wide.3$sample_id <- str_replace_all(ASV.wide.3$sample_id, "AF_","AP_")
ASV.wide.3 <- ASV.wide.3 %>% column_to_rownames(var = "sample_id")
ASV.wide.3 <- as.matrix(ASV.wide.3)

### Change name of Site "AF" to "AP" in metadata and add "Region" column
metadata <- as.data.frame(metadata) %>% rownames_to_column(var = "sample_id")  %>% as_tibble()

metadata <- as.data.frame(sapply(metadata,gsub,pattern="AF", replacement="AP")) %>% 
  select(sample_id, Country, `Location name`, Location.Code, Site, Latitude, Longitude, everything())

metadata <- metadata %>% column_to_rownames(var = "sample_id")

metadata$`No of samples`[metadata$Site == "AP_46"] <- 20

### Add 'Region' column to metadata
North.Thailand <- c('PT','PP','LN','HS','MK')
Central.Thailand <- 'PB'
South.Thailand <- c('RB','RN')
North.Malaysia <- c('AP','BA','US')
South.Malaysia <- c('KJ','SE','LA')
Singapore <- 'SW'

metadata <- metadata %>% 
  mutate(Region = case_when(Location.Code %in% North.Thailand ~ 'North.Thailand', Location.Code == 'PB' ~ 'Central.Thailand', Location.Code %in% South.Thailand ~ 'South.Thailand', Location.Code %in% North.Malaysia ~ 'North.Malaysia', Location.Code %in% South.Malaysia ~ 'South.Malaysia', Location.Code %in% 'SW' ~ 'Singapore'  ))

metadata <- metadata %>% select(Sample_id, Country, Region, Location.Code, Site, Location, everything())

### ensure numeric columns are numeric in the metadata. 
### 11 is used as an example. Specify any range of column numbers
metadata[11:ncol(metadata)] <- lapply(metadata[11:ncol(metadata)], as.numeric)

## non-numeric columns as factors
metadata$`Temp.adj (°C)` <- as.factor(metadata$`Temp.adj (°C)`)
metadata$Location.Code <- as.factor(metadata$Location.Code)

### Create new and clean phyloseqs 

ASV.wide.3.new = otu_table(ASV.wide.3, taxa_are_rows = FALSE)

samples = sample_data(metadata)

## Use same tax_table as clean phyloseq as no changes in tax table
ASV.tax.wide.3 = tax_table(ps.clean.3)

## new and clean phyloseq
new.ps.clean.3 <- phyloseq(ASV.wide.3.new, ASV.tax.wide.3, samples)

### Save clean and new phyloseq
saveRDS(new.ps.clean.3, "/path_to_output_directory/TMFZ/03_tabletax/new.ps.clean.3.rds")


### 12.6. Rarefaction to get rarefied dataset (To use for alpha diversity analyses) ####

### Sequencing depth in this case (amplicon sequencing) is the same as the total reads in each sample

### Different options to rarefy
min.seq.depth.3 <- min(sample_sums(new.ps.clean.3))
max.seq.depth.3 <- max(sample_sums(new.ps.clean.3))
mean.seq.depth.3 <- mean(sample_sums(new.ps.clean.3))
median.seq.depth.3 <- median(sample_sums(new.ps.clean.3))

### min.seq.depth.3 =  62868
### max.seq.depth.3 = 179080
### mean.seq.depth.3 = 115269.2
### median.seq.depth.3 = 112225

### Rarefying with median sequencing depth causes loss in samples as well as ASVs after random sampling
### Therefore, rarefy with minimum sequencing depth

rarefied.min.3 <- rarefy_even_depth(new.ps.clean.3,rngseed=123,sample.size=min(sample_sums(new.ps.clean.3)), replace=F)

rarefied.min.3

## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 50538 taxa and 395 samples ]
## sample_data() Sample Data:       [ 395 samples by 27 sample variables ]
## tax_table()   Taxonomy Table:    [ 50538 taxa by 7 taxonomic ranks ]
## refseq()      DNAStringSet:      [ 50538 reference sequences ]

summarize_phyloseq(rarefied.min.3)

sum(sample_sums(rarefied.min.3))
## Total number of reads = 24832860

### Save rarefied phyloseq
saveRDS(rarefied.min.3, "/path_to_output_directory/TMFZ/03_tabletax/rarefied.min.3.rds")


### 12.7. Get Final Working Dataset (Rarefied and Filtered) = ASVs with relative abundance > 1% from rarefied dataset. ####
### The final working dataset is used for all downstream analyses except for alpha diversity estimation.

### First Get phyloseq with abundance raw counts (integers) transformed to percentages
rarefied.min.prop.3 <- transform_sample_counts(rarefied.min.3, function(x) x / sum(x)*100)

total.3 = sample_sums(rarefied.min.prop.3)

## Rarefied and filtered (ASVs with relative abundance > 1%) phyloseq with abundance as percentages
rarefied.min.prop.exclude.1.3 = filter_taxa(rarefied.min.prop.3, function(x) sum(x > total.3*0.01) > 0,TRUE)

rarefied.min.prop.exclude.1.3

## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 572 taxa and 395 samples ]
## sample_data() Sample Data:       [ 395 samples by 27 sample variables ]
## tax_table()   Taxonomy Table:    [ 572 taxa by 7 taxonomic ranks ]
## refseq()      DNAStringSet:      [ 572 reference sequen


##  Rarefied and filtered (ASVs with relative abundance > 1%) phyloseq with abundance as integers
## Subset taxa with rel.ab > 1% from phyloseq object with counts as integers
subtaxa.rarefied.min.exclude.1.3 = taxa_names(rarefied.min.prop.exclude.1.3)
rarefied.min.int.exclude.1.3 <- prune_taxa(subtaxa.rarefied.min.exclude.1.3, rarefied.min.3)

rarefied.min.int.exclude.1.3

## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 572 taxa and 395 samples ]
## sample_data() Sample Data:       [ 395 samples by 27 sample variables ]
## tax_table()   Taxonomy Table:    [ 572 taxa by 7 taxonomic ranks ]
## refseq()      DNAStringSet:      [ 572 reference sequences ]

summarize_phyloseq(rarefied.min.int.exclude.1.3)
## Total number of reads = 20833623"

### Save rarefied and filtered phyloseqs (final working dataset)
saveRDS(rarefied.min.prop.3, "/path_to_output_directory/TMFZ/03_tabletax/rarefied.min.prop.3.rds")
saveRDS(rarefied.min.prop.exclude.1.3, "/path_to_output_directory/TMFZ/03_tabletax/rarefied.min.prop.exclude.1.3.rds")
saveRDS(rarefied.min.int.exclude.1.3, "/path_to_output_directory/TMFZ/03_tabletax/rarefied.min.int.exclude.1.3.rds")

### 12.8. Sanity check between ####

##  rarefied dataset as integers
rarefied.min.wide.3 <-  as.matrix(as.data.frame(rarefied.min.3@otu_table))

## rarefied dataset as percentages
rarefied.min.prop.wide.3 <- as.matrix(as.data.frame(rarefied.min.prop.3@otu_table))

## wide tables as percentages and excluded relative abundances of ASVs <= 1%
rarefied.min.prop.wide.exclude.1.3 <- as.data.frame(rarefied.min.prop.exclude.1.3@otu_table)

## wide tables as intgers and excluded relative abundances of ASVs <= 1%
rarefied.min.wide.exclude.1.3  <- as.data.frame(rarefied.min.int.exclude.1.3@otu_table)

##sanity check for integers and percentages
setequal(colnames(rarefied.min.wide.exclude.1.3),colnames(rarefied.min.wide.prop.exclude.1.3))
setequal(colnames(rarefied.min.wide.3),colnames(rarefied.min.prop.wide.3 ))
## TRUE


#############################################################################################

### 13. Detect % of reads from taxa considered as human contaminants #####

###############################################################################################

### Taxa considered as human contaminants were identified with the following genera
contaminants <- c("Bacteroides","Bifidobacterium","Corynebacterium","Cutibacterium","Escherichia","Faecalibacterium", "Haemophilus", "Klebsiella", "Lactobacillus", "Listeria", "Moraxella", "Neisseria", "Porphyromonas", "Prevotella", "Propionibacterium", "Salmonella", "Shigella", "Staphylococcus", "Streptococcus", "Veillonella")

## Tax_glom bacteria integers
rare.genus.exclude.1.3 <- tax_glom(rarefied.min.int.exclude.1.3, taxrank="Genus")

## melt 
rare.wide.exclude.1.melt.3 <- psmelt(rare.genus.exclude.1.3)

## 0.145%
contaminants.rarefied <- sum(rare.wide.exclude.1.melt.3[which(rare.wide.exclude.1.melt.3$Genus %in% contaminants), 3])/sum(rare.wide.exclude.1.melt.3$Abundance)*100

### sanity check wihtout tax_glom at genus level

## melt 
rare.ps.clean.wide.melt.exclude.1.3 <- psmelt(rarefied.min.int.exclude.1.3)

## 0.145%
contaminants.rarefied <- sum(rare.ps.clean.wide.melt.exclude.1.3[which(rare.ps.clean.wide.melt.exclude.1.3$Genus %in% contaminants), 3])/sum(rare.ps.clean.wide.melt.exclude.1.3$Abundance)*100

###################################################################################

### End of Section ###

