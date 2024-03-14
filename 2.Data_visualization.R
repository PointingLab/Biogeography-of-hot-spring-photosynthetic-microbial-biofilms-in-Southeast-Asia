###############################################################################################################################################

### NOTES: 

### R codes to process phyloseqs created with 1.RawSeqs_to_Phyloseq_short for data visualization
### Codes for all the figures (main text and supplementary) in the paper are presented sequentially.
### Codes for statistical tests are also included in the respective sections.

#################################################################################################################################################

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
library(plotly) # enables creation of interactive graphs, especially helpful for quality plots
library("phyloseq"); packageVersion("phyloseq")
library(colorRamp2)
library("Biostrings")
library(fantaxtic) # to assign NAs in taxa table as Unassigned and proper labels
library(phangorn)
library(microbiome)
library(ggpubr)
library(geosphere)
library(geodist)
library(ggrepel)
library(gtools)
library(ggforce)
library(pulsar)
library(SpiecEasi)
library(NetCoMi)
library(microbiomeutilities)
library("metacoder")
library(microeco)
library(magrittr)
library(ggcor)
library(microbiomeutilities)
library(DECIPHER)
library(phangorn)
library(microViz)

### Read all saved phyloseqs ####

### Rarefied
rarefied.min.3 <- readRDS("/path_to_output_directory/TMFZ/03_tabletax/rarefied.min.3.rds")
rarefied.min.prop.3 <- readRDS("/path_to_output_directory/TMFZ/03_tabletax/rarefied.min.prop.3.rds")

### Final working dataset (Rarefied and filtered)
rarefied.min.int.exclude.1.3 <- readRDS("/path_to_output_directory/TMFZ/03_tabletax/rarefied.min.int.exclude.1.3.rds")
rarefied.min.prop.exclude.1.3 <- readRDS("/path_to_output_directory/TMFZ/03_tabletax/rarefied.min.prop.exclude.1.3.rds")

### Extract dataframes from phyloseqs ####

## rarefied dataset as integers
rarefied.min.wide.3 <-  as.matrix(as.data.frame(rarefied.min.3@otu_table))

## rarefied dataset as percentages
rarefied.min.prop.wide.3 <- as.matrix(as.data.frame(rarefied.min.prop.3@otu_table))

## wide tables as intgers and relative abundances of ASVs > 1%
rarefied.min.wide.exclude.1.3  <- as.data.frame(rarefied.min.int.exclude.1.3@otu_table)

## wide tables as percentages and relative abundances of ASVs > 1%
rarefied.min.prop.wide.exclude.1.3 <- as.data.frame(rarefied.min.prop.exclude.1.3@otu_table)

## Get tax table from any of the final working daatset phyloseqs as both are the same
Tax <- as.matrix(as.data.frame(rarefied.min.int.exclude.1.3@tax_table)) 

## Get metadata. use any phyloseq as all have the same metadata
metadata <- as.matrix(as.data.frame(rarefied.min.int.exclude.1.3@sam_data))
  

####################################################################################################################

### Fig. S1. Sampling curves for 16S rRNA gene sequencing, shown for unrarefied data and for the sequencing data rarefied to the minimum sequencing depth (62,868 reads per sample ####

######################################################################################################################

## get otu table as matrix for rarefaction curves

## Function rarecurve draws a rarefaction curve for each row of the input data. The rarefaction curves are evaluated using the interval of step sample sizes, always including 1 and total sample size. 

#If sample is specified, a vertical line is drawn at sample with horizontal lines for the rarefied species richnesses.

### axis description: sample size - number od raw reads and Species - ASVs

# Number of INDIVIDULS per site (?)
raremax <- min(rowSums(rarefied.min.wide.3)) # = 62868; 
raremax

raremax.all <- min(rowSums(ASV.wide.3)) # = 62868; 
raremax.all

#total number of species at each site (row of data)
total.species.all <- specnumber(ASV.wide.3)
total.species.rarefy <- specnumber(rarefied.min.wide.3)

# rarefy, w/ raremax as input (?)
rarefied <- rarefy(ASV.wide.3, raremax)


plot(total.species.all,rarefied, xlab = "Observed No. of Species", 
     ylab = "Rarefied No. of Species",
     main = " plot(rarefy(ASV.wide.3, raremax))")


## without subsanpling
rarecurve.rarefied.nosub <- rarecurve(rarefied.min.wide.3, step = 250,
                                      col = "khaki3",
                                      cex = 0.6,label=F,
                                      main = "rarecurve()")

#saveRDS(rarecurve.rarefied.nosub, "rarecurve.rarefied.nosub.rds")

rarecurve.rarefied.nosub <- readRDS("rarecurve.rarefied.nosub.rds")


rarecurve.all.nosub <- rarecurve(ASV.wide.3, step = 250,
                                 col = "khaki3",
                                 cex = 0.6,label=F,
                                 main = "rarecurve.all.nosub")

#saveRDS(rarecurve.all.nosub, "rarecurve.all.nosub.rds")


## With subsampling
rarecurve.rarefied.sub <- rarecurve(rarefied.min.wide.3, step = 250, sample = raremax,
                                    col = "khaki3",
                                    cex = 0.6,label=F,
                                    main = "rarecurve()")

saveRDS(rarecurve.rarefied.sub, "rarecurve.rarefied.sub.rds")


rarecurve.all.sub <- rarecurve(ASV.wide.3, step = 250, sample = raremax.all,
                               col = "khaki3",
                               cex = 0.6,label=F,
                               main = "rarecurve.all.sub")

saveRDS(rarecurve.all.sub, "rarecurve.all.sub.rds")


###specaccum(): “Function specaccum finds species accumulation curves or the number of species for a certain number of sampled sites or individuals.”"

### Speiecesaccum Vs rarefaction

##Species accumulation “Species accumulation models and species poolmodels study collections of sites, and their species richness, or try to estimate the number of unseen species.” (Oksanen 2016)

###Rarefaction “Species richness increases with sample size, and differences in richness actually may be caused by differences in sample size. To solve this problem, we may try to rarefy species richness to the same number of individuals… Rarefaction is sometimes presented as an ecologically meaningful alternative to dubious diversity indices (Hurlbert, 1971), but the differences really seem to be small.” (Oksanen 2016)

sp2 <- specaccum(rarefied.min.wide.3[,1:100])
plot(sp2, main = "method = collector")

# method = "collector"
plot(specaccum(ASV.wide.3, method = "collector"),main = "method = collector")

# method = "rarefaction"
plot(specaccum(ASV.wide.3, method = "rarefaction"),main = "method = rarefaction")
legend("bottomleft", legend = "Note x-axis is # of plots")



##############################################################################################################

### Fig 1B : Distance decay plot illustrating the strong correlation between geographic and phylogenetic distance for hot spring biofilm communities (N = 395) ####
### Fig S3 : Distance decay plots of phylogenetic versus geographic distance for ecological groups of bacteria (N = 395) ####

##############################################################################################################

### Calculating UniFrac distances for distance decay requires phylogenetic tree ####

### Create phylogeentic tree and add to phyloseq

## first get reference sequences from phyloseq.
refseqs.rare.3 = refseq(rarefied.min.int.exclude.1.3)
seqs.rare.3 <- getSequences(refseqs.rare.3)

# This propagates to the tip labels of the tree
names(seqs.rare.3) <- taxa_names(refseqs.rare.3) 

# sequence alignment with 'DECIPHER'
alignment.rare.3 <- AlignSeqs(DNAStringSet(seqs.rare.3), anchor=NA, processors=NULL)

# converting into phyDat format
phang.align.rare.3 <- phyDat(as(alignment.rare.3, "matrix"), type="DNA")

# compute pairwaise distances
dm.rare.3 <- dist.ml(phang.align.rare.3)

# construct neighbor-joining tree
treeNJ.rare.3 <- NJ(dm.rare.3)

# computes likelihood of phylogenetic tree
fit.rare.3 <- pml(treeNJ.rare.3, data = phang.align.rare.3)

# re-fit a model
fitGTR.rare.3 <- update(fit.rare.3, k=4, inv=0.2)
fitGTR.rare.3 <- optim.pml(fitGTR.rare.3, model="GTR", optInv=TRUE, optGamma=TRUE, rearrangement = "stochastic", control = pml.control(trace = 0))

tree.rare.3 <- fitGTR.rare.3$tree

phytree.rare = phy_tree(tree.rare.3)

## This is an unrooted tree
is.rooted(phytree.rare)
## FALSE 

### Create phyloseqs with unrooted phylogentic trees (overall community and ecological groups) ####

com.rarefied.min.int.exclude.1.3 <- merge_phyloseq(rarefied.min.int.exclude.1.3, phytree.rare)

com.rarefied.min.prop.exclude.1.3 <- merge_phyloseq(rarefied.min.prop.exclude.1.3, phytree.rare)

## Subset phyloseqs for each of the ecological groups with proportions/percentages as abundance in otu_table ####

## Chloroflexi
rarefied.chloroflexi.prop.exclude.1.3 <- com.rarefied.min.prop.exclude.1.3 %>%
  subset_taxa(Phylum == "Chloroflexi" & Class  == "Chloroflexia")

# otu_table()   OTU Table:         [ 55 taxa and 395 samples ]
# sample_data() Sample Data:       [ 395 samples by 33 sample variables ]
# tax_table()   Taxonomy Table:    [ 55 taxa by 9 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 55 tips and 54 internal nodes ]
# refseq()      DNAStringSet:      [ 55 reference sequences ]


## Cyanobacteria
rarefied.cyanobacteria.prop.exclude.1.3 <- com.rarefied.min.prop.exclude.1.3 %>%
  subset_taxa(Phylum == "Cyanobacteria" & Class  == "Cyanobacteriia")

# otu_table()   OTU Table:         [ 74 taxa and 395 samples ]
# sample_data() Sample Data:       [ 395 samples by 33 sample variables ]
# tax_table()   Taxonomy Table:    [ 74 taxa by 9 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 74 tips and 73 internal nodes ]
# refseq()      DNAStringSet:      [ 74 reference sequences ]

## First list genera for photosynthetic Proteoacteria and Order for other photosynthetic groups
proteobacteria.genus <- c("DSSF69","Elioraea","Methylobacterium-Methylorubrum","Rhodomicrobium","Roseomonas","Sandaracinobacter","Tabrizicola","Unassigned Rhodobacteraceae (Family)","Unassigned Sphingomonadaceae (Family)","AAP99","Allochromatium","Caldimonas","Curvibacter","DSSD61","Thiolamprovum","Unassigned B1-7BS (Family)","Unassigned Burkholderiales (Order)","Unassigned Comamonadaceae (Family)","Unassigned Gammaproteobacteria (Class)","Unassigned Rhodocyclaceae (Family)","Unassigned Sutterellaceae (Family)","Z-35")

other.photo <- c("Chlorobiales","Chloracidobacteriales")

rarefied.proteobacteria.otherphoto.prop.exclude.1.3 <- com.rarefied.min.prop.exclude.1.3 %>%
  subset_taxa(Order %in% other.photo | Genus %in% proteobacteria.genus)

#otu_table()   OTU Table:         [ 58 taxa and 395 samples ]
#sample_data() Sample Data:       [ 395 samples by 33 sample variables ]
#tax_table()   Taxonomy Table:    [ 58 taxa by 9 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 58 tips and 57 internal nodes ]
#refseq()      DNAStringSet:      [ 58 reference sequences ]

## all photosynthetic
photosynthetic.class <- c("Cyanobacteriia","Chloroflexia")

rarefied.photosynthetic.prop.exclude.1.3 <- com.rarefied.min.prop.exclude.1.3 %>%
  subset_taxa( Class %in% photosynthetic.class | Order %in% other.photo | Genus %in% proteobacteria.genus )

#otu_table()   OTU Table:         [ 187 taxa and 395 samples ]
#sample_data() Sample Data:       [ 395 samples by 33 sample variables ]
#tax_table()   Taxonomy Table:    [ 187 taxa by 9 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 187 tips and 186 internal nodes ]
#refseq()      DNAStringSet:      [ 187 reference sequences ]

## Chemolitotrophs

#Phylum == 'Aquificota' & Class == 'Aquificae' ~ 'Chemolithoautotroph',
#Phylum == 'Calditrichota' & Class == 'Calditrichia' ~ 'Chemolithoautotroph',
#Phylum == 'Desulfobacterota' ~ 'Chemolithoautotroph',
#Phylum == 'Elusimicrobiota' ~ 'Chemolithoautotroph',
#Phylum == 'Nitrospirota' ~ 'Chemolithoautotroph',
#Phylum == 'Patescibacteria' ~ 'Chemolithoautotroph',
#Phylum == 'Sva0485' ~ 'Chemolithoautotroph',
#Class == 'Leptospirae' ~ 'Chemolithoautotroph'

chemolithotroph.phylum <- c("Aquificota","Calditrichota","Desulfobacterota","Elusimicrobiota","Nitrospirota","Patescibacteria","Sva0485")

chemolithotroph.class <- c("Aquificae","Calditrichia","Leptospirae")

rarefied.chemolithotrophs.prop.exclude.1.3 <- com.rarefied.min.prop.exclude.1.3 %>%
  subset_taxa(Phylum %in% chemolithotroph.phylum  |  Class %in% chemolithotroph.class)

#otu_table()   OTU Table:         [ 41 taxa and 395 samples ]
#sample_data() Sample Data:       [ 395 samples by 33 sample variables ]
#tax_table()   Taxonomy Table:    [ 41 taxa by 9 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 41 tips and 40 internal nodes ]
#refseq()      DNAStringSet:      [ 41 reference sequences ]

#### all  heterotrophs only
rarefied.only.heterotrophs.prop.exclude.1.3 <- com.rarefied.min.prop.exclude.1.3 %>%
  subset_taxa(!Phylum %in% chemolithotroph.phylum & !Class %in% chemolithotroph.class & !Class %in% photosynthetic.class & !Order %in% other.photo & !Genus %in% proteobacteria.genus )

#otu_table()   OTU Table:         [ 344 taxa and 395 samples ]
#sample_data() Sample Data:       [ 395 samples by 33 sample variables ]
#tax_table()   Taxonomy Table:    [ 344 taxa by 9 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 344 tips and 342 internal nodes ]
#refseq()      DNAStringSet:      [ 344 reference sequences ]

### save all phyloseqs

saveRDS(rarefied.chloroflexi.prop.exclude.1.3, "/path_to_output_directory/TMFZ/03_tabletax/rarefied.chloroflexi.prop.exclude.1.3.rds")
saveRDS(rarefied.cyanobacteria.prop.exclude.1.3, "/path_to_output_directory/TMFZ/03_tabletax/rarefied.cyanobacteria.prop.exclude.1.3.rds")
saveRDS(rarefied.proteobacteria.otherphoto.prop.exclude.1.3, "/path_to_output_directory/TMFZ/03_tabletax/rarefied.proteobacteria.otherphoto.prop.exclude.1.3.rds")
saveRDS(rarefied.photosynthetic.prop.exclude.1.3, "/path_to_output_directory/TMFZ/03_tabletax/rarefied.photosynthetic.prop.exclude.1.3.rds")
saveRDS(rarefied.chemolithotrophs.prop.exclude.1.3, "/path_to_output_directory/TMFZ/03_tabletax/rarefied.chemolithotrophs.prop.exclude.1.3.rds")
saveRDS(rarefied.only.heterotrophs.prop.exclude.1.3, "/path_to_output_directory/TMFZ/03_tabletax/rarefied.only.heterotrophs.prop.exclude.1.3.rds")

### Create uniFrac distance matrices for all phyloseqs #####

rarefied.unifrac.dist <- as.matrix(phyloseq::distance(com.rarefied.min.prop.exclude.1.3, method = "unifrac"))

rarefied.unifrac.chloroflexi.dist <- as.matrix(phyloseq::distance(rarefied.chloroflexi.prop.exclude.1.3, method = "unifrac"))

rarefied.unifrac.cyanobacteria.dist <- as.matrix(phyloseq::distance(rarefied.cyanobacteria.prop.exclude.1.3, method = "unifrac"))

rarefied.unifrac.proteobacteria.otherphoto.dist <- as.matrix(phyloseq::distance(rarefied.proteobacteria.otherphoto.prop.exclude.1.3, method = "unifrac"))

rarefied.unifrac.photosynthetic.dist <- as.matrix(phyloseq::distance(rarefied.photosynthetic.prop.exclude.1.3, method = "unifrac"))

rarefied.unifrac.chemolithotrophs.dist <- as.matrix(phyloseq::distance(rarefied.chemolithotrophs.prop.exclude.1.3, method = "unifrac"))

rarefied.unifrac.only.heterotrophs.dist <- as.matrix(phyloseq::distance(rarefied.only.heterotrophs.prop.exclude.1.3, method = "unifrac"))

### get gps co-ordinates ####
### use metadata from any phyloseq 
metadata.dd <- as.matrix(as.data.frame(com.rarefied.min.prop.exclude.1.3@sam_data)) 
metadata.dd <- as.data.frame(metadata.dd) %>% rownames_to_column(var = "sample_id")  %>% as_tibble() 
gps.all <- metadata.dd %>% select(sample_id,Longitude,Latittude) %>% column_to_rownames(var = "sample_id")
gps.all[1:ncol(gps.all)] <- lapply(gps.all[1:ncol(gps.all)], as.numeric)

#Vector of distances in km for all taxa
distgeo.all <- geodist(gps.all, measure = 'geodesic')/1000
rownames(distgeo.all) <- rownames(gps.all)
colnames(distgeo.all) <- rownames(gps.all)

### to get maximum distance
max(distgeo.all)

## melt distance matrices
rarefied.unifrac.dist.df <- melt(as.matrix(rarefied.unifrac.dist), varnames = c("row", "col"))

rarefied.unifrac.chloroflexi.dist.df <- melt(as.matrix(rarefied.unifrac.chloroflexi.dist), varnames = c("row", "col"))

rarefied.unifrac.cyanobacteria.dist.df <- melt(as.matrix(rarefied.unifrac.cyanobacteria.dist), varnames = c("row", "col"))

rarefied.unifrac.proteobacteria.otherphoto.dist.df <- melt(as.matrix(rarefied.unifrac.proteobacteria.otherphoto.dist), varnames = c("row", "col"))

rarefied.unifrac.photosynthetic.dist.df <- melt(as.matrix(rarefied.unifrac.photosynthetic.dist), varnames = c("row", "col"))

rarefied.unifrac.chemolithotrophs.dist.df <- melt(as.matrix(rarefied.unifrac.chemolithotrophs.dist), varnames = c("row", "col"))

rarefied.unifrac.only.heterotrophs.dist.df <- melt(as.matrix(rarefied.unifrac.only.heterotrophs.dist), varnames = c("row", "col"))

distgeo.all.df <- melt(as.matrix(distgeo.all), varnames = c("row", "col"))

### Add unique column names for distance matrix dataframes for better identification
colnames(rarefied.unifrac.dist.df)[3] <- "unweighted.unifrac.alltaxa"
colnames(rarefied.unifrac.chloroflexi.dist.df)[3] <- "unweighted.unifrac.chloroflexi"
colnames(rarefied.unifrac.cyanobacteria.dist.df)[3] <- "unweighted.unifrac.cyanobacteria"
colnames(rarefied.unifrac.proteobacteria.otherphoto.dist.df)[3] <- "unweighted.unifrac.proteobacteria.otherphoto"
colnames(rarefied.unifrac.photosynthetic.dist.df)[3] <- "unweighted.unifrac.all.photosynthetic"
colnames(rarefied.unifrac.chemolithotrophs.dist.df)[3] <- "unweighted.unifrac.chemolithotrophic"
colnames(rarefied.unifrac.only.heterotrophs.dist.df)[3] <- "unweighted.unifrac.heterotrophic"
colnames(distgeo.all.df)[3] <- "Kmdist"

### MERGE distance matrix and gps df
mergedist.unifrac <- merge(rarefied.unifrac.dist.df,distgeo.all.df)
mergedist.unifrac.chloroflexi <- merge(rarefied.unifrac.chloroflexi.dist.df,distgeo.all.df)
mergedist.unifrac.cyanobacteria <- merge(rarefied.unifrac.cyanobacteria.dist.df,distgeo.all.df)
mergedist.unifrac.proteobacteria.otherphoto <- merge(rarefied.unifrac.proteobacteria.otherphoto.dist.df,distgeo.all.df)
mergedist.unifrac.photosynthetic <- merge(rarefied.unifrac.photosynthetic.dist.df,distgeo.all.df)
mergedist.unifrac.chemolithotrophs <- merge(rarefied.unifrac.chemolithotrophs.dist.df,distgeo.all.df)
mergedist.unifrac.only.heterotrophs <- merge(rarefied.unifrac.only.heterotrophs.dist.df,distgeo.all.df)

### plot distance deacays ####

distdecay.unifrac <- mergedist.unifrac %>% ggplot(aes(x=Kmdist, y=unweighted.unifrac.alltaxa)) + geom_point(size = 1, alpha = 0.1, position = position_dodge(0.75)) + geom_smooth(formula = y ~ x, method = "lm", se = TRUE, color = "maroon") + stat_cor(label.x = 1000, label.y = 0.3) + stat_regline_equation(label.x = 1000, label.y = 0.25)  + theme(axis.text.x =element_text(size=14, face = "bold", angle = 0 , hjust = 0.5),axis.text.y =element_text(size=12, face = "bold"),axis.title.y=element_text(size = 12, face = "bold"),axis.title.x = element_text(size = 12, face = "bold")) +  theme(legend.position = "none") +labs(title="Distance deecay of Unweighted UniFrac distance (all taxa) against pairwise geographical distance in Km ", y = "Unweighted UniFrac dissimilarity", x = "Geographical distance (Km)")+ theme(strip.text.y = element_text(angle = 0, face = "bold", size = 12)) + theme(plot.title = element_text(hjust = 0.5, size = 11, face = "bold")) + ggeasy::easy_center_title() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) 

ggsave("distdecay.unifrac.pdf")
ggsave("distdecay.unifrac.png")


distdecay.unifrac.chloroflexi <- mergedist.unifrac.chloroflexi %>% ggplot(aes(x=Kmdist, y=unweighted.unifrac.chloroflexi)) + geom_point(size = 1, alpha = 0.1, position = position_dodge(0.75)) + geom_smooth(formula = y ~ x, method = "lm", se = TRUE, color = "maroon") + stat_cor(label.x = 1000, label.y = 0.7) + stat_regline_equation(label.x = 1000, label.y = 0.65)  + theme(axis.text.x =element_text(size=14, face = "bold", angle = 0 , hjust = 0.5),axis.text.y =element_text(size=12, face = "bold"),axis.title.y=element_text(size = 12, face = "bold"),axis.title.x = element_text(size = 12, face = "bold")) +  theme(legend.position = "none") +labs(title="Distance deecay of Unweighted UniFrac distance (chloroflexi) against pairwise geographical distance in Km ", y = "Unweighted UniFrac dissimilarity", x = "Geographical distance (Km)")+ theme(strip.text.y = element_text(angle = 0, face = "bold", size = 12)) + theme(plot.title = element_text(hjust = 0.5, size = 11, face = "bold")) + ggeasy::easy_center_title() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) 

ggsave("distdecay.unifrac.chloroflexi.pdf")
ggsave("distdecay.unifrac.chloroflexi.png")


distdecay.unifrac.cyanobacteria <- mergedist.unifrac.cyanobacteria %>% ggplot(aes(x=Kmdist, y=unweighted.unifrac.cyanobacteria)) + geom_point(size = 1, alpha = 0.1, position = position_dodge(0.75)) + geom_smooth(formula = y ~ x, method = "lm", se = TRUE, color = "maroon") + stat_cor(label.x = 1000, label.y = 0.15) + stat_regline_equation(label.x = 1000, label.y = 0.10)  + theme(axis.text.x =element_text(size=14, face = "bold", angle = 0 , hjust = 0.5),axis.text.y =element_text(size=12, face = "bold"),axis.title.y=element_text(size = 12, face = "bold"),axis.title.x = element_text(size = 12, face = "bold")) +  theme(legend.position = "none") +labs(title="Distance deecay of Unweighted UniFrac distance (cyanobacteria) against pairwise geographical distance in Km ", y = "Unweighted UniFrac dissimilarity", x = "Geographical distance (Km)")+ theme(strip.text.y = element_text(angle = 0, face = "bold", size = 12)) + theme(plot.title = element_text(hjust = 0.5, size = 11, face = "bold")) + ggeasy::easy_center_title() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) 

ggsave("distdecay.unifrac.cyanobacteria.pdf")
ggsave("distdecay.unifrac.cyanobacteria.png")


distdecay.unifrac.proteobacteria.otherphoto <- mergedist.unifrac.proteobacteria.otherphoto %>% ggplot(aes(x=Kmdist, y=unweighted.unifrac.proteobacteria.otherphoto)) + geom_point(size = 1, alpha = 0.1, position = position_dodge(0.75)) + geom_smooth(formula = y ~ x, method = "lm", se = TRUE, color = "maroon") + stat_cor(label.x = 1000, label.y = 0.78) + stat_regline_equation(label.x = 1000, label.y = 0.72)  + theme(axis.text.x =element_text(size=14, face = "bold", angle = 0 , hjust = 0.5),axis.text.y =element_text(size=12, face = "bold"),axis.title.y=element_text(size = 12, face = "bold"),axis.title.x = element_text(size = 12, face = "bold")) +  theme(legend.position = "none") +labs(title="Distance deecay of Unweighted UniFrac distance (proteobacteria and other photosynthetic taxa) against pairwise geographical distance in Km ", y = "Unweighted UniFrac dissimilarity", x = "Geographical distance (Km)")+ theme(strip.text.y = element_text(angle = 0, face = "bold", size = 12)) + theme(plot.title = element_text(hjust = 0.5, size = 11, face = "bold")) + ggeasy::easy_center_title() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) 

ggsave("distdecay.unifrac.proteobacteria.otherphoto.pdf")
ggsave("distdecay.unifrac.proteobacteria.otherphoto.png")


distdecay.unifrac.photosynthetic <- mergedist.unifrac.photosynthetic %>% ggplot(aes(x=Kmdist, y=unweighted.unifrac.all.photosynthetic)) + geom_point(size = 1, alpha = 0.1, position = position_dodge(0.75)) + geom_smooth(formula = y ~ x, method = "lm", se = TRUE, color = "maroon") + stat_cor(label.x = 1000, label.y = 0.12) + stat_regline_equation(label.x = 1000, label.y = 0.09)  + theme(axis.text.x =element_text(size=14, face = "bold", angle = 0 , hjust = 0.5),axis.text.y =element_text(size=12, face = "bold"),axis.title.y=element_text(size = 12, face = "bold"),axis.title.x = element_text(size = 12, face = "bold")) +  theme(legend.position = "none") +labs(title="Distance deecay of Unweighted UniFrac distance (all photosynthetic taxa) against pairwise geographical distance in Km ", y = "Unweighted UniFrac dissimilarity", x = "Geographical distance (Km)")+ theme(strip.text.y = element_text(angle = 0, face = "bold", size = 12)) + theme(plot.title = element_text(hjust = 0.5, size = 11, face = "bold")) + ggeasy::easy_center_title() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) 

ggsave("distdecay.unifrac.photosynthetic.pdf")
ggsave("distdecay.unifrac.photosynthetic.png")


distdecay.unifrac.chemolithotrophs <- mergedist.unifrac.chemolithotrophs %>% ggplot(aes(x=Kmdist, y=unweighted.unifrac.chemolithotrophic)) + geom_point(size = 1, alpha = 0.1, position = position_dodge(0.75)) + geom_smooth(formula = y ~ x, method = "lm", se = TRUE, color = "maroon") + stat_cor(label.x = 1200, label.y = 0.30) + stat_regline_equation(label.x = 1200, label.y = 0.25)  + theme(axis.text.x =element_text(size=14, face = "bold", angle = 0 , hjust = 0.5),axis.text.y =element_text(size=12, face = "bold"),axis.title.y=element_text(size = 12, face = "bold"),axis.title.x = element_text(size = 12, face = "bold")) +  theme(legend.position = "none") +labs(title="Distance deecay of Unweighted UniFrac distance (Chemolithotrophs) against pairwise geographical distance in Km ", y = "Unweighted UniFrac dissimilarity", x = "Geographical distance (Km)")+ theme(strip.text.y = element_text(angle = 0, face = "bold", size = 12)) + theme(plot.title = element_text(hjust = 0.5, size = 11, face = "bold")) + ggeasy::easy_center_title() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) 

ggsave("distdecay.unifrac.chemolithotrophs.pdf")
ggsave("distdecay.unifrac.chemolithotrophs.png")


distdecay.unifrac.only.heterotrophs <- mergedist.unifrac.only.heterotrophs %>% ggplot(aes(x=Kmdist, y=unweighted.unifrac.heterotrophic)) + geom_point(size = 1, alpha = 0.1, position = position_dodge(0.75)) + geom_smooth(formula = y ~ x, method = "lm", se = TRUE, color = "maroon") + stat_cor(label.x = 1000, label.y = 0.30) + stat_regline_equation(label.x = 1000, label.y = 0.25)  + theme(axis.text.x =element_text(size=14, face = "bold", angle = 0 , hjust = 0.5),axis.text.y =element_text(size=12, face = "bold"),axis.title.y=element_text(size = 12, face = "bold"),axis.title.x = element_text(size = 12, face = "bold")) +  theme(legend.position = "none") +labs(title="Distance deecay of Unweighted UniFrac distance (only heterotrophs) against pairwise geographical distance in Km ", y = "Unweighted UniFrac dissimilarity", x = "Geographical distance (Km)")+ theme(strip.text.y = element_text(angle = 0, face = "bold", size = 12)) + theme(plot.title = element_text(hjust = 0.5, size = 11, face = "bold")) + ggeasy::easy_center_title() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) 

ggsave("distdecay.unifrac.only.heterotrophs.pdf")
ggsave("distdecay.unifrac.only.heterotrophs.png")


####################################################################################################################

### Fig. 2A. Beta diversity patterns for hot spring photosynthetic biofilms indicated statistically supported biogeographic regions. ####
### Beta diversity was estimated using Principle Coordinate Analysis of unweighted UniFRAC distances

######################################################################################################################

### Set color and shape legends ####

# to get variables from phyloseq metadata
sample_variables(com.rarefied.min.prop.exclude.1.3)

#get variable levels
get_variable(com.rarefied.min.prop.exclude.1.3,"Temp.adj...C.")
get_variable(com.rarefied.min.prop.exclude.1.3,"Location.Code")

## Assign shapes to Location.Codes 
# Location codes : AP BA HS KJ LA LN MK PB PP PT RB RN SE SW US
loc.shape.names = c('AP', 'BA', 'HS', 'KJ', 'LA', 'LN', 'MK', 'PB', 'PP', 'PT', 'RB', 'RN', 'SE', 'SW', 'US')
loc.shape <- 0:(length(loc.shape.names)-1)
names(loc.shape) <- loc.shape.names
loc.shape["Taxa"] <- 16
loc.shape

## Assign colors to Temperatures
temp.col = c('38' = "mediumorchid1", '40' = "mediumpurple1", '42' = "blueviolet",'44' = "mediumpurple4", '45' = "lightskyblue2",  '46' = "deepskyblue" , '47'= "royalblue1",  '48' = "mediumblue", '49' = "cyan", '50' = "aquamarine2" , '51' = "springgreen2",  '53' = "yellowgreen", '55' = "forestgreen", '57' = "khaki3", '58' = "goldenrod1", '61' = "sienna2", '63' = "firebrick3", '66' = "darkred", 'Taxa' = "black")

##use temp breaks to specify temp legend first and then location
temp.breaks = c("38", "40", "42", "44" , "45", "46", "47", "48" , "49", "50", "51", "53", "55", "57", "58", "61" , "63" , "66","Taxa")

## Assign colors to Regions - for ellipses
region.col = c('North.Thailand' = "slateblue", 'Central.Thailand' = "turquoise3", 'South.Thailand' = "gold", 'North.Malaysia' = "chocolate1", 'South.Malaysia' = "lawngreen", 'Singapore'="palevioletred3")

#### Calcluate Unweighted UniFRAC ordiation ####

## In phyloseq : 
## "unifrac" :Original (unweighted) UniFrac distance, UniFrac
## "wunifrac" :weighted-UniFrac distance, UniFrac
## "dpcoa": sample-wise distance used in Double Principle Coordinate Analysis, DPCoA

## MDS/PCoA: Performs principal coordinate analysis (also called principle coordinate decomposition, multidimensional scaling (MDS), or classical scaling) of a distance matrix (Gower 1966), including two correction methods for negative eigenvalues. See pcoa for further details.
rare.unweighted.unifrac.PCoA = ordinate(com.rarefied.min.prop.exclude.1.3, method="PCoA", distance="unifrac", weighted=FALSE)

##Plot

clean_background <- theme(plot.background = element_rect("white"),
                          panel.background = element_rect(color = "black", fill = NA, linewidth = 2),
                          panel.grid = element_line("white"),
                          axis.line = element_line("gray25"),
                          axis.text = element_text(size = 12, color = "gray25"),
                          axis.title = element_text(color = "gray25"),
                          legend.text = element_text(size = 12),
                          legend.key = element_rect("white"))

rarefied.PCoA.UnweightedUniFrac.phyloseq.regionellipse.samples <- plot_ordination(com.rarefied.min.prop.exclude.1.3, rare.unweighted.unifrac.PCoA, type="samples", color = "Temp.adj...C.", shape = "Location.Code", title="PCoA on unweighted UniFrac distance") +   geom_point(size=7, position="jitter") + stat_ellipse(aes(fill = Region, group = Region), linetype = 0 ,type = "t",  level = 0.99, alpha = .2, geom = "polygon")   + clean_background  + scale_shape_manual(values= loc.shape)  + scale_color_manual(values = temp.col) + scale_fill_manual(values = region.col,breaks = c('North.Thailand','Central.Thailand','South.Thailand','North.Malaysia','South.Malaysia','Singapore')) + theme(legend.title = element_text(size = 15.5, face = "bold"), legend.text = element_text(size = 14),legend.position = "bottom", legend.direction = "horizontal",legend.box="vertical", legend.margin=margin()) + guides(colour = guide_legend(nrow=1,byrow=TRUE, override.aes = list(size = 7)), shape = guide_legend(nrow=1,byrow=TRUE,override.aes = list(size = 7)), fill = guide_legend(nrow=1,byrow=TRUE,override.aes = list(size = 7)) )

#### to remove double/open shapes during ordination which is only apparent with open lables and not closed or filed labels
rarefied.PCoA.UnweightedUniFrac.phyloseq.regionellipse.samples$layers <- rarefied.PCoA.UnweightedUniFrac.phyloseq.regionellipse.samples$layers[-1]

pdf("rarefied.PCoA.UnweightedUniFrac.phyloseq.regionellipse.samples.pdf", height = 12, width = 17)
rarefied.PCoA.UnweightedUniFrac.phyloseq.regionellipse.samples
dev.off()


## Significance testing on beta diversity ####

## get metadata
# make a data frame from the sample_data
sampledf <- data.frame(metadata)

## betadisper :: Homogeneity of  multivariavte dispersion test 
betadisp.region <- betadisper(rarefied.unifrac.dist, sampledf$Region)
anova(betadisp.region)

## Tukey's test on betadisper
betadisp.region.HSD <- TukeyHSD(betadisp.region)
betadisp.region.HSD

# Adonis test -  permanova
adonis.region <- adonis2(rarefied.unifrac.dist ~ Region, data = sampledf)


####################################################################################################################

### Fig. 2B. Phylogenetic diversity and relative abundance of major bacterial and archaeal lineages in biofilms (grey tree), and pairwise comparisons of lineages that were over-represented (blue branches and nodes) or under-represented (orange branches and nodes) between regions ####
### Fig. S5. Heat tree indicating the phylogenetic composition of photosynthetic biofilm community from Southeast Asian hot springs (N = 395). Branch thickness and node size denote relative abundance. #### 

######################################################################################################################

otu <- as.matrix(as.data.frame(otu_table(com.rarefied.min.int.exclude.1.3.rooted)))
otu <- t(as.data.frame(otu) %>% rownames_to_column(var = "sample_id")  %>% as_tibble())
otu <- janitor::row_to_names(otu, 1)
otu <- as.data.frame(otu) %>% rownames_to_column(var = "ASV")  %>% as_tibble()
sapply(otu,class)
otu[2:ncol(otu)] <- lapply(otu[2:ncol(otu)], as.numeric)

tax <- as.matrix(as.data.frame(tax_table(com.rarefied.min.int.exclude.1.3.rooted)))
tax <- as.data.frame(tax) %>% rownames_to_column(var = "ASV")  %>% as_tibble()
#tax <- tax %>% unite(Taxonomy, c(Kingdom, Phylum, Class, Order, Family, Genus, Species), sep = ";", remove = FALSE)

otu_tax <- left_join(otu, tax, by = "ASV")
sapply(otu_tax,class)

# Create tax_map class and class_cols = "Phylum" is used as taxon_names
#otu_tax_taxmap <- parse_tax_data(otu_tax,class_cols = "Taxonomy", class_sep = ";")
otu_tax_taxmap <- parse_tax_data(otu_tax,class_cols = c("Kingdom","Phylum","Class","Order","Family","Genus","Species"),named_by_rank = TRUE)

#changes names of the tble from default tax_data to otu counts.
names(otu_tax_taxmap$data) <- "otu_counts_taxa"

metadata <- (sample_data(com.rarefied.min.int.exclude.1.3.rooted))

#The combination of “reingold-tilford” and “davidson-harel” to be space-efficient for large plots.
#Supertaxa: The taxa a taxon is contained within. For example, Homo is a supertaxon of the species Homo sapiens

# Setting edge attributes: Just like the node size/color, the edge size/color can be used to display statistics. Edges can also have labels. This means you can plot 4 statistics in one graph, although more than 3 is typically overwhelmin

otu_tax_taxmap$data$tax_abund <- calc_taxon_abund(otu_tax_taxmap, "otu_counts_taxa",
                                                  cols = metadata$Sample_id,
                                                  groups = metadata$Region)

### Use taxmap2 to add read abundances to legend ######

otu_tax_taxmap2 <- otu_tax_taxmap 

otu_tax_taxmap2$data$tax_abund <- calc_taxon_abund(otu_tax_taxmap2, "otu_counts_taxa",
                                                   cols = metadata$Sample_id)

otu_tax_taxmap2$data$tax_occ <- calc_n_samples(otu_tax_taxmap2, "tax_abund", cols = metadata$Sample_id)

set.seed(3)
otu_tax_taxmap2 %>% 
  mutate_obs("tax_abund", abundance = rowSums(otu_tax_taxmap2$data$tax_abund[metadata[["Sample_id"]]])) %>%
  filter_taxa(taxon_ranks == "Species", supertaxa = TRUE) %>%
  heat_tree(node_label = taxon_names,
            node_size = n_obs,
            node_color = abundance,
            #edge_label = n_obs,
            initial_layout = "re",
            layout = "da",
            title = "ASV diversity",
            node_color_axis_label = "Read abundance",
            node_size_axis_label = "Number of OTUs",
            output_file = "Species.heattree.supplementary.labels.pdf")


#Comparing different regions 
#For each taxon, a Wilcoxon Rank Sum test was used to test for differences between the median abundances of samples in each treatment. We can use this information to create what we call a differential heat tree,

otu_tax_taxmap$data$diff_table <- compare_groups(otu_tax_taxmap, "tax_abund",
                                                 cols = metadata$Sample_id, # What columns of sample data to use
                                                 groups = metadata$Region) # What category each sample is assigned to


set.seed(1)
heat_tree_matrix(otu_tax_taxmap,
                 data = "diff_table",
                 node_size = n_obs, # n_obs is a function that calculates, in this case, the number of OTUs per taxon
                 #node_label = taxon_names,
                 node_color = log2_median_ratio, # A column from `obj$data$diff_table`
                 node_color_range = diverging_palette(), # The built-in palette for diverging data
                 node_color_trans = "linear", # The default is scaled by circle area
                 node_color_interval = c(-3, 3), # The range of `log2_median_ratio` to display
                 edge_color_interval = c(-3, 3), # The range of `log2_median_ratio` to display
                 node_size_axis_label = "Number of OTUs",
                 node_color_axis_label = "Log2 ratio median proportions",
                 layout = "davidson-harel", # The primary layout algorithm
                 initial_layout = "reingold-tilford", # The layout algorithm that initializes node locations
                 output_file = "differential_heat_tree.nolabels.pdf") # Saves the plot as a pdf file


####################################################################################################################

### Fig. S2.Alpha diversity metrics for the 40 hot spring locations sampled in this study. #####
### Samples are displayed in order of ascending temperature. ####

######################################################################################################################

## Get metadata 
rarefied.min.metadata.3 <- as.data.frame(metadata)

rarefied.min.metadata.3$Temp.adj...C. <- as.factor(rarefied.min.metadata.3$Temp.adj...C.)
rarefied.min.metadata.3$Longitude <- as.factor(rarefied.min.metadata.3$Longitude)
rarefied.min.metadata.3$Latittude <- as.factor(rarefied.min.metadata.3$Latittude)
rarefied.min.metadata.3$EC..mS. <- as.factor(rarefied.min.metadata.3$EC..mS.)
rarefied.min.metadata.3$pH <- as.factor(rarefied.min.metadata.3$pH)
rarefied.min.metadata.3$Carbonate..ppm. <- as.factor(rarefied.min.metadata.3$Carbonate..ppm.)

## change column names
names(rarefied.min.metadata.3)[names(rarefied.min.metadata.3) == "Temp.adj...C."] <- "Temp.adj (°C)"
names(rarefied.min.metadata.3)[names(rarefied.min.metadata.3) == "Temp...C."] <- "Temp.actual (°C)"

## Get variables for alpha diversity indices table
Site.rarefied.min.3 <- rarefied.min.metadata.3$Site
Sample_id.rarefied.min.3 <- rarefied.min.metadata.3$Sample_id
Temp.rarefied.min.3 <- rarefied.min.metadata.3$`Temp.adj (°C)`
Long.rarefied.min.3 <- rarefied.min.metadata.3$Longitude
Lat.rarefied.min.3 <- rarefied.min.metadata.3$Latittude
Cond.rarefied.min.3 <- rarefied.min.metadata.3$EC..mS. 
pH.rarefied.min.3 <- rarefied.min.metadata.3$pH 
carb.rarefied.min.3 <- rarefied.min.metadata.3$Carbonate..ppm.

### Calculate alpha diversity indices from otu table

### calculte Shannon Index (H)
rarefied.min.shannon.3 <- diversity(rarefied.min.wide.3, index="shannon")

## Pielou's eveness (J) for Shannon index
rarefied.min.pielou.3 <- diversity(rarefied.min.wide.3,index="shannon")/log(specnumber(rarefied.min.wide.3))

## Simpson's index (D,λ) and Gini-Simpson(GS)=1-λ
rarefied.min.simpson.3 <- diversity(rarefied.min.wide.3, index="simpson")

## Chao estimator (specpool and estimateR)
rarefied.min.chao.3 <- as.data.frame(estimateR(rarefied.min.wide.3))

### Create alpha diveristy indices table
diversity.indices.rarefied.min.3 <- data.frame(Site.rarefied.min.3,Sample_id.rarefied.min.3,Temp.rarefied.min.3, Lat.rarefied.min.3, Cond.rarefied.min.3, pH.rarefied.min.3,carb.rarefied.min.3)
diversity.indices.rarefied.min.3$rarefied.min.shannon.3 <- rarefied.min.shannon.3
diversity.indices.rarefied.min.3$rarefied.min.pielou.3 <- rarefied.min.pielou.3
diversity.indices.rarefied.min.3$rarefied.min.simpson.3 <- rarefied.min.simpson.3
diversity.indices.rarefied.min.3$rarefied.min.chao.3 <- t(rarefied.min.chao.3["S.chao1", ])

### Change x-axis variable to character for melt and later change them as numeric after getting long table

diversity.indices.rarefied.min.3$Temp.rarefied.min.3 <- as.character(diversity.indices.rarefied.min.3$Temp.rarefied.min.3)

diversity.indices.rarefied.min.3$Lat.rarefied.min.3 <- as.character(diversity.indices.rarefied.min.3$Lat.rarefied.min.3)

diversity.indices.rarefied.min.3$Cond.rarefied.min.3 <- as.character(diversity.indices.rarefied.min.3$Cond.rarefied.min.3)

diversity.indices.rarefied.min.3$pH.rarefied.min.3 <- as.character(diversity.indices.rarefied.min.3$pH.rarefied.min.3)

diversity.indices.rarefied.min.3$carb.rarefied.min.3 <- as.character(diversity.indices.rarefied.min.3$carb.rarefied.min.3)

rownames(diversity.indices.rarefied.min.3) <- NULL

setDT(diversity.indices.rarefied.min.3)

diversity.indices.rarefied.min.long.3  <- melt(diversity.indices.rarefied.min.3[1:nrow(diversity.indices.rarefied.min.3),])

### Create new facet label names for variables
new_var_names_rarefied_min_3 <- c(
  'rarefied.min.chao.3'="Chao1 estimator",
  'rarefied.min.shannon.3'="Shannon's index",
  'rarefied.min.pielou.3'="Pielou's evenness",
  'rarefied.min.simpson.3'="Gini-Simpson's index"
)

### Assign numeric variables and order table

## Temperature
diversity.indices.rarefied.min.long.3$Temp.rarefied.min.3 <- as.numeric(diversity.indices.rarefied.min.long.3$Temp.rarefied.min.3)

diversity.indices.rarefied.min.long.3 <- arrange(diversity.indices.rarefied.min.long.3, (Temp.rarefied.min.3), desc(value),)

## Latitude
diversity.indices.rarefied.min.long.3$Lat.rarefied.min.3 <- as.numeric(diversity.indices.rarefied.min.long.3$Lat.rarefied.min.3)

## Assign column number accordingly : "[,4]"
diversity.indices.rarefied.min.long.3[,4] <- round(diversity.indices.rarefied.min.long.3[,4], digits = 4)

diversity.indices.rarefied.min.long.3 <- arrange(diversity.indices.rarefied.min.long.3, desc(Lat.rarefied.min.3), (value),)


## Conductivity
diversity.indices.rarefied.min.long.3$Cond.rarefied.min.3 <- as.numeric(diversity.indices.rarefied.min.long.3$Cond.rarefied.min.3)

## Assign column number accordingly : "[,5]"
diversity.indices.rarefied.min.long.3[,5] <- round(diversity.indices.rarefied.min.long.3[,5], digits = 2)

diversity.indices.rarefied.min.long.3 <- arrange(diversity.indices.rarefied.min.long.3, desc(Cond.rarefied.min.3), (value),)

## pH
diversity.indices.rarefied.min.long.3$pH.rarefied.min.3 <- as.numeric(diversity.indices.rarefied.min.long.3$pH.rarefied.min.3)

## Assign column number accordingly : "[,6]"
diversity.indices.rarefied.min.long.3[,6] <- round(diversity.indices.rarefied.min.long.3[,6], digits = 2)

diversity.indices.rarefied.min.long.3 <- arrange(diversity.indices.rarefied.min.long.3, desc(pH.rarefied.min.3), (value),)

## Carbonate
diversity.indices.rarefied.min.long.3$carb.rarefied.min.3 <- as.numeric(diversity.indices.rarefied.min.long.3$carb.rarefied.min.3)

diversity.indices.rarefied.min.long.3 <- arrange(diversity.indices.rarefied.min.long.3, desc(carb.rarefied.min.3), (value),)

## plots 

## By Site and Temperature
rarefied.min.alpha.div.indices.loc.notempfacet.3 <- diversity.indices.rarefied.min.long.3  %>% filter(!variable %in% "rarefied.min.no.3") %>% mutate(variable = factor(variable, levels=c("rarefied.min.chao.3","rarefied.min.shannon.3","rarefied.min.simpson.3","rarefied.min.pielou.3"))) %>% mutate(Site.rarefied.min.3 = factor(Site.rarefied.min.3, levels= c("PB_38" , "RB_40" , "RN_42" , "HS_42" , "MK_42" , "MK_44" , "RN_45", "HS_45" , "RB_45", "SE_45", "SW_46" , "PB_46"  , "LN_46", "AP_46" , "LA_46" , "KJ_47" , "PP_48" , "SE_49", "PT_50", "MK_50"  , "BA_50"  , "RB_50" , "RN_51" , "HS_51" , "LN_51" , "SW_51" , "HS_53" , "PP_53" , "RB_55", "PT_55" , "SW_55" , "MK_55" , "US_57" , "PP_58" , "RN_61", "SW_61", "LN_63","US_63" , "PT_63" , "LN_66"))) %>% drop_na() %>% ggplot(aes(x=Site.rarefied.min.3, y=value, fill=variable)) + geom_boxplot(outlier.colour = "red") + stat_summary(fun=mean, geom="point", color="black", size=2, shape=23, position=position_dodge(0.75)) +  theme(axis.text.x =element_text(size=10, face = "bold", angle = 45 , hjust = 1),axis.text.y =element_text(size=12, face = "bold"),axis.title.y=element_text(size = 12, face = "bold"),axis.title.x = element_blank()) +  theme(legend.position = "none") +labs(title="Diversity indices of rarefied ASVs with error model 3", y = "Diversity index") + facet_grid(vars(variable), scales = "free",labeller = as_labeller(new_var_names_rarefied_min_3)) + theme(strip.text.y = element_text(angle = 0, face = "bold", size = 12)) + theme(plot.title = element_text(hjust = 0.5, size = 13, face = "bold")) + ggeasy::easy_center_title() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + geom_smooth(formula = y ~ x, method = "lm", se = TRUE, color = "grey40", aes(group = variable))  + stat_regline_equation(label.x = 1, label.y = 0.7,aes(group = variable)) + stat_cor(label.x = 1, label.y = 1, aes(group = variable)) 

ggsave("rarefied.min.alpha.div.indices.loc.notempfacet.3.pdf",height = 12, width = 20)
ggsave("rarefied.min.alpha.div.indices.loc.notempfacet.3.png", height = 12, width = 20)

## to get equations clearly for chao as it cannot be seen in plot with all indices
### Adjust accordingly for all variables
rarefied.min.chao.loc.notempfacet.3 <- diversity.indices.rarefied.min.long.3  %>% filter(variable == "rarefied.min.chao.3")  %>% mutate(Site.rarefied.min.3 = factor(Site.rarefied.min.3, levels= c("PB_38" , "RB_40" , "RN_42" , "HS_42" , "MK_42" , "MK_44" , "RN_45", "HS_45" , "RB_45", "SE_45", "SW_46" , "PB_46"  , "LN_46", "AP_46" , "LA_46" , "KJ_47" , "PP_48" , "SE_49", "PT_50", "MK_50"  , "BA_50"  , "RB_50" , "RN_51" , "HS_51" , "LN_51" , "SW_51" , "HS_53" , "PP_53" , "RB_55", "PT_55" , "SW_55" , "MK_55" , "US_57" , "PP_58" , "RN_61", "SW_61", "LN_63","US_63" , "PT_63" , "LN_66"))) %>% drop_na() %>% ggplot(aes(x=Site.rarefied.min.3, y=value, fill=variable)) + geom_boxplot(outlier.colour = "red") + stat_summary(fun=mean, geom="point", color="black", size=2, shape=23, position=position_dodge(0.75)) +  theme(axis.text.x =element_text(size=10, face = "bold", angle = 45 , hjust = 1),axis.text.y =element_text(size=12, face = "bold"),axis.title.y=element_text(size = 12, face = "bold"),axis.title.x = element_blank()) +  theme(legend.position = "none") +labs(title="Diversity indices of rarefied ASVs with error model 3", y = "Diversity index") + facet_grid(vars(variable), scales = "free",labeller = as_labeller(new_var_names_rarefied_min_3)) + theme(strip.text.y = element_text(angle = 0, face = "bold", size = 12)) + theme(plot.title = element_text(hjust = 0.5, size = 13, face = "bold")) + ggeasy::easy_center_title() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + geom_smooth(formula = y ~ x, method = "lm", se = TRUE, color = "grey40", aes(group = variable))  + stat_regline_equation(label.x = 1, label.y =2750 ,aes(group = variable)) + stat_cor(label.x = 1, label.y = 3000, aes(group = variable)) 

## By latitude
rarefied.min.alpha.div.indices.lat.3 <- diversity.indices.rarefied.min.long.3   %>% filter(!variable == "rarefied.min.no.3") %>% mutate(variable = factor(variable, levels=c("rarefied.min.chao.3","rarefied.min.shannon.3","rarefied.min.simpson.3","rarefied.min.pielou.3"))) %>% drop_na() %>% ggplot(aes(x=reorder(factor(as.numeric(Lat.rarefied.min.3)), -(as.numeric(Lat.rarefied.min.3))), y=value, fill=variable)) + geom_point(size=0.3) + stat_summary(fun=mean, geom="point", color="red", size=0.4, shape=23, position=position_dodge(0.75)) +  theme(axis.text.x =element_text(size=6, face = "bold", angle = 45 , hjust = 1),axis.text.y =element_text(size=6, face = "bold"),axis.title.y=element_blank(),axis.title.x = element_blank()) +  theme(legend.position = "none") + facet_wrap(vars(variable), scales = "free_y",labeller = as_labeller(new_var_names_rarefied_min_3),  ncol=1) + theme(strip.text.y = element_text(angle = 0, face = "bold", size = 8)) + theme(plot.title = element_blank()) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + geom_smooth(formula = y ~ x, method = "lm", se = TRUE, color = "grey40", aes(group = variable)) + stat_regline_equation(label.x = 1, label.y =0.7 ,aes(group = variable)) + stat_cor(label.x = 1, label.y = 1, aes(group = variable)) 

ggsave("rarefied.min.alpha.div.indices.lat.3.pdf",height = 11, width = 5)
ggsave("rarefied.min.alpha.div.indices.lat.3.png", height = 11, width = 3)

## By pH
rarefied.min.alpha.div.indices.pH.3 <- diversity.indices.rarefied.min.long.3  %>% filter(!variable == "rarefied.min.no.3") %>% mutate(variable = factor(variable, levels=c("rarefied.min.chao.3","rarefied.min.shannon.3","rarefied.min.simpson.3","rarefied.min.pielou.3"))) %>% drop_na() %>% ggplot(aes(x=reorder(factor(pH.rarefied.min.3), pH.rarefied.min.3), y=value, fill=variable)) + geom_point(size=0.3) + stat_summary(fun=mean, geom="point", color="red", size=0.4, shape=23, position=position_dodge(0.75)) +  theme(axis.text.x =element_text(size=6, face = "bold", angle = 45 , hjust = 1),axis.text.y =element_text(size=6, face = "bold"),axis.title.y=element_blank(),axis.title.x = element_blank()) +  theme(legend.position = "none") + facet_wrap(vars(variable), scales = "free_y",labeller = as_labeller(new_var_names_rarefied_min_3),  ncol=1) + theme(strip.text.y = element_text(angle = 0, face = "bold", size = 8)) + theme(plot.title = element_blank()) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + geom_smooth(formula = y ~ x, method = "lm", se = TRUE, color = "grey40", aes(group = variable))  + stat_regline_equation(label.x = 1, label.y =2000 ,aes(group = variable)) + stat_cor(label.x = 1, label.y = 1000, aes(group = variable)) 

ggsave("rarefied.min.alpha.div.indices.pH.3.pdf",height = 11, width = 3)
ggsave("rarefied.min.alpha.div.indices.pH.3.png", height = 11, width = 3)

## By conductivity
rarefied.min.alpha.div.indices.cond.3 <- diversity.indices.rarefied.min.long.3  %>% filter(!variable == "rarefied.min.no.3") %>% mutate(variable = factor(variable, levels=c("rarefied.min.chao.3","rarefied.min.shannon.3","rarefied.min.simpson.3","rarefied.min.pielou.3"))) %>% drop_na() %>% ggplot(aes(x=reorder(factor(Cond.rarefied.min.3), Cond.rarefied.min.3), y=value, fill=variable)) + geom_point(size=0.3) + stat_summary(fun=mean, geom="point", color="red", size=0.4, shape=23, position=position_dodge(0.75)) +  theme(axis.text.x =element_text(size=6, face = "bold", angle = 45 , hjust = 1),axis.text.y =element_text(size=6, face = "bold"),axis.title.y=element_blank(),axis.title.x = element_blank()) +  theme(legend.position = "none") + facet_wrap(vars(variable), scales = "free_y",labeller = as_labeller(new_var_names_rarefied_min_3),  ncol=1) + theme(strip.text.y = element_text(angle = 0, face = "bold", size = 8)) + theme(plot.title = element_blank()) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + geom_smooth(formula = y ~ x, method = "lm", se = TRUE, color = "grey40", aes(group = variable)) + stat_regline_equation(label.x = 1, label.y =0.7 ,aes(group = variable)) + stat_cor(label.x = 1, label.y = 1, aes(group = variable)) 


ggsave("rarefied.min.alpha.div.indices.cond.3.pdf",height = 11, width = 3)
ggsave("rarefied.min.alpha.div.indices.cond.3.png", height = 11, width = 3)


rarefied.min.alpha.div.indices.carb.3 <- diversity.indices.rarefied.min.long.3  %>% filter(!variable == "rarefied.min.no.3") %>% mutate(variable = factor(variable, levels=c("rarefied.min.chao.3","rarefied.min.shannon.3","rarefied.min.simpson.3","rarefied.min.pielou.3"))) %>% drop_na() %>% ggplot(aes(x=reorder(factor(carb.rarefied.min.3), carb.rarefied.min.3), y=value, fill=variable)) + geom_point(size=0.3) + stat_summary(fun=mean, geom="point", color="red", size=0.4, shape=23, position=position_dodge(0.75)) +  theme(axis.text.x =element_text(size=6, face = "bold", angle = 45 , hjust = 1),axis.text.y =element_text(size=6, face = "bold"),axis.title.y=element_blank(),axis.title.x = element_blank()) +  theme(legend.position = "none") + facet_wrap(vars(variable), scales = "free_y",labeller = as_labeller(new_var_names_rarefied_min_3),  ncol=1) + theme(strip.text.y = element_text(angle = 0, face = "bold", size = 8)) + theme(plot.title = element_blank()) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + geom_smooth(formula = y ~ x, method = "lm", se = TRUE, color = "grey40", aes(group = variable)) + stat_regline_equation(label.x = 1, label.y =0.7 ,aes(group = variable)) + stat_cor(label.x = 1, label.y = 1000, aes(group = variable)) 

ggsave("rarefied.min.alpha.div.indices.carb.3.pdf",height = 11, width = 3)
ggsave("rarefied.min.alpha.div.indices.carb.3.png", height = 11, width = 3)

### ANOVA on alpha diversity indices Vs all variables ####

## ANOVA between alpha diversity indices and site

chao.site.aov <- aov(as.numeric(rarefied.min.chao.3) ~ as.factor(Site.rarefied.min.3), data = diversity.indices.rarefied.min.3)
summary(chao.site.aov)

shannon.site.aov <- aov(as.numeric(rarefied.min.shannon.3) ~ as.factor(Site.rarefied.min.3), data = diversity.indices.rarefied.min.3)
summary(shannon.site.aov)

simpson.site.aov <- aov(as.numeric(rarefied.min.simpson.3) ~ as.factor(Site.rarefied.min.3), data = diversity.indices.rarefied.min.3)
summary(simpson.site.aov)

pielou.site.aov <- aov(as.numeric(rarefied.min.pielou.3) ~ as.factor(Site.rarefied.min.3), data = diversity.indices.rarefied.min.3)
summary(pielou.site.aov)


## ANOVA between alpha diversity indices and latitude

chao.lat.aov <- aov(as.numeric(rarefied.min.chao.3) ~ as.numeric(Lat.rarefied.min.3), data = diversity.indices.rarefied.min.3)
summary(chao.lat.aov)

shannon.lat.aov <- aov(as.numeric(rarefied.min.shannon.3) ~ as.numeric(Lat.rarefied.min.3), data = diversity.indices.rarefied.min.3)
summary(shannon.lat.aov)

simpson.lat.aov <- aov(as.numeric(rarefied.min.simpson.3) ~ as.numeric(Lat.rarefied.min.3), data = diversity.indices.rarefied.min.3)
summary(simpson.lat.aov)

pielou.lat.aov <- aov(as.numeric(rarefied.min.pielou.3) ~ as.numeric(Lat.rarefied.min.3), data = diversity.indices.rarefied.min.3)
summary(pielou.lat.aov)


## ANOVA between alpha diversity indices and pH

chao.ph.aov <- aov(as.numeric(rarefied.min.chao.3) ~ as.numeric(pH.rarefied.min.3), data = diversity.indices.rarefied.min.3)
summary(chao.ph.aov)

shannon.ph.aov <- aov(as.numeric(rarefied.min.shannon.3) ~ as.numeric(pH.rarefied.min.3), data = diversity.indices.rarefied.min.3)
summary(shannon.ph.aov)

simpson.ph.aov <- aov(as.numeric(rarefied.min.simpson.3) ~ as.numeric(pH.rarefied.min.3), data = diversity.indices.rarefied.min.3)
summary(simpson.ph.aov)

pielou.ph.aov <- aov(as.numeric(rarefied.min.pielou.3) ~ as.numeric(pH.rarefied.min.3), data = diversity.indices.rarefied.min.3)
summary(pielou.ph.aov)


## ANOVA between alpha diversity indices and conductivity

chao.cond.aov <- aov(as.numeric(rarefied.min.chao.3) ~ as.numeric(Cond.rarefied.min.3), data = diversity.indices.rarefied.min.3)
summary(chao.cond.aov)

shannon.cond.aov <- aov(as.numeric(rarefied.min.shannon.3) ~ as.numeric(Cond.rarefied.min.3), data = diversity.indices.rarefied.min.3)
summary(shannon.cond.aov)

simpson.cond.aov <- aov(as.numeric(rarefied.min.simpson.3) ~ as.numeric(Cond.rarefied.min.3), data = diversity.indices.rarefied.min.3)
summary(simpson.cond.aov)

pielou.cond.aov <- aov(as.numeric(rarefied.min.pielou.3) ~ as.numeric(Cond.rarefied.min.3), data = diversity.indices.rarefied.min.3)
summary(pielou.cond.aov)


## ANOVA between alpha diversity indices and carbonate

chao.carb.aov <- aov(as.numeric(rarefied.min.chao.3) ~ as.numeric(carb.rarefied.min.3), data = diversity.indices.rarefied.min.3)
summary(chao.carb.aov)

shannon.carb.aov <- aov(as.numeric(rarefied.min.shannon.3) ~ as.numeric(carb.rarefied.min.3), data = diversity.indices.rarefied.min.3)
summary(shannon.carb.aov)

simpson.carb.aov <- aov(as.numeric(rarefied.min.simpson.3) ~ as.numeric(carb.rarefied.min.3), data = diversity.indices.rarefied.min.3)
summary(simpson.carb.aov)

pielou.carb.aov <- aov(as.numeric(rarefied.min.pielou.3) ~ as.numeric(carb.rarefied.min.3), data = diversity.indices.rarefied.min.3)
summary(pielou.carb.aov)

## paired T-tests of alpha diversity indices Vs Site only as these were significant ####

chao.site.bonf <- pairwise.t.test(diversity.indices.rarefied.min.3$rarefied.min.chao.3 ,as.factor(diversity.indices.rarefied.min.3$Site.rarefied.min.3) , p.adj = "bonf")
chao.site.bonf.tb <- chao.site.bonf$p.value

shannon.site.bonf <- pairwise.t.test(diversity.indices.rarefied.min.3$rarefied.min.shannon.3 ,as.factor(diversity.indices.rarefied.min.3$Site.rarefied.min.3) , p.adj = "bonf")
shannon.site.bonf.tb <- shannon.site.bonf$p.value

simpson.site.bonf <- pairwise.t.test(diversity.indices.rarefied.min.3$rarefied.min.simpson.3 ,as.factor(diversity.indices.rarefied.min.3$Site.rarefied.min.3) , p.adj = "bonf")
simpson.site.bonf.tb <- simpson.site.bonf$p.value

pielou.site.bonf <- pairwise.t.test(diversity.indices.rarefied.min.3$rarefied.min.pielou.3 ,as.factor(diversity.indices.rarefied.min.3$Site.rarefied.min.3) , p.adj = "bonf")
pielou.site.bonf.tb <- pielou.site.bonf$p.value


### Post-hoc Tukey's test of alpha diversity indices Vs Site only as these were significant  ####

chao.tukey <- (TukeyHSD(chao.aov, "Site.rarefied.min.3", group = TRUE))
chao.tukey.tb <- as.data.frame(chao.tukey$Site.rarefied.min.3)

shannon.tukey <- (TukeyHSD(shannon.aov, "Site.rarefied.min.3", group = TRUE))
shannon.tukey.tb <- as.data.frame(shannon.tukey$Site.rarefied.min.3)

simpson.tukey <- (TukeyHSD(simpson.aov, "Site.rarefied.min.3", group = TRUE))
simpson.tukey.tb <- as.data.frame(simpson.tukey$Site.rarefied.min.3)

pielou.tukey <- (TukeyHSD(pielou.aov, "Site.rarefied.min.3", group = TRUE))
pielou.tukey.tb <- as.data.frame(pielou.tukey$Site.rarefied.min.3)



#################################################################################

### Fig. S4. Taxonomic composition of ASVs clustered at Phylum and Class level (N = 395) ####

#################################################################################

rarefied.class.prop.exclude.1.3 <- tax_glom(rarefied.min.prop.exclude.1.3, taxrank="Class")

rarefied.class.prop.wide.exclude.1.3 <- as.matrix(as.data.frame(rarefied.class.prop.exclude.1.3@otu_table))

colnames(rarefied.class.prop.wide.exclude.1.3) <- as.character(tax_table(rarefied.class.prop.exclude.1.3)[, "Class"])

meta <- as.data.frame(metadata) %>% select(Sample_id,Site,Temp.adj...C.)

## arrange all ASV tables
rarefied.class.prop.wide.exclude.1.3.reordered <- as.data.frame(rarefied.class.prop.wide.exclude.1.3) %>% rownames_to_column(var = 'Sample_id') 

## merge meta and reordered tables
merge.rare.class.ex1.meta <- merge(meta, rarefied.class.prop.wide.exclude.1.3.reordered) %>% as_tibble()

## order sample_ids alphanumerically
merge.rare.class.ex1.meta <- merge.rare.class.ex1.meta [gtools::mixedorder(merge.rare.class.ex1.meta$Sample_id), ] %>% column_to_rownames(var = "Sample_id")

## now create reordered matrix from lowest to highest temperature
rarefied.class.prop.wide.exclude.1.3.reordered <- as.data.frame((merge.rare.class.ex1.meta[order(merge.rare.class.ex1.meta[,2],decreasing=FALSE),])) %>% select(-Site,-Temp.adj...C.)

#order columns (Class) in alphabetical order and Move Unassigneds to the last column
rarefied.class.prop.wide.exclude.1.3.reordered <- rarefied.class.prop.wide.exclude.1.3.reordered[,order(colnames(rarefied.class.prop.wide.exclude.1.3.reordered))] %>% select(-contains("Unassigned"), everything())

### order meta according to sample_ids from lowest to highest temp to split matrix in complexheatmap
meta <- meta[gtools::mixedorder(meta$Sample_id), ]  
meta <- as.data.frame((meta[order(meta[,3],decreasing=FALSE),]))

# Define the annotation color for columns and rows
## column annotations

rare_annotation_col_1.3 = data.frame(
  Phylum = as.factor(tax_table(rarefied.class.prop.exclude.1.3 )[, "Phylum"]) )

rownames(rare_annotation_col_1.3) = colnames(rarefied.class.prop.wide.exclude.1.3)

rare_annotation_col_1.3 <- rare_annotation_col_1.3 %>% rownames_to_column
rare_annotation_col_1.3 <- rare_annotation_col_1.3 %>% rename(Class = rowname)
rare_annotation_col_1.3 <- rare_annotation_col_1.3 [order(rare_annotation_col_1.3$Class),]

#move unassigned to the end in column annotations - option 1
rare_match_1.3 <- grep('^Unassigned', strsplit(rare_annotation_col_1.3$Class, '\n'), value = TRUE)
Class <- as.character(rare_annotation_col_1.3$Class) %in% rare_match_1.3
rare_annotation_col_1.3.reordered <- rbind(rare_annotation_col_1.3[!Class,], rare_annotation_col_1.3[Class,])

# get colors for phylum
rare.Phylum.fac.1 <-   data.frame(Phylum = as.factor(rare_annotation_col_1.3.reordered$Phylum))
rare_phylum_levels_1 <- length(levels(rare.Phylum.fac.1$Phylum))
rare_phylum_col_1 <- colorRampPalette(brewer.pal(12, "Paired"))(rare_phylum_levels_1)
names(rare_phylum_col_1) <- levels(rare.Phylum.fac.1$Phylum)

#colors as list
rare.col.phylum.1 = list(Phylum = rare_phylum_col_1)

rare.class.column.annotations.ex1 <- HeatmapAnnotation(
  Phylum = rare_annotation_col_1.3.reordered$Phylum,
  col = rare.col.phylum.1,annotation_name_gp = gpar(fontsize = 10, fontface = "bold"),
  annotation_legend_param = list(title_gp = gpar(fontsize = 9,fontface = "bold"),
                                 labels_gp = gpar(fontsize = 7.5)),show_annotation_name = FALSE
)


## row annotation
col.env = list(Site = c(PB_38 = "#5881F9", RB_40 = "#A7A3FE", HS_42 = "#816AA7", MK_42 = "#8CCBD3", RN_42 = "#0DD7FD", MK_44 = "#0D8CA8", HS_45 = "#2AB2FC" , RB_45 = "#5881F8" , RN_45 = "#4B3BA7", SE_45 = "#4B00FD", AP_46 = "#0DDDD0" , LA_46 = "#0DE36D", LN_46 = "#B2C9B2", PB_46 = "#99D584", SW_46 = "#9CD700", KJ_47 = "#848F22", PP_48 = "#009458", SE_49 = "#2A7F72", BA_50 = "#3B4D16", MK_50 = "#E2BD7D" , PT_50 = "#DFC500" , RB_50 = "#FEAB2E", HS_51 = "#967D63", LN_51 = "#F76A4B", RN_51 = "#C8682A", SW_51 = "#8B5500", HS_53 = "#EAB1DF", PP_53 = "#E19AFE", MK_55 = "#DE0DFD", PT_55 = "#A655ED", RB_55 = "#BB00B8", SW_55 = "#7E1C7C", US_57 = "#FD85CA", PP_58 = "#F9168D", RN_61 = "#FD0DD1", SW_61 = "#FD8A91", LN_63 = "#990D6A",PT_63 = "#CC003D", US_63 = "#F81626", LN_66 = "#950D3D"))

scalebr5 = colorRamp2(c(0,10,25,40,55), c("seashell4","yellow","yellowgreen","springgreen3","darkgreen"))

rarefied.class.prop.wide.exclude.1.3.avg <- merge.rare.class.ex1.meta %>%
  group_by(Site,Temp.adj...C.) %>% 
  summarise(across(everything(), mean),.groups = 'keep') %>% as_tibble()

## no more sample_ids so arrange Site as rownames and order from lowest to highest temp
rarefied.class.prop.wide.exclude.1.3.avg.reordered <- as.data.frame((rarefied.class.prop.wide.exclude.1.3.avg[order(rarefied.class.prop.wide.exclude.1.3.avg[,2],decreasing=FALSE),]))  %>% column_to_rownames(var = "Site") %>% select(-Temp.adj...C.)

#order columns (Class) in alphabetical order and Move Unassigneds to the last column
rarefied.class.prop.wide.exclude.1.3.avg.reordered <- rarefied.class.prop.wide.exclude.1.3.avg.reordered[,order(colnames(rarefied.class.prop.wide.exclude.1.3.avg.reordered))] %>% select(-contains("Unassigned"), everything())

## row annotation
avg.row.annotations.left.1 <-  rowAnnotation(Site = meta.avg$Site,annotation_name_rot = 45,annotation_name_gp = gpar(fontsize = 7, fontface = "bold"),annotation_legend_param = list(title_gp = gpar(fontsize = 7, fontface = "bold"),labels_gp = gpar(fontsize = 6)),col = col.env, show_legend = FALSE)

meta.avg <- class.prop.wide.exclude.1.3.avg %>% select(Site,Temp.adj...C.) 
meta.avg <- as.data.frame((meta.avg[order(meta.avg[,2],decreasing=FALSE),]))
rownames(meta.avg) <- meta.avg$Site

## execute lines 983-987 together
pdf(file="rarefied.class.exclude.1.3.avg.compheatmap.nosplit.pdf", width = 15, height = 11)
rarefied.class.exclude.1.3.avg.compheatmap.nosplit <- Heatmap((as.matrix(rarefied.class.prop.wide.exclude.1.3.avg.reordered)), 
                                                              name = "Relative abundance",heatmap_legend_param = list(title = "Relative abundance (%)", title_gp = gpar(fontsize = 10, fontface = 'bold'), labels_gp = gpar(fontsize = 9),ncol = 3, legend_height = unit(7, "cm"),legend_width = unit(4, "cm"), legend_direction = "horizontal"),
                                                              column_title = "Rarefied heatmap of average relative abundances at the Class level", column_title_gp = gpar(fontsize = 13, fontface = 'bold'), row_title = "Samples",row_title_gp = gpar(fontsize = 10,fontface = 'bold'),row_names_gp = gpar(fontsize = 9), column_names_gp = gpar(fontsize = 8.9, fontface = "italic"), col = scalebr5,cluster_rows = FALSE, cluster_columns = FALSE , left_annotation = avg.row.annotations.left.1, bottom_annotation = rare.class.column.annotations.ex1, row_names_side = "left", row_names_centered = TRUE,  column_names_rot = 45, use_raster = TRUE, raster_by_magick = TRUE,width = ncol(rarefied.class.prop.wide.exclude.1.3.avg.reordered)*unit(5, "mm"), height = nrow(rarefied.class.prop.wide.exclude.1.3.avg.reordered)*unit(5, "mm"))
rarefied.class.exclude.1.3.avg.compheatmap.nosplit <- draw(rarefied.class.exclude.1.3.avg.compheatmap.nosplit , heatmap_legend_side="bottom", annotation_legend_side="right",legend_grouping = "original")



####################################################################################3

### Fig. 3A Relative abundance of the 25 most abundant ASVs in each sample (N = 395) shown clustered at genus level (bar plots). ####
### Fig. S6 Relative abundance of the 10 most abundant ASVs in ecological groups for each sample (N = 395). ####

###################################################################################

### use microviz to calculate PCoA ordination with unweighted unifrac and rooted tree as microviz needs a rooted tree and does not internally root tree like phyloseq.
### microviz prefers counts instead of proportions or percentages

## Root tree with ape ####
##https://github.com/joey711/phyloseq/issues/597

## function by joey to pick outgroup with longest terminal branch - It returns the string of the tip-label for the new proposed outgroup
pick_new_outgroup <- function(tree.unrooted){
  require("magrittr")
  require("data.table")
  require("ape") # ape::Ntip
  # tablify parts of tree that we need.
  treeDT <- 
    cbind(
      data.table(tree.unrooted$edge),
      data.table(length = tree.unrooted$edge.length)
    )[1:Ntip(tree.unrooted)] %>% 
    cbind(data.table(id = tree.unrooted$tip.label))
  # Take the longest terminal branch as outgroup
  new.outgroup <- treeDT[which.max(length)]$id
  return(new.outgroup)
}


## now to root with function
new.outgroup.rare = pick_new_outgroup(phytree.rare)
rootedTree.rare = ape::root(phytree.rare, outgroup=new.outgroup.rare, resolve.root=TRUE)

is.rooted(rootedTree.rare)


### get rooted phyloseqs using both integers and proportions

com.rarefied.min.int.exclude.1.3.rooted <- merge_phyloseq(rarefied.min.int.exclude.1.3, rootedTree.rare)
  
com.rarefied.min.prop.exclude.1.3.rooted <- merge_phyloseq(rarefied.min.prop.exclude.1.3, rootedTree.rare)

saveRDS(com.rarefied.min.int.exclude.1.3.rooted, "com.rarefied.min.int.exclude.1.3.rooted.rds")
saveRDS(com.rarefied.min.prop.exclude.1.3.rooted, "com.rarefied.min.prop.exclude.1.3.rooted.rds")


### rooted phyloseqs with integers of photosynthetic and heterotrophic taxa ####

## Chloroflexi
rarefied.chloroflexi.int.exclude.1.3.rooted <- com.rarefied.min.int.exclude.1.3.rooted %>%
  subset_taxa(Phylum == "Chloroflexi" & Class  == "Chloroflexia")

## Cyanobacteria
rarefied.cyanobacteria.int.exclude.1.3.rooted <- com.rarefied.min.int.exclude.1.3.rooted %>%
  subset_taxa(Phylum == "Cyanobacteria" & Class  == "Cyanobacteriia")

##  Phototynthetic Proteoacteria and other photosynthetic groups
rarefied.proteobacteria.otherphoto.int.exclude.1.3.rooted <- com.rarefied.min.int.exclude.1.3.rooted %>%
  subset_taxa(Order %in% other.photo | Genus %in% proteobacteria.genus)

## all photosynthetic
rarefied.photosynthetic.int.exclude.1.3.rooted <- com.rarefied.min.int.exclude.1.3.rooted %>%
  subset_taxa( Class %in% photosynthetic.class | Order %in% other.photo | Genus %in% proteobacteria.genus )

## Chemolitotrophs
rarefied.chemolithotrophs.int.exclude.1.3.rooted <- com.rarefied.min.int.exclude.1.3.rooted %>%
  subset_taxa(Phylum %in% chemolithotroph.phylum  |  Class %in% chemolithotroph.class)

#### all  heterotrophs only
rarefied.only.heterotrophs.int.exclude.1.3.rooted <- com.rarefied.min.int.exclude.1.3.rooted %>%
  subset_taxa(!Phylum %in% chemolithotroph.phylum & !Class %in% chemolithotroph.class & !Class %in% photosynthetic.class & !Order %in% other.photo & !Genus %in% proteobacteria.genus )

saveRDS(rarefied.chloroflexi.int.exclude.1.3.rooted, "rarefied.chloroflexi.int.exclude.1.3.rooted.rds")
saveRDS(rarefied.cyanobacteria.int.exclude.1.3.rooted, "rarefied.cyanobacteria.int.exclude.1.3.rooted.rds")
saveRDS(rarefied.proteobacteria.otherphoto.int.exclude.1.3.rooted, "rarefied.proteobacteria.otherphoto.int.exclude.1.3.rooted.rds")
saveRDS(rarefied.photosynthetic.int.exclude.1.3.rooted, "rarefied.photosynthetic.int.exclude.1.3.rooted.rds")
saveRDS(rarefied.chemolithotrophs.int.exclude.1.3.rooted, "rarefied.chemolithotrophs.int.exclude.1.3.rooted.rds")
saveRDS(rarefied.only.heterotrophs.int.exclude.1.3.rooted, "rarefied.only.heterotrophs.int.exclude.1.3.rooted.rds")

### Use microviz to calculate PCoA ordination with unweighted unifrac and rooted tree as microviz needs a rooted tree and does not internally root tree like phyloseq.

rarefied.PCoA.unweighted.unifrac.microviz <- com.rarefied.min.int.exclude.1.3.rooted %>%
  phyloseq_validate(verbose = TRUE) %>%
  tax_transform("identity", rank = "unique") %>%
  dist_calc("unifrac") %>%
  ord_calc("PCoA") 

rarefied.PCoA.unweighted.unifrac.microviz.chloroflexi <- rarefied.chloroflexi.int.exclude.1.3.rooted %>%
  phyloseq_validate(verbose = TRUE) %>%
  tax_transform("identity", rank = "unique") %>%
  dist_calc("unifrac") %>%
  ord_calc("PCoA") 

rarefied.PCoA.unweighted.unifrac.microviz.cyanobacteria <- rarefied.cyanobacteria.int.exclude.1.3.rooted %>%
  phyloseq_validate(verbose = TRUE) %>%
  tax_transform("identity", rank = "unique") %>%
  dist_calc("unifrac") %>%
  ord_calc("PCoA") 


rarefied.PCoA.unweighted.unifrac.microviz.proteobacteria.otherphoto <- rarefied.proteobacteria.otherphoto.int.exclude.1.3.rooted %>%
  phyloseq_validate(verbose = TRUE) %>%
  tax_transform("identity", rank = "unique") %>%
  dist_calc("unifrac") %>%
  ord_calc("PCoA") 

rarefied.PCoA.unweighted.unifrac.microviz.photosynthetic <- rarefied.photosynthetic.int.exclude.1.3.rooted %>%
  phyloseq_validate(verbose = TRUE) %>%
  tax_transform("identity", rank = "unique") %>%
  dist_calc("unifrac") %>%
  ord_calc("PCoA") 

rarefied.PCoA.unweighted.unifrac.microviz.chemotrophs <- rarefied.chemolithotrophs.int.exclude.1.3.rooted %>%
  phyloseq_validate(verbose = TRUE) %>%
  tax_transform("identity", rank = "unique") %>%
  dist_calc("unifrac") %>%
  ord_calc("PCoA") 

rarefied.PCoA.unweighted.unifrac.microviz.only.heterotrophs <- rarefied.only.heterotrophs.int.exclude.1.3.rooted %>%
  phyloseq_validate(verbose = TRUE) %>%
  tax_transform("identity", rank = "unique") %>%
  dist_calc("unifrac") %>%
  ord_calc("PCoA") 


region.col = c('North.Thailand' = "slateblue", 'Central.Thailand' = "turquoise3", 'South.Thailand' = "gold", 'North.Malaysia' = "chocolate1", 'South.Malaysia' = "lawngreen", 'Singapore'="palevioletred3")

loc.col = c('PB' = "turquoise3",'RB' = "gold", 'HS' = "slateblue",'MK' = "slateblue",'RN' = "gold", 'SE'= "lawngreen",'AP' ="chocolate1" ,'LA' ="lawngreen", 'LN'="slateblue", 'SW' ="palevioletred3",  'KJ' ="lawngreen",'PP' ="slateblue",'BA' = "chocolate1",'PT'="slateblue", 'US' = "chocolate1")

rarefied.PCoA.unweighted.unifrac.microviz.irisplot.regionanno.genus.25 <- rarefied.PCoA.unweighted.unifrac.microviz %>% ord_plot_iris(tax_level = "Genus", n_taxa = 25, ord_plot = "none",anno_colour = "Region", anno_colour_style = list(size = 25)) + scale_colour_manual(values = region.col) + theme(legend.title = element_text(size = 15.5, face = "bold"), legend.text = element_text(size = 7),legend.position = "bottom", legend.direction = "horizontal",legend.box="vertical", legend.margin=margin()) + guides(colour = guide_legend(nrow=1,byrow=TRUE, override.aes = list(size = 7)), fill = guide_legend(nrow=4,byrow=TRUE,override.aes = list(size = 7)) )

pdf("rarefied.PCoA.unweighted.unifrac.microviz.irisplot.regionanno.genus.25.pdf", width = 15, height = 10)
rarefied.PCoA.unweighted.unifrac.microviz.irisplot.regionanno.genus.25
dev.off()

rarefied.PCoA.unweighted.unifrac.microviz.irisplot.regionanno.chloroflexi.genus <- rarefied.PCoA.unweighted.unifrac.microviz.chloroflexi %>% ord_plot_iris(tax_level = "Genus", n_taxa = 10, ord_plot = "none",anno_colour = "Region", anno_colour_style = list(size = 10)) + scale_colour_manual(values = region.col) + theme(legend.title = element_text(size = 15.5, face = "bold"), legend.text = element_text(size = 14),legend.position = "bottom", legend.direction = "horizontal",legend.box="vertical", legend.margin=margin()) + guides(colour = guide_legend(nrow=1,byrow=TRUE, override.aes = list(size = 7)), fill = guide_legend(nrow=2,byrow=TRUE,override.aes = list(size = 7)) )

pdf("rarefied.PCoA.unweighted.unifrac.microviz.irisplot.regionanno.chloroflexi.genus.pdf", width = 15, height = 10)
rarefied.PCoA.unweighted.unifrac.microviz.irisplot.regionanno.chloroflexi.genus
dev.off()

rarefied.PCoA.unweighted.unifrac.microviz.irisplot.regionanno.cyanobacteria.genus <- rarefied.PCoA.unweighted.unifrac.microviz.cyanobacteria %>% ord_plot_iris(tax_level = "Genus", n_taxa = 10, ord_plot = "none",anno_colour = "Region", anno_colour_style = list(size = 10)) + scale_colour_manual(values = region.col) + theme(legend.title = element_text(size = 15.5, face = "bold"), legend.text = element_text(size = 14),legend.position = "bottom", legend.direction = "horizontal",legend.box="vertical", legend.margin=margin()) + guides(colour = guide_legend(nrow=1,byrow=TRUE, override.aes = list(size = 7)), fill = guide_legend(nrow=4,byrow=TRUE,override.aes = list(size = 7)) )

pdf("rarefied.PCoA.unweighted.unifrac.microviz.irisplot.regionanno.cyanobacteria.genus.pdf", width = 15, height = 10)
rarefied.PCoA.unweighted.unifrac.microviz.irisplot.regionanno.cyanobacteria.genus
dev.off()

rarefied.PCoA.unweighted.unifrac.microviz.irisplot.regionanno.proteobacteria.otherphoto.genus <- rarefied.PCoA.unweighted.unifrac.microviz.proteobacteria.otherphoto %>% ord_plot_iris(tax_level = "Genus", n_taxa = 10, ord_plot = "none",anno_colour = "Region", anno_colour_style = list(size = 10)) + scale_colour_manual(values = region.col) + theme(legend.title = element_text(size = 15.5, face = "bold"), legend.text = element_text(size = 14),legend.position = "bottom", legend.direction = "horizontal",legend.box="vertical", legend.margin=margin()) + guides(colour = guide_legend(nrow=1,byrow=TRUE, override.aes = list(size = 7)), fill = guide_legend(nrow=2,byrow=TRUE,override.aes = list(size = 7)) )

pdf("rarefied.PCoA.unweighted.unifrac.microviz.irisplot.regionanno.proteobacteria.otherphoto.genus.pdf", width = 15, height = 10)
rarefied.PCoA.unweighted.unifrac.microviz.irisplot.regionanno.proteobacteria.otherphoto.genus
dev.off()

### Photosynthetic 
rarefied.PCoA.unweighted.unifrac.microviz.irisplot.regionanno.photosynthetic.genus <- rarefied.PCoA.unweighted.unifrac.microviz.photosynthetic %>% ord_plot_iris(tax_level = "Genus", n_taxa = 10, ord_plot = "none",anno_colour = "Region", anno_colour_style = list(size = 10)) + scale_colour_manual(values = region.col) + theme(legend.title = element_text(size = 15.5, face = "bold"), legend.text = element_text(size = 14),legend.position = "bottom", legend.direction = "horizontal",legend.box="vertical", legend.margin=margin()) + guides(colour = guide_legend(nrow=1,byrow=TRUE, override.aes = list(size = 7)), fill = guide_legend(nrow=4,byrow=TRUE,override.aes = list(size = 7)) )

pdf("rarefied.PCoA.unweighted.unifrac.microviz.irisplot.regionanno.photosynthetic.genus.pdf", width = 15, height = 10)
rarefied.PCoA.unweighted.unifrac.microviz.irisplot.regionanno.photosynthetic.genus
dev.off()

### Chemolithotrophs #####

rarefied.PCoA.unweighted.unifrac.microviz.irisplot.regionanno.chemolithotrophs.genus <- rarefied.PCoA.unweighted.unifrac.microviz.chemotrophs %>% ord_plot_iris(tax_level = "Genus", n_taxa = 10, ord_plot = "none",anno_colour = "Region", anno_colour_style = list(size = 10)) + scale_colour_manual(values = region.col) + theme(legend.title = element_text(size = 15.5, face = "bold"), legend.text = element_text(size = 14),legend.position = "bottom", legend.direction = "horizontal",legend.box="vertical", legend.margin=margin()) + guides(colour = guide_legend(nrow=1,byrow=TRUE, override.aes = list(size = 7)), fill = guide_legend(nrow=4,byrow=TRUE,override.aes = list(size = 7)) )

pdf("rarefied.PCoA.unweighted.unifrac.microviz.irisplot.regionanno.chemolithotrophs.genus.pdf", width = 15, height = 10)
rarefied.PCoA.unweighted.unifrac.microviz.irisplot.regionanno.chemolithotrophs.genus
dev.off()

### Only Heterotrophs #####
rarefied.PCoA.unweighted.unifrac.microviz.irisplot.regionanno.only.heterotrophs.genus <- rarefied.PCoA.unweighted.unifrac.microviz.only.heterotrophs %>% ord_plot_iris(tax_level = "Genus", n_taxa = 10, ord_plot = "none",anno_colour = "Region", anno_colour_style = list(size = 10)) + scale_colour_manual(values = region.col) + theme(legend.title = element_text(size = 15.5, face = "bold"), legend.text = element_text(size = 14),legend.position = "bottom", legend.direction = "horizontal",legend.box="vertical", legend.margin=margin()) + guides(colour = guide_legend(nrow=1,byrow=TRUE, override.aes = list(size = 7)), fill = guide_legend(nrow=4,byrow=TRUE,override.aes = list(size = 7)) )

pdf("rarefied.PCoA.unweighted.unifrac.microviz.irisplot.regionanno.only.heterotrophs.genus.pdf", width = 15, height = 10)
rarefied.PCoA.unweighted.unifrac.microviz.irisplot.regionanno.only.heterotrophs.genus
dev.off()


###################################################################################

### Fig. 3B. Composition of the core microbiome for biofilms in each biogeographic region (N = 395) ####

###################################################################################

### CORE MEMBERS FOR 6 DIFFERENT REGIONS #####

## cluster 1 :North.Thailand <- c('PT','PP','LN','HS','MK') #####

cluster1.loc <- c("PT","PP","LN","HS","MK")

cluster1.PCoA.unweighted.unifrac <- com.rarefied.min.int.exclude.1.3.rooted %>%
  subset_samples(Location.Code %in% cluster1.loc)

cluster1.PCoA.unweighted.unifrac <- 
  prune_taxa(taxa_sums(cluster1.PCoA.unweighted.unifrac) > 0, cluster1.PCoA.unweighted.unifrac)

any(taxa_sums(cluster1.PCoA.unweighted.unifrac) == 0)

cluster1.PCoA.unweighted.unifrac.propTo1 <- transform(cluster1.PCoA.unweighted.unifrac, 'compositional')

core.rarefied.exc1.cluster1.genus.95 <- aggregate_rare(cluster1.PCoA.unweighted.unifrac.propTo1, "Genus", detection = 0, prevalence = 95/100)
taxa_names(core.rarefied.exc1.cluster1.genus.95)

## Cluster 2 : North.Malaysia <- c('AP','BA','US') ######

cluster2.loc <- c("AP","BA","US")

cluster2.PCoA.unweighted.unifrac <- com.rarefied.min.int.exclude.1.3.rooted %>%
  subset_samples( Location.Code %in% cluster2.loc)

cluster2.PCoA.unweighted.unifrac.propTo1 <- transform(cluster2.PCoA.unweighted.unifrac, 'compositional')

core.rarefied.exc1.cluster2.genus.95 <- aggregate_rare(cluster2.PCoA.unweighted.unifrac.propTo1, "Genus", detection = 0, prevalence = 95/100)
taxa_names(core.rarefied.exc1.cluster2.genus.95)


## Cluster 3 : South.Malaysia <- c('KJ','SE','LA') ######

cluster3.loc <- c("SE","LA","KJ") 

cluster3.PCoA.unweighted.unifrac <- com.rarefied.min.int.exclude.1.3.rooted %>%
  subset_samples(Location.Code %in% cluster3.loc)

cluster3.PCoA.unweighted.unifrac.propTo1 <- transform(cluster3.PCoA.unweighted.unifrac, 'compositional')

core.rarefied.exc1.cluster3.genus.95 <- aggregate_rare(cluster3.PCoA.unweighted.unifrac.propTo1, "Genus", detection = 0, prevalence = 95/100)
taxa_names(core.rarefied.exc1.cluster3.genus.95)

## Cluster 4 : All Sembawang : Singapore <- 'SW' #####

cluster4.PCoA.unweighted.unifrac <- com.rarefied.min.int.exclude.1.3.rooted %>%
  subset_samples(Location.Code == "SW")

cluster4.PCoA.unweighted.unifrac.propTo1 <- transform(cluster4.PCoA.unweighted.unifrac, 'compositional')

core.rarefied.exc1.cluster4.genus.95 <- aggregate_rare(cluster4.PCoA.unweighted.unifrac.propTo1, "Genus", detection = 0, prevalence = 95/100)
taxa_names(core.rarefied.exc1.cluster4.genus.95)

## Cluster 5 : South.Thailand <- c('RB','RN') #####

cluster5.loc <- c("RN","RB")

cluster5.PCoA.unweighted.unifrac <- com.rarefied.min.int.exclude.1.3.rooted %>%
  subset_samples(Location.Code %in% cluster5.loc)

cluster5.PCoA.unweighted.unifrac.propTo1 <- transform(cluster5.PCoA.unweighted.unifrac, 'compositional')

core.rarefied.exc1.cluster5.genus.95 <- aggregate_rare(cluster5.PCoA.unweighted.unifrac.propTo1, "Genus", detection = 0, prevalence = 95/100)
taxa_names(core.rarefied.exc1.cluster5.genus.95)

## Cluster 6 : Central.Thailand <- 'PB' ######

cluster6.PCoA.unweighted.unifrac <- com.rarefied.min.int.exclude.1.3.rooted %>%
  subset_samples(Location.Code == "PB" )

cluster6.PCoA.unweighted.unifrac.propTo1 <- transform(cluster6.PCoA.unweighted.unifrac, 'compositional')

core.rarefied.exc1.cluster6.genus.95 <- aggregate_rare(cluster6.PCoA.unweighted.unifrac.propTo1, "Genus", detection = 0, prevalence = 95/100)
taxa_names(core.rarefied.exc1.cluster6.genus.95)

### save cluster phyloseqs
saveRDS(cluster1.PCoA.unweighted.unifrac, "cluster1.PCoA.unweighted.unifrac.rds")
saveRDS(cluster2.PCoA.unweighted.unifrac, "cluster2.PCoA.unweighted.unifrac.rds")
saveRDS(cluster3.PCoA.unweighted.unifrac, "cluster3.PCoA.unweighted.unifrac.rds")
saveRDS(cluster4.PCoA.unweighted.unifrac, "cluster4.PCoA.unweighted.unifrac.rds")
saveRDS(cluster5.PCoA.unweighted.unifrac, "cluster5.PCoA.unweighted.unifrac.rds")
saveRDS(cluster6.PCoA.unweighted.unifrac, "cluster6.PCoA.unweighted.unifrac.rds")


### get core genera for entire dataset (not by each cluster)

com.rarefied.min.propTo1.exclude.1.3 <- transform(com.rarefied.min.int.exclude.1.3.rooted, 'compositional')

## sanity check for compositional data
is_compositional(com.rarefied.min.propTo1.exclude.1.3)

rare.genus.exc1.95 <- aggregate_rare(com.rarefied.min.propTo1.exclude.1.3, "Genus", detection = 0, prevalence = 95/100)
taxa_names(rare.genus.exc1.95)

## entire dataset core to subset
entire.dataset.core.95.genus.tosub <- c((taxa_names(rare.genus.exc1.95))[ !(c(taxa_names(rare.genus.exc1.95))) == 'Other'])
entire.dataset.core.95.genus.tosub  <- sort(unique(entire.dataset.core.95.genus.tosub ))
entire.dataset.core.95.genus.tosub 

### get raw abundances for all genera (core and non-core genera) from entire dataset
entiredata <- com.rarefied.min.int.exclude.1.3.rooted
entiredata.genera <- tax_glom(entiredata, taxrank = "Genus")
entiredata.genera.perc <- transform_sample_counts(entiredata.genera, function(OTU) (OTU/sum(OTU))*100)
entiredata.genera.melt <- psmelt(entiredata.genera.perc)

#change to character for easy-adjusted level
entiredata.genera.melt$Genus <- as.character(entiredata.genera.melt$Genus)

# Group according to faceting variable to get mean and median for further filtering
entiredata.genera.melt <- entiredata.genera.melt %>%
  group_by(Genus) %>%
  mutate(median=median(Abundance), mean=mean(Abundance))


## get mean abundances of genus 
entiredata.genera.melt.2 <- entiredata.genera.melt %>%
  group_by(Genus) %>%
  summarise(Abundance = mean(Abundance)) 

# to label core and non-core genus for entire dataset ###

entiredata.genera.melt.3 <- entiredata.genera.melt.2

entiredata.genera.melt.3$Genus[!(entiredata.genera.melt.3$Genus %in% entire.dataset.core.95.genus.tosub)] <- "non-core"

entiredata.genera.melt.3$Genus[(entiredata.genera.melt.3$Genus %in% entire.dataset.core.95.genus.tosub)] <- "core"

entiredata.genera.melt.4 <- entiredata.genera.melt.3 %>%
  group_by(Genus) %>%
  summarise(Abundance = sum(Abundance))

entiredata.genera.melt.5 <- as.data.frame(entiredata.genera.melt.4)
rownames(entiredata.genera.melt.5) <- entiredata.genera.melt.5$Genus

## colors to use for core and non-core : 

display.brewer.pal(n = 8, name = 'Paired')
brewer.pal(n = 16, name = "Paired")

colors <- c('#A6CEE3', '#1F78B4', '#B2DF8A', '#33A02C', '#FB9A99','#E31A1C','#FDBF6F','#FF7F00')

core.noncore <- c('rosybrown', 'lightgrey')


### Plot core genera as pie charts - entire dataset ####

## labels + percent
entiredata.genera.pie.1 <- plot_ly(entiredata.genera.melt.5, labels = ~Genus, values = ~Abundance, type = 'pie',marker = list(colors = core.noncore,line = list(color = '#FFFFFF', width = 2)), textposition = 'inside', textinfo = 'label+percent', insidetextfont = list(color = 'grey32'),hoverinfo = 'text') %>% layout(title = 'Entire dataset core and non-core genera',
                                                                                                                                                                                                                                                                                                                            xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),                                                                                                                                                                                                                                                                                                                          yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))
## percent
entiredata.genera.pie.2 <- plot_ly(entiredata.genera.melt.5, labels = ~Genus, values = ~Abundance, type = 'pie',marker = list(colors = core.noncore,line = list(color = '#FFFFFF', width = 2)), textposition = 'inside', textinfo = 'percent', insidetextfont = list(color = 'grey32'),hoverinfo = 'text') %>% layout(title = 'Entire dataset core and non-core genera',
                                                                                                                                                                                                                                                                                                                      xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),                                                                                                                                                                                                                                                                                                                   yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))
## no labels
entiredata.genera.pie.3 <- plot_ly(entiredata.genera.melt.5, labels = ~Genus, values = ~Abundance, type = 'pie',marker = list(colors = core.noncore, line = list(color = '#FFFFFF', width = 2)), textposition = 'inside', textinfo = 'none', insidetextfont = list(color = 'grey32'),hoverinfo = 'text') %>% layout(title = 'Entire dataset core and non-core genera',
                                                                                                                                                                                                                                                                                                                    xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),                                                                                                                                                                                                                                                                                                              yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))
### save plots as images or screenshots
entiredata.genera.pie.1
entiredata.genera.pie.2
entiredata.genera.pie.3


## colors
colors <- c('#A6CEE3', '#1F78B4', '#B2DF8A', '#33A02C', '#FB9A99','#E31A1C','#FDBF6F','#FF7F00', "#68A6CD", "#2A7FB7", "#569EA4", "#99CD91","#8CCC6E", "#52AF43","#5C9E42", "#B89B74", "#F88A89", "#ED4F50", "#E4201F","#F06C45","#FBB86B","#FDA440", "#FE870D", "#ED8F47", "#D5A7A9"  ,"#B294C7","#865FAB", "#825D99","#C7B699", "#F8F18F", "#D4A55B", "#B15928")

brewer.pal(n = 14, name = "Dark2")

colourCountGenus = 100
getPalette = colorRampPalette(brewer.pal(14, "Paired"))
#assign order names to colors to ensure subsetted order have same colors as full plot
value.genus.100 = getPalette(colourCountGenus)
print(value.genus.100)


cluster1.core <- c('#33A02C','#FB9A99','#A6CEE3','#B2DF8A','#1F78B4','#7570B3',"#CDBE99")

cluster2.core <- c('#33A02C','#FB9A99','#A6CEE3', '#B2DF8A','cyan',"#CDBE99")

cluster3.core <- c('#E31A1C','#990066',"#1B9E77", '#FF7F00','cyan',"#CDBE99")

cluster4.core <- c('#8B6508','#B2DF8A', '#D95F02','#1F78B4',"#CDBE99")

cluster5.core <- c("#FFFF99",'#B4CDCD', "#B2DF8A",'#CAB206',"#CDBE99")

cluster6.core <- c('#551A8B', '#FB9A99','#E31A1C','#990066',"#CDBE99")


### Plot core genera as pie charts - by region ####

North.Thailand.95.core <- c(taxa_names(core.rarefied.exc1.cluster1.genus.95))[ !(c(taxa_names(core.rarefied.exc1.cluster1.genus.95))) == 'Other']

North.Malaysia.95.core <-  c(taxa_names(core.rarefied.exc1.cluster2.genus.95))[ !(c(taxa_names(core.rarefied.exc1.cluster2.genus.95))) == 'Other']

South.Malaysia.95.core  <- c(taxa_names(core.rarefied.exc1.cluster3.genus.95))[ !(c(taxa_names(core.rarefied.exc1.cluster3.genus.95))) == 'Other']

Singapore.95.core  <- c(taxa_names(core.rarefied.exc1.cluster4.genus.95))[ !(c(taxa_names(core.rarefied.exc1.cluster4.genus.95))) == 'Other']

South.Thailand.95.core  <- c(taxa_names(core.rarefied.exc1.cluster5.genus.95))[ !(c(taxa_names(core.rarefied.exc1.cluster5.genus.95))) == 'Other']

Central.Thailand.95.core  <- c(taxa_names(core.rarefied.exc1.cluster6.genus.95))[ !(c(taxa_names(core.rarefied.exc1.cluster6.genus.95))) == 'Other']

## Cluster 1 ####

cluster1 <- 
  prune_taxa(taxa_sums(cluster1.PCoA.unweighted.unifrac) > 0, cluster1.PCoA.unweighted.unifrac)

cluster1.genera <- tax_glom(cluster1, taxrank = "Genus")


cluster1.genera.perc <- transform_sample_counts(cluster1.genera, function(OTU) (OTU/sum(OTU))*100)

cluster1.genera.melt <- psmelt(cluster1.genera.perc)

#change to character for easy-adjusted level
cluster1.genera.melt$Genus <- as.character(cluster1.genera.melt$Genus)


# Group according to faceting variable to get mean and median for further filtering
cluster1.genera.melt <- cluster1.genera.melt %>%
  group_by(Genus) %>%
  mutate(median=median(Abundance), mean=mean(Abundance))


## get mean abundances of genus identified at each cluster
cluster1.genera.melt.2 <- cluster1.genera.melt %>%
  group_by(Genus) %>%
  summarise(Abundance = mean(Abundance)) 

# to label core and non-core genus ###

cluster1.genera.melt.3 <- cluster1.genera.melt.2
cluster1.genera.melt.3$Genus[!(cluster1.genera.melt.3$Genus %in% North.Thailand.95.core)] <- "non-core"
cluster1.genera.melt.3$Genus[(cluster1.genera.melt.3$Genus %in% North.Thailand.95.core)] <- "core"

cluster1.genera.melt.5 <- cluster1.genera.melt.3 %>%
  group_by(Genus) %>%
  summarise(Abundance = sum(Abundance))

cluster1.genera.melt.5 <- as.data.frame(cluster1.genera.melt.5)
rownames(cluster1.genera.melt.5) <- cluster1.genera.melt.5$Genus


## labels + percent
cluster1.genera.pie.1 <- plot_ly(cluster1.genera.melt.5, labels = ~Genus, values = ~Abundance, type = 'pie',marker = list(colors = core.noncore,line = list(color = '#FFFFFF', width = 2)), textposition = 'inside', textinfo = 'label+percent', insidetextfont = list(color = '#FFFFFF'),hoverinfo = 'text') %>% layout(title = 'North Thailand (cluster1) core and non-core genera',
                                                                                                                                                                                                                                                                                                                         xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                                                                                                                                                                                                                                                                                                                         yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))
## percent
cluster1.genera.pie.2 <- plot_ly(cluster1.genera.melt.5, labels = ~Genus, values = ~Abundance, type = 'pie',marker = list(colors = core.noncore,line = list(color = '#FFFFFF', width = 2)), textposition = 'inside', textinfo = 'lpercent', insidetextfont = list(color = '#FFFFFF'),hoverinfo = 'text') %>% layout(title = 'North Thailand (cluster1) core and non-core genera',
                                                                                                                                                                                                                                                                                                                    xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                                                                                                                                                                                                                                                                                                                    yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))

## no labels
cluster1.genera.pie.3 <- plot_ly(cluster1.genera.melt.5, labels = ~Genus, values = ~Abundance, type = 'pie',marker = list(colors = core.noncore,line = list(color = '#FFFFFF', width = 2)), textinfo = 'none', insidetextfont = list(color = '#FFFFFF'),hoverinfo = 'text') %>% layout(title = 'North Thailand (cluster1) core and non-core genera',
                                                                                                                                                                                                                                                                                       xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                                                                                                                                                                                                                                                                                       yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))

## core genera pie chart ####
cluster1.core.genera <- cluster1.genera %>% subset_taxa(Genus %in% North.Thailand.95.core)

#sanity check
test <- as.data.frame(cluster1.core.genera@tax_table)
test$Genus
North.Thailand.95.core
identical(North.Thailand.95.core,sort(test$Genus))

cluster1.core.genera.perc <- transform_sample_counts(cluster1.core.genera, function(OTU) (OTU/sum(OTU))*100)

sample_sums(cluster1.core.genera.perc)

cluster1.core.genera.perc.melt <- psmelt(cluster1.core.genera.perc)

#change to character for easy-adjusted level
cluster1.core.genera.perc.melt$Genus <- as.character(cluster1.core.genera.perc.melt$Genus)

# Group according to faceting variable to get mean and median for further filtering
cluster1.core.genera.perc.melt <- cluster1.core.genera.perc.melt %>%
  group_by(Genus) %>%
  mutate(median=median(Abundance), mean=mean(Abundance))


### ## get mean abundances of genus identified at each cluster
cluster1.core.genera.perc.melt.2 <- cluster1.core.genera.perc.melt %>%
  group_by(Genus) %>%
  summarise(Abundance = mean(Abundance)) 


## Core genera with relative abundances < 5 are labelled as 'Others' in the final figures

cluster1.core.genera.perc.melt.3 <- cluster1.core.genera.perc.melt.2


keep.5 <- unique(cluster1.core.genera.perc.melt.3$Genus[cluster1.core.genera.perc.melt.3$Abundance >= 5]) 
cluster1.core.genera.perc.melt.3$Genus[!(cluster1.core.genera.perc.melt.3$Genus %in% keep.5)] <- "< 5%"

## to check if to include <5% as others in pie-chart for core-genera, but others make up 43%
cluster1.core.genera.perc.melt.4 <- cluster1.core.genera.perc.melt.3 %>%
  group_by(Genus) %>%
  summarise(Abundance = sum(Abundance))

# Add label position :  does not apply correctly on plot
cluster1.core.genera.perc.melt.5 <- cluster1.core.genera.perc.melt.4 %>%
  arrange(desc(Genus))

cluster1.core.genera.perc.melt.5 <- as.data.frame(cluster1.core.genera.perc.melt.5)
rownames(cluster1.core.genera.perc.melt.5) <- cluster1.core.genera.perc.melt.5$Genus

##label +  percentages
cluster1.core.genera.pie.1 <- plot_ly(cluster1.core.genera.perc.melt.5, labels = ~Genus, values = ~Abundance, type = 'pie', marker = list(colors = cluster1.core, line = list(color = '#FFFFFF', width = 2)), textposition = 'inside', textinfo = 'label+percent', insidetextfont = list(color = 'grey32'),hoverinfo = 'text') %>% layout(title = 'North Thailand (cluster1) core genera',
                                                                                                                                                                                                                                                                                                                                          xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                                                                                                                                                                                                                                                                                                                                          yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))

## only percentages
cluster1.core.genera.pie.2 <- plot_ly(cluster1.core.genera.perc.melt.5, labels = ~Genus, values = ~Abundance, type = 'pie', marker = list(colors = cluster1.core,line = list(color = '#FFFFFF', width = 2)), textposition = 'inside', textinfo = 'percent', insidetextfont = list(color = 'grey32')) %>% layout(title = 'North Thailand (cluster1) core genera',
                                                                                                                                                                                                                                                                                                                xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                                                                                                                                                                                                                                                                                                                yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))

##no labels
cluster1.core.genera.pie.3 <- plot_ly(cluster1.core.genera.perc.melt.5, values = ~Abundance, type = 'pie', marker = list(colors = cluster1.core, line = list(color = 'white', width = 2)),textinfo = 'none') %>% layout(title = 'North Thailand (cluster1) core genera', 
                                                                                                                                                                                                                        xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                                                                                                                                                                                                                        yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))

cluster1.core.genera.pie
cluster1.genera.pie


## Cluster 2 ####

cluster2 <- 
  prune_taxa(taxa_sums(cluster2.PCoA.unweighted.unifrac) > 0, cluster2.PCoA.unweighted.unifrac)

cluster2.genera <- tax_glom(cluster2, taxrank = "Genus")

test <- as.data.frame(cluster2.genera@tax_table)

cluster2.genera.perc <- transform_sample_counts(cluster2.genera, function(OTU) (OTU/sum(OTU))*100)

cluster2.genera.melt <- psmelt(cluster2.genera.perc)

#change to character for easy-adjusted level
cluster2.genera.melt$Genus <- as.character(cluster2.genera.melt$Genus)


# Group according to faceting variable to get mean and median for further filtering
cluster2.genera.melt <- cluster2.genera.melt %>%
  group_by(Genus) %>%
  mutate(median=median(Abundance), mean=mean(Abundance))


## get mean abundances of genus identified at each cluster
cluster2.genera.melt.2 <- cluster2.genera.melt %>%
  group_by(Genus) %>%
  summarise(Abundance = mean(Abundance)) 

# to label core and non-core genus ###
cluster2.genera.melt.3 <- cluster2.genera.melt.2

cluster2.genera.melt.3$Genus[!(cluster2.genera.melt.3$Genus %in% North.Malaysia.95.core)] <- "non-core"

cluster2.genera.melt.3$Genus[(cluster2.genera.melt.3$Genus %in% North.Malaysia.95.core)] <- "core"

cluster2.genera.melt.5 <- cluster2.genera.melt.3 %>%
  group_by(Genus) %>%
  summarise(Abundance = sum(Abundance))

cluster2.genera.melt.5 <- as.data.frame(cluster2.genera.melt.5)
rownames(cluster2.genera.melt.5) <- cluster2.genera.melt.5$Genus

## labels + percent
cluster2.genera.pie.1 <- plot_ly(cluster2.genera.melt.5, labels = ~Genus, values = ~Abundance, type = 'pie',marker = list(colors = core.noncore,line = list(color = '#FFFFFF', width = 2)), textposition = 'inside', textinfo = 'label+percent', insidetextfont = list(color = 'grey32'),hoverinfo = 'text') %>% layout(title = 'North Malayasia (cluster2) core and non-core genera',
                                                                                                                                                                                                                                                                                                                        xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                                                                                                                                                                                                                                                                                                                        yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))

## percent
cluster2.genera.pie.2 <- plot_ly(cluster2.genera.melt.5, labels = ~Genus, values = ~Abundance, type = 'pie',marker = list(colors = core.noncore,line = list(color = '#FFFFFF', width = 2)), textposition = 'inside', textinfo = 'percent', insidetextfont = list(color = 'grey32'),hoverinfo = 'text') %>% layout(title = 'North Malayasia (cluster2) core and non-core genera',
                                                                                                                                                                                                                                                                                                                  xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                                                                                                                                                                                                                                                                                                                  yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))


## no labels
cluster2.genera.pie.3 <- plot_ly(cluster2.genera.melt.5, labels = ~Genus, values = ~Abundance, type = 'pie',marker = list(colors = core.noncore,line = list(color = '#FFFFFF', width = 2)),  textinfo = 'none', insidetextfont = list(color = '#FFFFFF'),hoverinfo = 'text') %>% layout(title = 'North Malayasia (cluster2) core and non-core genera',
                                                                                                                                                                                                                                                                                        xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                                                                                                                                                                                                                                                                                        yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))


## core genera pie chart ####
cluster2.core.genera <- cluster2.genera %>% subset_taxa(Genus %in% North.Malaysia.95.core)

#sanity check
test <- as.data.frame(cluster2.core.genera@tax_table)
test$Genus
North.Malaysia.95.core
identical(North.Malaysia.95.core,sort(test$Genus))

cluster2.core.genera.perc <- transform_sample_counts(cluster2.core.genera, function(OTU) (OTU/sum(OTU))*100)

test <- as.data.frame(cluster2.core.genera.perc@tax_table)

sample_sums(cluster2.core.genera.perc)

cluster2.core.genera.perc.melt <- psmelt(cluster2.core.genera.perc)

#change to character for easy-adjusted level
cluster2.core.genera.perc.melt$Genus <- as.character(cluster2.core.genera.perc.melt$Genus)

# Group according to faceting variable to get mean and median for further filtering
cluster2.core.genera.perc.melt <- cluster2.core.genera.perc.melt %>%
  group_by(Genus) %>%
  mutate(median=median(Abundance), mean=mean(Abundance))


### ## get mean abundances of genus identified at each cluster
cluster2.core.genera.perc.melt.2 <- cluster2.core.genera.perc.melt %>%
  group_by(Genus) %>%
  summarise(Abundance = mean(Abundance)) 


cluster2.core.genera.perc.melt.3 <- cluster2.core.genera.perc.melt.2

keep.5 <- unique(cluster2.core.genera.perc.melt.3$Genus[cluster2.core.genera.perc.melt.3$Abundance >= 5]) 
cluster2.core.genera.perc.melt.3$Genus[!(cluster2.core.genera.perc.melt.3$Genus %in% keep.5)] <- "< 5%"

## to check if to include <5% as others in pie-chart for core-genera, but others make up 43%
cluster2.core.genera.perc.melt.4 <- cluster2.core.genera.perc.melt.3 %>%
  group_by(Genus) %>%
  summarise(Abundance = sum(Abundance))


# Add label position :  does not apply correctly on plot
cluster2.core.genera.perc.melt.5 <- cluster2.core.genera.perc.melt.4 %>%
  arrange(desc(Genus)) 

cluster2.core.genera.perc.melt.5 <- as.data.frame(cluster2.core.genera.perc.melt.5)
rownames(cluster2.core.genera.perc.melt.5) <- cluster2.core.genera.perc.melt.5$Genus


## with labels + percent
cluster2.core.genera.pie.1 <- plot_ly(cluster2.core.genera.perc.melt.5, labels = ~Genus, values = ~Abundance, type = 'pie', marker = list(colors = cluster2.core,line = list(color = '#FFFFFF', width = 2)), textposition = 'inside', textinfo = 'label+percent', insidetextfont = list(color = 'grey32'),hoverinfo = 'text') %>% layout(title = 'North Malaysia (cluster2) core genera',
                                                                                                                                                                                                                                                                                                                                         xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                                                                                                                                                                                                                                                                                                                                         yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))

## with percent
cluster2.core.genera.pie.2 <- plot_ly(cluster2.core.genera.perc.melt.5, labels = ~Genus, values = ~Abundance, type = 'pie', marker = list(colors = cluster2.core,line = list(color = '#FFFFFF', width = 2)), textposition = 'inside', textinfo = 'percent', insidetextfont = list(color = 'grey32')) %>% layout(title = 'North Malaysia (cluster2) core genera',
                                                                                                                                                                                                                                                                                                                xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                                                                                                                                                                                                                                                                                                                yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))


## no labels
cluster2.core.genera.pie.3 <- plot_ly(cluster2.core.genera.perc.melt.5, labels = ~Genus, values = ~Abundance, type = 'pie', marker = list(colors = cluster2.core,line = list(color = '#FFFFFF', width = 2)), textinfo = 'none', insidetextfont = list(color = '#FFFFFF')) %>% layout(title = 'North Malaysia (cluster2) core genera',
                                                                                                                                                                                                                                                                                     xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                                                                                                                                                                                                                                                                                     yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))

cluster2.core.genera.pie
cluster2.genera.pie

## Cluster 3 ####

cluster3 <- 
  prune_taxa(taxa_sums(cluster3.PCoA.unweighted.unifrac) > 0, cluster3.PCoA.unweighted.unifrac)

cluster3.genera <- tax_glom(cluster3, taxrank = "Genus")

test <- as.data.frame(cluster3.genera@tax_table)

cluster3.genera.perc <- transform_sample_counts(cluster3.genera, function(OTU) (OTU/sum(OTU))*100)

cluster3.genera.melt <- psmelt(cluster3.genera.perc)

#change to character for easy-adjusted level
cluster3.genera.melt$Genus <- as.character(cluster3.genera.melt$Genus)


# Group according to faceting variable to get mean and median for further filtering
cluster3.genera.melt <- cluster3.genera.melt %>%
  group_by(Genus) %>%
  mutate(median=median(Abundance), mean=mean(Abundance))


## get mean abundances of genus identified at each cluster
cluster3.genera.melt.2 <- cluster3.genera.melt %>%
  group_by(Genus) %>%
  summarise(Abundance = mean(Abundance)) 

# to label core and non-core genus ###
cluster3.genera.melt.3 <- cluster3.genera.melt.2

cluster3.genera.melt.3$Genus[!(cluster3.genera.melt.3$Genus %in% South.Malaysia.95.core)] <- "non-core"

cluster3.genera.melt.3$Genus[(cluster3.genera.melt.3$Genus %in% South.Malaysia.95.core)] <- "core"

cluster3.genera.melt.5 <- cluster3.genera.melt.3 %>%
  group_by(Genus) %>%
  summarise(Abundance = sum(Abundance))

cluster3.genera.melt.5 <- as.data.frame(cluster3.genera.melt.5)
rownames(cluster3.genera.melt.5) <- cluster3.genera.melt.5$Genus

## label + percent
cluster3.genera.pie.1 <- plot_ly(cluster3.genera.melt.5, labels = ~Genus, values = ~Abundance, type = 'pie',marker = list(colors = core.noncore,line = list(color = '#FFFFFF', width = 2)), textposition = 'inside', textinfo = 'label+percent', insidetextfont = list(color = 'grey32'),hoverinfo = 'text') %>% layout(title = 'South Malayasia (cluster3) core and non-core genera',
                                                                                                                                                                                                                                                                                                                        xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                                                                                                                                                                                                                                                                                                                        yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))

##  percent
cluster3.genera.pie.2 <- plot_ly(cluster3.genera.melt.5, labels = ~Genus, values = ~Abundance, type = 'pie',marker = list(colors = core.noncore,line = list(color = '#FFFFFF', width = 2)), textposition = 'inside', textinfo = 'percent', insidetextfont = list(color = 'grey32'),hoverinfo = 'text') %>% layout(title = 'South Malayasia (cluster3) core and non-core genera',
                                                                                                                                                                                                                                                                                                                  xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                                                                                                                                                                                                                                                                                                                  yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))

## no labels
cluster3.genera.pie.3 <- plot_ly(cluster3.genera.melt.5, labels = ~Genus, values = ~Abundance, type = 'pie',marker = list(colors = core.noncore,line = list(color = '#FFFFFF', width = 2)), textinfo = 'none', insidetextfont = list(color = '#FFFFFF'),hoverinfo = 'text') %>% layout(title = 'South Malayasia (cluster3) core and non-core genera',
                                                                                                                                                                                                                                                                                       xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                                                                                                                                                                                                                                                                                       yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))

cluster3.genera.pie


## core genera pie chart ####
cluster3.core.genera <- cluster3.genera %>% subset_taxa(Genus %in% South.Malaysia.95.core)

#sanity check
test <- as.data.frame(cluster3.core.genera@tax_table)
test$Genus
South.Malaysia.95.core
identical(South.Malaysia.95.core,sort(test$Genus))

cluster3.core.genera.perc <- transform_sample_counts(cluster3.core.genera, function(OTU) (OTU/sum(OTU))*100)

test <- as.data.frame(cluster3.core.genera.perc@tax_table)

sample_sums(cluster3.core.genera.perc)

cluster3.core.genera.perc.melt <- psmelt(cluster3.core.genera.perc)

#change to character for easy-adjusted level
cluster3.core.genera.perc.melt$Genus <- as.character(cluster3.core.genera.perc.melt$Genus)

# Group according to faceting variable to get mean and median for further filtering
cluster3.core.genera.perc.melt <- cluster3.core.genera.perc.melt %>%
  group_by(Genus) %>%
  mutate(median=median(Abundance), mean=mean(Abundance))


### ## get mean abundances of genus identified at each cluster
cluster3.core.genera.perc.melt.2 <- cluster3.core.genera.perc.melt %>%
  group_by(Genus) %>%
  summarise(Abundance = mean(Abundance)) 


cluster3.core.genera.perc.melt.3 <- cluster3.core.genera.perc.melt.2

keep.5 <- unique(cluster3.core.genera.perc.melt.3$Genus[cluster3.core.genera.perc.melt.3$Abundance >= 5]) 
cluster3.core.genera.perc.melt.3$Genus[!(cluster3.core.genera.perc.melt.3$Genus %in% keep.5)] <- "< 5%"

## to check if to include <5% as others in pie-chart for core-genera, but others make up 43%
cluster3.core.genera.perc.melt.4 <- cluster3.core.genera.perc.melt.3 %>%
  group_by(Genus) %>%
  summarise(Abundance = sum(Abundance))


# Add label position :  does not apply correctly on plot
cluster3.core.genera.perc.melt.5 <- cluster3.core.genera.perc.melt.4 %>%
  arrange(desc(Genus))

cluster3.core.genera.perc.melt.5 <- as.data.frame(cluster3.core.genera.perc.melt.5)
rownames(cluster3.core.genera.perc.melt.5) <- cluster3.core.genera.perc.melt.5$Genus


#label + percent
cluster3.core.genera.pie.1 <- plot_ly(cluster3.core.genera.perc.melt.5, labels = ~Genus, values = ~Abundance, type = 'pie', marker = list(colors = cluster3.core, line = list(color = '#FFFFFF', width = 2)), textposition = 'inside', textinfo = 'label+percent', insidetextfont = list(color = 'grey32'),hoverinfo = 'text') %>% layout(title = 'South Malaysia (cluster3) core genera',
                                                                                                                                                                                                                                                                                                                                          xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                                                                                                                                                                                                                                                                                                                                          yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))

## percent
cluster3.core.genera.pie.2 <- plot_ly(cluster3.core.genera.perc.melt.5, labels = ~Genus, values = ~Abundance, type = 'pie', marker = list(colors = cluster3.core,line = list(color = '#FFFFFF', width = 2)), textposition = 'inside', textinfo = 'percent', insidetextfont = list(color = 'grey32')) %>% layout(title = 'South Malaysia (cluster3) core genera',
                                                                                                                                                                                                                                                                                                                xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                                                                                                                                                                                                                                                                                                                yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))

## no label
cluster3.core.genera.pie.3 <- plot_ly(cluster3.core.genera.perc.melt.5, labels = ~Genus, values = ~Abundance, type = 'pie', marker = list(colors = cluster3.core,line = list(color = '#FFFFFF', width = 2)), textposition = 'inside', textinfo = 'none', insidetextfont = list(color = '#FFFFFF')) %>% layout(title = 'South Malaysia (cluster3) core genera',
                                                                                                                                                                                                                                                                                                              xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                                                                                                                                                                                                                                                                                                              yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))
cluster3.core.genera.pie
cluster3.genera.pie


## Cluster 4 ####

cluster4 <- 
  prune_taxa(taxa_sums(cluster4.PCoA.unweighted.unifrac) > 0, cluster4.PCoA.unweighted.unifrac)

cluster4.genera <- tax_glom(cluster4, taxrank = "Genus")

test <- as.data.frame(cluster4.genera@tax_table)

cluster4.genera.perc <- transform_sample_counts(cluster4.genera, function(OTU) (OTU/sum(OTU))*100)

cluster4.genera.melt <- psmelt(cluster4.genera.perc)

#change to character for easy-adjusted level
cluster4.genera.melt$Genus <- as.character(cluster4.genera.melt$Genus)


# Group according to faceting variable to get mean and median for further filtering
cluster4.genera.melt <- cluster4.genera.melt %>%
  group_by(Genus) %>%
  mutate(median=median(Abundance), mean=mean(Abundance))


## get mean abundances of genus identified at each cluster
cluster4.genera.melt.2 <- cluster4.genera.melt %>%
  group_by(Genus) %>%
  summarise(Abundance = mean(Abundance)) 

# to label core and non-core genus ###
cluster4.genera.melt.3 <- cluster4.genera.melt.2
cluster4.genera.melt.3$Genus[!(cluster4.genera.melt.3$Genus %in% Singapore.95.core)] <- "non-core"
cluster4.genera.melt.3$Genus[(cluster4.genera.melt.3$Genus %in% Singapore.95.core)] <- "core"

cluster4.genera.melt.5 <- cluster4.genera.melt.3 %>%
  group_by(Genus) %>%
  summarise(Abundance = sum(Abundance))

cluster4.genera.melt.5 <- as.data.frame(cluster4.genera.melt.5)
rownames(cluster4.genera.melt.5) <- cluster4.genera.melt.5$Genus

## label + percent
cluster4.genera.pie.1 <- plot_ly(cluster4.genera.melt.5, labels = ~Genus, values = ~Abundance, type = 'pie',marker = list(colors = core.noncore,line = list(color = '#FFFFFF', width = 2)), textposition = 'inside', textinfo = 'label+percent', insidetextfont = list(color = 'grey32'),hoverinfo = 'text') %>% layout(title = 'Singapore (cluster4) core and non-core genera',
                                                                                                                                                                                                                                                                                                                        xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                                                                                                                                                                                                                                                                                                                        yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))

##percent
cluster4.genera.pie.2 <- plot_ly(cluster4.genera.melt.5, labels = ~Genus, values = ~Abundance, type = 'pie',marker = list(colors = core.noncore,line = list(color = '#FFFFFF', width = 2)), textposition = 'inside', textinfo = 'percent', insidetextfont = list(color = 'grey32'),hoverinfo = 'text') %>% layout(title = 'Singapore (cluster4) core and non-core genera',
                                                                                                                                                                                                                                                                                                                  xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                                                                                                                                                                                                                                                                                                                  yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))

##no labels
cluster4.genera.pie.3 <- plot_ly(cluster4.genera.melt.5, labels = ~Genus, values = ~Abundance, type = 'pie',marker = list(colors = core.noncore,line = list(color = '#FFFFFF', width = 2)), textposition = 'inside', textinfo = 'none', insidetextfont = list(color = '#FFFFFF'),hoverinfo = 'text') %>% layout(title = 'Singapore (cluster4) core and non-core genera',
                                                                                                                                                                                                                                                                                                                xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                                                                                                                                                                                                                                                                                                                yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))

cluster4.genera.pie


## core genera pie chart ####
cluster4.core.genera <- cluster4.genera %>% subset_taxa(Genus %in% Singapore.95.core)

#sanity check
test <- as.data.frame(cluster4.core.genera@tax_table)
test$Genus
Singapore.95.core
identical(Singapore.95.core,sort(test$Genus))

cluster4.core.genera.perc <- transform_sample_counts(cluster4.core.genera, function(OTU) (OTU/sum(OTU))*100)

test <- as.data.frame(cluster4.core.genera.perc@tax_table)

sample_sums(cluster4.core.genera.perc)

cluster4.core.genera.perc.melt <- psmelt(cluster4.core.genera.perc)

#change to character for easy-adjusted level
cluster4.core.genera.perc.melt$Genus <- as.character(cluster4.core.genera.perc.melt$Genus)

# Group according to faceting variable to get mean and median for further filtering
cluster4.core.genera.perc.melt <- cluster4.core.genera.perc.melt %>%
  group_by(Genus) %>%
  mutate(median=median(Abundance), mean=mean(Abundance))


### ## get mean abundances of genus identified at each cluster
cluster4.core.genera.perc.melt.2 <- cluster4.core.genera.perc.melt %>%
  group_by(Genus) %>%
  summarise(Abundance = mean(Abundance)) 


cluster4.core.genera.perc.melt.3 <- cluster4.core.genera.perc.melt.2

keep.5 <- unique(cluster4.core.genera.perc.melt.3$Genus[cluster4.core.genera.perc.melt.3$Abundance >= 5]) 
cluster4.core.genera.perc.melt.3$Genus[!(cluster4.core.genera.perc.melt.3$Genus %in% keep.5)] <- "< 5%"

## to check if to include <5% as others in pie-chart for core-genera, but others make up 43%
cluster4.core.genera.perc.melt.4  <- cluster4.core.genera.perc.melt.3 %>%
  group_by(Genus) %>%
  summarise(Abundance = sum(Abundance))


# Add label position :  does not apply correctly on plot
cluster4.core.genera.perc.melt.5 <- cluster4.core.genera.perc.melt.4 %>%
  arrange(desc(Genus)) 

cluster4.core.genera.perc.melt.5 <- as.data.frame(cluster4.core.genera.perc.melt.5)
rownames(cluster4.core.genera.perc.melt.5) <- cluster4.core.genera.perc.melt.5$Genus


## label + percent
cluster4.core.genera.pie.1 <- plot_ly(cluster4.core.genera.perc.melt.5, labels = ~Genus, values = ~Abundance, type = 'pie', marker = list(colors = cluster4.core,line = list(color = '#FFFFFF', width = 2)), textposition = 'inside', textinfo = 'label+percent', insidetextfont = list(color = 'grey32'),hoverinfo = 'text') %>% layout(title = 'Singapore (cluster4) core genera with',
                                                                                                                                                                                                                                                                                                                                         xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                                                                                                                                                                                                                                                                                                                                         yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))

# percent
cluster4.core.genera.pie.2 <- plot_ly(cluster4.core.genera.perc.melt.5, labels = ~Genus, values = ~Abundance, type = 'pie', marker = list(colors = cluster4.core,line = list(color = '#FFFFFF', width = 2)), textposition = 'inside', textinfo = 'percent', insidetextfont = list(color = 'grey32')) %>% layout(title = 'Singapore (cluster4) core genera',
                                                                                                                                                                                                                                                                                                                xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                                                                                                                                                                                                                                                                                                                yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))

# no labels
cluster4.core.genera.pie.3 <- plot_ly(cluster4.core.genera.perc.melt.5, labels = ~Genus, values = ~Abundance, type = 'pie', marker = list(colors = cluster4.core,line = list(color = '#FFFFFF', width = 2)), textposition = 'inside', textinfo = 'none', insidetextfont = list(color = '#FFFFFF')) %>% layout(title = 'Singapore (cluster4) core genera',
                                                                                                                                                                                                                                                                                                              xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                                                                                                                                                                                                                                                                                                              yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))


cluster4.core.genera.pie
cluster4.genera.pie


## Cluster 5 ####

cluster5 <- 
  prune_taxa(taxa_sums(cluster5.PCoA.unweighted.unifrac) > 0, cluster5.PCoA.unweighted.unifrac)

cluster5.genera <- tax_glom(cluster5, taxrank = "Genus")

test <- as.data.frame(cluster5.genera@tax_table)

cluster5.genera.perc <- transform_sample_counts(cluster5.genera, function(OTU) (OTU/sum(OTU))*100)

cluster5.genera.melt <- psmelt(cluster5.genera.perc)

#change to character for easy-adjusted level
cluster5.genera.melt$Genus <- as.character(cluster5.genera.melt$Genus)


# Group according to faceting variable to get mean and median for further filtering
cluster5.genera.melt <- cluster5.genera.melt %>%
  group_by(Genus) %>%
  mutate(median=median(Abundance), mean=mean(Abundance))


## get mean abundances of genus identified at each cluster
cluster5.genera.melt.2 <- cluster5.genera.melt %>%
  group_by(Genus) %>%
  summarise(Abundance = mean(Abundance)) 

# to label core and non-core genus ###
South.Thailand.95.core
cluster5.genera.melt.3 <- cluster5.genera.melt.2
cluster5.genera.melt.3$Genus[!(cluster5.genera.melt.3$Genus %in% South.Thailand.95.core)] <- "non-core"
cluster5.genera.melt.3$Genus[(cluster5.genera.melt.3$Genus %in% South.Thailand.95.core)] <- "core"

cluster5.genera.melt.5 <- cluster5.genera.melt.3 %>%
  group_by(Genus) %>%
  summarise(Abundance = sum(Abundance))

cluster5.genera.melt.5 <- as.data.frame(cluster5.genera.melt.5)
rownames(cluster5.genera.melt.5) <- cluster5.genera.melt.5$Genus

# label + percent
cluster5.genera.pie.1 <- plot_ly(cluster5.genera.melt.5, labels = ~Genus, values = ~Abundance, type = 'pie',marker = list(colors = core.noncore,line = list(color = '#FFFFFF', width = 2)), textposition = 'inside', textinfo = 'label+percent', insidetextfont = list(color = 'grey32'),hoverinfo = 'text') %>% layout(title = 'South Thailand (cluster5) core and non-core genera',
                                                                                                                                                                                                                                                                                                                        xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                                                                                                                                                                                                                                                                                                                        yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))

# percent
cluster5.genera.pie.2 <- plot_ly(cluster5.genera.melt.5, labels = ~Genus, values = ~Abundance, type = 'pie',marker = list(colors = core.noncore,line = list(color = '#FFFFFF', width = 2)), textposition = 'inside', textinfo = 'percent', insidetextfont = list(color = 'grey32'),hoverinfo = 'text') %>% layout(title = 'South Thailand (cluster5) core and non-core genera',
                                                                                                                                                                                                                                                                                                                  xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                                                                                                                                                                                                                                                                                                                  yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))

# no labels
cluster5.genera.pie.3 <- plot_ly(cluster5.genera.melt.5, labels = ~Genus, values = ~Abundance, type = 'pie',marker = list(colors = core.noncore,line = list(color = '#FFFFFF', width = 2)), textposition = 'inside', textinfo = 'none', insidetextfont = list(color = '#FFFFFF'),hoverinfo = 'text') %>% layout(title = 'South Thailand (cluster5) core and non-core genera',
                                                                                                                                                                                                                                                                                                                xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                                                                                                                                                                                                                                                                                                                yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))

cluster5.genera.pie


## core genera pie chart ####
cluster5.core.genera <- cluster5.genera %>% subset_taxa(Genus %in% South.Thailand.95.core)

#sanity check
test <- as.data.frame(cluster5.core.genera@tax_table)
test$Genus
South.Thailand.95.core
identical(South.Thailand.95.core,sort(test$Genus))

cluster5.core.genera.perc <- transform_sample_counts(cluster5.core.genera, function(OTU) (OTU/sum(OTU))*100)

test <- as.data.frame(cluster5.core.genera.perc@tax_table)

sample_sums(cluster5.core.genera.perc)

cluster5.core.genera.perc.melt <- psmelt(cluster5.core.genera.perc)

#change to character for easy-adjusted level
cluster5.core.genera.perc.melt$Genus <- as.character(cluster5.core.genera.perc.melt$Genus)

# Group according to faceting variable to get mean and median for further filtering
cluster5.core.genera.perc.melt <- cluster5.core.genera.perc.melt %>%
  group_by(Genus) %>%
  mutate(median=median(Abundance), mean=mean(Abundance))


### ## get mean abundances of genus identified at each cluster
cluster5.core.genera.perc.melt.2 <- cluster5.core.genera.perc.melt %>%
  group_by(Genus) %>%
  summarise(Abundance = mean(Abundance)) 


cluster5.core.genera.perc.melt.3 <- cluster5.core.genera.perc.melt.2

keep.5 <- unique(cluster5.core.genera.perc.melt.3$Genus[cluster5.core.genera.perc.melt.3$Abundance >= 5]) 
cluster5.core.genera.perc.melt.3$Genus[!(cluster5.core.genera.perc.melt.3$Genus %in% keep.5)] <- "< 5%"

## to check if to include <5% as others in pie-chart for core-genera, but others make up 43%
cluster5.core.genera.perc.melt.4 <- cluster5.core.genera.perc.melt.3 %>%
  group_by(Genus) %>%
  summarise(Abundance = sum(Abundance))


# Add label position :  does not apply correctly on plot
cluster5.core.genera.perc.melt.5 <- cluster5.core.genera.perc.melt.4 %>%
  arrange(desc(Genus)) 

cluster5.core.genera.perc.melt.5 <- as.data.frame(cluster5.core.genera.perc.melt.5)
rownames(cluster5.core.genera.perc.melt.5) <- cluster5.core.genera.perc.melt.5$Genus


## label + percent
cluster5.core.genera.pie.1 <- plot_ly(cluster5.core.genera.perc.melt.5, labels = ~Genus, values = ~Abundance, type = 'pie', marker = list(colors = cluster5.core,line = list(color = '#FFFFFF', width = 2)), textposition = 'inside', textinfo = 'label+percent', insidetextfont = list(color = 'grey32'),hoverinfo = 'text') %>% layout(title = 'South Thailand (cluster5) core genera',
                                                                                                                                                                                                                                                                                                                                         xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                                                                                                                                                                                                                                                                                                                                         yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))


## percent
cluster5.core.genera.pie.2 <- plot_ly(cluster5.core.genera.perc.melt.5, labels = ~Genus, values = ~Abundance, type = 'pie', marker = list(colors = cluster5.core,line = list(color = '#FFFFFF', width = 2)), textposition = 'inside', textinfo = 'percent', insidetextfont = list(color = 'grey32')) %>% layout(title = 'South Thailand (cluster5) core genera',
                                                                                                                                                                                                                                                                                                                xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                                                                                                                                                                                                                                                                                                                yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))


## no labels
cluster5.core.genera.pie.3 <- plot_ly(cluster5.core.genera.perc.melt.5, labels = ~Genus, values = ~Abundance, type = 'pie', marker = list(colors = cluster5.core,line = list(color = '#FFFFFF', width = 2)), textposition = 'inside', textinfo = 'none', insidetextfont = list(color = '#FFFFFF')) %>% layout(title = 'South Thailand (cluster5) core genera',
                                                                                                                                                                                                                                                                                                              xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                                                                                                                                                                                                                                                                                                              yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))

cluster5.core.genera.pie
cluster5.genera.pie


## Cluster 6 ####

cluster6 <- 
  prune_taxa(taxa_sums(cluster6.PCoA.unweighted.unifrac) > 0, cluster6.PCoA.unweighted.unifrac)

cluster6.genera <- tax_glom(cluster6, taxrank = "Genus")

test <- as.data.frame(cluster6.genera@tax_table)

cluster6.genera.perc <- transform_sample_counts(cluster6.genera, function(OTU) (OTU/sum(OTU))*100)

cluster6.genera.melt <- psmelt(cluster6.genera.perc)

#change to character for easy-adjusted level
cluster6.genera.melt$Genus <- as.character(cluster6.genera.melt$Genus)


# Group according to faceting variable to get mean and median for further filtering
cluster6.genera.melt <- cluster6.genera.melt %>%
  group_by(Genus) %>%
  mutate(median=median(Abundance), mean=mean(Abundance))


## get mean abundances of genus identified at each cluster
cluster6.genera.melt.2 <- cluster6.genera.melt %>%
  group_by(Genus) %>%
  summarise(Abundance = mean(Abundance)) 

# to label core and non-core genus ###
cluster6.genera.melt.3 <- cluster6.genera.melt.2
cluster6.genera.melt.3$Genus[!(cluster6.genera.melt.3$Genus %in% Central.Thailand.95.core)] <- "non-core"
cluster6.genera.melt.3$Genus[(cluster6.genera.melt.3$Genus %in% Central.Thailand.95.core)] <- "core"

cluster6.genera.melt.5 <- cluster6.genera.melt.3 %>%
  group_by(Genus) %>%
  summarise(Abundance = sum(Abundance))

cluster6.genera.melt.5 <- as.data.frame(cluster6.genera.melt.5)
rownames(cluster6.genera.melt.5) <- cluster6.genera.melt.5$Genus

## label + percent
cluster6.genera.pie.1 <- plot_ly(cluster6.genera.melt.5, labels = ~Genus, values = ~Abundance, type = 'pie',marker = list(colors = core.noncore,line = list(color = '#FFFFFF', width = 2)), textposition = 'inside', textinfo = 'label+percent', insidetextfont = list(color = 'grey32'),hoverinfo = 'text') %>% layout(title = 'Central Thailand (cluster6) core and non-core genera',
                                                                                                                                                                                                                                                                                                                        xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                                                                                                                                                                                                                                                                                                                        yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))

## percent
cluster6.genera.pie.2 <- plot_ly(cluster6.genera.melt.5, labels = ~Genus, values = ~Abundance, type = 'pie',marker = list(colors = core.noncore,line = list(color = '#FFFFFF', width = 2)), textposition = 'inside', textinfo = 'percent', insidetextfont = list(color = 'grey32'),hoverinfo = 'text') %>% layout(title = 'Central Thailand (cluster6) core and non-core genera',
                                                                                                                                                                                                                                                                                                                  xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                                                                                                                                                                                                                                                                                                                  yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))


## no labels
cluster6.genera.pie.3 <- plot_ly(cluster6.genera.melt.5, labels = ~Genus, values = ~Abundance, type = 'pie',marker = list(colors = core.noncore,line = list(color = '#FFFFFF', width = 2)), textposition = 'inside', textinfo = 'none', insidetextfont = list(color = '#FFFFFF'),hoverinfo = 'text') %>% layout(title = 'Central Thailand (cluster6) core and non-core genera',
                                                                                                                                                                                                                                                                                                                xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                                                                                                                                                                                                                                                                                                                yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))

cluster6.genera.pie


## core genera pie chart ####
cluster6.core.genera <- cluster6.genera %>% subset_taxa(Genus %in% Central.Thailand.95.core)

#sanity check
test <- as.data.frame(cluster6.core.genera@tax_table)
test$Genus
Central.Thailand.95.core
identical(Central.Thailand.95.core,sort(test$Genus))

cluster6.core.genera.perc <- transform_sample_counts(cluster6.core.genera, function(OTU) (OTU/sum(OTU))*100)

test <- as.data.frame(cluster6.core.genera.perc@tax_table)

sample_sums(cluster6.core.genera.perc)

cluster6.core.genera.perc.melt <- psmelt(cluster6.core.genera.perc)

#change to character for easy-adjusted level
cluster6.core.genera.perc.melt$Genus <- as.character(cluster6.core.genera.perc.melt$Genus)

# Group according to faceting variable to get mean and median for further filtering
cluster6.core.genera.perc.melt <- cluster6.core.genera.perc.melt %>%
  group_by(Genus) %>%
  mutate(median=median(Abundance), mean=mean(Abundance))


### ## get mean abundances of genus identified at each cluster
cluster6.core.genera.perc.melt.2 <- cluster6.core.genera.perc.melt %>%
  group_by(Genus) %>%
  summarise(Abundance = mean(Abundance)) 


cluster6.core.genera.perc.melt.3 <- cluster6.core.genera.perc.melt.2

keep.5 <- unique(cluster6.core.genera.perc.melt.3$Genus[cluster6.core.genera.perc.melt.3$Abundance >= 5]) 
cluster6.core.genera.perc.melt.3$Genus[!(cluster6.core.genera.perc.melt.3$Genus %in% keep.5)] <- "< 5%"

## to check if to include <5% as others in pie-chart for core-genera, but others make up 43%
cluster6.core.genera.perc.melt.4  <- cluster6.core.genera.perc.melt.3 %>%
  group_by(Genus) %>%
  summarise(Abundance = sum(Abundance))


# Add label position :  does not apply correctly on plot
cluster6.core.genera.perc.melt.5 <- cluster6.core.genera.perc.melt.4 %>%
  arrange(desc(Genus)) 

cluster6.core.genera.perc.melt.5 <- as.data.frame(cluster6.core.genera.perc.melt.5)
rownames(cluster6.core.genera.perc.melt.5) <- cluster6.core.genera.perc.melt.5$Genus


## label + percent
cluster6.core.genera.pie.1 <- plot_ly(cluster6.core.genera.perc.melt.5, labels = ~Genus, values = ~Abundance, type = 'pie', marker = list(colors = cluster6.core,line = list(color = '#FFFFFF', width = 2)), textposition = 'inside', textinfo = 'label+percent', insidetextfont = list(color = 'grey32'),hoverinfo = 'text') %>% layout(title = 'Central Thailand (cluster6) core genera',
                                                                                                                                                                                                                                                                                                                                         xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                                                                                                                                                                                                                                                                                                                                         yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))

## percent
cluster6.core.genera.pie.2 <- plot_ly(cluster6.core.genera.perc.melt.5, labels = ~Genus, values = ~Abundance, type = 'pie', marker = list(colors = cluster6.core,line = list(color = '#FFFFFF', width = 2)), textposition = 'inside', textinfo = 'percent', insidetextfont = list(color = 'grey32')) %>% layout(title = 'Central Thailand (cluster6) core genera',
                                                                                                                                                                                                                                                                                                                xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                                                                                                                                                                                                                                                                                                                yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))


## no labels
cluster6.core.genera.pie.3 <- plot_ly(cluster6.core.genera.perc.melt.5, labels = ~Genus, values = ~Abundance, type = 'pie', marker = list(colors = cluster6.core,line = list(color = '#FFFFFF', width = 2)), textposition = 'inside', textinfo = 'none', insidetextfont = list(color = '#FFFFFF')) %>% layout(title = 'Central Thailand (cluster6) core genera',
                                                                                                                                                                                                                                                                                                              xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                                                                                                                                                                                                                                                                                                              yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))
cluster6.core.genera.pie
cluster6.genera.pie


###################################################################################

### Fig. 5A. Putative biotic interactions were estimated using co-occurrence network analysis. Modules of interaction are shown by circles sharing the same colour fill ####
### Fig. S7. Co-occurrence network analysis indicating ecological groups of interacting ASVs. ####
### Fig. 5B. Chord diagram illustrating major interactions between abundant ASVs clustered at genus level ####

######################################################################################

### Network analyses with SpiecEASI and NETCOMI ####

com.rarefied.min.int.wide.exclude.1.3.rooted <- as.matrix(as.data.frame(com.rarefied.min.int.exclude.1.3.rooted@otu_table))

### best to run both netConstruct and netAnalyze on high performance computers
ASV_net_spieceasi_mb <- netConstruct(com.rarefied.min.int.wide.exclude.1.3.rooted, 
                     measure = "spieceasi",
                     measurePar = list(method = "mb", sel.criterion='stars', 
                      lambda.min.ratio=1e-2, nlambda=20,
                        pulsar.params=list(rep.num=100, ncores=24),
                                       symBetaMode = "ave"),
                     sparsMethod = "none",
                     normMethod = "none",
                     verbose = 3)

ASV_anlayze_spieceasi_mb <- netAnalyze(ASV_net_spieceasi_mb, clustMethod = "cluster_fast_greedy",
                            hubPar = "eigenvector", weightDeg = FALSE, normDeg = FALSE)

summary(ASV_anlayze_spieceasi_mb, showCompSize= TRUE, showGlobal = TRUE, showCluster = TRUE, showHubs = TRUE, showCentr = "all", numbNodes = 30L)

pdf(file="ASV.spieceasi.mb.netcomi.nolabels.pdf", width = 25, height = 16)
plot(ASV_anlayze_spieceasi_mb,
     edgeTranspLow = 0,
     edgeTranspHigh = 40,
     edgeInvisFilter = "threshold",
     charToRm = "ASV",
     edgeInvisPar = 0.1,
     nodeColor = "cluster",
     colorVec = c("mediumaquamarine", "thistle3", "hotpink4", "darkgoldenrod4", "yellow3", "steelblue"),
     sameClustCol = TRUE,
     nodeSize = "eigenvector",
     repulsion = 0.8,
     rmSingles = TRUE,
     labelFont = 1,
     labelScale = FALSE, 
     cexLabels = 0,
     nodeSizeSpread = 3,
     cexNodes = 2,
     hubBorderCol = "firebrick2",
     hubBorderWidth = 3,
     highlightHubs = TRUE,
     title1 = "Network analyses on ASV level with spieceasi using method 'mb' (single nodes removed & edge filtered)", 
     showTitle = TRUE,
     cexTitle = 1.2)

legend(0.8, 1.2, cex = 2.5, title = "estimated correlation:",
       legend = c("+","-"), lty = 1, lwd = 4, col = c("#009900","red"),
       bty = "n", horiz = TRUE)

dev.off()


### add additional varibales to separate taxtable to add 2 additional columns

# Declare photosynthetoc and heterotrophic groups
## 1 Chloroflexi : Phylum == "Chloroflexi" & Class  == "Chloroflexia"

## 2 Cyanobacteria : Phylum == "Cyanobacteria" & Class  == "Cyanobacteriia"

## 3 Proteobacteria Genus:
proteobacteria.order <- c("DSSF69","Elioraea","Methylobacterium-Methylorubrum","Rhodomicrobium","Roseomonas","Sandaracinobacter","Tabrizicola","Unassigned Rhodobacteraceae (Family)","Unassigned Sphingomonadaceae (Family)","AAP99","Allochromatium","Caldimonas","Curvibacter","DSSD61","Thiolamprovum","Unassigned B1-7BS (Family)","Unassigned Burkholderiales (Order)","Unassigned Comamonadaceae (Family)","Unassigned Gammaproteobacteria (Class)","Unassigned Rhodocyclaceae (Family)","Unassigned Sutterellaceae (Family)","Z-35","Unassigned Hydrogenophilaceae (Family)")

## 4 Other photosynthetic : 
other.photo <- c("Chlorobiales","Chloracidobacteriales")

## All photosynthetic : ## photosynthetic.class <- c("Cyanobacteriia","Chloroflexia")
# ( Class %in% photosynthetic.class | Order %in% other.photo | Genus %in% proteobacteria.genus )

# Phylum == 'Cyanobacteria' & !Class == 'Cyanobacteriia' ~ 'Non-photosynthetic Cyanobacteria'

## Chemolithoautotrophs

net_spieceasi_taxtable <- as.data.frame(com.rarefied.min.int.exclude.1.3.rooted@tax_table) %>% rownames_to_column(var = "ASV")  %>% as_tibble()

#rownames(net_spieceasi_taxtable) <- net_spieceasi_taxtable$ASV

net_spieceasi_taxtable <- net_spieceasi_taxtable %>%
  mutate(PhotoHetero1 = case_when(Phylum == 'Aquificota' & Class == 'Aquificae' ~ 'Chemolithoautotroph',
                                  Phylum == 'Calditrichota' & Class == 'Calditrichia' ~ 'Chemolithoautotroph',
                                  Phylum == 'Desulfobacterota' ~ 'Chemolithoautotroph',
                                  Phylum == 'Elusimicrobiota' ~ 'Chemolithoautotroph',
                                  Phylum == 'Nitrospirota' ~ 'Chemolithoautotroph',
                                  Phylum == 'Patescibacteria' ~ 'Chemolithoautotroph',
                                  Phylum == 'Sva0485' ~ 'Chemolithoautotroph',
                                  Class == 'Leptospirae' ~ 'Chemolithoautotroph',
                                  Class == 'Chloroflexia' ~ 'Chloroflexia',
                                  Class == 'Cyanobacteriia' ~ 'Cyanobacteriia', 
                                  Genus %in% proteobacteria.order ~ 'Proteobacteria', 
                                  Order == 'Chlorobiales' ~ 'Chlorobiales', 
                                  Order == 'Chloracidobacteriales' ~ 'Chloracidobacteriales', 
                                  .default = "Heterotrophs" ))

net_spieceasi_taxtable <- net_spieceasi_taxtable %>%
  mutate(PhotoHetero2 = case_when(PhotoHetero1 == 'Chemolithoautotroph' ~ 'Others',
                                  PhotoHetero1 == 'Heterotrophs' ~ 'Others', .default = "Photosynthetic" ))

saveRDS(net_spieceasi_taxtable,"net_spieceasi_taxtable.rds")

# Get photosynthetic names
spieceasi_taxtable <- as.matrix(net_spieceasi_taxtable)
#rownames(spieceasi_taxtable) <- net_spieceasi_taxtable$ASV
photo <- as.factor(spieceasi_taxtable[, "PhotoHetero1"])
photo2 <- as.factor(spieceasi_taxtable[, "PhotoHetero2"])
names(photo) <- spieceasi_taxtable[, "ASV"]
names(photo2) <- spieceasi_taxtable[, "ASV"]
table(photo)
table(photo2)

# Define phylum colors
#Non-photosynthetic/Heterotrophs/Others                  Photosynthetic 
#                            385                             187
#
#
# Chemolithoautotroph Chloracidobacteriales          Chlorobiales          Chloroflexia        Cyanobacteriia 
#                   41                     5                    11                    55                    74 
#         Heterotrophs        Proteobacteria 
#                  344

photocol <- c("goldenrod1","cyan", "darkblue", "coral3", "forestgreen", "azure4", "purple3")
photocol2 <- c( "azure4", "forestgreen")


pdf(file="ASV.spieceasi.mb.netcomi.1.nolabels.pdf", width = 25, height = 16)
plot(ASV_anlayze_spieceasi_mb,
     edgeTranspLow = 0,
     edgeTranspHigh = 40,
     edgeInvisFilter = "threshold",
     edgeInvisPar = 0.1,
     charToRm = "ASV",
     nodeSize = "eigenvector",
     nodeColor = "feature", 
     featVecCol = photo, 
     colorVec =  photocol,
     repulsion = 0.8,
     rmSingles = TRUE,
     labelFont = 1,
     labelScale = FALSE, 
     cexLabels = 0,
     nodeSizeSpread = 3,
     cexNodes = 2,
     hubBorderCol = "firebrick2",
     hubBorderWidth = 3,
     highlightHubs = TRUE,
     title1 = "Network analyses on ASV level with spieceasi using method 'mb' (single nodes removed & edge filtered)", 
     showTitle = TRUE,
     cexTitle = 1.2)

# Colors used in the legend should be equally transparent as in the plot
photocol_transp <- colToTransp(photocol, 20)

legend(-1.2, 1.2, cex = 1.8, pt.cex = 2.3, title = "Photosynthetic & other taxa:", 
       legend=levels(photo), col = photocol_transp, bty = "n", pch = 16) 


legend(0.8, 1.2, cex = 2.5, title = "estimated correlation:",
       legend = c("+","-"), lty = 1, lwd = 4, col = c("#009900","red"),
       bty = "n", horiz = TRUE)

dev.off()

### Chord Diagram ####

### all tables as data.frame first before combining into microtable object
microeco.otu <- t(as.matrix(as.data.frame(otu_table(com.rarefied.min.int.exclude.1.3.rooted))))
microeco.otu <- as.data.frame(microeco.otu)
sapply(microeco.otu,class)
class(microeco.otu)

tax <- (as.matrix(as.data.frame(tax_table(com.rarefied.min.int.exclude.1.3.rooted))))
tax <- as.data.frame(tax)

metadata <- as.data.frame(sample_data(com.rarefied.min.int.exclude.1.3.rooted)) %>% as_tibble()
colnames(metadata)

### change column headers accordingly
metadata.micronet <- metadata %>% select("Sample_id","Country","Region","Location.Code","Site","Location","Latittude","Longitude","Human.usage..Y.N.","Water.flow..Pool.Flowing.","Nitrate..ppm.","Nitrite..ppm.","Temp.adj...C.","Carbonate..ppm.","pH","Total.alkalinity..ppm.","Phosphate..ppm.","H2S..ppm.","Temp...C.","EC..mS.")

rownames(metadata.micronet) <- metadata.micronet$Sample_id

metadata.micronet <- as.data.frame(as.matrix(metadata.micronet))

metadata.micronet[11:ncol(metadata.micronet)] <- lapply(metadata.micronet[11:ncol(metadata.micronet)], as.numeric)


# Let's create a microtable object with more information
microeco.dataset <- microtable$new(sample_table = metadata.micronet, otu_table = microeco.otu, tax_table = tax, phylo_tree = rootedTree.rare)

saveRDS(microeco.dataset,"microeco.dataset.rds")


microeco_transnetwork_spieceasi  <- trans_network$new(dataset = microeco.dataset, cor_method = NULL,taxa_level = "OTU")

saveRDS(microeco_transnetwork_spieceasi,"microeco_transnetwork_spieceasi.rds")

# Run on hpc
# with list providing the same parameters as used to create spieceasi mb plot with spieceasi package
microeco_transnetwork_spieceasi$cal_network(network_method = "SpiecEasi", SpiecEasi_method = "mb", ...list(sel.criterion='stars' ,lambda.min.ratio=1e-2,nlambda=20, pulsar.params=list(rep.num=100, ncores =24)))

### read into R after importing from HPC
microeco_transnetwork_spieceasi <- readRDS("microeco_transnetwork_spieceasi.rds")


# invoke igraph cluster_fast_greedy function for this undirected network 
## 16 modules identified.
microeco_transnetwork_spieceasi$cal_module(method = "cluster_fast_greedy",
                                           module_name_prefix = "M")

# calculate network attributes
microeco_transnetwork_spieceasi$cal_network_attr()
microeco_transnetwork_spieceasi$res_network_attr

# get node properties, calculate node roles, i.e. Module hubs, Network hubs, Connectors and Peripherals
microeco_transnetwork_spieceasi$get_node_table(node_roles = TRUE)
microeco_transnetwork_spieceasi$res_node_table 

# Get the edge property table, including connected nodes, label and weight.
microeco_transnetwork_spieceasi$get_edge_table()
microeco_transnetwork_spieceasi$res_edge_table 

#Get the adjacency matrix from the network graph
microeco_transnetwork_spieceasi$get_adjacency_matrix()
microeco_transnetwork_spieceasi$res_adjacency_matrix

#plot the node classification in terms of the within-module connectivity and among-module connectivity.
# add_label = TRUE can be used to directly add text label for points
microeco_transnetwork_spieceasi$plot_taxa_roles(use_type = 1, add_label = TRUE)

# plot node roles with phylum information
microeco_transnetwork_spieceasi$plot_taxa_roles(use_type = 2)

#Now, we show the eigengene analysis of modules. The eigengene of a module, i.e. the first principal component of PCA, represents the main variance of the abundance in the species of the module.
microeco_transnetwork_spieceasi$cal_eigen()
microeco_transnetwork_spieceasi$res_eigen

#perform correlation heatmap to show the associations between eigengenes and environmental factors.
# create trans_env object
microeco_transnetwork_spieceasi_correlation <- trans_env$new(dataset = microeco.dataset, add_data = metadata.micronet[, 24:31])

# calculate correlations
microeco_transnetwork_spieceasi_correlation$cal_cor(add_abund_table = microeco_transnetwork_spieceasi$res_eigen)


# default parameter represents using igraph plot.igraph function
microeco_transnetwork_spieceasi$plot_network(method = "igraph", layout = layout_with_kk,node_color = "module")

## The function cal_sum_links can sum the links (edge) number from one taxa to another or within the same taxa. The function plot_sum_links is used to show the result from the function cal_sum_links. This is very useful to fast see how many nodes are connected between different taxa or within one taxa. In terms of ‘Phylum’ level in the tutorial, the function cal_sum_links() sum the linkages number from one Phylum to another Phylum or the linkages in the same Phylum. So the numbers along the outside of the circular plot represent how many edges or linkages are related with the Phylum. For example, in terms of Proteobacteria, there are roughly total 900 edges associated with the OTUs in Proteobacteria, in which roughly 200 edges connect both OTUs in Proteobacteria and roughly 150 edges connect the OTUs from Proteobacteria with the OTUs from Chloroflexi.
microeco_transnetwork_spieceasi$cal_sum_links(taxa_level = "Family")
microeco_transnetwork_spieceasi$res_sum_links_pos
#devtools::install_github("mattflor/chorddiag")

library(chorddiag)

color_values = RColorBrewer::brewer.pal(12, "Paired")

color_values_2 = c("#3B4D16","#950D3D", "#848F22", "#A7A3FE","#B2C9B2", "#8CCBD3", "#0DD7FD", "#99D584", "#0D8CA8","#DFC500","#2A7F72","#FEAB2E", "#4B00FD", "#26E0A9" ,"#CC003D")

microeco_transnetwork_spieceasi$plot_sum_links(plot_pos = TRUE, plot_num = 15, color_values = color_values_2, groupnameFontsize = 8, showTicks = FALSE, groupnamePadding = 40)

microeco_transnetwork_spieceasi$plot_sum_links(plot_pos = FALSE, plot_num = 15, color_values = color_values_2)

#########################################################################################

### Fig. 4A. Mantel's test on the influence of abiotic variables on observed taxonomic composition for all taxa plus photosynthetic, chemoautotrophic, and chemoheterotrophic fractions of the community. ####
### Fig. 4B. Variance partioning ####

######################################################################################

## Mantel's test and correlation heatmap #####

# prepare data
tax <- as.data.frame(net_spieceasi_taxtable) 
tax <- as.data.frame(tax) %>% column_to_rownames(var = "ASV")


# Let's create a microtable object with more information
microeco.dataset.3 <- microtable$new(sample_table = metadata.micronet, otu_table = microeco.otu, tax_table = tax, phylo_tree = rootedTree.rare)

#### Mantel's test for all the data ####

#microeco.dataset.2 <- clone(microeco.dataset)
microeco.dataset.3$tidy_dataset()

#unifrac default FALSE; whether UniFrac index should be calculated,binary default FALSE; TRUE is used for jaccard and unweighted unifrac;
#microeco.dataset.2$cal_betadiv(unifrac = TRUE, binary = TRUE)
microeco.dataset.3$cal_betadiv(unifrac = TRUE, binary = TRUE)


# first perform mantel test 
# Nitrate and Nitrites excluded as they were zero
# all other env variables were considered as indication of human activity so not included.
## env_cols = 23:31: change column numbers accordingly
microeco.mantel <- trans_env$new(dataset = microeco.dataset.3, env_cols = 23:31)


#Mantel test between beta diversity matrix and environmental data.
microeco.mantel$cal_mantel(use_measure = "unwei_unifrac", method = "pearson" , partial_mantel = FALSE)

### use microeco.mantel.all to sibset for required columns
microeco.mantel.all <- data.frame(microeco.mantel$res_mantel)
microeco.mantel.2 <- data.frame(microeco.mantel$res_mantel) %>% .[, c(1, 2, 5, 7)]


# rename columns
colnames(microeco.mantel.2) <- c("spec", "env", "r", "p.value")

# generate interval data
microeco.mantel.2 %<>% dplyr::mutate(rd = cut(r, breaks = c(-Inf, 0.3, 0.6, Inf), labels = c("< 0.3", "0.3 - 0.6", ">= 0.6")), pd = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf), labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))

microeco.mantel$data_env


library(ggcor)
set_scale()

##### Mantel's with 2 groups : Photosynthetic and non-photosynthetic ####

# extract two phyla to show the steps
photosynthetic <- clone(microeco.dataset.3)
photosynthetic$tax_table <- photosynthetic$tax_table[photosynthetic$tax_table$PhotoHetero2 == "Photosynthetic", ]
#check.photo <- photosynthetic$tax_table 
photosynthetic$tidy_dataset()
photosynthetic$cal_betadiv(unifrac = TRUE, binary = TRUE)
#photosynthetic$beta_diversity

chemolithoautotroph <- clone(microeco.dataset.3)
chemolithoautotroph$tax_table <- chemolithoautotroph$tax_table[chemolithoautotroph$tax_table$PhotoHetero1 == "Chemolithoautotroph", ]
#check.chemo <- chemolithoautotroph$tax_table 
chemolithoautotroph$tidy_dataset()
chemolithoautotroph$cal_betadiv(unifrac = TRUE, binary = TRUE)

heterotroph <- clone(microeco.dataset.3)
heterotroph$tax_table <- heterotroph$tax_table[heterotroph$tax_table$PhotoHetero1 == "Heterotrophs", ]
check.hetero <- heterotroph$tax_table
heterotroph$tidy_dataset()
heterotroph$cal_betadiv(unifrac = TRUE, binary = TRUE)

#("DSSF69","Elioraea","Methylobacterium-Methylorubrum","Rhodomicrobium","Roseomonas","Sandaracinobacter","Tabrizicola","Unassigned Rhodobacteraceae (Family)","Unassigned Sphingomonadaceae (Family)","AAP99","Allochromatium","Caldimonas","Curvibacter","DSSD61","Thiolamprovum","Unassigned B1-7BS (Family)","Unassigned Burkholderiales (Order)","Unassigned Comamonadaceae (Family)","Unassigned Gammaproteobacteria (Class)","Unassigned Rhodocyclaceae (Family)","Unassigned Sutterellaceae (Family)","Z-35")

unique(check.chemo$Genus)

# first perform mantel test
photosynthetic.mantel <- trans_env$new(dataset = photosynthetic, env_cols = 23:31)
photosynthetic.mantel$cal_mantel(use_measure = "unwei_unifrac", method = "pearson" , partial_mantel = FALSE)

chemolithoautotroph.mantel <- trans_env$new(dataset = chemolithoautotroph, env_cols = 23:31)
chemolithoautotroph.mantel$cal_mantel(use_measure = "unwei_unifrac", method = "pearson" , partial_mantel = FALSE)

heterotroph.mantel <- trans_env$new(dataset = heterotroph, env_cols = 23:31)
heterotroph.mantel$cal_mantel(use_measure = "unwei_unifrac", method = "pearson" , partial_mantel = FALSE)

### use microeco.mantel.all to sibset for required columns
photosynthetic.mantel.all <- data.frame(photosynthetic.mantel$res_mantel)
photosynthetic.mantel.2 <- data.frame(photosynthetic.mantel$res_mantel) %>% .[, c(1, 2, 5, 7)]

chemolithoautotroph.mantel.all <- data.frame(chemolithoautotroph.mantel$res_mantel)
chemolithoautotroph.mantel.2 <- data.frame(chemolithoautotroph.mantel$res_mantel) %>% .[, c(1, 2, 5, 7)]

heterotroph.mantel.all <- data.frame(heterotroph.mantel$res_mantel)
heterotroph.mantel.2 <- data.frame(heterotroph.mantel$res_mantel) %>% .[, c(1, 2, 5, 7)]
heterotroph.mantel$res_mantel

# rename columns
colnames(photosynthetic.mantel.2) <- c("spec", "env", "r", "p.value")
colnames(chemolithoautotroph.mantel.2) <- c("spec", "env", "r", "p.value")
colnames(heterotroph.mantel.2) <- c("spec", "env", "r", "p.value")

# generate interval data
photosynthetic.mantel.2 %<>% dplyr::mutate(rd = cut(r, breaks = c(-Inf, 0.3, 0.6, Inf), labels = c("< 0.3", "0.3 - 0.6", ">= 0.6")), pd = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf), labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))

chemolithoautotroph.mantel.2 %<>% dplyr::mutate(rd = cut(r, breaks = c(-Inf, 0.3, 0.6, Inf), labels = c("< 0.3", "0.3 - 0.6", ">= 0.6")), pd = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf), labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))

heterotroph.mantel.2 %<>% dplyr::mutate(rd = cut(r, breaks = c(-Inf, 0.3, 0.6, Inf), labels = c("< 0.3", "0.3 - 0.6", ">= 0.6")), pd = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf), labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))


### Combine all 3 groups together

check.all <- microeco.mantel$dataset$tax_table

# extract a part of the results 
All <- microeco.mantel.2
Photosynthetic <- data.frame(spec = "Photosynthetic", photosynthetic.mantel$res_mantel) %>% .[, c(1, 3, 6, 8)]
Chemolithoautotroph <- data.frame(spec = "Chemolithoautotroph
", chemolithoautotroph.mantel$res_mantel) %>% .[, c(1, 3, 6, 8)]
Heterotroph <- data.frame(spec = "Heterotrophs", heterotroph.mantel$res_mantel) %>% .[, c(1, 3, 6, 8)]


# rename columns
colnames(Photosynthetic) <- colnames(Chemolithoautotroph) <- colnames(Heterotroph) <- c("spec", "env", "r", "p.value")

# generate interval data

All %<>% dplyr::mutate(rd = cut(r, breaks = c(-Inf, 0.3, 0.6, Inf), labels = c("< 0.3", "0.3 - 0.6", ">= 0.6")), pd = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf), labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))

Photosynthetic %<>% dplyr::mutate(rd = cut(r, breaks = c(-Inf, 0.3, 0.6, Inf), labels = c("< 0.3", "0.3 - 0.6", ">= 0.6")),
                                  pd = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf), labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))

Chemolithoautotroph %<>% dplyr::mutate(rd = cut(r, breaks = c(-Inf, 0.3, 0.6, Inf), labels = c("< 0.3", "0.3 - 0.6", ">= 0.6")),
                                       pd = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf), labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))

Heterotroph %<>% dplyr::mutate(rd = cut(r, breaks = c(-Inf, 0.3, 0.6, Inf), labels = c("< 0.3", "0.3 - 0.6", ">= 0.6")),
                               pd = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf), labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))

# cobine 3 tables
plot_table <- rbind(All, Photosynthetic , Chemolithoautotroph, Heterotroph )

library(ggcor)
set_scale()

main.groups <- quickcor(microeco.mantel$data_env, type = "upper", cor.test = TRUE, show.diag = FALSE) +
  geom_square() +
  geom_mark(sig.thres = 0.05, markonly = TRUE, color = "white") +
  anno_link(aes(colour = pd, size = rd), data = plot_table) +
  scale_size_manual(values = c(0.5, 1.5, 3)) +
  scale_colour_manual(values = c("#D95F02", "#1B9E77", "#A2A2A288")) +
  guides(size = guide_legend(title = "Mantel's r", override.aes = list(colour = "grey35"), order = 2),
         colour = guide_legend(title = "Mantel's p", override.aes = list(size = 3), order = 1),
         fill = guide_colorbar(title = "Pearson's r", order = 3)) 

pdf("main.groups.mantel.pdf", width = 10)
main.groups
dev.off()

## #Mantel statistic based on Pearson's correlation 
## Pearson's |r| is correlation between abiotic variables only.
pearsons<- data.frame(main.groups$data)
colnames(pearsons)

chosen5 <- c("Carbonate..ppm.","EC..mS.","Latittude","Longitude","pH")

pearsons.5 <- pearsons %>% filter(.row.names %in% chosen5 & .col.names %in% chosen5)

#### varpart ####

## Reducing the weight of rare species use hellinger transformation
otu <- as.matrix(as.data.frame(com.rarefied.min.int.exclude.1.3.rooted@otu_table))

metadata.varpart <- metadata.micronet

rownames(metadata.varpart) <- NULL

metadata.varpart <- data.frame(metadata.varpart)

colnames(metadata)

otu.hellinger.transform <- decostand(otu, method = "hellinger")

varp.chosen.2 <- varpart(otu.hellinger.transform, ~  EC..mS., ~ Latittude, ~ Carbonate..ppm., ~ pH, data = metadata.varpart)

#get partition table
summary(varp.chosen.2)

## indicates testable fractions
varp.chosen.2$part$fract

# plot the results into Venn's diagram (argument digits influences number of decimal digits shown in the diagram, Xnames the displayed names of variables, and bg background color of the diagram fractions; see ?varpart for details):

pdf("varpart.2.pdf", width = 8)
plot (varp.chosen.2, digits = 2, Xnames = c('EC (mS)','Latittude (D)', 'Carbonate (ppm)', 'pH'), bg = c('turquoise','khaki3', 'greenyellow', "mediumslateblue"))
dev.off()

## As stated in the varpart fucntion: "Use function ‘rda’ to test significance of fractions of interest"
otu.hellinger.chosen.2 <- rda(otu.hellinger.transform ~ Latittude + Carbonate..ppm. + EC..mS. + pH, data = metadata.varpart)

## ## The global model (fractions [a+b+c]):
anova(otu.hellinger.chosen.2)

###################################################################################

### Fig. 6A,B,C. Null models were employed to predict the influence of various evolutionary drivers on community assembly. Net relatedness index (betaNRI) estimates are shown by biogeographic region for all taxa (A), versus photosynthetic (B) and chemoheterotrophic (C) taxa ####
### Fig. S8. Net relatedness index (betaNRI) measures of the mean phylogenetic distance to the nearest taxon in the community for ecological groups of bacteria by biogeographic region.

######################################################################################

### change enviornmental variable columns accordingly

microeco.nullmodel <- trans_nullmodel$new(microeco.dataset.3, env_cols = 23:31)
saveRDS(microeco.nullmodel, "microeco.nullmodel.rds")

photosynthetic.nullmodel <- trans_nullmodel$new(photosynthetic, env_cols = 23:31)
saveRDS(photosynthetic.nullmodel, "photosynthetic.nullmodel.rds")

### create additional  null model objects for each of the photosynthetic groups:

chloroflexi <- clone(microeco.dataset.3)
chloroflexi$tax_table <- chloroflexi$tax_table[chloroflexi$tax_table$PhotoHetero1 == "Chloroflexia", ]
chloroflexi$tidy_dataset()
check.chloro <- chloroflexi$tax_table 
check.chloro <- chloroflexi$otu_table
chloroflexi$cal_betadiv(unifrac = TRUE, binary = TRUE)

chloroflexi.nullmodel <- trans_nullmodel$new(chloroflexi, env_cols = 23:31)
saveRDS(chloroflexi.nullmodel, "chloroflexi.nullmodel.rds")

cyanobacteria <- clone(microeco.dataset.3)
cyanobacteria$tax_table <- cyanobacteria$tax_table[cyanobacteria$tax_table$PhotoHetero1 == "Cyanobacteriia", ]
check.cyano <- cyanobacteria$tax_table 
cyanobacteria$tidy_dataset()
check.cyano <- cyanobacteria$otu_table 
cyanobacteria$cal_betadiv(unifrac = TRUE, binary = TRUE)

cyanobacteria.nullmodel <- trans_nullmodel$new(cyanobacteria, env_cols = 23:31)
saveRDS(cyanobacteria.nullmodel, "cyanobacteria.nullmodel.rds")


proteobacteria.genus <- c("DSSF69","Elioraea","Methylobacterium-Methylorubrum","Rhodomicrobium","Roseomonas","Sandaracinobacter","Tabrizicola","Unassigned Rhodobacteraceae (Family)","Unassigned Sphingomonadaceae (Family)","AAP99","Allochromatium","Caldimonas","Curvibacter","DSSD61","Thiolamprovum","Unassigned B1-7BS (Family)","Unassigned Burkholderiales (Order)","Unassigned Comamonadaceae (Family)","Unassigned Gammaproteobacteria (Class)","Unassigned Rhodocyclaceae (Family)","Unassigned Sutterellaceae (Family)","Z-35")

other.photo <- c("Chlorobiales","Chloracidobacteriales")

proteo.other <- clone(microeco.dataset.3)
proteo.other$tax_table <- subset(proteo.other$tax_table, Order %in% other.photo | Genus %in% proteobacteria.genus)
check.proteo.other <- proteo.other$tax_table 
proteo.other$tidy_dataset()
check.proteo.other <- proteo.other$otu_table
proteo.other$cal_betadiv(unifrac = TRUE, binary = TRUE)
names(proteo.other$beta_diversity)
proteo.other$beta_diversity

proteo.other.nullmodel <- trans_nullmodel$new(proteo.other, env_cols = 23:31)
saveRDS(proteo.other.nullmodel, "proteo.other.nullmodel.rds")

proteo.other.nullmodel$res_ses_betampd

chemolithoautotroph.nullmodel <- trans_nullmodel$new(chemolithoautotroph, env_cols = 23:31)
saveRDS(chemolithoautotroph.nullmodel, "chemolithoautotroph.nullmodel.rds")

heterotroph.nullmodel <- trans_nullmodel$new(heterotroph, env_cols = 23:31)
saveRDS(heterotroph.nullmodel, "heterotroph.nullmodel.rds")

### Run all the following on HPC and then import into R
### cal mantel's correlation
chloroflexi.nullmodel$cal_mantel_corr()
cyanobacteria.nullmodel$cal_mantel_corr()
proteo.other.nullmodel$cal_mantel_corr()
photosynthetic.nullmodel$cal_mantel_corr()
chemolithoautotroph.nullmodel$cal_mantel_corr()
heterotroph.nullmodel$cal_mantel_corr()
microeco.nullmodel$cal_mantel_corr()

### Using null model :  Calculate standardized effect size of betaMPD i.e. beta net relatedness index (betaNRI).beta mean pairwise phylogenetic distance (betaMPD)
chloroflexi.nullmodel$cal_ses_betampd(runs = 1000, abundance.weighted = TRUE)
cyanobacteria.nullmodel$cal_ses_betampd(runs = 1000, abundance.weighted = TRUE)
proteo.other.nullmodel$cal_ses_betampd(runs = 1000, abundance.weighted = TRUE)
photosynthetic.nullmodel$cal_ses_betampd(runs = 1000, abundance.weighted = TRUE)
chemolithoautotroph.nullmodel$cal_ses_betampd(runs = 1000, abundance.weighted = TRUE)
heterotroph.nullmodel$cal_ses_betampd(runs = 1000, abundance.weighted = TRUE)
microeco.nullmodel$cal_ses_betampd(runs = 1000, abundance.weighted = TRUE)

### Using null model : Calculate standardized effect size of betaMNTD, i.e. beta nearest taxon index (betaNTI). BetaNTI(ses.betamntd) can be used to indicate the phylogenetic terminal turnover
tmp.0 <- "/hpctmp/chrisg25/R_out/TMFZ/postDADA2/tmp.0"; dir.create(tmp.0)
microeco.nullmodel$cal_ses_betamntd(runs = 1000, abundance.weighted = TRUE, use_iCAMP = TRUE, iCAMP_tempdir = tmp.0,nworker = 24)

tmp.9 <- "/hpctmp/chrisg25/R_out/TMFZ/postDADA2/tmp.9"; dir.create(tmp.9)
photosynthetic.nullmodel$cal_ses_betamntd(runs = 1000, abundance.weighted = TRUE, use_iCAMP = TRUE, iCAMP_tempdir = tmp.9, nworker = 24)

tmp.10 <- "/hpctmp/chrisg25/R_out/TMFZ/postDADA2/tmp.10"; dir.create(tmp.10)
chemolithoautotroph.nullmodel$cal_ses_betamntd(runs = 1000, abundance.weighted = TRUE, use_iCAMP = TRUE, iCAMP_tempdir = tmp.10, nworker = 24)

tmp.11 <- "/hpctmp/chrisg25/R_out/TMFZ/postDADA2/tmp.11"; dir.create(tmp.11)
heterotroph.nullmodel$cal_ses_betamntd(runs = 1000, abundance.weighted = TRUE, use_iCAMP = TRUE, iCAMP_tempdir = tmp.11, nworker = 24)

tmp.12 <- "/hpctmp/chrisg25/R_out/TMFZ/postDADA2/tmp.12"; dir.create(tmp.12)
chloroflexi.nullmodel$cal_ses_betamntd(runs = 1000, abundance.weighted = TRUE, use_iCAMP = TRUE, iCAMP_tempdir = tmp.12, nworker = 24)

tmp.13 <- "/hpctmp/chrisg25/R_out/TMFZ/postDADA2/tmp.13"; dir.create(tmp.13)
cyanobacteria.nullmodel$cal_ses_betamntd(runs = 1000, abundance.weighted = TRUE, use_iCAMP = TRUE, iCAMP_tempdir = tmp.13, nworker = 24)

tmp.14 <- "/hpctmp/chrisg25/R_out/TMFZ/postDADA2/tmp.14"; dir.create(tmp.14)
proteo.other.nullmodel$cal_ses_betamntd(runs = 1000, abundance.weighted = TRUE, use_iCAMP = TRUE, iCAMP_tempdir = tmp.14, nworker = 24)

### Using null model : RCbray (Bray-Curtis-based Raup-Crick) can be calculated using function cal_rcbray() to assess whether the compositional turnover was governed primarily by drift (Chase et al. 2011). We applied null model to simulate species distribution by randomly sampling individuals from each species pool with preserving species occurrence frequency and sample species richness (Liu et al. 2017)
chloroflexi.nullmodel$cal_rcbray(runs = 1000)
cyanobacteria.nullmodel$cal_rcbray(runs = 1000)
proteo.other.nullmodel$cal_rcbray(runs = 1000)
photosynthetic.nullmodel$cal_rcbray(runs = 1000)
chemolithoautotroph.nullmodel$cal_rcbray(runs = 1000)
heterotroph.nullmodel$cal_rcbray(runs = 1000)
microeco.nullmodel$cal_rcbray(runs = 1000)


### Using null model : calculate the proportion of the inferred processes on the community assembly as shown in the references (Stegen et al. 2013; Liu et al. 2017).
chloroflexi.nullmodel$cal_process(use_betamntd = TRUE)
cyanobacteria.nullmodel$cal_process(use_betamntd = TRUE)
proteo.other.nullmodel$cal_process(use_betamntd = TRUE)
photosynthetic.nullmodel$cal_process(use_betamntd = TRUE)
chemolithoautotroph.nullmodel$cal_process(use_betamntd = TRUE)
heterotroph.nullmodel$cal_process(use_betamntd = TRUE)
microeco.nullmodel$cal_process(use_betamntd = TRUE)

saveRDS(chloroflexi.nullmodel, "chloroflexi.nullmodel.rds")
saveRDS(cyanobacteria.nullmodel, "cyanobacteria.nullmodel.rds")
saveRDS(proteo.other.nullmodel, "proteo.other.nullmodel.rds")
saveRDS(photosynthetic.nullmodel, "photosynthetic.nullmodel.rds")
saveRDS(chemolithoautotroph.nullmodel, "chemolithoautotroph.nullmodel.rds")
saveRDS(heterotroph.nullmodel, "heterotroph.nullmodel.rds")
saveRDS(microeco.nullmodel, "microeco.nullmodel.rds")

### Import the saved null model objects into R using 'readRDS' before running the following
microeco.nullmodel <- readRDS("microeco.nullmodel.rds")
photosynthetic.nullmodel <- readRDS("photosynthetic.nullmodel.rds")
chloroflexi.nullmodel <- readRDS("chloroflexi.nullmodel.rds")
cyanobacteria.nullmodel <- readRDS("cyanobacteria.nullmodel.rds")
proteo.other.nullmodel <- readRDS("proteo.other.nullmodel.rds")
chemolithoautotroph.nullmodel <- readRDS("chemolithoautotroph.nullmodel.rds")
heterotroph.nullmodel <- readRDS("heterotroph.nullmodel.rds")

#To plot the betaNRI, use plot_group_distance function in trans_beta class
# add betaNRI matrix to beta_diversity list
microeco.dataset.3$beta_diversity[["betaNRI"]] <- microeco.nullmodel$res_ses_betampd

photosynthetic$beta_diversity[["betaNRI"]] <- photosynthetic.nullmodel$res_ses_betampd

chemolithoautotroph$beta_diversity[["betaNRI"]] <- chemolithoautotroph.nullmodel$res_ses_betampd

heterotroph$beta_diversity[["betaNRI"]] <- heterotroph.nullmodel$res_ses_betampd

chloroflexi$beta_diversity[["betaNRI"]] <- chloroflexi.nullmodel$res_ses_betampd

cyanobacteria$beta_diversity[["betaNRI"]] <- cyanobacteria.nullmodel$res_ses_betampd

proteo.other$beta_diversity[["betaNRI"]] <- proteo.other.nullmodel$res_ses_betampd


microeco.nullmodel.2 <- trans_beta$new(dataset = microeco.dataset.3 , group = "Region", measure = "betaNRI")

photosynthetic.nullmodel.2 <- trans_beta$new(dataset = photosynthetic , group = "Region", measure = "betaNRI")

chemolithoautotroph.nullmodel.2 <- trans_beta$new(dataset = chemolithoautotroph , group = "Region", measure = "betaNRI")

heterotroph.nullmodel.2 <- trans_beta$new(dataset = heterotroph , group = "Region", measure = "betaNRI")

chloroflexi.nullmodel.2 <- trans_beta$new(dataset = chloroflexi , group = "Region", measure = "betaNRI")

cyanobacteria.nullmodel.2 <- trans_beta$new(dataset = cyanobacteria , group = "Region", measure = "betaNRI")

proteo.other.nullmodel.2 <- trans_beta$new(dataset = proteo.other , group = "Region", measure = "betaNRI")


# transform the distance for each group
microeco.nullmodel.2$cal_group_distance()
photosynthetic.nullmodel.2$cal_group_distance()
chemolithoautotroph.nullmodel.2$cal_group_distance()
heterotroph.nullmodel.2$cal_group_distance()
chloroflexi.nullmodel.2$cal_group_distance()
cyanobacteria.nullmodel.2$cal_group_distance()
proteo.other.nullmodel.2$cal_group_distance()


# within_group default TRUE; whether transform sample distance within groups, if FALSE,transform sample distance between any two groups.
microeco.nullmodel.2$cal_group_distance_diff(method = "wilcox")
photosynthetic.nullmodel.2$cal_group_distance_diff(method = "wilcox")
chemolithoautotroph.nullmodel.2$cal_group_distance_diff(method = "wilcox")
heterotroph.nullmodel.2$cal_group_distance_diff(method = "wilcox")
chloroflexi.nullmodel.2$cal_group_distance_diff(method = "wilcox")
cyanobacteria.nullmodel.2$cal_group_distance_diff(method = "wilcox")
proteo.other.nullmodel.2$cal_group_distance_diff(method = "wilcox")


# plot the results
betaNRI <- microeco.nullmodel.2$plot_group_distance(boxplot_add = "mean")
betaNRI.photo <- photosynthetic.nullmodel.2$plot_group_distance(boxplot_add = "mean")
betaNRI.hetero <- heterotroph.nullmodel.2$plot_group_distance(boxplot_add = "mean")
betaNRI.chemo <- chemolithoautotroph.nullmodel.2$plot_group_distance(boxplot_add = "mean")
betaNRI.chloroflexi <- chloroflexi.nullmodel.2$plot_group_distance(boxplot_add = "mean")
betaNRI.cyanobacteria <- cyanobacteria.nullmodel.2$plot_group_distance(boxplot_add = "mean")
betaNRI.proteo.other <- proteo.other.nullmodel.2$plot_group_distance(boxplot_add = "mean")

pdf("betaNRI.nullmodel.pdf", width = 15, height =10)
betaNRI + geom_hline(yintercept = -2, linetype = 2) + geom_hline(yintercept = 2, linetype = 2)
dev.off()

pdf("betaNRI.photo.nullmodel.pdf", width = 15, height =10)
betaNRI.photo + geom_hline(yintercept = -2, linetype = 2) + geom_hline(yintercept = 2, linetype = 2)
dev.off()

pdf("betaNRI.chemo.nullmodel.pdf", width = 15, height =10)
betaNRI.chemo + geom_hline(yintercept = -2, linetype = 2) + geom_hline(yintercept = 2, linetype = 2)
dev.off()

pdf("betaNRI.hetero.nullmodel.pdf", width = 15, height =10)
betaNRI.hetero + geom_hline(yintercept = -2, linetype = 2) + geom_hline(yintercept = 2, linetype = 2)
dev.off()

pdf("betaNRI.chloroflexi.nullmodel.pdf", width = 15, height =10)
betaNRI.chloroflexi + geom_hline(yintercept = -2, linetype = 2) + geom_hline(yintercept = 2, linetype = 2)
dev.off()

pdf("betaNRI.cyanobacteria.nullmodel.pdf", width = 15, height =10)
betaNRI.cyanobacteria + geom_hline(yintercept = -2, linetype = 2) + geom_hline(yintercept = 2, linetype = 2)
dev.off()

pdf("betaNRI.proteo.other.nullmodel.pdf", width = 15, height =10)
betaNRI.proteo.other + geom_hline(yintercept = -2, linetype = 2) + geom_hline(yintercept = 2, linetype = 2)
dev.off()


### Get fraction of evolutionary drivers

microeco.nullmodel$res_process
photosynthetic.nullmodel$res_process

# ## NO processes calculated for:
chloroflexi.nullmodel$res_process
cyanobacteria.nullmodel$res_process
proteo.other.nullmodel$res_process
chemolithoautotroph.nullmodel$res_process

heterotroph.nullmodel$res_process

### END OF SECTION ####

###############################################################################################################

