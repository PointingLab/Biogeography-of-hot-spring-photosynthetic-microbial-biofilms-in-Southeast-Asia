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

### Unrarefied
unrarefied.min.3 <- readRDS("/path_to_output_directory/TMFZ/03_tabletax/new.ps.clean.3.rds")

### Rarefied
rarefied.min.3 <- readRDS("/path_to_output_directory/TMFZ/03_tabletax/rarefied.min.3.rds")
rarefied.min.prop.3 <- readRDS("/path_to_output_directory/TMFZ/03_tabletax/rarefied.min.prop.3.rds")

### Final working dataset (Rarefied and filtered)
rarefied.min.int.exclude.1.3 <- readRDS("/path_to_output_directory/TMFZ/03_tabletax/rarefied.min.int.exclude.1.3.rds")
rarefied.min.prop.exclude.1.3 <- readRDS("/path_to_output_directory/TMFZ/03_tabletax/rarefied.min.prop.exclude.1.3.rds")

### Extract dataframes from phyloseqs ####

## unrarefied dataset as integers
unrarefied.min.wide.3 <-  as.matrix(as.data.frame(unrarefied.min.3@otu_table))

## rarefied dataset as integers
rarefied.min.wide.3 <-  as.matrix(as.data.frame(rarefied.min.3@otu_table))

## rarefied dataset as percentages
rarefied.min.prop.wide.3 <- as.matrix(as.data.frame(rarefied.min.prop.3@otu_table))

## wide tables as integers and relative abundances of ASVs > 1%
rarefied.min.wide.exclude.1.3  <- as.data.frame(rarefied.min.int.exclude.1.3@otu_table)

## wide tables as percentages and relative abundances of ASVs > 1%
rarefied.min.prop.wide.exclude.1.3 <- as.data.frame(rarefied.min.prop.exclude.1.3@otu_table)

## Get tax table from any of the final working daatset phyloseqs as both are the same
Tax <- as.matrix(as.data.frame(rarefied.min.int.exclude.1.3@tax_table)) 

## Get metadata. use any phyloseq as all have the same metadata
metadata <- as.matrix(as.data.frame(rarefied.min.int.exclude.1.3@sam_data))
  

##############################################################################################################

### Fig 1A : Map of the transect form North to SOuth of the Thailand-Malaysian-fault-zone (not shown here )
### Fig 1B : Distance decay plot illustrating the strong correlation between geographic and phylogenetic distance for hot spring biofilm communities (N = 395) ####
### Fig S5 : Distance decay plots of phylogenetic versus geographic distance for ecological groups of bacteria (N = 395) ####

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

### Fig. S1. Sampling curves for 16S rRNA gene sequencing, shown for unrarefied and rarefied data (rarefied to the minimum sequencing depth of 62,868 reads per sample) ####

######################################################################################################################

## get otu table as matrix for rarefaction curves

## Function rarecurve draws a rarefaction curve for each row of the input data. The rarefaction curves are evaluated using the interval of step sample sizes, always including 1 and total sample size. 

#If sample is specified, a vertical line is drawn at sample with horizontal lines for the rarefied species richness.
### axis description: sample size - number of raw reads and Species - ASVs

## Ensure extracted otu tables are as matrices
#class(unrarefied.min.wide.3) <- "matrix"
#class(rarefied.min.wide.3) <- "matrix"

# Number of INDIVIDULS per site (?)
raremax <- min(rowSums(rarefied.min.wide.3)) # = 62868; 
raremax

raremax.all <- min(rowSums(unrarefied.min.wide.3)) # = 62868; 
raremax.all


## without subsanpling
rarecurve.rarefied.nosub <- rarecurve(rarefied.min.wide.3, step = 250,
                                      col = "khaki3",
                                      cex = 0.6,label=F,
                                      main = "rarecurve()")



rarecurve.all.nosub <- rarecurve(unrarefied.min.wide.3, step = 250,
                                 col = "khaki3",
                                 cex = 0.6,label=F,
                                 main = "rarecurve.all.nosub")



#################################################################################

### Fig. S3. Taxonomic composition of ASVs clustered at Phylum and Class level (N = 395) ####

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


####################################################################################################################

### Fig. S4. Alpha diversity metrics for the 40 hot spring locations sampled in this study. #####
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

## Plots to check with which variable does alpha diversity vary significantly

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

## By carbonate
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

## either use ANOVA or permutest
anova(betadisp.region)
permutest(betadisp.region, pairwise = TRUE)

## Tukey's test on betadisper
betadisp.region.HSD <- TukeyHSD(betadisp.region)
betadisp.region.HSD

# Adonis test -  permanova
adonis.region <- adonis2(rarefied.unifrac.dist ~ Region, data = sampledf)


####################################################################################################################

### Fig. 2B. Phylogenetic diversity and relative abundance of major bacterial and archaeal lineages in biofilms (grey tree), and pairwise comparisons of lineages that were over-represented (blue branches and nodes) or under-represented (orange branches and nodes) between regions ####
### Fig. S6. Heat tree indicating the phylogenetic composition of photosynthetic biofilm community from Southeast Asian hot springs (N = 395). Branch thickness and node size denote relative abundance. #### 

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


# heat tree comparing different regions (pairwise) ####

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




### Use taxmap2 to add read abundances to legend ######

otu_tax_taxmap2 <- otu_tax_taxmap 

otu_tax_taxmap2$data$tax_abund <- calc_taxon_abund(otu_tax_taxmap2, "otu_counts_taxa",
                                                   cols = metadata$Sample_id)

otu_tax_taxmap2$data$tax_occ <- calc_n_samples(otu_tax_taxmap2, "tax_abund", cols = metadata$Sample_id)


## heat tree os taxonomic diversity ####

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
            output_file = "Species.heatree.supplementary.labels.pdf")




####################################################################################################################

### Fig. 2C. Functional profiling of hot spring photosynthetic biofilms using FAPROTAX ##### 
### Fig. 5C. Influence of abiotic variables on predicted metabolic function of biofilms was estimated using Pearson’s Correlation  ####

### CHECK ####

######################################################################################################################

# create object of trans_func
microeco.faprotax <- trans_func$new(dataset = microeco.dataset.3)

# mapping the taxonomy to the database
# this can recognize prokaryotes or fungi automatically if the names of taxonomic levels are standard.
# for fungi example, see : https://chiliubio.github.io/microeco_tutorial/other-dataset.html#fungi-data

# default database for prokaryotes is FAPROTAX database
microeco.faprotax$cal_spe_func(prok_database = "FAPROTAX")

# return t2$res_spe_func, 1 represent trait exists, 0 represent no or cannot confirmed.

microeco.faprotax$res_spe_func#[1:5, 1:2]

faprotax.func.binary <- as.data.frame(microeco.faprotax$res_spe_func)
colnames(faprotax.func.binary)

#faprotax.func.binary.all <- tibble::rownames_to_column(faprotax.func.binary, "ASVs")

## write ASV feature table to excel file ####
write.xlsx(faprotax.func.binary.all, file = "Taxonomy_faprotax.xlsx",
           sheetName = "faprotax.binary", append = TRUE, row.names = FALSE)

## The percentages of the OTUs having the same trait can reflect the functional redundancy of this function in the community.

# calculate the percentages for communities
# here do not consider the abundance - so based on binary?
microeco.faprotax$cal_spe_func_perc(abundance_weighted = FALSE)

faprotax.func.comm.perc.unwabund <- as.data.frame(microeco.faprotax$res_spe_func_perc)

## From v1.3.0, the trans_spe_func_perc function is implemented to get the long-format table for more flexible manipulation, e.g., filtering and grouping. The return res_spe_func_perc_trans in the object is the table for the following visualization. Note that this step is not necessary as plot_spe_func_perc function can automatically invoke this function if res_spe_func_perc_trans is not found.

microeco.faprotax$trans_spe_func_perc()
#Transformed long format table is stored in object$res_spe_func_perc_trans ...

faprotax.func.comm.perc.unwabund.long  <- as.data.frame(microeco.faprotax$res_spe_func_perc_trans)

## as it is byy smaples - 395- not possible to see clearly
pdf("faprotax.func.comm.pdf", width = 20)
microeco.faprotax$plot_spe_func_perc()
dev.off()






microeco.faprotax.2 <- trans_env$new(dataset = microeco.dataset.3, add_data = metadata.micronet.select.2)

microeco.faprotax.2$cal_cor(add_abund_table = microeco.faprotax$res_spe_func_perc, cor_method = "pearson")

pdf("faprotax.func.env.pdf")
microeco.faprotax.2$plot_cor(pheatmap = TRUE)
dev.off()


####################################################################################3

### Fig. 3A Relative abundance of the 25 most abundant ASVs in each sample (N = 395) shown clustered at genus level (bar plots). ####
### Fig. S7 Relative abundance of the 10 most abundant ASVs in ecological groups for each sample (N = 395). ####

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

### Photosynthetic ####
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

### Fig. 4A. Putative biotic interactions within (Zi) and among (Pi) modules for top 10 taxa by Class, for each region #####
### Fig. 4B. Chord diagram illustrating major interactions between abundant ASVs clustered at genus level ####
### Fig. S8. Putative biotic interactions within (Zi) and among (Pi) modules for all taxa by Class, for each region #####

### CHECK #####

######################################################################################

### Network analyses with SpiecEASI and NETCOMI ####

### Create taxtable with additional columns

net_spieceasi_taxtable <- as.data.frame(com.rarefied.min.int.exclude.1.3.rooted@tax_table) %>% rownames_to_column(var = "ASV")  %>% as_tibble()


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


## Spieceasi network plots in netcomi by regions ####

### Create phyloseq with taxtable containing functional groups in taxtable and ASVnos
add.taxtable <- net_spieceasi_taxtable #%>% column_to_rownames(var="ASV")
add.taxtable$ASVno <- add.taxtable$ASV
# netcomi will not plot plain numbers so leave the ASV in 
#add.taxtable <-  add.taxtable %>% mutate(ASVno = as.character(gsub("ASV", "", ASVno)))
add.taxtable <- add.taxtable %>% column_to_rownames(var="ASV")
add.taxtable <- as.matrix(add.taxtable)

com.rarefied.min.int.exclude.1.3.rooted.spnet <- com.rarefied.min.int.exclude.1.3.rooted
#`tax_table<-`(test, test2)
tax_table(com.rarefied.min.int.exclude.1.3.rooted.spnet) <- add.taxtable

test <- as.data.frame(com.rarefied.min.int.exclude.1.3.rooted.spnet@tax_table)

saveRDS(com.rarefied.min.int.exclude.1.3.rooted.spnet,"com.rarefied.min.int.exclude.1.3.rooted.spnet.rds")

## filter out ASVs with "Unassigned Bacteria (Kingdom) to run networks"

select.rarefied.min.int.exclude.1.3.rooted.spnet <- com.rarefied.min.int.exclude.1.3.rooted.spnet %>%
  subset_taxa(!Phylum == "Unassigned Bacteria (Kingdom)")


##phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 556 taxa and 395 samples ]
#sample_data() Sample Data:       [ 395 samples by 33 sample variables ]
#tax_table()   Taxonomy Table:    [ 556 taxa by 10 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 556 tips and 555 internal nodes ]
#refseq()      DNAStringSet:      [ 556 reference sequences ]

## 16 ASVs removed - Unassigned Bacteria (Kingdom)

test <- as.data.frame(as.matrix(select.rarefied.min.int.exclude.1.3.rooted.spnet@tax_table))

test.1 <- as.data.frame(as.matrix(com.rarefied.min.int.exclude.1.3.rooted.spnet@tax_table))

saveRDS(select.rarefied.min.int.exclude.1.3.rooted.spnet,"select.rarefied.min.int.exclude.1.3.rooted.spnet.rds")


### Get phyloseqs for each of the 6 regions 

## cluster 1 :North.Thailand <- c('PT','PP','LN','HS','MK') #####

NT.loc <- c("PT","PP","LN","HS","MK")

### Get phyloseqs for each of the 6 regions without ASVs with Unassigned Bacteria and use these

## phyloseqs without Unassigned Bacteria (Kingdom)
NT.select.rarefied.min.int.exclude.1.3.rooted.spnet <- select.rarefied.min.int.exclude.1.3.rooted.spnet %>%
  subset_samples(Location.Code %in% NT.loc)

NT.select.rarefied.min.int.exclude.1.3.rooted.spnet <- 
  prune_taxa(taxa_sums(NT.select.rarefied.min.int.exclude.1.3.rooted.spnet) > 0, NT.select.rarefied.min.int.exclude.1.3.rooted.spnet)

any(taxa_sums(NT.select.rarefied.min.int.exclude.1.3.rooted.spnet) == 0)

test.1 <- as.data.frame(NT.select.rarefied.min.int.exclude.1.3.rooted.spnet@tax_table)

test.2 <- as.data.frame(NT.select.rarefied.min.int.exclude.1.3.rooted.spnet@sam_data)

## Cluster 2 : North.Malaysia <- c('AP','BA','US') ######

NM.loc <- c("AP","BA","US")

## phyloseqs without Unassigned Bacteria (Kingdom)
NM.select.rarefied.min.int.exclude.1.3.rooted.spnet <- select.rarefied.min.int.exclude.1.3.rooted.spnet %>%
  subset_samples(Location.Code %in% NM.loc)

NM.select.rarefied.min.int.exclude.1.3.rooted.spnet <- 
  prune_taxa(taxa_sums(NM.select.rarefied.min.int.exclude.1.3.rooted.spnet) > 0, NM.select.rarefied.min.int.exclude.1.3.rooted.spnet)

any(taxa_sums(NM.select.rarefied.min.int.exclude.1.3.rooted.spnet) == 0)

test.1 <- as.data.frame(NM.select.rarefied.min.int.exclude.1.3.rooted.spnet@tax_table)

test.2 <- as.data.frame(NM.select.rarefied.min.int.exclude.1.3.rooted.spnet@sam_data)

## Cluster 3 : South.Malaysia <- c('KJ','SE','LA') ######

SM.loc <- c("SE","LA","KJ") 

# phyloseqs without Unassigned Bacteria (Kingdom)
SM.select.rarefied.min.int.exclude.1.3.rooted.spnet <- select.rarefied.min.int.exclude.1.3.rooted.spnet %>%
  subset_samples(Location.Code %in% SM.loc)

SM.select.rarefied.min.int.exclude.1.3.rooted.spnet <- 
  prune_taxa(taxa_sums(SM.select.rarefied.min.int.exclude.1.3.rooted.spnet) > 0, SM.select.rarefied.min.int.exclude.1.3.rooted.spnet)

any(taxa_sums(SM.select.rarefied.min.int.exclude.1.3.rooted.spnet) == 0)

test.1 <- as.data.frame(SM.select.rarefied.min.int.exclude.1.3.rooted.spnet@tax_table)

test.2 <- as.data.frame(SM.select.rarefied.min.int.exclude.1.3.rooted.spnet@sam_data)

## Cluster 4 : All Sembawang : Singapore <- 'SW' #####

# phyloseqs without Unassigned Bacteria (Kingdom)
SG.select.rarefied.min.int.exclude.1.3.rooted.spnet <- select.rarefied.min.int.exclude.1.3.rooted.spnet %>%
  subset_samples(Location.Code == "SW")

SG.select.rarefied.min.int.exclude.1.3.rooted.spnet <- 
  prune_taxa(taxa_sums(SG.select.rarefied.min.int.exclude.1.3.rooted.spnet) > 0, SG.select.rarefied.min.int.exclude.1.3.rooted.spnet)

any(taxa_sums(SG.select.rarefied.min.int.exclude.1.3.rooted.spnet) == 0)

test.1 <- as.data.frame(SG.select.rarefied.min.int.exclude.1.3.rooted.spnet@tax_table)

test.2 <- as.data.frame(SG.select.rarefied.min.int.exclude.1.3.rooted.spnet@sam_data)


## Cluster 5 : South.Thailand <- c('RB','RN') #####

ST <- c("RN","RB")

# phyloseqs without Unassigned Bacteria (Kingdom)
ST.select.rarefied.min.int.exclude.1.3.rooted.spnet <- select.rarefied.min.int.exclude.1.3.rooted.spnet %>%
  subset_samples(Location.Code %in% ST)

ST.select.rarefied.min.int.exclude.1.3.rooted.spnet <- 
  prune_taxa(taxa_sums(ST.select.rarefied.min.int.exclude.1.3.rooted.spnet) > 0, ST.select.rarefied.min.int.exclude.1.3.rooted.spnet)

any(taxa_sums(ST.select.rarefied.min.int.exclude.1.3.rooted.spnet) == 0)

test.1 <- as.data.frame(ST.select.rarefied.min.int.exclude.1.3.rooted.spnet@tax_table)

test.2 <- as.data.frame(ST.select.rarefied.min.int.exclude.1.3.rooted.spnet@sam_data)


## Cluster 6 : Central.Thailand <- 'PB' ######

# phyloseqs without Unassigned Bacteria (Kingdom)
CT.select.rarefied.min.int.exclude.1.3.rooted.spnet <- select.rarefied.min.int.exclude.1.3.rooted.spnet %>%
  subset_samples(Location.Code == "PB")

CT.select.rarefied.min.int.exclude.1.3.rooted.spnet <- 
  prune_taxa(taxa_sums(CT.select.rarefied.min.int.exclude.1.3.rooted.spnet) > 0, CT.select.rarefied.min.int.exclude.1.3.rooted.spnet)

any(taxa_sums(CT.select.rarefied.min.int.exclude.1.3.rooted.spnet) == 0)

test.1 <- as.data.frame(CT.select.rarefied.min.int.exclude.1.3.rooted.spnet@tax_table)

test.2 <- as.data.frame(CT.select.rarefied.min.int.exclude.1.3.rooted.spnet@sam_data)


### save cluster phyloseqs wihtout Unassigned Bacteria (Kingdom)
saveRDS(NT.select.rarefied.min.int.exclude.1.3.rooted.spnet, "NT.select.rarefied.min.int.exclude.1.3.rooted.spnet.rds")
saveRDS(NM.select.rarefied.min.int.exclude.1.3.rooted.spnet, "NM.select.rarefied.min.int.exclude.1.3.rooted.spnet.rds")
saveRDS(ST.select.rarefied.min.int.exclude.1.3.rooted.spnet, "ST.select.rarefied.min.int.exclude.1.3.rooted.spnet.rds")
saveRDS(SM.select.rarefied.min.int.exclude.1.3.rooted.spnet, "SM.select.rarefied.min.int.exclude.1.3.rooted.spnet.rds")
saveRDS(CT.select.rarefied.min.int.exclude.1.3.rooted.spnet, "CT.select.rarefied.min.int.exclude.1.3.rooted.spnet.rds")
saveRDS(SG.select.rarefied.min.int.exclude.1.3.rooted.spnet, "SG.select.rarefied.min.int.exclude.1.3.rooted.spnet.rds")


## Run on server ###

## NT
# Network construction
#NT_net_spieceasi_mb <- netConstruct(NT.rarefied.min.int.exclude.1.3.rooted.spnet, 
#                     measure = "spieceasi",measurePar = list(method = "mb",
#                     sel.criterion='stars', lambda.min.ratio=1e-2, nlambda=20,
##                      pulsar.params=list(rep.num=100),
#                     symBetaMode = "ave"), sparsMethod = "none",
#                      normMethod = "none",
#                     verbose = 3)

# Network analysis (see ?netAnalyze for details)
#NT_anlayze_spieceasi_mb <- netAnalyze(NT_net_spieceasi_mb, clustMethod = "cluster_fast_greedy", hubPar = "eigenvector", weightDeg = FALSE, normDeg = FALSE)
#

## NM
# Network construction
#NM_net_spieceasi_mb <- netConstruct(NM.rarefied.min.int.exclude.1.3.rooted.spnet, 
#                     measure = "spieceasi",measurePar = list(method = "mb",
#                     sel.criterion='stars', lambda.min.ratio=1e-2, nlambda=20,
#                      pulsar.params=list(rep.num=100),
#                     symBetaMode = "ave"), sparsMethod = "none",
#                      normMethod = "none",
#                     verbose = 3)

# Network analysis (see ?netAnalyze for details)
#NM_anlayze_spieceasi_mb <- netAnalyze(NM_net_spieceasi_mb, clustMethod = "cluster_fast_greedy", hubPar = "eigenvector", weightDeg = FALSE, normDeg = FALSE)
#

## SM
# Network construction
#SM_net_spieceasi_mb <- netConstruct(SM.rarefied.min.int.exclude.1.3.rooted.spnet, 
#                     measure = "spieceasi",measurePar = list(method = "mb",
#                     sel.criterion='stars', lambda.min.ratio=1e-2, nlambda=20,
#                      pulsar.params=list(rep.num=100),
#                     symBetaMode = "ave"), sparsMethod = "none",
#                      normMethod = "none",
#                    verbose = 3)

# Network analysis (see ?netAnalyze for details)
#SM_anlayze_spieceasi_mb <- netAnalyze(SM_net_spieceasi_mb, clustMethod = "cluster_fast_greedy", hubPar = "eigenvector", weightDeg = FALSE, normDeg = FALSE)


## SG
# Network construction
#SG_net_spieceasi_mb <- netConstruct(SG.rarefied.min.int.exclude.1.3.rooted.spnet, 
#                     measure = "spieceasi",measurePar = list(method = "mb",
#                     sel.criterion='stars', lambda.min.ratio=1e-2, nlambda=20,
#                      pulsar.params=list(rep.num=100),
#                     symBetaMode = "ave"), sparsMethod = "none",
#                      normMethod = "none",
#                     verbose = 3)

# Network analysis (see ?netAnalyze for details)
#SG_anlayze_spieceasi_mb <- netAnalyze(SG_net_spieceasi_mb, clustMethod = "cluster_fast_greedy", hubPar = "eigenvector", weightDeg = FALSE, normDeg = FALSE)


## ST
# Network construction
#ST_net_spieceasi_mb <- netConstruct(ST.rarefied.min.int.exclude.1.3.rooted.spnet, 
#                     measure = "spieceasi",measurePar = list(method = "mb",
#                     sel.criterion='stars', lambda.min.ratio=1e-2, nlambda=20,
#                      pulsar.params=list(rep.num=100),
#                     symBetaMode = "ave"), sparsMethod = "none",
#                      normMethod = "none",
#                     verbose = 3)

# Network analysis (see ?netAnalyze for details)
#ST_anlayze_spieceasi_mb <- netAnalyze(ST_net_spieceasi_mb, clustMethod = "cluster_fast_greedy", hubPar = "eigenvector", weightDeg = FALSE, normDeg = FALSE)


## CT
# Network construction
#CT_net_spieceasi_mb <- netConstruct(CT.rarefied.min.int.exclude.1.3.rooted.spnet, 
#                     measure = "spieceasi",measurePar = list(method = "mb",
#                     sel.criterion='stars', lambda.min.ratio=1e-2, nlambda=20,
#                      pulsar.params=list(rep.num=100),
#                     symBetaMode = "ave"), sparsMethod = "none",
#                      normMethod = "none",
#                     verbose = 3)


### Chord Diagram ####

select.rarefied.min.int.exclude.1.3.rooted.spnet <- readRDS("select.rarefied.min.int.exclude.1.3.rooted.spnet.rds")

## create microeco dataset wihtout "Unassigned Bacteria (Kingdom)" intaxa
select.microeco.dataset <- clone(microeco.dataset)

select.microeco.dataset$tax_table <- select.microeco.dataset$tax_table[select.microeco.dataset$tax_table$Phylum != "Unassigned Bacteria (Kingdom)", ]

## 556 ASVs
check.select <- select.microeco.dataset$tax_table 

check.select.phyloseq <- as.data.frame(select.rarefied.min.int.exclude.1.3.rooted.spnet@tax_table)

select.microeco.dataset$tidy_dataset()

check.select <- select.microeco.dataset$otu_table

check.select.phyloseq <- as.data.frame(select.rarefied.min.int.exclude.1.3.rooted.spnet@otu_table)

### Get microeco datasets for the 6 regions to use in chord diagram (excluding Unassigned bacteria at phykum level) and also for mantel's (all taxa)

## cluster 1 :North.Thailand <- c('PT','PP','LN','HS','MK') #####
NT.loc <- c("PT","PP","LN","HS","MK")

## convert phyloseq to microeco object - DO NOT USE as it changes tax_table
#NT.mircoeco.dataset <- phyloseq2meco(NT.rarefied.min.int.exclude.1.3.rooted.spnet)

## Instead subset based on samples

NT.microeco.dataset <- clone(microeco.dataset)

NT.microeco.dataset$sample_table <- subset(NT.microeco.dataset$sample_table, Region == "North.Thailand")

NT.microeco.dataset$tidy_dataset()

test <- as.data.frame(NT.rarefied.min.int.exclude.1.3.rooted.spnet@otu_table)
test.1 <- NT.microeco.dataset$otu_table


NT.select.microeco.dataset <- clone(select.microeco.dataset)

NT.select.microeco.dataset$sample_table <- subset(NT.microeco.dataset$sample_table, Region == "North.Thailand")

NT.select.microeco.dataset$tidy_dataset()

test <- as.data.frame(NT.select.rarefied.min.int.exclude.1.3.rooted.spnet@otu_table)
test.1 <- NT.select.microeco.dataset$otu_table

## Cluster 2 : North.Malaysia <- c('AP','BA','US') ######

NM.loc <- c("AP","BA","US")

NM.microeco.dataset <- clone(microeco.dataset)

NM.microeco.dataset$sample_table <- subset(NM.microeco.dataset$sample_table, Region == "North.Malaysia")

NM.microeco.dataset$tidy_dataset()

test <- as.data.frame(NM.rarefied.min.int.exclude.1.3.rooted.spnet@otu_table)
test.1 <- NM.microeco.dataset$otu_table


NM.select.microeco.dataset <- clone(select.microeco.dataset)

NM.select.microeco.dataset$sample_table <- subset(NM.microeco.dataset$sample_table, Region == "North.Malaysia")

NM.select.microeco.dataset$tidy_dataset()

test <- as.data.frame(NM.select.rarefied.min.int.exclude.1.3.rooted.spnet@otu_table)
test.1 <- NM.select.microeco.dataset$otu_table


## Cluster 3 : South.Malaysia <- c('KJ','SE','LA') ######

SM.loc <- c("SE","LA","KJ") 

SM.microeco.dataset <- clone(microeco.dataset)

SM.microeco.dataset$sample_table <- subset(SM.microeco.dataset$sample_table, Region == "South.Malaysia")

SM.microeco.dataset$tidy_dataset()

test <- as.data.frame(SM.rarefied.min.int.exclude.1.3.rooted.spnet@otu_table)
test.1 <- SM.microeco.dataset$otu_table


SM.select.microeco.dataset <- clone(select.microeco.dataset)

SM.select.microeco.dataset$sample_table <- subset(SM.microeco.dataset$sample_table, Region == "South.Malaysia")

SM.select.microeco.dataset$tidy_dataset()

test <- as.data.frame(SM.select.rarefied.min.int.exclude.1.3.rooted.spnet@otu_table)
test.1 <- SM.select.microeco.dataset$otu_table


## Cluster 4 : All Sembawang : Singapore <- 'SW' #####

SG.microeco.dataset <- clone(microeco.dataset)

SG.microeco.dataset$sample_table <- subset(SG.microeco.dataset$sample_table, Region == "Singapore")

SG.microeco.dataset$tidy_dataset()

test <- as.data.frame(SG.rarefied.min.int.exclude.1.3.rooted.spnet@otu_table)
test.1 <- SG.microeco.dataset$otu_table


SG.select.microeco.dataset <- clone(select.microeco.dataset)

SG.select.microeco.dataset$sample_table <- subset(SG.microeco.dataset$sample_table, Region == "Singapore")

SG.select.microeco.dataset$tidy_dataset()

test <- as.data.frame(SG.select.rarefied.min.int.exclude.1.3.rooted.spnet@otu_table)
test.1 <- SG.select.microeco.dataset$otu_table



## Cluster 5 : South.Thailand <- c('RB','RN') #####

ST <- c("RN","RB")

ST.microeco.dataset <- clone(microeco.dataset)

ST.microeco.dataset$sample_table <- subset(ST.microeco.dataset$sample_table, Region == "South.Thailand")

ST.microeco.dataset$tidy_dataset()

test <- as.data.frame(ST.rarefied.min.int.exclude.1.3.rooted.spnet@otu_table)
test.1 <- ST.microeco.dataset$otu_table


ST.select.microeco.dataset <- clone(select.microeco.dataset)

ST.select.microeco.dataset$sample_table <- subset(ST.microeco.dataset$sample_table, Region == "South.Thailand")

ST.select.microeco.dataset$tidy_dataset()

test <- as.data.frame(ST.select.rarefied.min.int.exclude.1.3.rooted.spnet@otu_table)
test.1 <- ST.select.microeco.dataset$otu_table


## Cluster 6 : Central.Thailand <- 'PB' ######

CT.microeco.dataset <- clone(microeco.dataset)

CT.microeco.dataset$sample_table <- subset(CT.microeco.dataset$sample_table, Region == "Central.Thailand")

CT.microeco.dataset$tidy_dataset()

test <- as.data.frame(CT.rarefied.min.int.exclude.1.3.rooted.spnet@otu_table)
test.1 <- CT.microeco.dataset$otu_table


CT.select.microeco.dataset <- clone(select.microeco.dataset)

CT.select.microeco.dataset$sample_table <- subset(CT.microeco.dataset$sample_table, Region == "Central.Thailand")

CT.select.microeco.dataset$tidy_dataset()

test <- as.data.frame(CT.select.rarefied.min.int.exclude.1.3.rooted.spnet@otu_table)
test.1 <- CT.select.microeco.dataset$otu_table


### save regions microeco for mantel's calculation in Rscript 5

saveRDS(CT.microeco.dataset,"CT.microeco.dataset.rds")
saveRDS(NT.microeco.dataset,"NT.microeco.dataset.rds")
saveRDS(ST.microeco.dataset,"ST.microeco.dataset.rds")
saveRDS(NM.microeco.dataset,"NM.microeco.dataset.rds")
saveRDS(SM.microeco.dataset,"SM.microeco.dataset.rds")
saveRDS(SG.microeco.dataset,"SG.microeco.dataset.rds")

saveRDS(CT.select.microeco.dataset,"CT.select.microeco.dataset.rds")
saveRDS(NT.select.microeco.dataset,"NT.select.microeco.dataset.rds")
saveRDS(ST.select.microeco.dataset,"ST.select.microeco.dataset.rds")
saveRDS(NM.select.microeco.dataset,"NM.select.microeco.dataset.rds")
saveRDS(SM.select.microeco.dataset,"SM.select.microeco.dataset.rds")
saveRDS(SG.select.microeco.dataset,"SG.select.microeco.dataset.rds")


## From microeco developer: network constructed with netcomi using spieceasi and create chord diagram 


### START here ####
# Compute layout
tmp_ASV_net_spieceasi_mb <- ASV_net_spieceasi_mb$adjaMat1
#tmp_NT_net_spieceasi_mb  <- NT_net_spieceasi_mb$adjaMat1
#tmp_NM_net_spieceasi_mb  <- NM_net_spieceasi_mb$adjaMat1
#tmp_ST_net_spieceasi_mb  <- ST_net_spieceasi_mb$adjaMat1
#tmp_SM_net_spieceasi_mb  <- SM_net_spieceasi_mb$adjaMat1
#tmp_CT_net_spieceasi_mb  <- CT_net_spieceasi_mb$adjaMat1
#tmp_SG_net_spieceasi_mb  <- SG_net_spieceasi_mb$adjaMat1

tmp_NT_select_net_spieceasi_mb  <- NT_select_net_spieceasi_mb$adjaMat1
tmp_NM_select_net_spieceasi_mb  <- NM_select_net_spieceasi_mb$adjaMat1
tmp_ST_select_net_spieceasi_mb  <- ST_select_net_spieceasi_mb$adjaMat1
tmp_SM_select_net_spieceasi_mb  <- SM_select_net_spieceasi_mb$adjaMat1
tmp_CT_select_net_spieceasi_mb  <- CT_select_net_spieceasi_mb$adjaMat1
tmp_SG_select_net_spieceasi_mb  <- SG_select_net_spieceasi_mb$adjaMat1


# 1 of diagonal should be removed
diag(tmp_ASV_net_spieceasi_mb) <- 0
#diag(tmp_NT_net_spieceasi_mb) <- 0
#diag(tmp_NM_net_spieceasi_mb) <- 0
#diag(tmp_ST_net_spieceasi_mb) <- 0
#diag(tmp_SM_net_spieceasi_mb) <- 0
#diag(tmp_CT_net_spieceasi_mb) <- 0
#diag(tmp_SG_net_spieceasi_mb) <- 0

diag(tmp_NT_select_net_spieceasi_mb) <- 0
diag(tmp_NM_select_net_spieceasi_mb) <- 0
diag(tmp_ST_select_net_spieceasi_mb) <- 0
diag(tmp_SM_select_net_spieceasi_mb) <- 0
diag(tmp_CT_select_net_spieceasi_mb) <- 0
diag(tmp_SG_select_net_spieceasi_mb) <- 0

## extract adjacency matrix

graph3.ASV.net.spiec.mb <- igraph::graph_from_adjacency_matrix(tmp_ASV_net_spieceasi_mb, weighted = TRUE)

graph3.NT.select.net.spiec.mb <- igraph::graph_from_adjacency_matrix(tmp_NT_select_net_spieceasi_mb, weighted = TRUE)

graph3.NM.select.net.spiec.mb <- igraph::graph_from_adjacency_matrix(tmp_NM_select_net_spieceasi_mb, weighted = TRUE)

graph3.ST.select.net.spiec.mb <- igraph::graph_from_adjacency_matrix(tmp_ST_select_net_spieceasi_mb, weighted = TRUE)

graph3.SM.select.net.spiec.mb <- igraph::graph_from_adjacency_matrix(tmp_SM_select_net_spieceasi_mb, weighted = TRUE)

graph3.CT.select.net.spiec.mb <- igraph::graph_from_adjacency_matrix(tmp_CT_select_net_spieceasi_mb, weighted = TRUE)

graph3.SG.select.net.spiec.mb <- igraph::graph_from_adjacency_matrix(tmp_SG_select_net_spieceasi_mb, weighted = TRUE)


# must assign names:
V(graph3.ASV.net.spiec.mb)$name <- as.character(as_ids(V(graph3.ASV.net.spiec.mb)))
V(graph3.NT.select.net.spiec.mb)$name <- as.character(as_ids(V(graph3.NT.select.net.spiec.mb)))
V(graph3.NM.select.net.spiec.mb)$name <- as.character(as_ids(V(graph3.NM.select.net.spiec.mb)))
V(graph3.ST.select.net.spiec.mb)$name <- as.character(as_ids(V(graph3.ST.select.net.spiec.mb)))
V(graph3.SM.select.net.spiec.mb)$name <- as.character(as_ids(V(graph3.SM.select.net.spiec.mb)))
V(graph3.CT.select.net.spiec.mb)$name <- as.character(as_ids(V(graph3.CT.select.net.spiec.mb)))
V(graph3.SG.select.net.spiec.mb)$name <- as.character(as_ids(V(graph3.SG.select.net.spiec.mb)))


E(graph3.ASV.net.spiec.mb)$label <- ifelse(E(graph3.ASV.net.spiec.mb)$weight > 0, '+', '-')
E(graph3.NT.select.net.spiec.mb)$label <- ifelse(E(graph3.NT.select.net.spiec.mb)$weight > 0, '+', '-')
E(graph3.NM.select.net.spiec.mb)$label <- ifelse(E(graph3.NM.select.net.spiec.mb)$weight > 0, '+', '-')
E(graph3.ST.select.net.spiec.mb)$label <- ifelse(E(graph3.ST.select.net.spiec.mb)$weight > 0, '+', '-')
E(graph3.SM.select.net.spiec.mb)$label <- ifelse(E(graph3.SM.select.net.spiec.mb)$weight > 0, '+', '-')
E(graph3.CT.select.net.spiec.mb)$label <- ifelse(E(graph3.CT.select.net.spiec.mb)$weight > 0, '+', '-')
E(graph3.SG.select.net.spiec.mb)$label <- ifelse(E(graph3.SG.select.net.spiec.mb)$weight > 0, '+', '-')



microeco.network <- trans_network$new(dataset = microeco.dataset,
                                      cor_method = NULL,taxa_level = "OTU")

NT.microeco.network <- trans_network$new(dataset = NT.select.microeco.dataset,
                                         cor_method = NULL,taxa_level = "OTU")

NM.microeco.network <- trans_network$new(dataset = NM.select.microeco.dataset,
                                         cor_method = NULL,taxa_level = "OTU")

ST.microeco.network <- trans_network$new(dataset = ST.select.microeco.dataset,
                                         cor_method = NULL,taxa_level = "OTU")

SM.microeco.network <- trans_network$new(dataset = SM.select.microeco.dataset,
                                         cor_method = NULL,taxa_level = "OTU")

CT.microeco.network <- trans_network$new(dataset = CT.select.microeco.dataset,
                                         cor_method = NULL,taxa_level = "OTU")

SG.microeco.network <- trans_network$new(dataset = SG.select.microeco.dataset,
                                         cor_method = NULL,taxa_level = "OTU")


microeco.network$cal_network(network_method = NULL)
NT.microeco.network$cal_network(network_method = NULL)
NM.microeco.network$cal_network(network_method = NULL)
ST.microeco.network$cal_network(network_method = NULL)
SM.microeco.network$cal_network(network_method = NULL)
CT.microeco.network$cal_network(network_method = NULL)
SG.microeco.network$cal_network(network_method = NULL)

test <- NT.microeco.dataset$tax_table
# tax_table is necessary for the summary when you use cal_sum_links function. Please manually assign it
microeco.network$tax_table <- microeco.dataset$tax_table
NT.microeco.network$tax_table <- NT.select.microeco.dataset$tax_table
NM.microeco.network$tax_table <- NM.select.microeco.dataset$tax_table
ST.microeco.network$tax_table <- ST.select.microeco.dataset$tax_table
SM.microeco.network$tax_table <- SM.select.microeco.dataset$tax_table
CT.microeco.network$tax_table <- CT.select.microeco.dataset$tax_table
SG.microeco.network$tax_table <- SG.select.microeco.dataset$tax_table

test <- NT.select.microeco.dataset$tax_table

# g1 to your trans_network object
microeco.network$res_network <- graph3.ASV.select.net.spiec.mb
NT.microeco.network$res_network <- graph3.NT.select.net.spiec.mb
NM.microeco.network$res_network <- graph3.NM.select.net.spiec.mb
ST.microeco.network$res_network <- graph3.ST.select.net.spiec.mb
SM.microeco.network$res_network <- graph3.SM.select.net.spiec.mb
CT.microeco.network$res_network <- graph3.CT.select.net.spiec.mb
SG.microeco.network$res_network <- graph3.SG.select.net.spiec.mb


## The default method "cluster_fast_greedy" can not be applied to directed network! Automatically switch to method "cluster_walktrap" ...Invoke cluster_walktrap function to find densely connected subgraphs ... Modules are assigned in network with attribute name -- module ...

## 4 modules
microeco.network$cal_module(method = "cluster_fast_greedy",
                            module_name_prefix = "M")

## 16 modules
NT.microeco.network$cal_module(method = "cluster_fast_greedy",
                               module_name_prefix = "M")

## 14 modules
NM.microeco.network$cal_module(method = "cluster_fast_greedy",
                               module_name_prefix = "M")

## 11 modules
ST.microeco.network$cal_module(method = "cluster_fast_greedy",
                               module_name_prefix = "M")

## 10 modules
SM.microeco.network$cal_module(method = "cluster_fast_greedy",
                               module_name_prefix = "M")

## 10 modules
CT.microeco.network$cal_module(method = "cluster_fast_greedy",
                               module_name_prefix = "M")

## 15 modules
SG.microeco.network$cal_module(method = "cluster_fast_greedy",
                               module_name_prefix = "M")


microeco.network$cal_network_attr()
#microeco.network$res_network_attr

NT.microeco.network$cal_network_attr()
#NT.microeco.network$res_network_attr

NM.microeco.network$cal_network_attr()
#NT.microeco.network$res_network_attr

ST.microeco.network$cal_network_attr()
#NT.microeco.network$res_network_attr

SM.microeco.network$cal_network_attr()
#NT.microeco.network$res_network_attr

CT.microeco.network$cal_network_attr()
#NT.microeco.network$res_network_attr

SG.microeco.network$cal_network_attr()
#NT.microeco.network$res_network_attr


microeco.network$get_node_table(node_roles = TRUE)
#microeco.network$res_node_table
# default TRUE; whether calculate node roles, i.e. Module hubs, Network hubs, Connectors and Peripheral

NT.microeco.network$get_node_table(node_roles = TRUE)
#NT.microeco.network$res_node_table 

NM.microeco.network$get_node_table(node_roles = TRUE)
#NT.microeco.network$res_node_table

ST.microeco.network$get_node_table(node_roles = TRUE)
#NT.microeco.network$res_node_table

CT.microeco.network$get_node_table(node_roles = TRUE)
#NT.microeco.network$res_node_table

SG.microeco.network$get_node_table(node_roles = TRUE)
#NT.microeco.network$res_node_table


### edges found in new microeco version

microeco.network$get_edge_table()
#microeco.network$res_edge_table 

NT.microeco.network$get_edge_table()
#NT.microeco.network$res_edge_table 

NM.microeco.network$get_edge_table()
ST.microeco.network$get_edge_table()
SM.microeco.network$get_edge_table()
CT.microeco.network$get_edge_table()
SG.microeco.network$get_edge_table()


microeco.network$get_adjacency_matrix()
#microeco.network$res_adjacency_matrix

NT.microeco.network$get_adjacency_matrix()
#NT.microeco.network$res_adjacency_matrix

NM.microeco.network$get_adjacency_matrix()
ST.microeco.network$get_adjacency_matrix()
SM.microeco.network$get_adjacency_matrix()
CT.microeco.network$get_adjacency_matrix()
SG.microeco.network$get_adjacency_matrix()


#plot the node classification in terms of the within-module connectivity and among-module connectivity.Plot the classification and importance of nodes

# add_label = TRUE can be used to directly add text label for points

#microeco.network$plot_taxa_roles(use_type = 1, add_label = TRUE)
#NT.microeco.network$plot_taxa_roles(use_type = 1, add_label = TRUE)

# plot node roles :  ####
# save as png 750, 600

## Get the node property table. The properties may include the node names, modules allocation, degree, betweenness, abundance, taxonomy, within-module connectivity and among-module connectivity <doi:10.1016/j.geoderma.2022.115866>.

color_values = RColorBrewer::brewer.pal(12, "Paired")

color_values_1 = c("Cyanobacteriia"="#3B4D16","Bacteroidia"="#950D3D", "Chloroflexia"="#848F22", "Blastocatellia"="#A7A3FE","Alphaproteobacteria" ="#B2C9B2",  "Unassigned Bacteria (Kingdom)" = "#8CCBD3", "Deinococci" = "#0DD7FD", "Planctomycetes"= "#99D584","Gammaproteobacteria" = "#0D8CA8","Anaerolineae" = "#DFC500")

color_values_2 = c("#3B4D16","#950D3D", "#848F22", "#A7A3FE","#B2C9B2", "#8CCBD3", "#0DD7FD", "#99D584", "#0D8CA8","#DFC500","#2A7F72","#FEAB2E", "#4B00FD", "#26E0A9" ,"#CC003D")


# top value is the p-value and bottom is the z-score
#Abundance expressed as a percentage; betweenness_centrality: betweenness centrality; betweenness_centrality: closeness centrality; eigenvector_centrality: eigenvector centrality; z: within-module connectivity; p: among-module connectivity.
microeco.network.node.table <- microeco.network$res_node_table

## plot at class level for plot_taxa_roles
## show_number : sorting according to the nodes number
## first identify all the class of the 572 nodes
## 56 class levels for overall
sort(unique(microeco.network$res_node_table$Class))

class.col = c("Acetothermiia" = "#BBCEA4", "Acidobacteriae" = "#EB8CBB", "Alphaproteobacteria" =  "#762AED", "Anaerolineae" =  "#A08666","Aquificae" = "#5D97B6", "Babeliae" = "#D7A9B9", "Bacilli" =  "#83B834", "Bacteroidia" = "#38AEC3", "Blastocatellia" =   "#CD58F3", "Brachyspirae" =  "#C3EF4B", "Brevinematia" = "#DDAB52", "Calditrichia" = "#B672AA", "Chitinivibrionia" =  "#5C97EC", "Chlorobia" =  "#4EE4BD", "Chloroflexia" = "#DF764F", "Chthonomonadetes" = "#AD9CB7", "Clostridia" = "#51669E", "Cyanobacteriia" = "#027742", "Deinococci" = "#CE877A", "Desulfovibrionia" =  "#fcbabd", "Desulfuromonadia" = "#C4EEE9" , "Dissulfuribacteria"= "#9457BE", "Elusimicrobia" = "#A692C7", "Fimbriimonadia" = "#B1A8A3", "Gammaproteobacteria" = "#e879fc", "Gemmatimonadetes" = "#4DECAC", "Ignavibacteria"  = "#53F0EF", "Kapabacteria" =  "#CF6785", "Kiritimatiellae" =  "#8C6072", "Kryptonia" =  "#C6BA97", "Leptospirae" =  "#60B6F0", "Microgenomatia" = "#EB64AB", "Myxococcia" = "#DCE196", "Nanoarchaeia" = "#BAB042", "Negativicutes" = "#58682E", "Nitrososphaeria" =  "#ACBED5", "Parcubacteria"= "#A156E2","Phycisphaerae" =  "#CF8AA4", "Planctomycetes"  ="#BAD0C2", "Polyangia"= "#E697E9", "Spirochaetia" = "#53F0EF", "Sumerlaeia" = "#ECD6A2",   "Synergistia" = "#56A096","Thermodesulfobacteria" =  "#EAC5AE", "Thermodesulfovibrionia" =  "#5EEB3C", "Thermotogae" = "#E32C97", "Unassigned Armatimonadota (Phylum)" =  "#9CB2A9",  "Unassigned Bacteria (Kingdom)" = "#758391", "Unassigned Desulfobacterota (Phylum)" = "#5C97EC", "Unassigned MBNT15 (Phylum)" = "#C8F0CF" , "Unassigned Planctomycetota (Phylum)" = "#913E84", "Unassigned Sva0485 (Phylum)"= "#ABAAED", "Unassigned TA06 (Phylum)" = "#55CDD2", "vadinHA49" = "#8974EE", "Vampirivibrionia" =  "#BAE3F4", "Verrucomicrobiae" =  "#EAE77A")  

?microeco
help(cluster_walktrap)/cluster_fast_greedy
library(randomcoloR)
n <- 100
palette <- distinctColorPalette(100)
palette <- randomColor(75)


microeco.network$plot_taxa_roles(use_type = 2,use_level = "Class",plot_color = "Class" , color_values = class.col, show_number = 1:56, label_text_size = 5) + stat_summary(fun=mean, geom="point",  color = "red", fill = "red", size=1, shape=23, position=position_dodge(0.75)) + stat_summary(aes(label=round(..y..,2)), fun=mean, geom="text", size=2, vjust = -0.5, fontface = "bold")

### NT taxa roles #####

NT.microeco.network.node.table <- NT.microeco.network$res_node_table
sort(unique(NT.microeco.network.node.table$Class)) # 46 levels

png("NT.microeco.network.taxa.roles.all.png", width = 1400, height = 600)
NT.microeco.network$plot_taxa_roles(use_type = 2,use_level = "Class",plot_color = "Class", color_values = class.col, show_number = 1:56, label_text_size = 5) + stat_summary(fun=mean, geom="point",  color = "red", fill = "red", size=1.5, shape=23, position=position_dodge(0.75)) 
dev.off()

png("NT.microeco.network.taxa.roles.10.png", width = 750, height = 600)
NT.microeco.network$plot_taxa_roles(use_type = 2,use_level = "Class",plot_color = "Class", color_values = class.col, show_number = 1:10, label_text_size = 5) + stat_summary(fun=mean, geom="point",  color = "red", fill = "red", size=2.5, shape=23, position=position_dodge(0.75)) 
dev.off()

pdf("NT.microeco.network.taxa.roles.all.pdf", width = 18, height = 8)
NT.microeco.network$plot_taxa_roles(use_type = 2,use_level = "Class",plot_color = "Class", color_values = class.col, show_number = 1:56, label_text_size = 5) + stat_summary(fun=mean, geom="point",  color = "red", fill = "red", size=1.5, shape=23, position=position_dodge(0.75)) 
dev.off()

pdf("NT.microeco.network.taxa.roles.10.pdf", width = 9, height = 7)
NT.microeco.network$plot_taxa_roles(use_type = 2,use_level = "Class",plot_color = "Class", color_values = class.col, show_number = 1:10, label_text_size = 4) + stat_summary(fun=mean, geom="point",  color = "red", fill = "red", size=2.5, shape=23, position=position_dodge(0.75)) 
dev.off()

### NM taxa roles #####

NM.microeco.network.node.table <- NM.microeco.network$res_node_table
sort(unique(NM.microeco.network.node.table$Class)) # 46 levels

png("NM.microeco.network.taxa.roles.all.png", width = 1400, height = 600)
NM.microeco.network$plot_taxa_roles(use_type = 2,use_level = "Class",plot_color = "Class", color_values = class.col, show_number = 1:56, label_text_size = 5) + stat_summary(fun=mean, geom="point",  color = "red", fill = "red", size=1.5, shape=23, position=position_dodge(0.75)) 
dev.off()

png("NM.microeco.network.taxa.roles.10.png", width = 750, height = 600)
NM.microeco.network$plot_taxa_roles(use_type = 2,use_level = "Class",plot_color = "Class", color_values = class.col, show_number = 1:10, label_text_size = 5) + stat_summary(fun=mean, geom="point",  color = "red", fill = "red", size=2.5, shape=23, position=position_dodge(0.75)) 
dev.off()

pdf("NM.microeco.network.taxa.roles.all.pdf", width = 18, height = 8)
NM.microeco.network$plot_taxa_roles(use_type = 2,use_level = "Class",plot_color = "Class", color_values = class.col, show_number = 1:56, label_text_size = 5) + stat_summary(fun=mean, geom="point",  color = "red", fill = "red", size=1.5, shape=23, position=position_dodge(0.75)) 
dev.off()

pdf("NM.microeco.network.taxa.roles.10.pdf", width = 9, height = 7)
NM.microeco.network$plot_taxa_roles(use_type = 2,use_level = "Class",plot_color = "Class", color_values = class.col, show_number = 1:10, label_text_size = 4) + stat_summary(fun=mean, geom="point",  color = "red", fill = "red", size=2.5, shape=23, position=position_dodge(0.75)) 
dev.off()


### ST taxa roles #####

ST.microeco.network.node.table <- ST.microeco.network$res_node_table
sort(unique(ST.microeco.network.node.table$Class)) # 46 levels

png("ST.microeco.network.taxa.roles.all.png", width = 1400, height = 600)
ST.microeco.network$plot_taxa_roles(use_type = 2,use_level = "Class",plot_color = "Class", color_values = class.col, show_number = 1:56, label_text_size = 5) + stat_summary(fun=mean, geom="point",  color = "red", fill = "red", size=1.5, shape=23, position=position_dodge(0.75)) 
dev.off()

png("ST.microeco.network.taxa.roles.10.png", width = 750, height = 600)
ST.microeco.network$plot_taxa_roles(use_type = 2,use_level = "Class",plot_color = "Class", color_values = class.col, show_number = 1:10, label_text_size = 5) + stat_summary(fun=mean, geom="point",  color = "red", fill = "red", size=2.5, shape=23, position=position_dodge(0.75)) 
dev.off()

pdf("ST.microeco.network.taxa.roles.all.pdf", width = 18, height = 8)
ST.microeco.network$plot_taxa_roles(use_type = 2,use_level = "Class",plot_color = "Class", color_values = class.col, show_number = 1:56, label_text_size = 5) + stat_summary(fun=mean, geom="point",  color = "red", fill = "red", size=1.5, shape=23, position=position_dodge(0.75)) 
dev.off()

pdf("ST.microeco.network.taxa.roles.10.pdf", width = 9, height = 7)
ST.microeco.network$plot_taxa_roles(use_type = 2,use_level = "Class",plot_color = "Class", color_values = class.col, show_number = 1:10, label_text_size = 4) + stat_summary(fun=mean, geom="point",  color = "red", fill = "red", size=2.5, shape=23, position=position_dodge(0.75)) 
dev.off()



### SM taxa roles #####

SM.microeco.network.node.table <- SM.microeco.network$res_node_table
sort(unique(SM.microeco.network.node.table$Class)) # 46 levels

png("SM.microeco.network.taxa.roles.all.png", width = 1400, height = 600)
SM.microeco.network$plot_taxa_roles(use_type = 2,use_level = "Class",plot_color = "Class", color_values = class.col, show_number = 1:56, label_text_size = 5) + stat_summary(fun=mean, geom="point",  color = "red", fill = "red", size=1.5, shape=23, position=position_dodge(0.75)) 
dev.off()

png("SM.microeco.network.taxa.roles.10.png", width = 750, height = 600)
SM.microeco.network$plot_taxa_roles(use_type = 2,use_level = "Class",plot_color = "Class", color_values = class.col, show_number = 1:10, label_text_size = 5) + stat_summary(fun=mean, geom="point",  color = "red", fill = "red", size=2.5, shape=23, position=position_dodge(0.75)) 
dev.off()

pdf("SM.microeco.network.taxa.roles.all.pdf", width = 18, height = 8)
SM.microeco.network$plot_taxa_roles(use_type = 2,use_level = "Class",plot_color = "Class", color_values = class.col, show_number = 1:56, label_text_size = 5) + stat_summary(fun=mean, geom="point",  color = "red", fill = "red", size=1.5, shape=23, position=position_dodge(0.75)) 
dev.off()

pdf("SM.microeco.network.taxa.roles.10.pdf", width = 9, height = 7)
SM.microeco.network$plot_taxa_roles(use_type = 2,use_level = "Class",plot_color = "Class", color_values = class.col, show_number = 1:10, label_text_size = 4) + stat_summary(fun=mean, geom="point",  color = "red", fill = "red", size=2.5, shape=23, position=position_dodge(0.75)) 
dev.off()


### CT taxa roles #####

CT.microeco.network.node.table <- CT.microeco.network$res_node_table
sort(unique(CT.microeco.network.node.table$Class)) # 46 levels

png("CT.microeco.network.taxa.roles.all.png", width = 1400, height = 600)
CT.microeco.network$plot_taxa_roles(use_type = 2,use_level = "Class",plot_color = "Class", color_values = class.col, show_number = 1:56, label_text_size = 5) + stat_summary(fun=mean, geom="point",  color = "red", fill = "red", size=1.5, shape=23, position=position_dodge(0.75)) 
dev.off()

png("CT.microeco.network.taxa.roles.10.png", width = 750, height = 600)
CT.microeco.network$plot_taxa_roles(use_type = 2,use_level = "Class",plot_color = "Class", color_values = class.col, show_number = 1:10, label_text_size = 5) + stat_summary(fun=mean, geom="point",  color = "red", fill = "red", size=2.5, shape=23, position=position_dodge(0.75)) 
dev.off()

pdf("CT.microeco.network.taxa.roles.all.pdf", width = 18, height = 8)
CT.microeco.network$plot_taxa_roles(use_type = 2,use_level = "Class",plot_color = "Class", color_values = class.col, show_number = 1:56, label_text_size = 5) + stat_summary(fun=mean, geom="point",  color = "red", fill = "red", size=1.5, shape=23, position=position_dodge(0.75)) 
dev.off()

pdf("CT.microeco.network.taxa.roles.10.pdf", width = 9, height = 7)
CT.microeco.network$plot_taxa_roles(use_type = 2,use_level = "Class",plot_color = "Class", color_values = class.col, show_number = 1:10, label_text_size = 4) + stat_summary(fun=mean, geom="point",  color = "red", fill = "red", size=2.5, shape=23, position=position_dodge(0.75)) 
dev.off()


### SG taxa roles #####

SG.microeco.network.node.table <- SG.microeco.network$res_node_table
sort(unique(SG.microeco.network.node.table$Class)) # 46 levels

png("SG.microeco.network.taxa.roles.all.png", width = 1400, height = 600)
SG.microeco.network$plot_taxa_roles(use_type = 2,use_level = "Class",plot_color = "Class", color_values = class.col, show_number = 1:56, label_text_size = 5) + stat_summary(fun=mean, geom="point",  color = "red", fill = "red", size=1.5, shape=23, position=position_dodge(0.75)) 
dev.off()

png("SG.microeco.network.taxa.roles.10.png", width = 750, height = 600)
SG.microeco.network$plot_taxa_roles(use_type = 2,use_level = "Class",plot_color = "Class", color_values = class.col, show_number = 1:10, label_text_size = 5) + stat_summary(fun=mean, geom="point",  color = "red", fill = "red", size=2.5, shape=23, position=position_dodge(0.75)) 
dev.off()

pdf("SG.microeco.network.taxa.roles.all.pdf", width = 18, height = 8)
SG.microeco.network$plot_taxa_roles(use_type = 2,use_level = "Class",plot_color = "Class", color_values = class.col, show_number = 1:56, label_text_size = 5) + stat_summary(fun=mean, geom="point",  color = "red", fill = "red", size=1.5, shape=23, position=position_dodge(0.75)) 
dev.off()

pdf("SG.microeco.network.taxa.roles.10.pdf", width = 9, height = 7)
SG.microeco.network$plot_taxa_roles(use_type = 2,use_level = "Class",plot_color = "Class", color_values = class.col, show_number = 1:10, label_text_size = 4) + stat_summary(fun=mean, geom="point",  color = "red", fill = "red", size=2.5, shape=23, position=position_dodge(0.75)) 
dev.off()




#Now, we show the eigengene analysis of modules. The eigengene of a module, i.e. the first principal component of PCA, represents the main variance of the abundance in the species of the module.
microeco.network$cal_eigen()
#microeco.network$res_eigen

NT.microeco.network$cal_eigen()
#microeco.network$res_eigen

NM.microeco.network$cal_eigen()
ST.microeco.network$cal_eigen()
SM.microeco.network$cal_eigen()
CT.microeco.network$cal_eigen()
SG.microeco.network$cal_eigen()


#perform correlation heatmap to show the associations between eigengenes and environmental factors.
# create trans_env object
microeco.network_correlation <- trans_env$new(dataset = microeco.dataset, add_data = metadata.micronet[, 23:29])

microeco.network_correlation$cal_cor(add_abund_table = microeco.network$res_eigen)

# plot the correlation heatmap
pdf("microeco.network.spieceasi.correlation.modules.pdf", height = 3)
microeco.network_correlation$plot_cor()
dev.off()

# default parameter represents using igraph plot.igraph function
pdf("microeco.network.pdf")
microeco.network$plot_network(method = "igraph", layout = layout_with_kk,node_color = "module")
dev.off()

## Chord diagrams : cal_sum_links ####
## This function is used to sum the links number from one taxa to another or in the same taxa, for example, at Phylum level. This is very useful to fast see how many nodes are connected between different taxa or within the taxa.

## use this chorddiagram as its better represents the original network when periperal nodes were compared
## use these to adjustnparameters : https://www.rdocumentation.org/packages/chorddiag/versions/0.1.3/topics/chorddiag

#microeco.network$cal_sum_links(taxa_level = "Family")

#microeco.network$plot_sum_links(plot_pos = TRUE, plot_num = 15, color_values = color_values_2, groupnameFontsize = 12, showTicks = FALSE, groupnamePadding = 2, margin = 175)

microeco.network$cal_sum_links(taxa_level = "Genus")
NT.microeco.network$cal_sum_links(taxa_level = "Genus")
NM.microeco.network$cal_sum_links(taxa_level = "Genus")
ST.microeco.network$cal_sum_links(taxa_level = "Genus")
SM.microeco.network$cal_sum_links(taxa_level = "Genus")
CT.microeco.network$cal_sum_links(taxa_level = "Genus")
SG.microeco.network$cal_sum_links(taxa_level = "Genus")


## To save : microeco.network$plot_sum_links(plot_pos = TRUE, plot_num = 15, color_values = color_values_2, groupnameFontsize = 7 or 5.5, showTicks = FALSE, groupnamePadding = 2, margin = 50,showGroupnames = TRUE)

# export as image png(750,900)

##Plot sum links ####
## Plot the summed linkages among taxa.

##########################################################################

microeco.network$plot_sum_links(plot_pos = TRUE, plot_num = 5, color_values = color_values_2, groupnameFontsize = 5.5, showTicks = FALSE, groupnamePadding = 2)


## compare this with the network constructed and analysed with microeco

microeco_transnetwork_spieceasi$plot_sum_links(plot_pos = TRUE ,plot_num = 15, color_values = color_values_2, groupnameFontsize = 5.5, showTicks = FALSE, groupnamePadding = 2)

## both are the same

microeco_transnetwork_spieceasi$plot_sum_links(plot_pos = FALSE ,plot_num = 15, color_values = color_values_2, groupnameFontsize = 5.5, showTicks = FALSE, groupnamePadding = 2)


########################################################

## get list of top 15  genus for all 6 regions

test <- NT.microeco.network$res_sum_links_pos

genus <- sort(unique(c(colnames(NT.microeco.network$res_sum_links_pos)[1:15], colnames(NM.microeco.network$res_sum_links_pos)[1:15], colnames(SM.microeco.network$res_sum_links_pos)[1:15], colnames(ST.microeco.network$res_sum_links_pos)[1:15], colnames(CT.microeco.network$res_sum_links_pos)[1:15], colnames(SG.microeco.network$res_sum_links_pos)[1:15] )))

genus.col.2 <- c("#1B9E77", "#B672AA" ,"#A6CEE3", "#CAB206", "#551A8B" ,"#95CDB3", "#FF7F00", "#1F78B4", "#CE877A", "#9BEAAF", "#7570B3", "#B15928", "#4DECAC", "#9457BE", "#92A362", "#D95F02", "#B2DF8A" , "#990066" ,"#33A02C" ,"#C4EEE9", "#E31A1C", "#B4CDCD", "#9B3741", "#8B6508", "#CDBE99", "#FB9A99", "#FFFF99", "goldenrod")


genus.col.NM = genus.col.2[rownames(NM.microeco.network$res_sum_links_pos)[1:15]] %>% unname
genus.col.NT = genus.col.2[rownames(NT.microeco.network$res_sum_links_pos)[1:15]] %>% unname
genus.col.SM = genus.col.2[rownames(SM.microeco.network$res_sum_links_pos)[1:15]] %>% unname
genus.col.ST = genus.col.2[rownames(ST.microeco.network$res_sum_links_pos)[1:15]] %>% unname
genus.col.CT = genus.col.2[rownames(CT.microeco.network$res_sum_links_pos)[1:15]] %>% unname
genus.col.SG.15 = genus.col.2[rownames(SG.microeco.network$res_sum_links_pos)[1:15]] %>% unname
genus.col.SG.16 = genus.col.2[rownames(SG.microeco.network$res_sum_links_pos)[1:16]] %>% unname


NT.microeco.network$plot_sum_links(plot_pos = TRUE, plot_num = 15, color_values = genus.col.NT , groupnameFontsize = 8, showTicks = FALSE, groupnamePadding = 2, margin = 140,showGroupnames = TRUE)

#NT.microeco.network$plot_sum_links(plot_pos = TRUE, plot_num = 15, color_values = genus.col.2 , groupnameFontsize = 11.7, showTicks = FALSE, groupnamePadding = 2, margin = 50,showGroupnames = TRUE)

NM.microeco.network$plot_sum_links(plot_pos = TRUE, plot_num = 15, color_values = genus.col.NM , groupnameFontsize = 8, showTicks = FALSE, groupnamePadding = 2, margin = 140,showGroupnames = TRUE)

SM.microeco.network$plot_sum_links(plot_pos = TRUE, plot_num = 15, color_values = genus.col.SM , groupnameFontsize = 8, showTicks = FALSE, groupnamePadding = 2, margin = 140,showGroupnames = TRUE)

ST.microeco.network$plot_sum_links(plot_pos = TRUE, plot_num = 15, color_values = genus.col.ST , groupnameFontsize = 8, showTicks = FALSE, groupnamePadding = 2, margin = 140,showGroupnames = TRUE)

CT.microeco.network$plot_sum_links(plot_pos = TRUE, plot_num = 15, color_values = genus.col.CT , groupnameFontsize = 8, showTicks = FALSE, groupnamePadding = 2, margin = 140,showGroupnames = TRUE)

SG.microeco.network$plot_sum_links(plot_pos = TRUE, plot_num = 15, color_values = genus.col.SG.15 , groupnameFontsize = 8, showTicks = FALSE, groupnamePadding = 2, margin = 140,showGroupnames = TRUE)

SG.microeco.network$plot_sum_links(plot_pos = TRUE, plot_num = 16, color_values = genus.col.SG.16 , groupnameFontsize = 8, showTicks = FALSE, groupnamePadding = 2, margin = 140,showGroupnames = TRUE)



#NT.microeco.network$plot_sum_links(TRUE, plot_num = 15, color_values = color_values_2,method = "circlize", transparency = 0.5, annotationTrackHeight = circlize::mm_h(c(5, 5)))

## with no lables #####

microeco.network$plot_sum_links(plot_pos = TRUE, plot_num = 15, color_values = color_values_2, groupnameFontsize = 11.7, showTicks = FALSE, groupnamePadding = 2, margin = 50, showGroupnames = FALSE)


NT.microeco.network$plot_sum_links(plot_pos = TRUE, plot_num = 15, color_values = genus.col.NT , groupnameFontsize = 11.7, showTicks = FALSE, groupnamePadding = 2, margin = 50, showGroupnames = FALSE)

NM.microeco.network$plot_sum_links(plot_pos = TRUE, plot_num = 15, color_values = genus.col.NM , groupnameFontsize = 11.7, showTicks = FALSE, groupnamePadding = 2, margin = 50, showGroupnames = FALSE)

ST.microeco.network$plot_sum_links(plot_pos = TRUE, plot_num = 15, color_values = genus.col.ST , groupnameFontsize = 11.7, showTicks = FALSE, groupnamePadding = 2, margin = 50, showGroupnames = FALSE)

SM.microeco.network$plot_sum_links(plot_pos = TRUE, plot_num = 15, color_values = genus.col.SM , groupnameFontsize = 11.7, showTicks = FALSE, groupnamePadding = 2, margin = 50, showGroupnames = FALSE)

CT.microeco.network$plot_sum_links(plot_pos = TRUE, plot_num = 15, color_values = genus.col.CT , groupnameFontsize = 11.7, showTicks = FALSE, groupnamePadding = 2, margin = 50, showGroupnames = FALSE)

SG.microeco.network$plot_sum_links(plot_pos = TRUE, plot_num = 15, color_values = genus.col.SG.15 , groupnameFontsize = 11.7, showTicks = FALSE, groupnamePadding = 2, margin = 50, showGroupnames = FALSE)

SG.microeco.network$plot_sum_links(plot_pos = TRUE, plot_num = 16, color_values = genus.col.SG.16 , groupnameFontsize = 11.7, showTicks = FALSE, groupnamePadding = 2, margin = 50, showGroupnames = FALSE)


#########################################################################################

### Fig. S9. Variation of abiotic variables between regions. 
### Tests for collinearity were performed for the full dataset using a Pearson’s correlation coefficient 

### CHECK #####

######################################################################################







#########################################################################################

### Fig. 5A)i)   Mantel's test on the influence of abiotic variables on observed taxonomic composition for all taxa, photosynthetic, chemoautotrophic, and heterotrophic fractions of the community.####
### Fig. 5A)ii)  Variance Partitioning for overall dataset (using transformed Hellinger distances versus abiotic variables indicating the relative contribution from the four most influential variables) (excluding collinear and latitude) ####
### Fig. 5A)iii) relative contribution of the four topmost variables to observed variation in community distribution at different spatial scales ####

### Fig. 5B. Random Forest ensemble modelling of the most influential variables on observed community distribution for the overall meta-community, each region, and ecological guilds ####
### Fig. S11B.Random Forest modeling of major abiotic influences on microbial communities with distance (latitude) included as a variable. ######

### CHECK #####

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



#### Random Forest ####

### Datasets for random forests : phyloseqs by regions and ecological guilds ####

### Get phyloseqs for each of the 6 regions with all ASVs ####

### Get phyloseqs with the 3 ecological guilds for each of the regions #####

## all photosynthetic

photosynthetic.class <- c("Cyanobacteriia","Chloroflexia")

proteobacteria.genus <- c("DSSF69","Elioraea","Methylobacterium-Methylorubrum","Rhodomicrobium","Roseomonas","Sandaracinobacter","Tabrizicola","Unassigned Rhodobacteraceae (Family)","Unassigned Sphingomonadaceae (Family)","AAP99","Allochromatium","Caldimonas","Curvibacter","DSSD61","Thiolamprovum","Unassigned B1-7BS (Family)","Unassigned Burkholderiales (Order)","Unassigned Comamonadaceae (Family)","Unassigned Gammaproteobacteria (Class)","Unassigned Rhodocyclaceae (Family)","Unassigned Sutterellaceae (Family)","Z-35")

other.photo <- c("Chlorobiales","Chloracidobacteriales")

## all photosynthetic phyloseq
photosynthetic.class <- c("Cyanobacteriia","Chloroflexia")

#rarefied.photosynthetic.prop.exclude.1.3 <- com.rarefied.min.prop.exclude.1.3 %>% subset_taxa( Class %in% photosynthetic.class | Order %in% other.photo | Genus %in% proteobacteria.genus )


## Chemolitotrophs

chemolithotroph.phylum <- c("Aquificota","Calditrichota","Desulfobacterota","Elusimicrobiota","Nitrospirota","Patescibacteria","Sva0485")

chemolithotroph.class <- c("Aquificae","Calditrichia","Leptospirae")

##rarefied.chemolithotrophs.prop.exclude.1.3 <- com.rarefied.min.prop.exclude.1.3 %>% subset_taxa(Phylum %in% chemolithotroph.phylum  |  Class %in% chemolithotroph.class)


#### all  heterotrophs only

#rarefied.only.heterotrophs.prop.exclude.1.3 <- com.rarefied.min.prop.exclude.1.3 %>% subset_taxa(!Phylum %in% chemolithotroph.phylum & !Class %in% chemolithotroph.class & !Class %in% photosynthetic.class & !Order %in% other.photo & !Genus %in% proteobacteria.genus )

#rarefied.only.heterotrophs.prop.wide.exclude.1.3 <- as.matrix(as.data.frame(rarefied.only.heterotrophs.prop.exclude.1.3@otu_table))

#test <- as.matrix(as.data.frame(rarefied.only.heterotrophs.prop.exclude.1.3@tax_table))


## cluster 1 :North.Thailand <- c('PT','PP','LN','HS','MK') #####

NT.loc <- c("PT","PP","LN","HS","MK")

test <- as.data.frame(NT.rarefied.min.int.exclude.1.3.rooted.spnet@tax_table)

## NT by ecological guilds

## photosynthetic
NT.rarefied.photosynthetic.min.int.exclude.1.3.rooted.spnet <- NT.rarefied.min.int.exclude.1.3.rooted.spnet %>% subset_taxa( Class %in% photosynthetic.class | Order %in% other.photo | Genus %in% proteobacteria.genus )

NT.rarefied.photosynthetic.min.int.exclude.1.3.rooted.spnet <- 
  prune_taxa(taxa_sums(NT.rarefied.photosynthetic.min.int.exclude.1.3.rooted.spnet) > 0, NT.rarefied.photosynthetic.min.int.exclude.1.3.rooted.spnet)

NT.test <- (as.data.frame(NT.rarefied.photosynthetic.min.int.exclude.1.3.rooted.spnet@otu_table))

NT.test.1 <- (as.data.frame(NT.rarefied.photosynthetic.min.int.exclude.1.3.rooted.spnet@tax_table))

## chemolithotrophic
NT.rarefied.chemolithotrophic.min.int.exclude.1.3.rooted.spnet <- NT.rarefied.min.int.exclude.1.3.rooted.spnet %>% subset_taxa(Phylum %in% chemolithotroph.phylum  |  Class %in% chemolithotroph.class)

NT.rarefied.chemolithotrophic.min.int.exclude.1.3.rooted.spnet <- 
  prune_taxa(taxa_sums(NT.rarefied.chemolithotrophic.min.int.exclude.1.3.rooted.spnet) > 0, NT.rarefied.chemolithotrophic.min.int.exclude.1.3.rooted.spnet)

NT.test <- (as.data.frame(NT.rarefied.chemolithotrophic.min.int.exclude.1.3.rooted.spnet@otu_table))

NT.test.1 <- (as.data.frame(NT.rarefied.chemolithotrophic.min.int.exclude.1.3.rooted.spnet@tax_table))


##heterotrophic

NT.rarefied.heterotrophic.min.int.exclude.1.3.rooted.spnet <- NT.rarefied.min.int.exclude.1.3.rooted.spnet %>% subset_taxa(!Phylum %in% chemolithotroph.phylum & !Class %in% chemolithotroph.class & !Class %in% photosynthetic.class & !Order %in% other.photo & !Genus %in% proteobacteria.genus )

NT.rarefied.heterotrophic.min.int.exclude.1.3.rooted.spnet <- 
  prune_taxa(taxa_sums(NT.rarefied.heterotrophic.min.int.exclude.1.3.rooted.spnet) > 0, NT.rarefied.heterotrophic.min.int.exclude.1.3.rooted.spnet)

test <- (as.data.frame(NT.rarefied.heterotrophic.min.int.exclude.1.3.rooted.spnet@otu_table))

test.1 <- (as.data.frame(NT.rarefied.heterotrophic.min.int.exclude.1.3.rooted.spnet@tax_table))


## Cluster 2 : North.Malaysia <- c('AP','BA','US') ######

NM.loc <- c("AP","BA","US")

test <- NM.rarefied.min.int.exclude.1.3.rooted.spnet@sam_data


## NM by ecological guilds

## photosynthetic
NM.rarefied.photosynthetic.min.int.exclude.1.3.rooted.spnet <- NM.rarefied.min.int.exclude.1.3.rooted.spnet %>% subset_taxa( Class %in% photosynthetic.class | Order %in% other.photo | Genus %in% proteobacteria.genus)

NM.rarefied.photosynthetic.min.int.exclude.1.3.rooted.spnet  <- 
  prune_taxa(taxa_sums(NM.rarefied.photosynthetic.min.int.exclude.1.3.rooted.spnet ) > 0, NM.rarefied.photosynthetic.min.int.exclude.1.3.rooted.spnet)


NM.test <- (as.data.frame(NM.rarefied.photosynthetic.min.int.exclude.1.3.rooted.spnet@otu_table))

NM.test.1 <- (as.data.frame(NM.rarefied.photosynthetic.min.int.exclude.1.3.rooted.spnet@tax_table))

## chemolithotrophic
NM.rarefied.chemolithotrophic.min.int.exclude.1.3.rooted.spnet <- NM.rarefied.min.int.exclude.1.3.rooted.spnet %>% subset_taxa(Phylum %in% chemolithotroph.phylum  |  Class %in% chemolithotroph.class)

NM.rarefied.chemolithotrophic.min.int.exclude.1.3.rooted.spnet  <- 
  prune_taxa(taxa_sums(NM.rarefied.chemolithotrophic.min.int.exclude.1.3.rooted.spnet) > 0, NM.rarefied.chemolithotrophic.min.int.exclude.1.3.rooted.spnet)

NM.test <- (as.data.frame(NM.rarefied.chemolithotrophic.min.int.exclude.1.3.rooted.spnet@otu_table))

NM.test.1 <- (as.data.frame(NM.rarefied.chemolithotrophic.min.int.exclude.1.3.rooted.spnet@tax_table))


##heterotrophic

NM.rarefied.heterotrophic.min.int.exclude.1.3.rooted.spnet <- NM.rarefied.min.int.exclude.1.3.rooted.spnet %>% subset_taxa(!Phylum %in% chemolithotroph.phylum & !Class %in% chemolithotroph.class & !Class %in% photosynthetic.class & !Order %in% other.photo & !Genus %in% proteobacteria.genus )

NM.rarefied.heterotrophic.min.int.exclude.1.3.rooted.spnet  <- 
  prune_taxa(taxa_sums(NM.rarefied.heterotrophic.min.int.exclude.1.3.rooted.spnet) > 0, NM.rarefied.heterotrophic.min.int.exclude.1.3.rooted.spnet)

test <- (as.data.frame(NM.rarefied.heterotrophic.min.int.exclude.1.3.rooted.spnet@otu_table))

test.1 <- (as.data.frame(NM.rarefied.heterotrophic.min.int.exclude.1.3.rooted.spnet@tax_table))

## Cluster 3 : South.Malaysia <- c('KJ','SE','LA') ######

SM.loc <- c("SE","LA","KJ") 

test <- SM.rarefied.min.int.exclude.1.3.rooted.spnet@sam_data

## SM by ecological guilds

## photosynthetic
SM.rarefied.photosynthetic.min.int.exclude.1.3.rooted.spnet <- SM.rarefied.min.int.exclude.1.3.rooted.spnet %>% subset_taxa( Class %in% photosynthetic.class | Order %in% other.photo | Genus %in% proteobacteria.genus )

SM.rarefied.photosynthetic.min.int.exclude.1.3.rooted.spnet  <- 
  prune_taxa(taxa_sums(SM.rarefied.photosynthetic.min.int.exclude.1.3.rooted.spnet) > 0, SM.rarefied.photosynthetic.min.int.exclude.1.3.rooted.spnet)

SM.test <- (as.data.frame(SM.rarefied.photosynthetic.min.int.exclude.1.3.rooted.spnet@otu_table))

SM.test.1 <- (as.data.frame(SM.rarefied.photosynthetic.min.int.exclude.1.3.rooted.spnet@tax_table))

## chemolithotrophic
SM.rarefied.chemolithotrophic.min.int.exclude.1.3.rooted.spnet <- SM.rarefied.min.int.exclude.1.3.rooted.spnet %>% subset_taxa(Phylum %in% chemolithotroph.phylum  |  Class %in% chemolithotroph.class)

SM.rarefied.chemolithotrophic.min.int.exclude.1.3.rooted.spnet <- 
  prune_taxa(taxa_sums(SM.rarefied.chemolithotrophic.min.int.exclude.1.3.rooted.spnet) > 0, SM.rarefied.chemolithotrophic.min.int.exclude.1.3.rooted.spnet)

test <- (as.data.frame(SM.rarefied.chemolithotrophic.min.int.exclude.1.3.rooted.spnet@otu_table))

test.1 <- (as.data.frame(SM.rarefied.chemolithotrophic.min.int.exclude.1.3.rooted.spnet@tax_table))


##heterotrophic

SM.rarefied.heterotrophic.min.int.exclude.1.3.rooted.spnet <- SM.rarefied.min.int.exclude.1.3.rooted.spnet %>% subset_taxa(!Phylum %in% chemolithotroph.phylum & !Class %in% chemolithotroph.class & !Class %in% photosynthetic.class & !Order %in% other.photo & !Genus %in% proteobacteria.genus )

SM.rarefied.heterotrophic.min.int.exclude.1.3.rooted.spnet <- 
  prune_taxa(taxa_sums(SM.rarefied.heterotrophic.min.int.exclude.1.3.rooted.spnet) > 0, SM.rarefied.heterotrophic.min.int.exclude.1.3.rooted.spnet)

test <- (as.data.frame(SM.rarefied.heterotrophic.min.int.exclude.1.3.rooted.spnet@otu_table))

test.1 <- (as.data.frame(SM.rarefied.heterotrophic.min.int.exclude.1.3.rooted.spnet@tax_table))

## Cluster 4 : All Sembawang : Singapore <- 'SW' #####

test <- (as.data.frame(SG.rarefied.min.int.exclude.1.3.rooted.spnet@tax_table))

## SG by ecological guilds

## photosynthetic
SG.rarefied.photosynthetic.min.int.exclude.1.3.rooted.spnet <- SG.rarefied.min.int.exclude.1.3.rooted.spnet %>% subset_taxa( Class %in% photosynthetic.class | Order %in% other.photo | Genus %in% proteobacteria.genus )

SG.rarefied.photosynthetic.min.int.exclude.1.3.rooted.spnet <- 
  prune_taxa(taxa_sums(SG.rarefied.photosynthetic.min.int.exclude.1.3.rooted.spnet) > 0, SG.rarefied.photosynthetic.min.int.exclude.1.3.rooted.spnet)

SG.test <- (as.data.frame(SG.rarefied.photosynthetic.min.int.exclude.1.3.rooted.spnet@otu_table))

SG.test.1 <- (as.data.frame(SG.rarefied.photosynthetic.min.int.exclude.1.3.rooted.spnet@tax_table))

## chemolithotrophic
SG.rarefied.chemolithotrophic.min.int.exclude.1.3.rooted.spnet <- SG.rarefied.min.int.exclude.1.3.rooted.spnet %>% subset_taxa(Phylum %in% chemolithotroph.phylum  |  Class %in% chemolithotroph.class)

SG.rarefied.chemolithotrophic.min.int.exclude.1.3.rooted.spnet <- 
  prune_taxa(taxa_sums(SG.rarefied.chemolithotrophic.min.int.exclude.1.3.rooted.spnet) > 0, SG.rarefied.chemolithotrophic.min.int.exclude.1.3.rooted.spnet)

test <- (as.matrix(SG.rarefied.chemolithotrophic.min.int.exclude.1.3.rooted.spnet@otu_table))

test.1 <- (as.matrix(SG.rarefied.chemolithotrophic.min.int.exclude.1.3.rooted.spnet@tax_table))


##heterotrophic

SG.rarefied.heterotrophic.min.int.exclude.1.3.rooted.spnet <- SG.rarefied.min.int.exclude.1.3.rooted.spnet %>% subset_taxa(!Phylum %in% chemolithotroph.phylum & !Class %in% chemolithotroph.class & !Class %in% photosynthetic.class & !Order %in% other.photo & !Genus %in% proteobacteria.genus )

SG.rarefied.heterotrophic.min.int.exclude.1.3.rooted.spnet <- 
  prune_taxa(taxa_sums(SG.rarefied.heterotrophic.min.int.exclude.1.3.rooted.spnet) > 0, SG.rarefied.heterotrophic.min.int.exclude.1.3.rooted.spnet)

test <- (as.data.frame(SG.rarefied.heterotrophic.min.int.exclude.1.3.rooted.spnet@otu_table))

test.1 <- (as.data.frame(SG.rarefied.heterotrophic.min.int.exclude.1.3.rooted.spnet@tax_table))

## Cluster 5 : South.Thailand <- c('RB','RN') #####

ST <- c("RN","RB")

## ST by ecological guilds

## photosynthetic
ST.rarefied.photosynthetic.min.int.exclude.1.3.rooted.spnet <- ST.rarefied.min.int.exclude.1.3.rooted.spnet %>% subset_taxa( Class %in% photosynthetic.class | Order %in% other.photo | Genus %in% proteobacteria.genus)

ST.rarefied.photosynthetic.min.int.exclude.1.3.rooted.spnet <- 
  prune_taxa(taxa_sums(ST.rarefied.photosynthetic.min.int.exclude.1.3.rooted.spnet) > 0, ST.rarefied.photosynthetic.min.int.exclude.1.3.rooted.spnet)

test <- (as.data.frame(ST.rarefied.photosynthetic.min.int.exclude.1.3.rooted.spnet@otu_table))

test.1 <- (as.data.frame(ST.rarefied.photosynthetic.min.int.exclude.1.3.rooted.spnet@tax_table))

## chemolithotrophic
ST.rarefied.chemolithotrophic.min.int.exclude.1.3.rooted.spnet <- ST.rarefied.min.int.exclude.1.3.rooted.spnet %>% subset_taxa(Phylum %in% chemolithotroph.phylum  |  Class %in% chemolithotroph.class)

ST.rarefied.chemolithotrophic.min.int.exclude.1.3.rooted.spnet <- 
  prune_taxa(taxa_sums(ST.rarefied.chemolithotrophic.min.int.exclude.1.3.rooted.spnet) > 0, ST.rarefied.chemolithotrophic.min.int.exclude.1.3.rooted.spnet)

test <- (as.data.frame(ST.rarefied.chemolithotrophic.min.int.exclude.1.3.rooted.spnet@otu_table))

test.1 <- (as.data.frame(ST.rarefied.chemolithotrophic.min.int.exclude.1.3.rooted.spnet@tax_table))


##heterotrophic

ST.rarefied.heterotrophic.min.int.exclude.1.3.rooted.spnet <- ST.rarefied.min.int.exclude.1.3.rooted.spnet %>% subset_taxa(!Phylum %in% chemolithotroph.phylum & !Class %in% chemolithotroph.class & !Class %in% photosynthetic.class & !Order %in% other.photo & !Genus %in% proteobacteria.genus )

ST.rarefied.heterotrophic.min.int.exclude.1.3.rooted.spnet <- 
  prune_taxa(taxa_sums(ST.rarefied.heterotrophic.min.int.exclude.1.3.rooted.spnet) > 0, ST.rarefied.heterotrophic.min.int.exclude.1.3.rooted.spnet)

test <- (as.data.frame(ST.rarefied.heterotrophic.min.int.exclude.1.3.rooted.spnet@otu_table))

test.1 <- (as.data.frame(ST.rarefied.heterotrophic.min.int.exclude.1.3.rooted.spnet@tax_table))

## Cluster 6 : Central.Thailand <- 'PB' ######

test <- CT.rarefied.min.int.exclude.1.3.rooted.spnet@tax_table

## CT by ecological guilds

## photosynthetic
CT.rarefied.photosynthetic.min.int.exclude.1.3.rooted.spnet <- CT.rarefied.min.int.exclude.1.3.rooted.spnet %>% subset_taxa( Class %in% photosynthetic.class | Order %in% other.photo | Genus %in% proteobacteria.genus )

CT.rarefied.photosynthetic.min.int.exclude.1.3.rooted.spnet <- 
  prune_taxa(taxa_sums(CT.rarefied.photosynthetic.min.int.exclude.1.3.rooted.spnet) > 0, CT.rarefied.photosynthetic.min.int.exclude.1.3.rooted.spnet)

test <- (as.data.frame(CT.rarefied.photosynthetic.min.int.exclude.1.3.rooted.spnet@otu_table))

test.1 <- (as.data.frame(CT.rarefied.photosynthetic.min.int.exclude.1.3.rooted.spnet@tax_table))

## chemolithotrophic
CT.rarefied.chemolithotrophic.min.int.exclude.1.3.rooted.spnet <- CT.rarefied.min.int.exclude.1.3.rooted.spnet %>% subset_taxa(Phylum %in% chemolithotroph.phylum  |  Class %in% chemolithotroph.class)

CT.rarefied.chemolithotrophic.min.int.exclude.1.3.rooted.spnet <- 
  prune_taxa(taxa_sums(CT.rarefied.chemolithotrophic.min.int.exclude.1.3.rooted.spnet) > 0, CT.rarefied.chemolithotrophic.min.int.exclude.1.3.rooted.spnet)

test <- (as.data.frame(CT.rarefied.chemolithotrophic.min.int.exclude.1.3.rooted.spnet@otu_table))

test.1 <- (as.data.frame(CT.rarefied.chemolithotrophic.min.int.exclude.1.3.rooted.spnet@tax_table))


##heterotrophic

CT.rarefied.heterotrophic.min.int.exclude.1.3.rooted.spnet <- CT.rarefied.min.int.exclude.1.3.rooted.spnet %>% subset_taxa(!Phylum %in% chemolithotroph.phylum & !Class %in% chemolithotroph.class & !Class %in% photosynthetic.class & !Order %in% other.photo & !Genus %in% proteobacteria.genus )

CT.rarefied.heterotrophic.min.int.exclude.1.3.rooted.spnet <- 
  prune_taxa(taxa_sums(CT.rarefied.heterotrophic.min.int.exclude.1.3.rooted.spnet) > 0, CT.rarefied.heterotrophic.min.int.exclude.1.3.rooted.spnet)

test <- (as.data.frame(CT.rarefied.heterotrophic.min.int.exclude.1.3.rooted.spnet@otu_table))

test.1 <- (as.data.frame(CT.rarefied.heterotrophic.min.int.exclude.1.3.rooted.spnet@tax_table))


### save cluster phyloseqs

saveRDS(NT.rarefied.photosynthetic.min.int.exclude.1.3.rooted.spnet, "NT.rarefied.photosynthetic.min.int.exclude.1.3.rooted.spnet.rds")
saveRDS(NT.rarefied.chemolithotrophic.min.int.exclude.1.3.rooted.spnet, "NT.rarefied.chemolithotrophic.min.int.exclude.1.3.rooted.spnet.rds")
saveRDS(NT.rarefied.heterotrophic.min.int.exclude.1.3.rooted.spnet, "NT.rarefied.heterotrophic.min.int.exclude.1.3.rooted.spnet.rds")


saveRDS(NM.rarefied.photosynthetic.min.int.exclude.1.3.rooted.spnet, "NM.rarefied.photosynthetic.min.int.exclude.1.3.rooted.spnet.rds")
saveRDS(NM.rarefied.chemolithotrophic.min.int.exclude.1.3.rooted.spnet, "NM.rarefied.chemolithotrophic.min.int.exclude.1.3.rooted.spnet.rds")
saveRDS(NM.rarefied.heterotrophic.min.int.exclude.1.3.rooted.spnet, "NM.rarefied.heterotrophic.min.int.exclude.1.3.rooted.spnet.rds")


saveRDS(SM.rarefied.photosynthetic.min.int.exclude.1.3.rooted.spnet, "SM.rarefied.photosynthetic.min.int.exclude.1.3.rooted.spnet.rds")
saveRDS(SM.rarefied.chemolithotrophic.min.int.exclude.1.3.rooted.spnet, "SM.rarefied.chemolithotrophic.min.int.exclude.1.3.rooted.spnet.rds")
saveRDS(SM.rarefied.heterotrophic.min.int.exclude.1.3.rooted.spnet, "SM.rarefied.heterotrophic.min.int.exclude.1.3.rooted.spnet.rds")


saveRDS(SG.rarefied.photosynthetic.min.int.exclude.1.3.rooted.spnet, "SG.rarefied.photosynthetic.min.int.exclude.1.3.rooted.spnet.rds")
saveRDS(SG.rarefied.chemolithotrophic.min.int.exclude.1.3.rooted.spnet, "SG.rarefied.chemolithotrophic.min.int.exclude.1.3.rooted.spnet.rds")
saveRDS(SG.rarefied.heterotrophic.min.int.exclude.1.3.rooted.spnet, "SG.rarefied.heterotrophic.min.int.exclude.1.3.rooted.spnet.rds")


saveRDS(ST.rarefied.photosynthetic.min.int.exclude.1.3.rooted.spnet, "ST.rarefied.photosynthetic.min.int.exclude.1.3.rooted.spnet.rds")
saveRDS(ST.rarefied.chemolithotrophic.min.int.exclude.1.3.rooted.spnet, "ST.rarefied.chemolithotrophic.min.int.exclude.1.3.rooted.spnet.rds")
saveRDS(ST.rarefied.heterotrophic.min.int.exclude.1.3.rooted.spnet, "ST.rarefied.heterotrophic.min.int.exclude.1.3.rooted.spnet.rds")


saveRDS(CT.rarefied.photosynthetic.min.int.exclude.1.3.rooted.spnet, "CT.rarefied.photosynthetic.min.int.exclude.1.3.rooted.spnet.rds")
saveRDS(CT.rarefied.chemolithotrophic.min.int.exclude.1.3.rooted.spnet, "CT.rarefied.chemolithotrophic.min.int.exclude.1.3.rooted.spnet.rds")
saveRDS(CT.rarefied.heterotrophic.min.int.exclude.1.3.rooted.spnet, "CT.rarefied.heterotrophic.min.int.exclude.1.3.rooted.spnet.rds")

### run unweighted PCOA ordinations for 6 regions as per Rscript 3 for random forest ######

## first Ordinate ####

## cluster 1 :North.Thailand <- c('PT','PP','LN','HS','MK') #####

NT.loc <- c("PT","PP","LN","HS","MK")

NT.unweighted.int.unifrac.PCoA = ordinate(NT.rarefied.min.int.exclude.1.3.rooted.spnet, method="PCoA", distance="unifrac", weighted=FALSE)

test <- as.data.frame(NT.rarefied.min.int.exclude.1.3.rooted.spnet@tax_table)
test <- as.data.frame(NT.unweighted.int.unifrac.PCoA$vectors)

NT.photosynthetic.unweighted.int.unifrac.PCoA = ordinate(NT.rarefied.photosynthetic.min.int.exclude.1.3.rooted.spnet, method="PCoA", distance="unifrac", weighted=FALSE)

NT.chemolithotrophic.unweighted.int.unifrac.PCoA = ordinate(NT.rarefied.chemolithotrophic.min.int.exclude.1.3.rooted.spnet, method="PCoA", distance="unifrac", weighted=FALSE)

NT.heterotrophic.unweighted.int.unifrac.PCoA = ordinate(NT.rarefied.heterotrophic.min.int.exclude.1.3.rooted.spnet , method="PCoA", distance="unifrac", weighted=FALSE)

## Cluster 2 : North.Malaysia <- c('AP','BA','US') ######

NM.loc <- c("AP","BA","US")

NM.unweighted.int.unifrac.PCoA = ordinate(NM.rarefied.min.int.exclude.1.3.rooted.spnet, method="PCoA", distance="unifrac", weighted=FALSE)

NM.photosynthetic.unweighted.int.unifrac.PCoA = ordinate(NM.rarefied.photosynthetic.min.int.exclude.1.3.rooted.spnet, method="PCoA", distance="unifrac", weighted=FALSE)

NM.chemolithotrophic.unweighted.int.unifrac.PCoA = ordinate(NM.rarefied.chemolithotrophic.min.int.exclude.1.3.rooted.spnet, method="PCoA", distance="unifrac", weighted=FALSE)

NM.heterotrophic.unweighted.int.unifrac.PCoA = ordinate(NM.rarefied.heterotrophic.min.int.exclude.1.3.rooted.spnet , method="PCoA", distance="unifrac", weighted=FALSE)




## Cluster 3 : South.Malaysia <- c('KJ','SE','LA') ######

SM.loc <- c("SE","LA","KJ") 

SM.unweighted.int.unifrac.PCoA = ordinate(SM.rarefied.min.int.exclude.1.3.rooted.spnet, method="PCoA", distance="unifrac", weighted=FALSE)

SM.photosynthetic.unweighted.int.unifrac.PCoA = ordinate(SM.rarefied.photosynthetic.min.int.exclude.1.3.rooted.spnet, method="PCoA", distance="unifrac", weighted=FALSE)

SM.chemolithotrophic.unweighted.int.unifrac.PCoA = ordinate(SM.rarefied.chemolithotrophic.min.int.exclude.1.3.rooted.spnet, method="PCoA", distance="unifrac", weighted=FALSE)

SM.heterotrophic.unweighted.int.unifrac.PCoA = ordinate(SM.rarefied.heterotrophic.min.int.exclude.1.3.rooted.spnet , method="PCoA", distance="unifrac", weighted=FALSE)

## Cluster 4 : All Sembawang : Singapore <- 'SG' #####

SG.unweighted.int.unifrac.PCoA = ordinate(SG.rarefied.min.int.exclude.1.3.rooted.spnet, method="PCoA", distance="unifrac", weighted=FALSE)

SG.photosynthetic.unweighted.int.unifrac.PCoA = ordinate(SG.rarefied.photosynthetic.min.int.exclude.1.3.rooted.spnet, method="PCoA", distance="unifrac", weighted=FALSE)

SG.chemolithotrophic.unweighted.int.unifrac.PCoA = ordinate(SG.rarefied.chemolithotrophic.min.int.exclude.1.3.rooted.spnet, method="PCoA", distance="unifrac", weighted=FALSE)

SG.heterotrophic.unweighted.int.unifrac.PCoA = ordinate(SG.rarefied.heterotrophic.min.int.exclude.1.3.rooted.spnet , method="PCoA", distance="unifrac", weighted=FALSE)

## Cluster 5 : South.Thailand <- c('RB','RN') #####

ST.unweighted.int.unifrac.PCoA = ordinate(ST.rarefied.min.int.exclude.1.3.rooted.spnet, method="PCoA", distance="unifrac", weighted=FALSE)

ST.photosynthetic.unweighted.int.unifrac.PCoA = ordinate(ST.rarefied.photosynthetic.min.int.exclude.1.3.rooted.spnet, method="PCoA", distance="unifrac", weighted=FALSE)

ST.chemolithotrophic.unweighted.int.unifrac.PCoA = ordinate(ST.rarefied.chemolithotrophic.min.int.exclude.1.3.rooted.spnet, method="PCoA", distance="unifrac", weighted=FALSE)

ST.heterotrophic.unweighted.int.unifrac.PCoA = ordinate(ST.rarefied.heterotrophic.min.int.exclude.1.3.rooted.spnet , method="PCoA", distance="unifrac", weighted=FALSE)

## Cluster 6 : Central.Thailand <- 'PB' ######

CT.unweighted.int.unifrac.PCoA = ordinate(CT.rarefied.min.int.exclude.1.3.rooted.spnet, method="PCoA", distance="unifrac", weighted=FALSE)

test <- CT.rarefied.min.int.exclude.1.3.rooted.spnet@tax_table

CT.photosynthetic.unweighted.int.unifrac.PCoA = ordinate(CT.rarefied.photosynthetic.min.int.exclude.1.3.rooted.spnet, method="PCoA", distance="unifrac", weighted=FALSE)

CT.chemolithotrophic.unweighted.int.unifrac.PCoA = ordinate(CT.rarefied.chemolithotrophic.min.int.exclude.1.3.rooted.spnet, method="PCoA", distance="unifrac", weighted=FALSE)

CT.heterotrophic.unweighted.int.unifrac.PCoA = ordinate(CT.rarefied.heterotrophic.min.int.exclude.1.3.rooted.spnet , method="PCoA", distance="unifrac", weighted=FALSE)

### Plots ####

#+ stat_ellipse(aes(fill = Location.Code , group = Location.Code), linetype = 0 ,type = "t",  level = 0.99, alpha = .2, geom = "polygon")  + scale_fill_manual(values = region.col,breaks = c('North.Thailand','Central.Thailand','South.Thailand','North.Malaysia','South.Malaysia','Singapore')) 

## +   geom_text_repel(aes(label = ASVno), size = 6, vjust = 1.5, max.overlaps = 200) + geom_point(size=3, position="jitter")

# , fill = guide_legend(nrow=1,byrow=TRUE,override.aes = list(size = 7))

loc.shape.names = c('AP', 'BA', 'HS', 'KJ', 'LA', 'LN', 'MK', 'PB', 'PP', 'PT', 'RB', 'RN', 'SE', 'SW', 'US')
loc.shape <- 0:(length(loc.shape.names)-1)
names(loc.shape) <- loc.shape.names
loc.shape["Taxa"] <- 16
loc.shape

temp.col = c('38' = "mediumorchid1", '40' = "mediumpurple1", '42' = "blueviolet",'44' = "mediumpurple4", '45' = "lightskyblue2",  '46' = "deepskyblue" , '47'= "royalblue1",  '48' = "mediumblue", '49' = "cyan", '50' = "aquamarine2" , '51' = "springgreen2",  '53' = "yellowgreen", '55' = "forestgreen", '57' = "khaki3", '58' = "goldenrod1", '61' = "sienna2", '63' = "firebrick3", '66' = "darkred", 'Taxa' = "black")


clean_background <- theme(plot.background = element_rect("white"),
                          panel.background = element_rect(color = "black", fill = NA, linewidth = 2),
                          panel.grid = element_line("white"),
                          axis.line = element_line("gray25"),
                          axis.text = element_text(size = 12, color = "gray25"),
                          axis.title = element_text(color = "gray25"),
                          legend.text = element_text(size = 12),
                          legend.key = element_rect("white"))


NT.PCoA.UnweightedUniFrac.phyloseq.region.samples <- plot_ordination(NT.rarefied.min.int.exclude.1.3.rooted.spnet, NT.unweighted.int.unifrac.PCoA, type="split", color = "Temp.adj...C.", shape = "Location.Code", title="PCoA on unweighted UniFrac distance for North Thailand") +   geom_point(size=5, position="jitter") +  geom_text(aes(label = ASVno), size = 6, vjust = 1.5, max.overlaps = 200) + geom_point(size=3, position="jitter")  + clean_background  + scale_shape_manual(values= loc.shape)  + scale_color_manual(values = temp.col) + theme(legend.title = element_text(size = 15.5, face = "bold"), legend.text = element_text(size = 14),legend.position = "bottom", legend.direction = "horizontal",legend.box="vertical", legend.margin=margin()) + guides(colour = guide_legend(nrow=1,byrow=TRUE, override.aes = list(size = 7)), shape = guide_legend(nrow=1,byrow=TRUE,override.aes = list(size = 7)) )


NT.PCoA.UnweightedUniFrac.phyloseq.region.samples$layers <- NT.PCoA.UnweightedUniFrac.phyloseq.region.samples$layers[-1]

pdf("NT.PCoA.UnweightedUniFrac.phyloseq.region.samples.pdf", height = 12, width = 17)
NT.PCoA.UnweightedUniFrac.phyloseq.region.samples
dev.off()

saveRDS(NT.PCoA.UnweightedUniFrac.phyloseq.region.samples, "NT.PCoA.UnweightedUniFrac.phyloseq.region.samples.rds")


NT.PCoA.photosynthetic.UnweightedUniFrac.phyloseq.region.samples <- plot_ordination(NT.rarefied.photosynthetic.min.int.exclude.1.3.rooted.spnet, NT.photosynthetic.unweighted.int.unifrac.PCoA, type="split", color = "Temp.adj...C.", shape = "Location.Code", title="PCoA on unweighted UniFrac distance for photosynthetic taxa in North Thailand") +   geom_point(size=5, position="jitter") +  geom_text(aes(label = ASVno), size = 6, vjust = 1.5, max.overlaps = 200) + geom_point(size=3, position="jitter")  + clean_background  + scale_shape_manual(values= loc.shape)  + scale_color_manual(values = temp.col) + theme(legend.title = element_text(size = 15.5, face = "bold"), legend.text = element_text(size = 14),legend.position = "bottom", legend.direction = "horizontal",legend.box="vertical", legend.margin=margin()) + guides(colour = guide_legend(nrow=1,byrow=TRUE, override.aes = list(size = 7)), shape = guide_legend(nrow=1,byrow=TRUE,override.aes = list(size = 7)) )


NT.PCoA.photosynthetic.UnweightedUniFrac.phyloseq.region.samples$layers <- NT.PCoA.photosynthetic.UnweightedUniFrac.phyloseq.region.samples$layers[-1]

pdf("NT.PCoA.photosynthetic.UnweightedUniFrac.phyloseq.region.samples.pdf", height = 12, width = 17)
NT.PCoA.photosynthetic.UnweightedUniFrac.phyloseq.region.samples
dev.off()

saveRDS(NT.PCoA.photosynthetic.UnweightedUniFrac.phyloseq.region.samples, "NT.PCoA.photosynthetic.UnweightedUniFrac.phyloseq.region.samples.rds")

NT.PCoA.chemolithotrophic.UnweightedUniFrac.phyloseq.region.samples <- plot_ordination(NT.rarefied.chemolithotrophic.min.int.exclude.1.3.rooted.spnet, NT.chemolithotrophic.unweighted.int.unifrac.PCoA, type="split", color = "Temp.adj...C.", shape = "Location.Code", title="PCoA on unweighted UniFrac distance for chemolithotrophic taxa in North Thailand") +   geom_point(size=5, position="jitter") +  geom_text(aes(label = ASVno), size = 6, vjust = 1.5, max.overlaps = 200) + geom_point(size=3, position="jitter")  + clean_background  + scale_shape_manual(values= loc.shape)  + scale_color_manual(values = temp.col) + theme(legend.title = element_text(size = 15.5, face = "bold"), legend.text = element_text(size = 14),legend.position = "bottom", legend.direction = "horizontal",legend.box="vertical", legend.margin=margin()) + guides(colour = guide_legend(nrow=1,byrow=TRUE, override.aes = list(size = 7)), shape = guide_legend(nrow=1,byrow=TRUE,override.aes = list(size = 7)) )

NT.PCoA.chemolithotrophic.UnweightedUniFrac.phyloseq.region.samples$layers <- NT.PCoA.chemolithotrophic.UnweightedUniFrac.phyloseq.region.samples$layers[-1]

pdf("NT.PCoA.chemolithotrophic.UnweightedUniFrac.phyloseq.region.samples.pdf", height = 12, width = 17)
NT.PCoA.chemolithotrophic.UnweightedUniFrac.phyloseq.region.samples
dev.off()

saveRDS(NT.PCoA.chemolithotrophic.UnweightedUniFrac.phyloseq.region.samples, "NT.PCoA.chemolithotrophic.UnweightedUniFrac.phyloseq.region.samples.rds")




NT.PCoA.heterotrophic.UnweightedUniFrac.phyloseq.region.samples <- plot_ordination(NT.rarefied.heterotrophic.min.int.exclude.1.3.rooted.spnet, NT.heterotrophic.unweighted.int.unifrac.PCoA, type="split", color = "Temp.adj...C.", shape = "Location.Code", title="PCoA on unweighted UniFrac distance for heterotrophic taxa in North Thailand") +   geom_point(size=5, position="jitter") +  geom_text(aes(label = ASVno), size = 6, vjust = 1.5, max.overlaps = 200) + geom_point(size=3, position="jitter")  + clean_background  + scale_shape_manual(values= loc.shape)  + scale_color_manual(values = temp.col) + theme(legend.title = element_text(size = 15.5, face = "bold"), legend.text = element_text(size = 14),legend.position = "bottom", legend.direction = "horizontal",legend.box="vertical", legend.margin=margin()) + guides(colour = guide_legend(nrow=1,byrow=TRUE, override.aes = list(size = 7)), shape = guide_legend(nrow=1,byrow=TRUE,override.aes = list(size = 7)) )

NT.PCoA.heterotrophic.UnweightedUniFrac.phyloseq.region.samples$layers <- NT.PCoA.heterotrophic.UnweightedUniFrac.phyloseq.region.samples$layers[-1]

pdf("NT.PCoA.heterotrophic.UnweightedUniFrac.phyloseq.region.samples.pdf", height = 12, width = 17)
NT.PCoA.heterotrophic.UnweightedUniFrac.phyloseq.region.samples
dev.off()

saveRDS(NT.PCoA.heterotrophic.UnweightedUniFrac.phyloseq.region.samples, "NT.PCoA.heterotrophic.UnweightedUniFrac.phyloseq.region.samples.rds")



NM.PCoA.UnweightedUniFrac.phyloseq.region.samples <- plot_ordination(NM.rarefied.min.int.exclude.1.3.rooted.spnet, NM.unweighted.int.unifrac.PCoA, type="split", color = "Temp.adj...C.", shape = "Location.Code", title="PCoA on unweighted UniFrac distance for North Malaysia") +   geom_point(size=5, position="jitter") +  geom_text(aes(label = ASVno), size = 6, vjust = 1.5, max.overlaps = 200) + geom_point(size=3, position="jitter")  + clean_background  + scale_shape_manual(values= loc.shape)  + scale_color_manual(values = temp.col) + theme(legend.title = element_text(size = 15.5, face = "bold"), legend.text = element_text(size = 14),legend.position = "bottom", legend.direction = "horizontal",legend.box="vertical", legend.margin=margin()) + guides(colour = guide_legend(nrow=1,byrow=TRUE, override.aes = list(size = 7)), shape = guide_legend(nrow=1,byrow=TRUE,override.aes = list(size = 7)) )


NM.PCoA.UnweightedUniFrac.phyloseq.region.samples$layers <- NM.PCoA.UnweightedUniFrac.phyloseq.region.samples$layers[-1]

pdf("NM.PCoA.UnweightedUniFrac.phyloseq.region.samples.pdf", height = 12, width = 17)
NM.PCoA.UnweightedUniFrac.phyloseq.region.samples
dev.off()

saveRDS(NM.PCoA.UnweightedUniFrac.phyloseq.region.samples, "NM.PCoA.UnweightedUniFrac.phyloseq.region.samples.rds")


NM.PCoA.photosynthetic.UnweightedUniFrac.phyloseq.region.samples <- plot_ordination(NM.rarefied.photosynthetic.min.int.exclude.1.3.rooted.spnet, NM.photosynthetic.unweighted.int.unifrac.PCoA, type="split", color = "Temp.adj...C.", shape = "Location.Code", title="PCoA on unweighted UniFrac distance for photosynthetic taxa in North Malaysia") +   geom_point(size=5, position="jitter") +  geom_text(aes(label = ASVno), size = 6, vjust = 1.5, max.overlaps = 200) + geom_point(size=3, position="jitter")  + clean_background  + scale_shape_manual(values= loc.shape)  + scale_color_manual(values = temp.col) + theme(legend.title = element_text(size = 15.5, face = "bold"), legend.text = element_text(size = 14),legend.position = "bottom", legend.direction = "horizontal",legend.box="vertical", legend.margin=margin()) + guides(colour = guide_legend(nrow=1,byrow=TRUE, override.aes = list(size = 7)), shape = guide_legend(nrow=1,byrow=TRUE,override.aes = list(size = 7)) )


NM.PCoA.photosynthetic.UnweightedUniFrac.phyloseq.region.samples$layers <- NM.PCoA.photosynthetic.UnweightedUniFrac.phyloseq.region.samples$layers[-1]

pdf("NM.PCoA.photosynthetic.UnweightedUniFrac.phyloseq.region.samples.pdf", height = 12, width = 17)
NM.PCoA.photosynthetic.UnweightedUniFrac.phyloseq.region.samples
dev.off()

saveRDS(NM.PCoA.photosynthetic.UnweightedUniFrac.phyloseq.region.samples, "NM.PCoA.photosynthetic.UnweightedUniFrac.phyloseq.region.samples.rds")

NM.PCoA.chemolithotrophic.UnweightedUniFrac.phyloseq.region.samples <- plot_ordination(NM.rarefied.chemolithotrophic.min.int.exclude.1.3.rooted.spnet, NM.chemolithotrophic.unweighted.int.unifrac.PCoA, type="split", color = "Temp.adj...C.", shape = "Location.Code", title="PCoA on unweighted UniFrac distance for chemolithotrophic taxa in North Malaysia") +   geom_point(size=5, position="jitter") +  geom_text(aes(label = ASVno), size = 6, vjust = 1.5, max.overlaps = 200) + geom_point(size=3, position="jitter")  + clean_background  + scale_shape_manual(values= loc.shape)  + scale_color_manual(values = temp.col) + theme(legend.title = element_text(size = 15.5, face = "bold"), legend.text = element_text(size = 14),legend.position = "bottom", legend.direction = "horizontal",legend.box="vertical", legend.margin=margin()) + guides(colour = guide_legend(nrow=1,byrow=TRUE, override.aes = list(size = 7)), shape = guide_legend(nrow=1,byrow=TRUE,override.aes = list(size = 7)) )

NM.PCoA.chemolithotrophic.UnweightedUniFrac.phyloseq.region.samples$layers <- NM.PCoA.chemolithotrophic.UnweightedUniFrac.phyloseq.region.samples$layers[-1]

pdf("NM.PCoA.chemolithotrophic.UnweightedUniFrac.phyloseq.region.samples.pdf", height = 12, width = 17)
NM.PCoA.chemolithotrophic.UnweightedUniFrac.phyloseq.region.samples
dev.off()

saveRDS(NM.PCoA.chemolithotrophic.UnweightedUniFrac.phyloseq.region.samples, "NM.PCoA.chemolithotrophic.UnweightedUniFrac.phyloseq.region.samples.rds")


NM.PCoA.heterotrophic.UnweightedUniFrac.phyloseq.region.samples <- plot_ordination(NM.rarefied.heterotrophic.min.int.exclude.1.3.rooted.spnet, NM.heterotrophic.unweighted.int.unifrac.PCoA, type="split", color = "Temp.adj...C.", shape = "Location.Code", title="PCoA on unweighted UniFrac distance for heterotrophic taxa in North Malaysia") +   geom_point(size=5, position="jitter") +  geom_text(aes(label = ASVno), size = 6, vjust = 1.5, max.overlaps = 200) + geom_point(size=3, position="jitter")  + clean_background  + scale_shape_manual(values= loc.shape)  + scale_color_manual(values = temp.col) + theme(legend.title = element_text(size = 15.5, face = "bold"), legend.text = element_text(size = 14),legend.position = "bottom", legend.direction = "horizontal",legend.box="vertical", legend.margin=margin()) + guides(colour = guide_legend(nrow=1,byrow=TRUE, override.aes = list(size = 7)), shape = guide_legend(nrow=1,byrow=TRUE,override.aes = list(size = 7)) )

NM.PCoA.heterotrophic.UnweightedUniFrac.phyloseq.region.samples$layers <- NM.PCoA.heterotrophic.UnweightedUniFrac.phyloseq.region.samples$layers[-1]

pdf("NM.PCoA.heterotrophic.UnweightedUniFrac.phyloseq.region.samples.pdf", height = 12, width = 17)
NM.PCoA.heterotrophic.UnweightedUniFrac.phyloseq.region.samples
dev.off()

saveRDS(NM.PCoA.heterotrophic.UnweightedUniFrac.phyloseq.region.samples, "NM.PCoA.heterotrophic.UnweightedUniFrac.phyloseq.region.samples.rds")


ST.PCoA.UnweightedUniFrac.phyloseq.region.samples <- plot_ordination(ST.rarefied.min.int.exclude.1.3.rooted.spnet, ST.unweighted.int.unifrac.PCoA, type="split", color = "Temp.adj...C.", shape = "Location.Code", title="PCoA on unweighted UniFrac distance for South Thailand") +   geom_point(size=5, position="jitter") +  geom_text(aes(label = ASVno), size = 6, vjust = 1.5, max.overlaps = 200) + geom_point(size=3, position="jitter")  + clean_background  + scale_shape_manual(values= loc.shape)  + scale_color_manual(values = temp.col) + theme(legend.title = element_text(size = 15.5, face = "bold"), legend.text = element_text(size = 14),legend.position = "bottom", legend.direction = "horizontal",legend.box="vertical", legend.margin=margin()) + guides(colour = guide_legend(nrow=1,byrow=TRUE, override.aes = list(size = 7)), shape = guide_legend(nrow=1,byrow=TRUE,override.aes = list(size = 7)) )


ST.PCoA.UnweightedUniFrac.phyloseq.region.samples$layers <- ST.PCoA.UnweightedUniFrac.phyloseq.region.samples$layers[-1]

pdf("ST.PCoA.UnweightedUniFrac.phyloseq.region.samples.pdf", height = 12, width = 17)
ST.PCoA.UnweightedUniFrac.phyloseq.region.samples
dev.off()

saveRDS(ST.PCoA.UnweightedUniFrac.phyloseq.region.samples, "ST.PCoA.UnweightedUniFrac.phyloseq.region.samples.rds")


ST.PCoA.photosynthetic.UnweightedUniFrac.phyloseq.region.samples <- plot_ordination(ST.rarefied.photosynthetic.min.int.exclude.1.3.rooted.spnet, ST.photosynthetic.unweighted.int.unifrac.PCoA, type="split", color = "Temp.adj...C.", shape = "Location.Code", title="PCoA on unweighted UniFrac distance for photosynthetic taxa in South Thailand") +   geom_point(size=5, position="jitter") +  geom_text(aes(label = ASVno), size = 6, vjust = 1.5, max.overlaps = 200) + geom_point(size=3, position="jitter")  + clean_background  + scale_shape_manual(values= loc.shape)  + scale_color_manual(values = temp.col) + theme(legend.title = element_text(size = 15.5, face = "bold"), legend.text = element_text(size = 14),legend.position = "bottom", legend.direction = "horizontal",legend.box="vertical", legend.margin=margin()) + guides(colour = guide_legend(nrow=1,byrow=TRUE, override.aes = list(size = 7)), shape = guide_legend(nrow=1,byrow=TRUE,override.aes = list(size = 7)) )


ST.PCoA.photosynthetic.UnweightedUniFrac.phyloseq.region.samples$layers <- ST.PCoA.photosynthetic.UnweightedUniFrac.phyloseq.region.samples$layers[-1]

pdf("ST.PCoA.photosynthetic.UnweightedUniFrac.phyloseq.region.samples.pdf", height = 12, width = 17)
ST.PCoA.photosynthetic.UnweightedUniFrac.phyloseq.region.samples
dev.off()

saveRDS(ST.PCoA.photosynthetic.UnweightedUniFrac.phyloseq.region.samples, "ST.PCoA.photosynthetic.UnweightedUniFrac.phyloseq.region.samples.rds")

ST.PCoA.chemolithotrophic.UnweightedUniFrac.phyloseq.region.samples <- plot_ordination(ST.rarefied.chemolithotrophic.min.int.exclude.1.3.rooted.spnet, ST.chemolithotrophic.unweighted.int.unifrac.PCoA, type="split", color = "Temp.adj...C.", shape = "Location.Code", title="PCoA on unweighted UniFrac distance for chemolithotrophic taxa in South Thailand") +   geom_point(size=5, position="jitter") +  geom_text(aes(label = ASVno), size = 6, vjust = 1.5, max.overlaps = 200) + geom_point(size=3, position="jitter")  + clean_background  + scale_shape_manual(values= loc.shape)  + scale_color_manual(values = temp.col) + theme(legend.title = element_text(size = 15.5, face = "bold"), legend.text = element_text(size = 14),legend.position = "bottom", legend.direction = "horizontal",legend.box="vertical", legend.margin=margin()) + guides(colour = guide_legend(nrow=1,byrow=TRUE, override.aes = list(size = 7)), shape = guide_legend(nrow=1,byrow=TRUE,override.aes = list(size = 7)) )

ST.PCoA.chemolithotrophic.UnweightedUniFrac.phyloseq.region.samples$layers <- ST.PCoA.chemolithotrophic.UnweightedUniFrac.phyloseq.region.samples$layers[-1]

pdf("ST.PCoA.chemolithotrophic.UnweightedUniFrac.phyloseq.region.samples.pdf", height = 12, width = 17)
ST.PCoA.chemolithotrophic.UnweightedUniFrac.phyloseq.region.samples
dev.off()

saveRDS(ST.PCoA.chemolithotrophic.UnweightedUniFrac.phyloseq.region.samples, "ST.PCoA.chemolithotrophic.UnweightedUniFrac.phyloseq.region.samples.rds")


ST.PCoA.heterotrophic.UnweightedUniFrac.phyloseq.region.samples <- plot_ordination(ST.rarefied.heterotrophic.min.int.exclude.1.3.rooted.spnet, ST.heterotrophic.unweighted.int.unifrac.PCoA, type="split", color = "Temp.adj...C.", shape = "Location.Code", title="PCoA on unweighted UniFrac distance for heterotrophic taxa in South Thailand") +   geom_point(size=5, position="jitter") +  geom_text(aes(label = ASVno), size = 6, vjust = 1.5, max.overlaps = 200) + geom_point(size=3, position="jitter")  + clean_background  + scale_shape_manual(values= loc.shape)  + scale_color_manual(values = temp.col) + theme(legend.title = element_text(size = 15.5, face = "bold"), legend.text = element_text(size = 14),legend.position = "bottom", legend.direction = "horizontal",legend.box="vertical", legend.margin=margin()) + guides(colour = guide_legend(nrow=1,byrow=TRUE, override.aes = list(size = 7)), shape = guide_legend(nrow=1,byrow=TRUE,override.aes = list(size = 7)) )

ST.PCoA.heterotrophic.UnweightedUniFrac.phyloseq.region.samples$layers <- ST.PCoA.heterotrophic.UnweightedUniFrac.phyloseq.region.samples$layers[-1]

pdf("ST.PCoA.heterotrophic.UnweightedUniFrac.phyloseq.region.samples.pdf", height = 12, width = 17)
ST.PCoA.heterotrophic.UnweightedUniFrac.phyloseq.region.samples
dev.off()

saveRDS(ST.PCoA.heterotrophic.UnweightedUniFrac.phyloseq.region.samples, "ST.PCoA.heterotrophic.UnweightedUniFrac.phyloseq.region.samples.rds")


SM.PCoA.UnweightedUniFrac.phyloseq.region.samples <- plot_ordination(SM.rarefied.min.int.exclude.1.3.rooted.spnet, SM.unweighted.int.unifrac.PCoA, type="split", color = "Temp.adj...C.", shape = "Location.Code", title="PCoA on unweighted UniFrac distance for SOuth Malaysia") +   geom_point(size=5, position="jitter") +  geom_text(aes(label = ASVno), size = 6, vjust = 1.5, max.overlaps = 200) + geom_point(size=3, position="jitter")  + clean_background  + scale_shape_manual(values= loc.shape)  + scale_color_manual(values = temp.col) + theme(legend.title = element_text(size = 15.5, face = "bold"), legend.text = element_text(size = 14),legend.position = "bottom", legend.direction = "horizontal",legend.box="vertical", legend.margin=margin()) + guides(colour = guide_legend(nrow=1,byrow=TRUE, override.aes = list(size = 7)), shape = guide_legend(nrow=1,byrow=TRUE,override.aes = list(size = 7)) )


SM.PCoA.UnweightedUniFrac.phyloseq.region.samples$layers <- SM.PCoA.UnweightedUniFrac.phyloseq.region.samples$layers[-1]

pdf("SM.PCoA.UnweightedUniFrac.phyloseq.region.samples.pdf", height = 12, width = 17)
SM.PCoA.UnweightedUniFrac.phyloseq.region.samples
dev.off()

saveRDS(SM.PCoA.UnweightedUniFrac.phyloseq.region.samples, "SM.PCoA.UnweightedUniFrac.phyloseq.region.samples.rds")


SM.PCoA.photosynthetic.UnweightedUniFrac.phyloseq.region.samples <- plot_ordination(SM.rarefied.photosynthetic.min.int.exclude.1.3.rooted.spnet, SM.photosynthetic.unweighted.int.unifrac.PCoA, type="split", color = "Temp.adj...C.", shape = "Location.Code", title="PCoA on unweighted UniFrac distance for photosynthetic taxa in South Malaysia") +   geom_point(size=5, position="jitter") +  geom_text(aes(label = ASVno), size = 6, vjust = 1.5, max.overlaps = 200) + geom_point(size=3, position="jitter")  + clean_background  + scale_shape_manual(values= loc.shape)  + scale_color_manual(values = temp.col) + theme(legend.title = element_text(size = 15.5, face = "bold"), legend.text = element_text(size = 14),legend.position = "bottom", legend.direction = "horizontal",legend.box="vertical", legend.margin=margin()) + guides(colour = guide_legend(nrow=1,byrow=TRUE, override.aes = list(size = 7)), shape = guide_legend(nrow=1,byrow=TRUE,override.aes = list(size = 7)) )


SM.PCoA.photosynthetic.UnweightedUniFrac.phyloseq.region.samples$layers <- SM.PCoA.photosynthetic.UnweightedUniFrac.phyloseq.region.samples$layers[-1]

pdf("SM.PCoA.photosynthetic.UnweightedUniFrac.phyloseq.region.samples.pdf", height = 12, width = 17)
SM.PCoA.photosynthetic.UnweightedUniFrac.phyloseq.region.samples
dev.off()

saveRDS(SM.PCoA.photosynthetic.UnweightedUniFrac.phyloseq.region.samples, "SM.PCoA.photosynthetic.UnweightedUniFrac.phyloseq.region.samples.rds")

SM.PCoA.chemolithotrophic.UnweightedUniFrac.phyloseq.region.samples <- plot_ordination(SM.rarefied.chemolithotrophic.min.int.exclude.1.3.rooted.spnet, SM.chemolithotrophic.unweighted.int.unifrac.PCoA, type="split", color = "Temp.adj...C.", shape = "Location.Code", title="PCoA on unweighted UniFrac distance for chemolithotrophic taxa in South Malaysia") +   geom_point(size=5, position="jitter") +  geom_text(aes(label = ASVno), size = 6, vjust = 1.5, max.overlaps = 200) + geom_point(size=3, position="jitter")  + clean_background  + scale_shape_manual(values= loc.shape)  + scale_color_manual(values = temp.col) + theme(legend.title = element_text(size = 15.5, face = "bold"), legend.text = element_text(size = 14),legend.position = "bottom", legend.direction = "horizontal",legend.box="vertical", legend.margin=margin()) + guides(colour = guide_legend(nrow=1,byrow=TRUE, override.aes = list(size = 7)), shape = guide_legend(nrow=1,byrow=TRUE,override.aes = list(size = 7)) )

SM.PCoA.chemolithotrophic.UnweightedUniFrac.phyloseq.region.samples$layers <- SM.PCoA.chemolithotrophic.UnweightedUniFrac.phyloseq.region.samples$layers[-1]

pdf("SM.PCoA.chemolithotrophic.UnweightedUniFrac.phyloseq.region.samples.pdf", height = 12, width = 17)
SM.PCoA.chemolithotrophic.UnweightedUniFrac.phyloseq.region.samples
dev.off()

saveRDS(SM.PCoA.chemolithotrophic.UnweightedUniFrac.phyloseq.region.samples, "SM.PCoA.chemolithotrophic.UnweightedUniFrac.phyloseq.region.samples.rds")


SM.PCoA.heterotrophic.UnweightedUniFrac.phyloseq.region.samples <- plot_ordination(SM.rarefied.heterotrophic.min.int.exclude.1.3.rooted.spnet, SM.heterotrophic.unweighted.int.unifrac.PCoA, type="split", color = "Temp.adj...C.", shape = "Location.Code", title="PCoA on unweighted UniFrac distance for heterotrophic taxa in South Malaysia") +   geom_point(size=5, position="jitter") +  geom_text(aes(label = ASVno), size = 6, vjust = 1.5, max.overlaps = 200) + geom_point(size=3, position="jitter")  + clean_background  + scale_shape_manual(values= loc.shape)  + scale_color_manual(values = temp.col) + theme(legend.title = element_text(size = 15.5, face = "bold"), legend.text = element_text(size = 14),legend.position = "bottom", legend.direction = "horizontal",legend.box="vertical", legend.margin=margin()) + guides(colour = guide_legend(nrow=1,byrow=TRUE, override.aes = list(size = 7)), shape = guide_legend(nrow=1,byrow=TRUE,override.aes = list(size = 7)) )

SM.PCoA.heterotrophic.UnweightedUniFrac.phyloseq.region.samples$layers <- SM.PCoA.heterotrophic.UnweightedUniFrac.phyloseq.region.samples$layers[-1]

pdf("SM.PCoA.heterotrophic.UnweightedUniFrac.phyloseq.region.samples.pdf", height = 12, width = 17)
SM.PCoA.heterotrophic.UnweightedUniFrac.phyloseq.region.samples
dev.off()

saveRDS(SM.PCoA.heterotrophic.UnweightedUniFrac.phyloseq.region.samples, "SM.PCoA.heterotrophic.UnweightedUniFrac.phyloseq.region.samples.rds")



SG.PCoA.UnweightedUniFrac.phyloseq.region.samples <- plot_ordination(SG.rarefied.min.int.exclude.1.3.rooted.spnet, SG.unweighted.int.unifrac.PCoA, type="split", color = "Temp.adj...C.", shape = "Location.Code", title="PCoA on unweighted UniFrac distance for Singapore") +   geom_point(size=5, position="jitter") +  geom_text(aes(label = ASVno), size = 6, vjust = 1.5, max.overlaps = 200) + geom_point(size=3, position="jitter")  + clean_background  + scale_shape_manual(values= loc.shape)  + scale_color_manual(values = temp.col) + theme(legend.title = element_text(size = 15.5, face = "bold"), legend.text = element_text(size = 14),legend.position = "bottom", legend.direction = "horizontal",legend.box="vertical", legend.margin=margin()) + guides(colour = guide_legend(nrow=1,byrow=TRUE, override.aes = list(size = 7)), shape = guide_legend(nrow=1,byrow=TRUE,override.aes = list(size = 7)) )


SG.PCoA.UnweightedUniFrac.phyloseq.region.samples$layers <- SG.PCoA.UnweightedUniFrac.phyloseq.region.samples$layers[-1]

pdf("SG.PCoA.UnweightedUniFrac.phyloseq.region.samples.pdf", height = 12, width = 17)
SG.PCoA.UnweightedUniFrac.phyloseq.region.samples
dev.off()

saveRDS(SG.PCoA.UnweightedUniFrac.phyloseq.region.samples, "SG.PCoA.UnweightedUniFrac.phyloseq.region.samples.rds")


SG.PCoA.photosynthetic.UnweightedUniFrac.phyloseq.region.samples <- plot_ordination(SG.rarefied.photosynthetic.min.int.exclude.1.3.rooted.spnet, SG.photosynthetic.unweighted.int.unifrac.PCoA, type="split", color = "Temp.adj...C.", shape = "Location.Code", title="PCoA on unweighted UniFrac distance for photosynthetic taxa in Singapore") +   geom_point(size=5, position="jitter") +  geom_text(aes(label = ASVno), size = 6, vjust = 1.5, max.overlaps = 200) + geom_point(size=3, position="jitter")  + clean_background  + scale_shape_manual(values= loc.shape)  + scale_color_manual(values = temp.col) + theme(legend.title = element_text(size = 15.5, face = "bold"), legend.text = element_text(size = 14),legend.position = "bottom", legend.direction = "horizontal",legend.box="vertical", legend.margin=margin()) + guides(colour = guide_legend(nrow=1,byrow=TRUE, override.aes = list(size = 7)), shape = guide_legend(nrow=1,byrow=TRUE,override.aes = list(size = 7)) )


SG.PCoA.photosynthetic.UnweightedUniFrac.phyloseq.region.samples$layers <- SG.PCoA.photosynthetic.UnweightedUniFrac.phyloseq.region.samples$layers[-1]

pdf("SG.PCoA.photosynthetic.UnweightedUniFrac.phyloseq.region.samples.pdf", height = 12, width = 17)
SG.PCoA.photosynthetic.UnweightedUniFrac.phyloseq.region.samples
dev.off()

saveRDS(SG.PCoA.photosynthetic.UnweightedUniFrac.phyloseq.region.samples, "SG.PCoA.photosynthetic.UnweightedUniFrac.phyloseq.region.samples.rds")

### NO ORDINATION FOR SG CHEMOLITHS> ONLY AXIS1 is computed ####
#SG.PCoA.chemolithotrophic.UnweightedUniFrac.phyloseq.region.samples <- plot_ordination(SG.rarefied.chemolithotrophic.min.int.exclude.1.3.rooted.spnet, SG.chemolithotrophic.unweighted.int.unifrac.PCoA, type="split", color = "Temp.adj...C.", shape = "Location.Code", title="PCoA on unweighted UniFrac distance for chemolithotrophic taxa in Singapore") +   geom_point(size=5, position="jitter") +  geom_text(aes(label = ASVno), size = 6, vjust = 1.5, max.overlaps = 200) + geom_point(size=3, position="jitter")  + clean_background  + scale_shape_manual(values= loc.shape)  + scale_color_manual(values = temp.col) + theme(legend.title = element_text(size = 15.5, face = "bold"), legend.text = element_text(size = 14),legend.position = "bottom", legend.direction = "horizontal",legend.box="vertical", legend.margin=margin()) + guides(colour = guide_legend(nrow=1,byrow=TRUE, override.aes = list(size = 7)), shape = guide_legend(nrow=1,byrow=TRUE,override.aes = list(size = 7)) )

#SG.PCoA.chemolithotrophic.UnweightedUniFrac.phyloseq.region.samples$layers <- SG.PCoA.chemolithotrophic.UnweightedUniFrac.phyloseq.region.samples$layers[-1]

#pdf("SG.PCoA.chemolithotrophic.UnweightedUniFrac.phyloseq.region.samples.pdf", height = 12, width = 17)
#SG.PCoA.chemolithotrophic.UnweightedUniFrac.phyloseq.region.samples
#dev.off()

#saveRDS(SG.PCoA.chemolithotrophic.UnweightedUniFrac.phyloseq.region.samples, "SG.PCoA.chemolithotrophic.UnweightedUniFrac.phyloseq.region.samples.rds")


SG.PCoA.heterotrophic.UnweightedUniFrac.phyloseq.region.samples <- plot_ordination(SG.rarefied.heterotrophic.min.int.exclude.1.3.rooted.spnet, SG.heterotrophic.unweighted.int.unifrac.PCoA, type="split", color = "Temp.adj...C.", shape = "Location.Code", title="PCoA on unweighted UniFrac distance for heterotrophic taxa in Singapore") +   geom_point(size=5, position="jitter") +  geom_text(aes(label = ASVno), size = 6, vjust = 1.5, max.overlaps = 200) + geom_point(size=3, position="jitter")  + clean_background  + scale_shape_manual(values= loc.shape)  + scale_color_manual(values = temp.col) + theme(legend.title = element_text(size = 15.5, face = "bold"), legend.text = element_text(size = 14),legend.position = "bottom", legend.direction = "horizontal",legend.box="vertical", legend.margin=margin()) + guides(colour = guide_legend(nrow=1,byrow=TRUE, override.aes = list(size = 7)), shape = guide_legend(nrow=1,byrow=TRUE,override.aes = list(size = 7)) )

SG.PCoA.heterotrophic.UnweightedUniFrac.phyloseq.region.samples$layers <- SG.PCoA.heterotrophic.UnweightedUniFrac.phyloseq.region.samples$layers[-1]

pdf("SG.PCoA.heterotrophic.UnweightedUniFrac.phyloseq.region.samples.pdf", height = 12, width = 17)
SG.PCoA.heterotrophic.UnweightedUniFrac.phyloseq.region.samples
dev.off()

saveRDS(SG.PCoA.heterotrophic.UnweightedUniFrac.phyloseq.region.samples, "SG.PCoA.heterotrophic.UnweightedUniFrac.phyloseq.region.samples.rds")



CT.PCoA.UnweightedUniFrac.phyloseq.region.samples <- plot_ordination(CT.rarefied.min.int.exclude.1.3.rooted.spnet, CT.unweighted.int.unifrac.PCoA, type="split", color = "Temp.adj...C.", shape = "Location.Code", title="PCoA on unweighted UniFrac distance for Central Thailand") +   geom_point(size=5, position="jitter") +  geom_text(aes(label = ASVno), size = 6, vjust = 1.5, max.overlaps = 200) + geom_point(size=3, position="jitter")  + clean_background  + scale_shape_manual(values= loc.shape)  + scale_color_manual(values = temp.col) + theme(legend.title = element_text(size = 15.5, face = "bold"), legend.text = element_text(size = 14),legend.position = "bottom", legend.direction = "horizontal",legend.box="vertical", legend.margin=margin()) + guides(colour = guide_legend(nrow=1,byrow=TRUE, override.aes = list(size = 7)), shape = guide_legend(nrow=1,byrow=TRUE,override.aes = list(size = 7)) )


CT.PCoA.UnweightedUniFrac.phyloseq.region.samples$layers <- CT.PCoA.UnweightedUniFrac.phyloseq.region.samples$layers[-1]

pdf("CT.PCoA.UnweightedUniFrac.phyloseq.region.samples.pdf", height = 12, width = 17)
CT.PCoA.UnweightedUniFrac.phyloseq.region.samples
dev.off()

saveRDS(CT.PCoA.UnweightedUniFrac.phyloseq.region.samples, "CT.PCoA.UnweightedUniFrac.phyloseq.region.samples.rds")


CT.PCoA.photosynthetic.UnweightedUniFrac.phyloseq.region.samples <- plot_ordination(CT.rarefied.photosynthetic.min.int.exclude.1.3.rooted.spnet, CT.photosynthetic.unweighted.int.unifrac.PCoA, type="split", color = "Temp.adj...C.", shape = "Location.Code", title="PCoA on unweighted UniFrac distance for photosynthetic taxa in Central Thailand") +   geom_point(size=5, position="jitter") +  geom_text(aes(label = ASVno), size = 6, vjust = 1.5, max.overlaps = 200) + geom_point(size=3, position="jitter")  + clean_background  + scale_shape_manual(values= loc.shape)  + scale_color_manual(values = temp.col) + theme(legend.title = element_text(size = 15.5, face = "bold"), legend.text = element_text(size = 14),legend.position = "bottom", legend.direction = "horizontal",legend.box="vertical", legend.margin=margin()) + guides(colour = guide_legend(nrow=1,byrow=TRUE, override.aes = list(size = 7)), shape = guide_legend(nrow=1,byrow=TRUE,override.aes = list(size = 7)) )


CT.PCoA.photosynthetic.UnweightedUniFrac.phyloseq.region.samples$layers <- CT.PCoA.photosynthetic.UnweightedUniFrac.phyloseq.region.samples$layers[-1]

pdf("CT.PCoA.photosynthetic.UnweightedUniFrac.phyloseq.region.samples.pdf", height = 12, width = 17)
CT.PCoA.photosynthetic.UnweightedUniFrac.phyloseq.region.samples
dev.off()

saveRDS(CT.PCoA.photosynthetic.UnweightedUniFrac.phyloseq.region.samples, "CT.PCoA.photosynthetic.UnweightedUniFrac.phyloseq.region.samples.rds")

CT.PCoA.chemolithotrophic.UnweightedUniFrac.phyloseq.region.samples <- plot_ordination(CT.rarefied.chemolithotrophic.min.int.exclude.1.3.rooted.spnet, CT.chemolithotrophic.unweighted.int.unifrac.PCoA, type="split", color = "Temp.adj...C.", shape = "Location.Code", title="PCoA on unweighted UniFrac distance for chemolithotrophic taxa in Central Thailand") +   geom_point(size=5, position="jitter") +  geom_text(aes(label = ASVno), size = 6, vjust = 1.5, max.overlaps = 200) + geom_point(size=3, position="jitter")  + clean_background  + scale_shape_manual(values= loc.shape)  + scale_color_manual(values = temp.col) + theme(legend.title = element_text(size = 15.5, face = "bold"), legend.text = element_text(size = 14),legend.position = "bottom", legend.direction = "horizontal",legend.box="vertical", legend.margin=margin()) + guides(colour = guide_legend(nrow=1,byrow=TRUE, override.aes = list(size = 7)), shape = guide_legend(nrow=1,byrow=TRUE,override.aes = list(size = 7)) )

CT.PCoA.chemolithotrophic.UnweightedUniFrac.phyloseq.region.samples$layers <- CT.PCoA.chemolithotrophic.UnweightedUniFrac.phyloseq.region.samples$layers[-1]

pdf("CT.PCoA.chemolithotrophic.UnweightedUniFrac.phyloseq.region.samples.pdf", height = 12, width = 17)
CT.PCoA.chemolithotrophic.UnweightedUniFrac.phyloseq.region.samples
dev.off()

saveRDS(CT.PCoA.chemolithotrophic.UnweightedUniFrac.phyloseq.region.samples, "CT.PCoA.chemolithotrophic.UnweightedUniFrac.phyloseq.region.samples.rds")


CT.PCoA.heterotrophic.UnweightedUniFrac.phyloseq.region.samples <- plot_ordination(CT.rarefied.heterotrophic.min.int.exclude.1.3.rooted.spnet, CT.heterotrophic.unweighted.int.unifrac.PCoA, type="split", color = "Temp.adj...C.", shape = "Location.Code", title="PCoA on unweighted UniFrac distance for heterotrophic taxa in Central Thailand") +   geom_point(size=5, position="jitter") +  geom_text(aes(label = ASVno), size = 6, vjust = 1.5, max.overlaps = 200) + geom_point(size=3, position="jitter")  + clean_background  + scale_shape_manual(values= loc.shape)  + scale_color_manual(values = temp.col) + theme(legend.title = element_text(size = 15.5, face = "bold"), legend.text = element_text(size = 14),legend.position = "bottom", legend.direction = "horizontal",legend.box="vertical", legend.margin=margin()) + guides(colour = guide_legend(nrow=1,byrow=TRUE, override.aes = list(size = 7)), shape = guide_legend(nrow=1,byrow=TRUE,override.aes = list(size = 7)) )

CT.PCoA.heterotrophic.UnweightedUniFrac.phyloseq.region.samples$layers <- CT.PCoA.heterotrophic.UnweightedUniFrac.phyloseq.region.samples$layers[-1]

pdf("CT.PCoA.heterotrophic.UnweightedUniFrac.phyloseq.region.samples.pdf", height = 12, width = 17)
CT.PCoA.heterotrophic.UnweightedUniFrac.phyloseq.region.samples
dev.off()

saveRDS(CT.PCoA.heterotrophic.UnweightedUniFrac.phyloseq.region.samples, "CT.PCoA.heterotrophic.UnweightedUniFrac.phyloseq.region.samples.rds")

### prepare dataframes.

rf.data.PCOA.all.reg <- rarefied.PCoA.UnweightedUniFrac.phyloseq.regionellipse.samples$data %>% select("Axis.1","Region","Carbonate..mg.L.", "pH", "Phosphate..mg.L.", "H2S..mg.L.", "Temp...C." , "Cond..mS.")

colnames(rf.data.PCOA.all.reg) <- c("Axis.1","Region","Carbonate..ppm.", "pH", "Phosphate..ppm.", "H2S..ppm.", "Temp...C." , "EC..mS.")

rf.data.PCOA.all.lat <- rarefied.PCoA.UnweightedUniFrac.phyloseq.regionellipse.samples$data %>% select("Axis.1","Latittude","Carbonate..mg.L.", "pH", "Phosphate..mg.L.", "H2S..mg.L.", "Temp...C." , "Cond..mS.")

colnames(rf.data.PCOA.all.lat) <- c("Axis.1","Latitude","Carbonate..ppm.", "pH", "Phosphate..ppm.", "H2S..ppm.", "Temp...C." , "EC..mS.")

rf.data.PCOA.select <- rarefied.PCoA.UnweightedUniFrac.phyloseq.regionellipse.samples$data %>% select("Axis.1","Carbonate..mg.L.", "pH", "Phosphate..mg.L.", "H2S..mg.L.", "Temp...C." , "Cond..mS.")

colnames(rf.data.PCOA.select) <- c("Axis.1","Carbonate..ppm.", "pH", "Phosphate..ppm.", "H2S..ppm.", "Temp...C." , "EC..mS.")


rf.data.PCOA.photo.reg <- rarefied.PCoA.UnweightedUniFrac.phyloseq.photosynthetic$data %>% filter(id.type == "Samples") %>% select("Axis.1","Region","Carbonate..mg.L.", "pH", "Phosphate..mg.L.", "H2S..mg.L.", "Temp...C." , "Cond..mS.")

colnames(rf.data.PCOA.photo.reg) <- c("Axis.1","Region","Carbonate..ppm.", "pH", "Phosphate..ppm.", "H2S..ppm.", "Temp...C." , "EC..mS.")

rf.data.PCOA.photo.lat <- rarefied.PCoA.UnweightedUniFrac.phyloseq.photosynthetic$data %>% filter(id.type == "Samples") %>% select("Axis.1","Latittude","Carbonate..mg.L.", "pH", "Phosphate..mg.L.", "H2S..mg.L.", "Temp...C." , "Cond..mS.")

colnames(rf.data.PCOA.photo.lat) <- c("Axis.1","Latitude","Carbonate..ppm.", "pH", "Phosphate..ppm.", "H2S..ppm.", "Temp...C." , "EC..mS.")


rf.data.PCOA.photo.select <- rarefied.PCoA.UnweightedUniFrac.phyloseq.photosynthetic$data %>% filter(id.type == "Samples") %>% select("Axis.1","Carbonate..mg.L.", "pH", "Phosphate..mg.L.", "H2S..mg.L.", "Temp...C." , "Cond..mS.")

colnames(rf.data.PCOA.photo.select) <- c("Axis.1","Carbonate..ppm.", "pH", "Phosphate..ppm.", "H2S..ppm.", "Temp...C." , "EC..mS.")


rf.data.PCOA.chemo.reg <- rarefied.PCoA.UnweightedUniFrac.phyloseq.chemolithotrophs$data %>% filter(id.type == "Samples") %>% select("Axis.1","Region","Carbonate..mg.L.", "pH", "Phosphate..mg.L.", "H2S..mg.L.", "Temp...C." , "Cond..mS.")

colnames(rf.data.PCOA.chemo.reg) <- c("Axis.1","Region","Carbonate..ppm.", "pH", "Phosphate..ppm.", "H2S..ppm.", "Temp...C." , "EC..mS.")

rf.data.PCOA.chemo.lat <- rarefied.PCoA.UnweightedUniFrac.phyloseq.chemolithotrophs$data %>% filter(id.type == "Samples") %>% select("Axis.1","Latittude","Carbonate..mg.L.", "pH", "Phosphate..mg.L.", "H2S..mg.L.", "Temp...C." , "Cond..mS.")

colnames(rf.data.PCOA.chemo.lat) <- c("Axis.1","Latitude","Carbonate..ppm.", "pH", "Phosphate..ppm.", "H2S..ppm.", "Temp...C." , "EC..mS.")

rf.data.PCOA.chemo.select <- rarefied.PCoA.UnweightedUniFrac.phyloseq.chemolithotrophs$data %>% filter(id.type == "Samples") %>% select("Axis.1","Carbonate..mg.L.", "pH", "Phosphate..mg.L.", "H2S..mg.L.", "Temp...C." , "Cond..mS.")

colnames(rf.data.PCOA.chemo.select) <- c("Axis.1","Carbonate..ppm.", "pH", "Phosphate..ppm.", "H2S..ppm.", "Temp...C." , "EC..mS.")



rf.data.PCOA.hetero.reg <- rarefied.PCoA.UnweightedUniFrac.phyloseq.only.heterotrophs$data %>% filter(id.type == "Samples") %>% select("Axis.1","Region","Carbonate..mg.L.", "pH", "Phosphate..mg.L.", "H2S..mg.L.", "Temp...C." , "Cond..mS.")

colnames(rf.data.PCOA.hetero.reg) <- c("Axis.1","Region","Carbonate..ppm.", "pH", "Phosphate..ppm.", "H2S..ppm.", "Temp...C." , "EC..mS.")


rf.data.PCOA.hetero.lat <- rarefied.PCoA.UnweightedUniFrac.phyloseq.only.heterotrophs$data %>% filter(id.type == "Samples") %>% select("Axis.1","Latittude","Carbonate..mg.L.", "pH", "Phosphate..mg.L.", "H2S..mg.L.", "Temp...C." , "Cond..mS.")

colnames(rf.data.PCOA.hetero.lat) <- c("Axis.1","Latitude","Carbonate..ppm.", "pH", "Phosphate..ppm.", "H2S..ppm.", "Temp...C." , "EC..mS.")

rf.data.PCOA.hetero.select <- rarefied.PCoA.UnweightedUniFrac.phyloseq.only.heterotrophs$data %>% filter(id.type == "Samples") %>% select("Axis.1","Carbonate..mg.L.", "pH", "Phosphate..mg.L.", "H2S..mg.L.", "Temp...C." , "Cond..mS.")

colnames(rf.data.PCOA.hetero.select) <- c("Axis.1","Carbonate..ppm.", "pH", "Phosphate..ppm.", "H2S..ppm.", "Temp...C." , "EC..mS.")



saveRDS(rf.data.PCOA.all.reg,"rf.data.PCOA.all.reg.rds")
saveRDS(rf.data.PCOA.all.lat,"rf.data.PCOA.all.lat.rds")
saveRDS(rf.data.PCOA.select, "rf.data.PCOA.select.rds")

saveRDS(rf.data.PCOA.photo.reg,"rf.data.PCOA.photo.reg.rds")
saveRDS(rf.data.PCOA.photo.lat,"rf.data.PCOA.photo.lat.rds")
saveRDS(rf.data.PCOA.photo.select,"rf.data.PCOA.photo.select.rds")

saveRDS(rf.data.PCOA.chemo.reg,"rf.data.PCOA.chemo.reg.rds")
saveRDS(rf.data.PCOA.chemo.lat,"rf.data.PCOA.chemo.lat.rds")
saveRDS(rf.data.PCOA.chemo.select,"rf.data.PCOA.chemo.select.rds")

saveRDS(rf.data.PCOA.hetero.reg,"rf.data.PCOA.hetero.reg.rds")
saveRDS(rf.data.PCOA.hetero.lat,"rf.data.PCOA.hetero.lat.rds")
saveRDS(rf.data.PCOA.hetero.select,"rf.data.PCOA.hetero.select.rds")


rf.data.PCOA.NT.all <- NT.PCoA.UnweightedUniFrac.phyloseq.region.samples$data %>% filter(id.type == "Samples") %>% select("Axis.1","Region","Latittude","Carbonate..ppm.", "pH", "Phosphate..ppm.", "H2S..ppm.", "Temp...C." , "EC..mS.")

rf.data.PCOA.NT.reg <- NT.PCoA.UnweightedUniFrac.phyloseq.region.samples$data %>% filter(id.type == "Samples") %>% select("Axis.1","Region","Carbonate..ppm.", "pH", "Phosphate..ppm.", "H2S..ppm.", "Temp...C." , "EC..mS.")

rf.data.PCOA.NT.lat <- NT.PCoA.UnweightedUniFrac.phyloseq.region.samples$data  %>% filter(id.type == "Samples") %>% select("Axis.1","Latittude","Carbonate..ppm.", "pH", "Phosphate..ppm.", "H2S..ppm.", "Temp...C." , "EC..mS.")

rf.data.PCOA.NT <- NT.PCoA.UnweightedUniFrac.phyloseq.region.samples$data  %>% filter(id.type == "Samples") %>% select("Axis.1","Carbonate..ppm.", "pH", "Phosphate..ppm.", "H2S..ppm.", "Temp...C." , "EC..mS.")

rf.data.PCOA.NT.photo <- NT.PCoA.photosynthetic.UnweightedUniFrac.phyloseq.region.samples$data  %>% filter(id.type == "Samples") %>% select("Axis.1","Carbonate..ppm.", "pH", "Phosphate..ppm.", "H2S..ppm.", "Temp...C." , "EC..mS.")

rf.data.PCOA.NT.chemo <- NT.PCoA.chemolithotrophic.UnweightedUniFrac.phyloseq.region.samples$data %>% filter(id.type == "Samples") %>% select("Axis.1","Carbonate..ppm.", "pH", "Phosphate..ppm.", "H2S..ppm.", "Temp...C." , "EC..mS.")

rf.data.PCOA.NT.hetero <- NT.PCoA.chemolithotrophic.UnweightedUniFrac.phyloseq.region.samples$data  %>% filter(id.type == "Samples") %>% select("Axis.1","Carbonate..ppm.", "pH", "Phosphate..ppm.", "H2S..ppm.", "Temp...C." , "EC..mS.")


rf.data.PCOA.NM <- NM.PCoA.UnweightedUniFrac.phyloseq.region.samples$data  %>% filter(id.type == "Samples") %>% select("Axis.1","Carbonate..ppm.", "pH", "Phosphate..ppm.", "H2S..ppm.", "Temp...C." , "EC..mS.")

rf.data.PCOA.NM.photo <- NM.PCoA.photosynthetic.UnweightedUniFrac.phyloseq.region.samples$data  %>% filter(id.type == "Samples") %>% select("Axis.1","Carbonate..ppm.", "pH", "Phosphate..ppm.", "H2S..ppm.", "Temp...C." , "EC..mS.")

rf.data.PCOA.NM.chemo <- NM.PCoA.chemolithotrophic.UnweightedUniFrac.phyloseq.region.samples$data %>% filter(id.type == "Samples") %>% select("Axis.1","Carbonate..ppm.", "pH", "Phosphate..ppm.", "H2S..ppm.", "Temp...C." , "EC..mS.")

rf.data.PCOA.NM.hetero <- NM.PCoA.chemolithotrophic.UnweightedUniFrac.phyloseq.region.samples$data  %>% filter(id.type == "Samples") %>% select("Axis.1","Carbonate..ppm.", "pH", "Phosphate..ppm.", "H2S..ppm.", "Temp...C." , "EC..mS.")


rf.data.PCOA.ST <- ST.PCoA.UnweightedUniFrac.phyloseq.region.samples$data  %>% filter(id.type == "Samples") %>% select("Axis.1","Carbonate..ppm.", "pH", "Phosphate..ppm.", "H2S..ppm.", "Temp...C." , "EC..mS.")

rf.data.PCOA.ST.photo <- ST.PCoA.photosynthetic.UnweightedUniFrac.phyloseq.region.samples$data  %>% filter(id.type == "Samples") %>% select("Axis.1","Carbonate..ppm.", "pH", "Phosphate..ppm.", "H2S..ppm.", "Temp...C." , "EC..mS.")

rf.data.PCOA.ST.chemo <- ST.PCoA.chemolithotrophic.UnweightedUniFrac.phyloseq.region.samples$data %>% filter(id.type == "Samples") %>% select("Axis.1","Carbonate..ppm.", "pH", "Phosphate..ppm.", "H2S..ppm.", "Temp...C." , "EC..mS.")

rf.data.PCOA.ST.hetero <- ST.PCoA.chemolithotrophic.UnweightedUniFrac.phyloseq.region.samples$data  %>% filter(id.type == "Samples") %>% select("Axis.1","Carbonate..ppm.", "pH", "Phosphate..ppm.", "H2S..ppm.", "Temp...C." , "EC..mS.")


rf.data.PCOA.SM <- SM.PCoA.UnweightedUniFrac.phyloseq.region.samples$data  %>% filter(id.type == "Samples") %>% select("Axis.1","Carbonate..ppm.", "pH", "Phosphate..ppm.", "H2S..ppm.", "Temp...C." , "EC..mS.")

rf.data.PCOA.SM.photo <- SM.PCoA.photosynthetic.UnweightedUniFrac.phyloseq.region.samples$data  %>% filter(id.type == "Samples") %>% select("Axis.1","Carbonate..ppm.", "pH", "Phosphate..ppm.", "H2S..ppm.", "Temp...C." , "EC..mS.")

rf.data.PCOA.SM.chemo <- SM.PCoA.chemolithotrophic.UnweightedUniFrac.phyloseq.region.samples$data %>% filter(id.type == "Samples") %>% select("Axis.1","Carbonate..ppm.", "pH", "Phosphate..ppm.", "H2S..ppm.", "Temp...C." , "EC..mS.")

rf.data.PCOA.SM.hetero <- SM.PCoA.chemolithotrophic.UnweightedUniFrac.phyloseq.region.samples$data  %>% filter(id.type == "Samples") %>% select("Axis.1","Carbonate..ppm.", "pH", "Phosphate..ppm.", "H2S..ppm.", "Temp...C." , "EC..mS.")


rf.data.PCOA.CT <- CT.PCoA.UnweightedUniFrac.phyloseq.region.samples$data  %>% filter(id.type == "Samples") %>% select("Axis.1","Carbonate..ppm.", "pH", "Phosphate..ppm.", "H2S..ppm.", "Temp...C." , "EC..mS.")

rf.data.PCOA.CT.photo <- CT.PCoA.photosynthetic.UnweightedUniFrac.phyloseq.region.samples$data  %>% filter(id.type == "Samples") %>% select("Axis.1","Carbonate..ppm.", "pH", "Phosphate..ppm.", "H2S..ppm.", "Temp...C." , "EC..mS.")

rf.data.PCOA.CT.chemo <- CT.PCoA.chemolithotrophic.UnweightedUniFrac.phyloseq.region.samples$data %>% filter(id.type == "Samples") %>% select("Axis.1","Carbonate..ppm.", "pH", "Phosphate..ppm.", "H2S..ppm.", "Temp...C." , "EC..mS.")

rf.data.PCOA.CT.hetero <- CT.PCoA.chemolithotrophic.UnweightedUniFrac.phyloseq.region.samples$data  %>% filter(id.type == "Samples") %>% select("Axis.1","Carbonate..ppm.", "pH", "Phosphate..ppm.", "H2S..ppm.", "Temp...C." , "EC..mS.")


rf.data.PCOA.SG <- SG.PCoA.UnweightedUniFrac.phyloseq.region.samples$data  %>% filter(id.type == "Samples") %>% select("Axis.1","Carbonate..ppm.", "pH", "Phosphate..ppm.", "H2S..ppm.", "Temp...C." , "EC..mS.")

rf.data.PCOA.SG.photo <- SG.PCoA.photosynthetic.UnweightedUniFrac.phyloseq.region.samples$data  %>% filter(id.type == "Samples") %>% select("Axis.1","Carbonate..ppm.", "pH", "Phosphate..ppm.", "H2S..ppm.", "Temp...C." , "EC..mS.")

rf.data.PCOA.SG.hetero <- SG.PCoA.heterotrophic.UnweightedUniFrac.phyloseq.region.samples$data  %>% filter(id.type == "Samples") %>% select("Axis.1","Carbonate..ppm.", "pH", "Phosphate..ppm.", "H2S..ppm.", "Temp...C." , "EC..mS.")


saveRDS(rf.data.PCOA.NT, "rf.data.PCOA.NT.rds")
saveRDS(rf.data.PCOA.NT.photo, "rf.data.PCOA.NT.photo.rds")
saveRDS(rf.data.PCOA.NT.chemo, "rf.data.PCOA.NT.chemo.rds")
saveRDS(rf.data.PCOA.NT.hetero, "rf.data.PCOA.NT.hetero.rds")

saveRDS(rf.data.PCOA.NM, "rf.data.PCOA.NM.rds")
saveRDS(rf.data.PCOA.NM.photo, "rf.data.PCOA.NM.photo.rds")
saveRDS(rf.data.PCOA.NM.chemo, "rf.data.PCOA.NM.chemo.rds")
saveRDS(rf.data.PCOA.NM.hetero, "rf.data.PCOA.NM.hetero.rds")

saveRDS(rf.data.PCOA.ST, "rf.data.PCOA.ST.rds")
saveRDS(rf.data.PCOA.ST.photo, "rf.data.PCOA.ST.photo.rds")
saveRDS(rf.data.PCOA.ST.chemo, "rf.data.PCOA.ST.chemo.rds")
saveRDS(rf.data.PCOA.ST.hetero, "rf.data.PCOA.ST.hetero.rds")

saveRDS(rf.data.PCOA.SM, "rf.data.PCOA.SM.rds")
saveRDS(rf.data.PCOA.SM.photo, "rf.data.PCOA.SM.photo.rds")
saveRDS(rf.data.PCOA.SM.chemo, "rf.data.PCOA.SM.chemo.rds")
saveRDS(rf.data.PCOA.SM.hetero, "rf.data.PCOA.SM.hetero.rds")

saveRDS(rf.data.PCOA.CT, "rf.data.PCOA.CT.rds")
saveRDS(rf.data.PCOA.CT.photo, "rf.data.PCOA.CT.photo.rds")
saveRDS(rf.data.PCOA.CT.chemo, "rf.data.PCOA.CT.chemo.rds")
saveRDS(rf.data.PCOA.CT.hetero, "rf.data.PCOA.CT.hetero.rds")

saveRDS(rf.data.PCOA.SG, "rf.data.PCOA.SG.rds")
saveRDS(rf.data.PCOA.SG.photo, "rf.data.PCOA.SG.photo.rds")
saveRDS(rf.data.PCOA.SG.hetero, "rf.data.PCOA.SG.hetero.rds")


response.PCOA.data <- rarefied.PCoA.UnweightedUniFrac.phyloseq.regionellipse.samples$data

response.PCOA.axis1 <- (response.PCOA.data %>% select(Axis.1))

response.PCOA.axis2 <- (response.PCOA.data %>% select(Axis.2))

predictor.PCOA <- (metadata) %>% select("Region","Latittude", "Carbonate..ppm.", "pH", "Total.alkalinity..ppm.", "Phosphate..ppm.", "H2S..ppm.", "Temp...C." , "EC..mS.")

predictor.PCOA.reg.noalk <- (metadata) %>% select("Region","Carbonate..ppm.", "pH", "Phosphate..ppm.", "H2S..ppm.", "Temp...C." , "EC..mS.")


predictor.PCOA.lat.noalk <- (metadata) %>% select("Latittude", "Carbonate..ppm.", "pH", "Phosphate..ppm.", "H2S..ppm.", "Temp...C." , "EC..mS.")


rf.data.PCOA <- data.frame(response.PCOA.axis1, predictor.PCOA)
rf.data.PCOA.axis2 <- data.frame(response.PCOA.axis2, predictor.PCOA)

rf.data.PCOA.noreg <- data.frame(response.PCOA.axis1, predictor.PCOA.noreg)
rf.data.PCOA.axis2.noreg <- data.frame(response.PCOA.axis2, predictor.PCOA.noreg)

set.seed(151)
rf.PCOA <- randomForest(Axis.1~., data = rf.data.PCOA, ntree = 10001,importance=TRUE, proximity = TRUE)
print(rf.PCOA)
names(rf.PCOA)

set.seed(151)
rf.PCOA.NT.all <- randomForest(Axis.1~., data = rf.data.PCOA.NT.all, ntree = 10001,importance=TRUE, proximity = TRUE)
print(rf.PCOA.NT.all)
names(rf.PCOA.NT.all)


set.seed(151)
rf.PCOA.NT.reg <- randomForest(Axis.1~., data = rf.data.PCOA.NT.reg, ntree = 10001,importance=TRUE, proximity = TRUE)
print(rf.PCOA.NT.reg)
names(rf.PCOA.NT.reg)


set.seed(151)
rf.PCOA.NT.lat <- randomForest(Axis.1~., data = rf.data.PCOA.NT.lat, ntree = 10001,importance=TRUE, proximity = TRUE)
print(rf.PCOA.NT.lat)
names(rf.PCOA.NT.lat)


set.seed(151)
rf.PCOA.NT <- randomForest(Axis.1~., data = rf.data.PCOA.NT, ntree = 10001,importance=TRUE, proximity = TRUE)
print(rf.PCOA.NT)
names(rf.PCOA.NT)

set.seed(151)
rf.PCOA.NM <- randomForest(Axis.1~., data = rf.data.PCOA.NM, ntree = 10001,importance=TRUE, proximity = TRUE)
print(rf.PCOA.NM)

set.seed(151)
rf.PCOA.ST <- randomForest(Axis.1~., data = rf.data.PCOA.ST, ntree = 10001,importance=TRUE, proximity = TRUE)
print(rf.PCOA.ST)

set.seed(151)
rf.PCOA.SM <- randomForest(Axis.1~., data = rf.data.PCOA.SM, ntree = 10001,importance=TRUE, proximity = TRUE)
print(rf.PCOA.SM)

set.seed(151)
rf.PCOA.CT <- randomForest(Axis.1~., data = rf.data.PCOA.CT, ntree = 10001,importance=TRUE, proximity = TRUE)
print(rf.PCOA.CT)

set.seed(151)
rf.PCOA.SG <- randomForest(Axis.1~., data = rf.data.PCOA.SG, ntree = 10001,importance=TRUE, proximity = TRUE)
print(rf.PCOA.SG)


imp.rf.PCOA.axis1.rfp <- as.data.frame(importance(rf.PCOA.rfp ))
imp.rf.PCOA.axis1.rfp <- imp.rf.PCOA.axis1.rfp %>% rownames_to_column(var="response")

colnames(imp.rf.PCOA.axis1.rfp) <- c("response","MSE","MSE.pval","IncNodePurity","IncNodePurity.pval")

variable.imp.PCOA.axis1.MSE.bubbleplot <- imp.rf.PCOA.axis1.rfp  %>% ggplot( aes(x = reorder(response,MSE), y = MSE)) + geom_point(aes(fill = MSE.pval), alpha = 0.5, shape = 21,size = 6) + coord_flip()  + scale_size_continuous(limits = c(0,0.15), range = c(0,0.15), breaks = c(0,0.05,0.10,0.15)) + labs( x= "Environemntal variables", y = "% Increase in mean squared error (%IncMSE)", fill = "Mean squared error p-value") + theme(legend.key=element_blank(), axis.title = element_text(colour = "black", face = "bold", size = 8) , axis.text.x = element_text(colour = "black", size = 12, face = "bold", angle = 0, hjust = 0.5),   axis.text.y = element_text(colour = "black", face = "bold", size = 12), legend.title = element_text(size = 10, face = "bold"), panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2), legend.position = "right", title = element_text(size = 8))  + ggtitle("Factors influencing microbial community structure (PCoA Axis.1)")  + scale_fill_viridis() + theme_bw()

ggsave("variable.imp.PCOA.axis1.MSE.bubbleplot.pdf", height = 7, width = 9)
ggsave("variable.imp.PCOA.axis1.MSE.bubbleplot.png", height = 8, width = 7.5)


## RUn rfp for regions on server ####

#rfp.PCOA.NT <- rfPermute(Axis.1 ~ ., rf.data.PCOA.NT, ntree = 10001, num.rep = 1000)

#rfp.PCOA.NM <- rfPermute(Axis.1 ~ ., rf.data.PCOA.NM, ntree = 10001, num.rep = 1000)

#rfp.PCOA.ST <- rfPermute(Axis.1 ~ ., rf.data.PCOA.ST, ntree = 10001, num.rep = 1000)

#rfp.PCOA.SM <- rfPermute(Axis.1 ~ ., rf.data.PCOA.SM, ntree = 10001, num.rep = 1000)

#rfp.PCOA.CT <- rfPermute(Axis.1 ~ ., rf.data.PCOA.CT, ntree = 10001, num.rep = 1000)

#rfp.PCOA.SG <- rfPermute(Axis.1 ~ ., rf.data.PCOA.SG, ntree = 10001, num.rep = 1000)

#rfp.PCOA.ST.photo <- rfPermute(Axis.1 ~ ., rf.data.PCOA.ST.photo, ntree = 10001, num.rep = 1000, num.cores = 12)

#aveRDS(rfp.PCOA.ST.photo,"rfp.PCOA.ST.photo.rds")


#rfp.PCOA.ST.chemo <- rfPermute(Axis.1 ~ ., rf.data.PCOA.ST.chemo, ntree = 10001, num.rep = 1000, num.cores = 12)

#saveRDS(rfp.PCOA.ST.chemo,"rfp.PCOA.ST.chemo.rds")


#rfp.PCOA.ST.hetero <- rfPermute(Axis.1 ~ ., rf.data.PCOA.ST.hetero, ntree = 10001, num.rep = 1000, num.cores = 12)

#saveRDS(rfp.PCOA.ST.hetero,"rfp.PCOA.ST.hetero.rds")


### Run model significance

#rfp.PCOA.ST.photo.sig <- rf.significance(x=rfp.PCOA.ST.photo$rf, xdata=rf.data.PCOA.ST.photo[,2:(ncol(rf.data.PCOA.ST.photo))] , nperm=1000 , ntree=10001 ,  q=0.99, p=0.05)

#saveRDS(rfp.PCOA.ST.photo.sig,"rfp.PCOA.ST.photo.sig.rds")


#rfp.PCOA.ST.chemo.sig <- rf.significance(x=rfp.PCOA.ST.chemo$rf, xdata=rf.data.PCOA.ST.chemo[,2:(ncol(rf.data.PCOA.ST.chemo))] , nperm=1000 , ntree=10001 ,  q=0.99, p=0.05)

#saveRDS(rfp.PCOA.ST.chemo.sig,"rfp.PCOA.ST.chemo.sig.rds")


#rfp.PCOA.ST.hetero.sig <- rf.significance(x=rfp.PCOA.ST.hetero$rf, xdata=rf.data.PCOA.ST.hetero[,2:(ncol(rf.data.PCOA.ST.hetero))] , nperm=1000 , ntree=10001 ,  q=0.99, p=0.05)

#saveRDS(rfp.PCOA.ST.hetero.sig,"rfp.PCOA.ST.hetero.sig.rds")


##############################################################

### Import rfp for regions and plot graphs

## IMport rfp  for groups and plot graphs #### 

## ALWAYS DETACH RANGER BEFORE RUNNING RFP

### NT ####

rfp.PCOA.NT <- readRDS("rfp.PCOA.NT.rds")
rfp.PCOA.NT.sig <- readRDS("rfp.PCOA.NT.sig.rds")

imp.rf.PCOA.NT.axis1.rfp <- as.data.frame(importance(rfp.PCOA.NT ))
imp.rf.PCOA.NT.axis1.rfp <- imp.rf.PCOA.NT.axis1.rfp %>% rownames_to_column(var="response")

colnames(imp.rf.PCOA.NT.axis1.rfp) <- c("response","MSE","MSE.pval","IncNodePurity","IncNodePurity.pval")

imp.rf.PCOA.NT.axis1.rfp$Region <- 'All taxa'


rfp.PCOA.NT.photo <- readRDS("rfp.PCOA.NT.photo.rds")
rfp.PCOA.NT.photo.sig <- readRDS("rfp.PCOA.NT.photo.sig.rds")

imp.rf.PCOA.NT.photo.axis1.rfp <- as.data.frame(importance(rfp.PCOA.NT.photo ))
imp.rf.PCOA.NT.photo.axis1.rfp <- imp.rf.PCOA.NT.photo.axis1.rfp %>% rownames_to_column(var="response")

colnames(imp.rf.PCOA.NT.photo.axis1.rfp) <- c("response","MSE","MSE.pval","IncNodePurity","IncNodePurity.pval")

imp.rf.PCOA.NT.photo.axis1.rfp$Region <- 'Photosynthetic'


rfp.PCOA.NT.chemo <- readRDS("rfp.PCOA.NT.chemo.rds")
rfp.PCOA.NT.chemo.sig <- readRDS("rfp.PCOA.NT.chemo.sig.rds")

imp.rf.PCOA.NT.chemo.axis1.rfp <- as.data.frame(importance(rfp.PCOA.NT.chemo ))
imp.rf.PCOA.NT.chemo.axis1.rfp <- imp.rf.PCOA.NT.chemo.axis1.rfp %>% rownames_to_column(var="response")

colnames(imp.rf.PCOA.NT.chemo.axis1.rfp) <- c("response","MSE","MSE.pval","IncNodePurity","IncNodePurity.pval")

imp.rf.PCOA.NT.chemo.axis1.rfp$Region <- 'Chemoautotrophic'


rfp.PCOA.NT.hetero <- readRDS("rfp.PCOA.NT.hetero.rds")
rfp.PCOA.NT.hetero.sig <- readRDS("rfp.PCOA.NT.hetero.sig.rds")

imp.rf.PCOA.NT.hetero.axis1.rfp <- as.data.frame(importance(rfp.PCOA.NT.hetero ))
imp.rf.PCOA.NT.hetero.axis1.rfp <- imp.rf.PCOA.NT.hetero.axis1.rfp %>% rownames_to_column(var="response")

colnames(imp.rf.PCOA.NT.hetero.axis1.rfp) <- c("response","MSE","MSE.pval","IncNodePurity","IncNodePurity.pval")

imp.rf.PCOA.NT.hetero.axis1.rfp$Region <- 'Heterotrophic'

### combine all data 

NT.rf <- rbind(imp.rf.PCOA.NT.axis1.rfp, imp.rf.PCOA.NT.photo.axis1.rfp, imp.rf.PCOA.NT.chemo.axis1.rfp, imp.rf.PCOA.NT.hetero.axis1.rfp)

## convert continuous to discrete values
NT.rf$MSE.pval <- as.factor(NT.rf$MSE.pval)

sort(unique(NT.rf$MSE.pval)) 
#[1] 0.000999000999000999 0.001998001998002    0.002997002997003   
#[4] 0.003996003996004    0.004995004995005    0.00699300699300699 
#[7] 0.00799200799200799  0.00999000999000999  0.034965034965035 

# Define the number of colors you want
NT.cols <- 9
NTcolors <- colorRampPalette(brewer.pal(8, "Dark2"))(NT.cols)

NT.PCOA.axis1.MSE.bubbleplot <- NT.rf  %>% ggplot( aes(x = reorder(response,MSE), y = MSE, shape = factor(Region), color = MSE.pval)) + geom_point(alpha = 0.5, size = 3) + coord_flip()   + labs( x= "Environemntal variables", y = "% Increase in mean squared error (%IncMSE)", fill = "Mean squared error p-value") +  theme_bw() + theme( axis.title = element_text(colour = "black", face = "bold", size = 10) , axis.text.x = element_text(colour = "black", size = 10, face = "bold", angle = 0, hjust = 0.5),   axis.text.y = element_text(colour = "black", face = "bold", size = 10), legend.title = element_text(size = 9, face = "bold"), panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2), legend.position = "right",plot.title = element_text(size=10))  + ggtitle("Factors influencing microbial community structure at NT (PCoA Axis.1)")  + scale_shape_manual(values=c(15,16,17,18)) +  scale_color_manual(values = NTcolors)

ggsave("NT.PCOA.axis1.MSE.bubbleplot.pdf", height = 5.5, width = 7)
ggsave("NT.PCOA.axis1.MSE.bubbleplot.png", height = 5.5, width = 7)

#### multiple 4 plots for each group with MSE.pval as factor and facet #######

NT.PCOA.axis1.MSE.facet.bubbleplot <- NT.rf %>% 
  mutate(response = tidytext::reorder_within(response, MSE, within = Region))   %>% ggplot( aes(x = MSE, y = response, color = MSE.pval)) + geom_point(alpha = 0.5, size = 3)  + labs( x= "Environemntal variables", y = "% Increase in mean squared error (%IncMSE)", fill = "Mean squared error p-value") +  theme_bw() + theme( axis.title = element_text(colour = "black", face = "bold", size = 10) , axis.text.x = element_text(colour = "black", size = 10, face = "bold", angle = 0, hjust = 0.5),   axis.text.y = element_text(colour = "black", face = "bold", size = 10), legend.title = element_text(size = 9, face = "bold"), panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2), legend.position = "bottom",plot.title = element_text(size=10))  + ggtitle("Factors influencing microbial community structure at NT (PCoA Axis.1)")  + scale_shape_manual(values=c(15,16,17,18))+  scale_color_manual(values = NTcolors) + tidytext::scale_y_reordered() + facet_wrap(vars(Region),nrow = 1,  scales = "free_y")

ggsave("NT.PCOA.axis1.MSE.facet.bubbleplot.pdf", height = 5.5, width = 12)
ggsave("NT.PCOA.axis1.MSE.facet.bubbleplot.png", height = 5.5, width = 12)


#########################################

variable.imp.PCOA.NT.axis1.MSE.bubbleplot <- imp.rf.PCOA.NT.axis1.rfp  %>% ggplot( aes(x = reorder(response,MSE), y = MSE)) + geom_point(aes(fill = MSE.pval), alpha = 0.5, shape = 21,size = 6) + coord_flip()  + scale_size_continuous(limits = c(0,0.15), range = c(0,0.15), breaks = c(0,0.05,0.10,0.15)) + labs( x= "Environemntal variables", y = "% Increase in mean squared error (%IncMSE)", fill = "Mean squared error p-value") + theme(legend.key=element_blank(), axis.title = element_text(colour = "black", face = "bold", size = 8) , axis.text.x = element_text(colour = "black", size = 12, face = "bold", angle = 0, hjust = 0.5),   axis.text.y = element_text(colour = "black", face = "bold", size = 12), legend.title = element_text(size = 10, face = "bold"), panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2), legend.position = "right", title = element_text(size = 8))  + ggtitle("Factors influencing microbial community structure at NT (PCoA Axis.1)")  + scale_fill_viridis() + theme_bw()

ggsave("variable.imp.PCOA.NT.axis1.MSE.bubbleplot.pdf", height = 7, width = 9)
ggsave("variable.imp.PCOA.NT.axis1.MSE.bubbleplot.png", height = 8, width = 7.5)

variable.imp.PCOA.NT.axis1.IncNodePurity.bubbleplot <- imp.rf.PCOA.NT.axis1.rfp  %>% ggplot( aes(x = reorder(response,IncNodePurity), y = IncNodePurity)) + geom_point(aes(fill = IncNodePurity.pval), alpha = 0.5, shape = 21,size = 6) + coord_flip()  + scale_size_continuous(limits = c(0,0.15), range = c(0,0.15), breaks = c(0,0.05,0.10,0.15)) + labs( x= "Environemntal variables", y = "Increase in node purity (IncNodePurity)", fill = "Increase in node purity p-value") + theme(legend.key=element_blank(), axis.title = element_text(colour = "black", face = "bold", size = 8) , axis.text.x = element_text(colour = "black", size = 12, face = "bold", angle = 0, hjust = 0.5),   axis.text.y = element_text(colour = "black", face = "bold", size = 12), legend.title = element_text(size = 10, face = "bold"), panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2), legend.position = "right", title = element_text(size = 8))  + ggtitle("Factors influencing microbial community structure at NT (PCoA Axis.1)")  + scale_fill_viridis() + theme_bw()

ggsave("variable.imp.PCOA.NT.axis1.IncNodePurity.bubbleplot.pdf", height = 7, width = 9)
ggsave("variable.imp.PCOA.NT.axis1.IncNodePurity.bubbleplot.png", height = 8, width = 7.5)


### NM ####

rfp.PCOA.NM <- readRDS("rfp.PCOA.NM.rds")
rfp.PCOA.NM.sig <- readRDS("rfp.PCOA.NM.sig.rds")

imp.rf.PCOA.NM.axis1.rfp <- as.data.frame(importance(rfp.PCOA.NM ))
imp.rf.PCOA.NM.axis1.rfp <- imp.rf.PCOA.NM.axis1.rfp %>% rownames_to_column(var="response")

colnames(imp.rf.PCOA.NM.axis1.rfp) <- c("response","MSE","MSE.pval","IncNodePurity","IncNodePurity.pval")

imp.rf.PCOA.NM.axis1.rfp$Region <- 'All taxa'


rfp.PCOA.NM.photo <- readRDS("rfp.PCOA.NM.photo.rds")
rfp.PCOA.NM.photo.sig <- readRDS("rfp.PCOA.NM.photo.sig.rds")

imp.rf.PCOA.NM.photo.axis1.rfp <- as.data.frame(importance(rfp.PCOA.NM.photo ))
imp.rf.PCOA.NM.photo.axis1.rfp <- imp.rf.PCOA.NM.photo.axis1.rfp %>% rownames_to_column(var="response")

colnames(imp.rf.PCOA.NM.photo.axis1.rfp) <- c("response","MSE","MSE.pval","IncNodePurity","IncNodePurity.pval")

imp.rf.PCOA.NM.photo.axis1.rfp$Region <- 'Photosynthetic'


rfp.PCOA.NM.chemo <- readRDS("rfp.PCOA.NM.chemo.rds")
rfp.PCOA.NM.chemo.sig <- readRDS("rfp.PCOA.NM.chemo.sig.rds")

imp.rf.PCOA.NM.chemo.axis1.rfp <- as.data.frame(importance(rfp.PCOA.NM.chemo ))
imp.rf.PCOA.NM.chemo.axis1.rfp <- imp.rf.PCOA.NM.chemo.axis1.rfp %>% rownames_to_column(var="response")

colnames(imp.rf.PCOA.NM.chemo.axis1.rfp) <- c("response","MSE","MSE.pval","IncNodePurity","IncNodePurity.pval")

imp.rf.PCOA.NM.chemo.axis1.rfp$Region <- 'Chemoautotrophic'


rfp.PCOA.NM.hetero <- readRDS("rfp.PCOA.NM.hetero.rds")
rfp.PCOA.NM.hetero.sig <- readRDS("rfp.PCOA.NM.hetero.sig.rds")

imp.rf.PCOA.NM.hetero.axis1.rfp <- as.data.frame(importance(rfp.PCOA.NM.hetero ))
imp.rf.PCOA.NM.hetero.axis1.rfp <- imp.rf.PCOA.NM.hetero.axis1.rfp %>% rownames_to_column(var="response")

colnames(imp.rf.PCOA.NM.hetero.axis1.rfp) <- c("response","MSE","MSE.pval","IncNodePurity","IncNodePurity.pval")

imp.rf.PCOA.NM.hetero.axis1.rfp$Region <- 'Heterotrophic'

### combine all data 

NM.rf <- rbind(imp.rf.PCOA.NM.axis1.rfp, imp.rf.PCOA.NM.photo.axis1.rfp, imp.rf.PCOA.NM.chemo.axis1.rfp, imp.rf.PCOA.NM.hetero.axis1.rfp)

## convert continuous to discrete values
NM.rf$MSE.pval <- as.factor(NM.rf$MSE.pval)

sort(unique(NM.rf$MSE.pval)) 
#[1] 0.000999000999000999 0.001998001998002    0.003996003996004   
# [4] 0.004995004995005    0.010989010989011    0.012987012987013   
## [7] 0.015984015984016    0.016983016983017    0.212787212787213   
#[10] 0.231768231768232    0.257742257742258    0.271728271728272   
#[13] 1  

# Define the number of colors you want
NM.cols <- 13
NMcolors <- colorRampPalette(brewer.pal(8, "Dark2"))(NM.cols)

NM.PCOA.axis1.MSE.bubbleplot <- NM.rf  %>% ggplot( aes(x = reorder(response,MSE), y = MSE, shape = factor(Region), color = MSE.pval)) + geom_point(alpha = 0.5, size = 3) + coord_flip()   + labs( x= "Environemntal variables", y = "% Increase in mean squared error (%IncMSE)", fill = "Mean squared error p-value") +  theme_bw() + theme( axis.title = element_text(colour = "black", face = "bold", size = 10) , axis.text.x = element_text(colour = "black", size = 10, face = "bold", angle = 0, hjust = 0.5),   axis.text.y = element_text(colour = "black", face = "bold", size = 10), legend.title = element_text(size = 9, face = "bold"), panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2), legend.position = "right",plot.title = element_text(size=10))  + ggtitle("Factors influencing microbial community structure at NM (PCoA Axis.1)")  + scale_shape_manual(values=c(15,16,17,18)) +  scale_color_manual(values = NMcolors)

ggsave("NM.PCOA.axis1.MSE.bubbleplot.pdf", height = 5.5, width = 7)
ggsave("NM.PCOA.axis1.MSE.bubbleplot.png", height = 5.5, width = 7)

#### multiple 4 plots for each group with MSE.pval as factor and facet #######

NM.PCOA.axis1.MSE.facet.bubbleplot <- NM.rf %>% 
  mutate(response = tidytext::reorder_within(response, MSE, within = Region))   %>% ggplot( aes(x = MSE, y = response, color = MSE.pval)) + geom_point(alpha = 0.5, size = 3)  + labs( x= "Environemntal variables", y = "% Increase in mean squared error (%IncMSE)", fill = "Mean squared error p-value") +  theme_bw() + theme( axis.title = element_text(colour = "black", face = "bold", size = 10) , axis.text.x = element_text(colour = "black", size = 10, face = "bold", angle = 0, hjust = 0.5),   axis.text.y = element_text(colour = "black", face = "bold", size = 10), legend.title = element_text(size = 9, face = "bold"), panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2), legend.position = "bottom",plot.title = element_text(size=10))  + ggtitle("Factors influencing microbial community structure at NM (PCoA Axis.1)")  + scale_shape_manual(values=c(15,16,17,18)) +  scale_color_manual(values = NMcolors) + tidytext::scale_y_reordered() + facet_wrap(vars(Region),nrow = 1,  scales = "free_y")

ggsave("NM.PCOA.axis1.MSE.facet.bubbleplot.pdf", height = 5.5, width = 12)
ggsave("NM.PCOA.axis1.MSE.facet.bubbleplot.png", height = 5.5, width = 12)


#################################


variable.imp.PCOA.NM.axis1.MSE.bubbleplot <- imp.rf.PCOA.NM.axis1.rfp  %>% ggplot( aes(x = reorder(response,MSE), y = MSE)) + geom_point(aes(fill = MSE.pval), alpha = 0.5, shape = 21,size = 6) + coord_flip()  + scale_size_continuous(limits = c(0,0.15), range = c(0,0.15), breaks = c(0,0.05,0.10,0.15)) + labs( x= "Environemntal variables", y = "% Increase in mean squared error (%IncMSE)", fill = "Mean squared error p-value") + theme(legend.key=element_blank(), axis.title = element_text(colour = "black", face = "bold", size = 8) , axis.text.x = element_text(colour = "black", size = 12, face = "bold", angle = 0, hjust = 0.5),   axis.text.y = element_text(colour = "black", face = "bold", size = 12), legend.title = element_text(size = 10, face = "bold"), panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2), legend.position = "right", title = element_text(size = 8))  + ggtitle("Factors influencing microbial community structure at NM (PCoA Axis.1)")  + scale_fill_viridis() + theme_bw()

ggsave("variable.imp.PCOA.NM.axis1.MSE.bubbleplot.pdf", height = 7, width = 9)
ggsave("variable.imp.PCOA.NM.axis1.MSE.bubbleplot.png", height = 8, width = 7.5)

variable.imp.PCOA.NM.axis1.IncNodePurity.bubbleplot <- imp.rf.PCOA.NM.axis1.rfp  %>% ggplot( aes(x = reorder(response,IncNodePurity), y = IncNodePurity)) + geom_point(aes(fill = IncNodePurity.pval), alpha = 0.5, shape = 21,size = 6) + coord_flip()  + scale_size_continuous(limits = c(0,0.15), range = c(0,0.15), breaks = c(0,0.05,0.10,0.15)) + labs( x= "Environemntal variables", y = "Increase in node purity (IncNodePurity)", fill = "Increase in node purity p-value") + theme(legend.key=element_blank(), axis.title = element_text(colour = "black", face = "bold", size = 8) , axis.text.x = element_text(colour = "black", size = 12, face = "bold", angle = 0, hjust = 0.5),   axis.text.y = element_text(colour = "black", face = "bold", size = 12), legend.title = element_text(size = 10, face = "bold"), panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2), legend.position = "right", title = element_text(size = 8))  + ggtitle("Factors influencing microbial community structure at NM (PCoA Axis.1)")  + scale_fill_viridis() + theme_bw()

ggsave("variable.imp.PCOA.NM.axis1.IncNodePurity.bubbleplot.pdf", height = 7, width = 9)
ggsave("variable.imp.PCOA.NM.axis1.IncNodePurity.bubbleplot.png", height = 8, width = 7.5)


### ST ####

rfp.PCOA.ST <- readRDS("rfp.PCOA.ST.rds")
rfp.PCOA.ST.sig <- readRDS("rfp.PCOA.ST.sig.rds")

imp.rf.PCOA.ST.axis1.rfp <- as.data.frame(importance(rfp.PCOA.ST ))
imp.rf.PCOA.ST.axis1.rfp <- imp.rf.PCOA.ST.axis1.rfp %>% rownames_to_column(var="response")

colnames(imp.rf.PCOA.ST.axis1.rfp) <- c("response","MSE","MSE.pval","IncNodePurity","IncNodePurity.pval")

imp.rf.PCOA.ST.axis1.rfp$Region <- 'All taxa'


rfp.PCOA.ST.photo <- readRDS("rfp.PCOA.ST.photo.rds")
rfp.PCOA.ST.photo.sig <- readRDS("rfp.PCOA.ST.photo.sig.rds")

imp.rf.PCOA.ST.photo.axis1.rfp <- as.data.frame(importance(rfp.PCOA.ST.photo ))
imp.rf.PCOA.ST.photo.axis1.rfp <- imp.rf.PCOA.ST.photo.axis1.rfp %>% rownames_to_column(var="response")

colnames(imp.rf.PCOA.ST.photo.axis1.rfp) <- c("response","MSE","MSE.pval","IncNodePurity","IncNodePurity.pval")

imp.rf.PCOA.ST.photo.axis1.rfp$Region <- 'Photosynthetic'


rfp.PCOA.ST.chemo <- readRDS("rfp.PCOA.ST.chemo.rds")
rfp.PCOA.ST.chemo.sig <- readRDS("rfp.PCOA.ST.chemo.sig.rds")

imp.rf.PCOA.ST.chemo.axis1.rfp <- as.data.frame(importance(rfp.PCOA.ST.chemo ))
imp.rf.PCOA.ST.chemo.axis1.rfp <- imp.rf.PCOA.ST.chemo.axis1.rfp %>% rownames_to_column(var="response")

colnames(imp.rf.PCOA.ST.chemo.axis1.rfp) <- c("response","MSE","MSE.pval","IncNodePurity","IncNodePurity.pval")

imp.rf.PCOA.ST.chemo.axis1.rfp$Region <- 'Chemoautotrophic'


rfp.PCOA.ST.hetero <- readRDS("rfp.PCOA.ST.hetero.rds")
rfp.PCOA.ST.hetero.sig <- readRDS("rfp.PCOA.ST.hetero.sig.rds")

imp.rf.PCOA.ST.hetero.axis1.rfp <- as.data.frame(importance(rfp.PCOA.ST.hetero ))
imp.rf.PCOA.ST.hetero.axis1.rfp <- imp.rf.PCOA.ST.hetero.axis1.rfp %>% rownames_to_column(var="response")

colnames(imp.rf.PCOA.ST.hetero.axis1.rfp) <- c("response","MSE","MSE.pval","IncNodePurity","IncNodePurity.pval")

imp.rf.PCOA.ST.hetero.axis1.rfp$Region <- 'Heterotrophic'

### combine all data 

ST.rf <- rbind(imp.rf.PCOA.ST.axis1.rfp, imp.rf.PCOA.ST.photo.axis1.rfp, imp.rf.PCOA.ST.chemo.axis1.rfp, imp.rf.PCOA.ST.hetero.axis1.rfp)

## convert continuous to discrete values
ST.rf$MSE.pval <- as.factor(ST.rf$MSE.pval)

sort(unique(ST.rf$MSE.pval)) 
#[1] 0.000999000999000999 0.001998001998002    0.002997002997003   
#[4] 0.003996003996004    0.00699300699300699  0.00799200799200799 
#[7] 0.051948051948052    0.0539460539460539   0.0889110889110889  
#[10] 0.0909090909090909   1   

# Define the number of colors you want
ST.cols <- 11
STcolors <- colorRampPalette(brewer.pal(8, "Dark2"))(ST.cols)

ST.PCOA.axis1.MSE.bubbleplot <- ST.rf  %>% ggplot( aes(x = reorder(response,MSE), y = MSE, shape = factor(Region), color = MSE.pval)) + geom_point(alpha = 0.5, size = 3) + coord_flip()   + labs( x= "Environemntal variables", y = "% Increase in mean squared error (%IncMSE)", fill = "Mean squared error p-value") +  theme_bw() + theme( axis.title = element_text(colour = "black", face = "bold", size = 10) , axis.text.x = element_text(colour = "black", size = 10, face = "bold", angle = 0, hjust = 0.5),   axis.text.y = element_text(colour = "black", face = "bold", size = 10), legend.title = element_text(size = 9, face = "bold"), panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2), legend.position = "right",plot.title = element_text(size=10))  + ggtitle("Factors influencing microbial community structure at ST (PCoA Axis.1)")  + scale_shape_manual(values=c(15,16,17,18)) +  scale_color_manual(values = STcolors)

ggsave("ST.PCOA.axis1.MSE.bubbleplot.pdf", height = 5.5, width = 7)
ggsave("ST.PCOA.axis1.MSE.bubbleplot.png", height = 5.5, width = 7)

#### multiple 4 plots for each group with MSE.pval as factor and facet #######

ST.PCOA.axis1.MSE.facet.bubbleplot <- ST.rf %>% 
  mutate(response = tidytext::reorder_within(response, MSE, within = Region))   %>% ggplot( aes(x = MSE, y = response, color = MSE.pval)) + geom_point(alpha = 0.5, size = 3)  + labs( x= "Environemntal variables", y = "% Increase in mean squared error (%IncMSE)", fill = "Mean squared error p-value") +  theme_bw() + theme( axis.title = element_text(colour = "black", face = "bold", size = 10) , axis.text.x = element_text(colour = "black", size = 10, face = "bold", angle = 0, hjust = 0.5),   axis.text.y = element_text(colour = "black", face = "bold", size = 10), legend.title = element_text(size = 9, face = "bold"), panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2), legend.position = "bottom",plot.title = element_text(size=10))  + ggtitle("Factors influencing microbial community structure at ST (PCoA Axis.1)")  + scale_shape_manual(values=c(15,16,17,18)) +  scale_color_manual(values = STcolors) + tidytext::scale_y_reordered() + facet_wrap(vars(Region),nrow = 1,  scales = "free_y")

ggsave("ST.PCOA.axis1.MSE.facet.bubbleplot.pdf", height = 5.5, width = 12)
ggsave("ST.PCOA.axis1.MSE.facet.bubbleplot.png", height = 5.5, width = 12)


###########
variable.imp.PCOA.ST.axis1.MSE.bubbleplot <- imp.rf.PCOA.ST.axis1.rfp  %>% ggplot( aes(x = reorder(response,MSE), y = MSE)) + geom_point(aes(fill = MSE.pval), alpha = 0.5, shape = 21,size = 6) + coord_flip()  + scale_size_continuous(limits = c(0,0.15), range = c(0,0.15), breaks = c(0,0.05,0.10,0.15)) + labs( x= "Environemntal variables", y = "% Increase in mean squared error (%IncMSE)", fill = "Mean squared error p-value") + theme(legend.key=element_blank(), axis.title = element_text(colour = "black", face = "bold", size = 8) , axis.text.x = element_text(colour = "black", size = 12, face = "bold", angle = 0, hjust = 0.5),   axis.text.y = element_text(colour = "black", face = "bold", size = 12), legend.title = element_text(size = 10, face = "bold"), panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2), legend.position = "right", title = element_text(size = 8))  + ggtitle("Factors influencing microbial community structure at ST (PCoA Axis.1)")  + scale_fill_viridis() + theme_bw()

ggsave("variable.imp.PCOA.ST.axis1.MSE.bubbleplot.pdf", height = 7, width = 9)
ggsave("variable.imp.PCOA.ST.axis1.MSE.bubbleplot.png", height = 8, width = 7.5)

variable.imp.PCOA.ST.axis1.IncNodePurity.bubbleplot <- imp.rf.PCOA.ST.axis1.rfp  %>% ggplot( aes(x = reorder(response,IncNodePurity), y = IncNodePurity)) + geom_point(aes(fill = IncNodePurity.pval), alpha = 0.5, shape = 21,size = 6) + coord_flip()  + scale_size_continuous(limits = c(0,0.15), range = c(0,0.15), breaks = c(0,0.05,0.10,0.15)) + labs( x= "Environemntal variables", y = "Increase in node purity (IncNodePurity)", fill = "Increase in node purity p-value") + theme(legend.key=element_blank(), axis.title = element_text(colour = "black", face = "bold", size = 8) , axis.text.x = element_text(colour = "black", size = 12, face = "bold", angle = 0, hjust = 0.5),   axis.text.y = element_text(colour = "black", face = "bold", size = 12), legend.title = element_text(size = 10, face = "bold"), panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2), legend.position = "right", title = element_text(size = 8))  + ggtitle("Factors influencing microbial community structure at ST (PCoA Axis.1)")  + scale_fill_viridis() + theme_bw()

ggsave("variable.imp.PCOA.ST.axis1.IncNodePurity.bubbleplot.pdf", height = 7, width = 9)
ggsave("variable.imp.PCOA.ST.axis1.IncNodePurity.bubbleplot.png", height = 8, width = 7.5)


### SM ####

rfp.PCOA.SM <- readRDS("rfp.PCOA.SM.rds")
rfp.PCOA.SM.sig <- readRDS("rfp.PCOA.SM.sig.rds")

imp.rf.PCOA.SM.axis1.rfp <- as.data.frame(importance(rfp.PCOA.SM ))
imp.rf.PCOA.SM.axis1.rfp <- imp.rf.PCOA.SM.axis1.rfp %>% rownames_to_column(var="response")

colnames(imp.rf.PCOA.SM.axis1.rfp) <- c("response","MSE","MSE.pval","IncNodePurity","IncNodePurity.pval")

imp.rf.PCOA.SM.axis1.rfp$Region <- 'All taxa'


rfp.PCOA.SM.photo <- readRDS("rfp.PCOA.SM.photo.rds")
rfp.PCOA.SM.photo.sig <- readRDS("rfp.PCOA.SM.photo.sig.rds")

imp.rf.PCOA.SM.photo.axis1.rfp <- as.data.frame(importance(rfp.PCOA.SM.photo ))
imp.rf.PCOA.SM.photo.axis1.rfp <- imp.rf.PCOA.SM.photo.axis1.rfp %>% rownames_to_column(var="response")

colnames(imp.rf.PCOA.SM.photo.axis1.rfp) <- c("response","MSE","MSE.pval","IncNodePurity","IncNodePurity.pval")

imp.rf.PCOA.SM.photo.axis1.rfp$Region <- 'Photosynthetic'


rfp.PCOA.SM.chemo <- readRDS("rfp.PCOA.SM.chemo.rds")
rfp.PCOA.SM.chemo.sig <- readRDS("rfp.PCOA.SM.chemo.sig.rds")

imp.rf.PCOA.SM.chemo.axis1.rfp <- as.data.frame(importance(rfp.PCOA.SM.chemo ))
imp.rf.PCOA.SM.chemo.axis1.rfp <- imp.rf.PCOA.SM.chemo.axis1.rfp %>% rownames_to_column(var="response")

colnames(imp.rf.PCOA.SM.chemo.axis1.rfp) <- c("response","MSE","MSE.pval","IncNodePurity","IncNodePurity.pval")

imp.rf.PCOA.SM.chemo.axis1.rfp$Region <- 'Chemoautotrophic'


rfp.PCOA.SM.hetero <- readRDS("rfp.PCOA.SM.hetero.rds")
rfp.PCOA.SM.hetero.sig <- readRDS("rfp.PCOA.SM.hetero.sig.rds")

imp.rf.PCOA.SM.hetero.axis1.rfp <- as.data.frame(importance(rfp.PCOA.SM.hetero ))
imp.rf.PCOA.SM.hetero.axis1.rfp <- imp.rf.PCOA.SM.hetero.axis1.rfp %>% rownames_to_column(var="response")

colnames(imp.rf.PCOA.SM.hetero.axis1.rfp) <- c("response","MSE","MSE.pval","IncNodePurity","IncNodePurity.pval")

imp.rf.PCOA.SM.hetero.axis1.rfp$Region <- 'Heterotrophic'

### combine all data 

SM.rf <- rbind(imp.rf.PCOA.SM.axis1.rfp, imp.rf.PCOA.SM.photo.axis1.rfp, imp.rf.PCOA.SM.chemo.axis1.rfp, imp.rf.PCOA.SM.hetero.axis1.rfp)

## convert continuous to discrete values
SM.rf$MSE.pval <- as.factor(SM.rf$MSE.pval)

sort(unique(SM.rf$MSE.pval)) 
#[1] 0.000999000999000999 0.00699300699300699  0.00799200799200799 
# [4] 0.00899100899100899  0.00999000999000999  0.015984015984016   
# [7] 0.021978021978022    0.027972027972028    0.030969030969031   
#[10] 0.040959040959041    0.0529470529470529   0.0539460539460539  
#[13] 0.0639360639360639   1    

# Define the number of colors you want
SM.cols <- 14
SMcolors <- colorRampPalette(brewer.pal(8, "Dark2"))(SM.cols)

SM.PCOA.axis1.MSE.bubbleplot <- SM.rf  %>% ggplot( aes(x = reorder(response,MSE), y = MSE, shape = factor(Region), color = MSE.pval)) + geom_point(alpha = 0.5, size = 3) + coord_flip()   + labs( x= "Environemntal variables", y = "% Increase in mean squared error (%IncMSE)", fill = "Mean squared error p-value") +  theme_bw() + theme( axis.title = element_text(colour = "black", face = "bold", size = 10) , axis.text.x = element_text(colour = "black", size = 10, face = "bold", angle = 0, hjust = 0.5),   axis.text.y = element_text(colour = "black", face = "bold", size = 10), legend.title = element_text(size = 9, face = "bold"), panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2), legend.position = "right",plot.title = element_text(size=10))  + ggtitle("Factors influencing microbial community structure at SM (PCoA Axis.1)")  + scale_shape_manual(values=c(15,16,17,18)) +  scale_color_manual(values = SMcolors)

ggsave("SM.PCOA.axis1.MSE.bubbleplot.pdf", height = 5.5, width = 7)
ggsave("SM.PCOA.axis1.MSE.bubbleplot.png", height = 5.5, width = 7)

#### multiple 4 plots for each group with MSE.pval as factor and facet #######

SM.PCOA.axis1.MSE.facet.bubbleplot <- SM.rf %>% 
  mutate(response = tidytext::reorder_within(response, MSE, within = Region))   %>% ggplot( aes(x = MSE, y = response, color = MSE.pval)) + geom_point(alpha = 0.5, size = 3)  + labs( x= "Environemntal variables", y = "% Increase in mean squared error (%IncMSE)", fill = "Mean squared error p-value") +  theme_bw() + theme( axis.title = element_text(colour = "black", face = "bold", size = 10) , axis.text.x = element_text(colour = "black", size = 10, face = "bold", angle = 0, hjust = 0.5),   axis.text.y = element_text(colour = "black", face = "bold", size = 10), legend.title = element_text(size = 9, face = "bold"), panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2), legend.position = "bottom",plot.title = element_text(size=10))  + ggtitle("Factors influencing microbial community structure at SM (PCoA Axis.1)")  + scale_shape_manual(values=c(15,16,17,18)) +  scale_color_manual(values = SMcolors) + tidytext::scale_y_reordered() + facet_wrap(vars(Region),nrow = 1,  scales = "free_y")

ggsave("SM.PCOA.axis1.MSE.facet.bubbleplot.pdf", height = 5.5, width = 12)
ggsave("SM.PCOA.axis1.MSE.facet.bubbleplot.png", height = 5.5, width = 12)


variable.imp.PCOA.SM.axis1.MSE.bubbleplot <- imp.rf.PCOA.SM.axis1.rfp  %>% ggplot( aes(x = reorder(response,MSE), y = MSE)) + geom_point(aes(fill = MSE.pval), alpha = 0.5, shape = 21,size = 6) + coord_flip()  + scale_size_continuous(limits = c(0,0.15), range = c(0,0.15), breaks = c(0,0.05,0.10,0.15)) + labs( x= "Environemntal variables", y = "% Increase in mean squared error (%IncMSE)", fill = "Mean squared error p-value") + theme(legend.key=element_blank(), axis.title = element_text(colour = "black", face = "bold", size = 8) , axis.text.x = element_text(colour = "black", size = 12, face = "bold", angle = 0, hjust = 0.5),   axis.text.y = element_text(colour = "black", face = "bold", size = 12), legend.title = element_text(size = 10, face = "bold"), panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2), legend.position = "right", title = element_text(size = 8))  + ggtitle("Factors influencing microbial community structure at SM (PCoA Axis.1)")  + scale_fill_viridis() + theme_bw()

ggsave("variable.imp.PCOA.SM.axis1.MSE.bubbleplot.pdf", height = 7, width = 9)
ggsave("variable.imp.PCOA.SM.axis1.MSE.bubbleplot.png", height = 8, width = 7.5)

variable.imp.PCOA.SM.axis1.IncNodePurity.bubbleplot <- imp.rf.PCOA.SM.axis1.rfp  %>% ggplot( aes(x = reorder(response,IncNodePurity), y = IncNodePurity)) + geom_point(aes(fill = IncNodePurity.pval), alpha = 0.5, shape = 21,size = 6) + coord_flip()  + scale_size_continuous(limits = c(0,0.15), range = c(0,0.15), breaks = c(0,0.05,0.10,0.15)) + labs( x= "Environemntal variables", y = "Increase in node purity (IncNodePurity)", fill = "Increase in node purity p-value") + theme(legend.key=element_blank(), axis.title = element_text(colour = "black", face = "bold", size = 8) , axis.text.x = element_text(colour = "black", size = 12, face = "bold", angle = 0, hjust = 0.5),   axis.text.y = element_text(colour = "black", face = "bold", size = 12), legend.title = element_text(size = 10, face = "bold"), panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2), legend.position = "right", title = element_text(size = 8))  + ggtitle("Factors influencing microbial community structure at SM (PCoA Axis.1)")  + scale_fill_viridis() + theme_bw()

ggsave("variable.imp.PCOA.SM.axis1.IncNodePurity.bubbleplot.pdf", height = 7, width = 9)
ggsave("variable.imp.PCOA.SM.axis1.IncNodePurity.bubbleplot.png", height = 8, width = 7.5)


### CT ####

rfp.PCOA.CT <- readRDS("rfp.PCOA.CT.rds")
rfp.PCOA.CT.sig <- readRDS("rfp.PCOA.CT.sig.rds")

imp.rf.PCOA.CT.axis1.rfp <- as.data.frame(importance(rfp.PCOA.CT ))
imp.rf.PCOA.CT.axis1.rfp <- imp.rf.PCOA.CT.axis1.rfp %>% rownames_to_column(var="response")

colnames(imp.rf.PCOA.CT.axis1.rfp) <- c("response","MSE","MSE.pval","IncNodePurity","IncNodePurity.pval")

imp.rf.PCOA.CT.axis1.rfp$Region <- 'All taxa'


rfp.PCOA.CT.photo <- readRDS("rfp.PCOA.CT.photo.rds")
rfp.PCOA.CT.photo.sig <- readRDS("rfp.PCOA.CT.photo.sig.rds")

imp.rf.PCOA.CT.photo.axis1.rfp <- as.data.frame(importance(rfp.PCOA.CT.photo ))
imp.rf.PCOA.CT.photo.axis1.rfp <- imp.rf.PCOA.CT.photo.axis1.rfp %>% rownames_to_column(var="response")

colnames(imp.rf.PCOA.CT.photo.axis1.rfp) <- c("response","MSE","MSE.pval","IncNodePurity","IncNodePurity.pval")

imp.rf.PCOA.CT.photo.axis1.rfp$Region <- 'Photosynthetic'


rfp.PCOA.CT.chemo <- readRDS("rfp.PCOA.CT.chemo.rds")
rfp.PCOA.CT.chemo.sig <- readRDS("rfp.PCOA.CT.chemo.sig.rds")

imp.rf.PCOA.CT.chemo.axis1.rfp <- as.data.frame(importance(rfp.PCOA.CT.chemo ))
imp.rf.PCOA.CT.chemo.axis1.rfp <- imp.rf.PCOA.CT.chemo.axis1.rfp %>% rownames_to_column(var="response")

colnames(imp.rf.PCOA.CT.chemo.axis1.rfp) <- c("response","MSE","MSE.pval","IncNodePurity","IncNodePurity.pval")

imp.rf.PCOA.CT.chemo.axis1.rfp$Region <- 'Chemoautotrophic'


rfp.PCOA.CT.hetero <- readRDS("rfp.PCOA.CT.hetero.rds")
rfp.PCOA.CT.hetero.sig <- readRDS("rfp.PCOA.CT.hetero.sig.rds")

imp.rf.PCOA.CT.hetero.axis1.rfp <- as.data.frame(importance(rfp.PCOA.CT.hetero ))
imp.rf.PCOA.CT.hetero.axis1.rfp <- imp.rf.PCOA.CT.hetero.axis1.rfp %>% rownames_to_column(var="response")

colnames(imp.rf.PCOA.CT.hetero.axis1.rfp) <- c("response","MSE","MSE.pval","IncNodePurity","IncNodePurity.pval")

imp.rf.PCOA.CT.hetero.axis1.rfp$Region <- 'Heterotrophic'

### combine all data 

CT.rf <- rbind(imp.rf.PCOA.CT.axis1.rfp, imp.rf.PCOA.CT.photo.axis1.rfp, imp.rf.PCOA.CT.chemo.axis1.rfp, imp.rf.PCOA.CT.hetero.axis1.rfp)

## convert continuous to discrete values
CT.rf$MSE.pval <- as.factor(CT.rf$MSE.pval)

sort(unique(CT.rf$MSE.pval)) 
# [1] 0.000999000999000999 0.001998001998002    0.003996003996004   
#[4] 0.004995004995005    0.00699300699300699  0.00799200799200799 
#[7] 0.01998001998002     0.027972027972028    1  

# Define the number of colors you want
CT.cols <- 9
CTcolors <- colorRampPalette(brewer.pal(8, "Dark2"))(CT.cols)

CT.PCOA.axis1.MSE.bubbleplot <- CT.rf  %>% ggplot( aes(x = reorder(response,MSE), y = MSE, shape = factor(Region), color = MSE.pval)) + geom_point(alpha = 0.5, size = 3) + coord_flip()   + labs( x= "Environemntal variables", y = "% Increase in mean squared error (%IncMSE)", fill = "Mean squared error p-value") +  theme_bw() + theme( axis.title = element_text(colour = "black", face = "bold", size = 10) , axis.text.x = element_text(colour = "black", size = 10, face = "bold", angle = 0, hjust = 0.5),   axis.text.y = element_text(colour = "black", face = "bold", size = 10), legend.title = element_text(size = 9, face = "bold"), panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2), legend.position = "right",plot.title = element_text(size=10))  + ggtitle("Factors influencing microbial community structure at CT (PCoA Axis.1)")  + scale_shape_manual(values=c(15,16,17,18)) +  scale_color_manual(values = CTcolors)

ggsave("CT.PCOA.axis1.MSE.bubbleplot.pdf", height = 5.5, width = 7)
ggsave("CT.PCOA.axis1.MSE.bubbleplot.png", height = 5.5, width = 7)

#### multiple 4 plots for each group with MSE.pval as factor and facet #######

CT.PCOA.axis1.MSE.facet.bubbleplot <- CT.rf %>% 
  mutate(response = tidytext::reorder_within(response, MSE, within = Region))   %>% ggplot( aes(x = MSE, y = response, color = MSE.pval)) + geom_point(alpha = 0.5, size = 3)  + labs( x= "Environemntal variables", y = "% Increase in mean squared error (%IncMSE)", fill = "Mean squared error p-value") +  theme_bw() + theme( axis.title = element_text(colour = "black", face = "bold", size = 10) , axis.text.x = element_text(colour = "black", size = 10, face = "bold", angle = 0, hjust = 0.5),   axis.text.y = element_text(colour = "black", face = "bold", size = 10), legend.title = element_text(size = 9, face = "bold"), panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2), legend.position = "bottom",plot.title = element_text(size=10))  + ggtitle("Factors influencing microbial community structure at CT (PCoA Axis.1)")  + scale_shape_manual(values=c(15,16,17,18)) +  scale_color_manual(values = CTcolors) + tidytext::scale_y_reordered() + facet_wrap(vars(Region),nrow = 1,  scales = "free_y")

ggsave("CT.PCOA.axis1.MSE.facet.bubbleplot.pdf", height = 5.5, width = 12)
ggsave("CT.PCOA.axis1.MSE.facet.bubbleplot.png", height = 5.5, width = 12)



variable.imp.PCOA.CT.axis1.MSE.bubbleplot <- imp.rf.PCOA.CT.axis1.rfp  %>% ggplot( aes(x = reorder(response,MSE), y = MSE)) + geom_point(aes(fill = MSE.pval), alpha = 0.5, shape = 21,size = 6) + coord_flip()  + scale_size_continuous(limits = c(0,0.15), range = c(0,0.15), breaks = c(0,0.05,0.10,0.15)) + labs( x= "Environemntal variables", y = "% Increase in mean squared error (%IncMSE)", fill = "Mean squared error p-value") + theme(legend.key=element_blank(), axis.title = element_text(colour = "black", face = "bold", size = 8) , axis.text.x = element_text(colour = "black", size = 12, face = "bold", angle = 0, hjust = 0.5),   axis.text.y = element_text(colour = "black", face = "bold", size = 12), legend.title = element_text(size = 10, face = "bold"), panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2), legend.position = "right", title = element_text(size = 8))  + ggtitle("Factors influencing microbial community structure at CT (PCoA Axis.1)")  + scale_fill_viridis() + theme_bw()

ggsave("variable.imp.PCOA.CT.axis1.MSE.bubbleplot.pdf", height = 7, width = 9)
ggsave("variable.imp.PCOA.CT.axis1.MSE.bubbleplot.png", height = 8, width = 7.5)

variable.imp.PCOA.CT.axis1.IncNodePurity.bubbleplot <- imp.rf.PCOA.CT.axis1.rfp  %>% ggplot( aes(x = reorder(response,IncNodePurity), y = IncNodePurity)) + geom_point(aes(fill = IncNodePurity.pval), alpha = 0.5, shape = 21,size = 6) + coord_flip()  + scale_size_continuous(limits = c(0,0.15), range = c(0,0.15), breaks = c(0,0.05,0.10,0.15)) + labs( x= "Environemntal variables", y = "Increase in node purity (IncNodePurity)", fill = "Increase in node purity p-value") + theme(legend.key=element_blank(), axis.title = element_text(colour = "black", face = "bold", size = 8) , axis.text.x = element_text(colour = "black", size = 12, face = "bold", angle = 0, hjust = 0.5),   axis.text.y = element_text(colour = "black", face = "bold", size = 12), legend.title = element_text(size = 10, face = "bold"), panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2), legend.position = "right", title = element_text(size = 8))  + ggtitle("Factors influencing microbial community structure at CT (PCoA Axis.1)")  + scale_fill_viridis() + theme_bw()

ggsave("variable.imp.PCOA.CT.axis1.IncNodePurity.bubbleplot.pdf", height = 7, width = 9)
ggsave("variable.imp.PCOA.CT.axis1.IncNodePurity.bubbleplot.png", height = 8, width = 7.5)


### SG ####

rfp.PCOA.SG <- readRDS("rfp.PCOA.SG.rds")
rfp.PCOA.SG.sig <- readRDS("rfp.PCOA.SG.sig.rds")

imp.rf.PCOA.SG.axis1.rfp <- as.data.frame(importance(rfp.PCOA.SG ))
imp.rf.PCOA.SG.axis1.rfp <- imp.rf.PCOA.SG.axis1.rfp %>% rownames_to_column(var="response")

colnames(imp.rf.PCOA.SG.axis1.rfp) <- c("response","MSE","MSE.pval","IncNodePurity","IncNodePurity.pval")

imp.rf.PCOA.SG.axis1.rfp$Region <- 'All taxa'


rfp.PCOA.SG.photo <- readRDS("rfp.PCOA.SG.photo.rds")
rfp.PCOA.SG.photo.sig <- readRDS("rfp.PCOA.SG.photo.sig.rds")

imp.rf.PCOA.SG.photo.axis1.rfp <- as.data.frame(importance(rfp.PCOA.SG.photo ))
imp.rf.PCOA.SG.photo.axis1.rfp <- imp.rf.PCOA.SG.photo.axis1.rfp %>% rownames_to_column(var="response")

colnames(imp.rf.PCOA.SG.photo.axis1.rfp) <- c("response","MSE","MSE.pval","IncNodePurity","IncNodePurity.pval")

imp.rf.PCOA.SG.photo.axis1.rfp$Region <- 'Photosynthetic'



rfp.PCOA.SG.hetero <- readRDS("rfp.PCOA.SG.hetero.rds")
rfp.PCOA.SG.hetero.sig <- readRDS("rfp.PCOA.SG.hetero.sig.rds")

imp.rf.PCOA.SG.hetero.axis1.rfp <- as.data.frame(importance(rfp.PCOA.SG.hetero ))
imp.rf.PCOA.SG.hetero.axis1.rfp <- imp.rf.PCOA.SG.hetero.axis1.rfp %>% rownames_to_column(var="response")

colnames(imp.rf.PCOA.SG.hetero.axis1.rfp) <- c("response","MSE","MSE.pval","IncNodePurity","IncNodePurity.pval")

imp.rf.PCOA.SG.hetero.axis1.rfp$Region <- 'Heterotrophic'

### combine all data 

SG.rf <- rbind(imp.rf.PCOA.SG.axis1.rfp, imp.rf.PCOA.SG.photo.axis1.rfp, imp.rf.PCOA.SG.hetero.axis1.rfp)

## convert continuous to discrete values
SG.rf$MSE.pval <- as.factor(SG.rf$MSE.pval)

sort(unique(SG.rf$MSE.pval)) 
#[1] 0.000999000999000999 0.001998001998002    0.172827172827173   
#[4] 0.244755244755245    0.253746253746254    1   

# Define the number of colors you want
SG.cols <- 6
SGcolors <- colorRampPalette(brewer.pal(6, "Dark2"))(SG.cols)

SG.PCOA.axis1.MSE.bubbleplot <- SG.rf  %>% ggplot( aes(x = reorder(response,MSE), y = MSE, shape = factor(Region), color = MSE.pval)) + geom_point(alpha = 0.5, size = 3) + coord_flip()   + labs( x= "Environemntal variables", y = "% Increase in mean squared error (%IncMSE)", fill = "Mean squared error p-value") +  theme_bw() + theme( axis.title = element_text(colour = "black", face = "bold", size = 10) , axis.text.x = element_text(colour = "black", size = 10, face = "bold", angle = 0, hjust = 0.5),   axis.text.y = element_text(colour = "black", face = "bold", size = 10), legend.title = element_text(size = 9, face = "bold"), panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2), legend.position = "right",plot.title = element_text(size=10))  + ggtitle("Factors influencing microbial community structure at SG (PCoA Axis.1)")  + scale_shape_manual(values=c(15,16,17,18)) +  scale_color_manual(values = SGcolors)

ggsave("SG.PCOA.axis1.MSE.bubbleplot.pdf", height = 5.5, width = 7)
ggsave("SG.PCOA.axis1.MSE.bubbleplot.png", height = 5.5, width = 7)

#### multiple 4 plots for each group with MSE.pval as factor and facet #######

SG.PCOA.axis1.MSE.facet.bubbleplot <- SG.rf %>% 
  mutate(response = tidytext::reorder_within(response, MSE, within = Region))   %>% ggplot( aes(x = MSE, y = response, color = MSE.pval)) + geom_point(alpha = 0.5, size = 3)  + labs( x= "Environemntal variables", y = "% Increase in mean squared error (%IncMSE)", fill = "Mean squared error p-value") +  theme_bw() + theme( axis.title = element_text(colour = "black", face = "bold", size = 10) , axis.text.x = element_text(colour = "black", size = 10, face = "bold", angle = 0, hjust = 0.5),   axis.text.y = element_text(colour = "black", face = "bold", size = 10), legend.title = element_text(size = 9, face = "bold"), panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2), legend.position = "bottom",plot.title = element_text(size=10))  + ggtitle("Factors influencing microbial community structure at SG (PCoA Axis.1)")  + scale_shape_manual(values=c(15,16,17,18)) +  scale_color_manual(values = SGcolors) + tidytext::scale_y_reordered() + facet_wrap(vars(Region),nrow = 1,  scales = "free_y")

ggsave("SG.PCOA.axis1.MSE.facet.bubbleplot.pdf", height = 5.5, width = 12)
ggsave("SG.PCOA.axis1.MSE.facet.bubbleplot.png", height = 5.5, width = 12)



variable.imp.PCOA.SG.axis1.MSE.bubbleplot <- imp.rf.PCOA.SG.axis1.rfp  %>% ggplot( aes(x = reorder(response,MSE), y = MSE)) + geom_point(aes(fill = MSE.pval), alpha = 0.5, shape = 21,size = 6) + coord_flip()  + scale_size_continuous(limits = c(0,0.15), range = c(0,0.15), breaks = c(0,0.05,0.10,0.15)) + labs( x= "Environemntal variables", y = "% Increase in mean squared error (%IncMSE)", fill = "Mean squared error p-value") + theme(legend.key=element_blank(), axis.title = element_text(colour = "black", face = "bold", size = 8) , axis.text.x = element_text(colour = "black", size = 12, face = "bold", angle = 0, hjust = 0.5),   axis.text.y = element_text(colour = "black", face = "bold", size = 12), legend.title = element_text(size = 10, face = "bold"), panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2), legend.position = "right", title = element_text(size = 8))  + ggtitle("Factors influencing microbial community structure at SG (PCoA Axis.1)")  + scale_fill_viridis() + theme_bw()

ggsave("variable.imp.PCOA.SG.axis1.MSE.bubbleplot.pdf", height = 7, width = 9)
ggsave("variable.imp.PCOA.SG.axis1.MSE.bubbleplot.png", height = 8, width = 7.5)

variable.imp.PCOA.SG.axis1.IncNodePurity.bubbleplot <- imp.rf.PCOA.SG.axis1.rfp  %>% ggplot( aes(x = reorder(response,IncNodePurity), y = IncNodePurity)) + geom_point(aes(fill = IncNodePurity.pval), alpha = 0.5, shape = 21,size = 6) + coord_flip()  + scale_size_continuous(limits = c(0,0.15), range = c(0,0.15), breaks = c(0,0.05,0.10,0.15)) + labs( x= "Environemntal variables", y = "Increase in node purity (IncNodePurity)", fill = "Increase in node purity p-value") + theme(legend.key=element_blank(), axis.title = element_text(colour = "black", face = "bold", size = 8) , axis.text.x = element_text(colour = "black", size = 12, face = "bold", angle = 0, hjust = 0.5),   axis.text.y = element_text(colour = "black", face = "bold", size = 12), legend.title = element_text(size = 10, face = "bold"), panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2), legend.position = "right", title = element_text(size = 8))  + ggtitle("Factors influencing microbial community structure at SG (PCoA Axis.1)")  + scale_fill_viridis() + theme_bw()

ggsave("variable.imp.PCOA.SG.axis1.IncNodePurity.bubbleplot.pdf", height = 7, width = 9)
ggsave("variable.imp.PCOA.SG.axis1.IncNodePurity.bubbleplot.png", height = 8, width = 7.5)


### ALL ####

## wihtout lat,long,tot.alk.

rfp.PCOA.select <- readRDS("rfp.PCOA.select.rds")
rfp.PCOA.select.sig <- readRDS("rfp.PCOA.select.sig.rds")

imp.rf.PCOA.select.axis1.rfp <- as.data.frame(importance(rfp.PCOA.select ))
imp.rf.PCOA.select.axis1.rfp <- imp.rf.PCOA.select.axis1.rfp %>% rownames_to_column(var="response")

colnames(imp.rf.PCOA.select.axis1.rfp) <- c("response","MSE","MSE.pval","IncNodePurity","IncNodePurity.pval")

imp.rf.PCOA.select.axis1.rfp$Region <- 'All taxa'


rfp.PCOA.photo.select <- readRDS("rfp.PCOA.photo.select.rds")
rfp.PCOA.photo.select.sig <- readRDS("rfp.PCOA.photo.select.sig.rds")

imp.rf.PCOA.photo.select.axis1.rfp <- as.data.frame(importance(rfp.PCOA.photo.select ))
imp.rf.PCOA.photo.select.axis1.rfp <- imp.rf.PCOA.photo.select.axis1.rfp %>% rownames_to_column(var="response")

colnames(imp.rf.PCOA.photo.select.axis1.rfp) <- c("response","MSE","MSE.pval","IncNodePurity","IncNodePurity.pval")

imp.rf.PCOA.photo.select.axis1.rfp$Region <- 'Photosynthetic'


rfp.PCOA.chemo.select <- readRDS("rfp.PCOA.chemo.select.rds")
rfp.PCOA.chemo.select.sig <- readRDS("rfp.PCOA.chemo.select.sig.rds")

imp.rf.PCOA.chemo.select.axis1.rfp <- as.data.frame(importance(rfp.PCOA.chemo.select ))
imp.rf.PCOA.chemo.select.axis1.rfp <- imp.rf.PCOA.chemo.select.axis1.rfp %>% rownames_to_column(var="response")

colnames(imp.rf.PCOA.chemo.select.axis1.rfp) <- c("response","MSE","MSE.pval","IncNodePurity","IncNodePurity.pval")

imp.rf.PCOA.chemo.select.axis1.rfp$Region <- 'Chemoautotrophic'


rfp.PCOA.hetero.select <- readRDS("rfp.PCOA.hetero.select.rds")
rfp.PCOA.hetero.select.sig <- readRDS("rfp.PCOA.hetero.select.sig.rds")

imp.rf.PCOA.hetero.select.axis1.rfp <- as.data.frame(importance(rfp.PCOA.hetero.select ))
imp.rf.PCOA.hetero.select.axis1.rfp <- imp.rf.PCOA.hetero.select.axis1.rfp %>% rownames_to_column(var="response")

colnames(imp.rf.PCOA.hetero.select.axis1.rfp) <- c("response","MSE","MSE.pval","IncNodePurity","IncNodePurity.pval")

imp.rf.PCOA.hetero.select.axis1.rfp$Region <- 'Heterotrophic'

### combine all data 

select.rf <- rbind(imp.rf.PCOA.select.axis1.rfp, imp.rf.PCOA.photo.select.axis1.rfp, imp.rf.PCOA.chemo.select.axis1.rfp, imp.rf.PCOA.hetero.select.axis1.rfp)

## convert continuous to discrete values
select.rf$MSE.pval <- as.factor(select.rf$MSE.pval)

sort(unique(select.rf$MSE.pval)) 
#[1] 0.000999000999000999 0.001998001998002    0.002997002997003   
#[4] 0.004995004995005    0.00699300699300699  0.0579420579420579   


# Define the number of colors you want
select.cols <- 6
selectcolors <- colorRampPalette(brewer.pal(8, "Dark2"))(select.cols)

select.PCOA.axis1.MSE.bubbleplot <- select.rf  %>% ggplot( aes(x = reorder(response,MSE), y = MSE, shape = factor(Region), color = MSE.pval)) + geom_point(alpha = 0.5, size = 3) + coord_flip()   + labs( x= "Environemntal variables", y = "% Increase in mean squared error (%IncMSE)", fill = "Mean squared error p-value") +  theme_bw() + theme( axis.title = element_text(colour = "black", face = "bold", size = 10) , axis.text.x = element_text(colour = "black", size = 10, face = "bold", angle = 0, hjust = 0.5),   axis.text.y = element_text(colour = "black", face = "bold", size = 10), legend.title = element_text(size = 9, face = "bold"), panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2), legend.position = "right",plot.title = element_text(size=10))  + ggtitle("Factors influencing microbial community structure at select (PCoA Axis.1)")  + scale_shape_manual(values=c(15,16,17,18)) +  scale_color_manual(values = selectcolors)

ggsave("select.PCOA.axis1.MSE.bubbleplot.pdf", height = 5.5, width = 7)
ggsave("select.PCOA.axis1.MSE.bubbleplot.png", height = 5.5, width = 7)

#### multiple 4 plots for each group with MSE.pval as factor and facet #######

select.PCOA.axis1.MSE.facet.bubbleplot <- select.rf %>% 
  mutate(response = tidytext::reorder_within(response, MSE, within = Region))   %>% ggplot( aes(x = MSE, y = response, color = MSE.pval)) + geom_point(alpha = 0.5, size = 3)  + labs( x= "Environemntal variables", y = "% Increase in mean squared error (%IncMSE)", fill = "Mean squared error p-value") +  theme_bw() + theme( axis.title = element_text(colour = "black", face = "bold", size = 10) , axis.text.x = element_text(colour = "black", size = 10, face = "bold", angle = 0, hjust = 0.5),   axis.text.y = element_text(colour = "black", face = "bold", size = 10), legend.title = element_text(size = 9, face = "bold"), panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2), legend.position = "bottom",plot.title = element_text(size=10))  + ggtitle("Factors influencing microbial community structure at select (PCoA Axis.1)")  + scale_shape_manual(values=c(15,16,17,18)) +  scale_color_manual(values = selectcolors) + tidytext::scale_y_reordered() + facet_wrap(vars(Region),nrow = 1,  scales = "free_y")

ggsave("select.PCOA.axis1.MSE.facet.bubbleplot.pdf", height = 5.5, width = 12)
ggsave("select.PCOA.axis1.MSE.facet.bubbleplot.png", height = 5.5, width = 12)


########################

### Combine all regions and all grouops and Overall #####

overall <- select.rf
overall$Region2 <- 'Overall'

NT <- NT.rf
NT$Region2 <- 'North.Thailand'

NM <- NM.rf
NM$Region2 <- 'North.Malaysia'

ST <- ST.rf
ST$Region2 <- 'South.Thailand'

SM <- SM.rf
SM$Region2 <- 'South.Malaysia'

CT <- CT.rf
CT$Region2 <- 'Central.Thailand'

SG <- SG.rf
SG$Region2 <- 'Singapore'


### combine all data 

all.rf <- rbind(overall, NT, NM, ST, SM, CT, SG)

## convert continuous to discrete values
#all.rf$MSE.pval <- as.factor(all.rf$MSE.pval)

unique(sort(str_sort(all.rf$MSE.pval, numeric = TRUE)))
## 34 unique

# Define the number of colors you want
all.cols <- 34
allcolors <- colorRampPalette(brewer.pal(8, "Dark2"))(all.cols)

temp$size_f = factor(temp$size, levels=c('50%','100%','150%','200%'))

all.rf$Region2 = factor(all.rf$Region2, levels=c('Overall','North.Thailand','Central.Thailand','South.Thailand', 'North.Malaysia', 'South.Malaysia', 'Singapore'))

all.rf$Region = factor(all.rf$Region, levels=c('All taxa','Photosynthetic','Chemoautotrophic','Heterotrophic'))

all.PCOA.axis1.MSE.facet.bubbleplot <- all.rf %>%
  mutate(response = tidytext::reorder_within(response, MSE, within = Region2))   %>% ggplot( aes(x = MSE, y = response, color = MSE.pval)) + geom_point(alpha = 0.5, size = 3)  + labs( x= "% Increase in mean squared error (%IncMSE)", y = "Environmental variables", fill = "Mean squared error p-value") +  theme_bw() + theme( axis.title = element_text(colour = "black", face = "bold", size = 10) , axis.text.x = element_text(colour = "black", size = 10, face = "bold", angle = 0, hjust = 0.5),   axis.text.y = element_text(colour = "black", face = "bold", size = 10), legend.title = element_text(size = 9, face = "bold"), panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2), legend.position = "right",plot.title = element_text(size=10))  + ggtitle("Factors influencing microbial community structure (PCoA Axis.1)")  + scale_shape_manual(values=c(15,16,17,18)) +  scale_color_manual(values = allcolors, breaks = unique(sort(str_sort(all.rf$MSE.pval, numeric = TRUE)))) + tidytext::scale_y_reordered() + facet_grid(Region2 ~ Region ,scales = "free_y") 

#facet_grid(vars(factor(Region2, levels=c('Overall','North Thailand','Central Thailand','South Thailand', 'North Malaysia', 'South Malaysia', 'Singapore'))), vars(factor(Region, levels=c('All taxa','Photosynthetic','Chemoautotrophic','Heterotrophic'))) ,scales = "free_y") 

#breaks = unique(sort(str_sort(all.rf$MSE.pval, numeric = TRUE)))

ggsave("all.PCOA.axis1.MSE.facet.bubbleplot.pdf", height = 10, width = 16)
ggsave("all.PCOA.axis1.MSE.facet.bubbleplot.png", height = 10, width = 16)


all.rf2 <- all.rf
all.rf2$MSE.pval <- as.numeric(as.character(all.rf2$MSE.pval))

all.rf2 %<>% dplyr::mutate(MSE.pval2 = cut(MSE.pval, breaks = c(-Inf, 0.001, 0.0025,0.005,0.0075,0.01,0.025,0.05,0.075,0.1,0.25, 0.5, Inf), labels = c("< 0.001", "0.001 - 0.0025", "0.0025 - 0.005", "0.005 - 0.0075", "0.0075 - 0.01", "0.01 - 0.025", "0.025 - 0.050", "0.050 - 0.075", "0.075 - 0.1", "0.1 - 0.25", "0.25 - 0.50",">= 0.5")))

unique(sort(str_sort(all.rf2$MSE.pval2, numeric = TRUE)))
## 12 unique

# Define the number of colors you want
all.cols.2 <- 12
allcolors2 <- colorRampPalette(brewer.pal(12, "Dark2"))(all.cols.2)

all.rf2$Region2 = factor(all.rf$Region2, levels=c('Overall','North.Thailand','Central.Thailand','South.Thailand', 'North.Malaysia', 'South.Malaysia', 'Singapore'))

all.rf2$Region = factor(all.rf$Region, levels=c('All taxa','Photosynthetic','Chemoautotrophic','Heterotrophic'))

all.PCOA.axis1.MSE.facet.bubbleplot.2 <- all.rf2 %>%
  mutate(response = tidytext::reorder_within(response, MSE, within = Region2))   %>% ggplot( aes(x = MSE, y = response, color = MSE.pval2)) + geom_point(alpha = 0.5, size = 3)  + labs( x= "% Increase in mean squared error (%IncMSE)", y = "Environmental variables", fill = "Mean squared error p-value") +  theme_bw() + theme( axis.title = element_text(colour = "black", face = "bold", size = 10) , axis.text.x = element_text(colour = "black", size = 10, face = "bold", angle = 0, hjust = 0.5),   axis.text.y = element_text(colour = "black", face = "bold", size = 10), legend.title = element_text(size = 9, face = "bold"), panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2), legend.position = "right",plot.title = element_text(size=10))  + ggtitle("Factors influencing microbial community structure (PCoA Axis.1)")  + scale_shape_manual(values=c(15,16,17,18)) +  scale_color_manual(values = allcolors2) + tidytext::scale_y_reordered() + facet_grid(Region2 ~ Region ,scales = "free_y") 

ggsave("all.PCOA.axis1.MSE.facet.bubbleplot.2.pdf", height = 10, width = 16)
ggsave("all.PCOA.axis1.MSE.facet.bubbleplot.2.png", height = 10, width = 16)




##############################

rfp.PCOA.all.reg <- readRDS("rfp.PCOA.all.reg.rds")
rfp.PCOA.all.reg.sig <- readRDS("rfp.PCOA.all.reg.sig.rds")

imp.rfp.PCOA.all.reg <- as.data.frame(importance(rfp.PCOA.all.reg ))
imp.rfp.PCOA.all.reg <- imp.rfp.PCOA.all.reg %>% rownames_to_column(var="response")

colnames(imp.rfp.PCOA.all.reg) <- c("response","MSE","MSE.pval","IncNodePurity","IncNodePurity.pval")

variable.imp.PCOA.axis1.all.reg.MSE.bubbleplot <- imp.rfp.PCOA.all.reg  %>% ggplot( aes(x = reorder(response,MSE), y = MSE)) + geom_point(aes(fill = MSE.pval), alpha = 0.5, shape = 21,size = 6) + coord_flip()  + scale_size_continuous(limits = c(0,0.15), range = c(0,0.15), breaks = c(0,0.05,0.10,0.15)) + labs( x= "Environemntal variables", y = "% Increase in mean squared error (%IncMSE)", fill = "Mean squared error p-value") + theme(legend.key=element_blank(), axis.title = element_text(colour = "black", face = "bold", size = 8) , axis.text.x = element_text(colour = "black", size = 12, face = "bold", angle = 0, hjust = 0.5),   axis.text.y = element_text(colour = "black", face = "bold", size = 12), legend.title = element_text(size = 10, face = "bold"), panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2), legend.position = "right", title = element_text(size = 8))  + ggtitle("Factors influencing all taxa (PCoA Axis.1) with region")  + scale_fill_viridis() + theme_bw()

ggsave("variable.imp.PCOA.axis1.all.reg.MSE.bubbleplot.pdf", height = 7, width = 9)
ggsave("variable.imp.PCOA.axis1.all.reg.MSE.bubbleplot.png", height = 8, width = 7.5)


variable.imp.PCOA.axis1.all.reg.IncNodePurity.bubbleplot <- imp.rfp.PCOA.all.reg  %>% ggplot( aes(x = reorder(response,IncNodePurity), y = IncNodePurity)) + geom_point(aes(fill = IncNodePurity.pval), alpha = 0.5, shape = 21,size = 6) + coord_flip()  + scale_size_continuous(limits = c(0,0.15), range = c(0,0.15), breaks = c(0,0.05,0.10,0.15)) + labs( x= "Environemntal variables", y = "Increase in node purity (IncNodePurity)", fill = "Increase in node purity p-value") + theme(legend.key=element_blank(), axis.title = element_text(colour = "black", face = "bold", size = 8) , axis.text.x = element_text(colour = "black", size = 12, face = "bold", angle = 0, hjust = 0.5),   axis.text.y = element_text(colour = "black", face = "bold", size = 12), legend.title = element_text(size = 10, face = "bold"), panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2), legend.position = "right", title = element_text(size = 8))  + ggtitle("Factors influencing all taxa (PCoA Axis.1) with region")  + scale_fill_viridis() + theme_bw()

ggsave("variable.imp.PCOA.axis1.all.reg.IncNodePurity.bubbleplot.pdf", height = 7, width = 9)
ggsave("variable.imp.PCOA.axis1.all.reg.IncNodePurity.bubbleplot.png", height = 8, width = 7.5)


###

rfp.PCOA.all.lat <- readRDS("rfp.PCOA.all.lat.rds")
rfp.PCOA.all.lat.sig <- readRDS("rfp.PCOA.all.lat.sig.rds")

imp.rfp.PCOA.all.lat <- as.data.frame(importance(rfp.PCOA.all.lat ))
imp.rfp.PCOA.all.lat <- imp.rfp.PCOA.all.lat %>% rownames_to_column(var="response")

colnames(imp.rfp.PCOA.all.lat) <- c("response","MSE","MSE.pval","IncNodePurity","IncNodePurity.pval")

variable.imp.PCOA.axis1.all.lat.MSE.bubbleplot <- imp.rfp.PCOA.all.lat  %>% ggplot( aes(x = reorder(response,MSE), y = MSE)) + geom_point(aes(fill = MSE.pval), alpha = 0.5, shape = 21,size = 6) + coord_flip()  + scale_size_continuous(limits = c(0,0.15), range = c(0,0.15), breaks = c(0,0.05,0.10,0.15)) + labs( x= "Environemntal variables", y = "% Increase in mean squared error (%IncMSE)", fill = "Mean squared error p-value") + theme(legend.key=element_blank(), axis.title = element_text(colour = "black", face = "bold", size = 8) , axis.text.x = element_text(colour = "black", size = 12, face = "bold", angle = 0, hjust = 0.5),   axis.text.y = element_text(colour = "black", face = "bold", size = 12), legend.title = element_text(size = 10, face = "bold"), panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2), legend.position = "right", title = element_text(size = 8))  + ggtitle("Factors influencing all taxa (PCoA Axis.1) with latitude")  + scale_fill_viridis() + theme_bw()

ggsave("variable.imp.PCOA.axis1.all.lat.MSE.bubbleplot.pdf", height = 7, width = 9)
ggsave("variable.imp.PCOA.axis1.all.lat.MSE.bubbleplot.png", height = 8, width = 7.5)


variable.imp.PCOA.axis1.all.lat.IncNodePurity.bubbleplot <- imp.rfp.PCOA.all.lat  %>% ggplot( aes(x = reorder(response,IncNodePurity), y = IncNodePurity)) + geom_point(aes(fill = IncNodePurity.pval), alpha = 0.5, shape = 21,size = 6) + coord_flip()  + scale_size_continuous(limits = c(0,0.15), range = c(0,0.15), breaks = c(0,0.05,0.10,0.15)) + labs( x= "Environemntal variables", y = "Increase in node purity (IncNodePurity)", fill = "Increase in node purity p-value") + theme(legend.key=element_blank(), axis.title = element_text(colour = "black", face = "bold", size = 8) , axis.text.x = element_text(colour = "black", size = 12, face = "bold", angle = 0, hjust = 0.5),   axis.text.y = element_text(colour = "black", face = "bold", size = 12), legend.title = element_text(size = 10, face = "bold"), panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2), legend.position = "right", title = element_text(size = 8))  + ggtitle("Factors influencing all taxa (PCoA Axis.1) with latitude")  + scale_fill_viridis() + theme_bw()

ggsave("variable.imp.PCOA.axis1.all.lat.IncNodePurity.bubbleplot.pdf", height = 7, width = 9)
ggsave("variable.imp.PCOA.axis1.all.lat.IncNodePurity.bubbleplot.png", height = 8, width = 7.5)


###

rfp.PCOA.photo.reg <- readRDS("rfp.PCOA.photo.reg.rds")
rfp.PCOA.photo.reg.sig <- readRDS("rfp.PCOA.photo.reg.sig.rds")

imp.rfp.PCOA.photo.reg <- as.data.frame(importance(rfp.PCOA.photo.reg ))
imp.rfp.PCOA.photo.reg <- imp.rfp.PCOA.photo.reg %>% rownames_to_column(var="response")

colnames(imp.rfp.PCOA.photo.reg) <- c("response","MSE","MSE.pval","IncNodePurity","IncNodePurity.pval")

variable.imp.PCOA.axis1.photo.reg.MSE.bubbleplot <- imp.rfp.PCOA.photo.reg  %>% ggplot( aes(x = reorder(response,MSE), y = MSE)) + geom_point(aes(fill = MSE.pval), alpha = 0.5, shape = 21,size = 6) + coord_flip()  + scale_size_continuous(limits = c(0,0.15), range = c(0,0.15), breaks = c(0,0.05,0.10,0.15)) + labs( x= "Environemntal variables", y = "% Increase in mean squared error (%IncMSE)", fill = "Mean squared error p-value") + theme(legend.key=element_blank(), axis.title = element_text(colour = "black", face = "bold", size = 8) , axis.text.x = element_text(colour = "black", size = 12, face = "bold", angle = 0, hjust = 0.5),   axis.text.y = element_text(colour = "black", face = "bold", size = 12), legend.title = element_text(size = 10, face = "bold"), panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2), legend.position = "right", title = element_text(size = 8))  + ggtitle("Factors influencing photosynthetic taxa (PCoA Axis.1) with region")  + scale_fill_viridis() + theme_bw()

ggsave("variable.imp.PCOA.axis1.photo.reg.MSE.bubbleplot.pdf", height = 7, width = 9)
ggsave("variable.imp.PCOA.axis1.photo.reg.MSE.bubbleplot.png", height = 8, width = 7.5)


variable.imp.PCOA.axis1.photo.reg.IncNodePurity.bubbleplot <- imp.rfp.PCOA.photo.reg  %>% ggplot( aes(x = reorder(response,IncNodePurity), y = IncNodePurity)) + geom_point(aes(fill = IncNodePurity.pval), alpha = 0.5, shape = 21,size = 6) + coord_flip()  + scale_size_continuous(limits = c(0,0.15), range = c(0,0.15), breaks = c(0,0.05,0.10,0.15)) + labs( x= "Environemntal variables", y = "Increase in node purity (IncNodePurity)", fill = "Increase in node purity p-value") + theme(legend.key=element_blank(), axis.title = element_text(colour = "black", face = "bold", size = 8) , axis.text.x = element_text(colour = "black", size = 12, face = "bold", angle = 0, hjust = 0.5),   axis.text.y = element_text(colour = "black", face = "bold", size = 12), legend.title = element_text(size = 10, face = "bold"), panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2), legend.position = "right", title = element_text(size = 8))  + ggtitle("Factors influencing photosynthetic taxa (PCoA Axis.1) with region")  + scale_fill_viridis() + theme_bw()

ggsave("variable.imp.PCOA.axis1.photo.reg.IncNodePurity.bubbleplot.pdf", height = 7, width = 9)
ggsave("variable.imp.PCOA.axis1.photo.reg.IncNodePurity.bubbleplot.png", height = 8, width = 7.5)

###

rfp.PCOA.photo.lat <- readRDS("rfp.PCOA.photo.lat.rds")
rfp.PCOA.photo.lat.sig <- readRDS("rfp.PCOA.photo.lat.sig.rds")

imp.rfp.PCOA.photo.lat <- as.data.frame(importance(rfp.PCOA.photo.lat ))
imp.rfp.PCOA.photo.lat <- imp.rfp.PCOA.photo.lat %>% rownames_to_column(var="response")

colnames(imp.rfp.PCOA.photo.lat) <- c("response","MSE","MSE.pval","IncNodePurity","IncNodePurity.pval")

variable.imp.PCOA.axis1.photo.lat.MSE.bubbleplot <- imp.rfp.PCOA.photo.lat  %>% ggplot( aes(x = reorder(response,MSE), y = MSE)) + geom_point(aes(fill = MSE.pval), alpha = 0.5, shape = 21,size = 6) + coord_flip()  + scale_size_continuous(limits = c(0,0.15), range = c(0,0.15), breaks = c(0,0.05,0.10,0.15)) + labs( x= "Environemntal variables", y = "% Increase in mean squared error (%IncMSE)", fill = "Mean squared error p-value") + theme(legend.key=element_blank(), axis.title = element_text(colour = "black", face = "bold", size = 8) , axis.text.x = element_text(colour = "black", size = 12, face = "bold", angle = 0, hjust = 0.5),   axis.text.y = element_text(colour = "black", face = "bold", size = 12), legend.title = element_text(size = 10, face = "bold"), panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2), legend.position = "right", title = element_text(size = 8))  + ggtitle("Factors influencing photosynthetic taxa (PCoA Axis.1) with latitude")  + scale_fill_viridis() + theme_bw()

ggsave("variable.imp.PCOA.axis1.photo.lat.MSE.bubbleplot.pdf", height = 7, width = 9)
ggsave("variable.imp.PCOA.axis1.photo.lat.MSE.bubbleplot.png", height = 8, width = 7.5)


variable.imp.PCOA.axis1.photo.lat.IncNodePurity.bubbleplot <- imp.rfp.PCOA.photo.lat  %>% ggplot( aes(x = reorder(response,IncNodePurity), y = IncNodePurity)) + geom_point(aes(fill = IncNodePurity.pval), alpha = 0.5, shape = 21,size = 6) + coord_flip()  + scale_size_continuous(limits = c(0,0.15), range = c(0,0.15), breaks = c(0,0.05,0.10,0.15)) + labs( x= "Environemntal variables", y = "Increase in node purity (IncNodePurity)", fill = "Increase in node purity p-value") + theme(legend.key=element_blank(), axis.title = element_text(colour = "black", face = "bold", size = 8) , axis.text.x = element_text(colour = "black", size = 12, face = "bold", angle = 0, hjust = 0.5),   axis.text.y = element_text(colour = "black", face = "bold", size = 12), legend.title = element_text(size = 10, face = "bold"), panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2), legend.position = "right", title = element_text(size = 8))  + ggtitle("Factors influencing photosynhetic taxa (PCoA Axis.1) with latitude")  + scale_fill_viridis() + theme_bw()

ggsave("variable.imp.PCOA.axis1.photo.lat.IncNodePurity.bubbleplot.pdf", height = 7, width = 9)
ggsave("variable.imp.PCOA.axis1.photo.lat.IncNodePurity.bubbleplot.png", height = 8, width = 7.5)


###

rfp.PCOA.chemo.reg <- readRDS("rfp.PCOA.chemo.reg.rds")
rfp.PCOA.chemo.reg.sig <- readRDS("rfp.PCOA.chemo.reg.sig.rds")

imp.rfp.PCOA.chemo.reg <- as.data.frame(importance(rfp.PCOA.chemo.reg ))
imp.rfp.PCOA.chemo.reg <- imp.rfp.PCOA.chemo.reg %>% rownames_to_column(var="response")

colnames(imp.rfp.PCOA.chemo.reg) <- c("response","MSE","MSE.pval","IncNodePurity","IncNodePurity.pval")

variable.imp.PCOA.axis1.chemo.reg.MSE.bubbleplot <- imp.rfp.PCOA.chemo.reg  %>% ggplot( aes(x = reorder(response,MSE), y = MSE)) + geom_point(aes(fill = MSE.pval), alpha = 0.5, shape = 21,size = 6) + coord_flip()  + scale_size_continuous(limits = c(0,0.15), range = c(0,0.15), breaks = c(0,0.05,0.10,0.15)) + labs( x= "Environemntal variables", y = "% Increase in mean squared error (%IncMSE)", fill = "Mean squared error p-value") + theme(legend.key=element_blank(), axis.title = element_text(colour = "black", face = "bold", size = 8) , axis.text.x = element_text(colour = "black", size = 12, face = "bold", angle = 0, hjust = 0.5),   axis.text.y = element_text(colour = "black", face = "bold", size = 12), legend.title = element_text(size = 10, face = "bold"), panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2), legend.position = "right", title = element_text(size = 8))  + ggtitle("Factors influencing chemoautotrophs (PCoA Axis.1) with region")  + scale_fill_viridis() + theme_bw()

ggsave("variable.imp.PCOA.axis1.chemo.reg.MSE.bubbleplot.pdf", height = 7, width = 9)
ggsave("variable.imp.PCOA.axis1.chemo.reg.MSE.bubbleplot.png", height = 8, width = 7.5)


variable.imp.PCOA.axis1.chemo.reg.IncNodePurity.bubbleplot <- imp.rfp.PCOA.chemo.reg  %>% ggplot( aes(x = reorder(response,IncNodePurity), y = IncNodePurity)) + geom_point(aes(fill = IncNodePurity.pval), alpha = 0.5, shape = 21,size = 6) + coord_flip()  + scale_size_continuous(limits = c(0,0.15), range = c(0,0.15), breaks = c(0,0.05,0.10,0.15)) + labs( x= "Environemntal variables", y = "Increase in node purity (IncNodePurity)", fill = "Increase in node purity p-value") + theme(legend.key=element_blank(), axis.title = element_text(colour = "black", face = "bold", size = 8) , axis.text.x = element_text(colour = "black", size = 12, face = "bold", angle = 0, hjust = 0.5),   axis.text.y = element_text(colour = "black", face = "bold", size = 12), legend.title = element_text(size = 10, face = "bold"), panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2), legend.position = "right", title = element_text(size = 8))  + ggtitle("Factors influencing chemoautotrophs (PCoA Axis.1) with region")  + scale_fill_viridis() + theme_bw()

ggsave("variable.imp.PCOA.axis1.chemo.reg.IncNodePurity.bubbleplot.pdf", height = 7, width = 9)
ggsave("variable.imp.PCOA.axis1.chemo.reg.IncNodePurity.bubbleplot.png", height = 8, width = 7.5)

###

rfp.PCOA.chemo.lat <- readRDS("rfp.PCOA.chemo.lat.rds")
rfp.PCOA.chemo.lat.sig <- readRDS("rfp.PCOA.chemo.lat.sig.rds")

imp.rfp.PCOA.chemo.lat <- as.data.frame(importance(rfp.PCOA.chemo.lat ))
imp.rfp.PCOA.chemo.lat <- imp.rfp.PCOA.chemo.lat %>% rownames_to_column(var="response")

colnames(imp.rfp.PCOA.chemo.lat) <- c("response","MSE","MSE.pval","IncNodePurity","IncNodePurity.pval")

variable.imp.PCOA.axis1.chemo.lat.MSE.bubbleplot <- imp.rfp.PCOA.chemo.lat  %>% ggplot( aes(x = reorder(response,MSE), y = MSE)) + geom_point(aes(fill = MSE.pval), alpha = 0.5, shape = 21,size = 6) + coord_flip()  + scale_size_continuous(limits = c(0,0.15), range = c(0,0.15), breaks = c(0,0.05,0.10,0.15)) + labs( x= "Environemntal variables", y = "% Increase in mean squared error (%IncMSE)", fill = "Mean squared error p-value") + theme(legend.key=element_blank(), axis.title = element_text(colour = "black", face = "bold", size = 8) , axis.text.x = element_text(colour = "black", size = 12, face = "bold", angle = 0, hjust = 0.5),   axis.text.y = element_text(colour = "black", face = "bold", size = 12), legend.title = element_text(size = 10, face = "bold"), panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2), legend.position = "right", title = element_text(size = 8))  + ggtitle("Factors influencing chemoautotrophs (PCoA Axis.1) with latitude")  + scale_fill_viridis() + theme_bw()

ggsave("variable.imp.PCOA.axis1.chemo.lat.MSE.bubbleplot.pdf", height = 7, width = 9)
ggsave("variable.imp.PCOA.axis1.chemo.lat.MSE.bubbleplot.png", height = 8, width = 7.5)


variable.imp.PCOA.axis1.chemo.lat.IncNodePurity.bubbleplot <- imp.rfp.PCOA.chemo.lat  %>% ggplot( aes(x = reorder(response,IncNodePurity), y = IncNodePurity)) + geom_point(aes(fill = IncNodePurity.pval), alpha = 0.5, shape = 21,size = 6) + coord_flip()  + scale_size_continuous(limits = c(0,0.15), range = c(0,0.15), breaks = c(0,0.05,0.10,0.15)) + labs( x= "Environemntal variables", y = "Increase in node purity (IncNodePurity)", fill = "Increase in node purity p-value") + theme(legend.key=element_blank(), axis.title = element_text(colour = "black", face = "bold", size = 8) , axis.text.x = element_text(colour = "black", size = 12, face = "bold", angle = 0, hjust = 0.5),   axis.text.y = element_text(colour = "black", face = "bold", size = 12), legend.title = element_text(size = 10, face = "bold"), panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2), legend.position = "right", title = element_text(size = 8))  + ggtitle("Factors influencing chemoautotrophs (PCoA Axis.1) with latitude")  + scale_fill_viridis() + theme_bw()

ggsave("variable.imp.PCOA.axis1.chemo.lat.IncNodePurity.bubbleplot.pdf", height = 7, width = 9)
ggsave("variable.imp.PCOA.axis1.chemo.lat.IncNodePurity.bubbleplot.png", height = 8, width = 7.5)


###

rfp.PCOA.hetero.reg <- readRDS("rfp.PCOA.hetero.reg.rds")
rfp.PCOA.hetero.reg.sig <- readRDS("rfp.PCOA.hetero.reg.sig.rds")

imp.rfp.PCOA.hetero.reg <- as.data.frame(importance(rfp.PCOA.hetero.reg ))
imp.rfp.PCOA.hetero.reg <- imp.rfp.PCOA.hetero.reg %>% rownames_to_column(var="response")

colnames(imp.rfp.PCOA.hetero.reg) <- c("response","MSE","MSE.pval","IncNodePurity","IncNodePurity.pval")

variable.imp.PCOA.axis1.hetero.reg.MSE.bubbleplot <- imp.rfp.PCOA.hetero.reg  %>% ggplot( aes(x = reorder(response,MSE), y = MSE)) + geom_point(aes(fill = MSE.pval), alpha = 0.5, shape = 21,size = 6) + coord_flip()  + scale_size_continuous(limits = c(0,0.15), range = c(0,0.15), breaks = c(0,0.05,0.10,0.15)) + labs( x= "Environemntal variables", y = "% Increase in mean squared error (%IncMSE)", fill = "Mean squared error p-value") + theme(legend.key=element_blank(), axis.title = element_text(colour = "black", face = "bold", size = 8) , axis.text.x = element_text(colour = "black", size = 12, face = "bold", angle = 0, hjust = 0.5),   axis.text.y = element_text(colour = "black", face = "bold", size = 12), legend.title = element_text(size = 10, face = "bold"), panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2), legend.position = "right", title = element_text(size = 8))  + ggtitle("Factors influencing heterotrophs (PCoA Axis.1) with region")  + scale_fill_viridis() + theme_bw()

ggsave("variable.imp.PCOA.axis1.hetero.reg.MSE.bubbleplot.pdf", height = 7, width = 9)
ggsave("variable.imp.PCOA.axis1.hetero.reg.MSE.bubbleplot.png", height = 8, width = 7.5)


variable.imp.PCOA.axis1.hetero.reg.IncNodePurity.bubbleplot <- imp.rfp.PCOA.hetero.reg  %>% ggplot( aes(x = reorder(response,IncNodePurity), y = IncNodePurity)) + geom_point(aes(fill = IncNodePurity.pval), alpha = 0.5, shape = 21,size = 6) + coord_flip()  + scale_size_continuous(limits = c(0,0.15), range = c(0,0.15), breaks = c(0,0.05,0.10,0.15)) + labs( x= "Environemntal variables", y = "Increase in node purity (IncNodePurity)", fill = "Increase in node purity p-value") + theme(legend.key=element_blank(), axis.title = element_text(colour = "black", face = "bold", size = 8) , axis.text.x = element_text(colour = "black", size = 12, face = "bold", angle = 0, hjust = 0.5),   axis.text.y = element_text(colour = "black", face = "bold", size = 12), legend.title = element_text(size = 10, face = "bold"), panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2), legend.position = "right", title = element_text(size = 8))  + ggtitle("Factors influencing heterotrophs (PCoA Axis.1) with region")  + scale_fill_viridis() + theme_bw()

ggsave("variable.imp.PCOA.axis1.hetero.reg.IncNodePurity.bubbleplot.pdf", height = 7, width = 9)
ggsave("variable.imp.PCOA.axis1.hetero.reg.IncNodePurity.bubbleplot.png", height = 8, width = 7.5)

###

rfp.PCOA.hetero.lat <- readRDS("rfp.PCOA.hetero.lat.rds")
rfp.PCOA.hetero.lat.sig <- readRDS("rfp.PCOA.hetero.lat.sig.rds")

imp.rfp.PCOA.hetero.lat <- as.data.frame(importance(rfp.PCOA.hetero.lat ))
imp.rfp.PCOA.hetero.lat <- imp.rfp.PCOA.hetero.lat %>% rownames_to_column(var="response")

colnames(imp.rfp.PCOA.hetero.lat) <- c("response","MSE","MSE.pval","IncNodePurity","IncNodePurity.pval")

variable.imp.PCOA.axis1.hetero.lat.MSE.bubbleplot <- imp.rfp.PCOA.hetero.lat  %>% ggplot( aes(x = reorder(response,MSE), y = MSE)) + geom_point(aes(fill = MSE.pval), alpha = 0.5, shape = 21,size = 6) + coord_flip()  + scale_size_continuous(limits = c(0,0.15), range = c(0,0.15), breaks = c(0,0.05,0.10,0.15)) + labs( x= "Environemntal variables", y = "% Increase in mean squared error (%IncMSE)", fill = "Mean squared error p-value") + theme(legend.key=element_blank(), axis.title = element_text(colour = "black", face = "bold", size = 8) , axis.text.x = element_text(colour = "black", size = 12, face = "bold", angle = 0, hjust = 0.5),   axis.text.y = element_text(colour = "black", face = "bold", size = 12), legend.title = element_text(size = 10, face = "bold"), panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2), legend.position = "right", title = element_text(size = 8))  + ggtitle("Factors influencing heterotrophs (PCoA Axis.1) with latitude")  + scale_fill_viridis() + theme_bw()

ggsave("variable.imp.PCOA.axis1.hetero.lat.MSE.bubbleplot.pdf", height = 7, width = 9)
ggsave("variable.imp.PCOA.axis1.hetero.lat.MSE.bubbleplot.png", height = 8, width = 7.5)


variable.imp.PCOA.axis1.hetero.lat.IncNodePurity.bubbleplot <- imp.rfp.PCOA.hetero.lat  %>% ggplot( aes(x = reorder(response,IncNodePurity), y = IncNodePurity)) + geom_point(aes(fill = IncNodePurity.pval), alpha = 0.5, shape = 21,size = 6) + coord_flip()  + scale_size_continuous(limits = c(0,0.15), range = c(0,0.15), breaks = c(0,0.05,0.10,0.15)) + labs( x= "Environemntal variables", y = "Increase in node purity (IncNodePurity)", fill = "Increase in node purity p-value") + theme(legend.key=element_blank(), axis.title = element_text(colour = "black", face = "bold", size = 8) , axis.text.x = element_text(colour = "black", size = 12, face = "bold", angle = 0, hjust = 0.5),   axis.text.y = element_text(colour = "black", face = "bold", size = 12), legend.title = element_text(size = 10, face = "bold"), panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2), legend.position = "right", title = element_text(size = 8))  + ggtitle("Factors influencing heterotrophs (PCoA Axis.1) with latitude")  + scale_fill_viridis() + theme_bw()

ggsave("variable.imp.PCOA.axis1.hetero.lat.IncNodePurity.bubbleplot.pdf", height = 7, width = 9)
ggsave("variable.imp.PCOA.axis1.hetero.lat.IncNodePurity.bubbleplot.png", height = 8, width = 7.5)

####

rfp.PCOA.select <- readRDS("rfp.PCOA.select.rds")
rfp.PCOA.select.sig <- readRDS("rfp.PCOA.select.sig.rds")

imp.rfp.PCOA.select <- as.data.frame(importance(rfp.PCOA.select))
imp.rfp.PCOA.select <- imp.rfp.PCOA.select %>% rownames_to_column(var="response")

colnames(imp.rfp.PCOA.select) <- c("response","MSE","MSE.pval","IncNodePurity","IncNodePurity.pval")


rfp.PCOA.photo.select <- readRDS("rfp.PCOA.photo.select.rds")
rfp.PCOA.photo.select.sig <- readRDS("rfp.PCOA.photo.select.sig.rds")

imp.rfp.PCOA.photo.select <- as.data.frame(importance(rfp.PCOA.photo.select))
imp.rfp.PCOA.photo.select <- imp.rfp.PCOA.photo.select %>% rownames_to_column(var="response")

colnames(imp.rfp.PCOA.photo.select) <- c("response","MSE","MSE.pval","IncNodePurity","IncNodePurity.pval")


rfp.PCOA.chemo.select <- readRDS("rfp.PCOA.chemo.select.rds")
rfp.PCOA.chemo.select.sig <- readRDS("rfp.PCOA.chemo.select.sig.rds")

imp.rfp.PCOA.chemo.select <- as.data.frame(importance(rfp.PCOA.chemo.select))
imp.rfp.PCOA.chemo.select <- imp.rfp.PCOA.chemo.select %>% rownames_to_column(var="response")

colnames(imp.rfp.PCOA.chemo.select) <- c("response","MSE","MSE.pval","IncNodePurity","IncNodePurity.pval")


rfp.PCOA.hetero.select <- readRDS("rfp.PCOA.hetero.select.rds")
rfp.PCOA.hetero.select.sig <- readRDS("rfp.PCOA.hetero.select.sig.rds")

imp.rfp.PCOA.hetero.select <- as.data.frame(importance(rfp.PCOA.hetero.select))
imp.rfp.PCOA.hetero.select <- imp.rfp.PCOA.hetero.select %>% rownames_to_column(var="response")

colnames(imp.rfp.PCOA.hetero.select) <- c("response","MSE","MSE.pval","IncNodePurity","IncNodePurity.pval")



#########################################################################################

### Fig. S9. Autocorrelation of abiotic variables by biogeogrphic regions ######
### CHECK

######################################################################################

pdf("autocorr.plot.regions.pdf", height = 14 ,width = 22)
microeco.mantel$cal_autocor(group = "Region")
dev.off()


#########################################################################################

### Fig. S10. Univariate (Mantel’s Test) and multivariate (Variance Partitioning) correlations for abiotic variables with microbial community distribution for individual regions. ######
### CHECK

#### Mantels test by regions and groups ####

## import regions dataset-CHANGE !!!!!!!!!!!!!!!

NT.microeco.dataset <- readRDS("NT.microeco.dataset.rds")
NM.microeco.dataset <- readRDS("NM.microeco.dataset.rds")
ST.microeco.dataset <- readRDS("ST.microeco.dataset.rds")
SM.microeco.dataset <- readRDS("SM.microeco.dataset.rds")
CT.microeco.dataset <- readRDS("CT.microeco.dataset.rds")
SG.microeco.dataset <- readRDS("SG.microeco.dataset.rds")


#unifrac default FALSE; whether UniFrac index should be calculated,binary default FALSE; TRUE is used for jaccard and unweighted unifrac;
#microeco.dataset.2$cal_betadiv(unifrac = TRUE, binary = TRUE)
NT.microeco.dataset$cal_betadiv(unifrac = TRUE, binary = TRUE)
NM.microeco.dataset$cal_betadiv(unifrac = TRUE, binary = TRUE)
ST.microeco.dataset$cal_betadiv(unifrac = TRUE, binary = TRUE)
SM.microeco.dataset$cal_betadiv(unifrac = TRUE, binary = TRUE)
CT.microeco.dataset$cal_betadiv(unifrac = TRUE, binary = TRUE)
SG.microeco.dataset$cal_betadiv(unifrac = TRUE, binary = TRUE)

# first perform mantel test 
# Nitrate and Nitrites excluded as they were zero
# all other env variables were considered as indication of human activity so not included.
# Total alkalinity covaries with Carbonate so removed
# no lat and long

NT.microeco.mantel.select.noalk <- trans_env$new(dataset = NT.microeco.dataset, env_cols = 23:28)
NM.microeco.mantel.select.noalk <- trans_env$new(dataset = NM.microeco.dataset, env_cols = 23:28)
ST.microeco.mantel.select.noalk <- trans_env$new(dataset = ST.microeco.dataset, env_cols = 23:28)
SM.microeco.mantel.select.noalk <- trans_env$new(dataset = SM.microeco.dataset, env_cols = 23:28)
CT.microeco.mantel.select.noalk <- trans_env$new(dataset = CT.microeco.dataset, env_cols = 23:28)
SG.microeco.mantel.select.noalk <- trans_env$new(dataset = SG.microeco.dataset, env_cols = 23:28)

colnames(SG.microeco.mantel.select.noalk$data_env)


#microeco.mantel$data_env
#test <- NT.microeco.mantel.select.noalk$data_env
#lapply(test, class)

#http://qiime.org/tutorials/distance_matrix_comparison.html#:~:text=The%20partial%20Mantel%20test%20is,three%20distance%20(dissimilarit
#The partial Mantel test is used to estimate the correlation between two matrices, A and B, while controlling for the effect of a control matrix C. The partial Mantel test is a first-order correlation analysis that utilizes three distance (dissimilarity) matrices. This test builds on the simple Mantel test by adding a third “control” matrix. The goal is to test the correlation between matrices A and B while controlling the effect of a third matrix C, in order to remove spurious correlations. The first distance matrix is the one that is permuted so that the correlation structure between the first and second distance matrices is kept constant (Oksanen et al., 2011). A popular use of the partial Mantel test is to compare a community distance matrix with another distance matrix derived from an environmental parameter, using geographic distance as the third “control” distance matrix.

# Full Mantel's test between beta diversity matrix and environmental data.

NT.microeco.mantel.select.noalk$cal_mantel(use_measure = "unwei_unifrac", method = "pearson" , partial_mantel = FALSE)

NM.microeco.mantel.select.noalk$cal_mantel(use_measure = "unwei_unifrac", method = "pearson" , partial_mantel = FALSE)

ST.microeco.mantel.select.noalk$cal_mantel(use_measure = "unwei_unifrac", method = "pearson", partial_mantel = FALSE)

SM.microeco.mantel.select.noalk$cal_mantel(use_measure = "unwei_unifrac", method = "pearson" , partial_mantel = FALSE)

CT.microeco.mantel.select.noalk$cal_mantel(use_measure = "unwei_unifrac", method = "pearson" , partial_mantel = FALSE)

SG.microeco.mantel.select.noalk$cal_mantel(use_measure = "unwei_unifrac", method = "pearson" , partial_mantel = FALSE)

#test$cal_mantel(use_measure = "unwei_unifrac", method = "pearson" , partial_mantel = TRUE)

### use microeco.mantel.all to subset for required columns

test <- SG.microeco.mantel.select.noalk$data_env

## all variables with non-zero values
NT.microeco.mantel.select.noalk.all <- data.frame(NT.microeco.mantel.select.noalk$res_mantel)

## H2S values are 0 throughout and cannot run full Mantel's
NM.microeco.mantel.select.noalk.all <- data.frame(NM.microeco.mantel.select.noalk$res_mantel)

## Phosphate values are 0 throughout and cannot run full Mantel's
ST.microeco.mantel.select.noalk.all <- data.frame(ST.microeco.mantel.select.noalk$res_mantel)

## H2S values are 0 throughout, EC is same value at 0.63, and cannot run full Mantel's
SM.microeco.mantel.select.noalk.all <- data.frame(SM.microeco.mantel.select.noalk$res_mantel)

## H2S values are 0 throughout, pH is same value at 7.6, and cannot run full Mantel's
CT.microeco.mantel.select.noalk.all <- data.frame(CT.microeco.mantel.select.noalk$res_mantel)

## Phosphate values are 0.1, pH is 6.8 and EC is 1.9 throughout and cannot run full Mantel's
SG.microeco.mantel.select.noalk.all <- data.frame(SG.microeco.mantel.select.noalk$res_mantel)


#microeco.mantel.select.noalk.2 <- data.frame(microeco.mantel.select.noalk$res_mantel) %>% .[, c(1, 2, 5, 7)]

#test.all <-  data.frame(test$res_mantel)

# rename columns
#colnames(microeco.mantel.select.noalk.2) <- c("spec", "env", "r", "p.value")

# generate interval data
#microeco.mantel.select.noalk.2 %<>% dplyr::mutate(rd = cut(r, breaks = c(-Inf, 0.3, 0.6, Inf), labels = c("< 0.3", "0.3 - 0.6", ">= 0.6")), pd = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf), labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))

library(ggcor)
set_scale()


##### Mantel's with groups for all 6 regions ####

test <- NT.microeco.dataset$tax_table

## Photosynthetic 

NT.photosynthetic <- clone(NT.microeco.dataset)
NT.photosynthetic$tax_table <- NT.photosynthetic$tax_table[NT.photosynthetic$tax_table$PhotoHetero2 == "Photosynthetic", ]
#check.photo <- NT.photosynthetic$tax_table 
NT.photosynthetic$tidy_dataset()
NT.photosynthetic$cal_betadiv(unifrac = TRUE, binary = TRUE)
#photosynthetic$beta_diversity

NM.photosynthetic <- clone(NM.microeco.dataset)
NM.photosynthetic$tax_table <- NM.photosynthetic$tax_table[NM.photosynthetic$tax_table$PhotoHetero2 == "Photosynthetic", ]
check.photo <- NM.photosynthetic$tax_table 
#NM.photosynthetic$tidy_dataset()
NM.photosynthetic$cal_betadiv(unifrac = TRUE, binary = TRUE)

ST.photosynthetic <- clone(ST.microeco.dataset)
ST.photosynthetic$tax_table <- ST.photosynthetic$tax_table[ST.photosynthetic$tax_table$PhotoHetero2 == "Photosynthetic", ]
#check.photo <- ST.photosynthetic$tax_table 
ST.photosynthetic$tidy_dataset()
ST.photosynthetic$cal_betadiv(unifrac = TRUE, binary = TRUE)

SM.photosynthetic <- clone(SM.microeco.dataset)
SM.photosynthetic$tax_table <- SM.photosynthetic$tax_table[SM.photosynthetic$tax_table$PhotoHetero2 == "Photosynthetic", ]
#check.photo <- SM.photosynthetic$tax_table 
SM.photosynthetic$tidy_dataset()
SM.photosynthetic$cal_betadiv(unifrac = TRUE, binary = TRUE)

CT.photosynthetic <- clone(CT.microeco.dataset)
CT.photosynthetic$tax_table <- CT.photosynthetic$tax_table[CT.photosynthetic$tax_table$PhotoHetero2 == "Photosynthetic", ]
#check.photo <- CT.photosynthetic$tax_table 
CT.photosynthetic$tidy_dataset()
CT.photosynthetic$cal_betadiv(unifrac = TRUE, binary = TRUE)

SG.photosynthetic <- clone(SG.microeco.dataset)
SG.photosynthetic$tax_table <- SG.photosynthetic$tax_table[SG.photosynthetic$tax_table$PhotoHetero2 == "Photosynthetic", ]
#check.photo <- SG.photosynthetic$tax_table 
SG.photosynthetic$tidy_dataset()
SG.photosynthetic$cal_betadiv(unifrac = TRUE, binary = TRUE)


## Chemolithoautotrophs

NT.chemolithoautotroph <- clone(NT.microeco.dataset)
NT.chemolithoautotroph$tax_table <- NT.chemolithoautotroph$tax_table[NT.chemolithoautotroph$tax_table$PhotoHetero1 == "Chemolithoautotroph", ]
check.chemo <- NT.chemolithoautotroph$tax_table 
NT.chemolithoautotroph$tidy_dataset()
NT.chemolithoautotroph$cal_betadiv(unifrac = TRUE, binary = TRUE)

NM.chemolithoautotroph <- clone(NM.microeco.dataset)
NM.chemolithoautotroph$tax_table <- NM.chemolithoautotroph$tax_table[NM.chemolithoautotroph$tax_table$PhotoHetero1 == "Chemolithoautotroph", ]
check.chemo <- NM.chemolithoautotroph$tax_table 
NM.chemolithoautotroph$tidy_dataset()
NM.chemolithoautotroph$cal_betadiv(unifrac = TRUE, binary = TRUE)

ST.chemolithoautotroph <- clone(ST.microeco.dataset)
ST.chemolithoautotroph$tax_table <- ST.chemolithoautotroph$tax_table[ST.chemolithoautotroph$tax_table$PhotoHetero1 == "Chemolithoautotroph", ]
check.chemo <- ST.chemolithoautotroph$tax_table 
ST.chemolithoautotroph$tidy_dataset()
ST.chemolithoautotroph$cal_betadiv(unifrac = TRUE, binary = TRUE)

SM.chemolithoautotroph <- clone(SM.microeco.dataset)
SM.chemolithoautotroph$tax_table <- SM.chemolithoautotroph$tax_table[SM.chemolithoautotroph$tax_table$PhotoHetero1 == "Chemolithoautotroph", ]
check.chemo <- SM.chemolithoautotroph$tax_table 
SM.chemolithoautotroph$tidy_dataset()
SM.chemolithoautotroph$cal_betadiv(unifrac = TRUE, binary = TRUE)

CT.chemolithoautotroph <- clone(CT.microeco.dataset)
CT.chemolithoautotroph$tax_table <- CT.chemolithoautotroph$tax_table[CT.chemolithoautotroph$tax_table$PhotoHetero1 == "Chemolithoautotroph", ]
check.chemo <- CT.chemolithoautotroph$tax_table 
CT.chemolithoautotroph$tidy_dataset()
CT.chemolithoautotroph$cal_betadiv(unifrac = TRUE, binary = TRUE)

SG.chemolithoautotroph <- clone(SG.microeco.dataset)
SG.chemolithoautotroph$tax_table <- SG.chemolithoautotroph$tax_table[SG.chemolithoautotroph$tax_table$PhotoHetero1 == "Chemolithoautotroph", ]
check.chemo <- SG.chemolithoautotroph$tax_table 
SG.chemolithoautotroph$tidy_dataset()
SG.chemolithoautotroph$cal_betadiv(unifrac = TRUE, binary = TRUE)


## Heterotrophs

NT.heterotroph <- clone(NT.microeco.dataset)
NT.heterotroph$tax_table <- NT.heterotroph$tax_table[NT.heterotroph$tax_table$PhotoHetero1 == "Heterotrophs", ]
NT.check.hetero <- NT.heterotroph$tax_table
NT.heterotroph$tidy_dataset()
NT.heterotroph$cal_betadiv(unifrac = TRUE, binary = TRUE)

NM.heterotroph <- clone(NM.microeco.dataset)
NM.heterotroph$tax_table <- NM.heterotroph$tax_table[NM.heterotroph$tax_table$PhotoHetero1 == "Heterotrophs", ]
NM.check.hetero <- NM.heterotroph$tax_table
NM.heterotroph$tidy_dataset()
NM.heterotroph$cal_betadiv(unifrac = TRUE, binary = TRUE)

ST.heterotroph <- clone(ST.microeco.dataset)
ST.heterotroph$tax_table <- ST.heterotroph$tax_table[ST.heterotroph$tax_table$PhotoHetero1 == "Heterotrophs", ]
ST.check.hetero <- ST.heterotroph$tax_table
ST.heterotroph$tidy_dataset()
ST.heterotroph$cal_betadiv(unifrac = TRUE, binary = TRUE)

SM.heterotroph <- clone(SM.microeco.dataset)
SM.heterotroph$tax_table <- SM.heterotroph$tax_table[SM.heterotroph$tax_table$PhotoHetero1 == "Heterotrophs", ]
SM.check.hetero <- SM.heterotroph$tax_table
SM.heterotroph$tidy_dataset()
SM.heterotroph$cal_betadiv(unifrac = TRUE, binary = TRUE)

CT.heterotroph <- clone(CT.microeco.dataset)
CT.heterotroph$tax_table <- CT.heterotroph$tax_table[CT.heterotroph$tax_table$PhotoHetero1 == "Heterotrophs", ]
CT.check.hetero <- CT.heterotroph$tax_table
CT.heterotroph$tidy_dataset()
CT.heterotroph$cal_betadiv(unifrac = TRUE, binary = TRUE)

SG.heterotroph <- clone(SG.microeco.dataset)
SG.heterotroph$tax_table <- SG.heterotroph$tax_table[SG.heterotroph$tax_table$PhotoHetero1 == "Heterotrophs", ]
SG.check.hetero <- SG.heterotroph$tax_table
SG.heterotroph$tidy_dataset()
SG.heterotroph$cal_betadiv(unifrac = TRUE, binary = TRUE)

#("DSSF69","Elioraea","Methylobacterium-Methylorubrum","Rhodomicrobium","Roseomonas","Sandaracinobacter","Tabrizicola","Unassigned Rhodobacteraceae (Family)","Unassigned Sphingomonadaceae (Family)","AAP99","Allochromatium","Caldimonas","Curvibacter","DSSD61","Thiolamprovum","Unassigned B1-7BS (Family)","Unassigned Burkholderiales (Order)","Unassigned Comamonadaceae (Family)","Unassigned Gammaproteobacteria (Class)","Unassigned Rhodocyclaceae (Family)","Unassigned Sutterellaceae (Family)","Z-35")

unique(check.chemo$Genus)

# first perform mantel test on all groups 

NT.photosynthetic.mantel.select.noalk <- trans_env$new(dataset = NT.photosynthetic, env_cols = 23:28)
NM.photosynthetic.mantel.select.noalk <- trans_env$new(dataset = NM.photosynthetic, env_cols = 23:28)
ST.photosynthetic.mantel.select.noalk <- trans_env$new(dataset = ST.photosynthetic, env_cols = 23:28)
SM.photosynthetic.mantel.select.noalk <- trans_env$new(dataset = SM.photosynthetic, env_cols = 23:28)
CT.photosynthetic.mantel.select.noalk <- trans_env$new(dataset = CT.photosynthetic, env_cols = 23:28)
SG.photosynthetic.mantel.select.noalk <- trans_env$new(dataset = SG.photosynthetic, env_cols = 23:28)


NT.photosynthetic.mantel.select.noalk$cal_mantel(use_measure = "unwei_unifrac", method = "pearson" , partial_mantel = FALSE)
NM.photosynthetic.mantel.select.noalk$cal_mantel(use_measure = "unwei_unifrac", method = "pearson" , partial_mantel = FALSE)
ST.photosynthetic.mantel.select.noalk$cal_mantel(use_measure = "unwei_unifrac", method = "pearson" , partial_mantel = FALSE)
SM.photosynthetic.mantel.select.noalk$cal_mantel(use_measure = "unwei_unifrac", method = "pearson" , partial_mantel = FALSE)
CT.photosynthetic.mantel.select.noalk$cal_mantel(use_measure = "unwei_unifrac", method = "pearson" , partial_mantel = FALSE)
SG.photosynthetic.mantel.select.noalk$cal_mantel(use_measure = "unwei_unifrac", method = "pearson" , partial_mantel = FALSE)


NT.chemolithoautotroph.mantel.select.noalk <- trans_env$new(dataset = NT.chemolithoautotroph, env_cols = 23:28)
NM.chemolithoautotroph.mantel.select.noalk <- trans_env$new(dataset = NM.chemolithoautotroph, env_cols = 23:28)
ST.chemolithoautotroph.mantel.select.noalk <- trans_env$new(dataset = ST.chemolithoautotroph, env_cols = 23:28)
SM.chemolithoautotroph.mantel.select.noalk <- trans_env$new(dataset = SM.chemolithoautotroph, env_cols = 23:28)
CT.chemolithoautotroph.mantel.select.noalk <- trans_env$new(dataset = CT.chemolithoautotroph, env_cols = 23:28)
SG.chemolithoautotroph.mantel.select.noalk <- trans_env$new(dataset = SG.chemolithoautotroph, env_cols = 23:28)


NT.chemolithoautotroph.mantel.select.noalk$cal_mantel(use_measure = "unwei_unifrac", method = "pearson" , partial_mantel = FALSE)
NM.chemolithoautotroph.mantel.select.noalk$cal_mantel(use_measure = "unwei_unifrac", method = "pearson" , partial_mantel = FALSE)
ST.chemolithoautotroph.mantel.select.noalk$cal_mantel(use_measure = "unwei_unifrac", method = "pearson" , partial_mantel = FALSE)
SM.chemolithoautotroph.mantel.select.noalk$cal_mantel(use_measure = "unwei_unifrac", method = "pearson" , partial_mantel = FALSE)
CT.chemolithoautotroph.mantel.select.noalk$cal_mantel(use_measure = "unwei_unifrac", method = "pearson" , partial_mantel = FALSE)
SG.chemolithoautotroph.mantel.select.noalk$cal_mantel(use_measure = "unwei_unifrac", method = "pearson" , partial_mantel = FALSE)


NT.heterotroph.mantel.select.noalk <- trans_env$new(dataset = NT.heterotroph, env_cols = 23:28)
NM.heterotroph.mantel.select.noalk <- trans_env$new(dataset = NM.heterotroph, env_cols = 23:28)
ST.heterotroph.mantel.select.noalk <- trans_env$new(dataset = ST.heterotroph, env_cols = 23:28)
SM.heterotroph.mantel.select.noalk <- trans_env$new(dataset = SM.heterotroph, env_cols = 23:28)
CT.heterotroph.mantel.select.noalk <- trans_env$new(dataset = CT.heterotroph, env_cols = 23:28)
SG.heterotroph.mantel.select.noalk <- trans_env$new(dataset = SG.heterotroph, env_cols = 23:28)


NT.heterotroph.mantel.select.noalk$cal_mantel(use_measure = "unwei_unifrac", method = "pearson" , partial_mantel = FALSE)
NM.heterotroph.mantel.select.noalk$cal_mantel(use_measure = "unwei_unifrac", method = "pearson" , partial_mantel = FALSE)
ST.heterotroph.mantel.select.noalk$cal_mantel(use_measure = "unwei_unifrac", method = "pearson" , partial_mantel = FALSE)
SM.heterotroph.mantel.select.noalk$cal_mantel(use_measure = "unwei_unifrac", method = "pearson" , partial_mantel = FALSE)
CT.heterotroph.mantel.select.noalk$cal_mantel(use_measure = "unwei_unifrac", method = "pearson" , partial_mantel = FALSE)
SG.heterotroph.mantel.select.noalk$cal_mantel(use_measure = "unwei_unifrac", method = "pearson" , partial_mantel = FALSE)


### use microeco.mantel.all to subset for required columns

#NT.photosynthetic.mantel.select.noalk.all <- data.frame(NT.photosynthetic.mantel.select.noalk$res_mantel)
#photosynthetic.mantel.select.noalk.2 <- data.frame(photosynthetic.mantel.select.noalk$res_mantel) %>% .[, c(1, 2, 5, 7)]

#NT.chemolithoautotroph.mantel.select.noalk.all <- data.frame(NT.chemolithoautotroph.mantel.select.noalk$res_mantel)
#chemolithoautotroph.mantel.select.noalk.2 <- data.frame(chemolithoautotroph.mantel.select.noalk$res_mantel) %>% .[, c(1, 2, 5, 7)]

#NT.heterotroph.mantel.select.noalk.all <- data.frame(NT.heterotroph.mantel.select.noalk$res_mantel)
#heterotroph.mantel.select.noalk.2 <- data.frame(heterotroph.mantel.select.noalk$res_mantel) %>% .[, c(1, 2, 5, 7)]


library(ggcor)
set_scale()


### Combine all 3 groups together - for each region ####

## Selected noalk full mantel's test

#test.all <- data.frame(spec = "All", microeco.mantel.select.noalk$res_mantel) %>% .[, c(1, 3, 6, 8)]

#All.select.noalk <- microeco.mantel.select.noalk.2

## North Thailand - all main groups ####

NT.All.select.noalk <- data.frame(spec = "All", NT.microeco.mantel.select.noalk$res_mantel) %>% .[, c(1, 3, 6, 8)]

NT.Photosynthetic.select.noalk <- data.frame(spec = "Photosynthetic", NT.photosynthetic.mantel.select.noalk$res_mantel) %>% .[, c(1, 3, 6, 8)]

NT.Chemolithoautotroph.select.noalk <- data.frame(spec = "Chemolithoautotroph", NT.chemolithoautotroph.mantel.select.noalk$res_mantel) %>% .[, c(1, 3, 6, 8)]

NT.Heterotroph.select.noalk <- data.frame(spec = "Heterotrophs", NT.heterotroph.mantel.select.noalk$res_mantel) %>% .[, c(1, 3, 6, 8)]

# rename columns
colnames(NT.All.select.noalk) <- colnames(NT.Photosynthetic.select.noalk) <- colnames(NT.Chemolithoautotroph.select.noalk) <- colnames(NT.Heterotroph.select.noalk) <- c("spec", "env", "r", "p.value")

# generate interval data

NT.All.select.noalk %<>% dplyr::mutate(rd = cut(r, breaks = c(-Inf, 0.3, 0.6, Inf), labels = c("< 0.3", "0.3 - 0.6", ">= 0.6")), pd = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf), labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))

NT.Photosynthetic.select.noalk %<>% dplyr::mutate(rd = cut(r, breaks = c(-Inf, 0.3, 0.6, Inf), labels = c("< 0.3", "0.3 - 0.6", ">= 0.6")),
                                                  pd = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf), labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))

NT.Chemolithoautotroph.select.noalk %<>% dplyr::mutate(rd = cut(r, breaks = c(-Inf, 0.3, 0.6, Inf), labels = c("< 0.3", "0.3 - 0.6", ">= 0.6")),
                                                       pd = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf), labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))

NT.Heterotroph.select.noalk %<>% dplyr::mutate(rd = cut(r, breaks = c(-Inf, 0.3, 0.6, Inf), labels = c("< 0.3", "0.3 - 0.6", ">= 0.6")),
                                               pd = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf), labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))

# cobine 3 tables
plot_table_select_noalk_NT <- rbind(NT.All.select.noalk, NT.Photosynthetic.select.noalk , NT.Chemolithoautotroph.select.noalk, NT.Heterotroph.select.noalk )

library(ggcor)
set_scale()

NT.main.groups.select.noalk <- quickcor(NT.microeco.mantel.select.noalk$data_env, type = "upper", cor.test = TRUE, show.diag = FALSE) +
  geom_square() +
  geom_mark(sig.thres = 0.05, markonly = TRUE, color = "white") +
  anno_link(aes(colour = pd, size = rd), data = plot_table_select_noalk_NT) +
  scale_size_manual(values = c(0.5, 1.5, 3)) +
  scale_colour_manual(values = c("#D95F02", "#1B9E77", "#A2A2A288")) +
  guides(size = guide_legend(title = "Mantel's r", override.aes = list(colour = "grey35"), order = 2),
         colour = guide_legend(title = "Mantel's p", override.aes = list(size = 3), order = 1),
         fill = guide_colorbar(title = "Pearson's r", order = 3)) 

pdf("NT.main.groups.mantel.select.noalk.pdf", width = 10)
NT.main.groups.select.noalk
dev.off()

ggsave("NT.main.groups.select.noalk.png", height = 10, width =10)


## North Malaysia - all main groups ####

NM.All.select.noalk <- data.frame(spec = "All", NM.microeco.mantel.select.noalk$res_mantel) %>% .[, c(1, 3, 6, 8)]

NM.Photosynthetic.select.noalk <- data.frame(spec = "Photosynthetic", NM.photosynthetic.mantel.select.noalk$res_mantel) %>% .[, c(1, 3, 6, 8)]

NM.Chemolithoautotroph.select.noalk <- data.frame(spec = "Chemolithoautotroph", NM.chemolithoautotroph.mantel.select.noalk$res_mantel) %>% .[, c(1, 3, 6, 8)]

NM.Heterotroph.select.noalk <- data.frame(spec = "Heterotrophs", NM.heterotroph.mantel.select.noalk$res_mantel) %>% .[, c(1, 3, 6, 8)]

# rename columns
colnames(NM.All.select.noalk) <- colnames(NM.Photosynthetic.select.noalk) <- colnames(NM.Chemolithoautotroph.select.noalk) <- colnames(NM.Heterotroph.select.noalk) <- c("spec", "env", "r", "p.value")

# generate interval data

NM.All.select.noalk %<>% dplyr::mutate(rd = cut(r, breaks = c(-Inf, 0.3, 0.6, Inf), labels = c("< 0.3", "0.3 - 0.6", ">= 0.6")), pd = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf), labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))

NM.Photosynthetic.select.noalk %<>% dplyr::mutate(rd = cut(r, breaks = c(-Inf, 0.3, 0.6, Inf), labels = c("< 0.3", "0.3 - 0.6", ">= 0.6")),
                                                  pd = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf), labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))

NM.Chemolithoautotroph.select.noalk %<>% dplyr::mutate(rd = cut(r, breaks = c(-Inf, 0.3, 0.6, Inf), labels = c("< 0.3", "0.3 - 0.6", ">= 0.6")),
                                                       pd = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf), labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))

NM.Heterotroph.select.noalk %<>% dplyr::mutate(rd = cut(r, breaks = c(-Inf, 0.3, 0.6, Inf), labels = c("< 0.3", "0.3 - 0.6", ">= 0.6")),
                                               pd = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf), labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))

# cobine 3 tables
plot_table_select_noalk_NM <- rbind(NM.All.select.noalk, NM.Photosynthetic.select.noalk , NM.Chemolithoautotroph.select.noalk, NM.Heterotroph.select.noalk )

library(ggcor)
set_scale()

NM.main.groups.select.noalk <- quickcor(NM.microeco.mantel.select.noalk$data_env, type = "upper", cor.test = TRUE, show.diag = FALSE) +
  geom_square() +
  geom_mark(sig.thres = 0.05, markonly = TRUE, color = "white") +
  anno_link(aes(colour = pd, size = rd), data = plot_table_select_noalk_NM) +
  scale_size_manual(values = c(0.5, 1.5, 3)) +
  scale_colour_manual(values = c("#D95F02", "#1B9E77", "#A2A2A288")) +
  guides(size = guide_legend(title = "Mantel's r", override.aes = list(colour = "grey35"), order = 2),
         colour = guide_legend(title = "Mantel's p", override.aes = list(size = 3), order = 1),
         fill = guide_colorbar(title = "Pearson's r", order = 3)) 

pdf("NM.main.groups.mantel.select.noalk.pdf", width = 10)
NM.main.groups.select.noalk
dev.off()

ggsave("NM.main.groups.select.noalk.png", height = 10, width =10)


## South Thailand - all main groups ####

ST.All.select.noalk <- data.frame(spec = "All", ST.microeco.mantel.select.noalk$res_mantel) %>% .[, c(1, 3, 6, 8)]

ST.Photosynthetic.select.noalk <- data.frame(spec = "Photosynthetic", ST.photosynthetic.mantel.select.noalk$res_mantel) %>% .[, c(1, 3, 6, 8)]

ST.Chemolithoautotroph.select.noalk <- data.frame(spec = "Chemolithoautotroph", ST.chemolithoautotroph.mantel.select.noalk$res_mantel) %>% .[, c(1, 3, 6, 8)]

ST.Heterotroph.select.noalk <- data.frame(spec = "Heterotrophs", ST.heterotroph.mantel.select.noalk$res_mantel) %>% .[, c(1, 3, 6, 8)]

# rename columns
colnames(ST.All.select.noalk) <- colnames(ST.Photosynthetic.select.noalk) <- colnames(ST.Chemolithoautotroph.select.noalk) <- colnames(ST.Heterotroph.select.noalk) <- c("spec", "env", "r", "p.value")

# generate interval data

ST.All.select.noalk %<>% dplyr::mutate(rd = cut(r, breaks = c(-Inf, 0.3, 0.6, Inf), labels = c("< 0.3", "0.3 - 0.6", ">= 0.6")), pd = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf), labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))

ST.Photosynthetic.select.noalk %<>% dplyr::mutate(rd = cut(r, breaks = c(-Inf, 0.3, 0.6, Inf), labels = c("< 0.3", "0.3 - 0.6", ">= 0.6")),
                                                  pd = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf), labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))

ST.Chemolithoautotroph.select.noalk %<>% dplyr::mutate(rd = cut(r, breaks = c(-Inf, 0.3, 0.6, Inf), labels = c("< 0.3", "0.3 - 0.6", ">= 0.6")),
                                                       pd = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf), labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))

ST.Heterotroph.select.noalk %<>% dplyr::mutate(rd = cut(r, breaks = c(-Inf, 0.3, 0.6, Inf), labels = c("< 0.3", "0.3 - 0.6", ">= 0.6")),
                                               pd = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf), labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))

# cobine 3 tables
plot_table_select_noalk_ST <- rbind(ST.All.select.noalk, ST.Photosynthetic.select.noalk , ST.Chemolithoautotroph.select.noalk, ST.Heterotroph.select.noalk )

library(ggcor)
set_scale()

ST.main.groups.select.noalk <- quickcor(ST.microeco.mantel.select.noalk$data_env, type = "upper", cor.test = TRUE, show.diag = FALSE) +
  geom_square() +
  geom_mark(sig.thres = 0.05, markonly = TRUE, color = "white") +
  anno_link(aes(colour = pd, size = rd), data = plot_table_select_noalk_ST) +
  scale_size_manual(values = c(0.5, 1.5, 3)) +
  scale_colour_manual(values = c("#D95F02", "#1B9E77", "#A2A2A288")) +
  guides(size = guide_legend(title = "Mantel's r", override.aes = list(colour = "grey35"), order = 2),
         colour = guide_legend(title = "Mantel's p", override.aes = list(size = 3), order = 1),
         fill = guide_colorbar(title = "Pearson's r", order = 3)) 

pdf("ST.main.groups.mantel.select.noalk.pdf", width = 10)
ST.main.groups.select.noalk
dev.off()

ggsave("ST.main.groups.select.noalk.png", height = 10, width =10)


## South Malaysia - all main groups ####

SM.All.select.noalk <- data.frame(spec = "All", SM.microeco.mantel.select.noalk$res_mantel) %>% .[, c(1, 3, 6, 8)]

SM.Photosynthetic.select.noalk <- data.frame(spec = "Photosynthetic", SM.photosynthetic.mantel.select.noalk$res_mantel) %>% .[, c(1, 3, 6, 8)]

SM.Chemolithoautotroph.select.noalk <- data.frame(spec = "Chemolithoautotroph", SM.chemolithoautotroph.mantel.select.noalk$res_mantel) %>% .[, c(1, 3, 6, 8)]

SM.Heterotroph.select.noalk <- data.frame(spec = "Heterotrophs", SM.heterotroph.mantel.select.noalk$res_mantel) %>% .[, c(1, 3, 6, 8)]

# rename columns
colnames(SM.All.select.noalk) <- colnames(SM.Photosynthetic.select.noalk) <- colnames(SM.Chemolithoautotroph.select.noalk) <- colnames(SM.Heterotroph.select.noalk) <- c("spec", "env", "r", "p.value")

# generate interval data

SM.All.select.noalk %<>% dplyr::mutate(rd = cut(r, breaks = c(-Inf, 0.3, 0.6, Inf), labels = c("< 0.3", "0.3 - 0.6", ">= 0.6")), pd = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf), labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))

SM.Photosynthetic.select.noalk %<>% dplyr::mutate(rd = cut(r, breaks = c(-Inf, 0.3, 0.6, Inf), labels = c("< 0.3", "0.3 - 0.6", ">= 0.6")),
                                                  pd = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf), labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))

SM.Chemolithoautotroph.select.noalk %<>% dplyr::mutate(rd = cut(r, breaks = c(-Inf, 0.3, 0.6, Inf), labels = c("< 0.3", "0.3 - 0.6", ">= 0.6")),
                                                       pd = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf), labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))

SM.Heterotroph.select.noalk %<>% dplyr::mutate(rd = cut(r, breaks = c(-Inf, 0.3, 0.6, Inf), labels = c("< 0.3", "0.3 - 0.6", ">= 0.6")),
                                               pd = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf), labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))

# cobine 3 tables
plot_table_select_noalk_SM <- rbind(SM.All.select.noalk, SM.Photosynthetic.select.noalk , SM.Chemolithoautotroph.select.noalk, SM.Heterotroph.select.noalk )

library(ggcor)
set_scale()

SM.main.groups.select.noalk <- quickcor(SM.microeco.mantel.select.noalk$data_env, type = "upper", cor.test = TRUE, show.diag = FALSE) +
  geom_square() +
  geom_mark(sig.thres = 0.05, markonly = TRUE, color = "white") +
  anno_link(aes(colour = pd, size = rd), data = plot_table_select_noalk_SM) +
  scale_size_manual(values = c(0.5, 1.5, 3)) +
  scale_colour_manual(values = c("#D95F02", "#1B9E77", "#A2A2A288")) +
  guides(size = guide_legend(title = "Mantel's r", override.aes = list(colour = "grey35"), order = 2),
         colour = guide_legend(title = "Mantel's p", override.aes = list(size = 3), order = 1),
         fill = guide_colorbar(title = "Pearson's r", order = 3)) 

pdf("SM.main.groups.mantel.select.noalk.pdf", width = 10)
SM.main.groups.select.noalk
dev.off()

ggsave("SM.main.groups.select.noalk.png", height = 10, width =10)


## Central Thailand - all main groups ####

CT.All.select.noalk <- data.frame(spec = "All", CT.microeco.mantel.select.noalk$res_mantel) %>% .[, c(1, 3, 6, 8)]

CT.Photosynthetic.select.noalk <- data.frame(spec = "Photosynthetic", CT.photosynthetic.mantel.select.noalk$res_mantel) %>% .[, c(1, 3, 6, 8)]

CT.Chemolithoautotroph.select.noalk <- data.frame(spec = "Chemolithoautotroph", CT.chemolithoautotroph.mantel.select.noalk$res_mantel) %>% .[, c(1, 3, 6, 8)]

CT.Heterotroph.select.noalk <- data.frame(spec = "Heterotrophs", CT.heterotroph.mantel.select.noalk$res_mantel) %>% .[, c(1, 3, 6, 8)]

# rename columns
colnames(CT.All.select.noalk) <- colnames(CT.Photosynthetic.select.noalk) <- colnames(CT.Chemolithoautotroph.select.noalk) <- colnames(CT.Heterotroph.select.noalk) <- c("spec", "env", "r", "p.value")

# generate interval data

CT.All.select.noalk %<>% dplyr::mutate(rd = cut(r, breaks = c(-Inf, 0.3, 0.6, Inf), labels = c("< 0.3", "0.3 - 0.6", ">= 0.6")), pd = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf), labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))

CT.Photosynthetic.select.noalk %<>% dplyr::mutate(rd = cut(r, breaks = c(-Inf, 0.3, 0.6, Inf), labels = c("< 0.3", "0.3 - 0.6", ">= 0.6")),
                                                  pd = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf), labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))

CT.Chemolithoautotroph.select.noalk %<>% dplyr::mutate(rd = cut(r, breaks = c(-Inf, 0.3, 0.6, Inf), labels = c("< 0.3", "0.3 - 0.6", ">= 0.6")),
                                                       pd = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf), labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))

CT.Heterotroph.select.noalk %<>% dplyr::mutate(rd = cut(r, breaks = c(-Inf, 0.3, 0.6, Inf), labels = c("< 0.3", "0.3 - 0.6", ">= 0.6")),
                                               pd = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf), labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))

# cobine 3 tables
plot_table_select_noalk_CT <- rbind(CT.All.select.noalk, CT.Photosynthetic.select.noalk , CT.Chemolithoautotroph.select.noalk, CT.Heterotroph.select.noalk )

library(ggcor)
set_scale()

CT.main.groups.select.noalk <- quickcor(CT.microeco.mantel.select.noalk$data_env, type = "upper", cor.test = TRUE, show.diag = FALSE) +
  geom_square() +
  geom_mark(sig.thres = 0.05, markonly = TRUE, color = "white") +
  anno_link(aes(colour = pd, size = rd), data = plot_table_select_noalk_CT) +
  scale_size_manual(values = c(0.5, 1.5, 3)) +
  scale_colour_manual(values = c("#D95F02", "#1B9E77", "#A2A2A288")) +
  guides(size = guide_legend(title = "Mantel's r", override.aes = list(colour = "grey35"), order = 2),
         colour = guide_legend(title = "Mantel's p", override.aes = list(size = 3), order = 1),
         fill = guide_colorbar(title = "Pearson's r", order = 3)) 

pdf("CT.main.groups.mantel.select.noalk.pdf", width = 10)
CT.main.groups.select.noalk
dev.off()

ggsave("CT.main.groups.select.noalk.png", height = 10, width =10)


## Singapore - all main groups ####

SG.All.select.noalk <- data.frame(spec = "All", SG.microeco.mantel.select.noalk$res_mantel) %>% .[, c(1, 3, 6, 8)]

SG.Photosynthetic.select.noalk <- data.frame(spec = "Photosynthetic", SG.photosynthetic.mantel.select.noalk$res_mantel) %>% .[, c(1, 3, 6, 8)]

SG.Chemolithoautotroph.select.noalk <- data.frame(spec = "Chemolithoautotroph", SG.chemolithoautotroph.mantel.select.noalk$res_mantel) %>% .[, c(1, 3, 6, 8)]

SG.Heterotroph.select.noalk <- data.frame(spec = "Heterotrophs", SG.heterotroph.mantel.select.noalk$res_mantel) %>% .[, c(1, 3, 6, 8)]

# rename columns
colnames(SG.All.select.noalk) <- colnames(SG.Photosynthetic.select.noalk) <- colnames(SG.Chemolithoautotroph.select.noalk) <- colnames(SG.Heterotroph.select.noalk) <- c("spec", "env", "r", "p.value")

# generate interval data

SG.All.select.noalk %<>% dplyr::mutate(rd = cut(r, breaks = c(-Inf, 0.3, 0.6, Inf), labels = c("< 0.3", "0.3 - 0.6", ">= 0.6")), pd = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf), labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))

SG.Photosynthetic.select.noalk %<>% dplyr::mutate(rd = cut(r, breaks = c(-Inf, 0.3, 0.6, Inf), labels = c("< 0.3", "0.3 - 0.6", ">= 0.6")),
                                                  pd = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf), labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))

SG.Chemolithoautotroph.select.noalk %<>% dplyr::mutate(rd = cut(r, breaks = c(-Inf, 0.3, 0.6, Inf), labels = c("< 0.3", "0.3 - 0.6", ">= 0.6")),
                                                       pd = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf), labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))

SG.Heterotroph.select.noalk %<>% dplyr::mutate(rd = cut(r, breaks = c(-Inf, 0.3, 0.6, Inf), labels = c("< 0.3", "0.3 - 0.6", ">= 0.6")),
                                               pd = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf), labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))

# cobine 3 tables
plot_table_select_noalk_SG <- rbind(SG.All.select.noalk, SG.Photosynthetic.select.noalk , SG.Chemolithoautotroph.select.noalk, SG.Heterotroph.select.noalk )

library(ggcor)
set_scale()

SG.main.groups.select.noalk <- quickcor(SG.microeco.mantel.select.noalk$data_env, type = "upper", cor.test = TRUE, show.diag = FALSE) +
  geom_square() +
  geom_mark(sig.thres = 0.05, markonly = TRUE, color = "white") +
  anno_link(aes(colour = pd, size = rd), data = plot_table_select_noalk_SG) +
  scale_size_manual(values = c(0.5, 1.5, 3)) +
  scale_colour_manual(values = c("#D95F02", "#1B9E77", "#A2A2A288")) +
  guides(size = guide_legend(title = "Mantel's r", override.aes = list(colour = "grey35"), order = 2),
         colour = guide_legend(title = "Mantel's p", override.aes = list(size = 3), order = 1),
         fill = guide_colorbar(title = "Pearson's r", order = 3)) 

pdf("SG.main.groups.mantel.select.noalk.pdf", width = 10)
SG.main.groups.select.noalk
dev.off()

ggsave("SG.main.groups.select.noalk.png", height = 10, width =10)



########################################################################

### all groups across overall and 6 regions #####

#NOTREQ
All.select.noalk <- data.frame(spec = "All", microeco.mantel.select.noalk$res_mantel) %>% .[, c(1, 3, 6, 8)]

NT.All.select.noalk <- data.frame(spec = "NT", NT.microeco.mantel.select.noalk$res_mantel) %>% .[, c(1, 3, 6, 8)]

NM.All.select.noalk <- data.frame(spec = "NM", NM.microeco.mantel.select.noalk$res_mantel) %>% .[, c(1, 3, 6, 8)]

ST.All.select.noalk <- data.frame(spec = "ST", ST.microeco.mantel.select.noalk$res_mantel) %>% .[, c(1, 3, 6, 8)]

SM.All.select.noalk <- data.frame(spec = "SM", SM.microeco.mantel.select.noalk$res_mantel) %>% .[, c(1, 3, 6, 8)]

CT.All.select.noalk <- data.frame(spec = "CT", CT.microeco.mantel.select.noalk$res_mantel) %>% .[, c(1, 3, 6, 8)]

SG.All.select.noalk <- data.frame(spec = "SG", SG.microeco.mantel.select.noalk$res_mantel) %>% .[, c(1, 3, 6, 8)]

#NOTREQ
# rename columns
colnames(All.select.noalk) <- colnames(NT.All.select.noalk) <- colnames(NM.All.select.noalk) <- colnames(ST.All.select.noalk) <- colnames(SM.All.select.noalk) <- colnames(CT.All.select.noalk) <- colnames(SG.All.select.noalk) <- c("spec", "env", "r", "p.value")

#NOTREQ
# generate interval data

All.select.noalk %<>% dplyr::mutate(rd = cut(r, breaks = c(-Inf, 0.3, 0.6, Inf), labels = c("< 0.3", "0.3 - 0.6", ">= 0.6")), pd = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf), labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))

NT.All.select.noalk %<>% dplyr::mutate(rd = cut(r, breaks = c(-Inf, 0.3, 0.6, Inf), labels = c("< 0.3", "0.3 - 0.6", ">= 0.6")), pd = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf), labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))

NM.All.select.noalk %<>% dplyr::mutate(rd = cut(r, breaks = c(-Inf, 0.3, 0.6, Inf), labels = c("< 0.3", "0.3 - 0.6", ">= 0.6")), pd = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf), labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))

ST.All.select.noalk %<>% dplyr::mutate(rd = cut(r, breaks = c(-Inf, 0.3, 0.6, Inf), labels = c("< 0.3", "0.3 - 0.6", ">= 0.6")), pd = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf), labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))

SM.All.select.noalk %<>% dplyr::mutate(rd = cut(r, breaks = c(-Inf, 0.3, 0.6, Inf), labels = c("< 0.3", "0.3 - 0.6", ">= 0.6")), pd = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf), labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))

CT.All.select.noalk %<>% dplyr::mutate(rd = cut(r, breaks = c(-Inf, 0.3, 0.6, Inf), labels = c("< 0.3", "0.3 - 0.6", ">= 0.6")), pd = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf), labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))

SG.All.select.noalk %<>% dplyr::mutate(rd = cut(r, breaks = c(-Inf, 0.3, 0.6, Inf), labels = c("< 0.3", "0.3 - 0.6", ">= 0.6")), pd = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf), labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))


# cobine 3 tables
plot_table_select_noalk_all_reg <- rbind(All.select.noalk, NT.All.select.noalk, ST.All.select.noalk, NM.All.select.noalk, SM.All.select.noalk, SG.All.select.noalk, CT.All.select.noalk)

saveRDS(plot_table_select_noalk_all_reg, "plot_table_select_noalk_all_reg.rds")

library(ggcor)
set_scale()


noalk.all.reg <- quickcor(microeco.mantel.select.noalk$data_env, type = "upper", cor.test = TRUE, show.diag = FALSE) +
  geom_square() +
  geom_mark(sig.thres = 0.05, markonly = TRUE, color = "white") +
  anno_link(aes(colour = pd, size = rd), data = plot_table_select_noalk_all_reg) +
  scale_size_manual(values = c(0.5, 1.5, 3)) +
  scale_colour_manual(values = c("#D95F02", "#1B9E77", "#A2A2A288")) +
  guides(size = guide_legend(title = "Mantel's r", override.aes = list(colour = "grey35"), order = 2),
         colour = guide_legend(title = "Mantel's p", override.aes = list(size = 3), order = 1),
         fill = guide_colorbar(title = "Pearson's r", order = 3)) 

pdf("noalk.all.reg.pdf", width = 10)
noalk.all.reg
dev.off()

ggsave("noalk.all.reg.png", height = 10, width =10)

save.image("TMFZ_pseudopool_workspace_5.RData")


```





```{r}

#### Mantels test by locations ####

location <- unique(metadata.micronet$Location.Code)

# [1] "AP" "BA" "HS" "KJ" "LA" "LN" "MK" "PB" "PP" "PT" "RB" "RN" "SE" "SW" "US"


## create dataset by locations

## AP ####

AP.microeco.dataset <- clone(microeco.dataset.3)

AP.microeco.dataset$sample_table <- subset(AP.microeco.dataset$sample_table, Location.Code == "AP")

AP.microeco.dataset$tidy_dataset()

a <- as.data.frame(AP.microeco.dataset$otu_table)
b <- as.data.frame(AP.microeco.dataset$sample_table)


## BA ####

BA.microeco.dataset <- clone(microeco.dataset.3)

BA.microeco.dataset$sample_table <- subset(BA.microeco.dataset$sample_table, Location.Code == "BA")

BA.microeco.dataset$tidy_dataset()

a <- as.data.frame(BA.microeco.dataset$otu_table)
b <- as.data.frame(BA.microeco.dataset$sample_table)


## HS ####

HS.microeco.dataset <- clone(microeco.dataset.3)

HS.microeco.dataset$sample_table <- subset(HS.microeco.dataset$sample_table, Location.Code == "HS")

HS.microeco.dataset$tidy_dataset()

a <- as.data.frame(HS.microeco.dataset$otu_table)
b <- as.data.frame(HS.microeco.dataset$sample_table)


## KJ ####

KJ.microeco.dataset <- clone(microeco.dataset.3)

KJ.microeco.dataset$sample_table <- subset(KJ.microeco.dataset$sample_table, Location.Code == "KJ")

KJ.microeco.dataset$tidy_dataset()

a <- as.data.frame(KJ.microeco.dataset$otu_table)
b <- as.data.frame(KJ.microeco.dataset$sample_table)


## LA ####

LA.microeco.dataset <- clone(microeco.dataset.3)

LA.microeco.dataset$sample_table <- subset(LA.microeco.dataset$sample_table, Location.Code == "LA")

LA.microeco.dataset$tidy_dataset()

a <- as.data.frame(LA.microeco.dataset$otu_table)
b <- as.data.frame(LA.microeco.dataset$sample_table)


## LN ####

LN.microeco.dataset <- clone(microeco.dataset.3)

LN.microeco.dataset$sample_table <- subset(LN.microeco.dataset$sample_table, Location.Code == "LN")

LN.microeco.dataset$tidy_dataset()

a <- as.data.frame(LN.microeco.dataset$otu_table)
b <- as.data.frame(LN.microeco.dataset$sample_table)


## MK ####

MK.microeco.dataset <- clone(microeco.dataset.3)

MK.microeco.dataset$sample_table <- subset(MK.microeco.dataset$sample_table, Location.Code == "MK")

MK.microeco.dataset$tidy_dataset()

a <- as.data.frame(MK.microeco.dataset$otu_table)
b <- as.data.frame(MK.microeco.dataset$sample_table)


## PB ####

PB.microeco.dataset <- clone(microeco.dataset.3)

PB.microeco.dataset$sample_table <- subset(PB.microeco.dataset$sample_table, Location.Code == "PB")

PB.microeco.dataset$tidy_dataset()

a <- as.data.frame(PB.microeco.dataset$otu_table)
b <- as.data.frame(PB.microeco.dataset$sample_table)


## PP ####

PP.microeco.dataset <- clone(microeco.dataset.3)

PP.microeco.dataset$sample_table <- subset(PP.microeco.dataset$sample_table, Location.Code == "PP")

PP.microeco.dataset$tidy_dataset()

a <- as.data.frame(PP.microeco.dataset$otu_table)
b <- as.data.frame(PP.microeco.dataset$sample_table)


## PT ####

PT.microeco.dataset <- clone(microeco.dataset.3)

PT.microeco.dataset$sample_table <- subset(PT.microeco.dataset$sample_table, Location.Code == "PT")

PT.microeco.dataset$tidy_dataset()

a <- as.data.frame(PT.microeco.dataset$otu_table)
b <- as.data.frame(PT.microeco.dataset$sample_table)


## RB ####

RB.microeco.dataset <- clone(microeco.dataset.3)

RB.microeco.dataset$sample_table <- subset(RB.microeco.dataset$sample_table, Location.Code == "RB")

RB.microeco.dataset$tidy_dataset()

a <- as.data.frame(RB.microeco.dataset$otu_table)
b <- as.data.frame(RB.microeco.dataset$sample_table)



## RN ####

RN.microeco.dataset <- clone(microeco.dataset.3)

RN.microeco.dataset$sample_table <- subset(RN.microeco.dataset$sample_table, Location.Code == "RN")

RN.microeco.dataset$tidy_dataset()

a <- as.data.frame(RN.microeco.dataset$otu_table)
b <- as.data.frame(RN.microeco.dataset$sample_table)


## SE ####

SE.microeco.dataset <- clone(microeco.dataset.3)

SE.microeco.dataset$sample_table <- subset(SE.microeco.dataset$sample_table, Location.Code == "SE")

SE.microeco.dataset$tidy_dataset()

a <- as.data.frame(SE.microeco.dataset$otu_table)
b <- as.data.frame(SE.microeco.dataset$sample_table)


## SW ####

SW.microeco.dataset <- clone(microeco.dataset.3)

SW.microeco.dataset$sample_table <- subset(SW.microeco.dataset$sample_table, Location.Code == "SW")

SW.microeco.dataset$tidy_dataset()

a <- as.data.frame(SW.microeco.dataset$otu_table)
b <- as.data.frame(SW.microeco.dataset$sample_table)


## US ####

US.microeco.dataset <- clone(microeco.dataset.3)

US.microeco.dataset$sample_table <- subset(US.microeco.dataset$sample_table, Location.Code == "US")

US.microeco.dataset$tidy_dataset()

a <- as.data.frame(US.microeco.dataset$otu_table)
b <- as.data.frame(US.microeco.dataset$sample_table)



#unifrac default FALSE; whether UniFrac index should be calculated,binary default FALSE; TRUE is used for jaccard and unweighted unifrac;
#microeco.dataset.2$cal_betadiv(unifrac = TRUE, binary = TRUE)
AP.microeco.dataset$cal_betadiv(unifrac = TRUE, binary = TRUE)
BA.microeco.dataset$cal_betadiv(unifrac = TRUE, binary = TRUE)
HS.microeco.dataset$cal_betadiv(unifrac = TRUE, binary = TRUE)
KJ.microeco.dataset$cal_betadiv(unifrac = TRUE, binary = TRUE)
LA.microeco.dataset$cal_betadiv(unifrac = TRUE, binary = TRUE)
LN.microeco.dataset$cal_betadiv(unifrac = TRUE, binary = TRUE)
MK.microeco.dataset$cal_betadiv(unifrac = TRUE, binary = TRUE)
PB.microeco.dataset$cal_betadiv(unifrac = TRUE, binary = TRUE)
PP.microeco.dataset$cal_betadiv(unifrac = TRUE, binary = TRUE)
PT.microeco.dataset$cal_betadiv(unifrac = TRUE, binary = TRUE)
RB.microeco.dataset$cal_betadiv(unifrac = TRUE, binary = TRUE)
RN.microeco.dataset$cal_betadiv(unifrac = TRUE, binary = TRUE)
SE.microeco.dataset$cal_betadiv(unifrac = TRUE, binary = TRUE)
SW.microeco.dataset$cal_betadiv(unifrac = TRUE, binary = TRUE)
US.microeco.dataset$cal_betadiv(unifrac = TRUE, binary = TRUE)



# first perform mantel test 
# Nitrate and Nitrites excluded as they were zero
# all other env variables were considered as indication of human activity so not included.
# Total alkalinity covaries with Carbonate so removed
# no lat and long

AP.microeco.mantel.select.noalk <- trans_env$new(dataset = AP.microeco.dataset, env_cols = 23:28)
BA.microeco.mantel.select.noalk <- trans_env$new(dataset = BA.microeco.dataset, env_cols = 23:28)
HS.microeco.mantel.select.noalk <- trans_env$new(dataset = HS.microeco.dataset, env_cols = 23:28)
KJ.microeco.mantel.select.noalk <- trans_env$new(dataset = KJ.microeco.dataset, env_cols = 23:28)
LA.microeco.mantel.select.noalk <- trans_env$new(dataset = LA.microeco.dataset, env_cols = 23:28)
LN.microeco.mantel.select.noalk <- trans_env$new(dataset = LN.microeco.dataset, env_cols = 23:28)
MK.microeco.mantel.select.noalk <- trans_env$new(dataset = MK.microeco.dataset, env_cols = 23:28)
PB.microeco.mantel.select.noalk <- trans_env$new(dataset = PB.microeco.dataset, env_cols = 23:28)
PP.microeco.mantel.select.noalk <- trans_env$new(dataset = PP.microeco.dataset, env_cols = 23:28)
PT.microeco.mantel.select.noalk <- trans_env$new(dataset = PT.microeco.dataset, env_cols = 23:28)
RB.microeco.mantel.select.noalk <- trans_env$new(dataset = RB.microeco.dataset, env_cols = 23:28)
RN.microeco.mantel.select.noalk <- trans_env$new(dataset = RN.microeco.dataset, env_cols = 23:28)
SE.microeco.mantel.select.noalk <- trans_env$new(dataset = SE.microeco.dataset, env_cols = 23:28)
SW.microeco.mantel.select.noalk <- trans_env$new(dataset = SW.microeco.dataset, env_cols = 23:28)
US.microeco.mantel.select.noalk <- trans_env$new(dataset = US.microeco.dataset, env_cols = 23:28)


#a <- AP.microeco.mantel.select.noalk$data_env
#lapply(test, class)

#http://qiime.org/tutorials/distance_matrix_comparison.html#:~:text=The%20partial%20Mantel%20test%20is,three%20distance%20(dissimilarit
#The partial Mantel test is used to estimate the correlation between two matrices, A and B, while controlling for the effect of a control matrix C. The partial Mantel test is a first-order correlation analysis that utilizes three distance (dissimilarity) matrices. This test builds on the simple Mantel test by adding a third “control” matrix. The goal is to test the correlation between matrices A and B while controlling the effect of a third matrix C, in order to remove spurious correlations. The first distance matrix is the one that is permuted so that the correlation structure between the first and second distance matrices is kept constant (Oksanen et al., 2011). A popular use of the partial Mantel test is to compare a community distance matrix with another distance matrix derived from an environmental parameter, using geographic distance as the third “control” distance matrix.

# Full Mantel's test between beta diversity matrix and environmental data.
## No mantel's for locations with no temperature variability

#AP.microeco.mantel.select.noalk$cal_mantel(use_measure = "unwei_unifrac", method = "pearson" , partial_mantel = FALSE)

#BA.microeco.mantel.select.noalk$cal_mantel(use_measure = "unwei_unifrac", method = "pearson" , partial_mantel = FALSE)

HS.microeco.mantel.select.noalk$cal_mantel(use_measure = "unwei_unifrac", method = "pearson", partial_mantel = FALSE)

#KJ.microeco.mantel.select.noalk$cal_mantel(use_measure = "unwei_unifrac", method = "pearson" , partial_mantel = FALSE)

#LA.microeco.mantel.select.noalk$cal_mantel(use_measure = "unwei_unifrac", method = "pearson" , partial_mantel = FALSE)

LN.microeco.mantel.select.noalk$cal_mantel(use_measure = "unwei_unifrac", method = "pearson" , partial_mantel = FALSE)

MK.microeco.mantel.select.noalk$cal_mantel(use_measure = "unwei_unifrac", method = "pearson" , partial_mantel = FALSE)

PB.microeco.mantel.select.noalk$cal_mantel(use_measure = "unwei_unifrac", method = "pearson" , partial_mantel = FALSE)

PP.microeco.mantel.select.noalk$cal_mantel(use_measure = "unwei_unifrac", method = "pearson" , partial_mantel = FALSE)

PT.microeco.mantel.select.noalk$cal_mantel(use_measure = "unwei_unifrac", method = "pearson" , partial_mantel = FALSE)

RB.microeco.mantel.select.noalk$cal_mantel(use_measure = "unwei_unifrac", method = "pearson" , partial_mantel = FALSE)

RN.microeco.mantel.select.noalk$cal_mantel(use_measure = "unwei_unifrac", method = "pearson" , partial_mantel = FALSE)

SE.microeco.mantel.select.noalk$cal_mantel(use_measure = "unwei_unifrac", method = "pearson" , partial_mantel = FALSE)

SW.microeco.mantel.select.noalk$cal_mantel(use_measure = "unwei_unifrac", method = "pearson" , partial_mantel = FALSE)

US.microeco.mantel.select.noalk$cal_mantel(use_measure = "unwei_unifrac", method = "pearson" , partial_mantel = FALSE)


### use microeco.mantel.all to subset for required columns

## no H2S and phosphate , no variability in conductivity, carbonate 
HS.microeco.mantel.select.noalk.all <- data.frame(HS.microeco.mantel.select.noalk$res_mantel)

LN.microeco.mantel.select.noalk.all <- data.frame(LN.microeco.mantel.select.noalk$res_mantel)

MK.microeco.mantel.select.noalk.all <- data.frame(MK.microeco.mantel.select.noalk$res_mantel)

PB.microeco.mantel.select.noalk.all <- data.frame(PB.microeco.mantel.select.noalk$res_mantel)

PP.microeco.mantel.select.noalk.all <- data.frame(PP.microeco.mantel.select.noalk$res_mantel)

PT.microeco.mantel.select.noalk.all <- data.frame(PT.microeco.mantel.select.noalk$res_mantel)

RB.microeco.mantel.select.noalk.all <- data.frame(RB.microeco.mantel.select.noalk$res_mantel)

RN.microeco.mantel.select.noalk.all <- data.frame(RN.microeco.mantel.select.noalk$res_mantel)

SE.microeco.mantel.select.noalk.all <- data.frame(SE.microeco.mantel.select.noalk$res_mantel)

SW.microeco.mantel.select.noalk.all <- data.frame(SW.microeco.mantel.select.noalk$res_mantel)

US.microeco.mantel.select.noalk.all <- data.frame(US.microeco.mantel.select.noalk$res_mantel)


### mantel's variables in a dataframe

HS.All.select.noalk <- data.frame(spec = "HS", HS.microeco.mantel.select.noalk$res_mantel) %>% .[, c(1, 3, 6, 8)]
LN.All.select.noalk <- data.frame(spec = "LN", LN.microeco.mantel.select.noalk$res_mantel) %>% .[, c(1, 3, 6, 8)]
MK.All.select.noalk <- data.frame(spec = "MK", MK.microeco.mantel.select.noalk$res_mantel) %>% .[, c(1, 3, 6, 8)]
PB.All.select.noalk <- data.frame(spec = "PB", PB.microeco.mantel.select.noalk$res_mantel) %>% .[, c(1, 3, 6, 8)]
PP.All.select.noalk <- data.frame(spec = "PP", PP.microeco.mantel.select.noalk$res_mantel) %>% .[, c(1, 3, 6, 8)]
PT.All.select.noalk <- data.frame(spec = "PT", PT.microeco.mantel.select.noalk$res_mantel) %>% .[, c(1, 3, 6, 8)]
RB.All.select.noalk <- data.frame(spec = "RB", RB.microeco.mantel.select.noalk$res_mantel) %>% .[, c(1, 3, 6, 8)]
RN.All.select.noalk <- data.frame(spec = "RN", RN.microeco.mantel.select.noalk$res_mantel) %>% .[, c(1, 3, 6, 8)]
SE.All.select.noalk <- data.frame(spec = "SE", SE.microeco.mantel.select.noalk$res_mantel) %>% .[, c(1, 3, 6, 8)]
SW.All.select.noalk <- data.frame(spec = "SW", SW.microeco.mantel.select.noalk$res_mantel) %>% .[, c(1, 3, 6, 8)]
US.All.select.noalk <- data.frame(spec = "US", US.microeco.mantel.select.noalk$res_mantel) %>% .[, c(1, 3, 6, 8)]


# rename columns
colnames(HS.All.select.noalk) <- colnames(LN.All.select.noalk) <- colnames(MK.All.select.noalk) <- colnames(PB.All.select.noalk) <- colnames(PP.All.select.noalk) <- colnames(PT.All.select.noalk) <-  colnames(RB.All.select.noalk) <- colnames(RN.All.select.noalk) <- colnames(SE.All.select.noalk) <- colnames(SW.All.select.noalk) <- colnames(US.All.select.noalk) <-  c("spec", "env", "r", "p.value")

# generate interval data

HS.All.select.noalk %<>% dplyr::mutate(rd = cut(r, breaks = c(-Inf, 0.3, 0.6, Inf), labels = c("< 0.3", "0.3 - 0.6", ">= 0.6")), pd = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf), labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))

LN.All.select.noalk %<>% dplyr::mutate(rd = cut(r, breaks = c(-Inf, 0.3, 0.6, Inf), labels = c("< 0.3", "0.3 - 0.6", ">= 0.6")), pd = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf), labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))

MK.All.select.noalk %<>% dplyr::mutate(rd = cut(r, breaks = c(-Inf, 0.3, 0.6, Inf), labels = c("< 0.3", "0.3 - 0.6", ">= 0.6")), pd = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf), labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))

PB.All.select.noalk %<>% dplyr::mutate(rd = cut(r, breaks = c(-Inf, 0.3, 0.6, Inf), labels = c("< 0.3", "0.3 - 0.6", ">= 0.6")), pd = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf), labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))

PP.All.select.noalk %<>% dplyr::mutate(rd = cut(r, breaks = c(-Inf, 0.3, 0.6, Inf), labels = c("< 0.3", "0.3 - 0.6", ">= 0.6")), pd = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf), labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))

PT.All.select.noalk %<>% dplyr::mutate(rd = cut(r, breaks = c(-Inf, 0.3, 0.6, Inf), labels = c("< 0.3", "0.3 - 0.6", ">= 0.6")), pd = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf), labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))

RB.All.select.noalk %<>% dplyr::mutate(rd = cut(r, breaks = c(-Inf, 0.3, 0.6, Inf), labels = c("< 0.3", "0.3 - 0.6", ">= 0.6")), pd = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf), labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))

RN.All.select.noalk %<>% dplyr::mutate(rd = cut(r, breaks = c(-Inf, 0.3, 0.6, Inf), labels = c("< 0.3", "0.3 - 0.6", ">= 0.6")), pd = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf), labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))

SE.All.select.noalk %<>% dplyr::mutate(rd = cut(r, breaks = c(-Inf, 0.3, 0.6, Inf), labels = c("< 0.3", "0.3 - 0.6", ">= 0.6")), pd = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf), labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))

SW.All.select.noalk %<>% dplyr::mutate(rd = cut(r, breaks = c(-Inf, 0.3, 0.6, Inf), labels = c("< 0.3", "0.3 - 0.6", ">= 0.6")), pd = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf), labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))

US.All.select.noalk %<>% dplyr::mutate(rd = cut(r, breaks = c(-Inf, 0.3, 0.6, Inf), labels = c("< 0.3", "0.3 - 0.6", ">= 0.6")), pd = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf), labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))


# combine all locations
plot_table_select_noalk_locations <- rbind(HS.All.select.noalk, LN.All.select.noalk, MK.All.select.noalk, PB.All.select.noalk, PP.All.select.noalk, PT.All.select.noalk, RB.All.select.noalk, RN.All.select.noalk, SE.All.select.noalk, SW.All.select.noalk,  US.All.select.noalk )


#########################################################################################

### Fig. S11A. Variance partitioning of major abiotic influences on microbial communities with distance (latitude) included as a variable. #####

######################################################################################

### Create locations dataset for locations that have a partial Mantel's values:

## HS, LN, MK, PB, PP, PT, RB, RN, SE, SW,US

HS.rarefied.min.int.exclude.1.3.rooted <- com.rarefied.min.int.exclude.1.3.rooted %>%
  subset_samples(Location.Code == "HS")

HS.rarefied.min.int.exclude.1.3.rooted <- 
  prune_taxa(taxa_sums(HS.rarefied.min.int.exclude.1.3.rooted) > 0, HS.rarefied.min.int.exclude.1.3.rooted)

any(taxa_sums(HS.rarefied.min.int.exclude.1.3.rooted) == 0)

test <- as.data.frame(HS.rarefied.min.int.exclude.1.3.rooted@otu_table)


LN.rarefied.min.int.exclude.1.3.rooted <- com.rarefied.min.int.exclude.1.3.rooted %>%
  subset_samples(Location.Code == "LN")

LN.rarefied.min.int.exclude.1.3.rooted <- 
  prune_taxa(taxa_sums(LN.rarefied.min.int.exclude.1.3.rooted) > 0, LN.rarefied.min.int.exclude.1.3.rooted)

any(taxa_sums(LN.rarefied.min.int.exclude.1.3.rooted) == 0)

test <- as.data.frame(LN.rarefied.min.int.exclude.1.3.rooted@otu_table)


MK.rarefied.min.int.exclude.1.3.rooted <- com.rarefied.min.int.exclude.1.3.rooted %>%
  subset_samples(Location.Code == "MK")

MK.rarefied.min.int.exclude.1.3.rooted <- 
  prune_taxa(taxa_sums(MK.rarefied.min.int.exclude.1.3.rooted) > 0, MK.rarefied.min.int.exclude.1.3.rooted)

any(taxa_sums(MK.rarefied.min.int.exclude.1.3.rooted) == 0)

test <- as.data.frame(MK.rarefied.min.int.exclude.1.3.rooted@otu_table)


PB.rarefied.min.int.exclude.1.3.rooted <- com.rarefied.min.int.exclude.1.3.rooted %>%
  subset_samples(Location.Code == "PB")

PB.rarefied.min.int.exclude.1.3.rooted <- 
  prune_taxa(taxa_sums(PB.rarefied.min.int.exclude.1.3.rooted) > 0, PB.rarefied.min.int.exclude.1.3.rooted)

any(taxa_sums(PB.rarefied.min.int.exclude.1.3.rooted) == 0)

test <- as.data.frame(PB.rarefied.min.int.exclude.1.3.rooted@otu_table)


PP.rarefied.min.int.exclude.1.3.rooted <- com.rarefied.min.int.exclude.1.3.rooted %>%
  subset_samples(Location.Code == "PP")

PP.rarefied.min.int.exclude.1.3.rooted <- 
  prune_taxa(taxa_sums(PP.rarefied.min.int.exclude.1.3.rooted) > 0, PP.rarefied.min.int.exclude.1.3.rooted)

any(taxa_sums(PP.rarefied.min.int.exclude.1.3.rooted) == 0)

test <- as.data.frame(PP.rarefied.min.int.exclude.1.3.rooted@otu_table)


PT.rarefied.min.int.exclude.1.3.rooted <- com.rarefied.min.int.exclude.1.3.rooted %>%
  subset_samples(Location.Code == "PT")

PT.rarefied.min.int.exclude.1.3.rooted <- 
  prune_taxa(taxa_sums(PT.rarefied.min.int.exclude.1.3.rooted) > 0, PT.rarefied.min.int.exclude.1.3.rooted)

any(taxa_sums(PT.rarefied.min.int.exclude.1.3.rooted) == 0)

test <- as.data.frame(PT.rarefied.min.int.exclude.1.3.rooted@otu_table)


RB.rarefied.min.int.exclude.1.3.rooted <- com.rarefied.min.int.exclude.1.3.rooted %>%
  subset_samples(Location.Code == "RB")

RB.rarefied.min.int.exclude.1.3.rooted <- 
  prune_taxa(taxa_sums(RB.rarefied.min.int.exclude.1.3.rooted) > 0, RB.rarefied.min.int.exclude.1.3.rooted)

any(taxa_sums(RB.rarefied.min.int.exclude.1.3.rooted) == 0)

test <- as.data.frame(RB.rarefied.min.int.exclude.1.3.rooted@otu_table)


RN.rarefied.min.int.exclude.1.3.rooted <- com.rarefied.min.int.exclude.1.3.rooted %>%
  subset_samples(Location.Code == "RN")

RN.rarefied.min.int.exclude.1.3.rooted <- 
  prune_taxa(taxa_sums(RN.rarefied.min.int.exclude.1.3.rooted) > 0, RN.rarefied.min.int.exclude.1.3.rooted)

any(taxa_sums(RN.rarefied.min.int.exclude.1.3.rooted) == 0)

test <- as.data.frame(RN.rarefied.min.int.exclude.1.3.rooted@otu_table)


SE.rarefied.min.int.exclude.1.3.rooted <- com.rarefied.min.int.exclude.1.3.rooted %>%
  subset_samples(Location.Code == "SE")

SE.rarefied.min.int.exclude.1.3.rooted <- 
  prune_taxa(taxa_sums(SE.rarefied.min.int.exclude.1.3.rooted) > 0, SE.rarefied.min.int.exclude.1.3.rooted)

any(taxa_sums(SE.rarefied.min.int.exclude.1.3.rooted) == 0)

test <- as.data.frame(SE.rarefied.min.int.exclude.1.3.rooted@otu_table)


SW.rarefied.min.int.exclude.1.3.rooted <- com.rarefied.min.int.exclude.1.3.rooted %>%
  subset_samples(Location.Code == "SW")

SW.rarefied.min.int.exclude.1.3.rooted <- 
  prune_taxa(taxa_sums(SW.rarefied.min.int.exclude.1.3.rooted) > 0, SW.rarefied.min.int.exclude.1.3.rooted)

any(taxa_sums(SW.rarefied.min.int.exclude.1.3.rooted) == 0)

test <- as.data.frame(SW.rarefied.min.int.exclude.1.3.rooted@otu_table)


US.rarefied.min.int.exclude.1.3.rooted <- com.rarefied.min.int.exclude.1.3.rooted %>%
  subset_samples(Location.Code == "US")

US.rarefied.min.int.exclude.1.3.rooted <- 
  prune_taxa(taxa_sums(US.rarefied.min.int.exclude.1.3.rooted) > 0, US.rarefied.min.int.exclude.1.3.rooted)

any(taxa_sums(US.rarefied.min.int.exclude.1.3.rooted) == 0)

test <- as.data.frame(US.rarefied.min.int.exclude.1.3.rooted@otu_table)


## Reducing the weight of rare species use hellinger transformation ####
otu <- as.matrix(as.data.frame(com.rarefied.min.int.exclude.1.3.rooted@otu_table))

NT.otu <- as.matrix(as.data.frame(NT.rarefied.min.int.exclude.1.3.rooted.spnet@otu_table))
NM.otu <- as.matrix(as.data.frame(NM.rarefied.min.int.exclude.1.3.rooted.spnet@otu_table))
ST.otu <- as.matrix(as.data.frame(ST.rarefied.min.int.exclude.1.3.rooted.spnet@otu_table))
SM.otu <- as.matrix(as.data.frame(SM.rarefied.min.int.exclude.1.3.rooted.spnet@otu_table))
CT.otu <- as.matrix(as.data.frame(CT.rarefied.min.int.exclude.1.3.rooted.spnet@otu_table))
SG.otu <- as.matrix(as.data.frame(SG.rarefied.min.int.exclude.1.3.rooted.spnet@otu_table))

HS.otu <- as.matrix(as.data.frame(HS.rarefied.min.int.exclude.1.3.rooted@otu_table))
LN.otu <- as.matrix(as.data.frame(LN.rarefied.min.int.exclude.1.3.rooted@otu_table))
MK.otu <- as.matrix(as.data.frame(MK.rarefied.min.int.exclude.1.3.rooted@otu_table))
PB.otu <- as.matrix(as.data.frame(PB.rarefied.min.int.exclude.1.3.rooted@otu_table))
PP.otu <- as.matrix(as.data.frame(PP.rarefied.min.int.exclude.1.3.rooted@otu_table))
PT.otu <- as.matrix(as.data.frame(PT.rarefied.min.int.exclude.1.3.rooted@otu_table))
RB.otu <- as.matrix(as.data.frame(RB.rarefied.min.int.exclude.1.3.rooted@otu_table))
RN.otu <- as.matrix(as.data.frame(RN.rarefied.min.int.exclude.1.3.rooted@otu_table))
SE.otu <- as.matrix(as.data.frame(SE.rarefied.min.int.exclude.1.3.rooted@otu_table))
SW.otu <- as.matrix(as.data.frame(SW.rarefied.min.int.exclude.1.3.rooted@otu_table))
US.otu <- as.matrix(as.data.frame(US.rarefied.min.int.exclude.1.3.rooted@otu_table))


## get metadata #####
metadata <- as.data.frame(sample_data(com.rarefied.min.int.exclude.1.3.rooted)) %>% as_tibble()
colnames(metadata)
metadata.varpart <- metadata %>% select("Sample_id","Country","Region","Location.Code","Site","Location","Human.usage..Y.N.","Water.flow..Pool.Flowing.","Copper..mg.L.","Cyanuric.acid..mg.L.","Nitrate..ppm.","Nitrite..ppm.","Temp.adj...C.","Chloride..mg.L.","Chlorine.dioxide..mg.L.","Flouride..mg.L.","Hardness..mg.L.","Lead..mg.L.","MPS..mg.L.","QUAT.QAC..mg.L.","Total.Chlorine..mg.l.","Iron..ppm.","Latittude","Longitude","Carbonate..ppm.","pH","Total.alkalinity..ppm.","Phosphate..ppm.","H2S..ppm.","Temp...C.","EC..mS.")

rownames(metadata.varpart) <- metadata.varpart$Sample_id

metadata.varpart <- as.data.frame(as.matrix(metadata.varpart))

metadata.varpart[11:ncol(metadata.varpart)] <- lapply(metadata.varpart[11:ncol(metadata.varpart)], as.numeric)

sapply(metadata.varpart,class)

rownames(metadata.varpart) <- NULL

metadata.varpart <- data.frame(metadata.varpart)

colnames(metadata)

sapply(metadata.varpart,class)

NT.metadata <-  as.data.frame(sample_data(NT.rarefied.min.int.exclude.1.3.rooted.spnet))
rownames(NT.metadata) <- NULL
NT.metadata <- data.frame(NT.metadata)

NM.metadata <-  as.data.frame(sample_data(NM.rarefied.min.int.exclude.1.3.rooted.spnet))
rownames(NM.metadata) <- NULL
NM.metadata <- data.frame(NM.metadata)

ST.metadata <-  as.data.frame(sample_data(ST.rarefied.min.int.exclude.1.3.rooted.spnet))
rownames(ST.metadata) <- NULL
ST.metadata <- data.frame(ST.metadata)

SM.metadata <-  as.data.frame(sample_data(SM.rarefied.min.int.exclude.1.3.rooted.spnet))
rownames(SM.metadata) <- NULL
SM.metadata <- data.frame(SM.metadata)

CT.metadata <-  as.data.frame(sample_data(CT.rarefied.min.int.exclude.1.3.rooted.spnet))
rownames(CT.metadata) <- NULL
CT.metadata <- data.frame(CT.metadata)

SG.metadata <-  as.data.frame(sample_data(SG.rarefied.min.int.exclude.1.3.rooted.spnet))
rownames(SG.metadata) <- NULL
SG.metadata <- data.frame(SG.metadata)



HS.metadata <-  as.data.frame(sample_data(HS.rarefied.min.int.exclude.1.3.rooted))
rownames(HS.metadata) <- NULL
HS.metadata <- data.frame(HS.metadata)

LN.metadata <-  as.data.frame(sample_data(LN.rarefied.min.int.exclude.1.3.rooted))
rownames(LN.metadata) <- NULL
LN.metadata <- data.frame(LN.metadata)

MK.metadata <-  as.data.frame(sample_data(MK.rarefied.min.int.exclude.1.3.rooted))
rownames(MK.metadata) <- NULL
MK.metadata <- data.frame(MK.metadata)

PB.metadata <-  as.data.frame(sample_data(PB.rarefied.min.int.exclude.1.3.rooted))
rownames(PB.metadata) <- NULL
PB.metadata <- data.frame(PB.metadata)

PP.metadata <-  as.data.frame(sample_data(PP.rarefied.min.int.exclude.1.3.rooted))
rownames(PP.metadata) <- NULL
PP.metadata <- data.frame(PP.metadata)

PT.metadata <-  as.data.frame(sample_data(PT.rarefied.min.int.exclude.1.3.rooted))
rownames(PT.metadata) <- NULL
PT.metadata <- data.frame(PT.metadata)

RB.metadata <-  as.data.frame(sample_data(RB.rarefied.min.int.exclude.1.3.rooted))
rownames(RB.metadata) <- NULL
RB.metadata <- data.frame(RB.metadata)

RN.metadata <-  as.data.frame(sample_data(RN.rarefied.min.int.exclude.1.3.rooted))
rownames(RN.metadata) <- NULL
RN.metadata <- data.frame(RN.metadata)

SE.metadata <-  as.data.frame(sample_data(SE.rarefied.min.int.exclude.1.3.rooted))
rownames(SE.metadata) <- NULL
SE.metadata <- data.frame(SE.metadata)

SW.metadata <-  as.data.frame(sample_data(SW.rarefied.min.int.exclude.1.3.rooted))
rownames(SW.metadata) <- NULL
SW.metadata <- data.frame(SW.metadata)

US.metadata <-  as.data.frame(sample_data(US.rarefied.min.int.exclude.1.3.rooted))
rownames(US.metadata) <- NULL
US.metadata <- data.frame(US.metadata)



### hellinger transform #####
otu.hellinger.transform <- decostand(otu, method = "hellinger")
NT.hellinger.transform <- decostand(NT.otu, method = "hellinger")
NM.hellinger.transform <- decostand(NM.otu, method = "hellinger")
ST.hellinger.transform <- decostand(ST.otu, method = "hellinger")
SM.hellinger.transform <- decostand(SM.otu, method = "hellinger")
CT.hellinger.transform <- decostand(CT.otu, method = "hellinger")
SG.hellinger.transform <- decostand(SG.otu, method = "hellinger")

HS.hellinger.transform <- decostand(HS.otu, method = "hellinger")
LN.hellinger.transform <- decostand(LN.otu, method = "hellinger")
MK.hellinger.transform <- decostand(MK.otu, method = "hellinger")
PB.hellinger.transform <- decostand(PB.otu, method = "hellinger")
PP.hellinger.transform <- decostand(PP.otu, method = "hellinger")
PT.hellinger.transform <- decostand(PT.otu, method = "hellinger")
RB.hellinger.transform <- decostand(RB.otu, method = "hellinger")
RN.hellinger.transform <- decostand(RN.otu, method = "hellinger")
SE.hellinger.transform <- decostand(SE.otu, method = "hellinger")
SW.hellinger.transform <- decostand(SW.otu, method = "hellinger")
US.hellinger.transform <- decostand(US.otu, method = "hellinger")



# Run a PCA on the Hellinger-transformed sh data and extract the results ####
# Copper Cyanuric acid Fluoride is 0

## 1st fraction :
otu.hellinger.rDA.1 <- rda(otu.hellinger.transform ~ Temp...C. + H2S..ppm. + EC..mS. + Carbonate..ppm., data = metadata.varpart)

otu.hellinger.rDA.2 <- rda(otu.hellinger.transform ~ pH + Chloride..mg.L. + Chlorine.dioxide..mg.L. + Hardness..mg.L., data = metadata.varpart)

otu.hellinger.rDA.3 <- rda(otu.hellinger.transform ~ Iron..ppm. + Lead..mg.L. + MPS..mg.L. + Phosphate..ppm. , data = metadata.varpart)

otu.hellinger.rDA.4 <- rda(otu.hellinger.transform ~ QUAT.QAC..mg.L. + Total.alkalinity..ppm. + Total.Chlorine..mg.l., data = metadata.varpart)

otu.hellinger.rDA.all <- rda(otu.hellinger.transform ~ Temp...C. + H2S..ppm. + EC..mS. + Carbonate..ppm. + pH + Chloride..mg.L. + Chlorine.dioxide..mg.L. + Hardness..mg.L. + Iron..ppm. + Lead..mg.L. + MPS..mg.L. + Phosphate..ppm. +  QUAT.QAC..mg.L. + Total.alkalinity..ppm. + Total.Chlorine..mg.l., data = metadata.varpart)

otu.hellinger.chosen <- rda(otu.hellinger.transform ~ Temp...C. + H2S..ppm. + EC..mS. + pH, data = metadata.varpart)

otu.hellinger.chosen

otu.hellinger.chosen.2 <- rda(otu.hellinger.transform ~ Latittude + Carbonate..ppm. + EC..mS. + pH, data = metadata.varpart)

otu.hellinger.chosen.2

otu.hellinger.chosen.3 <- rda(otu.hellinger.transform ~ Latittude + Condition(Carbonate..ppm.) + Condition(EC..mS.) + Condition(pH), data = metadata.varpart)

otu.hellinger.chosen.3
anova(otu.hellinger.chosen.3)

otu.hellinger.chosen.4 <- rda(otu.hellinger.transform ~ Condition(Latittude) + (Carbonate..ppm.) + Condition(EC..mS.) + Condition(pH), data = metadata.varpart)

otu.hellinger.chosen.4
anova(otu.hellinger.chosen.4)

otu.hellinger.chosen.5 <- rda(otu.hellinger.transform ~ Condition(Latittude) + Condition(Carbonate..ppm.) + (EC..mS.) + Condition(pH), data = metadata.varpart)

otu.hellinger.chosen.5
anova(otu.hellinger.chosen.5)

otu.hellinger.chosen.6 <- rda(otu.hellinger.transform ~ Condition(Latittude) + Condition(Carbonate..ppm.) + Condition(EC..mS.) + (pH), data = metadata.varpart)

otu.hellinger.chosen.6
anova(otu.hellinger.chosen.6)

otu.hellinger.chosen.7 <- rda(otu.hellinger.transform ~ H2S..ppm. + EC..mS. + Carbonate..ppm. + pH , data = metadata.varpart)

otu.hellinger.chosen.8 <- rda(otu.hellinger.transform ~ Condition(H2S..ppm.) + Condition(EC..mS.) + Carbonate..ppm. + pH , data = metadata.varpart)

##
otu.hellinger.rDA.temp <- rda(otu.hellinger.transform ~ Temp...C., data = metadata.varpart)

otu.hellinger.rDA.Sul <- rda(otu.hellinger.transform ~ H2S..ppm., data = metadata.varpart)
otu.hellinger.rDA.EC <- rda(otu.hellinger.transform ~ EC..mS., data = metadata.varpart)
##otu.hellinger.rDA.Carb <- rda(otu.hellinger.transform ~ Carbonate..ppm., data = metadata.varpart)
otu.hellinger.rDA.pH <- rda(otu.hellinger.transform ~ pH, data = metadata.varpart)
otu.hellinger.rDA.carb <- rda(otu.hellinger.transform ~ Carbonate..ppm., data = metadata.varpart)

#otu.hellinger.rDA.Chloride <- rda(otu.hellinger.transform ~ Chloride..mg.L., data = metadata.varpart)
#otu.hellinger.rDA.ChloDio <- rda(otu.hellinger.transform ~ Chlorine.dioxide..mg.L., data = metadata.varpart)
#otu.hellinger.rDA.Hard <- rda(otu.hellinger.transform ~ Hardness..mg.L., data = metadata.varpart)
#otu.hellinger.rDA.Iron <- rda(otu.hellinger.transform ~ Iron..ppm., data = metadata.varpart)
#otu.hellinger.rDA.Lead <- rda(otu.hellinger.transform ~ Lead..mg.L., data = metadata.varpart)
#otu.hellinger.rDA.MPS <- rda(otu.hellinger.transform ~ MPS..mg.L., data = metadata.varpart)
#otu.hellinger.rDA.Phos <- rda(otu.hellinger.transform ~ Phosphate..ppm. , data = metadata.varpart)
#otu.hellinger.rDA.Q <- rda(otu.hellinger.transform ~ QUAT.QAC..mg.L., data = metadata.varpart)
#otu.hellinger.rDA.alk <- rda(otu.hellinger.transform ~ Total.alkalinity..ppm. , data = metadata.varpart)
#otu.hellinger.rDA.4 <- rda(otu.hellinger.transform ~  Total.Chlorine..mg.l., data = metadata.varpart)

#summary(otu.hellinger.rDA)


## to plot using https://fukamilab.github.io/BIO202/06-B-constrained-ordination.html

## RUN varpart #####
##Either
#varp.1 <- varpart (otu, ~ `Temp (oC)`, ~ `H2S (ng/l)`, ~ `Carbonate (mg/l)` , data = env.samples.wide,transfo = "hel")

## OR
varp.1 <- varpart (otu.hellinger.transform, ~ Temp...C., ~ H2S..ppm., ~ EC..mS., ~ Carbonate..ppm.,data = metadata.varpart)

varp.1

varp.2 <- varpart (otu.hellinger.transform, ~ pH, ~ Chloride..mg.L., ~ Chlorine.dioxide..mg.L., ~ Hardness..mg.L.,data = metadata.varpart)

varp.3 <- varpart (otu.hellinger.transform, ~ Iron..ppm., ~ Lead..mg.L., ~ MPS..mg.L., ~ Phosphate..ppm. ,data = metadata.varpart)

varp.4 <- varpart (otu.hellinger.transform, ~ QUAT.QAC..mg.L., ~ Total.alkalinity..ppm., ~ Total.Chlorine..mg.l. ,data = metadata.varpart)

varp.chosen <- varpart (otu.hellinger.transform, ~ EC..mS., ~Temp...C., ~ H2S..ppm., ~ pH, data = metadata.varpart)

varp.chosen

varp.chosen.2 <- varpart (otu.hellinger.transform, ~  EC..mS., ~ Latittude, ~ Carbonate..ppm., ~ pH, data = metadata.varpart)

varp.chosen.2$part$fract

#get partition table
summary(varp.chosen.2)
## indicates testable fractions
varp.chosen.2$part$fract

## overall ######
varp.chosen.3 <- varpart (otu.hellinger.transform, ~  EC..mS., ~ pH,  ~ Carbonate..ppm., ~ H2S..ppm., data = metadata.varpart)

#get partition table
summary(varp.chosen.3)
## indicates testable fractions
varp.chosen.3$part$fract


# We may plot the results into Venn's diagram (argument digits influences number of decimal digits shown in the diagram, Xnames the displayed names of variables, and bg background color of the diagram fractions; see ?varpart for details):

plot (varp.1, digits = 2, Xnames = c('`Temp (oC)`', '`H2S (ppm)`', 'EC (mS)', 'Carbonate (ppm)'), bg = c('maroon', 'darkgoldenrod1', "turquoise", "mediumslateblue"))

plot (varp.2, digits = 2, Xnames = c('`pH`', '`Chloride (mg/L)`', 'Chlorine.dioxide (mg/L)', 'Hardness (mg/L)'), bg = c('maroon', 'darkgoldenrod1', "turquoise", "mediumslateblue"))

plot (varp.3, digits = 2, Xnames = c('`Iron (ppm)`', '`Lead (mg/L)`', 'MPS (mg/L)', 'Phosphate (ppm)'), bg = c('maroon', 'darkgoldenrod1', "turquoise", "mediumslateblue"))

plot (varp.4, digits = 2, Xnames = c('` QUAT.QAC (mg/L)`', '`Total alkalinity (ppm)`', 'Total Chlorine (mg/L)', 'Phosphate (ppm)'), bg = c('maroon', 'darkgoldenrod1', "turquoise"))

pdf("varpart.pdf", width = 8)
plot (varp.chosen, digits = 2, Xnames = c('EC (mS)', 'Temp (°C)', 'H2S (ppm)','pH'), bg = c('turquoise', 'darkgoldenrod1', "maroon", "mediumslateblue"))
dev.off()



pdf("varpart.2.pdf", width = 8)
plot (varp.chosen.2, digits = 2, Xnames = c('EC (mS)','Latittude (D)', 'Carbonate (ppm)', 'pH'), bg = c('turquoise','khaki3', 'greenyellow', "mediumslateblue"))
dev.off()

pdf("varpart.3.pdf", width = 8)
plot (varp.chosen.3, digits = 2, Xnames = c('EC (mS)','pH', 'Carbonate (ppm)', 'Sulfide (ppm)'), bg = c('turquoise','khaki3', 'greenyellow', "mediumslateblue"))
dev.off()



## import mantel's table for all

plot_table_select_noalk_all_reg <- readRDS("plot_table_select_noalk_all_reg.rds")
plot_table_select_noalk_locations <- readRDS("plot_table_select_noalk_locations.rds")



### Varpart by regions ####


## NT - all varibales with variation ####
varp.chosen.NT <- varpart (NT.hellinger.transform, ~  EC..mS., ~ pH,  ~ Carbonate..ppm., ~ H2S..ppm., data = NT.metadata)

#get partition table
summary(varp.chosen.NT)
## indicates testable fractions
varp.chosen.NT$part$fract

pdf("varpart.NT.pdf", width = 8)
plot (varp.chosen.NT, digits = 2, Xnames = c('EC (mS)', 'pH','Carbonate (ppm)','Sulfide (ppm)' ), bg = c('turquoise','khaki3','greenyellow',"mediumslateblue"))
dev.off()


## NM - no H2S variation #####
## BOTH ph and phosphate with same r value in Mantel's

# WARNING MESSAGE - Warning: collinearity detected in cbind(X1,X3: mm = 2, m = 1Warning: collinearity detected in cbind(X1,X2,X3): mm = 3, m = 2Warning: collinearity detected in cbind(X1,X3,X4): mm = 3, m = 2Warning: collinearity detected in cbind(X1,X2,X3,X4): mm = 4, m = 3
varp.chosen.NM.ph <- varpart (NM.hellinger.transform, ~  EC..mS., ~ pH,  ~ Carbonate..ppm., ~ Temp...C., data = NM.metadata)

## WARNING MESSAGE - # WARNING MESSAGE - Warning: collinearity detected in cbind(X1,X3: mm = 2, m = 1Warning: collinearity detected in cbind(X1,X2,X3): mm = 3, m = 2Warning: collinearity detected in cbind(X1,X3,X4): mm = 3, m = 2Warning: collinearity detected in cbind(X1,X2,X3,X4): mm = 4, m = 3
varp.chosen.NM.phos <- varpart (NM.hellinger.transform, ~  EC..mS., ~ Phosphate..ppm.,  ~ Carbonate..ppm., ~ Temp...C., data = NM.metadata)

#get partition table
summary(varp.chosen.NM.ph)
summary(varp.chosen.NM.phos)
## indicates testable fractions
varp.chosen.NM.ph$part$fract
varp.chosen.NM.phos$part$fract

pdf("varpart.NM.ph.pdf", width = 8)
plot (varp.chosen.NM.ph, digits = 2, Xnames = c('EC (mS)', 'pH','Carbonate (ppm)','Temperature (°C)' ), bg = c('turquoise','khaki3','greenyellow',"mediumslateblue"))
dev.off()

pdf("varpart.NM.phos.pdf", width = 8)
plot (varp.chosen.NM.phos, digits = 2, Xnames = c('EC (mS)', 'Phosphate (ppm)' ,'Carbonate (ppm)','Temperature (°C)' ), bg = c('turquoise','khaki3','greenyellow',"mediumslateblue"))
dev.off()

## ST
varp.chosen.ST <- varpart (ST.hellinger.transform, ~  EC..mS., ~ pH,  ~ Carbonate..ppm., ~ Temp...C., data = ST.metadata)

#get partition table
summary(varp.chosen.ST)
## indicates testable fractions
varp.chosen.ST$part$fract

pdf("varpart.ST.pdf", width = 8)
plot (varp.chosen.ST, digits = 2, Xnames = c('EC (mS)', 'pH' ,'Carbonate (ppm)','Temperature (°C)' ), bg = c('turquoise','khaki3','greenyellow',"mediumslateblue"))
dev.off()

## SM
varp.chosen.SM <- varpart (SM.hellinger.transform, ~  Phosphate..ppm., ~ pH,  ~ Carbonate..ppm., ~ Temp...C., data = SM.metadata)

#get partition table
summary(varp.chosen.SM)
## indicates testable fractions
varp.chosen.SM$part$fract

pdf("varpart.SM.pdf", width = 8)
plot (varp.chosen.SM, digits = 2, Xnames = c( 'Phosphate (ppm)' , 'pH' ,'Carbonate (ppm)','Temperature (°C)' ), bg = c('turquoise','khaki3','greenyellow',"mediumslateblue"))
dev.off()



## SG
varp.chosen.SG <- varpart (SG.hellinger.transform, ~ H2S..ppm., ~ Carbonate..ppm., ~ Temp...C., data = SG.metadata)

#get partition table
summary(varp.chosen.SG)
## indicates testable fractions
varp.chosen.SG$part$fract


pdf("varpart.SG.pdf", width = 8)
plot (varp.chosen.SG, digits = 2, Xnames = c( 'Sulfide (ppm)' , 'pH' ,'Carbonate (ppm)','Temperature (°C)' ), bg = c('turquoise','khaki3','greenyellow'))
dev.off()


## CT
varp.chosen.CT <- varpart (CT.hellinger.transform, ~ EC..mS., ~ Phosphate..ppm. , ~ Carbonate..ppm., ~ Temp...C., data = CT.metadata)

#get partition table
summary(varp.chosen.CT)
## indicates testable fractions
varp.chosen.CT$part$fract


pdf("varpart.CT.pdf", width = 8)
plot (varp.chosen.CT, digits = 2, Xnames = c('EC (mS)', 'Phosphate (ppm)','Carbonate (ppm)', 'Temperature (°C)'), bg = c('turquoise','khaki3','greenyellow',"mediumslateblue"))
dev.off()



#### Locations varpart #####

## HS
varp.chosen.HS <- varpart (HS.hellinger.transform, ~ pH, ~ Temp...C., data = HS.metadata)

#get partition table
summary(varp.chosen.HS)
## indicates testable fractions
varp.chosen.HS$part$fract


pdf("varpart.HS.pdf", width = 8)
plot (varp.chosen.HS, digits = 2, Xnames = c('pH', 'Temp (°C)'), bg = c('turquoise','khaki3'))
dev.off()

## LN
varp.chosen.LN <- varpart (LN.hellinger.transform, ~ pH, ~ Temp...C., ~ Carbonate..ppm.,  data = LN.metadata)

#get partition table
summary(varp.chosen.LN)
## indicates testable fractions
varp.chosen.LN$part$fract


pdf("varpart.LN.pdf", width = 8)
plot (varp.chosen.LN, digits = 2, Xnames = c('pH', 'Temp (°C)','Carbonate (ppm)'), bg = c('turquoise','khaki3','greenyellow'))
dev.off()


## PB
varp.chosen.PB <- varpart (PB.hellinger.transform, ~ EC..mS., ~ Phosphate..ppm. ,  ~ Temp...C., ~ Carbonate..ppm.,  data = PB.metadata)

#get partition table
summary(varp.chosen.PB)
## indicates testable fractions
varp.chosen.PB$part$fract


pdf("varpart.PB.pdf", width = 8)
plot (varp.chosen.PB, digits = 2, Xnames = c('EC (mS)', 'Phosphate (ppm)', 'Temp (°C)','Carbonate (ppm)'), bg = c('turquoise','khaki3','greenyellow',"mediumslateblue"))
dev.off()


## PP
varp.chosen.PP <- varpart (PP.hellinger.transform, ~ Temp...C., ~ Carbonate..ppm.,  data = PP.metadata)

#get partition table
summary(varp.chosen.PP)
## indicates testable fractions
varp.chosen.PP$part$fract


pdf("varpart.PP.pdf", width = 8)
plot (varp.chosen.PP, digits = 2, Xnames = c('Temp (°C)','Carbonate (ppm)'), bg = c('turquoise','khaki3'))
dev.off()


## PT
varp.chosen.PT <- varpart (PT.hellinger.transform, ~ Temp...C., ~ Carbonate..ppm.,  data = PT.metadata)

#get partition table
summary(varp.chosen.PT)
## indicates testable fractions
varp.chosen.PT$part$fract


pdf("varpart.PT.pdf", width = 8)
plot (varp.chosen.PT, digits = 2, Xnames = c('Temp (°C)','Carbonate (ppm)'), bg = c('turquoise','khaki3'))
dev.off()


## RB
varp.chosen.RB <- varpart (RB.hellinger.transform, ~ pH, ~ Temp...C., ~ EC..mS.,  data = RB.metadata)

#get partition table
summary(varp.chosen.RB)
## indicates testable fractions
varp.chosen.RB$part$fract


pdf("varpart.RB.pdf", width = 8)
plot (varp.chosen.RB, digits = 2, Xnames = c('pH', 'Temp (°C)','EC (mS)'), bg = c('turquoise','khaki3','greenyellow'))
dev.off()


## RN
varp.chosen.RN <- varpart (RN.hellinger.transform, ~ Temp...C., ~ EC..mS.,  data = RN.metadata)

#get partition table
summary(varp.chosen.RN)
## indicates testable fractions
varp.chosen.RN$part$fract


pdf("varpart.RN.pdf", width = 8)
plot (varp.chosen.RN, digits = 2, Xnames = c( 'Temp (°C)','EC (mS)'), bg = c('turquoise','khaki3'))
dev.off()


## SW
varp.chosen.SW <- varpart (SW.hellinger.transform, ~ H2S..ppm., ~ Carbonate..ppm., ~ Temp...C.,  data = SW.metadata)

#get partition table
summary(varp.chosen.SW)
## indicates testable fractions
varp.chosen.SW$part$fract


pdf("varpart.SW.pdf", width = 8)
plot (varp.chosen.SW, digits = 2, Xnames = c( '`H2S (ppm)`','Carbonate (ppm)','Temp (°C)'), bg = c('turquoise','khaki3'))
dev.off()




#### Boxplot for varpart for % explained

## create dataframe

Scale <- c('Overall', 'Region', 'Region', 'Region', 'Region', 'Region', 'Region', 'Location', 'Location', 'Location', 'Location', 'Location', 'Location','Location','Location')

Scale.1 <- c('overall', 'North.Thailand', 'Central.Thailand', 'South.Thailand', 'North.Malaysia', 'South.Malaysia', 'Singpore', 'HS', 'LN', 'PB', 'PP', 'PT', 'RB','RN','SW')

### 100 - residuals
percent.explained <- c(26.8, 48.3, 62.3, 55.7, 49.7, 66.4, 47.4, 34.6, 84.8, 62.3, 78.2,66.3, 72.3, 55.0, 47.4)

varpart.box <-  data.frame(Scale,Scale.1, percent.explained)

varpart.box %>% group_by(Scale) %>%
  summarise(Mean = mean(percent.explained)) 


varpart.boxplot <- varpart.box %>% ggplot( aes(x = reorder(Scale,percent.explained), y = percent.explained, fill = Scale)) + geom_boxplot() +
  scale_fill_manual(values=c("#69b3a2","lightskyblue2", "sienna2")) + labs( y = "% variation explained") + theme(legend.key=element_blank(), axis.title = element_text(colour = "black", face = "bold", size = 8) , axis.text.x = element_text(colour = "black", size = 12, face = "bold", angle = 0, hjust = 0.5),   axis.text.y = element_text(colour = "black", face = "bold", size = 12), legend.title = element_text(size = 10, face = "bold"), panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2), legend.position = "right", title = element_text(size = 8))   + theme_bw()

ggsave("varpart.boxplot.pdf", height = 6, width = 5)
ggsave("varpart.boxplot.png", height = 6, width = 5)






###################################################################################

### Fig. 6A Net relatedness index (betaNRI) measures of the mean phylogenetic distance to the nearest taxon in the community by sites. ####
### Fig. 6B. Null model employed to predict the influence of various evolutionary drivers on community assembly for different spatial scales and ecological guilds#####

######################################################################################

### change environmental variable columns accordingly

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

pdf("betaNRI.nullmodel.pdf", width = 15, height =10)
betaNRI + geom_hline(yintercept = -2, linetype = 2) + geom_hline(yintercept = 2, linetype = 2)
dev.off()

### frequency ######
microeco.region.frequency.nullmodel <- readRDS("microeco.nullmodel.rds")

microeco.region.frequency.nullmodel$cal_mantel_corr()
saveRDS(microeco.region.frequency.nullmodel, "microeco.region.frequency.nullmodel.rds")

microeco.region.frequency.nullmodel$cal_ses_betampd(runs = 1000, abundance.weighted = TRUE, null.model = "frequency")
saveRDS(microeco.region.frequency.nullmodel, "microeco.region.frequency.nullmodel.rds")

tmp.3 <- "/hpctmp/chrisg25/R_out/TMFZ/postDADA2/tmp.3" ; dir.create(tmp.3)
microeco.region.frequency.nullmodel$cal_ses_betamntd(runs = 1000, abundance.weighted = TRUE, use_iCAMP = TRUE, iCAMP_tempdir = tmp.3,nworker = 24 , null.model = "frequency")
saveRDS(microeco.region.frequency.nullmodel, "microeco.region.frequency.nullmodel.rds")

microeco.region.frequency.nullmodel$cal_rcbray(runs = 1000)
saveRDS(microeco.region.frequency.nullmodel, "microeco.region.frequency.nullmodel.rds")

microeco.region.frequency.nullmodel$cal_process(use_betamntd = TRUE, group = "Region")
saveRDS(microeco.region.frequency.nullmodel, "microeco.region.frequency.nullmodel.rds")

## overall dataset for= null model - frequency
microeco.frequency.nullmodel <- readRDS("microeco.frequency.nullmodel.rds")
photosynthetic.frequency.nullmodel <- readRDS("photosynthetic.frequency.nullmodel.rds")
chemolithoautotroph.frequency.nullmodel <- readRDS("chemolithoautotroph.frequency.nullmodel.rds")
heterotroph.frequency.nullmodel <- readRDS("heterotroph.frequency.nullmodel.rds")

## inbuilt region binning (overall dataset) for all null models
microeco.region.frequency.nullmodel <- readRDS("microeco.region.frequency.nullmodel.rds")

## inbuilt location binning (overall dataset) for all null models
microeco.location.frequency.nullmodel <- readRDS("microeco.location.frequency.nullmodel.rds")

## inbuilt site binning (overall dataset) for all null models
microeco.site.frequency.nullmodel <- readRDS("microeco.site.frequency.nullmodel.rds")


### Get all process in onne dataframe for frequency #####

overall.proc.frequency <- as.data.frame(microeco.frequency.nullmodel$res_process)
overall.proc.frequency$scale <- "overall"
overall.proc.frequency$null.model <- "frequency"

photo.proc.frequency <- as.data.frame(photosynthetic.frequency.nullmodel$res_process)
photo.proc.frequency$scale <- "photosynthetic"
photo.proc.frequency$null.model <- "frequency"

chemo.proc.frequency <- as.data.frame(chemolithoautotroph.frequency.nullmodel$res_process)
chemo.proc.frequency$scale <- "chemoautotrophic"
chemo.proc.frequency$null.model <- "frequency"

hetero.proc.frequency <- as.data.frame(heterotroph.frequency.nullmodel$res_process)
hetero.proc.frequency$scale <- "heterotrophic"
hetero.proc.frequency$null.model <- "frequency"

microeco.region.proc.frequency <- as.data.frame(microeco.region.frequency.nullmodel$res_process)

region.inbuilt.frequency.process.avg <- microeco.region.proc.frequency %>% 
  group_by(process) %>%
  summarise(across(percentage, mean, na.rm = TRUE))

region.inbuilt.frequency.process.avg$scale <- "region"
region.inbuilt.frequency.process.avg$null.model <- "frequency"


microeco.location.proc.frequency <- as.data.frame(microeco.location.frequency.nullmodel$res_process)

location.inbuilt.frequency.process.avg <- microeco.location.proc.frequency %>% 
  group_by(process) %>%
  summarise(across(percentage, mean, na.rm = TRUE))

location.inbuilt.frequency.process.avg$scale <- "location"
location.inbuilt.frequency.process.avg$null.model <- "frequency"


microeco.site.proc.frequency <- as.data.frame(microeco.site.frequency.nullmodel$res_process)

site.inbuilt.frequency.process.avg <- microeco.site.proc.frequency %>% 
  group_by(process) %>%
  summarise(across(percentage, mean, na.rm = TRUE))

site.inbuilt.frequency.process.avg$scale <- "site"
site.inbuilt.frequency.process.avg$null.model <- "frequency"

frequency.inbuilt.proc <- rbind(overall.proc.frequency, photo.proc.frequency, chemo.proc.frequency , hetero.proc.frequency, region.inbuilt.frequency.process.avg, location.inbuilt.frequency.process.avg, site.inbuilt.frequency.process.avg) 

frequency.nullmodel.facet.bubbleplot <- all.nullmodel %>% filter(null.model == "frequency") %>% mutate(scale = factor(scale, levels=c("overall","region","location","site","photosynthetic","chemoautotrophic","heterotrophic")))  %>% ggplot( aes(x = percentage, y = reorder(process,percentage), color = process)) + geom_point(alpha = 0.75, size = 9)  + labs( x= "% of process", y = "process") +  theme_bw() + theme( axis.title = element_text(colour = "black", face = "bold", size = 12) , axis.text.x = element_text(colour = "black", size = 15, face = "bold", angle = 0, hjust = 0.5),   axis.text.y = element_text(colour = "black", face = "bold", size = 15), legend.title = element_text(size = 9, face = "bold"), panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2), legend.position = "bottom", legend.text = element_text(size = 9, face = "bold") ,strip.text.x = element_text(size = 18, colour = "black", face = "bold" ))  + ggtitle("Processes influencing microbial community structure across various scales and 'frequency' null model")  +  scale_color_manual(values = c("#AA4371", "#E7B800", "#FC4E07","forestgreen", "purple")) + tidytext::scale_y_reordered() + facet_wrap(. ~ scale, scales = "free_y")

ggsave("frequency.nullmodel.facet.bubbleplot.pdf", height = 16, width = 25)
ggsave("frequency.nullmodel.facet.bubbleplot.png", height = 14, width = 25)



### END OF SECTION ####

###############################################################################################################

