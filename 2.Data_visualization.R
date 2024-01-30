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
library(microViz)
library(DECIPHER)
library("metacoder")
library(microeco)
library(magrittr)
library(ggcor)

### Load R environment saved at the end of 1.RawSeqs_to_Phyloseq_short to get objects unsaved as .rds
## if not then there is no need to reload the R environment
load(file = "TMFZ_Renv.RData")

### Read all saved phyloseqs ####

### Rarefied
rarefied.min.3 <- readRDS("/path_to_output_directory/TMFZ/03_tabletax/rarefied.min.3")

### Final working dataset (Rarefied and filtered)
rarefied.min.prop.exclude.1.3 <- readRDS("/path_to_output_directory/TMFZ/03_tabletax/rarefied.min.prop.exclude.1.3")
rarefied.min.int.exclude.1.3 <- readRDS("/path_to_output_directory/TMFZ/03_tabletax/rarefied.min.int.exclude.1.3")


####################################################################################################################

### Fig. S1. Sampling curves for 16S rRNA gene sequencing, shown for unrarefied data and for the sequencing data rarefied to the minimum sequencing depth (62,868 reads per sample ####

######################################################################################################################

## get otu table as matrix for rarefaction curves

## Function rarecurve draws a rarefaction curve for each row of the input data. The rarefaction curves are evaluated using the interval of step sample sizes, always including 1 and total sample size. 

#If sample is specified, a vertical line is drawn at sample with horizontal lines for the rarefied species richnesses.

### axis description: sample size - number od raw reads and Species - ASVs

ps.clean.3 <- readRDS("ps.clean.3.rds")

ASV.wide.3 <- as.matrix(as.data.frame(ps.clean.3@otu_table))



rarefied.min.3 <- readRDS("rarefied.min.3.rds")

rarefied.min.wide.3 <- as.matrix(rarefied.min.3@otu_table)



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
### Fig. S3. Distance decay plots of phylogenetic versus geographic distance for ecological groups of bacteria (N = 395) ####

##############################################################################################################

gps <- (metadata %>% select(Longitude,Latittude))

rarefied.unifrac.dist <- as.matrix(phyloseq::distance(com.rarefied.min.int.exclude.1.3, method = "unifrac"))

rarefied.unifrac.chloroflexi.dist <- as.matrix(phyloseq::distance(rarefied.chloroflexi.prop.exclude.1.3, method = "unifrac"))

rarefied.unifrac.cyanobacteria.dist <- as.matrix(phyloseq::distance(rarefied.cyanobacteria.prop.exclude.1.3, method = "unifrac"))

rarefied.unifrac.proteobacteria.otherphoto.dist <- as.matrix(phyloseq::distance(rarefied.proteobacteria.otherphoto.prop.exclude.1.3, method = "unifrac"))

rarefied.unifrac.photosynthetic.dist <- as.matrix(phyloseq::distance(rarefied.photosynthetic.prop.exclude.1.3, method = "unifrac"))

rarefied.unifrac.chemolithotrophs.dist <- as.matrix(phyloseq::distance(rarefied.chemolithotrophs.prop.exclude.1.3, method = "unifrac"))

rarefied.unifrac.only.heterotrophs.dist <- as.matrix(phyloseq::distance(rarefied.only.heterotrophs.prop.exclude.1.3, method = "unifrac"))

test <- as.data.frame(rarefied.only.heterotrophs.prop.exclude.1.3@otu_table)

### get gps co-ordinates

## all taxa
metadata.dd <- as.matrix(as.data.frame(com.rarefied.min.int.exclude.1.3@sam_data)) 
metadata.dd <- as.data.frame(metadata.dd) %>% rownames_to_column(var = "sample_id")  %>% as_tibble() 
gps.all <- metadata.dd %>% select(sample_id,Longitude,Latittude) %>% column_to_rownames(var = "sample_id")
gps.all[1:ncol(gps.all)] <- lapply(gps.all[1:ncol(gps.all)], as.numeric)

#gps for chloroflexi
metadata.chloroflexi.dd <- as.matrix(as.data.frame(rarefied.chloroflexi.prop.exclude.1.3@sam_data)) 
metadata.chloroflexi.dd <- as.data.frame(metadata.chloroflexi.dd) %>% rownames_to_column(var = "sample_id")  %>% as_tibble() 
gps.chloroflexi <- metadata.chloroflexi.dd %>% select(sample_id,Longitude,Latittude) %>% column_to_rownames(var = "sample_id")
gps.chloroflexi[1:ncol(gps.chloroflexi)] <- lapply(gps.chloroflexi[1:ncol(gps.chloroflexi)], as.numeric)


identical(gps.all,gps.chloroflexi)
### all gps data is the same for all photosynthetic and non-photosynhetic groups.



#Vector of distances in km for all taxa
distgeo.all <- geodist(gps.all, measure = 'geodesic')/1000
rownames(distgeo.all) <- rownames(gps.all)
colnames(distgeo.all) <- rownames(gps.all)

citation("geodist")

#library(geosphere)
#distgeo3 <- as.matrix(distm(gps[,1:2],fun = distHaversine)[,1]/1000)

max(distgeo.all)

## melt distance matrices
ASV.bray.dist.df <- melt(as.matrix(ASV.bray.dist), varnames = c("row", "col"))

rarefied.unifrac.dist.df <- melt(as.matrix(rarefied.unifrac.dist), varnames = c("row", "col"))

rarefied.unifrac.chloroflexi.dist.df <- melt(as.matrix(rarefied.unifrac.chloroflexi.dist), varnames = c("row", "col"))

rarefied.unifrac.cyanobacteria.dist.df <- melt(as.matrix(rarefied.unifrac.cyanobacteria.dist), varnames = c("row", "col"))

rarefied.unifrac.proteobacteria.otherphoto.dist.df <- melt(as.matrix(rarefied.unifrac.proteobacteria.otherphoto.dist), varnames = c("row", "col"))

rarefied.unifrac.photosynthetic.dist.df <- melt(as.matrix(rarefied.unifrac.photosynthetic.dist), varnames = c("row", "col"))

rarefied.unifrac.chemolithotrophs.dist.df <- melt(as.matrix(rarefied.unifrac.chemolithotrophs.dist), varnames = c("row", "col"))

rarefied.unifrac.only.heterotrophs.dist.df <- melt(as.matrix(rarefied.unifrac.only.heterotrophs.dist), varnames = c("row", "col"))

distgeo.all.df <- melt(as.matrix(distgeo.all), varnames = c("row", "col"))

colnames(ASV.bray.dist.df)[3] <- "Braydist"
colnames(rarefied.unifrac.dist.df)[3] <- "unweighted.unifrac.alltaxa"
colnames(rarefied.unifrac.chloroflexi.dist.df)[3] <- "unweighted.unifrac.chloroflexi"
colnames(rarefied.unifrac.cyanobacteria.dist.df)[3] <- "unweighted.unifrac.cyanobacteria"
colnames(rarefied.unifrac.proteobacteria.otherphoto.dist.df)[3] <- "unweighted.unifrac.proteobacteria.otherphoto"
colnames(rarefied.unifrac.photosynthetic.dist.df)[3] <- "unweighted.unifrac.all.photosynthetic"
colnames(rarefied.unifrac.chemolithotrophs.dist.df)[3] <- "unweighted.unifrac.chemolithotrophic"
colnames(rarefied.unifrac.only.heterotrophs.dist.df)[3] <- "unweighted.unifrac.heterotrophic"
colnames(distgeo.all.df)[3] <- "Kmdist"

### MERGE distance matrix and gps df
mergedist.bray <- merge(ASV.bray.dist.df,distgeo.all.df)
mergedist.unifrac <- merge(rarefied.unifrac.dist.df,distgeo.all.df)
mergedist.unifrac.chloroflexi <- merge(rarefied.unifrac.chloroflexi.dist.df,distgeo.all.df)
mergedist.unifrac.cyanobacteria <- merge(rarefied.unifrac.cyanobacteria.dist.df,distgeo.all.df)
mergedist.unifrac.proteobacteria.otherphoto <- merge(rarefied.unifrac.proteobacteria.otherphoto.dist.df,distgeo.all.df)
mergedist.unifrac.photosynthetic <- merge(rarefied.unifrac.photosynthetic.dist.df,distgeo.all.df)
mergedist.unifrac.chemolithotrophs <- merge(rarefied.unifrac.chemolithotrophs.dist.df,distgeo.all.df)
mergedist.unifrac.only.heterotrophs <- merge(rarefied.unifrac.only.heterotrophs.dist.df,distgeo.all.df)

### plot distance deacays ####

distdecay.bray <- mergedist.bray %>% ggplot(aes(x=Kmdist, y=Braydist)) + geom_point(size = 1, alpha = 0.1, position = position_dodge(0.75)) + geom_smooth(formula = y ~ x, method = "lm", se = TRUE, color = "maroon") + stat_cor(label.x = 1000, label.y = 0.5) + stat_regline_equation(label.x = 1000, label.y = 0.45)  + theme(axis.text.x =element_text(size=14, face = "bold", angle = 0 , hjust = 0.5),axis.text.y =element_text(size=12, face = "bold"),axis.title.y=element_text(size = 12, face = "bold"),axis.title.x = element_text(size = 12, face = "bold")) +  theme(legend.position = "none") +labs(title="Distance deecay of Bray-curtis distance against pairwise geographical distance in Km ", y = "Bray-Curtis dissimilarity", x = "Geographical distance (Km)")+ theme(strip.text.y = element_text(angle = 0, face = "bold", size = 12)) + theme(plot.title = element_text(hjust = 0.5, size = 11, face = "bold")) + ggeasy::easy_center_title() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) 

ggsave("distdecay.bray.pdf")

model.bray <- lm(Braydist ~ Kmdist , data= mergedist.bray)
summary(model.bray) 

distdecay.unifrac <- mergedist.unifrac %>% ggplot(aes(x=Kmdist, y=unweighted.unifrac.alltaxa)) + geom_point(size = 1, alpha = 0.1, position = position_dodge(0.75)) + geom_smooth(formula = y ~ x, method = "lm", se = TRUE, color = "maroon") + stat_cor(label.x = 1000, label.y = 0.3) + stat_regline_equation(label.x = 1000, label.y = 0.25)  + theme(axis.text.x =element_text(size=14, face = "bold", angle = 0 , hjust = 0.5),axis.text.y =element_text(size=12, face = "bold"),axis.title.y=element_text(size = 12, face = "bold"),axis.title.x = element_text(size = 12, face = "bold")) +  theme(legend.position = "none") +labs(title="Distance deecay of Unweighted UniFrac distance (all taxa) against pairwise geographical distance in Km ", y = "Unweighted UniFrac dissimilarity", x = "Geographical distance (Km)")+ theme(strip.text.y = element_text(angle = 0, face = "bold", size = 12)) + theme(plot.title = element_text(hjust = 0.5, size = 11, face = "bold")) + ggeasy::easy_center_title() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) 

ggsave("distdecay.unifrac.pdf")
ggsave("distdecay.unifrac.png")

model.unifrac <- lm(unweighted.unifrac.alltaxa ~ Kmdist , data= mergedist.unifrac)
summary(model.unifrac)
anova(model.unifrac)

distdecay.unifrac.chloroflexi <- mergedist.unifrac.chloroflexi %>% ggplot(aes(x=Kmdist, y=unweighted.unifrac.chloroflexi)) + geom_point(size = 1, alpha = 0.1, position = position_dodge(0.75)) + geom_smooth(formula = y ~ x, method = "lm", se = TRUE, color = "maroon") + stat_cor(label.x = 1000, label.y = 0.7) + stat_regline_equation(label.x = 1000, label.y = 0.65)  + theme(axis.text.x =element_text(size=14, face = "bold", angle = 0 , hjust = 0.5),axis.text.y =element_text(size=12, face = "bold"),axis.title.y=element_text(size = 12, face = "bold"),axis.title.x = element_text(size = 12, face = "bold")) +  theme(legend.position = "none") +labs(title="Distance deecay of Unweighted UniFrac distance (chloroflexi) against pairwise geographical distance in Km ", y = "Unweighted UniFrac dissimilarity", x = "Geographical distance (Km)")+ theme(strip.text.y = element_text(angle = 0, face = "bold", size = 12)) + theme(plot.title = element_text(hjust = 0.5, size = 11, face = "bold")) + ggeasy::easy_center_title() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) 

ggsave("distdecay.unifrac.chloroflexi.pdf")
ggsave("distdecay.unifrac.chloroflexi.png")

model.unifrac.chloroflexi <- lm(unweighted.unifrac.chloroflexi ~ Kmdist , data= mergedist.unifrac.chloroflexi)
summary(model.unifrac.chloroflexi)
anova(model.unifrac.chloroflexi)


distdecay.unifrac.cyanobacteria <- mergedist.unifrac.cyanobacteria %>% ggplot(aes(x=Kmdist, y=unweighted.unifrac.cyanobacteria)) + geom_point(size = 1, alpha = 0.1, position = position_dodge(0.75)) + geom_smooth(formula = y ~ x, method = "lm", se = TRUE, color = "maroon") + stat_cor(label.x = 1000, label.y = 0.15) + stat_regline_equation(label.x = 1000, label.y = 0.10)  + theme(axis.text.x =element_text(size=14, face = "bold", angle = 0 , hjust = 0.5),axis.text.y =element_text(size=12, face = "bold"),axis.title.y=element_text(size = 12, face = "bold"),axis.title.x = element_text(size = 12, face = "bold")) +  theme(legend.position = "none") +labs(title="Distance deecay of Unweighted UniFrac distance (cyanobacteria) against pairwise geographical distance in Km ", y = "Unweighted UniFrac dissimilarity", x = "Geographical distance (Km)")+ theme(strip.text.y = element_text(angle = 0, face = "bold", size = 12)) + theme(plot.title = element_text(hjust = 0.5, size = 11, face = "bold")) + ggeasy::easy_center_title() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) 

ggsave("distdecay.unifrac.cyanobacteria.pdf")
ggsave("distdecay.unifrac.cyanobacteria.png")

model.unifrac.cyanobacteria <- lm(unweighted.unifrac.cyanobacteria ~ Kmdist , data= mergedist.unifrac.cyanobacteria)
summary(model.unifrac.cyanobacteria)
anova(model.unifrac.cyanobacteria)

distdecay.unifrac.proteobacteria.otherphoto <- mergedist.unifrac.proteobacteria.otherphoto %>% ggplot(aes(x=Kmdist, y=unweighted.unifrac.proteobacteria.otherphoto)) + geom_point(size = 1, alpha = 0.1, position = position_dodge(0.75)) + geom_smooth(formula = y ~ x, method = "lm", se = TRUE, color = "maroon") + stat_cor(label.x = 1000, label.y = 0.78) + stat_regline_equation(label.x = 1000, label.y = 0.72)  + theme(axis.text.x =element_text(size=14, face = "bold", angle = 0 , hjust = 0.5),axis.text.y =element_text(size=12, face = "bold"),axis.title.y=element_text(size = 12, face = "bold"),axis.title.x = element_text(size = 12, face = "bold")) +  theme(legend.position = "none") +labs(title="Distance deecay of Unweighted UniFrac distance (proteobacteria and other photosynthetic taxa) against pairwise geographical distance in Km ", y = "Unweighted UniFrac dissimilarity", x = "Geographical distance (Km)")+ theme(strip.text.y = element_text(angle = 0, face = "bold", size = 12)) + theme(plot.title = element_text(hjust = 0.5, size = 11, face = "bold")) + ggeasy::easy_center_title() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) 

ggsave("distdecay.unifrac.proteobacteria.otherphoto.pdf")
ggsave("distdecay.unifrac.proteobacteria.otherphoto.png")

model.unifrac.proteobacteria.otherphoto <- lm(unweighted.unifrac.proteobacteria.otherphoto ~ Kmdist , data= mergedist.unifrac.proteobacteria.otherphoto)
summary(model.unifrac.proteobacteria.otherphoto)
anova(model.unifrac.proteobacteria.otherphoto)

distdecay.unifrac.photosynthetic <- mergedist.unifrac.photosynthetic %>% ggplot(aes(x=Kmdist, y=unweighted.unifrac.all.photosynthetic)) + geom_point(size = 1, alpha = 0.1, position = position_dodge(0.75)) + geom_smooth(formula = y ~ x, method = "lm", se = TRUE, color = "maroon") + stat_cor(label.x = 1000, label.y = 0.12) + stat_regline_equation(label.x = 1000, label.y = 0.09)  + theme(axis.text.x =element_text(size=14, face = "bold", angle = 0 , hjust = 0.5),axis.text.y =element_text(size=12, face = "bold"),axis.title.y=element_text(size = 12, face = "bold"),axis.title.x = element_text(size = 12, face = "bold")) +  theme(legend.position = "none") +labs(title="Distance deecay of Unweighted UniFrac distance (all photosynthetic taxa) against pairwise geographical distance in Km ", y = "Unweighted UniFrac dissimilarity", x = "Geographical distance (Km)")+ theme(strip.text.y = element_text(angle = 0, face = "bold", size = 12)) + theme(plot.title = element_text(hjust = 0.5, size = 11, face = "bold")) + ggeasy::easy_center_title() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) 

ggsave("distdecay.unifrac.photosynthetic.pdf")
ggsave("distdecay.unifrac.photosynthetic.png")


model.unifrac.photosynthetic <- lm(unweighted.unifrac.all.photosynthetic ~ Kmdist , data= mergedist.unifrac.photosynthetic)
summary(model.unifrac.photosynthetic)
anova(model.unifrac.photosynthetic)


distdecay.unifrac.chemolithotrophs <- mergedist.unifrac.chemolithotrophs %>% ggplot(aes(x=Kmdist, y=unweighted.unifrac.chemolithotrophic)) + geom_point(size = 1, alpha = 0.1, position = position_dodge(0.75)) + geom_smooth(formula = y ~ x, method = "lm", se = TRUE, color = "maroon") + stat_cor(label.x = 1200, label.y = 0.30) + stat_regline_equation(label.x = 1200, label.y = 0.25)  + theme(axis.text.x =element_text(size=14, face = "bold", angle = 0 , hjust = 0.5),axis.text.y =element_text(size=12, face = "bold"),axis.title.y=element_text(size = 12, face = "bold"),axis.title.x = element_text(size = 12, face = "bold")) +  theme(legend.position = "none") +labs(title="Distance deecay of Unweighted UniFrac distance (Chemolithotrophs) against pairwise geographical distance in Km ", y = "Unweighted UniFrac dissimilarity", x = "Geographical distance (Km)")+ theme(strip.text.y = element_text(angle = 0, face = "bold", size = 12)) + theme(plot.title = element_text(hjust = 0.5, size = 11, face = "bold")) + ggeasy::easy_center_title() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) 

ggsave("distdecay.unifrac.chemolithotrophs.pdf")
ggsave("distdecay.unifrac.chemolithotrophs.png")

model.unifrac.chemolithotrophs <- lm(unweighted.unifrac.chemolithotrophic ~ Kmdist , data= mergedist.unifrac.chemolithotrophs)
summary(model.unifrac.chemolithotrophs)
anova(model.unifrac.chemolithotrophs)

distdecay.unifrac.only.heterotrophs <- mergedist.unifrac.only.heterotrophs %>% ggplot(aes(x=Kmdist, y=unweighted.unifrac.heterotrophic)) + geom_point(size = 1, alpha = 0.1, position = position_dodge(0.75)) + geom_smooth(formula = y ~ x, method = "lm", se = TRUE, color = "maroon") + stat_cor(label.x = 1000, label.y = 0.30) + stat_regline_equation(label.x = 1000, label.y = 0.25)  + theme(axis.text.x =element_text(size=14, face = "bold", angle = 0 , hjust = 0.5),axis.text.y =element_text(size=12, face = "bold"),axis.title.y=element_text(size = 12, face = "bold"),axis.title.x = element_text(size = 12, face = "bold")) +  theme(legend.position = "none") +labs(title="Distance deecay of Unweighted UniFrac distance (only heterotrophs) against pairwise geographical distance in Km ", y = "Unweighted UniFrac dissimilarity", x = "Geographical distance (Km)")+ theme(strip.text.y = element_text(angle = 0, face = "bold", size = 12)) + theme(plot.title = element_text(hjust = 0.5, size = 11, face = "bold")) + ggeasy::easy_center_title() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) 

ggsave("distdecay.unifrac.only.heterotrophs.pdf")
ggsave("distdecay.unifrac.only.heterotrophs.png")

model.unifrac.only.heterotrophs <- lm(unweighted.unifrac.heterotrophic ~ Kmdist , data= mergedist.unifrac.only.heterotrophs)
summary(model.unifrac.only.heterotrophs)
anova(model.unifrac.only.heterotrophs)


####################################################################################################################

### Fig. 2A. Beta diversity patterns for hot spring photosynthetic biofilms indicated statistically supported biogeographic regions. ####
### Beta diversity was estimated using Principle Coordinate Analysis of unweighted UniFRAC distances

######################################################################################################################

### Set color and shape legends

# to get variable names from phyloseq metadata
sample_variables(com.rarefied.min.prop.exclude.1.3)
get_variable(com.rarefied.min.prop.exclude.1.3,"Temp.adj...C.")
get_variable(com.rarefied.min.prop.exclude.1.3,"Location.Code")

## shape = Location

#AP BA HS KJ LA LN MK PB PP PT RB RN SE SW US

loc.shape.names = c('AP', 'BA', 'HS', 'KJ', 'LA', 'LN', 'MK', 'PB', 'PP', 'PT', 'RB', 'RN', 'SE', 'SW', 'US')
loc.shape <- 0:(length(loc.shape.names)-1)
names(loc.shape) <- loc.shape.names
loc.shape["Taxa"] <- 16
loc.shape

## Color = Temperature

temp.col = c('38' = "mediumorchid1", '40' = "mediumpurple1", '42' = "blueviolet",'44' = "mediumpurple4", '45' = "lightskyblue2",  '46' = "deepskyblue" , '47'= "royalblue1",  '48' = "mediumblue", '49' = "cyan", '50' = "aquamarine2" , '51' = "springgreen2",  '53' = "yellowgreen", '55' = "forestgreen", '57' = "khaki3", '58' = "goldenrod1", '61' = "sienna2", '63' = "firebrick3", '66' = "darkred", 'Taxa' = "black")

temp.col

##use temp breaks to specify temp legend first and then location
temp.breaks = c("38", "40", "42", "44" , "45", "46", "47", "48" , "49", "50", "51", "53", "55", "57", "58", "61" , "63" , "66","Taxa")

region.col = c('North.Thailand' = "slateblue", 'Central.Thailand' = "turquoise3", 'South.Thailand' = "gold", 'North.Malaysia' = "chocolate1", 'South.Malaysia' = "lawngreen", 'Singapore'="palevioletred3")

#### Calcluate Unweighted UniFRAC ordiation ####

## "unifrac" :Original (unweighted) UniFrac distance, UniFrac
## "wunifrac" :weighted-UniFrac distance, UniFrac
## "dpcoa": sample-wise distance used in Double Principle Coordinate Analysis, DPCoA

## MDS/PCoA: Performs principal coordinate analysis (also called principle coordinate decomposition, multidimensional scaling (MDS), or classical scaling) of a distance matrix (Gower 1966), including two correction methods for negative eigenvalues. See pcoa for further details.
rare.unweighted.unifrac.PCoA = ordinate(com.rarefied.min.prop.exclude.1.3, method="PCoA", distance="unifrac", weighted=FALSE)

rare.unweighted.int.unifrac.PCoA = ordinate(com.rarefied.min.int.exclude.1.3, method="PCoA", distance="unifrac", weighted=FALSE)

##Plot
rarefied.PCoA.UnweightedUniFrac.phyloseq.regionellipse.samples <- plot_ordination(com.rarefied.min.prop.exclude.1.3, rare.unweighted.unifrac.PCoA, type="samples", color = "Temp.adj...C.", shape = "Location.Code", title="PCoA on unweighted UniFrac distance") +   geom_point(size=7, position="jitter") + stat_ellipse(aes(fill = Region, group = Region), linetype = 0 ,type = "t",  level = 0.99, alpha = .2, geom = "polygon")   + clean_background  + scale_shape_manual(values= loc.shape)  + scale_color_manual(values = temp.col) + scale_fill_manual(values = region.col,breaks = c('North.Thailand','Central.Thailand','South.Thailand','North.Malaysia','South.Malaysia','Singapore')) + theme(legend.title = element_text(size = 15.5, face = "bold"), legend.text = element_text(size = 14),legend.position = "bottom", legend.direction = "horizontal",legend.box="vertical", legend.margin=margin()) + guides(colour = guide_legend(nrow=1,byrow=TRUE, override.aes = list(size = 7)), shape = guide_legend(nrow=1,byrow=TRUE,override.aes = list(size = 7)), fill = guide_legend(nrow=1,byrow=TRUE,override.aes = list(size = 7)) )

#scale_fill_manual(values = my_values, breaks = nmds.breaks) 

rarefied.PCoA.UnweightedUniFrac.phyloseq.regionellipse.samples$layers <- rarefied.PCoA.UnweightedUniFrac.phyloseq.regionellipse.samples$layers[-1]

pdf("rarefied.PCoA.UnweightedUniFrac.phyloseq.regionellipse.samples.pdf", height = 12, width = 17)
rarefied.PCoA.UnweightedUniFrac.phyloseq.regionellipse.samples
dev.off()

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
            #node_color_range = c("skyblue", "yellowgreen","bisque", "coral", "firebrick4"),
            #edge_label = n_obs,
            initial_layout = "re",
            layout = "da",
            title = "ASV diversity",
            node_color_axis_label = "Read abundance",
            node_size_axis_label = "Number of OTUs",
            output_file = "Species.heattree.supplementary.labels.pdf")


#Comparing any number of treatments/groups
#For each taxon, a Wilcoxon Rank Sum test was used to test for differences between the median abundances of samples in each treatment. We can use this information to create what we call a differential heat tree,

otu_tax_taxmap$data$diff_table <- compare_groups(otu_tax_taxmap, "tax_abund",
                                                 cols = metadata$Sample_id, # What columns of sample data to use
                                                 groups = metadata$Region) # What category each sample is assigned to

set.seed(1)
heat_tree_matrix(otu_tax_taxmap,
                 data = "diff_table",
                 node_size = n_obs, # n_obs is a function that calculates, in this case, the number of OTUs per taxon
                 node_label = taxon_names,
                 node_color = log2_median_ratio, # A column from `obj$data$diff_table`
                 node_color_range = diverging_palette(), # The built-in palette for diverging data
                 node_color_trans = "linear", # The default is scaled by circle area
                 node_color_interval = c(-3, 3), # The range of `log2_median_ratio` to display
                 edge_color_interval = c(-3, 3), # The range of `log2_median_ratio` to display
                 node_size_axis_label = "Number of OTUs",
                 node_color_axis_label = "Log2 ratio median proportions",
                 layout = "davidson-harel", # The primary layout algorithm
                 initial_layout = "reingold-tilford", # The layout algorithm that initializes node locations
                 output_file = "differential_heat_tree.pdf") # Saves the plot as a pdf file

####################################################################################################################

### Fig. S2.Alpha diversity metrics for the 40 hot spring locations sampled in this study.Samples are displayed in order of ascending temperature. ####

######################################################################################################################

rarefied.min.metadata.3 <- as.data.frame(new.rarefied.min.3@sam_data)

rarefied.min.metadata.3$Temp.adj...C. <- as.factor(rarefied.min.metadata.3$Temp.adj...C.)
rarefied.min.metadata.3$Longitude <- as.factor(rarefied.min.metadata.3$Longitude)
rarefied.min.metadata.3$Latittude <- as.factor(rarefied.min.metadata.3$Latittude)
rarefied.min.metadata.3$EC..mS. <- as.factor(rarefied.min.metadata.3$EC..mS.)
rarefied.min.metadata.3$pH <- as.factor(rarefied.min.metadata.3$pH)
rarefied.min.metadata.3$Carbonate..ppm. <- as.factor(rarefied.min.metadata.3$Carbonate..ppm.)

names(rarefied.min.metadata.3)[names(rarefied.min.metadata.3) == "Temp.adj...C."] <- "Temp.adj (°C)"
names(rarefied.min.metadata.3)[names(rarefied.min.metadata.3) == "Temp...C."] <- "Temp.actual (°C)"

Site.rarefied.min.3 <- rarefied.min.metadata.3$Site
Sample_id.rarefied.min.3 <- rarefied.min.metadata.3$Sample_id
Temp.rarefied.min.3 <- rarefied.min.metadata.3$`Temp.adj (°C)`
Long.rarefied.min.3 <- rarefied.min.metadata.3$Longitude
Lat.rarefied.min.3 <- rarefied.min.metadata.3$Latittude
Cond.rarefied.min.3 <- rarefied.min.metadata.3$EC..mS. 
pH.rarefied.min.3 <- rarefied.min.metadata.3$pH 
carb.rarefied.min.3 <- rarefied.min.metadata.3$Carbonate..ppm.


rarefied.min.wide.3 <- as.matrix(as.data.frame(new.rarefied.min.3@otu_table))


### calculte Shannon Index (H)
rarefied.min.shannon.3 <- diversity(rarefied.min.wide.3, index="shannon")

## Pielou's eveness (J) for Shannon index
rarefied.min.pielou.3 <- diversity(rarefied.min.wide.3,index="shannon")/log(specnumber(rarefied.min.wide.3))

## Simpson's index (D,λ) and Gini-Simpson(GS)=1-λ
rarefied.min.simpson.3 <- diversity(rarefied.min.wide.3, index="simpson")

## Chao estimator (specpool and estimateR)
rarefied.min.chao.3 <- as.data.frame(estimateR(rarefied.min.wide.3))

diversity.indices.rarefied.min.3 <- data.frame(Site.rarefied.min.3,Sample_id.rarefied.min.3,Temp.rarefied.min.3, Lat.rarefied.min.3, Cond.rarefied.min.3, pH.rarefied.min.3,carb.rarefied.min.3)
diversity.indices.rarefied.min.3$rarefied.min.shannon.3 <- rarefied.min.shannon.3
diversity.indices.rarefied.min.3$rarefied.min.pielou.3 <- rarefied.min.pielou.3
diversity.indices.rarefied.min.3$rarefied.min.simpson.3 <- rarefied.min.simpson.3
diversity.indices.rarefied.min.3$rarefied.min.chao.3 <- t(rarefied.min.chao.3["S.chao1", ])

### imp to change x-axis variable to character to later change them as numeric after getting long table

diversity.indices.rarefied.min.3$Temp.rarefied.min.3 <- as.character(diversity.indices.rarefied.min.3$Temp.rarefied.min.3)

diversity.indices.rarefied.min.3$Lat.rarefied.min.3 <- as.character(diversity.indices.rarefied.min.3$Lat.rarefied.min.3)

diversity.indices.rarefied.min.3$Cond.rarefied.min.3 <- as.character(diversity.indices.rarefied.min.3$Cond.rarefied.min.3)

diversity.indices.rarefied.min.3$pH.rarefied.min.3 <- as.character(diversity.indices.rarefied.min.3$pH.rarefied.min.3)

diversity.indices.rarefied.min.3$carb.rarefied.min.3 <- as.character(diversity.indices.rarefied.min.3$carb.rarefied.min.3)

rownames(diversity.indices.rarefied.min.3) <- NULL

setDT(diversity.indices.rarefied.min.3)

diversity.indices.rarefied.min.long.3  <- melt(diversity.indices.rarefied.min.3[1:nrow(diversity.indices.rarefied.min.3),])

new_var_names_rarefied_min_3 <- c(
  'rarefied.min.chao.3'="Chao1 estimator",
  'rarefied.min.shannon.3'="Shannon's index",
  'rarefied.min.pielou.3'="Pielou's evenness",
  'rarefied.min.simpson.3'="Gini-Simpson's index"
)

## Temperature
diversity.indices.rarefied.min.long.3$Temp.rarefied.min.3 <- as.numeric(diversity.indices.rarefied.min.long.3$Temp.rarefied.min.3)

diversity.indices.rarefied.min.long.3 <- arrange(diversity.indices.rarefied.min.long.3, (Temp.rarefied.min.3), desc(value),)

## Latitude
diversity.indices.rarefied.min.long.3$Lat.rarefied.min.3 <- as.numeric(diversity.indices.rarefied.min.long.3$Lat.rarefied.min.3)

diversity.indices.rarefied.min.long.3[,4] <- round(diversity.indices.rarefied.min.long.3[,4], digits = 4)

diversity.indices.rarefied.min.long.3 <- arrange(diversity.indices.rarefied.min.long.3, desc(Lat.rarefied.min.3), (value),)


## Conductivity
diversity.indices.rarefied.min.long.3$Cond.rarefied.min.3 <- as.numeric(diversity.indices.rarefied.min.long.3$Cond.rarefied.min.3)

diversity.indices.rarefied.min.long.3[,5] <- round(diversity.indices.rarefied.min.long.3[,5], digits = 2)

diversity.indices.rarefied.min.long.3 <- arrange(diversity.indices.rarefied.min.long.3, desc(Cond.rarefied.min.3), (value),)

## pH
diversity.indices.rarefied.min.long.3$pH.rarefied.min.3 <- as.numeric(diversity.indices.rarefied.min.long.3$pH.rarefied.min.3)

diversity.indices.rarefied.min.long.3[,6] <- round(diversity.indices.rarefied.min.long.3[,6], digits = 2)

diversity.indices.rarefied.min.long.3 <- arrange(diversity.indices.rarefied.min.long.3, desc(pH.rarefied.min.3), (value),)

## Carbonate
diversity.indices.rarefied.min.long.3$carb.rarefied.min.3 <- as.numeric(diversity.indices.rarefied.min.long.3$carb.rarefied.min.3)

diversity.indices.rarefied.min.long.3 <- arrange(diversity.indices.rarefied.min.long.3, desc(carb.rarefied.min.3), (value),)


## plot

rarefied.min.alpha.div.indices.loc.notempfacet.3 <- diversity.indices.rarefied.min.long.3  %>% filter(!variable %in% "rarefied.min.no.3") %>% mutate(variable = factor(variable, levels=c("rarefied.min.chao.3","rarefied.min.shannon.3","rarefied.min.simpson.3","rarefied.min.pielou.3"))) %>% mutate(Site.rarefied.min.3 = factor(Site.rarefied.min.3, levels= c("PB_38" , "RB_40" , "RN_42" , "HS_42" , "MK_42" , "MK_44" , "RN_45", "HS_45" , "RB_45", "SE_45", "SW_46" , "PB_46"  , "LN_46", "AP_46" , "LA_46" , "KJ_47" , "PP_48" , "SE_49", "PT_50", "MK_50"  , "BA_50"  , "RB_50" , "RN_51" , "HS_51" , "LN_51" , "SW_51" , "HS_53" , "PP_53" , "RB_55", "PT_55" , "SW_55" , "MK_55" , "US_57" , "PP_58" , "RN_61", "SW_61", "LN_63","US_63" , "PT_63" , "LN_66"))) %>% drop_na() %>% ggplot(aes(x=Site.rarefied.min.3, y=value, fill=variable)) + geom_boxplot(outlier.colour = "red") + stat_summary(fun=mean, geom="point", color="black", size=2, shape=23, position=position_dodge(0.75)) +  theme(axis.text.x =element_text(size=10, face = "bold", angle = 45 , hjust = 1),axis.text.y =element_text(size=12, face = "bold"),axis.title.y=element_text(size = 12, face = "bold"),axis.title.x = element_blank()) +  theme(legend.position = "none") +labs(title="Diversity indices of rarefied ASVs with error model 3", y = "Diversity index") + facet_grid(vars(variable), scales = "free",labeller = as_labeller(new_var_names_rarefied_min_3)) + theme(strip.text.y = element_text(angle = 0, face = "bold", size = 12)) + theme(plot.title = element_text(hjust = 0.5, size = 13, face = "bold")) + ggeasy::easy_center_title() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + geom_smooth(formula = y ~ x, method = "lm", se = TRUE, color = "grey40", aes(group = variable))  + stat_regline_equation(label.x = 1, label.y = 0.7,aes(group = variable)) + stat_cor(label.x = 1, label.y = 1, aes(group = variable)) 


#+ scale_x_discrete(labels=c("A" = "A:Lower channel at 45.7°C", "B" = "B:Mid channel at 51°C","C" = "C:Upper channel at 55.3°C", "D" = "D:Source at 61.4°C"))

rarefied.min.alpha.div.indices.loc.notempfacet.3

ggsave("rarefied.min.alpha.div.indices.loc.notempfacet.3.pdf",height = 12, width = 20)
ggsave("rarefied.min.alpha.div.indices.loc.notempfacet.3.png", height = 12, width = 20)

## to get equations clearly for chao as it cannot be seen in plot with all indices
rarefied.min.chao.loc.notempfacet.3 <- diversity.indices.rarefied.min.long.3  %>% filter(variable == "rarefied.min.chao.3")  %>% mutate(Site.rarefied.min.3 = factor(Site.rarefied.min.3, levels= c("PB_38" , "RB_40" , "RN_42" , "HS_42" , "MK_42" , "MK_44" , "RN_45", "HS_45" , "RB_45", "SE_45", "SW_46" , "PB_46"  , "LN_46", "AP_46" , "LA_46" , "KJ_47" , "PP_48" , "SE_49", "PT_50", "MK_50"  , "BA_50"  , "RB_50" , "RN_51" , "HS_51" , "LN_51" , "SW_51" , "HS_53" , "PP_53" , "RB_55", "PT_55" , "SW_55" , "MK_55" , "US_57" , "PP_58" , "RN_61", "SW_61", "LN_63","US_63" , "PT_63" , "LN_66"))) %>% drop_na() %>% ggplot(aes(x=Site.rarefied.min.3, y=value, fill=variable)) + geom_boxplot(outlier.colour = "red") + stat_summary(fun=mean, geom="point", color="black", size=2, shape=23, position=position_dodge(0.75)) +  theme(axis.text.x =element_text(size=10, face = "bold", angle = 45 , hjust = 1),axis.text.y =element_text(size=12, face = "bold"),axis.title.y=element_text(size = 12, face = "bold"),axis.title.x = element_blank()) +  theme(legend.position = "none") +labs(title="Diversity indices of rarefied ASVs with error model 3", y = "Diversity index") + facet_grid(vars(variable), scales = "free",labeller = as_labeller(new_var_names_rarefied_min_3)) + theme(strip.text.y = element_text(angle = 0, face = "bold", size = 12)) + theme(plot.title = element_text(hjust = 0.5, size = 13, face = "bold")) + ggeasy::easy_center_title() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + geom_smooth(formula = y ~ x, method = "lm", se = TRUE, color = "grey40", aes(group = variable))  + stat_regline_equation(label.x = 1, label.y =2750 ,aes(group = variable)) + stat_cor(label.x = 1, label.y = 3000, aes(group = variable)) 

rarefied.min.chao.loc.notempfacet.3

rarefied.min.alpha.div.indices.lat.3 <- diversity.indices.rarefied.min.long.3   %>% filter(!variable == "rarefied.min.no.3") %>% mutate(variable = factor(variable, levels=c("rarefied.min.chao.3","rarefied.min.shannon.3","rarefied.min.simpson.3","rarefied.min.pielou.3"))) %>% drop_na() %>% ggplot(aes(x=reorder(factor(as.numeric(Lat.rarefied.min.3)), -(as.numeric(Lat.rarefied.min.3))), y=value, fill=variable)) + geom_point(size=0.3) + stat_summary(fun=mean, geom="point", color="red", size=0.4, shape=23, position=position_dodge(0.75)) +  theme(axis.text.x =element_text(size=6, face = "bold", angle = 45 , hjust = 1),axis.text.y =element_text(size=6, face = "bold"),axis.title.y=element_blank(),axis.title.x = element_blank()) +  theme(legend.position = "none") + facet_wrap(vars(variable), scales = "free_y",labeller = as_labeller(new_var_names_rarefied_min_3),  ncol=1) + theme(strip.text.y = element_text(angle = 0, face = "bold", size = 8)) + theme(plot.title = element_blank()) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + geom_smooth(formula = y ~ x, method = "lm", se = TRUE, color = "grey40", aes(group = variable)) + stat_regline_equation(label.x = 1, label.y =0.7 ,aes(group = variable)) + stat_cor(label.x = 1, label.y = 1, aes(group = variable)) 


model <- lm(value ~  as.numeric(Lat.rarefied.min.3) , data= diversity.indices.rarefied.min.long.3)
summary(model)


ggsave("rarefied.min.alpha.div.indices.lat.3.pdf",height = 11, width = 5)
ggsave("rarefied.min.alpha.div.indices.lat.3.png", height = 11, width = 3)


rarefied.min.alpha.div.indices.pH.3 <- diversity.indices.rarefied.min.long.3  %>% filter(!variable == "rarefied.min.no.3") %>% mutate(variable = factor(variable, levels=c("rarefied.min.chao.3","rarefied.min.shannon.3","rarefied.min.simpson.3","rarefied.min.pielou.3"))) %>% drop_na() %>% ggplot(aes(x=reorder(factor(pH.rarefied.min.3), -pH.rarefied.min.3), y=value, fill=variable)) + geom_point(size=0.3) + stat_summary(fun=mean, geom="point", color="red", size=0.4, shape=23, position=position_dodge(0.75)) +  theme(axis.text.x =element_text(size=6, face = "bold", angle = 45 , hjust = 1),axis.text.y =element_text(size=6, face = "bold"),axis.title.y=element_blank(),axis.title.x = element_blank()) +  theme(legend.position = "none") + facet_wrap(vars(variable), scales = "free_y",labeller = as_labeller(new_var_names_rarefied_min_3),  ncol=1) + theme(strip.text.y = element_text(angle = 0, face = "bold", size = 8)) + theme(plot.title = element_blank()) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + geom_smooth(formula = y ~ x, method = "lm", se = TRUE, color = "grey40", aes(group = variable))  + stat_regline_equation(label.x = 1, label.y =2000 ,aes(group = variable)) + stat_cor(label.x = 1, label.y = 1000, aes(group = variable)) 

rarefied.min.alpha.div.indices.pH.3


ggsave("rarefied.min.alpha.div.indices.pH.3.pdf",height = 11, width = 3)
ggsave("rarefied.min.alpha.div.indices.pH.3.png", height = 11, width = 3)


rarefied.min.alpha.div.indices.cond.3 <- diversity.indices.rarefied.min.long.3  %>% filter(!variable == "rarefied.min.no.3") %>% mutate(variable = factor(variable, levels=c("rarefied.min.chao.3","rarefied.min.shannon.3","rarefied.min.simpson.3","rarefied.min.pielou.3"))) %>% drop_na() %>% ggplot(aes(x=reorder(factor(Cond.rarefied.min.3), -Cond.rarefied.min.3), y=value, fill=variable)) + geom_point(size=0.3) + stat_summary(fun=mean, geom="point", color="red", size=0.4, shape=23, position=position_dodge(0.75)) +  theme(axis.text.x =element_text(size=6, face = "bold", angle = 45 , hjust = 1),axis.text.y =element_text(size=6, face = "bold"),axis.title.y=element_blank(),axis.title.x = element_blank()) +  theme(legend.position = "none") + facet_wrap(vars(variable), scales = "free_y",labeller = as_labeller(new_var_names_rarefied_min_3),  ncol=1) + theme(strip.text.y = element_text(angle = 0, face = "bold", size = 8)) + theme(plot.title = element_blank()) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + geom_smooth(formula = y ~ x, method = "lm", se = TRUE, color = "grey40", aes(group = variable)) + stat_regline_equation(label.x = 1, label.y =0.7 ,aes(group = variable)) + stat_cor(label.x = 1, label.y = 1, aes(group = variable)) 


ggsave("rarefied.min.alpha.div.indices.cond.3.pdf",height = 11, width = 3)
ggsave("rarefied.min.alpha.div.indices.cond.3.png", height = 11, width = 3)


rarefied.min.alpha.div.indices.carb.3 <- diversity.indices.rarefied.min.long.3  %>% filter(!variable == "rarefied.min.no.3") %>% mutate(variable = factor(variable, levels=c("rarefied.min.chao.3","rarefied.min.shannon.3","rarefied.min.simpson.3","rarefied.min.pielou.3"))) %>% drop_na() %>% ggplot(aes(x=reorder(factor(Cond.rarefied.min.3), -Cond.rarefied.min.3), y=value, fill=variable)) + geom_point(size=0.3) + stat_summary(fun=mean, geom="point", color="red", size=0.4, shape=23, position=position_dodge(0.75)) +  theme(axis.text.x =element_text(size=6, face = "bold", angle = 45 , hjust = 1),axis.text.y =element_text(size=6, face = "bold"),axis.title.y=element_blank(),axis.title.x = element_blank()) +  theme(legend.position = "none") + facet_wrap(vars(variable), scales = "free_y",labeller = as_labeller(new_var_names_rarefied_min_3),  ncol=1) + theme(strip.text.y = element_text(angle = 0, face = "bold", size = 8)) + theme(plot.title = element_blank()) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + geom_smooth(formula = y ~ x, method = "lm", se = TRUE, color = "grey40", aes(group = variable)) + stat_regline_equation(label.x = 1, label.y =0.7 ,aes(group = variable)) + stat_cor(label.x = 1, label.y = 1000, aes(group = variable)) 

ggsave("rarefied.min.alpha.div.indices.carb.3.pdf",height = 11, width = 3)
ggsave("rarefied.min.alpha.div.indices.carb.3.png", height = 11, width = 3)

#################################################################################

### Fig. S4. Taxonomic composition of ASVs clustered at Phylum and Class level (N = 395) ####

#################################################################################

rarefied.class.prop.exclude.1.3 <- tax_glom(rarefied.min.prop.exclude.1.3, taxrank="Class")

rarefied.class.prop.wide.exclude.1.3 <- as.matrix(as.data.frame(rarefied.class.prop.exclude.1.3@otu_table))

colnames(rarefied.class.prop.wide.exclude.1.3) <- as.character(tax_table(rarefied.class.prop.exclude.1.3)[, "Class"])

meta <- rarefied.class.prop.wide.exclude.1.3 %>% select(Sample_id,Site,Temp.adj...C.)

## arrange all ASV tables
rarefied.class.prop.wide.exclude.1.3.reordered <- as.data.frame(rarefied.class.prop.wide.exclude.1.3) %>% rownames_to_column(var = 'Sample_id') 

## merge meta and reordered
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

class.row.annotations.left.2 <- rowAnnotation(Site = anno_block(gp = gpar(
  fill = c("#5881F9","#A7A3FE",  "#816AA7",  "#8CCBD3", "#0DD7FD",  "#0D8CA8",  "#2AB2FC" ,  "#0000FF" ,  "#4B3BA7",  "#4B00FD", "#0DDDD0" ,  "#0DE36D",  "#B2C9B2", "#99D584", "#9CD700", "#848F22",  "#009458",  "#2A7F72",  "#3B4D16",  "#E2BD7D" ,  "#DFC500" ,  "#FEAB2E", "#967D63", "#F76A4B", "#C8682A", "#8B5500",  "#EAB1DF",  "#E19AFE",  "#DE0DFD", "#A655ED",  "#BB00B8",  "#7E1C7C", "#FD85CA",  "#F9168D", "#FD0DD1",  "#FD8A91",  "#990D6A","#CC003D",  "#F81626",  "#950D3D")),
  labels = c("PB_38" , "RB_40" , "HS_42" , "MK_42" , "RN_42" , "MK_44" , "HS_45", "RB_45" , "RN_45", "SE_45", "AP_46" , "LA_46", "LN_46" , "PB_46" , "SW_46" , "KJ_47" , "PP_48" , "SE_49", "BA_50", "MK_50"  , "PT_50"  , "RB_50" , "HS_51" , "LN_51" , "RN_51" , "SW_51" , "HS_53" , "PP_53" , "MK_55", "PT_55" , "RB_55" , "SW_55" , "US_57" , "PP_58" , "RN_61", "SW_61", "LN_63","PT_63" , "US_63" , "LN_66"),
  labels_gp = gpar(col = "black", fontsize = 8, fontface = "bold")))

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

## get metadata for split
meta.avg <- as.data.frame((meta.avg[order(meta.avg[,2],decreasing=FALSE),]))
rownames(meta.avg) <- meta.avg$Site

class.row.annotations.left.2 <- rowAnnotation(Site = anno_block(gp = gpar(
  fill = c("#5881F9","#A7A3FE",  "#816AA7",  "#8CCBD3", "#0DD7FD",  "#0D8CA8",  "#2AB2FC" ,  "#0000FF" ,  "#4B3BA7",  "#4B00FD", "#0DDDD0" ,  "#0DE36D",  "#B2C9B2", "#99D584", "#9CD700", "#848F22",  "#009458",  "#2A7F72",  "#3B4D16",  "#E2BD7D" ,  "#DFC500" ,  "#FEAB2E", "#967D63", "#F76A4B", "#C8682A", "#8B5500",  "#EAB1DF",  "#E19AFE",  "#DE0DFD", "#A655ED",  "#BB00B8",  "#7E1C7C", "#FD85CA",  "#F9168D", "#FD0DD1",  "#FD8A91",  "#990D6A","#CC003D",  "#F81626",  "#950D3D")),
  labels = c("PB_38" , "RB_40" , "HS_42" , "MK_42" , "RN_42" , "MK_44" , "HS_45", "RB_45" , "RN_45", "SE_45", "AP_46" , "LA_46", "LN_46" , "PB_46" , "SW_46" , "KJ_47" , "PP_48" , "SE_49", "BA_50", "MK_50"  , "PT_50"  , "RB_50" , "HS_51" , "LN_51" , "RN_51" , "SW_51" , "HS_53" , "PP_53" , "MK_55", "PT_55" , "RB_55" , "SW_55" , "US_57" , "PP_58" , "RN_61", "SW_61", "LN_63","PT_63" , "US_63" , "LN_66"),
  labels_gp = gpar(col = "black", fontsize = 8, fontface = "bold"), labels_rot = 0))

pdf(file="rarefied.class.exclude.1.3.avg.compheatmap.nosplit.pdf", width = 15, height = 11)

rarefied.class.exclude.1.3.avg.compheatmap.nosplit <- Heatmap((as.matrix(rarefied.class.prop.wide.exclude.1.3.avg.reordered)), 
                                                              name = "Relative abundance",heatmap_legend_param = list(title = "Relative abundance (%)", title_gp = gpar(fontsize = 10, fontface = 'bold'), labels_gp = gpar(fontsize = 9),ncol = 3, legend_height = unit(7, "cm"),legend_width = unit(4, "cm"), legend_direction = "horizontal"),
                                                              column_title = "Rarefied heatmap of average relative abundances at the Class level", column_title_gp = gpar(fontsize = 13, fontface = 'bold'), row_title = "Samples",row_title_gp = gpar(fontsize = 10,fontface = 'bold'),row_names_gp = gpar(fontsize = 9), column_names_gp = gpar(fontsize = 8.9, fontface = "italic"), col = scalebr5,cluster_rows = FALSE, cluster_columns = FALSE , left_annotation = avg.row.annotations.left.1, bottom_annotation = rare.class.column.annotations.ex1,row_names_side = "left", row_names_centered = TRUE,  column_names_rot = 45, use_raster = TRUE, raster_by_magick = TRUE,width = ncol(class.prop.wide.exclude.1.3.avg.reordered)*unit(5, "mm"), height = nrow(class.prop.wide.exclude.1.3.avg.reordered)*unit(5, "mm"))

#left_annotation = class.row.annotations.left.2, split = data.frame(Site = test$Site),

rarefied.class.exclude.1.3.avg.compheatmap.nosplit <- draw(rarefied.class.exclude.1.3.avg.compheatmap.nosplit , heatmap_legend_side="bottom", annotation_legend_side="right",legend_grouping = "original")

pdf(file="rarefied.chloroflexi.exclude.1.3.avg.compheatmap.nosplit.pdf", width = 15, height = 11)


####################################################################################3

### Fig. 3A Relative abundance of the 25 most abundant ASVs in each sample (N = 395) shown clustered at genus level (bar plots). ####
### Fig. S6 Relative abundance of the 10 most abundant ASVs in ecological groups for each sample (N = 395). ####

###################################################################################

### use microviz to calculate PCoA ordination with unweighted unifrac and rooted tree as microviz needs a rooted tree and does not internally root tree like phyloseq.
### microviz prefers counts instead of proportions or percentages

#new tax tables
Tax.rare.wide.new.og = tax_table(new.rarefied.min.int.exclude.1.3@tax_table)

samples = sample_data(metadata)

phytree.rare = phy_tree(tree.rare.3)

## check if initial tree is unrooted or rooted
is.rooted(phytree.rare)

## root tree with ape ####
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

## check if correct number of tip labels are present by clicking on tree object

### new complete rarefied phyloseqs with original taxtable,rooted  and unrooted phylogenetic trees #####

### phyloseqs at 1% excluded ####

com.rarefied.min.int.exclude.1.3.unrooted <- phyloseq(rare.ASV.wide.new,Tax.rare.wide.new.og,samples,phytree.rare,refseqs.rare.3)

com.rarefied.min.prop.exclude.1.3.unrooted <- phyloseq(rare.prop.wide.new,Tax.rare.wide.new.og,samples,phytree.rare,refseqs.rare.3)

com.rarefied.min.int.exclude.1.3.rooted <- phyloseq(rare.ASV.wide.new,Tax.rare.wide.new.og,samples,rootedTree.rare,refseqs.rare.3)

com.rarefied.min.prop.exclude.1.3.rooted <- phyloseq(rare.prop.wide.new,Tax.rare.wide.new.og,samples,rootedTree.rare,refseqs.rare.3)

test <- as.data.frame(com.rarefied.min.prop.exclude.1.3.rooted@otu_table)

### Save phyloseqs with 1% excluded taxa, metadata with region for core microbiome

saveRDS(com.rarefied.min.int.exclude.1.3.rooted, "com.rarefied.min.int.exclude.1.3.rooted.rds")
saveRDS(com.rarefied.min.prop.exclude.1.3.rooted, "com.rarefied.min.prop.exclude.1.3.rooted.rds")


sample_variables(com.rarefied.min.int.exclude.1.3.rooted)
get_variable(com.rarefied.min.int.exclude.1.3.rooted,"Site")


### rooted phyloseqs with integers of photosynthetic and heterotrophic taxa ####
### unlike phyloseq, microviz requires integers instead of proportions

## get complete phyloseqs with subsets

## Chloroflexi
rarefied.chloroflexi.int.exclude.1.3.rooted <- com.rarefied.min.int.exclude.1.3.rooted %>%
  subset_taxa(Phylum == "Chloroflexi" & Class  == "Chloroflexia")

## Cyanobacteria
rarefied.cyanobacteria.int.exclude.1.3.rooted <- com.rarefied.min.int.exclude.1.3.rooted %>%
  subset_taxa(Phylum == "Cyanobacteria" & Class  == "Cyanobacteriia")

##  Phototynthetic Proteoacteria and other photosynthetic groups

proteobacteria.genus <- c("DSSF69","Elioraea","Methylobacterium-Methylorubrum","Rhodomicrobium","Roseomonas","Sandaracinobacter","Tabrizicola","Unassigned Rhodobacteraceae (Family)","Unassigned Sphingomonadaceae (Family)","AAP99","Allochromatium","Caldimonas","Curvibacter","DSSD61","Thiolamprovum","Unassigned B1-7BS (Family)","Unassigned Burkholderiales (Order)","Unassigned Comamonadaceae (Family)","Unassigned Gammaproteobacteria (Class)","Unassigned Rhodocyclaceae (Family)","Unassigned Sutterellaceae (Family)","Z-35")

other.photo <- c("Chlorobiales","Chloracidobacteriales")

rarefied.proteobacteria.otherphoto.int.exclude.1.3.rooted <- com.rarefied.min.int.exclude.1.3.rooted %>%
  subset_taxa(Order %in% other.photo | Genus %in% proteobacteria.genus)

## all photosynthetic
photosynthetic.class <- c("Cyanobacteriia","Chloroflexia")

rarefied.photosynthetic.int.exclude.1.3.rooted <- com.rarefied.min.int.exclude.1.3.rooted %>%
  subset_taxa( Class %in% photosynthetic.class | Order %in% other.photo | Genus %in% proteobacteria.genus )

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

rarefied.chemolithotrophs.int.exclude.1.3.rooted <- com.rarefied.min.int.exclude.1.3.rooted %>%
  subset_taxa(Phylum %in% chemolithotroph.phylum  |  Class %in% chemolithotroph.class)

#### all  heterotrophs only
rarefied.only.heterotrophs.int.exclude.1.3.rooted <- com.rarefied.min.int.exclude.1.3.rooted %>%
  subset_taxa(!Phylum %in% chemolithotroph.phylum & !Class %in% chemolithotroph.class & !Class %in% photosynthetic.class & !Order %in% other.photo & !Genus %in% proteobacteria.genus )

### save all phyloseqs

saveRDS(rarefied.chloroflexi.int.exclude.1.3.rooted, "rarefied.chloroflexi.int.exclude.1.3.rooted.rds")
saveRDS(rarefied.cyanobacteria.int.exclude.1.3.rooted, "rarefied.cyanobacteria.int.exclude.1.3.rooted.rds")
saveRDS(rarefied.proteobacteria.otherphoto.int.exclude.1.3.rooted, "rarefied.proteobacteria.otherphoto.int.exclude.1.3.rooted.rds")
saveRDS(rarefied.photosynthetic.int.exclude.1.3.rooted, "rarefied.photosynthetic.int.exclude.1.3.rooted.rds")
saveRDS(rarefied.chemolithotrophs.int.exclude.1.3.rooted, "rarefied.chemolithotrophs.int.exclude.1.3.rooted.rds")
saveRDS(rarefied.only.heterotrophs.int.exclude.1.3.rooted, "rarefied.only.heterotrophs.int.exclude.1.3.rooted.rds")

### Use microviz to calculate PCoA ordination with unweighted unifrac and rooted tree as microviz needs a rooted tree and does not internally root tree like phyloseq.
### microviz prefers counts instead of proportions or percentages

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

com.rarefied.min.propTo1.exclude.1.3 <- transform(com.rarefied.min.int.exclude.1.3.rooted, 'compositional')

is_compositional(com.rarefied.min.propTo1.exclude.1.3)

### get core microbiome for each PCoA cluster
metadata.ex.1.mb <- as.matrix(as.data.frame(com.rarefied.min.propTo1.exclude.1.3@sam_data))
sample_variables(com.rarefied.min.propTo1.exclude.1.3)

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

core.genus.all.clusters.95 <- list( 'North.Thailand'= c(taxa(core.rarefied.exc1.cluster1.genus.95))[ !(c(taxa(core.rarefied.exc1.cluster1.genus.95))) == 'Other'], 'North.Malaysia' = c(taxa(core.rarefied.exc1.cluster2.genus.95))[ !(c(taxa(core.rarefied.exc1.cluster2.genus.95))) == 'Other'], 'South.Malaysia' = c(taxa(core.rarefied.exc1.cluster3.genus.95))[ !(c(taxa(core.rarefied.exc1.cluster3.genus.95))) == 'Other'], 'Singapore' = c(taxa(core.rarefied.exc1.cluster4.genus.95))[ !(c(taxa(core.rarefied.exc1.cluster4.genus.95))) == 'Other'], 'South.Thailand' = c(taxa(core.rarefied.exc1.cluster5.genus.95))[ !(c(taxa(core.rarefied.exc1.cluster5.genus.95))) == 'Other'], 'Central.Thailand' = c(taxa(core.rarefied.exc1.cluster6.genus.95))[ !(c(taxa(core.rarefied.exc1.cluster6.genus.95))) == 'Other'])

### use to subset from phyloseq
core.genus.all.clusters.95.tosub <- c((taxa(core.rarefied.exc1.cluster1.genus.95))[ !(c(taxa(core.rarefied.exc1.cluster1.genus.95))) == 'Other'],c(taxa(core.rarefied.exc1.cluster2.genus.95))[ !(c(taxa(core.rarefied.exc1.cluster2.genus.95))) == 'Other'], c(taxa(core.rarefied.exc1.cluster3.genus.95))[ !(c(taxa(core.rarefied.exc1.cluster3.genus.95))) == 'Other'], c(taxa(core.rarefied.exc1.cluster4.genus.95))[ !(c(taxa(core.rarefied.exc1.cluster4.genus.95))) == 'Other'], c(taxa(core.rarefied.exc1.cluster5.genus.95))[ !(c(taxa(core.rarefied.exc1.cluster5.genus.95))) == 'Other'], c(taxa(core.rarefied.exc1.cluster6.genus.95))[ !(c(taxa(core.rarefied.exc1.cluster6.genus.95))) == 'Other'])

core.genus.all.clusters.95.tosub <- sort(unique(core.genus.all.clusters.95.tosub))


### Plot core genera as pie charts - entire dataset ####

entiredata <- com.rarefied.min.int.exclude.1.3.rooted
test <- as.data.frame(entiredata@otu_table)

entiredata.genera <- tax_glom(entiredata, taxrank = "Genus")

test <- as.data.frame(entiredata.genera@tax_table)

entiredata.genera.perc <- transform_sample_counts(entiredata.genera, function(OTU) (OTU/sum(OTU))*100)

entiredata.genera.melt <- psmelt(entiredata.genera.perc)

#change to character for easy-adjusted level
entiredata.genera.melt$Genus <- as.character(entiredata.genera.melt$Genus)


# Group according to faceting variable to get mean and median for further filtering
entiredata.genera.melt <- entiredata.genera.melt %>%
  group_by(Genus) %>%
  mutate(median=median(Abundance), mean=mean(Abundance))


## get mean abundances of genus identified at each cluster
entiredata.genera.melt.2 <- entiredata.genera.melt %>%
  group_by(Genus) %>%
  summarise(Abundance = mean(Abundance)) 

# to label core and non-core genus ###

entire.dataset.core.95.genus.tosub

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





## labels + percent
entiredata.genera.pie.1 <- plot_ly(entiredata.genera.melt.5, labels = ~Genus, values = ~Abundance, type = 'pie',marker = list(colors = core.noncore,line = list(color = '#FFFFFF', width = 2)), textposition = 'inside', textinfo = 'label+percent', insidetextfont = list(color = 'grey32'),hoverinfo = 'text') %>% layout(title = 'Entire dataset core and non-core genera',
                                                                                                                                                                                                                                                                                                                            xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                                                                                                                                                                                                                                                                                                                            yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))

## percent
entiredata.genera.pie.2 <- plot_ly(entiredata.genera.melt.5, labels = ~Genus, values = ~Abundance, type = 'pie',marker = list(colors = core.noncore,line = list(color = '#FFFFFF', width = 2)), textposition = 'inside', textinfo = 'percent', insidetextfont = list(color = 'grey32'),hoverinfo = 'text') %>% layout(title = 'Entire dataset core and non-core genera',
                                                                                                                                                                                                                                                                                                                      xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                                                                                                                                                                                                                                                                                                                      yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))

## no labels
entiredata.genera.pie.3 <- plot_ly(entiredata.genera.melt.5, labels = ~Genus, values = ~Abundance, type = 'pie',marker = list(colors = core.noncore, line = list(color = '#FFFFFF', width = 2)), textposition = 'inside', textinfo = 'none', insidetextfont = list(color = 'grey32'),hoverinfo = 'text') %>% layout(title = 'Entire dataset core and non-core genera',
                                                                                                                                                                                                                                                                                                                    xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                                                                                                                                                                                                                                                                                                                    yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))


entiredata.genera.pie.1
entiredata.genera.pie.2
entiredata.genera.pie.3

## core genera pie chart ####

entiredata.core.95.genus <- entiredata.genera  %>%
  subset_taxa(Genus %in% entire.dataset.core.95.genus.tosub )

entiredata.core.95.genus <- 
  prune_taxa(taxa_sums(entiredata.core.95.genus) > 0, entiredata.core.95.genus)

#sanity check
test <- as.data.frame(entiredata.core.95.genus@tax_table)
test$Genus
entire.dataset.core.95.genus.tosub
identical(North.Thailand.95.core,sort(test$Genus))

entiredata.core.95.perc <- transform_sample_counts(entiredata.core.95.genus, function(OTU) (OTU/sum(OTU))*100)

test <- as.data.frame(entiredata.core.95.perc@tax_table)

sample_sums(entiredata.core.95.perc)

entiredata.core.95.perc.melt <- psmelt(entiredata.core.95.perc)

#change to character for easy-adjusted level
entiredata.core.95.perc.melt$Genus <- as.character(entiredata.core.95.perc.melt$Genus)

# Group according to faceting variable to get mean and median for further filtering
entiredata.core.95.perc.melt <- entiredata.core.95.perc.melt %>%
  group_by(Genus) %>%
  mutate(median=median(Abundance), mean=mean(Abundance))


### ## get mean abundances of genus identified at each cluster
entiredata.core.95.perc.melt.2 <- entiredata.core.95.perc.melt %>%
  group_by(Genus) %>%
  summarise(Abundance = mean(Abundance)) 

## without excluding any genus
entiredata.core.95.perc.melt.3 <- entiredata.core.95.perc.melt.2 %>%
  group_by(Genus) %>%
  summarise(Abundance = sum(Abundance))

## Genus with relative abundances > 5
## Genus with relativ abundance > 10 only resulted in 2 genera
entiredata.core.95.perc.melt.4 <- entiredata.core.95.perc.melt.2

keep.5 <- unique(entiredata.core.95.perc.melt.4$Genus[entiredata.core.95.perc.melt.4$Abundance >= 5]) 
entiredata.core.95.perc.melt.4$Genus[!(entiredata.core.95.perc.melt.4$Genus %in% keep.5)] <- "< 5%"

## to check if to include <5% as others in pie-chart for core-genera, but others make up 43%
entiredata.core.95.perc.melt.5 <- entiredata.core.95.perc.melt.4 %>%
  group_by(Genus) %>%
  summarise(Abundance = sum(Abundance))

entiredata.core.95.perc.melt.3 <- as.data.frame(entiredata.core.95.perc.melt.3)
rownames(entiredata.core.95.perc.melt.3) <- entiredata.core.95.perc.melt.3$Genus

entiredata.core.95.perc.melt.5 <- as.data.frame(entiredata.core.95.perc.melt.5)
rownames(entiredata.core.95.perc.melt.5) <- entiredata.core.95.perc.melt.5$Genus

## colors
colors <- c('#A6CEE3', '#1F78B4', '#B2DF8A', '#33A02C', '#FB9A99','#E31A1C','#FDBF6F','#FF7F00', "#68A6CD", "#2A7FB7", "#569EA4", "#99CD91","#8CCC6E", "#52AF43","#5C9E42", "#B89B74", "#F88A89", "#ED4F50", "#E4201F","#F06C45","#FBB86B","#FDA440", "#FE870D", "#ED8F47", "#D5A7A9"  ,"#B294C7","#865FAB", "#825D99","#C7B699", "#F8F18F", "#D4A55B", "#B15928")

[1] "#66C2A5" "#FC8D62" "#8DA0CB" "#E78AC3" "#A6D854" "#FFD92F" "#E5C494"
[8] "#B3B3B3"

[1] "#1B9E77" "#D95F02" "#7570B3" "#E7298A" "#66A61E" "#E6AB02" "#A6761D"
[8] "#666666"

[1] "#A6CEE3" "#97C4DD" "#88BAD8" "#79B1D3" "#69A7CE" "#5B9EC8" "#4C94C3" "#3C8BBE" "#2D81B9" "#1F78B4" "#2F83AF" "#3F8EAA" "#509AA6"
[14] "#60A5A1" "#70B19C" "#81BC98" "#91C893" "#A1D38E" "#B2DF8A" "#A3D87F" "#95D175" "#87CA6A" "#79C360" "#6BBB55" "#5DB44B" "#4FAD40"
[27] "#41A636" "#33A02C" "#499F38" "#5F9E44" "#759E50" "#8B9D5C" "#A29C68" "#B89B74" "#CE9B80" "#E49A8C" "#FB9A99" "#F88B8B" "#F57D7D"
[40] "#F36F6F" "#F06161" "#ED5253" "#EB4445" "#E83637" "#E52829" "#E31A1C" "#E52C25" "#E83E2E" "#EB5137" "#EE6340" "#F1754A" "#F48853"
[53] "#F79A5C" "#FAAC65" "#FDBE6E" "#FDB762" "#FDB056" "#FDA949" "#FDA23D" "#FE9B31" "#FE9424" "#FE8D18" "#FE860C" "#FE7F00" "#F98417"
[66] "#F38A2F" "#ED9047" "#E7955F" "#E19B76" "#DBA18E" "#D5A6A6" "#CFACBE" "#CAB2D6" "#BFA4CF" "#B497C8" "#A98AC1" "#9F7DBB" "#9471B4"
[79] "#8963AD" "#7F57A7" "#7449A0" "#6A3D9A" "#7A5299" "#8B6899" "#9B7D99" "#AC9399" "#BCA899" "#CDBE99" "#DDD399" "#EEE999" "#FFFF99"
[92] "#F6EC8C" "#EDDA7F" "#E5C773" "#DCB566" "#D3A25A" "#CB904D" "#C27D41" "#B96B34" "#B15928"

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

entire.dataset.core <- c("#CDBE99",'#7570B3','#1F78B4','#B2DF8A','#E31A1C','#B15928',"#FFFF99",'#551A8B')


##label +  percentages
entiredata.core.genera.pie.1 <- plot_ly(entiredata.core.95.perc.melt.5, labels = ~Genus, values = ~Abundance, type = 'pie', marker = list(colors = entire.dataset.core, line = list(color = '#FFFFFF', width = 2)), textposition = 'inside', textinfo = 'label+percent', insidetextfont = list(color = 'grey32'),hoverinfo = 'text') %>% layout(title = 'Entire dataset core genera',
                                                                                                                                                                                                                                                                                                                                                xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                                                                                                                                                                                                                                                                                                                                                yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))


##percent
entiredata.core.genera.pie.2 <- plot_ly(entiredata.core.95.perc.melt.5, labels = ~Genus, values = ~Abundance, type = 'pie', marker = list(colors = entire.dataset.core, line = list(color = '#FFFFFF', width = 2)), textposition = 'inside', textinfo = 'percent', insidetextfont = list(color = 'grey32'),hoverinfo = 'text') %>% layout(title = 'Entire dataset core genera',
                                                                                                                                                                                                                                                                                                                                          xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                                                                                                                                                                                                                                                                                                                                          yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))


##no labels
entiredata.core.genera.pie.3 <- plot_ly(entiredata.core.95.perc.melt.5, values = ~Abundance, type = 'pie', marker = list(colors = entire.dataset.core ,line = list(color = 'white', width = 2)),textinfo = 'none') %>% layout(title = 'Entire dataset core genera', 
                                                                                                                                                                                                                              xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                                                                                                                                                                                                                              yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))


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

test <- as.data.frame(cluster1.genera@tax_table)

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
North.Thailand.95.core
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


cluster1.genera.pie

##arrange order of sample_id's for long table in alphanumerical order
#core.95prev.perc.melt.1 <- core.95prev.perc.melt.1 [gtools::mixedorder(core.95prev.perc.melt.1$Sample_id), ]
#cluster1.genera.melt.5$Genus <- factor(cluster1.genera.melt.5$Genus, #levels=unique(mixedsort(as.character(cluster1.genera.melt.5$Genus))))

#choose number and ramp up colors
#length(unique(cluster1.genera.melt.5$Genus))
#colourCountGenus2 = length(unique(cluster1.genera.melt.5$Genus))
#getPalette = colorRampPalette(brewer.pal(12, "Paired"))

#assign order names to colors to ensure subsetted order have same colors as full plot
#value.genus.2 = getPalette(colourCountGenus2)
#names(value.genus.2) = unique(cluster1.genera.melt.5$Genus)
#print(value.genus.2)

# Add label position
#cluster1.genera.melt.5 <- cluster1.genera.melt.5 %>%
#  arrange(desc(Genus)) %>%
#  mutate(lab.ypos = cumsum(Abundance) - 0.5*Abundance)

#cluster1.genera.pie <- ggplot(cluster1.genera.melt.5 , aes(x = "", y = Abundance, fill = Genus)) + geom_bar(stat = "identity", position="stack", width = 1, color = "white") + scale_fill_manual(values = value.genus.2) + theme_void() + theme(legend.position = "bottom", strip.background = element_blank(),legend.text = element_text(size=23), legend.key.size =  unit(1, "cm")) + coord_polar("y", start=0) + geom_text(aes(y = lab.ypos,label = percent(Abundance/100)), position = position_stack(vjust = 0.5), size=20,  color = "white") #+ guides(fill = guide_legend(nrow = 5)) 

#cluster1.genera.pie

#pdf("cluster1.genera.pie.pdf", width = 28, height = 15)
#cluster1.genera.pie
#dev.off()

## core genera pie chart ####
cluster1.core.genera <- cluster1.genera %>% subset_taxa(Genus %in% North.Thailand.95.core)

#sanity check
test <- as.data.frame(cluster1.core.genera@tax_table)
test$Genus
North.Thailand.95.core
identical(North.Thailand.95.core,sort(test$Genus))

cluster1.core.genera.perc <- transform_sample_counts(cluster1.core.genera, function(OTU) (OTU/sum(OTU))*100)

test <- as.data.frame(cluster1.core.genera.perc@tax_table)

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


## Genus with relative abundances > 5
## Genus with relativ abundance > 10 only resulted in 2 genera
cluster1.core.genera.perc.melt.3 <- cluster1.core.genera.perc.melt.2

keep.5 <- unique(cluster1.core.genera.perc.melt.3$Genus[cluster1.core.genera.perc.melt.3$Abundance >= 5]) 
cluster1.core.genera.perc.melt.3$Genus[!(cluster1.core.genera.perc.melt.3$Genus %in% keep.5)] <- "< 5%"

## to check if to include <5% as others in pie-chart for core-genera, but others make up 43%
cluster1.core.genera.perc.melt.4 <- cluster1.core.genera.perc.melt.3 %>%
  group_by(Genus) %>%
  summarise(Abundance = sum(Abundance))

##### filter out rare taxa = < 10% from dataset ####
#cluster1.core.genera.perc.melt.4  <- cluster1.core.genera.perc.melt.3 %>% filter(Genus != "< 5%")

# Add label position :  does not apply correctly on plot
cluster1.core.genera.perc.melt.5 <- cluster1.core.genera.perc.melt.4 %>%
  arrange(desc(Genus)) #%>%
#mutate(norm.abundance = (Abundance/sum(Abundance))*100) #%>%
#  mutate(lab.ypos = cumsum(Abundance) - 0.5*Abundance)


#choose number and ramp up colors
#length(unique(cluster1.core.genera.perc.melt.5$Genus))
#colourCountGenus6 = length(unique(cluster1.core.genera.perc.melt.5$Genus))
#getPalette = colorRampPalette(brewer.pal(12, "Paired"))

#assign order names to colors to ensure subsetted order have same colors as full plot
#value.genus.6 = getPalette(colourCountGenus6)
#names(value.genus.6) = unique(cluster1.core.genera.perc.melt.5$Genus)
#print(value.genus.6)

cluster1.core.genera.perc.melt.5 <- as.data.frame(cluster1.core.genera.perc.melt.5)
rownames(cluster1.core.genera.perc.melt.5) <- cluster1.core.genera.perc.melt.5$Genus

## colors
colors <- c('#A6CEE3', '#1F78B4', '#B2DF8A', '#33A02C', '#FB9A99','#E31A1C','#FDBF6F','#FF7F00', "#68A6CD", "#2A7FB7", "#569EA4", "#99CD91","#8CCC6E", "#52AF43","#5C9E42", "#B89B74", "#F88A89", "#ED4F50", "#E4201F","#F06C45","#FBB86B","#FDA440", "#FE870D", "#ED8F47", "#D5A7A9"  ,"#B294C7","#865FAB", "#825D99","#C7B699", "#F8F18F", "#D4A55B", "#B15928")


brewer.pal(n = 14, name = "Paired")

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

entire.dataset.core <- c("#CDBE99",'#7570B3','#1F78B4','#B2DF8A','#E31A1C','#B15928',"#FFFF99",'#551A8B')

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

#cluster1.core.genera.pie <- ggplot(cluster1.core.genera.perc.melt.5 , aes(x = "", y = norm.abundance, fill = Genus)) + geom_bar(stat = "identity", position="stack", width = 1, color = "white", aes(fill=fct_reorder(Genus,(norm.abundance)))) + scale_fill_manual(values = value.genus.6) + theme_void() + theme(legend.position = "bottom", strip.background = element_blank(),legend.text = element_text(size=5), legend.key.size =  unit(1, "cm")) + coord_polar("y", start=0) + geom_text(aes(y = lab.ypos,label = percent(norm.abundance/100)), position = position_stack(vjust = 0.5), size=5,  color = "white") #+ guides(fill = guide_legend(nrow = 5)) 

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


cluster2.genera.pie


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


## Genus with relative abundances > 5
## Genus with relativ abundance > 10 only resulted in 2 genera
cluster2.core.genera.perc.melt.3 <- cluster2.core.genera.perc.melt.2

keep.5 <- unique(cluster2.core.genera.perc.melt.3$Genus[cluster2.core.genera.perc.melt.3$Abundance >= 5]) 
cluster2.core.genera.perc.melt.3$Genus[!(cluster2.core.genera.perc.melt.3$Genus %in% keep.5)] <- "< 5%"

## to check if to include <5% as others in pie-chart for core-genera, but others make up 43%
cluster2.core.genera.perc.melt.4 <- cluster2.core.genera.perc.melt.3 %>%
  group_by(Genus) %>%
  summarise(Abundance = sum(Abundance))

##### filter out rare taxa = < 10% from dataset ####
#cluster2.core.genera.perc.melt.4  <- cluster2.core.genera.perc.melt.3 %>% filter(Genus != "< 5%")

# Add label position :  does not apply correctly on plot
cluster2.core.genera.perc.melt.5 <- cluster2.core.genera.perc.melt.4 %>%
  arrange(desc(Genus)) #%>%
#mutate(norm.abundance = (Abundance/sum(Abundance))*100) #%>%
#  mutate(lab.ypos = cumsum(Abundance) - 0.5*Abundance)

cluster2.core.genera.perc.melt.5 <- as.data.frame(cluster2.core.genera.perc.melt.5)
rownames(cluster2.core.genera.perc.melt.5) <- cluster2.core.genera.perc.melt.5$Genus

## colors
colors <- c('#A6CEE3', '#1F78B4', '#B2DF8A', '#33A02C', '#FB9A99','#E31A1C','#FDBF6F','#FF7F00', "#68A6CD", "#2A7FB7", "#569EA4", "#99CD91","#8CCC6E", "#52AF43","#5C9E42", "#B89B74", "#F88A89", "#ED4F50", "#E4201F","#F06C45","#FBB86B","#FDA440", "#FE870D", "#ED8F47", "#D5A7A9"  ,"#B294C7","#865FAB", "#825D99","#C7B699", "#F8F18F", "#D4A55B", "#B15928")


brewer.pal(n = 14, name = "Paired")

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

entire.dataset.core <- c("#CDBE99",'#7570B3','#1F78B4','#B2DF8A','#E31A1C','#B15928',"#FFFF99",'#551A8B')

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


## Genus with relative abundances > 5
## Genus with relativ abundance > 10 only resulted in 2 genera
cluster3.core.genera.perc.melt.3 <- cluster3.core.genera.perc.melt.2

keep.5 <- unique(cluster3.core.genera.perc.melt.3$Genus[cluster3.core.genera.perc.melt.3$Abundance >= 5]) 
cluster3.core.genera.perc.melt.3$Genus[!(cluster3.core.genera.perc.melt.3$Genus %in% keep.5)] <- "< 5%"

## to check if to include <5% as others in pie-chart for core-genera, but others make up 43%
cluster3.core.genera.perc.melt.4 <- cluster3.core.genera.perc.melt.3 %>%
  group_by(Genus) %>%
  summarise(Abundance = sum(Abundance))

##### filter out rare taxa = < 10% from dataset ####
#cluster3.core.genera.perc.melt.4  <- cluster3.core.genera.perc.melt.3 %>% filter(Genus != "< 5%")

# Add label position :  does not apply correctly on plot
cluster3.core.genera.perc.melt.5 <- cluster3.core.genera.perc.melt.4 %>%
  arrange(desc(Genus)) #%>%
#mutate(norm.abundance = (Abundance/sum(Abundance))*100) #%>%
#  mutate(lab.ypos = cumsum(Abundance) - 0.5*Abundance)

cluster3.core.genera.perc.melt.5 <- as.data.frame(cluster3.core.genera.perc.melt.5)
rownames(cluster3.core.genera.perc.melt.5) <- cluster3.core.genera.perc.melt.5$Genus

## colors
colors <- c('#A6CEE3', '#1F78B4', '#B2DF8A', '#33A02C', '#FB9A99','#E31A1C','#FDBF6F','#FF7F00', "#68A6CD", "#2A7FB7", "#569EA4", "#99CD91","#8CCC6E", "#52AF43","#5C9E42", "#B89B74", "#F88A89", "#ED4F50", "#E4201F","#F06C45","#FBB86B","#FDA440", "#FE870D", "#ED8F47", "#D5A7A9"  ,"#B294C7","#865FAB", "#825D99","#C7B699", "#F8F18F", "#D4A55B", "#B15928")

[1] "#66C2A5" "#FC8D62" "#8DA0CB" "#E78AC3" "#A6D854" "#FFD92F" "#E5C494"
[8] "#B3B3B3"

[1] "#1B9E77" "#D95F02" "#7570B3" "#E7298A" "#66A61E" "#E6AB02" "#A6761D"
[8] "#666666"


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

entire.dataset.core <- c("#CDBE99",'#7570B3','#1F78B4','#B2DF8A','#E31A1C','#B15928',"#FFFF99",'#551A8B')

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


## Genus with relative abundances > 5
## Genus with relativ abundance > 10 only resulted in 2 genera
cluster4.core.genera.perc.melt.3 <- cluster4.core.genera.perc.melt.2

keep.5 <- unique(cluster4.core.genera.perc.melt.3$Genus[cluster4.core.genera.perc.melt.3$Abundance >= 5]) 
cluster4.core.genera.perc.melt.3$Genus[!(cluster4.core.genera.perc.melt.3$Genus %in% keep.5)] <- "< 5%"

## to check if to include <5% as others in pie-chart for core-genera, but others make up 43%
cluster4.core.genera.perc.melt.4  <- cluster4.core.genera.perc.melt.3 %>%
  group_by(Genus) %>%
  summarise(Abundance = sum(Abundance))

##### filter out rare taxa = < 10% from dataset ####
#cluster4.core.genera.perc.melt.4  <- cluster4.core.genera.perc.melt.3 %>% filter(Genus != "< 5%")

# Add label position :  does not apply correctly on plot
cluster4.core.genera.perc.melt.5 <- cluster4.core.genera.perc.melt.4 %>%
  arrange(desc(Genus)) #%>%
#mutate(norm.abundance = (Abundance/sum(Abundance))*100) #%>%
#  mutate(lab.ypos = cumsum(Abundance) - 0.5*Abundance)

cluster4.core.genera.perc.melt.5 <- as.data.frame(cluster4.core.genera.perc.melt.5)
rownames(cluster4.core.genera.perc.melt.5) <- cluster4.core.genera.perc.melt.5$Genus

## colors
colors <- c('#A6CEE3', '#1F78B4', '#B2DF8A', '#33A02C', '#FB9A99','#E31A1C','#FDBF6F','#FF7F00', "#68A6CD", "#2A7FB7", "#569EA4", "#99CD91","#8CCC6E", "#52AF43","#5C9E42", "#B89B74", "#F88A89", "#ED4F50", "#E4201F","#F06C45","#FBB86B","#FDA440", "#FE870D", "#ED8F47", "#D5A7A9"  ,"#B294C7","#865FAB", "#825D99","#C7B699", "#F8F18F", "#D4A55B", "#B15928")

[1] "#66C2A5" "#FC8D62" "#8DA0CB" "#E78AC3" "#A6D854" "#FFD92F" "#E5C494"
[8] "#B3B3B3"

[1] "#1B9E77" "#D95F02" "#7570B3" "#E7298A" "#66A61E" "#E6AB02" "#A6761D"
[8] "#666666"


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

entire.dataset.core <- c("#CDBE99",'#7570B3','#1F78B4','#B2DF8A','#E31A1C','#B15928',"#FFFF99",'#551A8B')

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


## Genus with relative abundances > 5
## Genus with relativ abundance > 10 only resulted in 2 genera
cluster5.core.genera.perc.melt.3 <- cluster5.core.genera.perc.melt.2

keep.5 <- unique(cluster5.core.genera.perc.melt.3$Genus[cluster5.core.genera.perc.melt.3$Abundance >= 5]) 
cluster5.core.genera.perc.melt.3$Genus[!(cluster5.core.genera.perc.melt.3$Genus %in% keep.5)] <- "< 5%"

## to check if to include <5% as others in pie-chart for core-genera, but others make up 43%
cluster5.core.genera.perc.melt.4 <- cluster5.core.genera.perc.melt.3 %>%
  group_by(Genus) %>%
  summarise(Abundance = sum(Abundance))

##### filter out rare taxa = < 10% from dataset ####
#cluster5.core.genera.perc.melt.4  <- cluster5.core.genera.perc.melt.3 %>% filter(Genus != "< 5%")

# Add label position :  does not apply correctly on plot
cluster5.core.genera.perc.melt.5 <- cluster5.core.genera.perc.melt.4 %>%
  arrange(desc(Genus)) #%>%
#mutate(norm.abundance = (Abundance/sum(Abundance))*100) #%>%
#  mutate(lab.ypos = cumsum(Abundance) - 0.5*Abundance)

cluster5.core.genera.perc.melt.5 <- as.data.frame(cluster5.core.genera.perc.melt.5)
rownames(cluster5.core.genera.perc.melt.5) <- cluster5.core.genera.perc.melt.5$Genus

## colors
colors <- c('#A6CEE3', '#1F78B4', '#B2DF8A', '#33A02C', '#FB9A99','#E31A1C','#FDBF6F','#FF7F00', "#68A6CD", "#2A7FB7", "#569EA4", "#99CD91","#8CCC6E", "#52AF43","#5C9E42", "#B89B74", "#F88A89", "#ED4F50", "#E4201F","#F06C45","#FBB86B","#FDA440", "#FE870D", "#ED8F47", "#D5A7A9"  ,"#B294C7","#865FAB", "#825D99","#C7B699", "#F8F18F", "#D4A55B", "#B15928")

[1] "#66C2A5" "#FC8D62" "#8DA0CB" "#E78AC3" "#A6D854" "#FFD92F" "#E5C494"
[8] "#B3B3B3"

[1] "#1B9E77" "#D95F02" "#7570B3" "#E7298A" "#66A61E" "#E6AB02" "#A6761D"
[8] "#666666"


> 
  
  
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

entire.dataset.core <- c("#CDBE99",'#7570B3','#1F78B4','#B2DF8A','#E31A1C','#B15928',"#FFFF99",'#551A8B')

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


## Genus with relative abundances > 5
## Genus with relativ abundance > 10 only resulted in 2 genera
cluster6.core.genera.perc.melt.3 <- cluster6.core.genera.perc.melt.2

keep.5 <- unique(cluster6.core.genera.perc.melt.3$Genus[cluster6.core.genera.perc.melt.3$Abundance >= 5]) 
cluster6.core.genera.perc.melt.3$Genus[!(cluster6.core.genera.perc.melt.3$Genus %in% keep.5)] <- "< 5%"

## to check if to include <5% as others in pie-chart for core-genera, but others make up 43%
cluster6.core.genera.perc.melt.4  <- cluster6.core.genera.perc.melt.3 %>%
  group_by(Genus) %>%
  summarise(Abundance = sum(Abundance))

##### filter out rare taxa = < 10% from dataset ####
#cluster6.core.genera.perc.melt.4  <- cluster6.core.genera.perc.melt.3 %>% filter(Genus != "< 5%")

# Add label position :  does not apply correctly on plot
cluster6.core.genera.perc.melt.5 <- cluster6.core.genera.perc.melt.4 %>%
  arrange(desc(Genus)) #%>%
#mutate(norm.abundance = (Abundance/sum(Abundance))*100) #%>%
#  mutate(lab.ypos = cumsum(Abundance) - 0.5*Abundance)

cluster6.core.genera.perc.melt.5 <- as.data.frame(cluster6.core.genera.perc.melt.5)
rownames(cluster6.core.genera.perc.melt.5) <- cluster6.core.genera.perc.melt.5$Genus

## colors
colors <- c('#A6CEE3', '#1F78B4', '#B2DF8A', '#33A02C', '#FB9A99','#E31A1C','#FDBF6F','#FF7F00', "#68A6CD", "#2A7FB7", "#569EA4", "#99CD91","#8CCC6E", "#52AF43","#5C9E42", "#B89B74", "#F88A89", "#ED4F50", "#E4201F","#F06C45","#FBB86B","#FDA440", "#FE870D", "#ED8F47", "#D5A7A9"  ,"#B294C7","#865FAB", "#825D99","#C7B699", "#F8F18F", "#D4A55B", "#B15928")

[1] "#66C2A5" "#FC8D62" "#8DA0CB" "#E78AC3" "#A6D854" "#FFD92F" "#E5C494"
[8] "#B3B3B3"

[1] "#1B9E77" "#D95F02" "#7570B3" "#E7298A" "#66A61E" "#E6AB02" "#A6761D"
[8] "#666666"


> 
  
  
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

entire.dataset.core <- c("#CDBE99",'#7570B3','#1F78B4','#B2DF8A','#E31A1C','#B15928',"#FFFF99",'#551A8B')


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

### Fig. 4A. Mantel's test on the influence of abiotic variables on observed taxonomic composition for all taxa plus photosynthetic, chemoautotrophic, and chemoheterotrophic fractions of the community. ####
### Fig. 4B. Variance partioning ####

######################################################################################

com.rarefied.min.int.exclude.1.3.rooted <- readRDS("com.rarefied.min.int.exclude.1.3.rooted.rds")
net_spieceasi_taxtable <- readRDS("net_spieceasi_taxtable.rds")
microeco.otu <- readRDS("microeco.otu.rds")
microeco.dataset <- readRDS("microeco.dataset.rds")

### for info on Venenivibrio
## A genus in the bacterial phylum Aquificota appears to be endemic to Aotearoa- New Zealand : A specific environmental niche that increases habitat isolation was identified, with maximal read abundance of Venenivibrio occurring at pH 4-6, 50-70 °C, and low oxidation-reduction potentials.

complete.dataset.melt <- psmelt(com.rarefied.min.int.exclude.1.3.rooted)
venenivibrio <- complete.dataset.melt %>% filter(Genus == "Venenivibrio")

NT.ASVs <-  com.rarefied.min.int.exclude.1.3.rooted 

## Mantel's test and correlation heatmap #####

## https://chiliubio.github.io/microeco_tutorial/other-examples-1.html#mantel-test-correlation-heatmap


# prepare data
tax <- as.data.frame(net_spieceasi_taxtable) 
tax <- as.data.frame(tax) %>% column_to_rownames(var = "ASV")

metadata <- as.data.frame(sample_data(com.rarefied.min.int.exclude.1.3.rooted)) %>% as_tibble()
colnames(metadata)
metadata.micronet <- metadata %>% select("Sample_id","Country","Region","Location.Code","Site","Location","Human.usage..Y.N.","Water.flow..Pool.Flowing.","Copper..mg.L.","Cyanuric.acid..mg.L.","Nitrate..ppm.","Nitrite..ppm.","Temp.adj...C.","Chloride..mg.L.","Chlorine.dioxide..mg.L.","Flouride..mg.L.","Hardness..mg.L.","Lead..mg.L.","MPS..mg.L.","QUAT.QAC..mg.L.","Total.Chlorine..mg.l.","Iron..ppm.","Latittude","Longitude","Carbonate..ppm.","pH","Total.alkalinity..ppm.","Phosphate..ppm.","H2S..ppm.","Temp...C.","EC..mS.")

rownames(metadata.micronet) <- metadata.micronet$Sample_id

metadata.micronet <- as.data.frame(as.matrix(metadata.micronet))

metadata.micronet[11:ncol(metadata.micronet)] <- lapply(metadata.micronet[11:ncol(metadata.micronet)], as.numeric)

sapply(metadata.micronet,class)

## import rooted tree
rootedTree.rare <- ape::read.tree("rootedTree.rare.tree")
class(rootedTree.rare)

# Let's create a microtable object with more information
microeco.dataset.3 <- microtable$new(sample_table = metadata.micronet, otu_table = microeco.otu, tax_table = tax, phylo_tree = rootedTree.rare)

microeco.dataset.test <- microtable$new(sample_table = metadata.micronet, otu_table = microeco.otu,tax_table = tax, phylo_tree = rootedTree.rare)

#### Mantel's test for all the data ####

#microeco.dataset.2 <- clone(microeco.dataset)
microeco.dataset.3$tidy_dataset()

#unifrac default FALSE; whether UniFrac index should be calculated,binary default FALSE; TRUE is used for jaccard and unweighted unifrac;
#microeco.dataset.2$cal_betadiv(unifrac = TRUE, binary = TRUE)
microeco.dataset.3$cal_betadiv(unifrac = TRUE, binary = TRUE)
print(microeco.dataset.3$beta_diversity)
names(microeco.dataset.3$beta_diversity)

#default
microeco.dataset.test$cal_betadiv()
names(microeco.dataset.test$beta_diversity)
print(microeco.dataset.test$beta_diversity)

# first perform mantel test 
# Nitrate and Nitrites excluded as they were zero
# all other env variables were considered as indication of human activity so not included.
microeco.mantel <- trans_env$new(dataset = microeco.dataset.3, env_cols = 23:31)
microeco.mantel$datase

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

#### varpart

## Reducing the weight of rare species use hellinger transformation
otu <- as.matrix(as.data.frame(com.rarefied.min.int.exclude.1.3.rooted@otu_table))

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

otu.hellinger.transform <- decostand(otu, method = "hellinger")

varp.chosen.2 <- varpart (otu.hellinger.transform, ~  EC..mS., ~ Latittude, ~ Carbonate..ppm., ~ pH, data = metadata.varpart)

varp.chosen.2$part$fract

#get partition table
summary(varp.chosen.2)

## indicates testable fractions
varp.chosen.2$part$fract

# We may plot the results into Venn's diagram (argument digits influences number of decimal digits shown in the diagram, Xnames the displayed names of variables, and bg background color of the diagram fractions; see ?varpart for details):

pdf("varpart.2.pdf", width = 8)
plot (varp.chosen.2, digits = 2, Xnames = c('EC (mS)','Latittude (D)', 'Carbonate (ppm)', 'pH'), bg = c('turquoise','khaki3', 'greenyellow', "mediumslateblue"))
dev.off()

## As stated in the varpart fucntion: "Use function ‘rda’ to test significance of fractions of interest"

fractions.varpart.rda <- rda(otu.hellinger.transform ~ `H2S..ppm.` + `EC..mS.` + `pH` + Condition (`Temp...C.`), data =  metadata.varpart)


## ## The global model (fractions [a+b+c]):
anova(otu.hellinger.chosen.2)

## Conditional (partial) effect of Temperature
anova(fractions.varpart.rda)
###################################################################################

### Fig. 5A. Putative biotic interactions were estimated using co-occurrence network analysis. Modules of interaction are shown by circles sharing the same colour fill ####
### Fig. S7. Co-occurrence network analysis indicating ecological groups of interacting ASVs. ####
### Fig. 5B. Chord diagram illustrating major interactions between abundant ASVs clustered at genus level ####

######################################################################################

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
#photocol3 <- c("slateblue2", "khaki3", "darkseagreen1", "pink4", "cadetblue", "darkgoldenrod1")

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

### Chord Diagram

### all tables as data.frame first before combining into microtable object
microeco.otu <- t(as.matrix(as.data.frame(otu_table(com.rarefied.min.int.exclude.1.3.rooted))))
microeco.otu <- as.data.frame(microeco.otu)
sapply(microeco.otu,class)
class(microeco.otu)

tax <- (as.matrix(as.data.frame(tax_table(com.rarefied.min.int.exclude.1.3.rooted))))
tax <- as.data.frame(tax)

metadata <- as.data.frame(sample_data(com.rarefied.min.int.exclude.1.3.rooted)) %>% as_tibble()
colnames(metadata)
metadata.micronet <- metadata %>% select("Sample_id","Country","Region","Location.Code","Site","Location","Latittude","Longitude","Human.usage..Y.N.","Water.flow..Pool.Flowing.","Copper..mg.L.","Cyanuric.acid..mg.L.","Nitrate..ppm.","Nitrite..ppm.","Temp.adj...C.","Chloride..mg.L.","Chlorine.dioxide..mg.L.","Flouride..mg.L.","Hardness..mg.L.","Lead..mg.L.","MPS..mg.L.","QUAT.QAC..mg.L.","Total.Chlorine..mg.l.","Carbonate..ppm.","pH","Total.alkalinity..ppm.","Phosphate..ppm.","H2S..ppm.","Temp...C.","Iron..ppm.","EC..mS.")

rownames(metadata.micronet) <- metadata.micronet$Sample_id

metadata.micronet <- as.data.frame(as.matrix(metadata.micronet))

metadata.micronet[11:ncol(metadata.micronet)] <- lapply(metadata.micronet[11:ncol(metadata.micronet)], as.numeric)

sapply(metadata.micronet,class)

## import rooted tree
rootedTree.rare <- ape::read.tree("rootedTree.rare.tree")
class(rootedTree.rare)

# Let's create a microtable object with more information
microeco.dataset <- microtable$new(sample_table = metadata.micronet, otu_table = microeco.otu, tax_table = tax, phylo_tree = rootedTree.rare)

microeco.dataset$sample_table

saveRDS(microeco.dataset,"microeco.dataset.rds")

#microeco_transnetwork_spieceasi  <- trans_network$new(
#  dataset = microeco.dataset,
#  cor_method = NULL,taxa_level = "OTU")

#saveRDS(microeco_transnetwork_spieceasi,"microeco_transnetwork_spieceasi.rds")

# run on hpc
# with list providing the same parameters as used to create spieceasi mb plot with spieceasi package
#microeco_transnetwork_spieceasi$cal_network(network_method = "SpiecEasi", SpiecEasi_method = "mb", ...list(sel.criterion='stars' ,lambda.min.ratio=1e-2,nlambda=20, pulsar.params=list(rep.num=100, ncores =24)))

microeco_transnetwork_spieceasi <- readRDS("microeco_transnetwork_spieceasi.rds")

save(microeco_transnetwork_spieceasi, file = "microeco_transnetwork_spieceasi.RData")
zip("microeco_transnetwork_spieceas.1.zip", "microeco_transnetwork_spieceasi.RData")


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

# plot the correlation heatmap
pdf("microeco.spieceasi.correlation.modules.pdf")
microeco_transnetwork_spieceasi_correlation$plot_cor()
dev.off()

# default parameter represents using igraph plot.igraph function
microeco_transnetwork_spieceasi$plot_network(method = "igraph", layout = layout_with_kk,node_color = "module")

# use ggraph method; require ggraph package
# If ggraph is not installed; first install it with command: install.packages("ggraph")
microeco_transnetwork_spieceasi$plot_network(method = "ggraph", node_color = "Phylum")

pdf("microeco.spieceasi.ggraph.pdf", width = 18, height = 12)
microeco_transnetwork_spieceasi$plot_network(method = "ggraph", node_color = "module")
dev.off()

# For dynamic network plots :use networkD3 package method for the dynamic network visualization in R
# If networkD3 is not installed; first install it with command: 
install.packages("networkD3")
library(networkD3)
#microeco_transnetwork_spieceasi$plot_network(method = "networkD3", node_color = "module")
#microeco_transnetwork_spieceasi$plot_network(method = "networkD3", node_color = "Phylum")

## circos plot
## The function cal_sum_links can sum the links (edge) number from one taxa to another or within the same taxa. The function plot_sum_links is used to show the result from the function cal_sum_links. This is very useful to fast see how many nodes are connected between different taxa or within one taxa. In terms of ‘Phylum’ level in the tutorial, the function cal_sum_links() sum the linkages number from one Phylum to another Phylum or the linkages in the same Phylum. So the numbers along the outside of the circular plot represent how many edges or linkages are related with the Phylum. For example, in terms of Proteobacteria, there are roughly total 900 edges associated with the OTUs in Proteobacteria, in which roughly 200 edges connect both OTUs in Proteobacteria and roughly 150 edges connect the OTUs from Proteobacteria with the OTUs from Chloroflexi.
microeco_transnetwork_spieceasi$cal_sum_links(taxa_level = "Family")
microeco_transnetwork_spieceasi$res_sum_links_pos
#devtools::install_github("mattflor/chorddiag")

library(chorddiag)

color_values = RColorBrewer::brewer.pal(12, "Paired")


color_values_2 = c("#3B4D16","#950D3D", "#848F22", "#A7A3FE","#B2C9B2", "#8CCBD3", "#0DD7FD", "#99D584", "#0D8CA8","#DFC500","#2A7F72","#FEAB2E", "#4B00FD", "#26E0A9" ,"#CC003D")

microeco_transnetwork_spieceasi$plot_sum_links(plot_pos = TRUE, plot_num = 15, color_values = color_values_2, groupnameFontsize = 8, showTicks = FALSE, groupnamePadding = 40)

microeco_transnetwork_spieceasi$plot_sum_links(plot_pos = FALSE, plot_num = 15, color_values = color_values_2)

# From v1.2.0, method = "circlize" is available for conveniently saving the static plot
#if(!require("circlize")) install.packages("circlize")
library("circlize")
microeco_transnetwork_spieceasi$plot_sum_links(method = "circlize", transparency = 0.2, plot_num = 15, color_values = color_values_2, annotationTrackHeight = circlize::mm_h(c(5, 5)))


###################################################################################

### Fig. 6A,B,C. Null models were employed to predict the influence of various evolutionary drivers on community assembly. Net relatedness index (betaNRI) estimates are shown by biogeographic region for all taxa (A), versus photosynthetic (B) and chemoheterotrophic (C) taxa ####
### Fig. S8. Net relatedness index (betaNRI) measures of the mean phylogenetic distance to the nearest taxon in the community for ecological groups of bacteria by biogeographic region.

######################################################################################

microeco.nullmodel <- trans_nullmodel$new(microeco.dataset.3, env_cols = 23:31)
saveRDS(microeco.nullmodel, "microeco.nullmodel.rds")

##test
microeco.dataset.test.nullmodel <- trans_nullmodel$new(microeco.dataset.test, env_cols = 23:31)
microeco.dataset.test.nullmodel$cal_ses_betampd(runs = 1, abundance.weighted = TRUE, null.model = c("taxa.labels", "richness", "frequency", "sample.pool", "phylogeny.pool",
                                                                                                    "independentswap", "trialswap")[1])
data(sample_info_16S)
data(env_data_16S)

photosynthetic.nullmodel <- trans_nullmodel$new(photosynthetic, env_cols = 23:31)
saveRDS(photosynthetic.nullmodel, "photosynthetic.nullmodel.rds")

### create additional pbjects for each of the photosynthetic groups:

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

