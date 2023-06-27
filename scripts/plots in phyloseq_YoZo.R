## investigate the final 16S library


remotes::install_github("microbiome/microbiome")
library("microbiome")
library("vegan")
library("tidyr")
library("ggplot2")
library("phyloseq")


# ps.prune is the final phyloseq object from pipeline, can load it in here...
### or load in the raw data and re-assemble

ps.prune<-readRDS("output/TSCC_output/ps.prune.RDS")


######### if reassembling
# read in metadata and re-format so it has the str necessary
M<-read.csv("output/TSCC_output/sam_data.csv")

# set row names
row.names(M)<-M$X 

# remove junk column
M<-M[-1]

# format metadata (have to redo this if loading back in)
make.fac<-c("Year", "Location", "Lake", "Sample.type", "Organism", "Functional.group", "sample_control", "Miseq.ANL")

M[make.fac]<-lapply(M[make.fac], factor) # make all these factors

sample_data(ps.prune)<-M


#### assemble the phyloseq object from raw files
PS.fin<-
  read_phyloseq(
    otu.file = "output/TSCC_output/otu_table.csv",
    taxonomy.file = "output/TSCC_output/tax_table.csv",
    metadata.file = "output/TSCC_output/sam_data.csv",
    sep = ","
  )

#replace sample data with reformatted metadata
sample_data(PS.fin)<-M

# check that levels exist
levels(get_variable(PS.fin, "Lake"))
levels(get_variable(PS.fin, "Organism"))




################ Now on to plotting

PS.fin<-ps.prune

###### some richness plots
#richness by Lake
richness.plot_Lake<-plot_richness(PS.fin, x="Lake", measures=c("Observed", "Shannon"))  + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
richness.plot_Lake
dev.copy(pdf, "figures/richness.plot.location.pdf", height=4, width=10)
dev.off() 

#richness by organism (or sample type i.e, water)
richness.plot_organism<-plot_richness(PS.fin, x="Organism", measures=c("Observed", "Shannon")) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
richness.plot_organism
dev.copy(pdf, "figures/richness.plot.organism.pdf", height=4, width=12)
dev.off() 

#richness by Time.point
richness.plot_Time<-plot_richness(PS.fin, x="Year", measures=c("Observed", "Shannon")) + theme_bw()
richness.plot_Time
dev.copy(pdf, "figures/richness.plot.Year.pdf", height=4, width=10)
dev.off() 
###### 



##### Read counts and rarefaction curves
# Make a data frame with a column for the read counts of each sample
sample_sum_df<-as.data.frame(sample_sums(PS.fin))
colnames(sample_sum_df)<-"read.sum"
sample_sum_df$sampleNames<-rownames(sample_sum_df)

# merge in the reads
run.metaD<-merge(M, sample_sum_df, by="sampleNames", all.y=TRUE)
write.csv(run.metaD, "output/run.metaD.final.csv")


# Histogram of sample read counts
hist.depth<-ggplot(sample_sum_df, aes(x = read.sum)) + 
  geom_histogram(color = "black", fill = "gold2", binwidth = 2500) +
  ggtitle("Distribution of sample sequencing depth") + 
  xlab("Read counts") +
  theme(axis.title.y = element_blank()) + theme_classic() + geom_vline(xintercept=1000, lty=2)

hist.depth
dev.copy(pdf, "figures/hist.depth.pdf", height=4, width=5)
dev.off() 

##
### need to do some some smithing to get to work

OTU <- otu_table(PS.fin)
class(OTU) <- "matrix" # as.matrix() will do nothing
## you get a warning here, but this is what we need to have

OTU <- t(OTU) # transpose observations to rows

# version 1
pdf(file="figures/rare.raw.pdf", height=6, width=10)
rarecurve(OTU, step=50, cex=0.5, xlim=c(0,10000), ylim=(c(0,500)), label=FALSE)
#abline(v = 5000, lty = "dotted", col="red", lwd=2)
dev.off() 



#### more plots
pdf(file="figures/read.by.species.pdf", height=4, width=12)
boxplot(run.metaD$read.sum~run.metaD$Organism, cex.axis=0.5, cex.lab=0.8)
dev.off() 

pdf(file="figures/read.by.sample.pdf", height=6, width=5)
boxplot(run.metaD$read.sum~run.metaD$Sample.type, cex.axis=0.5, cex.lab=0.8)
dev.off() 

pdf(file="figures/reads.sample.pdf", height=4, width=7)
ggplot(run.metaD, aes(x=sample_control, y=read.sum, color=Organism)) + geom_boxplot()
dev.off() 

pdf(file="figures/reads.year.pdf", height=4, width=7)
ggplot(run.metaD, aes(x=Year, y=read.sum, color=Year)) + geom_boxplot()
dev.off() 
##### 





### PCoA
############ PCoA
sample_variables(PS.fin)

# make colors for sites
library(RColorBrewer)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))


###
T1<- subset_samples(PS.fin, Year=="2021" | Year=="2022")
ORD <- ordinate(T1, method='MDS', distance='bray')

NMDS.ord.organism<-plot_ordination(
  physeq = T1,                                                   
  ordination = ORD) +                                                
  geom_point(aes(color = Organism, shape=Year), size = 3) +    
  stat_ellipse(level=0.9, linetype = 2, aes(color=Organism)) +
  #scale_color_brewer(palette = "Dark2") +
  #scale_fill_brewer(palette = "Dark2") + 
  ggtitle("2021-22 Zoop") +
  theme_classic()   

NMDS.ord.organism
dev.copy(pdf, "figures/NMDS.ord.zoop.pdf", height=6, width=9)
dev.off() 


########## 

NMDS.ord.functional<-plot_ordination(
  physeq = T1,                                                   
  ordination = ORD) +                                                
  geom_point(aes(color = Functional.group), size = 3) +    
  stat_ellipse(level=0.9, linetype = 2, aes(color=Functional.group)) +
  #scale_color_brewer(palette = "Dark2") +
  #scale_fill_brewer(palette = "Dark2") + ggtitle("Time1") +
  theme_classic()   

NMDS.ord.functional
dev.copy(pdf, "figures/NMDS.ord.func.group.pdf", height=6, width=9)
dev.off() 

########## 

NMDS.ord.types<-plot_ordination(
  physeq = PS.fin,                                                   
  ordination = ORD) +                                                
  geom_point(aes(color = Sample.type, shape=Year), size = 3) +    
  stat_ellipse(level=0.9, linetype = 2, aes(color=Sample.type)) +
  #scale_color_brewer(palette = "Dark2") +
  #scale_fill_brewer(palette = "Dark2") + ggtitle("Time1") +
  theme_classic()   

NMDS.ord.types
dev.copy(pdf, "figures/NMDS.ord.types.pdf", height=6, width=9)
dev.off() 

########## 

NMDS.ord.lake<-plot_ordination(
  physeq = PS.fin,                                                   
  ordination = ORD) +                                                
  geom_point(aes(color = Lake, shape=Year), size = 3) +    
  #scale_color_brewer(palette = "Dark2") +
  #scale_fill_brewer(palette = "Dark2") + ggtitle("Time1") +
  theme_classic()   

NMDS.ord.lake
dev.copy(pdf, "figures/NMDS.ord.lake.pdf", height=6, width=9)
dev.off() 

### stacked bar plot
pdf(file="figures/stack.bar1.pdf", height=5, width=12)
plot_bar(PS.fin, x="Sample.type", fill="Phylum")
dev.off() 


## other
pdf(file="figures/stack.bar1.pdf", height=5, width=12)
p = plot_bar(PS.fin, "Phylum", fill="Phylum", facet_grid=Year~Organism)
dev.off() 