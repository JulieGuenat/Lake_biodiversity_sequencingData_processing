###################### DATA CLEANING #######################
################# Julie Guenat ############################

# Packages ####
library(metabaR)
library(tidyverse)
library(reshape2)
library(RColorBrewer)
library(scales)
library(cowplot)

#set directory ####

setwd("C:/Users/jguenat/Documents/PhD_Analyses/Cleaning_data")

# Import data ####
## raw OBITools output ####
euka03_P1_DNA_raw <- read.table("Obitools_output/14.tab_P1DNA_Euka03.txt", header =T, sep="\t")
euka03_P1_DNA_raw <- euka03_P1_DNA_raw %>% select("id", "cluster_weight", "count", starts_with("sample"), 
                                                                     "obitag_bestid", "obitag_bestmatch",  "obitag_match_count", "obitag_rank",
                                                                     "obitag_similarity_method", "scientific_name", "taxid", "taxonomic_path",
                                                                     "sequence")

### read table ####
reads_P1_euka03_DNA <-t(euka03_P1_DNA_raw[,grep("sample\\.", colnames(euka03_P1_DNA_raw))])
rownames(reads_P1_euka03_DNA) = gsub("sample\\.", "", rownames(reads_P1_euka03_DNA))
colnames(reads_P1_euka03_DNA) = euka03_P1_DNA_raw$id
rownames(reads_P1_euka03_DNA) = gsub("\\.", "_", rownames(reads_P1_euka03_DNA))
reads_P1_euka03_DNA <- as.matrix(reads_P1_euka03_DNA)
dim(reads_P1_euka03_DNA)
#view(reads_P1_euka03_DNA)

### motus table ####
motus_P1_euka03_DNA = euka03_P1_DNA_raw[,grep("sample\\.", colnames(euka03_P1_DNA_raw), invert = T)]
rownames(motus_P1_euka03_DNA)=motus_P1_euka03_DNA$id
#View(motus_P1_euka03_DNA)

## Samples dataset
sample_P1<-read.csv("sample_files/Samplefile_P1.csv", header= T, sep= ";")
#View(sample_P1)
#We relabel samples file rownames to be the Sample ID so we can merge by Sample_ID with the reads table.
sample_P1 = sample_P1[order(sample_P1$Sample_ID),]
sample_P1 = sample_P1[!duplicated(sample_P1$Sample_ID),]
sample_P1$Sample_ID <- gsub("-", "_", sample_P1$Sample_ID)
rownames(sample_P1)=sample_P1$Sample_ID
#str(samples)
samp_P1 = sample_P1
#[match(rownames(reads_P1_euka03_DNA), rownames(sample_P1)),]
#View(samp_P1)

## NGS filterfile ####
ngs_P1 <- read.table("ngs_filters/P1DNA_Euka03_ngsfilter.txt", sep="\t")
ngs_P1$V2 <- gsub("-", "_", ngs_P1$V2)

## PCRs GPS ####
pcr_gps_P1 <- read.csv("PCR_GPS/PCR_GPS_P1.csv", header = T, sep= ";") 
pcr_gps_P1$Sample_ID <- gsub("-", "_", pcr_gps_P1$Sample_ID)

pcrs_P1 <- merge(samp_P1, pcr_gps_P1, by.y="Sample_ID")

### merge PCR and GPS ####
pcrs_P1 <- merge(pcrs_P1, ngs_P1[,c(2:5)], by.x<-"Sample_ID", by.y = "V2", all=T)

rownames(pcrs_P1) <- pcrs_P1$Sample_ID
names(pcrs_P1)[names(pcrs_P1) == "Sample_ID"] <- "sample_id"
names(pcrs_P1)[names(pcrs_P1) == "V3"] <- "tags"
names(pcrs_P1)[names(pcrs_P1) == "V4"] <- "primer_fwd"
names(pcrs_P1)[names(pcrs_P1) == "V5"] <- "primer_rev"

pcrs_P1 <- pcrs_P1 %>% separate(tags, c("tag_fwd", "tag_rev"), ":")
pcrs_P1 <- as.data.frame(pcrs_P1)
pcrs_P1 <- subset(pcrs_P1, !is.na(pcrs_P1$type))
row.names(pcrs_P1)<- pcrs_P1$sample_id
#view(pcrs_P1)

pcrs_P1 <- pcrs_P1[,-1]
names(pcrs_P1)[names(pcrs_P1) == "Sample"] <- "sample_id"

## Metadata ####
metadata<-read.csv("Samples_info_Planaqua.csv", header= T, sep=";")

#For the sample info, we need to select the samples corresponding at the library. 
metadata_P1 <- merge(metadata, samp_P1, by.x = "Samples_ID", by.y = "Sample_pla")

#DOUBLE CHECK HERE
metadata_P1 = metadata_P1[!duplicated(metadata_P1[,c('Sample')]),]
rownames(metadata_P1)=metadata_P1$Sample
metadata_P1 = as.data.frame(metadata_P1)
#View(metadata_P1)

# METABAR LIST ####
P1_euka03_DNA <- metabarlist_generator(reads = reads_P1_euka03_DNA, motus = motus_P1_euka03_DNA, pcrs = pcrs_P1, samples = metadata_P1)
summary_metabarlist(P1_euka03_DNA)

# we can remove all the useless dataset from the env:
rm(euka03_P1_DNA_raw, metadata, metadata_P1, motus_P1_euka03_DNA, ngs_P1, pcr_gps_P1, pcrs_P1, samp_P1, sample_P1, reads_P1_euka03_DNA)

# OVERVIEW OF CONTAMINANTS ####
# Compute the number of reads per pcr                                   
P1_euka03_DNA$pcrs$nb_reads <- rowSums(P1_euka03_DNA$reads)
# Compute the number of motus per pcr
P1_euka03_DNA$pcrs$nb_motus <- rowSums(P1_euka03_DNA$reads>0)

# Create an input table (named check1) for ggplot of 3 columns: 
#  (i) control type 
#  (ii) a vector indicated whether it corresponds to nb_reads or nb_motus, 
#  (iii) the corresponding values.

check1 <- reshape2::melt(P1_euka03_DNA$pcrs[,c("control_type", "nb_reads", "nb_motus")])

p1<-ggplot(data <- check1, aes(x=control_type, y=value, color=control_type)) + 
  geom_boxplot() + theme_bw() + 
  geom_jitter(alpha=0.2) + 
  scale_color_manual(values = c("brown", "red", "cyan4","pink"), na.value = "darkgrey") +
  facet_wrap(~variable, scales = "free_y") + 
  theme(axis.text.x = element_text(angle=45, h=1))
p1
#ggsave(filename="1_P1_Euka03_DNA_barplot.contaminants.pdf", plot=p1, width = 297, height = 210, units = "mm", device="pdf", path="Output_metabar_data/P1_euka03_DNA/figures_cleaning")

p2<-ggplot(P1_euka03_DNA$pcrs, aes(x=nb_reads, y=nb_motus, color = control_type)) + 
  geom_point() + theme_bw() + 
  scale_y_log10() + scale_x_log10() + 
  scale_color_manual(values = c("brown", "red", "cyan4","pink"), na.value = "darkgrey")
p2
#ggsave(filename="2_P1_Euka03_DNA_plot.motus.vs.reads.pdf", plot=p2, width = 297, height = 210, units = "mm", device="pdf", path="Output_metabar_data/P1_euka03_DNA/figures_cleaning")

p3<-ggpcrplate(P1_euka03_DNA, FUN = function(m){rowSums(m$reads)}, legend_title = "# of reads per PCR")
p3
#ggsave(filename="3_P1_Euka03_DNA_pcr.plots.pdf", plot=p3, width = 297, height = 210, units = "mm", device="pdf", path="Output_metabar_data/P1_euka03_DNA/figures_cleaning")

tag.list <- as.character(unique(P1_euka03_DNA$pcrs$tag_rev))
p4<-ggpcrtag(P1_euka03_DNA, 
             legend_title = "# of reads per PCR", 
             FUN = function(m) {rowSums(m$reads)},
             taglist = tag.list) 
p4
#ggsave(filename="4_P1_Euka03_DNA_pcr.plots.by.tags.pdf", plot=p4, width = 297, height = 210, units = "mm", device="pdf", path="Output_metabar_data/P1_euka03_DNA/figures_cleaning")

# FLAGGING SPURIOUS SIGNAL ####
## Extraction control ####
P1_euka03_DNA <- contaslayer(P1_euka03_DNA, 
                             control_types = c("extraction"),
                             output_col = "not_an_extraction_conta")

table(P1_euka03_DNA$motus$not_an_extraction_conta)

id <- !P1_euka03_DNA$motus$not_an_extraction_conta
max.conta <- rownames(P1_euka03_DNA$motus[id,])[which.max(P1_euka03_DNA$motus[id, "count"])]

#... and its distribution and relative abundance in each pcr
# the "#reads of most abundant contaminat" is probably refferring to the % of reads that are labelled contaminant (check function above)
ggpcrplate(P1_euka03_DNA, legend_title = "#reads of most \nabundant contaminant",
           FUN = function(m) {m$reads[, max.conta]/rowSums(m$reads)})

# Compute relative abundance of all pcr contaminants together
a <- data.frame(conta.relab = rowSums(P1_euka03_DNA$reads[,!P1_euka03_DNA$motus$not_an_extraction_conta]) /
                  rowSums(P1_euka03_DNA$reads))
# Add information on control types
a$control_type <- P1_euka03_DNA$pcrs$control_type[match(rownames(a), rownames(P1_euka03_DNA$pcrs))]

p6<-ggplot(a, aes(x=control_type, y=conta.relab, color=control_type)) +
  geom_boxplot() + geom_jitter(alpha=0.5) +
  scale_color_manual(values = c("brown", "red", "cyan4","pink"), na.value = "darkgrey") +
  labs(x=NULL, y="Prop. Reads (log10)") +
  theme_bw() +
  scale_y_log10()
p6
#ggsave(filename="6_P1_Euka03_DNA_boxplot_contamination_extraction.pdf", plot=p6, width = 297, height = 210, units = "mm", device="pdf", path="Output_metabar_data/P1_euka03_DNA/figures_cleaning")

# flag pcrs with total contaminant relative abundance > 10% of reads)
P1_euka03_DNA$pcrs$low_contamination_level <- 
  ifelse(a$conta.relab[match(rownames(P1_euka03_DNA$pcrs), rownames(a))]>1e-1,  F, T)
# Proportion of potentially functional (TRUE) vs. failed (FALSE) pcrs
# (controls included) based on this criterion
table(P1_euka03_DNA$pcrs$low_contamination_level) / nrow(P1_euka03_DNA$pcrs)

## PCR control ####
P1_euka03_DNA <- contaslayer(P1_euka03_DNA, 
                             control_types = c("pcr"),
                             output_col = "not_an_pcr_conta")
table(P1_euka03_DNA$motus$not_an_pcr_conta)

# Identify the most common contaminant
# get contaminant ids
id <- !P1_euka03_DNA$motus$not_an_pcr_conta
max.conta <- rownames(P1_euka03_DNA$motus[id,])[which.max(P1_euka03_DNA$motus[id, "count"])]

#... and its distribution and relative abundance in each pcr
# the "#reads of most abundant contaminat" is probably refferring to the % of reads that are labelled contaminant (check function above)
ggpcrplate(P1_euka03_DNA, legend_title = "#reads of most \nabundant contaminant",
           FUN = function(m) {m$reads[, max.conta]/rowSums(m$reads)})

# Compute relative abundance of all pcr contaminants together
a <- data.frame(conta.relab = rowSums(P1_euka03_DNA$reads[,!P1_euka03_DNA$motus$not_an_pcr_conta]) /
                  rowSums(P1_euka03_DNA$reads))
# Add information on control types
a$control_type <- P1_euka03_DNA$pcrs$control_type[match(rownames(a), rownames(P1_euka03_DNA$pcrs))]

p7<-ggplot(a, aes(x=control_type, y=conta.relab, color=control_type)) +
  geom_boxplot() + geom_jitter(alpha=0.5) +
  scale_color_manual(values = c("brown", "red", "cyan4","pink"), na.value = "darkgrey") +
  labs(x=NULL, y="Prop. Reads (log10)") +
  theme_bw() +
  scale_y_log10()

p7
#ggsave(filename="7_P1_Euka03_DNA_boxplot_contamination_pcr.pdf", plot=p7, width = 297, height = 210, units = "mm", device="pdf", path="Output_metabar_data/P1_euka03_DNA/figures_cleaning")

# flag pcrs with total contaminant relative abundance > 10% of reads)
P1_euka03_DNA$pcrs$low_contamination_level <- 
  ifelse(a$conta.relab[match(rownames(P1_euka03_DNA$pcrs), rownames(a))]>1e-1,  F, T)
# Proportion of potentially functional (TRUE) vs. failed (FALSE) pcrs
# (controls included) based on this criterion
table(P1_euka03_DNA$pcrs$low_contamination_level) / nrow(P1_euka03_DNA$pcrs)

## Sequencing control ####
P1_euka03_DNA <- contaslayer(P1_euka03_DNA, 
                             control_types = c("sequencing"),
                             output_col = "not_an_seq_conta")
table(P1_euka03_DNA$motus$not_an_seq_conta)

# Identify the most common contaminant
# get contaminant ids
id <- !P1_euka03_DNA$motus$not_an_seq_conta
max.conta <- rownames(P1_euka03_DNA$motus[id,])[which.max(P1_euka03_DNA$motus[id, "count"])]

#... and its distribution and relative abundance in each pcr
# the "#reads of most abundant contaminat" is probably refferring to the % of reads that are labelled contaminant (check function above)
ggpcrplate(P1_euka03_DNA, legend_title = "#reads of most \nabundant contaminant",
           FUN = function(m) {m$reads[, max.conta]/rowSums(m$reads)})

# Compute relative abundance of all pcr contaminants together
a <- data.frame(conta.relab = rowSums(P1_euka03_DNA$reads[,!P1_euka03_DNA$motus$not_an_seq_conta]) /
                  rowSums(P1_euka03_DNA$reads))
# Add information on control types
a$control_type <- P1_euka03_DNA$pcrs$control_type[match(rownames(a), rownames(P1_euka03_DNA$pcrs))]

p8<-ggplot(a, aes(x=control_type, y=conta.relab, color=control_type)) +
  geom_boxplot() + geom_jitter(alpha=0.5) +
  scale_color_manual(values = c("brown", "red", "cyan4","pink"), na.value = "darkgrey") +
  labs(x=NULL, y="Prop. Reads (log10)") +
  theme_bw() +
  scale_y_log10()
p8
#ggsave(filename="8_P1_Euka03_DNA_boxplot_contamination_sequencing.pdf", plot=p8, width = 297, height = 210, units = "mm", device="pdf", path="Output_metabar_data/P1_euka03_DNA/figures_cleaning")

# flag pcrs with total contaminant relative abundance > 10% of reads)
P1_euka03_DNA$pcrs$low_contamination_level <- 
  ifelse(a$conta.relab[match(rownames(P1_euka03_DNA$pcrs), rownames(a))]>1e-1,  F, T)
# Proportion of potentially functional (TRUE) vs. failed (FALSE) pcrs
# (controls included) based on this criterion
table(P1_euka03_DNA$pcrs$low_contamination_level) / nrow(P1_euka03_DNA$pcrs)

## Non-target MOTUs ####
#Flag MOTUs corresponding to target (TRUE) vs. non-target (FALSE) taxa
P1_euka03_DNA$motus$target_taxon <- grepl("Eukaryota", P1_euka03_DNA$motus$taxonomic_path)

# Proportion of each of these over total number of MOTUs
table(P1_euka03_DNA$motus$target_taxon) / nrow(P1_euka03_DNA$motus)

# Plot the unweighted distribution of MOTUs similarity scores 
a <- 
  ggplot(P1_euka03_DNA$motus, aes(x=obitag_bestid)) + 
  geom_histogram(color="grey", fill="white", bins=20) + 
  geom_vline(xintercept = 0.80, col="orange", lty=2) + 
  theme_bw() + 
  theme(panel.grid = element_blank()) + 
  labs(x="% similarity against best match", y="# MOTUs")

# Same for the weighted distribution
b <- 
  ggplot(P1_euka03_DNA$motus, 
         aes(x=obitag_bestid, y = ..count.., weight = count)) + 
  geom_histogram(color="grey", fill="white", bins=20) + 
  geom_vline(xintercept = 0.80, col="orange", lty=2) + 
  theme_bw() + 
  theme(panel.grid = element_blank()) + 
  labs(x="% similarity against best match", y="# Reads")

# Combine plots into one
p9<-ggdraw() + 
  draw_plot(a, x=0, y=0, width = 0.5) + 
  draw_plot(b, x=0.5, y=0, width = 0.5)
p9
#ggsave(filename="9_P1_Euka03_DNA_barplot_db_best_match.pdf", plot=p9, width = 297, height = 210, units = "mm", device="pdf", path="Output_metabar_data/P1_euka03_DNA/figures_cleaning")

# Flag not degraded (TRUE) vs. potentially degraded sequences (FALSE)
P1_euka03_DNA$motus$not_degraded <-
  ifelse(P1_euka03_DNA$motus$obitag_bestid < 0.80, F, T)

# Proportion of each of these over total number of MOTUs
table(P1_euka03_DNA$motus$not_degraded) / nrow(P1_euka03_DNA$motus)

## Detecting PCR outliers ####
#sequencing depth

ggplot(P1_euka03_DNA$pcrs, aes(nb_reads)) +
  geom_histogram(bins=40, color="grey", fill="white") + 
  geom_vline(xintercept = 1e2, lty=2, color="orange") + # threshold
  scale_x_log10() + 
  labs(x="# Reads (with all MOTUs and PCRs)", 
       y="# PCRs") +
  theme_bw() + 
  theme(panel.grid = element_blank())

# Flag pcrs with an acceptable sequencing depth (TRUE) or inacceptable one (FALSE)
P1_euka03_DNA$pcrs$seqdepth_ok <- ifelse(P1_euka03_DNA$pcrs$nb_reads < 100, F, T)

# Proportion of each of these over total number of pcrs, control excluded
table(P1_euka03_DNA$pcrs$seqdepth_ok[P1_euka03_DNA$pcrs$type=="sample"]) /
  nrow(P1_euka03_DNA$pcrs[P1_euka03_DNA$pcrs$type=="sample",])

# second way == reproducibility between replicates
# Subsetting the metabarlist
P1_euka03_DNA_sub <- subset_metabarlist(P1_euka03_DNA,
                                        table = "pcrs",
                                        indices = P1_euka03_DNA$pcrs$nb_reads>0 & (
                                          is.na(P1_euka03_DNA$pcrs$control_type) |
                                            P1_euka03_DNA$pcrs$control_type=="positive"))

# First visualization
comp1 = pcr_within_between(P1_euka03_DNA_sub)
p10<- check_pcr_thresh(comp1)
p10
#ggsave(filename="10_P1_Euka03_DNA_density_plot_withinbetween_pcr_rep.pdf", plot=p10, width = 297, height = 210, units = "mm", device="pdf", path="Output_metabar_data/P1_euka03_DNA/figures_cleaning")

# Flagging
P1_euka03_DNA_sub <- pcrslayer(P1_euka03_DNA_sub, output_col = "replicating_pcr", plot = F)

# Proportion of replicating pcrs (TRUE)
table(P1_euka03_DNA_sub$pcrs$replicating_pcr) /
  nrow(P1_euka03_DNA_sub$pcrs)

# Intersection with the sequencing depth criterion
table(P1_euka03_DNA_sub$pcrs$seqdepth_ok,
      P1_euka03_DNA_sub$pcrs$replicating_pcr)

# Distinguish between pcrs obtained from samples from positive controls
mds = check_pcr_repl(P1_euka03_DNA_sub,
                     groups = P1_euka03_DNA_sub$pcrs$type,
                     funcpcr = P1_euka03_DNA_sub$pcrs$replicating_pcr)
mds + labs(color = "pcr type") + scale_color_manual(values = c("cyan4", "gray"))

P1_euka03_DNA$pcrs$replicating_pcr <- NA
P1_euka03_DNA$pcrs[rownames(P1_euka03_DNA_sub$pcrs),"replicating_pcr"] <- P1_euka03_DNA_sub$pcrs$replicating_pcr

## lowering tag jump ####
# Define a vector of thresholds to test
thresholds <- c(0,1e-4,1e-3, 1e-2, 3e-2, 5e-2)

# Run the tests and stores the results in a list
tests <- lapply(thresholds, function(x) tagjumpslayer(P1_euka03_DNA,x))
names(tests) <- paste("t_", thresholds, sep="")

# Format the data for ggplot with amount of reads at each threshold
tmp <- melt(as.matrix(do.call("rbind", lapply(tests, function(x) rowSums(x$reads)))))
colnames(tmp) <- c("threshold", "sample", "abundance")

# Add richness in MOTUs at each threshold
tmp$richness <-
  melt(as.matrix(do.call("rbind", lapply(tests, function(x) {
    rowSums(x$reads > 0)
  }))))$value

# Add control type information on pcrs and make data curation threshold numeric
tmp$controls <- P1_euka03_DNA$pcrs$control_type[match(tmp$sample, rownames(P1_euka03_DNA$pcrs))]
tmp$threshold <- as.numeric(gsub("t_", "", tmp$threshold))

# New table formatting for ggplot
tmp2 <- melt(tmp, id.vars=colnames(tmp)[-grep("abundance|richness", colnames(tmp))])

p11<- ggplot(tmp2, aes(x=as.factor(threshold), y=value)) +
  geom_boxplot(color="grey40") +
  geom_vline(xintercept = which(levels(as.factor(tmp2$threshold)) == "0.01"), col="orange", lty=2) +
  geom_jitter(aes(color=controls), width = 0.2, alpha=0.5) +
  scale_color_manual(values = c("brown", "red", "cyan4","pink"), na.value = "darkgrey") +
  facet_wrap(~variable+controls, scale="free_y", ncol=5) +
  theme_bw() +
  scale_y_log10() +
  labs(x="MOTU pcr : total abundance filtering threshold", y="# Reads/MOTUs") +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        axis.text.x = element_text(angle=40, h=1),
        legend.position = "none")
p11

#ggsave(filename="11_P1_Euka03_DNA_tagjump.pdf", plot=p11, width = 297, height = 210, units = "mm", device="pdf", path="Output_metabar_data/P1_euka03_DNA/figures_cleaning")

# Create a table of MOTUs quality criteria
# noise is identified as FALSE in soil_euk, the "!" transforms it to TRUE
motus.qual <- !P1_euka03_DNA$motus[,c("not_an_extraction_conta", "not_an_pcr_conta", "not_an_seq_conta", "target_taxon", "not_degraded")]

colnames(motus.qual) <- c("extraction_conta", "neg_conta", "seq_conta", "untargeted_taxon", "degraded_seq")

# Proportion of MOTUs potentially artifactual (TRUE) based on the criteria used
prop.table(table(apply(motus.qual, 1, sum) > 0))

# Corresponding proportion of artifactual reads (TRUE)
prop.table(xtabs(P1_euka03_DNA$motus$count~apply(motus.qual, 1, sum) > 0))

# Proportion of MOTUs and reads potentially artifactual for each criterion
apply(motus.qual, 2, sum) / nrow(motus.qual)
apply(motus.qual, 2, function(x) sum(P1_euka03_DNA$motus$count[x])/sum(P1_euka03_DNA$motus$count))

tmp.motus <-
  apply(sapply(1:ncol(motus.qual), function(x) {
    ifelse(motus.qual[,x]==T, colnames(motus.qual)[x], NA)}), 1, function(x) {
      paste(sort(unique(x)), collapse = "|")
    })
tmp.motus <- as.data.frame(gsub("^$", "not_artefactual", tmp.motus))
colnames(tmp.motus) <-  "artefact_type"

p12<-ggplot(tmp.motus, aes(x=1, fill=artefact_type)) +
  geom_bar() +  xlim(0, 2) +
  labs(fill="Artifact type") +
  coord_polar(theta="y") + theme_void() +
  scale_fill_brewer(palette = "Set3") +
  theme(legend.direction = "vertical") +
  ggtitle("MOTUs artefacts overview")

p12
ggsave(filename="12_P1_Euka03_DNA_MOTUs_Artefacts_overview.pdf", plot=p12, width = 297, height = 210, units = "mm", device="pdf", path="Output_metabar_data/P1_euka03_DNA/figures_cleaning")

# Create a table of pcrs quality criteria
# noise is identified as FALSE in soil_euk, the "!" transforms it to TRUE
pcrs.qual <- !P1_euka03_DNA$pcrs[,c("low_contamination_level", "seqdepth_ok")]
colnames(pcrs.qual) <- c("high_contamination_level", "low_seqdepth")

# Proportion of pcrs potentially artifactual (TRUE) based on the criteria used
# excluding controls
prop.table(table(apply(pcrs.qual[P1_euka03_DNA$pcrs$type=="sample",], 1, sum) > 0))

# Proportion of MOTUs and reads potentially artifactual for each criterion
apply(pcrs.qual[P1_euka03_DNA$pcrs$type=="sample",], 2, sum) / nrow(pcrs.qual[P1_euka03_DNA$pcrs$type=="sample",])

tmp.pcrs <-
  apply(sapply(1:ncol(pcrs.qual), function(x) {
    ifelse(pcrs.qual[P1_euka03_DNA$pcrs$type=="sample",x]==T,
           colnames(pcrs.qual)[x], NA)}), 1, function(x) {
             paste(sort(unique(x)), collapse = "|")
           })
tmp.pcrs <- as.data.frame(gsub("^$", "not_artefactual", tmp.pcrs))

colnames(tmp.pcrs) <- "artefact_type"

ggplot(tmp.pcrs, aes(x=1, fill=artefact_type)) +
  geom_bar() +  xlim(0, 2) +
  labs(fill="Artifact type") +
  coord_polar(theta="y") + theme_void() +
  scale_fill_brewer(palette = "Set3") +
  theme(legend.direction = "vertical") +
  ggtitle("PCR artefacts overview")

# DATA CLEANING ####
## Removing spurious signal ####
# Use tag-jump corrected metabarlist with the threshold identified above
tmp <- tests[["t_0.01"]]

# Subset on MOTUs: we keep motus that are defined as TRUE following the
# three criteria below (sum of three TRUE is equal to 3 with the rowSums function)
tmp <- subset_metabarlist(tmp, "motus",
                          indices = rowSums(tmp$motus[,c("not_an_extraction_conta", "not_an_pcr_conta", "not_an_seq_conta", "target_taxon", "not_degraded")]) == 5)
summary_metabarlist(tmp)

# Same with PCRs
#I have NA in the column "low_contamination_level" and in "replicating_pcr"
tmp$pcrs$low_contamination_level[is.na(tmp$pcrs$low_contamination_level)] <- FALSE
tmp$pcrs$replicating_pcr[is.na(tmp$pcrs$replicating_pcr)] <- FALSE

# Subset on pcrs and keep only controls
P1_euka03_DNA_clean <- subset_metabarlist(tmp, "pcrs",
                                          indices = rowSums(tmp$pcrs[,c("low_contamination_level",
                                                                        "seqdepth_ok", "replicating_pcr")]) == 3 &
                                            tmp$pcrs$type == "sample")
summary_metabarlist(P1_euka03_DNA_clean)

#Checking if there is still some empty PCR
if(sum(colSums(P1_euka03_DNA_clean$reads)==0)>0){print("empty motus present")}
if(sum(rowSums(P1_euka03_DNA_clean$reads)==0)>0){print("empty pcrs present")}

#Since we have now removed certain MOTUs or reduced their read counts, we need to update some parameters in the metabarlist (e.g. counts, etc.).
P1_euka03_DNA_clean$motus$count = colSums(P1_euka03_DNA_clean$reads)
P1_euka03_DNA_clean$pcrs$nb_reads_postmetabaR = rowSums(P1_euka03_DNA_clean$reads)
P1_euka03_DNA_clean$pcrs$nb_motus_postmetabaR = rowSums(ifelse(P1_euka03_DNA_clean$reads>0, T, F))

#We can now compare some basic characteristics before and after data curation with metabaR.
check <- melt(P1_euka03_DNA_clean$pcrs[,c("nb_reads", "nb_reads_postmetabaR",
                                          "nb_motus", "nb_motus_postmetabaR")])
check$type <- ifelse(grepl("motus", check$variable), "richness", "abundance")

p13<- ggplot(data = check, aes(x = variable, y = value)) +
  geom_boxplot( color = "darkgrey") +
  geom_jitter(alpha=0.1, color = "darkgrey") +
  theme_bw() +
  facet_wrap(~type, scales = "free", ncol = 5) +
  theme(axis.text.x = element_text(angle=45, h=1))
p13

#ggsave(filename="13_P1_Euka03_DNA_effects_cleaning_on_samples.pdf", plot=p13, width = 297, height = 210, units = "mm", device="pdf", path="Output_metabar_data/P1_euka03_DNA/figures_cleaning")


# DATA AGGREGATION ####
P1_euka03_DNA_agg_mean <- aggregate_pcrs(P1_euka03_DNA_clean, FUN = FUN_agg_pcrs_mean)
summary_metabarlist(P1_euka03_DNA_agg_mean)

write.csv(P1_euka03_DNA_agg_mean$motus,"Output_metabar_data/P1_euka03_DNA/files/P1_euka03_DNA_motus_final.csv", row.names = T)
write.csv(P1_euka03_DNA_agg_mean$pcrs,"Output_metabar_data/P1_euka03_DNA/files/P1_euka03_DNA_pcrs_final.csv", row.names = F)
write.csv(P1_euka03_DNA_agg_mean$samples,"Output_metabar_data/P1_euka03_DNA/files/P1_euka03_DNA_samples_final.csv", row.names = F)
write.csv(P1_euka03_DNA_agg_mean$reads,"Output_metabar_data/P1_euka03_DNA/files/P1_euka03_DNA_reads_mean_final.csv", row.names = T)

write.csv(P1_euka03_DNA_agg_prob$reads,"Output_metabar_data/P1_euka03_DNA/files/P1_euka03_DNA_reads_prob_final.csv", row.names = T)


P1_euka03_DNA_agg_prob <- aggregate_pcrs(P1_euka03_DNA_clean, FUN = FUN_agg_pcrs_prob)
summary_metabarlist(P1_euka03_DNA_agg_prob)

P1_euka03_DNA_agg_sum <- aggregate_pcrs(P1_euka03_DNA_clean, FUN = FUN_agg_pcrs_sum)
summary_metabarlist(P1_euka03_DNA_agg_sum)

write.csv(P1_euka03_DNA_agg_sum$reads,"Output_metabar_data/P1_euka03_DNA/files/P1_euka03_DNA_reads_sum_final.csv", row.names = T)

