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
euka01_P3_DNA_raw <- read.table("Obitools_output/14.tab_P3DNA_euka01.txt", header =T, sep="\t")
euka01_P3_DNA_raw <- euka01_P3_DNA_raw %>% select("id", "cluster_weight", "count", starts_with("sample"), 
                                                    "obitag_bestid", "obitag_bestmatch",  "obitag_match_count", "obitag_rank",
                                                    "obitag_similarity_method", "scientific_name", "taxid", "taxonomic_path",
                                                    "sequence")

### read table ####
reads_P3_euka01_DNA <-t(euka01_P3_DNA_raw[,grep("sample\\.", colnames(euka01_P3_DNA_raw))])
rownames(reads_P3_euka01_DNA) = gsub("sample\\.", "", rownames(reads_P3_euka01_DNA))
colnames(reads_P3_euka01_DNA) = euka01_P3_DNA_raw$id
rownames(reads_P3_euka01_DNA) = gsub("\\.", "_", rownames(reads_P3_euka01_DNA))
reads_P3_euka01_DNA <- as.matrix(reads_P3_euka01_DNA)
dim(reads_P3_euka01_DNA)
#view(reads_P3_euka01_DNA)

### motus table ####
motus_P3_euka01_DNA = euka01_P3_DNA_raw[,grep("sample\\.", colnames(euka01_P3_DNA_raw), invert = T)]
rownames(motus_P3_euka01_DNA)=motus_P3_euka01_DNA$id
#View(motus_P3_euka01_DNA)

## Samples dataset
sample_P3<-read.csv("sample_files/Samplefile_P3.csv", header= T, sep= ";")
#View(sample_P3)
#We relabel samples file rownames to be the Sample ID so we can merge by Sample_ID with the reads table.
sample_P3 = sample_P3[order(sample_P3$Sample_ID),]
sample_P3 = sample_P3[!duplicated(sample_P3$Sample_ID),]
sample_P3$Sample_ID <- gsub("-", "_", sample_P3$Sample_ID)
rownames(sample_P3)=sample_P3$Sample_ID
#str(samples)
samp_P3 = sample_P3
#[match(rownames(reads_P3_euka01_DNA), rownames(sample_P3)),]
#View(samp_P3)

## NGS filterfile ####
ngs_P3 <- read.table("ngs_filters/P3DNA_euka01_ngsfilter.txt", sep="\t")
ngs_P3$V2 <- gsub("-", "_", ngs_P3$V2)

## PCRs GPS ####
pcr_gps_P3 <- read.csv("PCR_GPS/PCR_GPS_P3.csv", header = T, sep= ";") 
pcr_gps_P3$Sample_ID <- gsub("-", "_", pcr_gps_P3$Sample_ID)

pcrs_P3 <- merge(samp_P3, pcr_gps_P3, by.y="Sample_ID")

### merge PCR and GPS ####
pcrs_P3 <- merge(pcrs_P3, ngs_P3[,c(2:5)], by.x<-"Sample_ID", by.y = "V2", all=T)

rownames(pcrs_P3) <- pcrs_P3$Sample_ID
names(pcrs_P3)[names(pcrs_P3) == "Sample_ID"] <- "sample_id"
names(pcrs_P3)[names(pcrs_P3) == "V3"] <- "tags"
names(pcrs_P3)[names(pcrs_P3) == "V4"] <- "primer_fwd"
names(pcrs_P3)[names(pcrs_P3) == "V5"] <- "primer_rev"

pcrs_P3 <- pcrs_P3 %>% separate(tags, c("tag_fwd", "tag_rev"), ":")
pcrs_P3 <- as.data.frame(pcrs_P3)
pcrs_P3 <- subset(pcrs_P3, !is.na(pcrs_P3$type))
row.names(pcrs_P3)<- pcrs_P3$sample_id
#view(pcrs_P3)

pcrs_P3 <- pcrs_P3[,-1]
names(pcrs_P3)[names(pcrs_P3) == "Sample"] <- "sample_id"

## Metadata ####
metadata<-read.csv("Samples_info_Planaqua.csv", header= T, sep=";")

#For the sample info, we need to select the samples corresponding at the library. 
metadata_P3 <- merge(metadata, samp_P3, by.x = "Samples_ID", by.y = "Sample_pla")

#DOUBLE CHECK HERE
metadata_P3 = metadata_P3[!duplicated(metadata_P3[,c('Sample')]),]
rownames(metadata_P3)=metadata_P3$Sample
metadata_P3 = as.data.frame(metadata_P3)
#View(metadata_P3)

# METABAR LIST ####
P3_euka01_DNA <- metabarlist_generator(reads = reads_P3_euka01_DNA, motus = motus_P3_euka01_DNA, pcrs = pcrs_P3, samples = metadata_P3)
summary_metabarlist(P3_euka01_DNA)

# we can remove all the useless dataset from the env:
rm(euka01_P3_DNA_raw, metadata, metadata_P3, motus_P3_euka01_DNA, ngs_P3, pcr_gps_P3, pcrs_P3, samp_P3, sample_P3, reads_P3_euka01_DNA)

# OVERVIEW OF CONTAMINANTS ####
# Compute the number of reads per pcr                                   
P3_euka01_DNA$pcrs$nb_reads <- rowSums(P3_euka01_DNA$reads)
# Compute the number of motus per pcr
P3_euka01_DNA$pcrs$nb_motus <- rowSums(P3_euka01_DNA$reads>0)

# Create an input table (named check1) for ggplot of 3 columns: 
#  (i) control type 
#  (ii) a vector indicated whether it corresponds to nb_reads or nb_motus, 
#  (iii) the corresponding values.

check1 <- reshape2::melt(P3_euka01_DNA$pcrs[,c("control_type", "nb_reads", "nb_motus")])

P1<-ggplot(data <- check1, aes(x=control_type, y=value, color=control_type)) + 
  geom_boxplot() + theme_bw() + 
  geom_jitter(alpha=0.2) + 
  scale_color_manual(values = c("brown", "red", "cyan4","pink"), na.value = "darkgrey") +
  facet_wrap(~variable, scales = "free_y") + 
  theme(axis.text.x = element_text(angle=45, h=1))
P1
#ggsave(filename="1_P3_euka01_DNA_barplot.contaminants.pdf", plot=P1, width = 297, height = 210, units = "mm", device="pdf", path="Output_metabar_data/P3_euka01_DNA/figures_cleaning")

P2<-ggplot(P3_euka01_DNA$pcrs, aes(x=nb_reads, y=nb_motus, color = control_type)) + 
  geom_point() + theme_bw() + 
  scale_y_log10() + scale_x_log10() + 
  scale_color_manual(values = c("brown", "red", "cyan4","pink"), na.value = "darkgrey")
P2
#ggsave(filename="2_P3_euka01_DNA_plot.motus.vs.reads.pdf", plot=P2, width = 297, height = 210, units = "mm", device="pdf", path="Output_metabar_data/P3_euka01_DNA/figures_cleaning")

p3<-ggpcrplate(P3_euka01_DNA, FUN = function(m){rowSums(m$reads)}, legend_title = "# of reads per PCR")
p3
#ggsave(filename="3_P3_euka01_DNA_pcr.plots.pdf", plot=p3, width = 297, height = 210, units = "mm", device="pdf", path="Output_metabar_data/P3_euka01_DNA/figures_cleaning")

tag.list <- as.character(unique(P3_euka01_DNA$pcrs$tag_rev))
p4<-ggpcrtag(P3_euka01_DNA, 
             legend_title = "# of reads per PCR", 
             FUN = function(m) {rowSums(m$reads)},
             taglist = tag.list) 
p4
#ggsave(filename="4_P3_euka01_DNA_pcr.plots.by.tags.pdf", plot=p4, width = 297, height = 210, units = "mm", device="pdf", path="Output_metabar_data/P3_euka01_DNA/figures_cleaning")

# FLAGGING SPURIOUS SIGNAL ####
## Extraction control ####
P3_euka01_DNA <- contaslayer(P3_euka01_DNA, 
                              control_types = c("extraction"),
                              output_col = "not_an_extraction_conta")

table(P3_euka01_DNA$motus$not_an_extraction_conta)

id <- !P3_euka01_DNA$motus$not_an_extraction_conta
max.conta <- rownames(P3_euka01_DNA$motus[id,])[which.max(P3_euka01_DNA$motus[id, "count"])]

#... and its distribution and relative abundance in each pcr
# the "#reads of most abundant contaminat" is probably refferring to the % of reads that are labelled contaminant (check function above)
ggpcrplate(P3_euka01_DNA, legend_title = "#reads of most \nabundant contaminant",
           FUN = function(m) {m$reads[, max.conta]/rowSums(m$reads)})

# Compute relative abundance of all pcr contaminants together
a <- data.frame(conta.relab = rowSums(P3_euka01_DNA$reads[,!P3_euka01_DNA$motus$not_an_extraction_conta]) /
                  rowSums(P3_euka01_DNA$reads))
# Add information on control types
a$control_type <- P3_euka01_DNA$pcrs$control_type[match(rownames(a), rownames(P3_euka01_DNA$pcrs))]

p6<-ggplot(a, aes(x=control_type, y=conta.relab, color=control_type)) +
  geom_boxplot() + geom_jitter(alpha=0.5) +
  scale_color_manual(values = c("brown", "red", "cyan4","pink"), na.value = "darkgrey") +
  labs(x=NULL, y="Prop. Reads (log10)") +
  theme_bw() +
  scale_y_log10()
p6
#ggsave(filename="6_P3_euka01_DNA_boxplot_contamination_extraction.pdf", plot=p6, width = 297, height = 210, units = "mm", device="pdf", path="Output_metabar_data/P3_euka01_DNA/figures_cleaning")

# flag pcrs with total contaminant relative abundance > 10% of reads)
P3_euka01_DNA$pcrs$low_contamination_level <- 
  ifelse(a$conta.relab[match(rownames(P3_euka01_DNA$pcrs), rownames(a))]>1e-1,  F, T)
# Proportion of potentially functional (TRUE) vs. failed (FALSE) pcrs
# (controls included) based on this criterion
table(P3_euka01_DNA$pcrs$low_contamination_level) / nrow(P3_euka01_DNA$pcrs)

## PCR control ####
P3_euka01_DNA <- contaslayer(P3_euka01_DNA, 
                              control_types = c("pcr"),
                              output_col = "not_an_pcr_conta")
table(P3_euka01_DNA$motus$not_an_pcr_conta)

# Identify the most common contaminant
# get contaminant ids
id <- !P3_euka01_DNA$motus$not_an_pcr_conta
max.conta <- rownames(P3_euka01_DNA$motus[id,])[which.max(P3_euka01_DNA$motus[id, "count"])]

#... and its distribution and relative abundance in each pcr
# the "#reads of most abundant contaminat" is probably refferring to the % of reads that are labelled contaminant (check function above)
ggpcrplate(P3_euka01_DNA, legend_title = "#reads of most \nabundant contaminant",
           FUN = function(m) {m$reads[, max.conta]/rowSums(m$reads)})

# Compute relative abundance of all pcr contaminants together
a <- data.frame(conta.relab = rowSums(P3_euka01_DNA$reads[,!P3_euka01_DNA$motus$not_an_pcr_conta]) /
                  rowSums(P3_euka01_DNA$reads))
# Add information on control types
a$control_type <- P3_euka01_DNA$pcrs$control_type[match(rownames(a), rownames(P3_euka01_DNA$pcrs))]

p7<-ggplot(a, aes(x=control_type, y=conta.relab, color=control_type)) +
  geom_boxplot() + geom_jitter(alpha=0.5) +
  scale_color_manual(values = c("brown", "red", "cyan4","pink"), na.value = "darkgrey") +
  labs(x=NULL, y="Prop. Reads (log10)") +
  theme_bw() +
  scale_y_log10()

p7
#ggsave(filename="7_P3_euka01_DNA_boxplot_contamination_pcr.pdf", plot=p7, width = 297, height = 210, units = "mm", device="pdf", path="Output_metabar_data/P3_euka01_DNA/figures_cleaning")

# flag pcrs with total contaminant relative abundance > 10% of reads)
P3_euka01_DNA$pcrs$low_contamination_level <- 
  ifelse(a$conta.relab[match(rownames(P3_euka01_DNA$pcrs), rownames(a))]>1e-1,  F, T)
# Proportion of potentially functional (TRUE) vs. failed (FALSE) pcrs
# (controls included) based on this criterion
table(P3_euka01_DNA$pcrs$low_contamination_level) / nrow(P3_euka01_DNA$pcrs)

## Sequencing control ####
P3_euka01_DNA <- contaslayer(P3_euka01_DNA, 
                              control_types = c("sequencing"),
                              output_col = "not_an_seq_conta")
table(P3_euka01_DNA$motus$not_an_seq_conta)

# Identify the most common contaminant
# get contaminant ids
id <- !P3_euka01_DNA$motus$not_an_seq_conta
max.conta <- rownames(P3_euka01_DNA$motus[id,])[which.max(P3_euka01_DNA$motus[id, "count"])]

#... and its distribution and relative abundance in each pcr
# the "#reads of most abundant contaminat" is probably refferring to the % of reads that are labelled contaminant (check function above)
ggpcrplate(P3_euka01_DNA, legend_title = "#reads of most \nabundant contaminant",
           FUN = function(m) {m$reads[, max.conta]/rowSums(m$reads)})

# Compute relative abundance of all pcr contaminants together
a <- data.frame(conta.relab = rowSums(P3_euka01_DNA$reads[,!P3_euka01_DNA$motus$not_an_seq_conta]) /
                  rowSums(P3_euka01_DNA$reads))
# Add information on control types
a$control_type <- P3_euka01_DNA$pcrs$control_type[match(rownames(a), rownames(P3_euka01_DNA$pcrs))]

p8<-ggplot(a, aes(x=control_type, y=conta.relab, color=control_type)) +
  geom_boxplot() + geom_jitter(alpha=0.5) +
  scale_color_manual(values = c("brown", "red", "cyan4","pink"), na.value = "darkgrey") +
  labs(x=NULL, y="Prop. Reads (log10)") +
  theme_bw() +
  scale_y_log10()
p8
#ggsave(filename="8_P3_euka01_DNA_boxplot_contamination_sequencing.pdf", plot=p8, width = 297, height = 210, units = "mm", device="pdf", path="Output_metabar_data/P3_euka01_DNA/figures_cleaning")

# flag pcrs with total contaminant relative abundance > 10% of reads)
P3_euka01_DNA$pcrs$low_contamination_level <- 
  ifelse(a$conta.relab[match(rownames(P3_euka01_DNA$pcrs), rownames(a))]>1e-1,  F, T)
# Proportion of potentially functional (TRUE) vs. failed (FALSE) pcrs
# (controls included) based on this criterion
table(P3_euka01_DNA$pcrs$low_contamination_level) / nrow(P3_euka01_DNA$pcrs)

## Non-target MOTUs ####
#Flag MOTUs corresponding to target (TRUE) vs. non-target (FALSE) taxa
P3_euka01_DNA$motus$target_taxon <- grepl("Eukaryota", P3_euka01_DNA$motus$taxonomic_path)

# Proportion of each of these over total number of MOTUs
table(P3_euka01_DNA$motus$target_taxon) / nrow(P3_euka01_DNA$motus)

# Plot the unweighted distribution of MOTUs similarity scores 
a <- 
  ggplot(P3_euka01_DNA$motus, aes(x=obitag_bestid)) + 
  geom_histogram(color="grey", fill="white", bins=20) + 
  geom_vline(xintercept = 0.80, col="orange", lty=2) + 
  theme_bw() + 
  theme(panel.grid = element_blank()) + 
  labs(x="% similarity against best match", y="# MOTUs")

# Same for the weighted distribution
b <- 
  ggplot(P3_euka01_DNA$motus, 
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
#ggsave(filename="9_P3_euka01_DNA_barplot_db_best_match.pdf", plot=p9, width = 297, height = 210, units = "mm", device="pdf", path="Output_metabar_data/P3_euka01_DNA/figures_cleaning")

# Flag not degraded (TRUE) vs. potentially degraded sequences (FALSE)
P3_euka01_DNA$motus$not_degraded <-
  ifelse(P3_euka01_DNA$motus$obitag_bestid < 0.80, F, T)

# Proportion of each of these over total number of MOTUs
table(P3_euka01_DNA$motus$not_degraded) / nrow(P3_euka01_DNA$motus)

## Detecting PCR outliers ####
#sequencing depth

ggplot(P3_euka01_DNA$pcrs, aes(nb_reads)) +
  geom_histogram(bins=40, color="grey", fill="white") + 
  geom_vline(xintercept = 1e2, lty=2, color="orange") + # threshold
  scale_x_log10() + 
  labs(x="# Reads (with all MOTUs and PCRs)", 
       y="# PCRs") +
  theme_bw() + 
  theme(panel.grid = element_blank())

# Flag pcrs with an acceptable sequencing depth (TRUE) or inacceptable one (FALSE)
P3_euka01_DNA$pcrs$seqdepth_ok <- ifelse(P3_euka01_DNA$pcrs$nb_reads < 100, F, T)

# Proportion of each of these over total number of pcrs, control excluded
table(P3_euka01_DNA$pcrs$seqdepth_ok[P3_euka01_DNA$pcrs$type=="sample"]) /
  nrow(P3_euka01_DNA$pcrs[P3_euka01_DNA$pcrs$type=="sample",])

# second way == reproducibility between replicates
# Subsetting the metabarlist
P3_euka01_DNA_sub <- subset_metabarlist(P3_euka01_DNA,
                                         table = "pcrs",
                                         indices = P3_euka01_DNA$pcrs$nb_reads>0 & (
                                           is.na(P3_euka01_DNA$pcrs$control_type) |
                                             P3_euka01_DNA$pcrs$control_type=="positive"))

# First visualization
comP3 = pcr_within_between(P3_euka01_DNA_sub)
P10<- check_pcr_thresh(comP3)
P10
#ggsave(filename="10_P3_euka01_DNA_density_plot_withinbetween_pcr_rep.pdf", plot=P10, width = 297, height = 210, units = "mm", device="pdf", path="Output_metabar_data/P3_euka01_DNA/figures_cleaning")

# Flagging
P3_euka01_DNA_sub <- pcrslayer(P3_euka01_DNA_sub, output_col = "replicating_pcr", plot = F)

# Proportion of replicating pcrs (TRUE)
table(P3_euka01_DNA_sub$pcrs$replicating_pcr) /
  nrow(P3_euka01_DNA_sub$pcrs)

# Intersection with the sequencing depth criterion
table(P3_euka01_DNA_sub$pcrs$seqdepth_ok,
      P3_euka01_DNA_sub$pcrs$replicating_pcr)

# Distinguish between pcrs obtained from samples from positive controls
mds = check_pcr_repl(P3_euka01_DNA_sub,
                     groups = P3_euka01_DNA_sub$pcrs$type,
                     funcpcr = P3_euka01_DNA_sub$pcrs$replicating_pcr)
mds + labs(color = "pcr type") + scale_color_manual(values = c("cyan4", "gray"))

P3_euka01_DNA$pcrs$replicating_pcr <- NA
P3_euka01_DNA$pcrs[rownames(P3_euka01_DNA_sub$pcrs),"replicating_pcr"] <- P3_euka01_DNA_sub$pcrs$replicating_pcr

## lowering tag jump ####
# Define a vector of thresholds to test
thresholds <- c(0,1e-4,1e-3, 1e-2, 3e-2, 5e-2)

# Run the tests and stores the results in a list
tests <- lapply(thresholds, function(x) tagjumpslayer(P3_euka01_DNA,x))
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
tmp$controls <- P3_euka01_DNA$pcrs$control_type[match(tmp$sample, rownames(P3_euka01_DNA$pcrs))]
tmp$threshold <- as.numeric(gsub("t_", "", tmp$threshold))

# New table formatting for ggplot
tmP3 <- melt(tmp, id.vars=colnames(tmp)[-grep("abundance|richness", colnames(tmp))])

P11<- ggplot(tmP3, aes(x=as.factor(threshold), y=value)) +
  geom_boxplot(color="grey40") +
  geom_vline(xintercept = which(levels(as.factor(tmP3$threshold)) == "0.01"), col="orange", lty=2) +
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
P11

ggsave(filename="11_P3_euka01_DNA_tagjump.pdf", plot=P11, width = 297, height = 210, units = "mm", device="pdf", path="Output_metabar_data/P3_euka01_DNA/figures_cleaning")

# Create a table of MOTUs quality criteria
# noise is identified as FALSE in soil_euk, the "!" transforms it to TRUE
motus.qual <- !P3_euka01_DNA$motus[,c("not_an_extraction_conta", "not_an_pcr_conta", "not_an_seq_conta", "target_taxon", "not_degraded")]

colnames(motus.qual) <- c("extraction_conta", "neg_conta", "seq_conta", "untargeted_taxon", "degraded_seq")

# Proportion of MOTUs potentially artifactual (TRUE) based on the criteria used
prop.table(table(apply(motus.qual, 1, sum) > 0))

# Corresponding proportion of artifactual reads (TRUE)
prop.table(xtabs(P3_euka01_DNA$motus$count~apply(motus.qual, 1, sum) > 0))

# Proportion of MOTUs and reads potentially artifactual for each criterion
apply(motus.qual, 2, sum) / nrow(motus.qual)
apply(motus.qual, 2, function(x) sum(P3_euka01_DNA$motus$count[x])/sum(P3_euka01_DNA$motus$count))

tmp.motus <-
  apply(sapply(1:ncol(motus.qual), function(x) {
    ifelse(motus.qual[,x]==T, colnames(motus.qual)[x], NA)}), 1, function(x) {
      paste(sort(unique(x)), collapse = "|")
    })
tmp.motus <- as.data.frame(gsub("^$", "not_artefactual", tmp.motus))
colnames(tmp.motus) <-  "artefact_type"

P12<-ggplot(tmp.motus, aes(x=1, fill=artefact_type)) +
  geom_bar() +  xlim(0, 2) +
  labs(fill="Artifact type") +
  coord_polar(theta="y") + theme_void() +
  scale_fill_brewer(palette = "Set3") +
  theme(legend.direction = "vertical") +
  ggtitle("MOTUs artefacts overview")

P12
ggsave(filename="12_P3_euka01_DNA_MOTUs_Artefacts_overview.pdf", plot=P12, width = 297, height = 210, units = "mm", device="pdf", path="Output_metabar_data/P3_euka01_DNA/figures_cleaning")

# Create a table of pcrs quality criteria
# noise is identified as FALSE in soil_euk, the "!" transforms it to TRUE
pcrs.qual <- !P3_euka01_DNA$pcrs[,c("low_contamination_level", "seqdepth_ok")]
colnames(pcrs.qual) <- c("high_contamination_level", "low_seqdepth")

# Proportion of pcrs potentially artifactual (TRUE) based on the criteria used
# excluding controls
prop.table(table(apply(pcrs.qual[P3_euka01_DNA$pcrs$type=="sample",], 1, sum) > 0))

# Proportion of MOTUs and reads potentially artifactual for each criterion
apply(pcrs.qual[P3_euka01_DNA$pcrs$type=="sample",], 2, sum) / nrow(pcrs.qual[P3_euka01_DNA$pcrs$type=="sample",])

tmp.pcrs <-
  apply(sapply(1:ncol(pcrs.qual), function(x) {
    ifelse(pcrs.qual[P3_euka01_DNA$pcrs$type=="sample",x]==T,
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
P3_euka01_DNA_clean <- subset_metabarlist(tmp, "pcrs",
                                           indices = rowSums(tmp$pcrs[,c("low_contamination_level",
                                                                         "seqdepth_ok", "replicating_pcr")]) == 3 &
                                             tmp$pcrs$type == "sample")
summary_metabarlist(P3_euka01_DNA_clean)

#Checking if there is still some empty PCR
if(sum(colSums(P3_euka01_DNA_clean$reads)==0)>0){print("empty motus present")}
if(sum(rowSums(P3_euka01_DNA_clean$reads)==0)>0){print("empty pcrs present")}

#Since we have now removed certain MOTUs or reduced their read counts, we need to update some parameters in the metabarlist (e.g. counts, etc.).
P3_euka01_DNA_clean$motus$count = colSums(P3_euka01_DNA_clean$reads)
P3_euka01_DNA_clean$pcrs$nb_reads_postmetabaR = rowSums(P3_euka01_DNA_clean$reads)
P3_euka01_DNA_clean$pcrs$nb_motus_postmetabaR = rowSums(ifelse(P3_euka01_DNA_clean$reads>0, T, F))

#We can now compare some basic characteristics before and after data curation with metabaR.
check <- melt(P3_euka01_DNA_clean$pcrs[,c("nb_reads", "nb_reads_postmetabaR",
                                           "nb_motus", "nb_motus_postmetabaR")])
check$type <- ifelse(grepl("motus", check$variable), "richness", "abundance")

P13<- ggplot(data = check, aes(x = variable, y = value)) +
  geom_boxplot( color = "darkgrey") +
  geom_jitter(alpha=0.1, color = "darkgrey") +
  theme_bw() +
  facet_wrap(~type, scales = "free", ncol = 5) +
  theme(axis.text.x = element_text(angle=45, h=1))
P13

ggsave(filename="13_P3_euka01_DNA_effects_cleaning_on_samples.pdf", plot=P13, width = 297, height = 210, units = "mm", device="pdf", path="Output_metabar_data/P3_euka01_DNA/figures_cleaning")


# DATA AGGREGATION ####
P3_euka01_DNA_agg_mean <- aggregate_pcrs(P3_euka01_DNA_clean, FUN = FUN_agg_pcrs_mean)
summary_metabarlist(P3_euka01_DNA_agg_mean)

P3_euka01_DNA_agg_prob <- aggregate_pcrs(P3_euka01_DNA_clean, FUN = FUN_agg_pcrs_prob)
summary_metabarlist(P3_euka01_DNA_agg_prob)

write.csv(P3_euka01_DNA_agg_mean$motus,"Output_metabar_data/P3_euka01_DNA/files/P3_euka01_DNA_motus_final.csv", row.names = T)
write.csv(P3_euka01_DNA_agg_mean$pcrs,"Output_metabar_data/P3_euka01_DNA/files/P3_euka01_DNA_pcrs_final.csv", row.names = F)
write.csv(P3_euka01_DNA_agg_mean$samples,"Output_metabar_data/P3_euka01_DNA/files/P3_euka01_DNA_samples_final.csv", row.names = F)
write.csv(P3_euka01_DNA_agg_mean$reads,"Output_metabar_data/P3_euka01_DNA/files/P3_euka01_DNA_reads_mean_final.csv", row.names = T)

P3_euka01_DNA_agg_sum <- aggregate_pcrs(P3_euka01_DNA_clean, FUN = FUN_agg_pcrs_sum)
summary_metabarlist(P3_euka01_DNA_agg_sum)

write.csv(P3_euka01_DNA_agg_sum$reads,"Output_metabar_data/P3_euka01_DNA/files/P3_euka01_DNA_reads_sum_final.csv", row.names = T)


write.csv(P3_euka01_DNA_agg_prob$reads,"Output_metabar_data/P3_euka01_DNA/files/P3_euka01_DNA_reads_prob_final.csv", row.names = T)
