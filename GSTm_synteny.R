library(gggenomes)
library(tidyverse)
library(gggenes)
library(patchwork)

#colors
c25 <- c(
  "dodgerblue2", "antiquewhite", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "black", "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "#E31A1C",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown", "salmon2", "palegreen4", "grey33",
  "lightblue", "antiquewhite3","cornflowerblue")


#file with coordinates
f <- "/mnt/md0/mouse_detox/synteny/gene_flank_locs_final.tsv"

#load in data
ncbi <- read.table(f, header = T, sep = "\t")


#flip taxa so all have same orientation
p <- ncbi[which(ncbi$Gene == "Ampd2" & ncbi$Orientation == "+"),]$Species_name
f <- ncbi[which(ncbi$Species_name %in% p),]
fg <- ncbi[which(!ncbi$Species_name %in% p),]
mx <- f %>% group_by(Species_name) %>% summarise(max = max(End))

for(x in 1:nrow(mx)){
  f[which(f$Species_name %in% mx[x,1]),3:2] <- abs(f[which(f$Species_name %in% mx[x,1]),2:3] - mx[x,2]$max)
}

f$Orientation <- case_when(f$Orientation == "+" ~ as.factor("-"), 
                                                    f$Orientation == "-" ~ as.factor("+"))

#combine newly oriented taxa
oriented <- rbind(f, fg)
cg <- which(oriented$Species_name == "Cricetulus_griseus")
or_sub <-oriented[-cg,]
or_sub <- subset(or_sub, Species_name != "Peromyscus_maniculatus")
or_sub <- subset(or_sub, Species_name != "Onychomys_torridus")
or_sub <- subset(or_sub, Gene != "Csf1")
or_sub <- subset(or_sub, Gene != "Gnat2")
or_sub <- subset(or_sub, Gene != "Gm12500")
or_sub <- subset(or_sub, Gene != "Unch")
or_sub <- subset(or_sub, Gene != "Rpl13")
or_sub$seq_id <- paste(or_sub$Species_name, or_sub$Chromosome, sep = "_")
or_sub <- subset(or_sub, seq_id != "Peromyscus_leucopus_Chr22")
or_sub <- subset(or_sub, seq_id != "Phyllotis_vacarrum_Scaf1")
or_sub[which(or_sub$Species_name == "Mastomys_coucha" & or_sub$Gene == "Eps8l3"), 2] <- or_sub[which(or_sub$Species_name == "Mastomys_coucha" & or_sub$Gene == "Eps8l3"), 2] + 185000
or_sub$seq_id <- gsub("_", " ", or_sub$seq_id)

#get ready for plotting
plt <- or_sub[,c(4,2,3,8,11,7)]
colnames(plt) <- c("strand", "start", "end", "Common_name", "seq_id", "Gene")
fctrs <- unique(plt$seq_id)[c(3,5,4,2,6,7,1)]
plt$seq_id <- factor(plt$seq_id, levels = fctrs)

#make linking dataframe
lnk <- merge(plt, plt, by = "Gene")
names(lnk) <- gsub(".x","", names(lnk))
names(lnk) <- gsub(".y","2", names(lnk))

#

#make plot
fp <- gggenomes(genes = plt, links = lnk) + 
  geom_seq() +
  geom_bin_label(size = 2.55) +
  geom_gene()+
  geom_link_line(linetype = "dotted", color = "grey", alpha = 0.7, linewidth = 0.5)+
  #geom_feat() + 
  geom_gene(aes(fill=Gene)) +
  geom_gene_tag(aes(label=Gene), nudge_y=0.1, check_overlap = T) + 
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))+
  ggtitle("GST\u03bc Clusters") + scale_fill_manual(values = c25)
fp

cairo_pdf("/mnt/md0/mouse_detox/synteny/synteny_plot.pdf", width = 8.45, height = 6.5)
fp
dev.off()
