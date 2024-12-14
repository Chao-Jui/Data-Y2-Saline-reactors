# library(devtools)


# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("Biostrings")
# BiocManager::install("S4Vectors")

#install_github("poyuliu/KTU")

library(Biostrings)
library(KTU)
library(tidyr) 
library(dplyr)
library(ggplot2)
library(ggpubr)

wd="Your/PATH/KTU/UASB_Daughter1_110" #  where the input files are saved
setwd(wd)


#import ASV table
asv <- read.csv("ASV_table_normalized_UASB_Daughter1_110.csv") # import QIIME feature table with taxonomy annotation
tasv=as.data.frame(t(asv))
colnames(asv)=tasv$V1
asv <- asv[-1, ]
#asv <- asv[,-ncol(asv)] # remove taxonomy annotation column
#import represent sequenge fasta
sequences <- readDNAStringSet("dna-sequences.fasta")
#import the metadata (Feature.ID to taxa) 
#the metadata TSV file from 138_taxonomy.qzv in qiime2 view
metadata=read.csv("metadata.csv")
metadata <- metadata[-1, ] # remove the first row which is usless information
#metadata <- metadata[!names(metadata) %in% "Confidence"]
metadata <- separate(metadata, Taxon, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"), sep = "; ")

# Filter the rows with the Phylum column 
metadata_p__Chloroflexi <- metadata %>% filter(Phylum == "p__Chloroflexi")
metadata_p__Firmicutes <- metadata %>% filter(Phylum == "p__Firmicutes")
metadata_p__Proteobacteria <- metadata %>% filter(Phylum == "p__Proteobacteria")
metadata_p__Planctomycetota <- metadata %>% filter(Phylum == "p__Planctomycetota")
metadata_p__Bacteroidota <- metadata %>% filter(Phylum == "p__Bacteroidota")
#remove the ID in the file abrove 
ids_to_remove <- c("p__Chloroflexi","p__Firmicutes","p__Proteobacteria","p__Planctomycetota","p__Bacteroidota")
metadata_p__resudual <- metadata %>% filter(!Phylum %in% ids_to_remove)

#Subset of the ASV table
asv_p__Chloroflexi=asv %>% filter(`#OTU ID` %in% metadata_p__Chloroflexi$Feature.ID)
asv_p__Firmicutes=asv %>% filter(`#OTU ID` %in% metadata_p__Firmicutes$Feature.ID)
asv_p__Proteobacteria=asv %>% filter(`#OTU ID` %in% metadata_p__Proteobacteria$Feature.ID)
asv_p__Planctomycetota=asv %>% filter(`#OTU ID` %in% metadata_p__Planctomycetota$Feature.ID)
asv_p__Bacteroidota=asv %>% filter(`#OTU ID` %in% metadata_p__Bacteroidota$Feature.ID)
asv_p__resudual=asv %>% filter(`#OTU ID` %in% metadata_p__resudual$Feature.ID)

# To convert all character columns in a dataframe to numeric, except for a specific column (#OTU ID)
# Identify character columns excluding `#OTUID`
character_cols <- colnames(asv_p__Chloroflexi)[sapply(asv_p__Chloroflexi, is.character)]
character_cols <- setdiff(character_cols, "#OTU ID")
# Convert character columns to numeric
asv_p__Chloroflexi[character_cols] <- lapply(asv_p__Chloroflexi[character_cols], function(x) as.numeric(as.character(x)))
#
character_cols <- colnames(asv_p__Firmicutes)[sapply(asv_p__Firmicutes, is.character)]
character_cols <- setdiff(character_cols, "#OTU ID")
asv_p__Firmicutes[character_cols] <- lapply(asv_p__Firmicutes[character_cols], function(x) as.numeric(as.character(x)))
#
character_cols <- colnames(asv_p__Proteobacteria)[sapply(asv_p__Proteobacteria, is.character)]
character_cols <- setdiff(character_cols, "#OTU ID")
asv_p__Proteobacteria[character_cols] <- lapply(asv_p__Proteobacteria[character_cols], function(x) as.numeric(as.character(x)))
#
character_cols <- colnames(asv_p__Planctomycetota)[sapply(asv_p__Planctomycetota, is.character)]
character_cols <- setdiff(character_cols, "#OTU ID")
asv_p__Planctomycetota[character_cols] <- lapply(asv_p__Planctomycetota[character_cols], function(x) as.numeric(as.character(x)))
#
character_cols <- colnames(asv_p__Bacteroidota)[sapply(asv_p__Bacteroidota, is.character)]
character_cols <- setdiff(character_cols, "#OTU ID")
asv_p__Bacteroidota[character_cols] <- lapply(asv_p__Bacteroidota[character_cols], function(x) as.numeric(as.character(x)))
#
character_cols <- colnames(asv_p__resudual)[sapply(asv_p__resudual, is.character)]
character_cols <- setdiff(character_cols, "#OTU ID")
asv_p__resudual[character_cols] <- lapply(asv_p__resudual[character_cols], function(x) as.numeric(as.character(x)))

#Check column type
sapply(asv_p__Chloroflexi, class)


# Select sequences by their IDs
sequences_p__Chloroflexi<- sequences[names(sequences) %in% metadata_p__Chloroflexi$Feature.ID]
sequences_p__Firmicutes<- sequences[names(sequences) %in% metadata_p__Firmicutes$Feature.ID]
sequences_p__Proteobacteria<- sequences[names(sequences) %in% metadata_p__Proteobacteria$Feature.ID]
sequences_p__Planctomycetota<- sequences[names(sequences) %in% metadata_p__Planctomycetota$Feature.ID]
sequences_p__Bacteroidota<- sequences[names(sequences) %in% metadata_p__Bacteroidota$Feature.ID]
sequences_p__resudual <- sequences[names(sequences) %in% metadata_p__resudual$Feature.ID]

#save the file as new fasta
writeXStringSet(sequences_p__Chloroflexi, "dna-sequences-Chloroflexi.fasta")
writeXStringSet(sequences_p__Firmicutes, "dna-sequences-Firmicutes.fasta")
writeXStringSet(sequences_p__Proteobacteria, "dna-sequences-Proteobacteria.fasta")
writeXStringSet(sequences_p__Planctomycetota, "dna-sequences-Planctomycetota.fasta")
writeXStringSet(sequences_p__Bacteroidota, "dna-sequences-Bacteroidota.fasta")
writeXStringSet(sequences_p__resudual, "dna-sequences-resudual.fasta")

#Run KTU

#Chloroflexi
kluster.Chloroflexi <- klustering(repseq = "dna-sequences-Chloroflexi.fasta",
                                  feature.table = asv_p__Chloroflexi,
                                  write.fasta = TRUE,cores = 10)
saveRDS(kluster.Chloroflexi,file="kluster-Chloroflexi.RDS") 
ktu_eval.Chloroflexi <- KTUsim.eval(klusterRDS = "kluster-Chloroflexi.RDS", ASVfasta = "dna-sequences-Chloroflexi.fasta")
write.csv(kluster.Chloroflexi$KTU.table,file="KTU.table.Chloroflexi.csv")
write.csv(kluster.Chloroflexi$ReqSeq ,file="rep-seq.Chloroflexi.csv")
write.csv(kluster.Chloroflexi$clusters ,file="cluster.Chloroflexi.csv")

#Firmicutes
kluster.Firmicutes <- klustering(repseq = "dna-sequences-Firmicutes.fasta",
                                 feature.table = asv_p__Firmicutes,
                                 write.fasta = TRUE,cores = 10)
saveRDS(kluster.Firmicutes,file="kluster-Firmicutes.RDS") 
ktu_eval.Firmicutes <- KTUsim.eval(klusterRDS = "kluster-Firmicutes.RDS", ASVfasta = "dna-sequences-Firmicutes.fasta")
write.csv(kluster.Firmicutes$KTU.table,file="KTU.table.Firmicutes.csv")
write.csv(kluster.Firmicutes$ReqSeq ,file="rep-seq.Firmicutes.csv")
write.csv(kluster.Firmicutes$clusters ,file="cluster.Firmicutes.csv")

#Proteobacteria
kluster.Proteobacteria <- klustering(repseq = "dna-sequences-Proteobacteria.fasta",
                                     feature.table = asv_p__Proteobacteria,
                                     write.fasta = TRUE,cores = 10)
saveRDS(kluster.Proteobacteria,file="kluster-Proteobacteria.RDS") 
ktu_eval.Proteobacteria <- KTUsim.eval(klusterRDS = "kluster-Proteobacteria.RDS", ASVfasta = "dna-sequences-Proteobacteria.fasta")
write.csv(kluster.Proteobacteria$KTU.table,file="KTU.table.Proteobacteria.csv")
write.csv(kluster.Proteobacteria$ReqSeq ,file="rep-seq.Proteobacteria.csv")
write.csv(kluster.Proteobacteria$clusters ,file="cluster.Proteobacteria.csv")

#Planctomycetota
kluster.Planctomycetota <- klustering(repseq = "dna-sequences-Planctomycetota.fasta",
                                      feature.table = asv_p__Planctomycetota,
                                      write.fasta = TRUE,cores = 10)
saveRDS(kluster.Planctomycetota,file="kluster-Planctomycetota.RDS") 
ktu_eval.Planctomycetota <- KTUsim.eval(klusterRDS = "kluster-Planctomycetota.RDS", ASVfasta = "dna-sequences-Planctomycetota.fasta")
write.csv(kluster.Planctomycetota$KTU.table,file="KTU.table.Planctomycetota.csv")
write.csv(kluster.Planctomycetota$ReqSeq ,file="rep-seq.Planctomycetota.csv")
write.csv(kluster.Planctomycetota$clusters ,file="cluster.Planctomycetota.csv")

#Bacteroidota
kluster.Bacteroidota <- klustering(repseq = "dna-sequences-Bacteroidota.fasta",
                                   feature.table = asv_p__Bacteroidota,
                                   write.fasta = TRUE,cores = 10)
saveRDS(kluster.Bacteroidota,file="kluster-Bacteroidota.RDS") 
ktu_eval.Bacteroidota <- KTUsim.eval(klusterRDS = "kluster-Bacteroidota.RDS", ASVfasta = "dna-sequences-Bacteroidota.fasta")
write.csv(kluster.Bacteroidota$KTU.table,file="KTU.table.Bacteroidota.csv")
write.csv(kluster.Bacteroidota$ReqSeq ,file="rep-seq.Bacteroidota.csv")
write.csv(kluster.Bacteroidota$clusters ,file="cluster.Bacteroidota.csv")

#resudual
kluster.resudual <- klustering(repseq = "dna-sequences-resudual.fasta",
                               feature.table = asv_p__resudual,
                               write.fasta = TRUE,cores = 10)
saveRDS(kluster.resudual,file="kluster-resudual.RDS") 
ktu_eval.resudual <- KTUsim.eval(klusterRDS = "kluster-resudual.RDS", ASVfasta = "dna-sequences-resudual.fasta")
write.csv(kluster.resudual$KTU.table,file="KTU.table.resudual.csv")
write.csv(kluster.resudual$ReqSeq ,file="rep-seq.resudual.csv")
write.csv(kluster.resudual$clusters ,file="cluster.resudual.csv")




#Merge the sub result of the KTU 
KTU.table.summary=rbind(kluster.Chloroflexi$KTU.table,kluster.Firmicutes$KTU.table,kluster.Proteobacteria$KTU.table,
                    kluster.Planctomycetota$KTU.table,kluster.Bacteroidota$KTU.table,kluster.resudual$KTU.table)

kmer.table.summary=cbind(kluster.Chloroflexi$kmer.table,kluster.Firmicutes$kmer.table,kluster.Proteobacteria$kmer.table,
                     kluster.Planctomycetota$kmer.table,kluster.Bacteroidota$kmer.table,kluster.resudual$kmer.table)

ktu_eval.summary=rbind(ktu_eval.Chloroflexi$eachmean,ktu_eval.Firmicutes$eachmean,ktu_eval.Proteobacteria$eachmean,
                         ktu_eval.Planctomycetota$eachmean,ktu_eval.Bacteroidota$eachmean,ktu_eval.resudual$eachmean)

write.csv(KTU.table.summary,file="KTU.table.summary.csv")
write.csv(kmer.table.summary,file="kmer.table.summary.csv")
write.csv(ktu_eval.summary,file="ktu_evaluate.summary.csv")


# Convert count KTU table to relative abundance table
KTU.table.summary_relative_abundance <- KTU.table.summary %>%
  mutate(across(everything(), ~ . / sum(.)))
write.csv(KTU.table.summary_relative_abundance,file="KTU.table.summary_relative_abundance.csv")

#histogram
#how many ASVs are clustered in to a KTU
p.n.feature.histogram=
  ggplot(ktu_eval.summary, aes(x = n.feature)) +
  geom_histogram(binwidth = 1, fill = "navyblue", color = "black") +
  labs(x = "# of ASVs / KTU", y = "Frequency") +
  #scale_x_continuous(limits = c(0, 15),breaks=seq(0,15,1))+
  theme(axis.text.y = element_text(colour = "black", size = 18), #, 1face = "bold"
        axis.text.x = element_text(colour = "black", size = 18), 
        legend.text = element_text(size = 20, colour ="black"), 
        #legend.position = c(0.6,0.88),
        legend.position = "bottom",
        legend.direction = "horizontal", ##Vertical  #horizontal
        legend.box = "vertical", ##Vertical  #horizontal
        axis.title.y = element_text( size = 22), 
        axis.title.x = element_text( size = 20, colour = "black"),
        #axis.title.x = element_blank(),#hide label of x axis
        legend.title = element_text(size = 20, colour = "black"), 
        axis.ticks=element_line(colour="black",size=1,linetype=1), #Modify the size and color of the line label on the axis
        axis.ticks.length=unit(0.4,"lines"),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 2)#框線
  )

#mean ASV similarity within a KTU
p.similarity.pct.histogram=
  ggplot(ktu_eval.summary, aes(x = similarity.pct)) +
  geom_histogram(binwidth = 0.5, fill = "navyblue", color = "black") +
  labs(x = "% of ASVs within KTU", y = "Frequency") +
  theme(axis.text.y = element_text(colour = "black", size = 18), #, 1face = "bold"
        axis.text.x = element_text(colour = "black", size = 18), 
        legend.text = element_text(size = 20, colour ="black"), 
        #legend.position = c(0.6,0.88),
        legend.position = "bottom",
        legend.direction = "horizontal", ##Vertical  #horizontal
        legend.box = "vertical", ##Vertical  #horizontal
        axis.title.y = element_text( size = 22), 
        axis.title.x = element_text( size = 20, colour = "black"),
        #axis.title.x = element_blank(),#hide label of x axis
        legend.title = element_text(size = 20, colour = "black"), 
        axis.ticks=element_line(colour="black",size=1,linetype=1), #Modify the size and color of the line label on the axis
        axis.ticks.length=unit(0.4,"lines"),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 2)#框線
  )

#Divergence of ASVs within KTU 
p.divergence.histogram=
  ggplot(ktu_eval.summary, aes(x = divergence)) +
  geom_histogram(binwidth = 0.002, fill = "navyblue", color = "black") +
  labs(x = "Divergence of ASVs within KTU", y = "Frequency") +
  theme(axis.text.y = element_text(colour = "black", size = 18), #, 1face = "bold"
        axis.text.x = element_text(colour = "black", size = 18), 
        legend.text = element_text(size = 20, colour ="black"), 
        #legend.position = c(0.6,0.88),
        legend.position = "bottom",
        legend.direction = "horizontal", ##Vertical  #horizontal
        legend.box = "vertical", ##Vertical  #horizontal
        axis.title.y = element_text(size = 22), 
        axis.title.x = element_text(size = 20, colour = "black"),
        #axis.title.x = element_blank(),#hide label of x axis
        legend.title = element_text(size = 20, colour = "black"), 
        axis.ticks=element_line(colour="black",size=1,linetype=1), #Modify the size and color of the line label on the axis
        axis.ticks.length=unit(0.4,"lines"),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 2)#框線
  )


P.merge.KTU.histogram=
  ggarrange(p.n.feature.histogram,p.similarity.pct.histogram,p.divergence.histogram,
            align = "hv",
            labels = c("A", "B", "C", "D"),
            ncol=1, nrow=3,
            #common.legend = TRUE, legend="bottom",
            font.label = list(size = 30))

ggsave(plot = P.merge.KTU.histogram,"01_P.merge.KTU.histogram.png", width = 5.5, height = 11)



# #how many ASVs are clustered in to a KTU
# histNdistri(histdata = ktu_eval$eachmean$n.feature,las=1,breaks=30,
#             coll = colset.d.5[1],
#             xlab.text = "# of ASVs / KTU",
#             legend.posit = "topright",legend.text = paste0(round(mean(ktu_eval$eachmean$n.feature),2)," ASVs/KTU on average"))
# #mean ASV similarity within a KTU
# histNdistri(histdata = ktu_eval$eachmean$similarity.pct,las=1,breaks=30,
#             coll = colset.d.5[2],
#             xlab.text = "% of ASVs within KTU",
#             legend.posit = "topleft",legend.text = paste0("Mean similarity: ",round(ktu_eval$globalmean,2),"%"))
# #Divergence of ASVs within KTU  
# histNdistri(histdata = ktu_eval$eachmean$divergence[!is.na(ktu_eval$eachmean$divergence)],las=1,breaks=30,
#             coll = colset.d.5[3],
#             xlab.text = "Divergence of ASVs within KTU",
#             legend.posit = "topright",legend.text = paste0("Mean divergence: ",round(ktu_eval$globaldivergence,2)))

library(Biostrings)
#make a KTU-based dna-sequences.fasta 
fasta <- readDNAStringSet("dna-sequences.fasta")
cluster1 <- read.csv("KTU.table.Bacteroidota.csv")
cluster2 <- read.csv("KTU.table.Chloroflexi.csv")
cluster3 <- read.csv("KTU.table.Firmicutes.csv")
cluster4 <- read.csv("KTU.table.Planctomycetota.csv")
cluster5 <- read.csv("KTU.table.Proteobacteria.csv")
cluster6 <- read.csv("KTU.table.resudual.csv")

name_list <- c(cluster1$X,cluster2$X,cluster3$X,cluster4$X,cluster5$X,cluster6$X)
length(name_list)

fasta_KTU <- fasta[names(fasta) %in% name_list]

writeXStringSet(fasta_KTU, "dna-sequences-KTU.fasta")

