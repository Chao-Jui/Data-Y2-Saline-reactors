rm(list=ls())

set.seed(123321)

library(ggpubr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(scales)
library(zoo)



# seting the path
wd="D:/R_project/UASB_Daughter_Y2_1_110_scientific_data" #  where the input files are saved
setwd(wd)

# Figure2
UASB_Saline_data <- read.csv("UASB_Saline_data.csv", fileEncoding = "ISO-8859-1",check.names = FALSE)

colors1 <- brewer.pal(9, "Set1")

plot.temp=
  ggplot(UASB_Saline_data, aes(x = Day, y = `Temp (¢XC)`, color = Reactor)) +
  geom_line(size = 1.2) +
  #facet_grid(. ~ Reactor) +
  scale_fill_manual(values = colors1)+
  labs(title = "", x = "Operation day", y = "Temperature (°C)", fill = "") +
  scale_x_continuous(breaks=seq(0,110,20))+
  scale_y_continuous(limits=c(30,35), breaks=seq(30,36,1))+
  theme(axis.text.y = element_text(colour = "black", size = 24), #, 1face = "bold"
        axis.text.x = element_text(colour = "black", size = 24), 
        strip.text = element_text(colour = "black", size = 24),
        legend.text = element_text(size = 24, colour ="black"), 
        legend.title = element_blank(),
        legend.position = "bottom",
        legend.direction = "horizontal", ##Vertical  #horizontal
        axis.title.y = element_text(size = 24 , colour = "black"), 
        axis.title.x = element_text(size = 24, colour = "black"),
        axis.ticks=element_line(colour="black",size=1,linetype=1), #Modify the size and color of the line label on the axis
        axis.ticks.length=unit(0.4,"lines"),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 2),
        strip.background=element_rect(colour=NA,fill=NA), #the 
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key=element_blank())

plot.cond=
  ggplot(UASB_Saline_data, aes(x = Day, y = `EC (mS/cm)`, color = Reactor)) +
  geom_line(size = 1.2) +
  #facet_grid(. ~ Reactor) +
  scale_fill_manual(values = colors1)+
  labs(title = "", x = "Operation day", y = "EC (mS/cm)", fill = "") +
  scale_x_continuous(breaks=seq(0,110,20))+
  scale_y_continuous(limits=c(0,60), breaks=seq(0,60,10))+
  theme(axis.text.y = element_text(colour = "black", size = 24), #, 1face = "bold"
        axis.text.x = element_text(colour = "black", size = 24), 
        strip.text = element_text(colour = "black", size = 24),
        legend.text = element_text(size = 24, colour ="black"), 
        legend.title = element_blank(),
        legend.position = "bottom",
        legend.direction = "horizontal", ##Vertical  #horizontal
        axis.title.y = element_text(size = 24 , colour = "black"), 
        axis.title.x = element_text(size = 24, colour = "black"),
        axis.ticks=element_line(colour="black",size=1,linetype=1), #Modify the size and color of the line label on the axis
        axis.ticks.length=unit(0.4,"lines"),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 2),
        strip.background=element_rect(colour=NA,fill=NA), #the 
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key=element_blank())

P.merge=
  ggarrange(
    plot.temp + theme (plot.margin = margin(t = 0.5, r = 0.5, b = 0, l = 0.5, "cm")),
    plot.cond + theme (plot.margin = margin(t = 0.5, r = 0.5, b = 0, l = 0.5, "cm")),
    align = "hv",
    labels = c("a", "b"),
    ncol=2, nrow=1,
    common.legend = TRUE, legend="bottom",
    font.label = list(size = 30))

ggsave(plot = P.merge,"Figure2.jpg", path=wd, width = 15, height = 6)


# Figure 3

#input the data
KTU.table <- read.csv("KTU.table.summary_relative_abundance.csv")
metadata=read.csv("metadata.csv")

#data cleaning and merge
colnames(KTU.table)[colnames(KTU.table) == "X"] <- "Feature.ID"
metadata <- metadata[-1, ]
metadata <- metadata[ ,-3]

KTU.table.annotate=merge(KTU.table, metadata, by = "Feature.ID")

# adjusting data structure
long_KTU.table.annotate <- pivot_longer(
  KTU.table.annotate,
  cols = -c(Feature.ID, Taxon),  # Exclude these columns
  names_to = "Sample_ID",         # Name for the new column containing variable names
  values_to = "abundance"            # Name for the new column containing values
)

long_KTU.table.annotate$Sample_ID2=long_KTU.table.annotate$Sample_ID
long_KTU.table.annotate <- separate(long_KTU.table.annotate, Sample_ID2, into = c("Reactor", "Day"), sep = "_D")
long_KTU.table.annotate$Day=as.numeric(long_KTU.table.annotate$Day)
long_KTU.table.annotate$Reactor <- gsub("Daughter", "D", long_KTU.table.annotate$Reactor)
long_KTU.table.annotate$abundance=long_KTU.table.annotate$abundance*100
long_KTU.table.annotate$abundance=as.numeric(long_KTU.table.annotate$abundance)

# color code
# Generate a large palette (up to the number of unique Feature.ID values)
# color code# seting 1
num_colors <- length(unique(long_KTU.table.annotate$Feature.ID))
custom_palette <- colorRampPalette(brewer.pal(12, "Set3"))(num_colors)
# color code# seting 2
colors <- brewer.pal(12, "Set3")
color_code1 <- rep(colors, length.out = length(unique(long_KTU.table.annotate$Feature.ID)))
# color code# seting 3
hex <- hue_pal()(16)
color_code2 <- rep(hex, length.out = length(unique(long_KTU.table.annotate$Feature.ID)))
# color code# seting 4
random_colors <- sample(colors(), size = length(unique(long_KTU.table.annotate$Feature.ID)), replace = TRUE)


# Create the taxonomic plot
plot.taxa.bar=
  ggplot(long_KTU.table.annotate, aes(x = Day, y = abundance, fill = Feature.ID)) +
  geom_bar(stat = "identity", position = "stack", width = 1) +
  facet_grid(. ~ Reactor) +
  scale_fill_manual(values = color_code1)+
  labs(title = "", x = "Operation day", y = "Relative Abundance (%)", fill = "Taxa") +
  scale_x_continuous(breaks=seq(0,110,20))+
  theme(axis.text.y = element_text(colour = "black", size = 24), #, 1face = "bold"
        axis.text.x = element_text(colour = "black", size = 24), 
        strip.text = element_text(colour = "black", size = 24),
        legend.text = element_text(size = 24, colour ="black"), 
        legend.title = element_blank(),
        legend.position = "none",
        legend.direction = "horizontal", ##Vertical  #horizontal
        axis.title.y = element_text(size = 24 , colour = "black"), 
        axis.title.x = element_text(size = 24, colour = "black"),
        axis.ticks=element_line(colour="black",size=1,linetype=1), #Modify the size and color of the line label on the axis
        axis.ticks.length=unit(0.4,"lines"),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 2),
        strip.background=element_rect(colour=NA,fill=NA), #the 
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key=element_blank())

plot.taxa.bar
# save the data 
ggsave(plot = plot.taxa.bar,"Figure3.jpg", path=wd, width = 15, height = 6)


# Figure 4

#Impute Missing Value (linear)
acid=UASB_Saline_data[,c('Eff Formate (mg/L)','Eff Acetate (mg/L)','Eff Propionate (mg/L)','Eff Isobutyrate (mg/L)','Eff Butyrate (mg/L)','Eff Isovalerate (mg/L)')]

acid=as.data.frame(na.approx(acid))

acid$Day=UASB_Saline_data$Day
acid$Reactor=UASB_Saline_data$Reactor


# adjusting data structure
long_acid <- pivot_longer(
  acid,
  cols = -c(Day, Reactor),  # select the column
  names_to = "VFA",       # Name for the new column containing variable names
  values_to = "concertation"       # Name for the new column containing values
)

#long_acid$VFA <- gsub(" (mg/L)", "", long_acid$VFA)
long_acid$VFA <- gsub("\\s*\\(mg/L\\)", "", long_acid$VFA)
long_acid$VFA <- gsub("Eff\\s*", "", long_acid$VFA)

long_acid$VFA <- factor(long_acid$VFA, levels = c("Formate","Acetate","Propionate","Isobutyrate","Butyrate","Isovalerate"))

colors1 <- brewer.pal(9, "Set1")

plot.VFA.line=
  ggplot(long_acid, aes(x = Day, y = concertation, color = VFA)) +
  geom_line(size = 1) +
  facet_grid(. ~ Reactor) +
  scale_fill_manual(values = colors1)+
  labs(title = "", x = "Operation day", y = "Concertation (mg/L)", fill = "") +
  scale_x_continuous(breaks=seq(0,110,20))+
  scale_y_continuous(breaks=seq(0,500,100))+
  theme(axis.text.y = element_text(colour = "black", size = 24), #, 1face = "bold"
        axis.text.x = element_text(colour = "black", size = 24), 
        strip.text = element_text(colour = "black", size = 24),
        legend.text = element_text(size = 24, colour ="black"), 
        legend.title = element_blank(),
        legend.position = "bottom",
        legend.direction = "horizontal", ##vertical  #horizontal
        axis.title.y = element_text(size = 24 , colour = "black"), 
        axis.title.x = element_text(size = 24, colour = "black"),
        axis.ticks=element_line(colour="black",size=1,linetype=1), #Modify the size and color of the line label on the axis
        axis.ticks.length=unit(0.4,"lines"),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 2),
        strip.background=element_rect(colour=NA,fill=NA), #the 
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key=element_blank())+
  guides(colour = guide_legend(nrow = 1, byrow = T))

plot.VFA.line
# save the data 
ggsave(plot = plot.VFA.line,"Figure4.jpg", path=wd, width = 15, height = 6)
