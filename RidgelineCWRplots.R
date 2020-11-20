###------
#Infer K2P distances between CWRs and crop species
#and display these interspecific distances using ridgeline plots
###------

###~~~
#Load packages
###~~~
library(ape)
library(stringr)
library(ggridges)
library(ggplot2)
library(forcats)
library(dplyr)
library(ggpubr)
library(grid)

####---
# Theobroma cacao
####---
###~~~
#Import aligned FASTA file
###~~~
fas <- read.FASTA("FASTA_files/genbank_query_theocwrs_its_edited_trimmed.fasta")

###~~~
#Infer Kimura's 2-parameters distances
###~~~
#This model is also known as K80
distDNA <- dist.dna(fas, model = "K80")
dist <- as.matrix(distDNA)

###~~~
#Convert pairwise distance matrix into 3 cols: Seq1, Seq2, Dist
###~~~
OUT <- NULL
for(i in 1:nrow(dist)){
  tmp <- dist[i,]
  tmp <- cbind(rep(rownames(dist)[i], length(tmp)), names(tmp), as.vector(tmp))
  OUT <- rbind(OUT, tmp)
}
distSimp <- as.data.frame(OUT)
colnames(distSimp) <- c("Seq1", "Seq2", "Dist")

###~~~
#Add 2 cols with species names (to infer intersp. dist.)
###~~~
distSimp$Sp1 <- paste(sapply(strsplit(as.vector(distSimp$Seq1), split="_"), "[[", 2), sapply(strsplit(as.vector(distSimp$Seq1), split="_"), "[[", 3), sep=" ") 
distSimp$Sp2 <- paste(sapply(strsplit(as.vector(distSimp$Seq2), split="_"), "[[", 2), sapply(strsplit(as.vector(distSimp$Seq2), split="_"), "[[", 3), sep=" ") 

###~~~
#Tidy matrix to infer intraspecific distances sp vs sp
###~~~

#List of species
sp <- unique(distSimp$Sp1)

#Intrasp
DistIntra <- NULL
for(i in 1:length(sp)){
  tmp <- subset(distSimp, distSimp$Sp1 == sp[i] & distSimp$Sp2 == sp[i])
  DistIntra <- rbind(DistIntra, tmp)
}

#Remove species with < 2 DNA sequences
NSeqSp <- table(DistIntra$Sp1)
ToRm <- names(NSeqSp)[which(NSeqSp <= 1)]

if(length(ToRm) > 0){
  DistIntra <- subset(DistIntra, !(DistIntra$Sp1 %in% ToRm))
}

###~~~
#Tidy matrix to infer interspecific distances between CWRs and crop
###~~~
#We only want each species (=CWR) vs. crop
#For the crop we do infraspecific distance (crop vs. crop)
crop <- "Theobroma cacao"

#Interspecific distances 
DistInter <- NULL
for(i in 1:length(sp)){
  tmp <- subset(distSimp, distSimp$Sp1 == sp[i] & distSimp$Sp2 == crop)
  DistInter <- rbind(DistInter, tmp)
}

#Remove crop against crop
DistInter <- subset(DistInter, !(DistInter$Sp1 %in% crop))

###~~~
#Order species for plots
###~~~
#We use mean distances to sort species (but could change that to mean or max) 

##
#Intra
orderPlotIntra <- aggregate(as.numeric(as.vector(Dist)) ~ Sp1, mean, data = DistIntra)
orderPlotIntra <- orderPlotIntra[order(orderPlotIntra[,2], decreasing = F),]

#Add species order in DistIntra
DistIntra$SpOrd <- rep("NA", nrow(DistIntra))
for(i in 1:length(orderPlotIntra$Sp1)){
  DistIntra$SpOrd[grep(orderPlotIntra$Sp1[i], DistIntra$Sp1)] <- i
}
#Pad tip numbers to allow sorting them (for ridgeline plot)
DistIntra$SpOrd <- str_pad(DistIntra$SpOrd, 2, pad = "0")

##
#Inter
orderPlotInter <- aggregate(as.numeric(as.vector(Dist)) ~ Sp1, mean, data = DistInter)
orderPlotInter <- orderPlotInter[order(orderPlotInter[,2], decreasing = F),]

#Add species order in DistInter
DistInter$SpOrd <- rep("NA", nrow(DistInter))
for(i in 1:length(orderPlotInter$Sp1)){
  DistInter$SpOrd[grep(orderPlotInter$Sp1[i], DistInter$Sp1)] <- i
}
#Pad tip numbers to allow sorting them (for ridgeline plot)
DistInter$SpOrd <- str_pad(DistInter$SpOrd, 2, pad = "0")

###~~~
#Draw plot: Ridgeline plots
###~~~

###
#Intra
IntraPlot <- DistIntra %>%
  mutate(text = fct_reorder(Sp1, as.numeric(SpOrd))) %>%
  ggplot(aes(x = as.numeric(as.vector(Dist)), y = text, fill = text)) +
  #scale_fill_manual(values = c("grey", "blue","pink")) +
  #geom_density_ridges(stat="binline", bins = 30) +
  geom_density_ridges() +
  theme_ridges() + 
  xlab("Kimura's 2-parameters genetic distance") +
  ylab("Intraspecific comparision          ") +
  #geom_segment(mapping=aes(x=-0.2, y=2, xend=-0.2, yend=9.5)) +
  theme(legend.position = "none", text = element_text(size = 8), axis.text.y = element_text(size = 7, face = "italic"), axis.text.x = element_text(size = 7))  

###
#Inter
InterPlot <- DistInter %>%
  mutate(text = fct_reorder(Sp1, as.numeric(SpOrd))) %>%
  ggplot(aes(x = as.numeric(as.vector(Dist)), y = text, fill = text)) +
  #scale_fill_manual(values = c("grey", "blue","pink")) +
  #geom_density_ridges(stat="binline", bins = 30) +
  geom_density_ridges() +
  theme_ridges() + 
  xlab("Kimura's 2-parameters genetic distance") +
  ylab("Interspecific comparision (CWRs vs. crop)          ") +
  #geom_segment(mapping=aes(x=-0.2, y=2, xend=-0.2, yend=9.5)) +
  theme(legend.position = "none", text = element_text(size = 8), axis.text.y = element_text(size = 7, face = "italic"), axis.text.x = element_text(size = 7))  

###~~~
#Arrange and save plots as pdf
###~~~
pdf("Theobroma_cacao_K2P.pdf")
grid.newpage()
# Create layout : nrow = 3, ncol = 3
pushViewport(viewport(layout = grid.layout(nrow = 1, ncol = 5)))
# A helper function to define a region on the layout
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
} 

# Arrange plots
print(IntraPlot, vp = define_region(row = 1, col = 1:2))
print(InterPlot, vp = define_region(row = 1, col = 3:5))
dev.off()

pdf("Theobroma_cacao_K2P.pdf")
ggarrange(IntraPlot, InterPlot, ncol = 2, nrow = 1, labels="auto")
dev.off()


####---
# Vanilla planifolia
####---
###~~~
#Import aligned FASTA file
###~~~
fas <- read.FASTA("FASTA_files/Vanilla_5.8S_2020_11_09_GenBank_edited_trimmed.fasta")

###~~~
#Infer Kimura's 2-parameters distances
###~~~
#This model is also known as K80
distDNA <- dist.dna(fas, model = "K80")
dist <- as.matrix(distDNA)

###~~~
#Convert pairwise distance matrix into 3 cols: Seq1, Seq2, Dist
###~~~
OUT <- NULL
for(i in 1:nrow(dist)){
  tmp <- dist[i,]
  tmp <- cbind(rep(rownames(dist)[i], length(tmp)), names(tmp), as.vector(tmp))
  OUT <- rbind(OUT, tmp)
}
distSimp <- as.data.frame(OUT)
colnames(distSimp) <- c("Seq1", "Seq2", "Dist")

#Find hybrids and replace names
distSimp$Seq1 <- gsub("Vanilla_planifolia_x_Vanilla_pompona", "Vanilla_planXpomp", distSimp$Seq1)
distSimp$Seq2 <- gsub("Vanilla_planifolia_x_Vanilla_pompona", "Vanilla_planXpomp", distSimp$Seq2)

###~~~
#Add 2 cols with species names (to infer intersp. dist.)
###~~~
distSimp$Sp1 <- paste(sapply(strsplit(as.vector(distSimp$Seq1), split="_"), "[[", 2), sapply(strsplit(as.vector(distSimp$Seq1), split="_"), "[[", 3), sep=" ") 
distSimp$Sp2 <- paste(sapply(strsplit(as.vector(distSimp$Seq2), split="_"), "[[", 2), sapply(strsplit(as.vector(distSimp$Seq2), split="_"), "[[", 3), sep=" ") 

###~~~
#Tidy matrix to infer intraspecific distances sp vs sp
###~~~

#List of species
sp <- unique(distSimp$Sp1)

#Intrasp
DistIntra <- NULL
for(i in 1:length(sp)){
  tmp <- subset(distSimp, distSimp$Sp1 == sp[i] & distSimp$Sp2 == sp[i])
  DistIntra <- rbind(DistIntra, tmp)
}

#Remove species with < 2 DNA sequences
NSeqSp <- table(DistIntra$Sp1)
ToRm <- names(NSeqSp)[which(NSeqSp <= 1)]

if(length(ToRm) > 0){
  DistIntra <- subset(DistIntra, !(DistIntra$Sp1 %in% ToRm))
}

###~~~
#Tidy matrix to infer interspecific distances between CWRs and crop
###~~~
#We only want each species (=CWR) vs. crop
#For the crop we do infraspecific distance (crop vs. crop)
crop <- "Vanilla planifolia"

#Interspecific distances 
DistInter <- NULL
for(i in 1:length(sp)){
  tmp <- subset(distSimp, distSimp$Sp1 == sp[i] & distSimp$Sp2 == crop)
  DistInter <- rbind(DistInter, tmp)
}

#Remove crop against crop
DistInter <- subset(DistInter, !(DistInter$Sp1 %in% crop))

###~~~
#Order species for plots
###~~~
#We use mean distances to sort species (but could change that to mean or max) 

##
#Intra
orderPlotIntra <- aggregate(as.numeric(as.vector(Dist)) ~ Sp1, mean, data = DistIntra)
orderPlotIntra <- orderPlotIntra[order(orderPlotIntra[,2], decreasing = F),]

#Add species order in DistIntra
DistIntra$SpOrd <- rep("NA", nrow(DistIntra))
for(i in 1:length(orderPlotIntra$Sp1)){
  DistIntra$SpOrd[grep(orderPlotIntra$Sp1[i], DistIntra$Sp1)] <- i
}
#Pad tip numbers to allow sorting them (for ridgeline plot)
DistIntra$SpOrd <- str_pad(DistIntra$SpOrd, 2, pad = "0")

##
#Inter
orderPlotInter <- aggregate(as.numeric(as.vector(Dist)) ~ Sp1, mean, data = DistInter)
orderPlotInter <- orderPlotInter[order(orderPlotInter[,2], decreasing = F),]

#Add species order in DistInter
DistInter$SpOrd <- rep("NA", nrow(DistInter))
for(i in 1:length(orderPlotInter$Sp1)){
  DistInter$SpOrd[grep(orderPlotInter$Sp1[i], DistInter$Sp1)] <- i
}
#Pad tip numbers to allow sorting them (for ridgeline plot)
DistInter$SpOrd <- str_pad(DistInter$SpOrd, 2, pad = "0")

###~~~
#Draw plot: Ridgeline plots
###~~~

###
#Intra
IntraPlot <- DistIntra %>%
  mutate(text = fct_reorder(Sp1, as.numeric(SpOrd))) %>%
  ggplot(aes(x = as.numeric(as.vector(Dist)), y = text, fill = text)) +
  #scale_fill_manual(values = c("grey", "blue","pink")) +
  #geom_density_ridges(stat="binline", bins = 30) +
  geom_density_ridges() +
  theme_ridges() + 
  xlab("Kimura's 2-parameters genetic distance") +
  ylab("Intraspecific comparision          ") +
  #geom_segment(mapping=aes(x=-0.2, y=2, xend=-0.2, yend=9.5)) +
  theme(legend.position = "none", text = element_text(size = 8), axis.text.y = element_text(size = 7, face = "italic"), axis.text.x = element_text(size = 7))  

###
#Inter
InterPlot <- DistInter %>%
  mutate(text = fct_reorder(Sp1, as.numeric(SpOrd))) %>%
  ggplot(aes(x = as.numeric(as.vector(Dist)), y = text, fill = text)) +
  #scale_fill_manual(values = c("grey", "blue","pink")) +
  #geom_density_ridges(stat="binline", bins = 30) +
  geom_density_ridges() +
  theme_ridges() + 
  xlab("Kimura's 2-parameters genetic distance") +
  ylab("Interspecific comparision (CWRs vs. crop)          ") +
  #geom_segment(mapping=aes(x=-0.2, y=2, xend=-0.2, yend=9.5)) +
  theme(legend.position = "none", text = element_text(size = 8), axis.text.y = element_text(size = 7, face = "italic"), axis.text.x = element_text(size = 7))  

###~~~
#Arrange and save plots as pdf
###~~~

pdf("Vanilla_planifolia_K2P.pdf")
grid.newpage()
# Create layout : nrow = 3, ncol = 3
pushViewport(viewport(layout = grid.layout(nrow = 1, ncol = 5)))
# A helper function to define a region on the layout
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
} 

# Arrange plots
print(IntraPlot, vp = define_region(row = 1, col = 1:2))
print(InterPlot, vp = define_region(row = 1, col = 3:5))
dev.off()

ggarrange(IntraPlot, InterPlot, ncol = 2, nrow = 1, labels="auto")





