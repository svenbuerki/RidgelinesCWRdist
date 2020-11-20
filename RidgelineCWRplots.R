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
#Tidy matrix to infer interspecific distances between CWRs and crop
###~~~
#We only want each species (=CWR) vs. crop
#For the crop we do infraspecific distance (crop vs. crop)
crop <- "Theobroma cacao"
sp <- unique(distSimp$Sp1)

DistInput <- NULL
for(i in 1:length(sp)){
  tmp <- subset(distSimp, distSimp$Sp1 == sp[i] & distSimp$Sp2 == crop)
  DistInput <- rbind(DistInput, tmp)
}

###~~~
#Order species (or CWRs) for plot
###~~~
#We use min distances to sort species (but could change that to mean or max) 
orderPlot <- aggregate(as.numeric(as.vector(Dist)) ~ Sp1, min, data = DistInput)
orderPlot <- orderPlot[order(orderPlot[,2], decreasing = F),]

#Add species order in DistInput
DistInput$SpOrd <- rep("NA", nrow(DistInput))
for(i in 1:length(orderPlot$Sp1)){
  DistInput$SpOrd[grep(orderPlot$Sp1[i], DistInput$Sp1)] <- i
}
#Pad tip numbers to allow sorting them (for ridgeline plot)
DistInput$SpOrd <- str_pad(DistInput$SpOrd, 2, pad = "0")

###~~~
#Draw plot: Ridgeline plots
###~~~
ylabtitle <- expression(paste("CWRs vs.", italic("Theobroma cacao"), sep=" "))

CWRsp2 <- DistInput %>%
  mutate(text = fct_reorder(Sp1, as.numeric(SpOrd))) %>%
  ggplot(aes(x = as.numeric(as.vector(Dist)), y = text, fill = text)) +
  #scale_fill_manual(values = c("grey", "blue","pink")) +
  #geom_density_ridges(stat="binline", bins = 30) +
  geom_density_ridges() +
  theme_ridges() + 
  xlab("Kimura's 2-parameters genetic distance") +
  ylab(ylabtitle) +
  #geom_segment(mapping=aes(x=-0.2, y=2, xend=-0.2, yend=9.5)) +
  theme(legend.position = "none", text = element_text(size = 8), axis.text.y = element_text(size = 7, face = "italic"), axis.text.x = element_text(size = 7))  

###~~~
#Export plot as pdf
###~~~
pdf(paste(gsub(" ", "_", crop), "_vs_CWRs.pdf", sep=''))
CWRsp2
dev.off()
