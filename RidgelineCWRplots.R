#Load packages
require(ape)
library(stringr)
library(ggridges)
library(ggplot2)
library(forcats)

#Import aligned FASTA file
fas <- read.FASTA("FASTA_files/genbank_query_theocwrs_its_edited_trimmed.fasta")

#Infer pairwise genetic distance
distDNA <- dist.dna(fas)
dist <- as.matrix(distDNA)

# Convert distance matrix into 3 cols: Seq1, Seq2, Dist
OUT <- NULL
for(i in 1:nrow(dist)){
  tmp <- dist[i,]
  tmp <- cbind(rep(rownames(dist)[i], length(tmp)), names(tmp), as.vector(tmp))
  OUT <- rbind(OUT, tmp)
}
distSimp <- as.data.frame(OUT)
colnames(distSimp) <- c("Seq1", "Seq2", "Dist")

#Add 2 cols with species
distSimp$Sp1 <- paste(sapply(strsplit(as.vector(distSimp$Seq1), split="_"), "[[", 2), sapply(strsplit(as.vector(distSimp$Seq1), split="_"), "[[", 3), sep="_") 
distSimp$Sp2 <- paste(sapply(strsplit(as.vector(distSimp$Seq2), split="_"), "[[", 2), sapply(strsplit(as.vector(distSimp$Seq2), split="_"), "[[", 3), sep="_") 

#Want a matrix with sp1 (sp) and cwr

cwr <- "Theobroma_cacao"
sp <- unique(distSimp$Sp1)

DistInput <- NULL
for(i in 1:length(sp)){
  tmp <- subset(distSimp, distSimp$Sp1 == sp[i] & distSimp$Sp2 == cwr)
  DistInput <- rbind(DistInput, tmp)
}

#Order sp for plot
orderPlot <- aggregate(as.numeric(Dist) ~ Sp1, min, data = DistInput)
orderPlot <- orderPlot[order(orderPlot[,2], decreasing = F),]

DistInput$SpOrd <- rep("NA", nrow(DistInput))
for(i in 1:length(orderPlot$Sp1)){
  DistInput$SpOrd[grep(orderPlot$Sp1[i], DistInput$Sp1)] <- i
}
#Pad tip numbers to allow sorting them (for ridgeline plot)
DistInput$SpOrd <- str_pad(DistInput$SpOrd, 2, pad = "0")


CWRsp2 <- DistInput %>%
  mutate(text = fct_reorder(Sp1, as.numeric(SpOrd))) %>%
  ggplot(aes(x = as.numeric(Dist), y = text, fill = text)) +
  #scale_fill_manual(values = c("grey", "blue","pink")) +
  #geom_density_ridges(stat="binline", bins = 30) +
  geom_density_ridges() +
  theme_ridges() + 
  xlab("Genetic distance") +
  ylab(paste("Species vs.", gsub("_", " ", cwr), sep=" ")) +
  theme(legend.position = "none", text = element_text(size = 7), axis.text.y = element_text(size = 7), axis.text.x = element_text(size = 7))  

pdf("Theobroma_cacao.pdf")
CWRsp2
dev.off()
