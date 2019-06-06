#! /usr/local/bin/Rscript
# Usage: ./102-boxplots.R salmo-salar

# Plot lastz output with error bars
library(tidyverse)
library(matrixStats)

args <- commandArgs(trailingOnly = TRUE)
args<-c("salmo-salar")

protos<-as_tibble(read.table(paste("./data/",args[1],"/",args[1], "-protokaryotype.txt", sep=""))) %>%
  rename(Protokaryotype = V1, Comparison = V2)
protos$Protokaryotype<-as.character(protos$Protokaryotype)
protos$Comparison<-as.character(protos$Comparison)

data<-as_tibble(read.table(paste("./outputs/101/",args[1],"-101-lastz.result", sep="")))
data$V2<-as.character(data$V2)
data$V7<-as.character(data$V7)
data$V12<-as.character(data$V12)

#Split 
data<-data %>% separate(V12, c("Top", "Bottom"))
data$Top<-as.numeric(data$Top)
data$Bottom<-as.numeric(data$Bottom)

#Calculate percent similarity
data<-data %>% mutate(Similarity = Top/Bottom*100)

#Create labels
data<-data %>% mutate(Comparison = paste(V2,V7, sep = "-"))

#Join up protokaryotype
data<- data %>% filter(Comparison %in% protos$Comparison)
data <-left_join(data, protos)

#filter for alignment lengths
data<-data %>% filter(Bottom >= 1000)

#filter for similarity?
#data<-data %>% filter(Similarity >= 80)

#Reduce comparisons to chromosomes of interest
#Sensitive to order of chromos
#data<-data %>% filter(Comparison %in% c("omy01-omy23", "omy02-omy03","omy07-omy18"))

#Calculate the medians to put on the plot later (just for now)
#medians <- aggregate(Similarity ~  Comparison, data, )

#Calculate the means and order data by high -> low
#data <- data %>% group_by(Comparison) %>% mutate(Mean = mean(Similarity)) %>%
#  arrange(desc(Mean))

#Calculate weighted medians
data <- data %>% separate(V14, sep="[/]", into = c("AlignmentLength","Length"))
data$AlignmentLength <- as.numeric(data$AlignmentLength)
data <- data %>% group_by(Comparison) %>% mutate(Median = weightedMedian(Similarity, w=AlignmentLength))
data$Comparison <- reorder(data$Comparison, desc(data$Median))
data$Protokaryotype <- reorder(data$Protokaryotype, desc(data$Median))

#Sample Sizes?
samplesize <- data %>% group_by(Comparison) %>% count_()
samplesize <- left_join(samplesize, protos)
#samplesize <- samplesize %>% filter(n >= 500)
#data <- data %>% filter(Comparison %in% samplesize$Comparison)

pdf(paste("./outputs/102/",args[1],"-boxplots.pdf", sep=""), width =11, height = 8.5/2)

ggplot(data)+geom_boxplot(aes(x=Protokaryotype, y=Similarity, weight=AlignmentLength),
                          outlier.size=0.1, outlier.alpha=0.5,
                          outlier.color="gray50", outlier.shape=15)+
  theme_classic()+
  ylab("Perecent Similarity")+
  theme(axis.text.x= element_text(angle=45,hjust=1, face = "bold"))+
  geom_text(data = samplesize, aes(label = Comparison,
                             x = Protokaryotype, y = 100 + 2), size=2.5, angle=45)+
  ylim(75,105)
  #geom_text(data = samplesize, aes(label = n, x=Protokaryotype,
  #                                 y = 100+1), size=3)

dev.off()

#Save output as .rda renamed by species
assign(paste(args[1]), data)

save(list=paste(args[1]), file=paste("./outputs/102/",args[1],".rda",sep=""))

