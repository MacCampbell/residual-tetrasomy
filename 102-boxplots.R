#! /usr/local/bin/Rscript
# Usage: ./102-boxplots.R salmo-salar

# Plot lastz output with error bars
library(tidyverse)
library(matrixStats)

args <- commandArgs(trailingOnly = TRUE)
#args<-c("salmo-salar")

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

#Sample Sizes?
samplesize <- data %>% group_by(Comparison) %>% count_()


#Example code
#ggplot(tetDf)+geom_boxplot(aes(Block,PID))+
#  labs(x="Tetrasomic Pairing",y="Percent ID")+
# coord_cartesian(ylim = c(80,100))+
#theme_classic()+
#theme(text = element_text(face="bold", size=18))+
#theme(axis.text.x= element_text(face="bold", size=12,angle=45,hjust=1))+
#theme(axis.text.y= element_text(face="bold", size=12))+
#scale_x_discrete(limits=meanDf$Block)

pdf(paste("./outputs/102/",args[1],"-boxplots.pdf", sep=""), width =11/2, height = 8.5/2)

ggplot(data)+geom_boxplot(aes(x=Comparison, y=Similarity, weight=AlignmentLength),
                          outlier.size=0.1, outlier.alpha=0.5,
                          outlier.color="gray50", outlier.shape=15)+
  theme_classic()+
  ylab("Perecent Similarity")+
  theme(axis.text.x= element_text(angle=45,hjust=1))+
  #geom_text(data = means, aes(label = round(Similarity,2),
   #                           x = Comparison, y = Similarity + 2)) +
  geom_text(data = samplesize, aes(label = n, x=Comparison,
                                   y = 100+1))

dev.off()

#Save output as .rda renamed by species
assign(paste(args[1]), data)

save(list=paste(args[1]), file=paste("./outputs/102/",args[1],".rda",sep=""))

