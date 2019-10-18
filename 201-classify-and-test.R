#! /usr/local/bin/RScript
# Usage
# ./201-classify-and-test.R salmo-salar

library(tidyverse)
library(ggvis)
library(class)
library(gmodels)
library(viridis)
library(matrixStats)
library(caret)

args <- commandArgs(trailingOnly = TRUE)

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
data<- data %>% dplyr::filter(Comparison %in% protos$Comparison)
data <-left_join(data, protos)

#filter for alignment lengths
data<-data %>% dplyr::filter(Bottom >= 1000)

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


#Create summary object
summary <- data %>% group_by(Protokaryotype, Median) %>% mutate(Mean=mean(Similarity)) %>% 
  select(Protokaryotype, Mean, Median) %>% ungroup() %>% group_by(Protokaryotype, Mean, Median) %>% 
  summarize()

#create a training set

head<-head(summary, 4)
tail<-tail(summary, 4)
training<-rbind(head,tail)

train.labels<-factor(c("Tetrasomic","Tetrasomic","Tetrasomic", "Tetrasomic",
                       "Disomic","Disomic","Disomic","Disomic"))
training$Type<-train.labels


#Identifying the best k

trControl <- trainControl(method = "repeatedcv",
                          number = 100,
                          repeats = 10)
fit <- train(Type ~ .,
             method     = "knn",
             tuneGrid   = expand.grid(k = 1:10),
             trControl  = trControl,
             metric     = "Accuracy",
             data       = training)

pred <- knn(train = training[, 3, drop=FALSE], test = summary[, 3, drop=FALSE], cl = train.labels, k=fit$finalModel$k, prob = TRUE)

#How did it work?
summary$Prediction <- pred
summary$Probability <-round(attr(pred, "prob"),2)


#Now to generate boxplots as well
total<-left_join(data, summary)

line<-total %>% select(Protokaryotype, Median) %>% unique() %>% ungroup() %>% mutate(Mean=mean(Median))
yint<-line$Mean  %>% unique()

#Now to test for differences
disomic<-filter(total, Prediction=="Disomic")
tetrasomic<-filter(total, Prediction=="Tetrasomic")
kw <- wilcox.test(disomic$Similarity, tetrasomic$Similarity)

#Now to get a name that is in line with the paper

taxon <- ifelse(args[1]=="t-thymallus", "T. thy",
         ifelse(args[1]=="salmo-salar", "S. sal",
         ifelse(args[1]=="s-alpinus", "S. alp",
         ifelse(args[1]=="o-mykiss", "O. myk",
         ifelse(args[1]=="o-tshaw", "O. tsh",
         ifelse(args[1]=="o-kisutch", "O. kis",
                "needs proper name"))))))
pdf(paste("./outputs/201/",args[1],"-accuracy-v-neighbors.pdf", sep=""), width=7, height=5)
ggplot(fit)
dev.off()

pdf(paste("./outputs/201/",args[1],"-predictions.pdf", sep=""), width=8.5, height=11/6)

ggplot(summary)+geom_bar(aes(x=Protokaryotype, y=Median, fill=Prediction), color="black", stat="identity", alpha=0.5) + 
  scale_fill_viridis_d(direction=-1) + 
  geom_text(aes(x=Protokaryotype, y=(Median+3)), label=summary$Probability, size=3)+
  coord_cartesian(ylim=c(75,100))+
  theme_classic()+
  theme(axis.text.x=element_text(angle=45, hjust=1, face="bold"))+
  theme(axis.title.x = element_blank())+
  theme(axis.title.y = element_blank())+
  ggtitle(taxon)+
  theme(plot.title = element_text(hjust = 0.5, size=12, face="bold"))
  



dev.off()

pdf(paste("./outputs/201/",args[1],"-boxplots.pdf", sep=""), width =8.5, height = 11/6)
ggplot(total)+geom_boxplot(aes(x=Protokaryotype, y=Similarity, weight=AlignmentLength, fill=Prediction),
                           outlier.size=0.05, outlier.alpha=0.5, outlier.shape=15,
                           outlier.stroke=0.25,
                           alpha=0.5)+
  theme_classic()+
  ylab("Perecent Similarity")+
  theme(axis.text.x= element_text(angle=45,hjust=1, face = "bold"))+
  geom_hline(yintercept = yint, alpha=0.75, linetype="dashed", color="grey50")+
  geom_text(data = samplesize, aes(label = n, x=Protokaryotype,
                                   y = 100+2), size=2)+
  #geom_text(data = samplesize, aes(label = Comparison,
  #                           x = Protokaryotype, y = 100 + 2), size=2.5, angle=45)+
  ylim(75,103)+
  scale_fill_viridis_d(direction=-1)+
  theme(legend.position = "none")+
  ggtitle(taxon)+
  #ggtitle(paste(args[1], "W =", round(kw$statistic,2), "p-value =", kw$p.value, sep=" "))+
  theme(plot.title = element_text(hjust = 0.5, size=12, face="bold"))+
  theme(axis.title.x = element_blank())+
  theme(axis.title.y = element_blank())

dev.off()

