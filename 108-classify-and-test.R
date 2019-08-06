#! /usr/local/bin/RScript
# Usage
# ./108-classify-and-test.R salmo-salar

library(tidyverse)
library(ggvis)
library(class)
library(gmodels)
library(viridis)
library(matrixStats)
library(caret)
### This is recycled
args <- commandArgs(trailingOnly = TRUE)
#args<-c("s-alpinus")

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

eight<-c(2,20,22,9,11,23,25,1)
tetra<-c(2,20,22,9,11,23,25)
data <- data %>% mutate(Type = ifelse(Protokaryotype %in% eight, "Tetrasomic", "Disomic"))


#Sample Sizes?
samplesize <- data %>% group_by(Comparison) %>% count_()
samplesize <- left_join(samplesize, protos)


#Handy normalization function if needed
#normalize <- function(x) {
# num <- x - min(x)
#denom <- max(x) - min(x)
#return (num/denom)
#}

#Using "summary" object from previous script
summary <- data %>% group_by(Protokaryotype, Median) %>% mutate(Mean=mean(Similarity)) %>% 
  select(Protokaryotype, Mean, Median) %>% ungroup() %>% group_by(Protokaryotype, Mean, Median) %>% 
  summarize()

#summary
#> summary
#Protokaryotype     Mean   Median
#1              11 92.17410 93.94942

#create a training set
#To do, create a way to generate it programmatically instead of
#training<-summary %>% dplyr::filter(Protokaryotype %in% c(16,21,24,23,11,20,2,9))

head<-head(summary, 4)
tail<-tail(summary, 4)
training<-rbind(head,tail)

train.labels<-factor(c("Tetrasomic","Tetrasomic","Tetrasomic", "Tetrasomic",
                       "Disomic","Disomic","Disomic","Disomic"))
training$Type<-train.labels


#Can we identify the best k?
#trControl <- trainControl(method  = "cv",
#                          number  = 10)

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

pdf(paste("./outputs/108/",args[1],"-predictions.pdf", sep=""), width=7, height=5)
ggplot(fit)

ggplot(summary)+geom_bar(aes(x=Protokaryotype, y=Median, fill=Prediction), color="black", stat="identity", alpha=0.5) + 
  scale_fill_viridis_d(direction=-1) + 
  geom_text(aes(x=Protokaryotype, y=(Median+5)), label=summary$Probability, size=2)+
  theme_classic()


dev.off()

pdf(paste("./outputs/108/",args[1],"-boxplots.pdf", sep=""), width =11, height = 8.5/2)
ggplot(total)+geom_boxplot(aes(x=Protokaryotype, y=Similarity, weight=AlignmentLength, fill=Prediction),
                           outlier.size=0.1, outlier.alpha=0.5, outlier.shape=15,
                           alpha=0.5)+
  theme_classic()+
  ylab("Perecent Similarity")+
  theme(axis.text.x= element_text(angle=45,hjust=1, face = "bold"))+
  geom_hline(yintercept = yint, alpha=0.75, linetype="dashed", color="grey50")+
  geom_text(data = samplesize, aes(label = n, x=Protokaryotype,
                                   y = 100+1), size=3)+
  #geom_text(data = samplesize, aes(label = Comparison,
  #                           x = Protokaryotype, y = 100 + 2), size=2.5, angle=45)+
  ylim(75,103)+
  scale_fill_viridis_d(direction=-1)+
  theme(legend.position = "none")+
  ggtitle(paste(args[1], "W =", round(kw$statistic,2), "p-value =", kw$p.value, sep=" "))+
  theme(plot.title = element_text(hjust = 0.5))

dev.off()

