# Don't have PK info for arctic charr/grayling
library(tidyverse)
library(ggvis)
library(class)
library(gmodels)
library(viridis)
library(matrixStats)
library(caret)

args<-c("t-thymallus")

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
data$Comparison<-gsub("\\D","", data$V2)
data$Comparison<-as.numeric(data$Comparison)

#Previously did this:
#data<-data %>% mutate(Comparison = paste(V2,V7, sep = "-"))


#filter for alignment lengths
data<-data %>% dplyr::filter(Bottom >= 1000)


#Calculate weighted medians
data <- data %>% separate(V14, sep="[/]", into = c("AlignmentLength","Length"))
data$AlignmentLength <- as.numeric(data$AlignmentLength)
data <- data %>% group_by(Comparison) %>% mutate(Median = weightedMedian(Similarity, w=AlignmentLength))
data$Comparison <- reorder(data$Comparison, desc(data$Median))


#Sample Sizes?
samplesize <- data %>% group_by(Comparison) %>% count_()


#Handy normalization function if needed
#normalize <- function(x) {
# num <- x - min(x)
#denom <- max(x) - min(x)
#return (num/denom)
#}

#Using "summary" object from previous script
summary <- data %>% group_by(Comparison, Median) %>% mutate(Mean=mean(Similarity)) %>% 
  select(Comparison, Mean, Median) %>% ungroup() %>% group_by(Comparison, Mean, Median) %>% 
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

line<-total %>% select(Comparison, Median) %>% unique() %>% ungroup() %>% mutate(Mean=mean(Median))
yint<-line$Mean  %>% unique()

#Now to test for differences

disomic<-filter(total, Prediction=="Disomic")
tetrasomic<-filter(total, Prediction=="Tetrasomic")
kw <- wilcox.test(disomic$Similarity, tetrasomic$Similarity)

pdf(paste("./outputs/108/", args[1], "-accuracy-v-neighbors.pdf",sep=""), width=7, height=5)
ggplot(fit)
dev.off()

pdf(paste("./outputs/108/",args[1],"-predictions.pdf", sep=""), width=7, height=5)

ggplot(summary)+geom_bar(aes(x=Comparison, y=Median, fill=Prediction), color="black", stat="identity", alpha=0.5) + 
  scale_fill_viridis_d(direction=-1) + 
  geom_text(aes(x=Comparison, y=(Median+5)), label=summary$Probability, size=2)+
  theme_classic()+
  theme(axis.text.x=element_text(angle=45, hjust=1))+
  xlab("Protokaryotype")



dev.off()

pdf(paste("./outputs/108/",args[1],"-boxplots.pdf", sep=""), width =8.5, height = 11/6)
ggplot(total)+geom_boxplot(aes(x=Comparison, y=Similarity, weight=AlignmentLength, fill=Prediction),
                           outlier.size=0.05, outlier.alpha=0.5, outlier.shape=15,
                           outlier.stroke=0.25,
                           alpha=0.5)+
  theme_classic()+
  ylab("Perecent Similarity")+
  theme(axis.text.x= element_text(angle=45,hjust=1, face = "bold"))+
  geom_hline(yintercept = yint, alpha=0.75, linetype="dashed", color="grey50")+
  geom_text(data = samplesize, aes(label = n, x=Comparison,
                                   y = 100+2), size=2)+
  #geom_text(data = samplesize, aes(label = Comparison,
  #                           x = Protokaryotype, y = 100 + 2), size=2.5, angle=45)+
  ylim(75,103)+
  scale_fill_viridis_d(direction=-1)+
  theme(legend.position = "none")+
  ggtitle(paste(args[1], "W =", round(kw$statistic,2), "p-value =", kw$p.value, sep=" "))+
  ggtitle("T. thy")+
  theme(plot.title = element_text(hjust = 0.5, size=8))+
  theme(axis.title.x = element_blank())+
  theme(axis.title.y = element_blank())
  #theme(axis.text.x=element_text(angle=45, hjust=1))+
  #xlab("Protokaryotype")

dev.off()