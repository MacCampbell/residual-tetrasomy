#! /usr/local/bin/Rscript
# Usage: 
# Goal: To make a series of pairwise comparisons to place like with like

library(tidyverse)
library(matrixStats)
library(viridis)
library(mixtools)
#args <- commandArgs(trailingOnly = TRUE)
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

eight<-c(2,20,22,9,11,23,25,1)

data <- data %>% mutate(Type = ifelse(Protokaryotype %in% eight, "Tetrasomic", "Disomic"))


#Sample Sizes?
samplesize <- data %>% group_by(Comparison) %>% count_()
samplesize <- left_join(samplesize, protos)

#Generate all combinations
combos<-combn(levels(data$Protokaryotype),2)


Compare <- function(combo) {
  v1 <- filter(data, Protokaryotype == combo[1])
  v2 <- filter(data, Protokaryotype == combo[2])
 # result <- wilcox.test(v1$Similarity, v2$Similarity)
 # return(result$p.value)
  df <- rbind(v1,v2)
  result<-kruskal.test(df$Similarity, df$Protokaryotype)
  return(result$p.value)
}

compared<-apply(combos, 2, Compare)

## Can I get this somehow into matrix form?
x <- as_tibble(setNames(data.frame(matrix(ncol=3, nrow=0)), c("Chrom1", "Chrom2", "pvalue")))
for (i in (1:length(compared)) ) {
  x<-add_row(x, Chrom1=combos[1,i], Chrom2=combos[2,i], pvalue=compared[i])
}
compared<-apply(combos, 2, Compare)
nonsig <- x %>% filter(pvalue >= 0.005)

summary <- data %>% group_by(Protokaryotype, Median) %>% mutate(Mean=mean(Similarity)) %>% 
  select(Protokaryotype, Mean, Median) %>% ungroup() %>% group_by(Protokaryotype, Mean, Median) %>% 
  summarize()

##A parametric type approach

Medians <- function(combo) {
  v1 <- filter(summary, Protokaryotype == combo[1])
  v2 <- filter(summary, Protokaryotype == combo[2])
  # result <- wilcox.test(v1$Similarity, v2$Similarity)
  # return(result$p.value)
  stat<-abs(v1$Median - v2$Median)
  return(stat)
}

summarized <- apply(combos, 2, Medians)

y <- as_tibble(setNames(data.frame(matrix(ncol=3, nrow=0)), c("Chrom1", "Chrom2", "Difference")))
for (i in (1:length(summarized)) ) {
  y<-add_row(y, Chrom1=combos[1,i], Chrom2=combos[2,i], Difference=summarized[i])
}

highdist<-filter(y, Difference>6)


ggplot(y, aes(x=Difference))+geom_histogram() + theme_bw() + ylab("Count")
model = normalmixEM(y$Difference, k=2)
estimators<-as.data.frame(model[c("lambda", "mu", "sigma")])

model1 = normalmixEM(y$Difference, k=1)
model2 = normalmixEM(y$Difference, k=2)
model3 = normalmixEM(y$Difference, k=3)
model4 = normalmixEM(y$Difference, k=4)
model5 = normalmixEM(y$Difference, k=5)
plot(model, which=2)


ggplot(y, aes(x=Difference))+geom_histogram((aes(y = (..count..)/sum(..count..)))) +
  stat_function(fun = dnorm, n = 101, args = list(mean = estimators$mu[1], sd = estimators$sigma[1])) +
  theme_bw() + ylab("Count")

#Convet to ggplot https://stackoverflow.com/questions/51043753/overlay-histogram-with-empirical-density-and-dnorm-function

n=nrow(y)
binwidth = 0.5
ggplot(y, aes(x=Difference))+geom_histogram(binwidth = binwidth) +
  stat_function(fun = function(x)
  {dnorm(x, mean = estimators$mu[1], sd = estimators$sigma[1]) *n*binwidth},
  aes(colour = "Distribution 1"), size = 1)+
  stat_function(fun = function(x)
  {dnorm(x, mean = estimators$mu[2], sd = estimators$sigma[2]) *n*binwidth},
  aes(colour = "Distribution 2"), size = 1)+
  theme_bw() + ylab("Count")

#http://tinyheero.github.io/2015/10/13/mixture-model.html
plot_mix_comps <- function(x, mu, sigma, lam) {
  lam * dnorm(x, mu, sigma)
}

mixmdl<-model2
estimators<-as.data.frame(mixmdl[c("lambda", "mu", "sigma")])

ggplot(y, aes(x=Difference)) +
  geom_histogram(aes(x=Difference, ..density..), binwidth = 0.25) +
  stat_function(geom = "line", fun = plot_mix_comps,
             args = list(mixmdl$mu[1], mixmdl$sigma[1], lam = mixmdl$lambda[1]),
             colour = "red", lwd = 1.5) +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(mixmdl$mu[2], mixmdl$sigma[2], lam = mixmdl$lambda[2]),
                colour = "blue", lwd = 1.5) +
 # stat_function(geom = "line", fun = plot_mix_comps,
  #              args = list(mixmdl$mu[3], mixmdl$sigma[3], lam = mixmdl$lambda[3]),
   #             colour = "yellow", lwd = 1.5) +
  theme_bw()