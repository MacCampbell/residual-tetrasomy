#! /usr/local/bin/Rscript
# Usage: ./104-plot-color salmo-salar
# Or $ cat taxon-list.txt | while read line; do ./104-plot-color.R $line; done;
# Plot lastz output with error bars
rm(list=ls())
library(tidyverse)
library(matrixStats)
library(viridis, data.table)
library(readxl)



#species <- read_excel("C:\\Users\\Wlarson\\Dropbox\\multiplot\\genomes.xlsx", col_names = FALSE)
df_total <- NULL
species=c("O.kis","O.myk","O.tsh","S.alp","S.sal","T.thy")
# O.kis
# O.myk
# O.tsh
# S.alp
# S.sal
# T.thy

#x <-1

for (i in species) {
#args <- commandArgs(trailingOnly = TRUE)  
#args<-c(i[x])
protos<-as_tibble(read.table(paste("C:\\Users\\Wlarson\\Dropbox\\multiplot\\protokaryo\\",i, "-protokaryotype.txt", sep=""))) %>%
  rename(Protokaryotype = V1, Comparison = V2)
protos$Protokaryotype<-as.character(protos$Protokaryotype)
protos$Comparison<-as.character(protos$Comparison)


data<-as_tibble(read.table(paste("C:\\Users\\Wlarson\\Dropbox\\multiplot\\lastz\\",i,"-101-lastz.txt", sep="")))
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
print(i)
data$species <-i

df_total <- rbind(df_total,data)

}

#perm test

coho=df_total[which(df_total$species=="O.kis"),]
mykiss=df_total[which(df_total$species=="O.myk"),]
chin=df_total[which(df_total$species=="O.tsh"),]
char=df_total[which(df_total$species=="S.alp"),]
atlantic=df_total[which(df_total$species=="S.sal"),]
grayling=df_total[which(df_total$species=="T.thy"),]
PKs=unique(coho$Protokaryotype)

head(coho)
#test to see if each PK has more sim than expected
#plan, get all aligns for each pk, get mean sim, randomly draw same number of aligns out of full dataset,
#get mean sim, compare mean similarity, count number of times test mean diff is more than actual observed diff between 
#mean of current chr and mean of whole dataset
#https://www.uvm.edu/~dhowell/StatPages/Randomization%20Tests/RandomizationTestsOverview.html
#it seemed too extreme when just sampling random alignments, everything was sig, i think ill sample random pks instead
#balls when sampling full PKs only the highest duplicates were outliers
#new plan, separate magic 8 and others, see if each chr is sig diff from magic 8 or from rest of chrs, go back to pulling
#by chr


perm_test=function(input_df=coho, pk="2",num_perm=1000)
{
  magic_8_pks=c(2,6,9,11,20,22,23,25)
  not_magic_pks=c(1,3,4,5,7,8,10,12,13,14,15,16,17,18,19,21,24)
  magic_8_dat=subset(input_df, subset = Protokaryotype %in% magic_8_pks)
  not_magic_dat=subset(input_df, subset = Protokaryotype %in% not_magic_pks)
  overall_mean=mean(input_df$Similarity)
  current_pk_dat=input_df[which(input_df$Protokaryotype==pk),]
  current_pk_mean=mean(current_pk_dat$Similarity)
  observed_diff_in_means_total=abs(overall_mean-current_pk_mean)
  
  naligns=nrow(current_pk_dat[,1])
  magic_8_counter=0
  not_magic_counter=0
  for(i in 1:num_perm)
  {
    #calc sim means
    rand_magic8_mean=mean(sample(magic_8_dat$Similarity, naligns, replace=FALSE))
    rand_not_magic_mean=mean(sample(not_magic_dat$Similarity, naligns, replace=FALSE))
    #calc diffs in means
    sim_mean_diff_magic_8=abs(rand_magic8_mean-current_pk_mean)
    if(sim_mean_diff_magic_8 > observed_diff_in_means_total){magic_8_counter=magic_8_counter+1}
    #not magic 8
    sim_mean_diff_not_magic=abs(rand_not_magic_mean-current_pk_mean)
    if(sim_mean_diff_not_magic > observed_diff_in_means_total){not_magic_counter=not_magic_counter+1}
  }  
  p_value_magic_8=magic_8_counter/num_perm
  p_value_not_magic=not_magic_counter/num_perm
  print(pk)
  print(c("p val magic 8=",p_value_magic_8))
  print(c("p val not magic=",p_value_not_magic))
  return(cbind(pk, p_value_magic_8,p_value_not_magic))
}

coho_sig_diffs=NULL
for(i in sort(as.numeric(PKs)))
{
  p_val_per_pk=perm_test(pk=as.character(i))
  #sig="no"
  #if(p_val_per_pk<0.05){sig="yes"}
  coho_sig_diffs=rbind(coho_sig_diffs, p_val_per_pk)
}
write.csv(coho_sig_diffs, row.names=FALSE, file="C:\\Users\\Wlarson\\Dropbox\\multiplot\\coho_sig_diffs.csv")

mykiss_sig_diffs=NULL
for(i in sort(as.numeric(PKs)))
{
  p_val_per_pk=perm_test(input_df=mykiss, pk=as.character(i))
  #sig="no"
  #if(p_val_per_pk<0.05){sig="yes"}
  mykiss_sig_diffs=rbind(mykiss_sig_diffs, p_val_per_pk)
}
write.csv(mykiss_sig_diffs, row.names=FALSE, file="C:\\Users\\Wlarson\\Dropbox\\multiplot\\mykiss_sig_diffs.csv")


chin_sig_diffs=NULL
for(i in sort(as.numeric(PKs)))
{
  p_val_per_pk=perm_test(input_df=chin, pk=as.character(i))
  #sig="no"
  #if(p_val_per_pk<0.05){sig="yes"}
  chin_sig_diffs=rbind(chin_sig_diffs, p_val_per_pk)
}
write.csv(chin_sig_diffs, row.names=FALSE, file="C:\\Users\\Wlarson\\Dropbox\\multiplot\\chin_sig_diffs.csv")


char_sig_diffs=NULL
for(i in sort(as.numeric(PKs)))
{
  p_val_per_pk=perm_test(input_df=char, pk=as.character(i))
  #sig="no"
  #if(p_val_per_pk<0.05){sig="yes"}
  char_sig_diffs=rbind(char_sig_diffs, p_val_per_pk)
}
write.csv(char_sig_diffs, row.names=FALSE, file="C:\\Users\\Wlarson\\Dropbox\\multiplot\\char_sig_diffs.csv")


atlantic_sig_diffs=NULL
for(i in sort(as.numeric(PKs)))
{
  p_val_per_pk=perm_test(input_df=atlantic, pk=as.character(i))
  #sig="no"
  #if(p_val_per_pk<0.05){sig="yes"}
  atlantic_sig_diffs=rbind(atlantic_sig_diffs, p_val_per_pk)
}
write.csv(atlantic_sig_diffs, row.names=FALSE, file="C:\\Users\\Wlarson\\Dropbox\\multiplot\\atlantic_sig_diffs.csv")


grayling_sig_diffs=NULL
for(i in sort(as.numeric(PKs)))
{
  p_val_per_pk=perm_test(input_df=grayling, pk=as.character(i))
  #sig="no"
  #if(p_val_per_pk<0.05){sig="yes"}
  grayling_sig_diffs=rbind(grayling_sig_diffs, p_val_per_pk)
}
write.csv(grayling_sig_diffs, row.names=FALSE, file="C:\\Users\\Wlarson\\Dropbox\\multiplot\\grayling_sig_diffs.csv")

