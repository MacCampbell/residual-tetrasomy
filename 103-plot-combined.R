# Plot lastz output of multiple species on one plot
library(tidyverse)
library(matrixStats)
#Loading premade dfs
load("./outputs/102/salmo-salar.rda")
load("./outputs/102/o-mykiss.rda")
load("./outputs/102/t-thymallus.rda")

#Adding sample size variable (already grouped by comparison)
`salmo-salar` <- `salmo-salar` %>% add_count(Comparison)
`o-mykiss`    <- `o-mykiss`    %>% add_count(Comparison)
`t-thymallus` <- `t-thymallus` %>% add_count(Comparison)


#Adding a Length based on splitting the coverage field:
#Fraction of the entire input sequence (target or query, whichever is shorter) that is covered by the alignment block (see Coverage). This is written as two fields. The first field is a fraction, written as <n>/<d>. The second field contains the same value, computed as a percentage.
`salmo-salar` <- `salmo-salar` %>% separate(V14, sep="[/]", into = c("AlignmentLength","Length"))
`o-mykiss`    <- `o-mykiss` %>%    separate(V14, sep="[/]", into = c("AlignmentLength","Length"))
`t-thymallus`<- `t-thymallus` %>%  separate(V14, sep="[/]", into = c("AlignmentLength","Length"))

`salmo-salar`$AlignmentLength <-as.numeric(`salmo-salar`$AlignmentLength)
`o-mykiss`$AlignmentLength    <-as.numeric(`o-mykiss`$AlignmentLength)
`t-thymallus`$AlignmentLength <-as.numeric(`t-thymallus`$AlignmentLength)

#Adding Median variable (already grouped by comparison)
`salmo-salar` <- `salmo-salar` %>% mutate(Median = weightedMedian(Similarity, w=AlignmentLength))
`o-mykiss`    <- `o-mykiss`    %>% mutate(Median = weightedMedian(Similarity, w=AlignmentLength))
`t-thymallus` <- `t-thymallus` %>% mutate(Median = weightedMedian(Similarity, w=AlignmentLength))


#Filter salmo-salar based on the comparisons we want
`salmo-salar` <- `salmo-salar` %>% filter(Comparison %in% c("ssa02-ssa12","ssa11-ssa26","ssa07-ssa17",
                                                           "ssa02-ssa05", "ssa03-ssa06","ssa01-ssa18",
                                                           "ssa16-ssa17","ssa04-ssa08"))
#Add another column for homologies that we can use to plot o-mykiss
`salmo-salar` <- `salmo-salar` %>% mutate(HomologousPair = Comparison)

`o-mykiss` <- `o-mykiss` %>% mutate(HomologousPair = ifelse(Comparison == "CM007947.1-CM007951.1", "ssa02-ssa12",
                                                     ifelse(Comparison == "CM007940.1-CM007960.1", "ssa11-ssa26",
                                                     ifelse(Comparison == "CM007949.1-CM007955.1", "ssa07-ssa17",
                                                     ifelse(Comparison == "CM007936.1-CM007937.1", "ssa02-ssa05",
                                                     ifelse(Comparison == "CM007946.1-CM007947.1", "ssa03-ssa06",
                                                     ifelse(Comparison == "CM007935.1-CM007957.1", "ssa01-ssa18",
                                                     ifelse(Comparison == "CM007941.1-CM007952.1", "ssa16-ssa17",
                                                     ifelse(Comparison == "CM007944.1-CM007953.1", "ssa04-ssa08",
                                                                          "NA")))))))))

`t-thymallus` <- `t-thymallus` %>% mutate(HomologousPair = ifelse(Comparison == "CM015039.1-CM015040.1", "ssa02-ssa12",
                                                           ifelse(Comparison == "CM015025.1-CM015026.1", "ssa11-ssa26",
                                                           ifelse(Comparison == "CM015019.1-CM015020.1", "ssa07-ssa17",
                                                           ifelse(Comparison == "CM015013.1-CM015014.1", "ssa02-ssa05",
                                                           ifelse(Comparison == "CM014992.1-CM014993.1", "ssa03-ssa06",
                                                           ifelse(Comparison == "CM015033.1-CM015034.1", "ssa01-ssa18",
                                                           ifelse(Comparison == "CM015017.1-CM015018.1", "ssa16-ssa17",
                                                           ifelse(Comparison == "CM015023.1-CM015024.1", "ssa04-ssa08",
                                                                                "NA")))))) )))


#Filtering the columns for plotting 

`salmo-salar` <- `salmo-salar` %>% select(HomologousPair, Similarity, Mean, Median, n, AlignmentLength) %>% mutate(Species = "salmo-salar")
`o-mykiss`    <- `o-mykiss`  %>% select(HomologousPair, Similarity, Mean, Median, n, AlignmentLength) %>% mutate(Species = "o-mykiss")
`t-thymallus` <- `t-thymallus`  %>% select(HomologousPair, Similarity, Mean, Median, n, AlignmentLength) %>% mutate(Species = "t-thymallus")

data<-bind_rows(`salmo-salar`, `o-mykiss`, `t-thymallus`)


#Reorder by median
order<-`salmo-salar` %>% select(HomologousPair, Median) %>% unique() %>% arrange(desc(Median))

#Plot in a specified order
data$Species<- factor(data$Species, levels = c("salmo-salar", "o-mykiss",
                                               "t-thymallus"))
data$HomologousPair <- factor(data$HomologousPair, levels = order$HomologousPair)
#                                                      levels = c("ssa03-ssa06",
#                                                             "ssa02-ssa05",
#                                                             "ssa11-ssa26",
#                                                            "ssa02-ssa12",
#                                                            "ssa16-ssa17",
#                                                            "ssa07-ssa17",
#                                                            "ssa04-ssa08",
#                                                            "ssa01-ssa18"))

#Get Medians for plotting
totOrder <- data %>% group_by(Species, Comparison, HomologousPair, Median) %>% summarize()
  
#Note: To change x-axis labels in differnt facets, you can use "scales="free_x" in facet_grid
# e.g. https://github.com/duttashi/visualizer/issues/47
ggplot(data)+geom_boxplot(aes(x=HomologousPair, y=Similarity, weight=AlignmentLength), 
                          outlier.size=0.1, outlier.alpha=0.5,
                          outlier.color="gray50", outlier.shape=15)+
  ylab("Perecent Similarity")+
  xlab("Homologous Pair")+
  facet_wrap(.~Species, scales="free_x", nrow=2)+
  geom_text(data = totOrder, aes(label = round(Median,2),
                              x = HomologousPair, y = 100 + 2),
            size=2)+
  theme_bw()+
  theme(axis.text.x= element_text(angle=45,hjust=1))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggsave("./outputs/103/combined-plot.pdf", width=11/2, height=6)
