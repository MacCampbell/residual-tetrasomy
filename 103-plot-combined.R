# Plot lastz output of multiple species on one plot
library(tidyverse)

#Loading premade dfs
load("./outputs/102/salmo-salar.rda")
load("./outputs/102/o-mykiss.rda")

#Adding sample size variable (already grouped by comparison)
`salmo-salar` <- `salmo-salar` %>% add_count(Comparison)
`o-mykiss`    <- `o-mykiss`    %>% add_count(Comparison)

#Filter salmo-salar based on the comparisons we want
`salmo-salar` <- `salmo-salar` %>% filter(Comparison %in% c("ssa02-ssa12","ssa11-ssa26","ssa07-ssa17",
                                                           "ssa02-ssa05", "ssa03-ssa06","ssa01-ssa18",
                                                           "ssa16-ssa17","ssa04-ssa08"))
#Add another column for homologies that we can use to plot o-mykiss
`salmo-salar` <- `salmo-salar` %>% mutate(HomologousPair = Comparison)

`o-mykiss` <- `o-mykiss` %>% mutate(HomologousPair = ifelse(Comparison == "CM007947.1-CM007951.1","ssa02-ssa12",
                                              ifelse(Comparison == "CM007940.1-CM007960.1","ssa11-ssa26",
                                                    ifelse(Comparison == "CM007949.1-CM007955.1", "ssa07-ssa17",
                                                           ifelse(Comparison == "CM007936.1-CM007937.1", "ssa02-ssa05",
                                              ifelse(Comparison == "CM007946.1-CM007947.1","ssa03-ssa06",
                                                     ifelse(Comparison == "CM007935.1-CM007957.1", "ssa01-ssa18",
                                                            ifelse(Comparison =="CM007941.1-CM007952.1", "ssa16-ssa17",
                                                                   ifelse(Comparison=="CM007944.1-CM007953.1", "ssa04-ssa08",
                                                                          "NA")))))) )))
#Filtering the columns for plotting 

`salmo-salar` <- `salmo-salar` %>% select(HomologousPair, Similarity, Mean, n) %>% mutate(Species = "salmo-salar")
`o-mykiss`    <- `o-mykiss`  %>% select(HomologousPair, Similarity, Mean, n) %>% mutate(Species = "o-mykiss")

data<-bind_rows(`salmo-salar`, `o-mykiss`)

#Plot in a specified order
data$Species<- factor(data$Species, levels = c("salmo-salar", "o-mykiss"))
data$HomologousPair <- factor(data$HomologousPair, levels = c("ssa03-ssa06",
                                                             "ssa02-ssa05",
                                                             "ssa11-ssa26",
                                                            "ssa02-ssa12",
                                                            "ssa16-ssa17",
                                                            "ssa07-ssa17",
                                                            "ssa04-ssa08",
                                                            "ssa01-ssa18"))



ggplot(data)+geom_boxplot(aes(x=HomologousPair, y=Similarity))+
  theme_classic()+
  ylab("Perecent Similarity")+
  theme(axis.text.x= element_text(angle=45,hjust=1))+
  facet_grid(.~Species)

ggsave("./outputs/103/combined-plot.pdf", width=11/2, height=8.5/2)
