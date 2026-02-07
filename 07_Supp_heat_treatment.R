## Purpose: quantify effectiveness of heat treatment for endophyte removal
## Author: Tom Miller

library(tidyverse)
library(readxl)
setwd("C:/Users/tm9/Dropbox/github/endo-range-limits")

##read in original plant data (we should include these as csv's in the repo)
datini <- read.csv("https://www.dropbox.com/scl/fi/b93bvocqltadc36xirak2/Initialdata.csv?rlkey=8hd3z4th35lqvtfvam83kb972&dl=1", stringsAsFactors = F)

##read in greenhouse inventory (all the plants we reared -- those in the experiment and many others)
greenhouse<-read_excel("Data/2022-23_Greenhouse_Inventory.xlsx",sheet="Range Limit Greenhouse")
str(greenhouse)

#why is real_endo a character vector?
unique(greenhouse$real_endo) ##because the excel function populates with "NA"

##how many genotypes from each soure population were used?
source_tab1<-datini %>% 
  group_by(Species,Population) %>% 
  summarise(number_plants=n(),number_genotypes=length(unique(GreenhouseID)))

##what was the average (and range of) clonal replication per population?
source_tab2<-datini %>% 
  group_by(Species,Population,GreenhouseID) %>% 
  summarise(clones=n()) %>% 
  group_by(Species,Population) %>% 
  summarise(mean_clones=mean(clones),min_clones=min(clones),max_clones=max(clones))

## combine into table for manuscript
left_join(source_tab1,source_tab2,by=c("Species","Population"))

##filter the source pops in the greenhouse inventory 
##to only include those used in the experiment
greenhouse %>% 
  filter(POPULATION %in% datini$Population) %>% 
  #filter only plants with scores
  filter(!is.na(Leaf_peel_liberal)|!is.na(Leaf_peel_conservative)|
           !is.na(seed_score_liberal)|!is.na(seed_score_conservative)|
           !is.na(agrinostics_score_liberal)|!is.na(agrinostics_score_conservative)) %>% 
  #derive "real endophyte" score as the average of all possible scores
  rowwise() %>%
  mutate(endo_score=mean(c_across(Leaf_peel_liberal:agrinostics_score_conservative),na.rm=T))->scored_greenhouse

#confirm that clones from the same genotype always share the same endo status
scored_greenhouse %>% 
  filter(!is.na(`Clone ID`)) %>% ##keep only genotypes that were cloned
  group_by(SPECIES,`INDIVIDUAL ID`) %>% 
  summarise(genotype_endo=mean(endo_score)) %>% 
  group_by(SPECIES,genotype_endo) %>% 
  summarise(n())
#there were 4 instances of scores diverging between clones of the same genotype
#3 in AGHY, 1 in POAU -- have a closer look at these instances
scored_greenhouse %>% 
  filter(!is.na(`Clone ID`)) %>% 
  group_by(SPECIES,`INDIVIDUAL ID`) %>% 
  summarise(genotype_endo=mean(endo_score)) %>% 
  filter(genotype_endo>0 & genotype_endo<1) %>% 
  select(SPECIES,`INDIVIDUAL ID`)->problem_genotypes
greenhouse %>% filter(`INDIVIDUAL ID` %in% problem_genotypes$`INDIVIDUAL ID`) %>% View
## these are false alarms -- all clones have the same genotype, they are just split score

## reduce database to genotypes, collapse clones
scored_greenhouse %>% 
  select(SPECIES,POPULATION,`INDIVIDUAL ID`,
         Leaf_peel_liberal,Leaf_peel_conservative,seed_score_liberal,
         seed_score_conservative,agrinostics_score_liberal,agrinostics_score_conservative) %>% 
  group_by(SPECIES,POPULATION,`INDIVIDUAL ID`) %>% 
  distinct() %>% 
  ##add indicator of which score type they had
  mutate(leaf_scored=!is.na(Leaf_peel_liberal) | !is.na(Leaf_peel_conservative),
         seed_scored=!is.na(seed_score_liberal) | !is.na(seed_score_conservative),
         tiller_scored=!is.na(agrinostics_score_liberal) | !is.na(agrinostics_score_conservative))->scored_greenhouse_genotypes

## I am sometimes averaging over different methods (eg leaf and agrinostics)
## what is their percentage agreement?
scored_greenhouse_genotypes %>% 
##subset to scores from multiple methods
  filter((leaf_scored+seed_scored+tiller_scored)>1) %>% 
  mutate(mean_endo=)

## calculate endo prevalence

  
  

