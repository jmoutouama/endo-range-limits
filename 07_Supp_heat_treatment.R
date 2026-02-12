## Purpose: quantify effectiveness of heat treatment for endophyte removal
## Author: Tom Miller

library(tidyverse)
library(magrittr)
library(readxl)
library(xtable)
setwd("C:/Users/tm9/Dropbox/github/endo-range-limits")

##read in original plant data (we should include these as csv's in the repo)
datini <- read.csv("https://www.dropbox.com/scl/fi/b93bvocqltadc36xirak2/Initialdata.csv?rlkey=8hd3z4th35lqvtfvam83kb972&dl=1", stringsAsFactors = F)
##combine SHS and HUNT pops
datini %<>% mutate(Population = recode_factor(Population, SHS = "HUNT"))

##read in greenhouse inventory (all the plants we reared -- those in the experiment and many others)
greenhouse<-read_excel("Data/2022-23_Greenhouse_Inventory.xlsx",sheet="Range Limit Greenhouse")
str(greenhouse)
#why is real_endo a character vector?
unique(greenhouse$real_endo) ##because the excel function populates with "NA"
##combine SHS and HUNT pops
greenhouse %<>% mutate(POPULATION = recode_factor(POPULATION, SHS = "HUNT"))

# source populations ------------------------------------------------------

##read in source info
sources<-read.csv("Data/endo_range_limits_experiment- Source Population.csv")
## drop SHS since it's a duplicate row with HUNT
sources %<>% filter(Population!="SHS")
## customize notes for ms table
sources[1,c("Population","Notes")]
sources$msNotes[1]<-"Private property, San Jacinto County, TX"
sources[2,c("Population","Notes")]
sources$msNotes[2]<-"Roadside collection, Columbus, TX"
sources[3,c("Population","Notes")]
sources$msNotes[3]<-"Stengl Lost Pines Field Station, Bastrop, TX"
sources[4,c("Population","Notes")]
sources$msNotes[4]<-"Roadside population, DeRidder, LA"
sources[5,c("Population","Notes")]
sources$msNotes[5]<-"Walter P. Jacobs Nature Center, Shreveport, LA"
sources[6,c("Population","Notes")]
sources$msNotes[6]<-"Kistachie National Forest, LA"
sources[7,c("Population","Notes")]
sources$msNotes[7]<-"Stephen F. Austin Experimental Forest, TX"
sources[8,c("Population","Notes")]
sources$msNotes[8]<-"Sam Houston National Forest, TX"
sources[9,c("Population","Notes")]
sources$msNotes[9]<-"Pineywoods Environmental Research Laboratory, Huntsville, TX"
sources[10,c("Population","Notes")]
sources$msNotes[10]<-"Llano River Field Station, Junction, TX"
sources[11,c("Population","Notes")]
sources$msNotes[11]<-"Palmetto State Park, TX"
sources[12,c("Population","Notes")]
sources$msNotes[12]<-"Jean Laffitte Park, LA"
##update species names
sources %<>% 
  mutate(Species2 = fct_recode(factor(Species),
                          "A. hyemalis" = "AGHY",
                          "E. virginicus" = "ELVI",
                          "P. autumnalis" = "POAU"))

##how many genotypes from each soure population were used?
source_tab1<-datini %>% 
  group_by(Species,Population) %>% 
  summarise(number_plants=n(),number_genotypes=length(unique(GreenhouseID)))

datini %>% 
  group_by(Species) %>% 
  summarise(number_plants=n(),number_genotypes=length(unique(GreenhouseID)))

##what was the average (and range of) clonal replication per population?
source_tab2<-datini %>% 
  group_by(Species,Population,GreenhouseID) %>% 
  summarise(clones=n()) %>% 
  group_by(Species,Population) %>% 
  summarise(mean_clones=mean(clones),min_clones=min(clones),max_clones=max(clones))

## how many sources per species?
sources %>% group_by(Species) %>% summarise(n())##but POAU WB not used
##what was the average number of genotypes per source?
source_tab1 %>% summarise(mean(number_genotypes))
##and the average number of clones per genotype?
source_tab2 %>% summarise(mean(mean_clones))
  
## combine into table for manuscript, print out as .tex file, import into manuscript
supp_table_sources<-left_join(sources %>% filter(Population!="WB"),#no WB in the experiment
          left_join(source_tab1,source_tab2,by=c("Species","Population")),
          by=c("Species","Population")) %>% 
  select(Species2,Population,Lat,Long,msNotes,number_plants,number_genotypes,mean_clones) %>% 
  rename(
    Species = Species2,
    Population = Population,
    Latitude = Lat,
    Longitude = Long,
    Notes = msNotes,
    `Total plants` = number_plants,
    `Unique genotypes` = number_genotypes,
    `Mean clones per genotype` = mean_clones) %>% 
  mutate(Species = paste0("\\textit{", Species, "}"))
## format latex table
tab<-xtable(supp_table_sources,digits = c(0,0,0,6,6,0,0,0,2))
align(tab)<-"lllrrp{0.2\\textwidth}p{2cm}p{2cm}p{2cm}" ## wrap notes
print(tab,file = "Manuscript/source_table.tex",include.rownames=F,
  sanitize.text.function=identity,floating=F)


# endophyte scores --------------------------------------------------------

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
  select(SPECIES,POPULATION,`INDIVIDUAL ID`,ENDOPHYTE,
         Leaf_peel_liberal,Leaf_peel_conservative,seed_score_liberal,
         seed_score_conservative,agrinostics_score_liberal,agrinostics_score_conservative) %>% 
  group_by(SPECIES,POPULATION,`INDIVIDUAL ID`,ENDOPHYTE) %>% 
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
  rowwise() %>%
  mutate(mean_liberal=mean(c_across(c(Leaf_peel_liberal,seed_score_liberal,agrinostics_score_liberal)),na.rm=T),
         mean_conservative=mean(c_across(c(Leaf_peel_conservative,seed_score_conservative,agrinostics_score_conservative)),na.rm=T)) %>% 
  ##find mismatrches (mean is not zero or one)
  mutate(mismatch_liberal=mean_liberal!=0 & mean_liberal!=1,
         mismatch_conservative=mean_conservative!=0 & mean_conservative!=1) %>% 
  group_by() %>% 
  summarise(n(),
            mean(mismatch_liberal==0),
            mean(mismatch_conservative==0))

## calculate endo prevalence conditional on expected status
scored_greenhouse_genotypes %>% 
  ## I will use max() here so that in the few cases where scores differed, 
  ## we will assume the positive score is correct
  mutate(mean_liberal=max(c_across(c(Leaf_peel_liberal,seed_score_liberal,agrinostics_score_liberal)),na.rm=T),
         mean_conservative=max(c_across(c(Leaf_peel_conservative,seed_score_conservative,agrinostics_score_conservative)),na.rm=T)) %>%
  group_by(SPECIES,POPULATION,ENDOPHYTE) %>% 
  summarise(count=n(),
            prevalence_liberal=mean(mean_liberal),
            prevalence_conservative=mean(mean_conservative))->endo_scores

##mean prevalence by species
endo_scores %>% group_by(SPECIES,ENDOPHYTE) %>% summarise(mean(prevalence_liberal))

treatment_table<-endo_scores %>% 
  ungroup() %>%
  mutate(Species = fct_recode(factor(SPECIES),
                               "A. hyemalis" = "AGHY",
                               "E. virginicus" = "ELVI",
                               "P. autumnalis" = "POAU"),
         Treatment = fct_recode(factor(ENDOPHYTE),"Heat"="E+/- Heat","Control"="E+ Con")) %>% 
  rename(Population=POPULATION,
        `\\%E+`=prevalence_liberal,
         N=count) %>% 
  select(Species,Population,Treatment,N,`\\%E+`)%>% 
  mutate(Species = paste0("\\textit{", Species, "}"))

print(xtable(treatment_table),file = "Manuscript/treatment_table.tex",include.rownames=F,
      sanitize.text.function=identity,floating = FALSE,)

##what was the mean E+ prevalence in controls?
treatment_table %>% filter(Treatment=="Control") %>% summarise(mean(`\\%E+`))