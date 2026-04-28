## This script reads in range limits census data and analyzes herbivory

# remove all objects and clear workspace
rm(list = ls(all=TRUE))

##tom's local wd
setwd("C:/Users/tm9/Dropbox/github/endo-range-limits")

##load packages
library(tidyverse)
library(googlesheets4)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
library(googledrive)
library(readxl)
library(bayesplot)
library(tidybayes)
library(posterior)
library(magrittr)

##read in endo range limits database from google drive
##read in data as text so that nothing gets coerced
file_id<-"11m_-0KYtLtp5J0EqGWZ8RzNM-JfsDB3B"
drive_download(as_id(file_id), path = "Data/endo_range_limits_experiment.xlsx", overwrite = TRUE)
##plot-level treatment assignments
plots<-read_excel("Data/endo_range_limits_experiment.xlsx",sheet="Initial Plot Data", col_types = "text")
##initial plant data (including source and endo status)
init23<-read_excel("Data/endo_range_limits_experiment.xlsx",sheet="Initial Plant Data", col_types = "text")
##2024 replant
init24<-read_excel("Data/endo_range_limits_experiment.xlsx",sheet="2024 Re-plant inital data", col_types = "text")
##2023-25 census
census23<-read_excel("Data/endo_range_limits_experiment.xlsx",sheet="2023 Data Census", col_types = "text")
census24<-read_excel("Data/endo_range_limits_experiment.xlsx",sheet="2024 Data Census", col_types = "text")
census25<-read_excel("Data/endo_range_limits_experiment.xlsx",sheet="2025 Data Census", col_types = "text")


# data prep ---------------------------------------------------------------


##there is one unfortunated underscored tag id: 731_2, which was the second use of 731
##give it a new unique tag not used by anything else
##it only appears once in the initial data, and once in the 23 census
init23[init23$Tag_ID=="731_2","Tag_ID"]<-"10001"
census23[census23$Tag_ID=="731_2","Tag_ID"]<-"10001"

##drop all the initial plant tags starting with D (these were dead on first census)
init23 %<>% filter(substr(Tag_ID,1,1)!="D") %>% 
  ##coerce remaining tags, plot, and endo to integer
  mutate(Tag_ID = as.integer(as.numeric(Tag_ID)),
         Plot = as.integer(as.numeric(Plot)),
         Endo = as.integer(as.numeric(Endo)))
##warning OK: there is 1 NA tag ID -- dead ELVI in LAF

plots %<>% mutate(Plot = as.integer(as.numeric(Plot)),
                  Herbivory = case_when(Herbivory == "0.0" ~ "Access",
                                        Herbivory == "1.0" ~ "Exclusion"))

##merge initial plant and plot treatment data
init_plot <- left_join(init23,plots %>% select(Site,Plot,Species,Herbivory),
                       by=c("Site","Species","Plot"))
##2024 is more complicated because of replant
init24 %<>% filter(substr(Tag_ID,1,1)!="D") %>% 
  ##coerce remaining tags, plot, and endo to integer
  mutate(Tag_ID = as.integer(as.numeric(Tag_ID)),
         Plot = as.integer(as.numeric(Plot)),
         Endo = as.integer(as.numeric(Endo)),
         Herbivory = case_when(Herbivory_treatment == "0.0" ~ "Access",
                               Herbivory_treatment == "1.0" ~ "Exclusion"))

init_combo<-bind_rows(init_plot,init24 %>% select(names(init_plot))) %>% 
  filter(!is.na(Tag_ID))

census23 %<>% 
  ##drop dead plants
  filter(substr(Tag_ID,1,1)!="D") %>% 
  mutate(Tag_ID = as.integer(as.numeric(Tag_ID)),
         Tiller = as.integer(as.numeric(Tiller_23)),
         Inflor = as.integer(as.numeric(Inf_23)),
         TillerHerb = as.integer(as.numeric(tiller_Herb)))
#warning OK: one row of NA data, not sure why
##are all the 2023 census tags in the initial plant data? -- YES
census23[which(!(census23$Tag_ID %in% init_plot$Tag_ID)),]

##merge 2023 censusu and initial plant-plot data
dat23<-left_join(census23,init_plot,by="Tag_ID")

census24 %<>% 
  ##drop dead plants
  filter(substr(Tag_ID,1,1)!="D") %>% 
  mutate(Tag_ID = as.integer(as.numeric(Tag_ID)),
         Tiller = as.integer(as.numeric(Tiller_24)),
         AttachedInf = as.integer(as.numeric(attachedInf_24)),
         BrokenInf = as.integer(as.numeric(brokenInf_24)),
         TillerHerb = as.integer(as.numeric(tiller_Herb)))
  #warning OK
##are all the 2024 census tags in either init23 or init24?
sum(census24$Tag_ID %in% init23$Tag_ID)+sum(census24$Tag_ID %in% init24$Tag_ID)
length(census24$Tag_ID)
##it looks like there were 3 tagIDs that matched into both init23 and init24:
census24[(census24$Tag_ID %in% init23$Tag_ID) & (census24$Tag_ID %in% init24$Tag_ID),]
##tags 2097, 2098, and 2099 -- these are actually the same data in both places, 
##because they are the sole 3 survivors at SON

dat24<-left_join(census24,init_combo,by="Tag_ID",multiple="first")
##multiple=first deals with two init_combo rows for tag 2097

census25 %<>% 
  ##drop dead plants
  filter(substr(Tag_ID,1,1)!="D") %>% 
  mutate(Tag_ID = as.integer(as.numeric(Tag_ID)),
         Tiller = as.integer(as.numeric(Tiller_25)),
         AttachedInf = as.integer(as.numeric(attachedInf_25)),
         BrokenInf = as.integer(as.numeric(brokenInf_25)),
         TillerHerb = as.integer(as.numeric(tiller_Herb)))

dat25<-left_join(census25,init_combo,by=c("Site","Species","Tag_ID"))

##combine years
bind_rows(dat23 %>% 
            mutate(Year=2023,
                   AttachedInf=NA,
                   BrokenInf=NA) %>% 
            select(Year,Site,Species,Plot,Herbivory,Tag_ID,Population,Endo,Tiller,AttachedInf,BrokenInf,TillerHerb),
          dat24 %>% 
            mutate(Year=2024) %>% 
            select(Year,Site,Species,Plot,Herbivory,Tag_ID,Population,Endo,Tiller,AttachedInf,BrokenInf,TillerHerb),
          dat25 %>% 
            mutate(Year=2025) %>% 
            select(Year,Site,Species,Plot,Herbivory,Tag_ID,Population,Endo,Tiller,AttachedInf,BrokenInf,TillerHerb)
) %>% 
  mutate(Site=factor(Site,levels=c("SON","KER","BFL","BAS","COL","HUN","LAF")))->herb_dat_tidy

#write.csv(herb_dat_tidy,"Data/herbivory_tidy_data.csv")


# analysis ----------------------------------------------------------------
herb_dat_tidy<-read.csv("Data/herbivory_tidy_data.csv")

#quick visual of proportion of damaged tillers by endo
herb_dat_tidy %>% 
  mutate(TillerHerb=ifelse(TillerHerb>Tiller,Tiller,TillerHerb),
         HerbProp=TillerHerb/Tiller) %>% 
  filter(Herbivory=="Access") %>% 
  ggplot()+
  geom_boxplot(aes(x=Site,y=HerbProp,fill=as.factor(Endo)))+
  facet_grid(cols=vars(Species))

#quick visual of proportion of damaged tillers by fency
herb_dat_tidy %>% 
  mutate(TillerHerb=ifelse(TillerHerb>Tiller,Tiller,TillerHerb),
         HerbProp=TillerHerb/Tiller) %>% 
  ggplot()+
  geom_boxplot(aes(x=Site,y=HerbProp,fill=as.factor(Herbivory)))+
  facet_grid(cols=vars(Species))

##combine endo and herbivory
herb_dat_tidy %>% 
  filter(Tiller>0 & !is.na(Species) & !is.na(Herbivory)) %>% 
  mutate(TillerHerb=ifelse(TillerHerb>Tiller,Tiller,TillerHerb),
         HerbProp=TillerHerb/Tiller) %>% 
  ggplot()+
  geom_boxplot(aes(x=Site,y=HerbProp,fill=as.factor(Endo)))+
  facet_grid(Herbivory~Species)

##inf damage
herb_dat_tidy %>% 
  filter(Year>2023) %>% 
  filter(Species!="AGHY") %>% 
  mutate(TotInf=BrokenInf+AttachedInf) %>% 
  filter(TotInf>0) %>% 
  mutate(InfDam = BrokenInf/TotInf) %>% 
  ggplot()+
  geom_boxplot(aes(x=Site,y=InfDam,fill=as.factor(Endo)))+
  facet_grid(Herbivory~Species)

##bundle data for stan

herb_dat<-herb_dat_tidy %>% 
  ##there are a few cases where damaged leaves were counted instead of tillers
  mutate(TillerHerb=ifelse(TillerHerb>Tiller,Tiller,TillerHerb)) %>% 
  ##keep only living plants
  filter(Tiller>0) %>% 
  select(Species,Site,Plot,Endo,Herbivory,Year,Tiller,TillerHerb) %>% 
  drop_na() %>% 
  ##convert species, sites, and years to integer
  mutate(species=as.integer(as.factor(Species)),
         fence=as.integer(as.factor(Herbivory))-1,
         site=as.integer(as.factor(Site)),
         year=Year-2022,
         unique_plot=interaction(Species,Site,Plot),
         plot_int=as.integer(as.factor(unique_plot))) 

stan_dat<-list(n=nrow(herb_dat),
               n_species=length(unique(herb_dat$species)),
               n_endo=length(unique(herb_dat$Endo)),
               n_fence=length(unique(herb_dat$fence)),
               n_site=length(unique(herb_dat$site)),
               n_year=length(unique(herb_dat$year)),
               n_plot=max(herb_dat$plot_int),
               y_tillers=herb_dat$Tiller,
               y_damaged=herb_dat$TillerHerb,
               species=herb_dat$species,
               endo=herb_dat$Endo,
               fence=herb_dat$fence,
               site=herb_dat$site,
               year=herb_dat$year,
               plot=herb_dat$plot_int)

##compile model
herb_model<-stan_model("stan/herbivory.stan")

##sample model
herb_fit<-sampling(herb_model,data=stan_dat,chains=3,iter=5000,
                   pars=c("beta0","beta_endo","beta_fence",
                          "eps_site","eps_year","eps_plot",
                          "sigma_site","sigma_year","sigma_plot",
                          "y_rep"),include=T)

##a few trace plots...
mcmc_trace(herb_fit,pars=c("beta0[1]","beta0[2]","beta0[3]"))

##posterior predictive check
y_rep<-rstan::extract(herb_fit,pars="y_rep")
ppc_dens_overlay(stan_dat$y_damaged,y_rep$y_rep[1:100,])

species_names <- c("AGHY", "ELVI", "POAU")
beta_draws <- herb_fit %>%
  spread_draws(
    beta0[species],
    beta_endo[species],
    beta_fence[species]
  ) %>%
  mutate(
    species_name = species_names[species]
  ) %>%
  pivot_longer(
    cols = c(beta0, beta_endo, beta_fence),
    names_to = "parameter",
    values_to = "estimate"
  ) %>%
  mutate(
    parameter = recode(
      parameter,
      beta0 = "Intercept",
      beta_endo = "Endophyte effect",
      beta_fence = "Fence effect"
    ),
    parameter = factor(parameter, levels = c("Intercept", "Endophyte effect", "Fence effect"))
  )

beta_summary <- beta_draws %>%
  group_by(parameter, species, species_name) %>%
  median_qi(estimate, .width = 0.9) %>%
  ungroup()

ggplot(beta_summary,
       aes(x = estimate, y = fct_reorder(species_name, estimate),
           color = species_name)) +
  geom_vline(xintercept = 0, linewidth = 0.4, linetype = "dashed") +
  geom_linerange(
    aes(xmin = .lower, xmax = .upper, linewidth = factor(.width)),
    position = position_dodge(width = 0.5)
  ) +
  geom_point(size = 4) +
  facet_wrap(~ parameter, scales = "free_x") +
  scale_linewidth_manual(
    values = c(`0.9` = 1.4),
    name = "Credible interval"
  ) +
  labs(
    x = "Posterior estimate",
    y = NULL,
    color = "Species"
  ) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )


## site effects
site_names <- c("SON","KER","BFL","BAS","COL","HUN","LAF")  
site_draws <- herb_fit %>%
  spread_draws(eps_site[site]) %>%
  mutate(
    site_name = site_names[site]
  )

site_summary <- site_draws %>%
  group_by(site, site_name) %>%
  median_qi(eps_site, .width = 0.95) %>%
  ungroup()

ggplot(site_summary,
       aes(x = eps_site,
           y = fct_reorder(site_name, eps_site))) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4) +
  geom_linerange(aes(xmin = .lower, xmax = .upper), linewidth = 1) +
  geom_point(size = 2) +
  labs(
    x = "Site random effect",
    y = NULL
  ) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank()
  )
