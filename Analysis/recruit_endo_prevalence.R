## This script fits models for recruit endo prevalence (from agrinostics scores)
library(tidyverse)
library(googlesheets4)
library(rstan)
library(googledrive)
library(readxl)
library(bayesplot)
library(tidybayes)
library(posterior)
setwd("C:/Users/tm9/Dropbox/github/ELVI-endophyte-density")

## read in agrinostics scores from lab drive
path1<-gs4_get("https://docs.google.com/spreadsheets/d/1vHibZHiy-vdkwIuLU1EAah9cBx4XqloSkyTFDXMtoYw/edit?gid=1882843135#gid=1882843135")
recruits<-read_sheet(path1,sheet="recruit_harvest")
## read in plot-level treatment data
## this is in the drive but it's an .xlsx file so I need to read it differently


file_id<-"11m_-0KYtLtp5J0EqGWZ8RzNM-JfsDB3B"
drive_download(as_id(file_id), path = "Data/endo_range_limits_experiment.xlsx", overwrite = TRUE)
plots <- readxl::read_excel("Data//endo_range_limits_experiment.xlsx",sheet="Initial Plot Data")
## plot data should have unique plot numbers that line up with recruit data
## also convert species, site, and fence treatments to positive integer for indexing
plots %>%
  ##drop sonora bc we have no recruit samples
  filter(Site!="SON") %>% 
  mutate(unique_plot = interaction(Species,Site,Plot),
         unique_plot_id = row_number(),
         species_int = as.integer(as.factor(Species)),
         site_int = as.integer(as.factor(Site)),
         fence_int = Herbivory+1)->plots

## proceeding here with endoScore_lib (leaf peels usually corresponded well where we have them)
## were there any NAs?
recruits[is.na(recruits$endoScore_lib),] %>% View
## drop NAs, and drop the bastrop natural population and the tagged census plants
recruits %<>% 
  filter(!is.na(endoScore_lib)) %>% 
  filter(individual=="NA") %>% 
  filter(subplot!="NATURAL_POPULATION") %>% 
  ## make a unique plot number
  mutate(unique_plot = interaction(species,site,plot)) %>% 
  ## make sure these plot numbers align with the plots df
  left_join(.,plots,by="unique_plot")

## make subplot data frame
recruits %>% 
  ## drop extras
  filter(subplot!="EXTRA") %>% 
  ##tally positives and samples by species, site, plot, subplot
  group_by(species,site,unique_plot,unique_plot_id,subplot) %>% 
  summarise(samples=n(),
            positives=sum(endoScore_lib))->subplot_recruits

## make a data frame for extras
recruits %>% 
  ## filter extras
  filter(subplot=="EXTRA") %>% 
  ##tally positives and samples by species, site, plot
  group_by(species,site,unique_plot,unique_plot_id) %>% 
  summarise(samples=n(),
            positives=sum(endoScore_lib))->extra_recruits

## tidy and bundle data
prevalence_data <- list(
  n_spp = max(plots$species_int),
  n_sites = max(plots$site_int), #6 sites, no Sonora data
  n_fence = max(plots$fence_int), #herbivory=0 is half fence
  n_plots = max(plots$unique_plot_id),
  initprev = plots$Endo_freq,
  species = plots$species_int,
  site = plots$site_int,
  fence = plots$fence_int,
  y1_n = nrow(subplot_recruits),
  y1_plot = subplot_recruits$unique_plot_id,
  y1_samples = subplot_recruits$samples,
  y1_positives = subplot_recruits$positives,
  y2_n = nrow(extra_recruits),
  y2_plot = extra_recruits$unique_plot_id,
  y2_samples = extra_recruits$samples,
  y2_positives = extra_recruits$positives  
)

## compile model
prevalence_model<-stan_model("Analysis/stan/recruit_prevalence.stan")
prevalence_samples<-sampling(prevalence_model, data = prevalence_data,
                             chains=3,iter=5000,thin=2,
                             cores = parallel::detectCores() - 1,
                       pars = c("beta0","beta1","rho1","rho2","p"),   
                       save_warmup=F)
#saveRDS(prevalence_samples,"Analysis/stan/prevalence_samples.rds")
prevalence_samples<-readRDS("Analysis/stan/prevalence_samples.rds")
## I got this code from CHatGPT, then modified
# Convert MCMC samples to tidy format
tidy_samples <- prevalence_samples %>%
  spread_draws(p[i])  # Extract all p[i] elements

# Compute posterior summaries
summary_p <- tidy_samples %>%
  filter(i %in% unique(recruits$unique_plot_id)) %>% #filter out plots that we have no data for
  group_by(i) %>%
  summarise(
    mean = mean(p),
    lower = quantile(p, 0.025),  # 2.5% credible interval
    upper = quantile(p, 0.975)   # 97.5% credible interval
  )

# Caterpillar plot
ggplot(summary_p, aes(x = mean, y = reorder(i, mean))) +  # Order by mean
  geom_point() +  # Posterior mean
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.3) +  # Credible interval
  theme_minimal() +
  labs(title = "Caterpillar Plot of p[i]",
       x = "Posterior Mean and 95% CI", y = "Index (i)")

ggplot(summary_p, aes(x = mean, y = reorder(i, mean), color = upper - lower)) +
  geom_point(size = 2) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.3) +
  scale_color_viridis_c() +  # Use a color scale to show uncertainty
  theme_minimal() +
  labs(title = "Caterpillar Plot with Uncertainty Highlighted",
       x = "Posterior Mean and 95% CI", y = "Index (i)", color = "CI Width")

## merge the posterior estimates with the plots df
summary_p %>% 
  left_join(.,recruits %>% 
              group_by(unique_plot_id,Endo_freq,Site,Species,Herbivory) %>% 
              summarise(n_samples=n()),
            by=c("i"="unique_plot_id"))->plot_init_final

## get mean regression parameters and generate prediction lines
spp_levels<-c("AGHY","ELVI","POAU")
site_levels<-c("BAS","BFL","COL","HUN","KER","LAF")
fence_levels<-c("Half","Full")

regression_params <- prevalence_samples %>%
  spread_draws(beta0[i,j,k], beta1[i,j,k]) %>%
  group_by(i, j, k) %>%
  summarise(
    beta0 = mean(beta0),
    beta1 = mean(beta1)
  ) %>%
  mutate(
    Species = factor(i, levels = 1:length(spp_levels), labels = spp_levels),
    Site = factor(j, levels = 1:length(site_levels), labels = site_levels),
    Fence = factor(k, levels = 1:length(fence_levels), labels = fence_levels)
  ) %>%
  select(Species, Site, Fence, beta0, beta1)

# Generate x values for the regression lines
x_vals <- seq(0, 1, length.out = 50)  # Match x-axis range
# Expand data to create a regression line for each facet
regression_lines <- expand.grid(Species=spp_levels,Site=site_levels,Fence=fence_levels,Endo_freq=x_vals) %>% 
  merge(.,regression_params,by=c("Species","Site","Fence"))
regression_lines$pred<-1/(1+exp(-(regression_lines$beta0+regression_lines$beta1*regression_lines$Endo_freq)))

##re-order sites from west to east
plot_init_final$Site <- factor(plot_init_final$Site, levels = c("KER", "BFL", "BAS", "COL", "HUN", "LAF"))
plot_init_final$Fence <- as.factor(ifelse(plot_init_final$Herbivory==1,"Full","Half"))

##drop site*spp*fence combos where we have no data
regression_lines<-semi_join(regression_lines,plot_init_final,by=c("Species","Site","Fence"))

ggplot(plot_init_final)+
  geom_abline(intercept=0,slope=1,linetype="dotted")+
  geom_line(data=regression_lines,aes(x=Endo_freq,y=pred,col=Fence),size=1) +
  geom_point(aes(x=Endo_freq,y=mean,col=Fence,size=n_samples),
             position = position_dodge(width = 0.3))+
  geom_errorbar(aes(x=Endo_freq,ymin=lower,ymax=upper,col=Fence),
                width=0,position = position_dodge(width = 0.3)) +
  facet_grid(rows=vars(Species),cols=vars(Site))+
  scale_size_continuous(breaks=c(10,30,50))+
  xlim(0,1)+ylim(0,1)+theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(x="Initial E+ prevalence",y="Final E+ prevaence")

##generate density histograms of predicted eqm prevalence
posterior_df <- prevalence_samples %>%
  spread_draws(beta0[i,j,k], beta1[i,j,k]) %>% 
  sample_n(1000)%>%
  mutate(
    Species = factor(i, levels = 1:length(spp_levels), labels = spp_levels),
    Site = factor(j, levels = 1:length(site_levels), labels = site_levels),
    Fence = factor(k, levels = 1:length(fence_levels), labels = fence_levels)
  )

##write function to generate eqm
logistic<-function(x){1/(1+exp(-x))}
find_eqm_prev <- function(beta0,beta1,start_prev,iter=500){
  for(t in 1:iter){
    if(t==1){prev<-start_prev}
    prev<-logistic(beta0+beta1*prev)
  }
  return(prev)
}

posterior_eqm <- posterior_df %>%
  mutate(eqm_prev05 = mapply(find_eqm_prev, beta0, beta1, start_prev=0.05),
         eqm_prev50 = mapply(find_eqm_prev, beta0, beta1, start_prev=0.50),
         eqm_prev95 = mapply(find_eqm_prev, beta0, beta1, start_prev=0.95)) 
##re-order sites from west to east
posterior_eqm$Site <- factor(posterior_eqm$Site, levels = c("KER", "BFL", "BAS", "COL", "HUN", "LAF"))
##drop site*spp*fence combos where we have no data
posterior_eqm<-semi_join(posterior_eqm,plot_init_final,by=c("Species","Site","Fence"))

ggplot(posterior_eqm)+
  geom_histogram(aes(x=eqm_prev95,fill=Fence))+
  facet_grid(rows=vars(Species),cols=vars(Site))+
  xlim(0,1)+theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(x="Equilibrium E+ prevalence")

##AGHY only
posterior_eqm %>% 
  filter(Species=="AGHY") %>% 
  ggplot()+
  geom_histogram(aes(x=eqm_prev95,fill=Fence))+
  facet_grid(cols=vars(Site))+
  xlim(0,1)+theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(x="Equilibrium E+ prevalence")

##Full fence only
posterior_eqm %>% 
  filter(Fence=="Full") %>% 
  ggplot()+
  geom_histogram(aes(x=eqm_prev95))+
  facet_grid(rows=vars(Species),cols=vaSrs(Site))+
  xlim(0,1)+theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(x="Equilibrium E+ prevalence")
