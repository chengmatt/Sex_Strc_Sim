# Purpose: To plot OM scenarios
# Creator: Matthew LH. Cheng (UAF - CFOS)
# Date: 10/19/23


# Set up ------------------------------------------------------------------

library(here)
library(tidyverse)
library(ggpubr)

# Set up quick theme
theme_tj = function() {
  theme_bw() +
    theme(legend.position = "none",
          axis.text = element_text(size = 13, color = "black"),
          axis.title = element_text(size = 15, color =  "black"),
          strip.text = element_text(size = 15),
          legend.title = element_text(size = 15),
          legend.text = element_text(size = 13),
          plot.title = element_text(size = 16))
}


# Experiment 1 ------------------------------------------------------------

# experiment 1 path
exp1_path = here("output", "Experiment 1")
# Read in scenarios
experiment_1_files = list.files(exp1_path)
experiment_1_files = experiment_1_files[!str_detect(experiment_1_files, "(30,10)|(10,30)")]
# Storage containers
vonB_all = data.frame()
waa_all = data.frame()
fishageselex_all = data.frame()
srvageselex_all = data.frame()
catch_sex_all = data.frame()
ssbs_all = data.frame()
M_all = data.frame()
total_biom_all = data.frame()
naa_store_all = data.frame()
sr_store_all = data.frame()

# Loop through to extract values
for(i in 1:length(experiment_1_files)) {
  
  # Load in OM object
  load(here(exp1_path, experiment_1_files[i], paste(experiment_1_files[i], ".RData", sep = "")))
  
  # Von bert function
  vonB_df = reshape2::melt(oms$vonB)
  names(vonB_df) = c("Age", "Sex", "Value")
  vonB_df$OM = experiment_1_files[i]
  vonB_all = rbind(vonB_all, vonB_df)
  
  # Weight at age
  waa_df = reshape2::melt(oms$waa)
  names(waa_df) = c("Age", "Sex", "Value")
  waa_df$OM = experiment_1_files[i]
  waa_all = rbind(waa_df, waa_all)
  
  # Age-based selectivity - fishery and survey
  fishageselex_df = reshape2::melt(oms$FishAge_Selex)
  names(fishageselex_df) = c("Age", "Sex", "Fleet", "Value")
  fishageselex_df$OM = experiment_1_files[i]
  fishageselex_all = rbind(fishageselex_all, fishageselex_df)
  
  # survey
  srvageselex_df = reshape2::melt(oms$SrvAge_Selex)
  names(srvageselex_df) = c("Age", "Sex", "Fleet", "Value")
  srvageselex_df$OM = experiment_1_files[i]
  srvageselex_all = rbind(srvageselex_df, srvageselex_all)
  
  # Catch
  catch_sex = reshape2::melt(oms$Total_Catch_Sex)
  names(catch_sex) = c("Years", "Sex", "Fleet", "Sim", "Catch")
  catch_sex$OM = experiment_1_files[i]
  catch_sex_all = rbind(catch_sex, catch_sex_all)
  
  # SSB
  ssbs = reshape2::melt(oms$SSB)
  names(ssbs) = c("Years", "Sim", "SSB")
  ssbs$OM = experiment_1_files[i]
  ssbs_all = rbind(ssbs, ssbs_all)
  
  # total biomass
  total_biom = reshape2::melt(oms$Total_Biom)
  names(total_biom) = c("Years", "Sim", "Biomass")
  total_biom$OM = experiment_1_files[i]
  total_biom_all = rbind(total_biom, total_biom_all)
  
  # Get natural mortality
  M_store = data.frame(M = rep(oms$M, (oms$n_years-1)), 
                 year = rep(1:(oms$n_years-1), rep(2, (oms$n_years-1))),
                 sex = c(1, 2), OM = experiment_1_files[i])
  
  M_all = rbind(M_store, M_all)
  
  # Get numbers at age
  naa_store = reshape2::melt(oms$NAA)
  names(naa_store) = c("Years", "Age", "Sex", "Sim", "Numbers")
  naa_store$OM = experiment_1_files[i]
  naa_store_all = rbind(naa_store, naa_store_all)
  
  # Get Sex ratio information over time
  sr_store = naa_store %>% 
    group_by(Years, Sim) %>% 
    mutate(Total_N = sum(Numbers)) %>% # get total numbers
    ungroup()
  
  # summarize across ages (need to take unique for some odd reason...)
  sr_store = sr_store %>% 
    group_by(Years, Sim, Sex, OM) %>% 
    summarize(Total_Sex = sum(Numbers) / Total_N) %>% 
    unique()
  sr_store_all = rbind(sr_store_all, sr_store)

} # end i loop


# Plot! -------------------------------------------------------------------

(vonB_plot = ggplot(vonB_all %>% 
                      mutate(Sex = ifelse(Sex == 1, "Female", 'Male')), 
                    aes(x = Age, y = Value, color = factor(Sex))) +
  geom_line(size = 2, alpha = 0.8) +
  facet_wrap(~OM, ncol = 1) +
  labs(x = "Age", y = "Value", color = "Sex", title = "Length-at-age") +
  theme_tj() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = c(0.65,0.9)))

waa_plot = ggplot(waa_all, aes(x = Age, y = Value, color = factor(Sex))) +
  geom_line(size = 2, alpha = 0.8) +
  facet_wrap(~OM, ncol = 1) +
  labs(x = "Age", y = "", color = "Sex", title = "Weight-at-age") +
  theme_tj() +
  theme(plot.title = element_text(hjust = 0.5))

fishageselex_plot = ggplot(fishageselex_all, aes(x = Age, y = Value, color = factor(Sex))) +
  geom_line(size = 2, alpha = 0.8) +
  facet_wrap(~OM, ncol = 1) +
  labs(x = "Age", y = " ", color = "Sex", title = "Fishery Selectivity") +
  theme_tj() +
  theme(plot.title = element_text(hjust = 0.5))

srvageselex_plot = ggplot(srvageselex_all, aes(x = Age, y = Value, color = factor(Sex))) +
  geom_line(size = 2, alpha = 0.8) +
  facet_wrap(~OM, ncol = 1) +
  labs(x = "Age", y = "", color = "Sex", title = "Survey Selectivity") +
  theme_tj() +
  theme(plot.title = element_text(hjust = 0.5))

catch_sex_all_sum = catch_sex_all %>% 
  filter(Catch != 0) %>% 
  group_by(Years, Sex, OM) %>% 
  summarize(Median = median(Catch),
            lwr_95 = quantile(Catch, 0.025),
            upr_95 = quantile(Catch, 0.975))

catch_plot = ggplot(catch_sex_all_sum,
       aes(x = Years, y = Median, color = factor(Sex), fill = factor(Sex),
           ymin = lwr_95, ymax = upr_95)) +
  geom_line(size = 1.5) +
  geom_ribbon(alpha = 0.35) +
  facet_wrap(~ OM, ncol = 1) +
  labs(x = "Years", y = "", color = "Sex", fill = "Sex", title = "Catch") +
  theme_tj() +
  theme(plot.title = element_text(hjust = 0.5))

ssb_sum = ssbs_all %>% 
  filter(Years != max(Years)) %>% 
  group_by(Years, OM) %>% 
  summarize(Median = median(SSB),
            lwr_95 = quantile(SSB, 0.025),
            upr_95 = quantile(SSB, 0.975))

(ssb_plot = ggplot(ssb_sum,
                    aes(x = Years, y = Median, ymin = lwr_95, ymax = upr_95)) +
  geom_line(size = 1.5) +
  geom_ribbon(alpha = 0.35) +
  facet_wrap(~OM, ncol = 1) +
  labs(x = "Years", y = "", color = "Sex", fill = "Sex", title = "Spawning Stock Biomass") +
  theme_tj() +
  theme(plot.title = element_text(hjust = 0.5)))

total_biom_sum = total_biom_all %>% 
  filter(Years != max(Years)) %>% 
  group_by(Years, OM) %>% 
  summarize(Median = median(Biomass),
            lwr_95 = quantile(Biomass, 0.025),
            upr_95 = quantile(Biomass, 0.975))

(biomass_plot = ggplot(total_biom_sum,
                   aes(x = Years, y = Median, ymin = lwr_95, ymax = upr_95)) +
    geom_line(size = 1.5) +
    geom_ribbon(alpha = 0.35) +
    facet_wrap(~OM, ncol = 1) +
    labs(x = "Years", y = "", color = "Sex", fill = "Sex", 
         title = "Total Biomass") +
    theme_tj() +
    theme(plot.title = element_text(hjust = 0.5)))

natmort_plot = ggplot(M_all, aes(x = year, y = M, color = factor(sex))) +
  geom_line(size = 1.5) +
  facet_wrap(~OM, ncol = 1) +
  labs(x = "Years", y = "", color = "Sex", fill = "Sex",
       title = "Natural Mortality") +
  theme_tj() +
  ylim(0.05, 0.15) +
  theme(plot.title = element_text(hjust = 0.5))

naa_sum = naa_store_all %>% 
  filter(Years != max(Years)) %>% 
  group_by(Years, Sex, OM, Sim) %>% 
  summarize(naa_sum = sum(Numbers)) %>% 
  group_by(Years, Sex, OM) %>% 
  summarize(Median = median(naa_sum),
            lwr_95 = quantile(naa_sum, 0.025),
            upr_95 = quantile(naa_sum, 0.975))

naa_plot = ggplot(naa_sum, aes(x = Years, y = Median,
                               ymin = lwr_95, ymax = upr_95, color = factor(Sex),
                               fill = factor(Sex))) +
  geom_line(size = 1.5) +
  geom_ribbon(alpha = 0.35) +
  facet_wrap(~OM, ncol = 1) +
  labs(x = "Years", y = "", color = "Sex", fill = "Sex", 
       title = "Total Numbers") +
  theme_tj() +
  theme(plot.title = element_text(hjust = 0.5))

# plot sex ratios changing over time
sr_sum = sr_store_all %>% 
  group_by(Years, Sex, OM) %>% 
  summarize(median = median(Total_Sex),
            lwr_95 = quantile(Total_Sex, 0.025),
            upr_95 = quantile(Total_Sex, 0.975))

# sex ratio across time summary
sr_sum_plot = ggplot(sr_sum %>% filter(!str_detect(OM, "No")), 
                     aes(x = Years, y = median, color = factor(Sex), fill = factor(Sex), 
                         ymin = lwr_95, ymax = upr_95)) +
  geom_line(size = 1.5) +
  geom_ribbon(alpha = 0.35) +
  facet_wrap(~OM, ncol = 1) +
  labs(x = "Years", y = "", color = "Sex", fill = "Sex", title = "Population Sex Ratio") +
  theme_tj() +
  theme(plot.title = element_text(hjust = 0.5))


pdf(here("figs", "OM_Exp1.pdf"), width = 30, height = 10)
ggarrange(vonB_plot, waa_plot, fishageselex_plot, 
          srvageselex_plot, natmort_plot, sr_sum_plot, 
          naa_plot, catch_plot, ssb_plot, biomass_plot,
          nrow = 1)
dev.off()



# Experiment 2 ------------------------------------------------------------

# experiment 2 path
exp2_path = here("output", "Experiment 2")
# Read in scenarios
experiment_2_files = list.files(exp2_path)

# Storage containers
fishageselex_all = data.frame()
catch_sex_all = data.frame()
ssbs_all = data.frame()
sr_all = data.frame()
total_biom_all = data.frame()
naa_store_all = data.frame()
vonB_all = data.frame()
waa_all = data.frame()
sr_store_all = data.frame()

# Loop through to extract values
for(i in 1:length(experiment_2_files)) {
  
  # Load in OM object
  load(here(exp2_path, experiment_2_files[i], paste(experiment_2_files[i], ".RData", sep = "")))
  
  # Von bert function
  vonB_df = reshape2::melt(oms$vonB)
  names(vonB_df) = c("Age", "Sex", "Value")
  vonB_df$OM = experiment_2_files[i]
  vonB_all = rbind(vonB_all, vonB_df)
  
  # Weight at age
  waa_df = reshape2::melt(oms$waa)
  names(waa_df) = c("Age", "Sex", "Value")
  waa_df$OM = experiment_2_files[i]
  waa_all = rbind(waa_df, waa_all)
  
  # Age-based selectivity - fishery and survey
  fishageselex_df = reshape2::melt(oms$FishAge_Selex)
  names(fishageselex_df) = c("Age", "Sex", "Fleet", "Value")
  fishageselex_df$OM = experiment_2_files[i]
  fishageselex_all = rbind(fishageselex_all, fishageselex_df)
  
  # Catch
  catch_sex = reshape2::melt(oms$Total_Catch_Sex)
  names(catch_sex) = c("Years", "Sex", "Fleet", "Sim", "Catch")
  catch_sex$OM = experiment_2_files[i]
  catch_sex_all = rbind(catch_sex, catch_sex_all)
  
  # SSB
  ssbs = reshape2::melt(oms$SSB)
  names(ssbs) = c("Years", "Sim", "SSB")
  ssbs$OM = experiment_2_files[i]
  ssbs_all = rbind(ssbs, ssbs_all)
  
  # total biomass
  total_biom = reshape2::melt(oms$Total_Biom)
  names(total_biom) = c("Years", "Sim", "Biomass")
  total_biom$OM = experiment_2_files[i]
  total_biom_all = rbind(total_biom, total_biom_all)
  
  # Get numbers at age
  naa_store = reshape2::melt(oms$NAA)
  names(naa_store) = c("Years", "Age", "Sex", "Sim", "Numbers")
  naa_store$OM = experiment_2_files[i]
  naa_store_all = rbind(naa_store, naa_store_all)

  # Get Sex ratio information over time
  sr_store = naa_store %>% 
    group_by(Years, Sim) %>% 
    mutate(Total_N = sum(Numbers)) %>% # get total numbers
    ungroup()
  
  # summarize across ages (need to take unique for some odd reason...)
  sr_store = sr_store %>% 
    group_by(Years, Sim, Sex, OM) %>% 
    summarize(Total_Sex = sum(Numbers) / Total_N) %>% 
    unique()
  sr_store_all = rbind(sr_store_all, sr_store)
  
} # end i loop

# length at age plot
(vonB_plot = ggplot(vonB_all %>% 
                      filter(!str_detect(OM, "No")) %>% 
                      mutate(Sex = ifelse(Sex == 1, "Female", 'Male')), 
                    aes(x = Age, y = Value, color = factor(Sex))) +
    geom_line(size = 2, alpha = 0.8) +
    facet_wrap(~OM, ncol = 1) +
    labs(x = "Age", y = "Value", color = "Sex", title = "Length-at-age") +
    theme_tj() +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = c(0.75,0.85)))

# weight at age plot
waa_plot = ggplot(waa_all %>% 
                    filter(!str_detect(OM, "No")), 
                  aes(x = Age, y = Value, color = factor(Sex))) +
  geom_line(size = 2, alpha = 0.8) +
  facet_wrap(~OM, ncol = 1) +
  labs(x = "Age", y = "", color = "Sex", title = "Weight-at-age") +
  theme_tj() +
  theme(plot.title = element_text(hjust = 0.5))

# fishery selex age plot
fishageselex_plot = ggplot(fishageselex_all %>% 
                             filter(!str_detect(OM, "No")) %>% 
                             mutate(Sex = ifelse(Sex == 1, "Female", 'Male')), 
                           aes(x = Age, y = Value, color = factor(Sex))) +
  geom_line(size = 2, alpha = 0.8) +
  facet_wrap(~OM, ncol = 1) +
  labs(x = "Age", y = "", color = "Sex", title = "Fishery Selectivity") +
  theme_tj() +
  theme(plot.title = element_text(hjust = 0.5))

# sumamrize catch
catch_sex_all_sum = catch_sex_all %>% 
  filter(Catch != 0) %>% 
  group_by(Years, Sex, OM) %>% 
  summarize(Median = median(Catch),
            lwr_95 = quantile(Catch, 0.025),
            upr_95 = quantile(Catch, 0.975))

# catch plot
catch_plot = ggplot(catch_sex_all_sum %>% 
                      filter(!str_detect(OM, "No")),
                    aes(x = Years, y = Median, color = factor(Sex), fill = factor(Sex),
                        ymin = lwr_95, ymax = upr_95)) +
  geom_line(size = 1.5) +
  geom_ribbon(alpha = 0.35) +
  facet_wrap(~ OM, ncol = 1) +
  labs(x = "Years", y = "", color = "Sex", fill = "Sex", title = "Catch") +
  theme_tj() +
  theme(plot.title = element_text(hjust = 0.5))

# summarize ssb
ssb_sum = ssbs_all %>% 
  filter(Years != max(Years)) %>% 
  group_by(Years, OM) %>% 
  summarize(Median = median(SSB),
            lwr_95 = quantile(SSB, 0.025),
            upr_95 = quantile(SSB, 0.975))

# ssb plot
(ssb_plot = ggplot(ssb_sum %>% 
                     filter(!str_detect(OM, "No")),
                   aes(x = Years, y = Median, ymin = lwr_95, ymax = upr_95)) +
    geom_line(size = 1.5) +
    geom_ribbon(alpha = 0.35) +
    facet_wrap(~OM, ncol = 1) +
    labs(x = "Years", y = "", color = "Sex", fill = "Sex", title = "Spawning Stock Biomass") +
    theme_tj() +
    theme(plot.title = element_text(hjust = 0.5)))

# summarize biomass
total_biom_sum = total_biom_all %>% 
  filter(Years != max(Years)) %>% 
  group_by(Years, OM) %>% 
  summarize(Median = median(Biomass),
            lwr_95 = quantile(Biomass, 0.025),
            upr_95 = quantile(Biomass, 0.975))

# biomass plot
(biomass_plot = ggplot(total_biom_sum %>% 
                         filter(!str_detect(OM, "No")),
                       aes(x = Years, y = Median, ymin = lwr_95, ymax = upr_95)) +
    geom_line(size = 1.5) +
    geom_ribbon(alpha = 0.35) +
    facet_wrap(~OM, ncol = 1) +
    labs(x = "Years", y = "", color = "Sex", fill = "Sex", 
         title = "Total Biomass") +
    theme_tj() +
    theme(plot.title = element_text(hjust = 0.5)))

# summarize total numbers
naa_sum = naa_store_all %>% 
  filter(Years != max(Years)) %>% 
  group_by(Years, Sex, OM, Sim) %>% 
  summarize(naa_sum = sum(Numbers)) %>% 
  group_by(Years, Sex, OM) %>% 
  summarize(Median = median(naa_sum),
            lwr_95 = quantile(naa_sum, 0.025),
            upr_95 = quantile(naa_sum, 0.975))

# total numbers plot
naa_plot = ggplot(naa_sum %>% 
                    filter(!str_detect(OM, "No")), aes(x = Years, y = Median,
                    ymin = lwr_95, ymax = upr_95, color = factor(Sex),
                    fill = factor(Sex))) +
  geom_line(size = 1.5) +
  geom_ribbon(alpha = 0.35) +
  facet_wrap(~OM, ncol = 1) +
  labs(x = "Years", y = "", color = "Sex", fill = "Sex", 
       title = "Total Numbers") +
  theme_tj() +
  theme(plot.title = element_text(hjust = 0.5))


# plot sex ratios changing over time
sr_sum = sr_store_all %>% 
  group_by(Years, Sex, OM) %>% 
  summarize(median = median(Total_Sex),
            lwr_95 = quantile(Total_Sex, 0.025),
            upr_95 = quantile(Total_Sex, 0.975))

# sex ratio across time summary
sr_sum_plot = ggplot(sr_sum %>% filter(!str_detect(OM, "No")), 
       aes(x = Years, y = median, color = factor(Sex), fill = factor(Sex), 
           ymin = lwr_95, ymax = upr_95)) +
  geom_line(size = 1.5) +
  geom_ribbon(alpha = 0.35) +
  facet_wrap(~OM, ncol = 1) +
  labs(x = "Years", y = "", color = "Sex", fill = "Sex", title = "Population Sex Ratio") +
  theme_tj() +
  theme(plot.title = element_text(hjust = 0.5))

pdf(here("figs", "OM_Exp2.pdf"), width = 32, height = 12)
ggarrange(vonB_plot, waa_plot, fishageselex_plot, sr_sum_plot,
          naa_plot, catch_plot, ssb_plot, biomass_plot, nrow = 1)
dev.off()

