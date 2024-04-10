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
          axis.title = element_text(size = 13, color =  "black"),
          strip.text = element_text(size = 13),
          legend.title = element_text(size = 13),
          legend.text = element_text(size = 13),
          plot.title = element_text(size = 13))
}


# Experiment 1 ------------------------------------------------------------

# experiment 1 path
exp1_path = here("output", "Experiment 1")
# Read in scenarios
experiment_1_files = list.files(exp1_path)

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

waa_plot = ggplot(waa_all %>% 
                    filter(OM == "Across_100") %>% 
                    mutate(Sex = ifelse(Sex == 1, "Female", 'Male')),
                  aes(x = Age, y = Value, color = factor(Sex))) +
  geom_line(size = 2, alpha = 0.8) +
  labs(x = "Age", y = "", color = "Sex", title = "Weight-at-age") +
  theme_tj() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = c(0.2,0.9))

fishageselex_plot = ggplot(fishageselex_all %>% 
                             filter(OM == "Across_100"),
                           aes(x = Age, y = Value, color = factor(Sex))) +
  geom_line(size = 2, alpha = 0.8) +
  labs(x = "Age", y = " ", color = "Sex", title = "Fishery Selectivity") +
  theme_tj() +
  theme(plot.title = element_text(hjust = 0.5))

catch_sex_all_sum = catch_sex_all %>% 
  filter(Catch != 0) %>% 
  group_by(Years, Sex, OM) %>% 
  summarize(Median = median(Catch),
            lwr_95 = quantile(Catch, 0.025),
            upr_95 = quantile(Catch, 0.975))

catch_plot = ggplot(catch_sex_all_sum %>% 
                      filter(OM == "Across_100"),
       aes(x = Years, y = Median, color = factor(Sex), fill = factor(Sex),
           ymin = lwr_95, ymax = upr_95)) +
  geom_line(size = 1.5) +
  geom_ribbon(alpha = 0.35, color = NA) +
  labs(x = "Years", y = "", color = "Sex", fill = "Sex", title = "Catch") +
  theme_tj() +
  theme(plot.title = element_text(hjust = 0.5))

ssb_sum = ssbs_all %>% 
  filter(Years != max(Years)) %>% 
  group_by(Years, OM) %>% 
  summarize(Median = median(SSB),
            lwr_95 = quantile(SSB, 0.025),
            upr_95 = quantile(SSB, 0.975))

(ssb_plot = ggplot(ssb_sum %>% 
                     filter(OM == "Across_100"),
                    aes(x = Years, y = Median, ymin = lwr_95, ymax = upr_95)) +
  geom_line(size = 1.5) +
  geom_ribbon(alpha = 0.35, color = NA) +
  labs(x = "Years", y = "", color = "Sex", fill = "Sex", title = "Spawning Stock Biomass") +
  theme_tj() +
  theme(plot.title = element_text(hjust = 0.5)))

natmort_plot = ggplot(M_all %>% 
                        filter(OM == "Across_100"), 
                      aes(x = year, y = M, color = factor(sex))) +
  geom_line(size = 1.5) +
  labs(x = "Years", y = "", color = "Sex", fill = "Sex",
       title = "Natural Mortality") +
  theme_tj() +
  ylim(0.05, 0.15) +
  theme(plot.title = element_text(hjust = 0.5))


# Get example of proportions within vs across
prop_within = naa_store_all %>% 
  filter(Years == 20, OM == "Across_100", Sim == 1) %>% 
  group_by(Years, Sex) %>% 
  mutate(Prop = Numbers/sum(Numbers),
         Type = "Within")

# Get example of proportions within vs across
prop_across = naa_store_all %>% 
  filter(Years == 20, OM == "Across_100", Sim == 1) %>% 
  group_by(Years) %>% 
  mutate(Prop = Numbers/sum(Numbers),
         Type = "Across")

# bind
prop_df <- rbind(prop_within, prop_across)
prop_plot <- ggplot(prop_df, aes(x = Age, y = Prop, color = factor(Sex))) +
  geom_line(size = 2, alpha = 0.85) +
  labs(x = "Age", y = "Proportion", color = "Sex", fill = "Sex",
       title = "Exp 1 (Treatment of Proportions)") +
  scale_color_manual(labels = c("Female", "Male"), 
                     values = c("#DC3220", "#005AB5")) +
  facet_wrap(~Type) +
  theme_tj() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = c(0.1, 0.8))

# plot sex ratios changing over time
sr_sum = sr_store_all %>% 
  group_by(Years, Sex, OM) %>% 
  summarize(median = median(Total_Sex),
            lwr_95 = quantile(Total_Sex, 0.025),
            upr_95 = quantile(Total_Sex, 0.975))

# sex ratio across time summary
sr_sum_plot = ggplot(sr_sum %>% filter(!str_detect(OM, "No")) %>% 
                       filter(OM == "Across_100"), 
                     aes(x = Years, y = median, color = factor(Sex), fill = factor(Sex), 
                         ymin = lwr_95, ymax = upr_95)) +
  geom_line(size = 1.5) +
  geom_ribbon(alpha = 0.35, color = NA) +
  labs(x = "Years", y = "", color = "Sex", fill = "Sex", title = "Population Sex Ratio") +
  theme_tj() +
  theme(plot.title = element_text(hjust = 0.5))


pdf(here("figs", "OM_Exp1.pdf"), width = 20, height = 12)
ggarrange(waa_plot, fishageselex_plot, 
          natmort_plot ,prop_within_plot,
          sr_sum_plot, catch_plot, ssb_plot,  prop_across_plot, nrow = 2, ncol = 4)
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
M_all = data.frame()

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
  
  
  # Get natural mortality
  M_store = data.frame(M = rep(oms$M, (oms$n_years-1)), 
                       year = rep(1:(oms$n_years-1), rep(2, (oms$n_years-1))),
                       sex = c(1, 2), OM = experiment_2_files[i])
  
  M_all = rbind(M_store, M_all)
  
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

# weight at age plot
waa_plot = ggplot(waa_all %>% 
                    filter(!str_detect(OM, "No")) %>% 
                    mutate(Sex = ifelse(Sex == 1, "Female", 'Male')), 
                  aes(x = Age, y = Value, color = factor(Sex))) +
  geom_line(size = 2, alpha = 0.8) +
  scale_color_manual(values = c("#DC3220", "#005AB5")) +
  facet_wrap(~OM, ncol = 1) +
  labs(x = "Age", y = "Weight-at-age", color = "Sex", title = "") +
  theme_tj() +
  theme(plot.title = element_text(hjust = 0.5))

# fishery selex age plot
fishageselex_plot = ggplot(fishageselex_all %>% 
                             filter(!str_detect(OM, "No")) %>% 
                             mutate(Sex = ifelse(Sex == 1, "Female", 'Male')), 
                           aes(x = Age, y = Value, color = factor(Sex))) +
  geom_line(size = 2, alpha = 0.8) +
  scale_color_manual(values = c("#DC3220", "#005AB5")) +
  facet_wrap(~OM, ncol = 1) +
  labs(x = "Age", y = "Fishery Selectivity", color = "Sex", 
       title = "Exp 2 (Sex Dimorphism & Sex-Specific Catch)") +
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
  scale_color_manual(values = c("#DC3220", "#005AB5")) +
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
    scale_color_manual(values = c("#DC3220", "#005AB5")) +
    facet_wrap(~OM, ncol = 1) +
    labs(x = "Years", y = "", color = "Sex", fill = "Sex", title = "Spawning Stock Biomass") +
    theme_tj() +
    theme(plot.title = element_text(hjust = 0.5)))

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
  geom_ribbon(alpha = 0.35, color = NA) +
  scale_color_manual(values = c("#DC3220", "#005AB5")) +
  facet_wrap(~OM, ncol = 1) +
  labs(x = "Years", y = "", color = "Sex", fill = "Sex", title = "Population Sex Ratio") +
  theme_tj() +
  theme(plot.title = element_text(hjust = 0.5))

natmort_plot = ggplot(M_all, 
                      aes(x = year, y = M, color = factor(sex))) +
  geom_line(size = 1.5) +
  facet_wrap(~OM, ncol = 1) +
  labs(x = "Years", y = "Natural Mortality", color = "Sex", fill = "Sex",
       title = "") +
  theme_tj() +
  scale_color_manual(values = c("#DC3220", "#005AB5")) +
  ylim(0.05, 0.15) +
  theme(plot.title = element_text(hjust = 0.5))

pdf(here("figs", "OM_Exp2.pdf"), width = 17, height = 12)
ggarrange(waa_plot, fishageselex_plot, natmort_plot, 
          sr_sum_plot, catch_plot, ssb_plot, nrow = 2, ncol = 3)
dev.off()

# Make experiment 2 plots
exp2_plots <- (waa_plot | fishageselex_plot | natmort_plot) +
  theme(plot.title = element_text(hjust = .5, size = 16))

# Experiment 3 ------------------------------------------------------------

# experiment 3 path
exp3_path = here("output", "Experiment 3")
# Read in scenarios
experiment_3_files = list.files(exp3_path)

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
M_all = data.frame()

# Loop through to extract values
for(i in 1:length(experiment_3_files)) {
  
  # Load in OM object
  load(here(exp3_path, experiment_3_files[i], paste(experiment_3_files[i], ".RData", sep = "")))
  
  # Von bert function
  vonB_df = reshape2::melt(oms$vonB)
  names(vonB_df) = c("Age", "Sex", "Value")
  vonB_df$OM = experiment_3_files[i]
  vonB_all = rbind(vonB_all, vonB_df)
  
  # Weight at age
  waa_df = reshape2::melt(oms$waa)
  names(waa_df) = c("Age", "Sex", "Value")
  waa_df$OM = experiment_3_files[i]
  waa_all = rbind(waa_df, waa_all)
  
  # Age-based selectivity - fishery and survey
  fishageselex_df = reshape2::melt(oms$FishAge_Selex)
  names(fishageselex_df) = c("Age", "Sex", "Fleet", "Value")
  fishageselex_df$OM = experiment_3_files[i]
  fishageselex_all = rbind(fishageselex_all, fishageselex_df)
  
  # Catch
  catch_sex = reshape2::melt(oms$Total_Catch_Sex)
  names(catch_sex) = c("Years", "Sex", "Fleet", "Sim", "Catch")
  catch_sex$OM = experiment_3_files[i]
  catch_sex_all = rbind(catch_sex, catch_sex_all)
  
  # SSB
  ssbs = reshape2::melt(oms$SSB)
  names(ssbs) = c("Years", "Sim", "SSB")
  ssbs$OM = experiment_3_files[i]
  ssbs_all = rbind(ssbs, ssbs_all)
  
  
  # Get natural mortality
  M_store = data.frame(M = rep(oms$M, (oms$n_years-1)), 
                       year = rep(1:(oms$n_years-1), rep(2, (oms$n_years-1))),
                       sex = c(1, 2), OM = experiment_3_files[i])
  
  M_all = rbind(M_store, M_all)
  
  # total biomass
  total_biom = reshape2::melt(oms$Total_Biom)
  names(total_biom) = c("Years", "Sim", "Biomass")
  total_biom$OM = experiment_3_files[i]
  total_biom_all = rbind(total_biom, total_biom_all)
  
  # Get numbers at age
  naa_store = reshape2::melt(oms$NAA)
  names(naa_store) = c("Years", "Age", "Sex", "Sim", "Numbers")
  naa_store$OM = experiment_3_files[i]
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
          legend.position = c(0.2,0.85)))

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
  geom_ribbon(alpha = 0.35, color = NA) +
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
    geom_ribbon(alpha = 0.35, color = NA) +
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
    geom_ribbon(alpha = 0.35, color = NA) +
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
  geom_ribbon(alpha = 0.35, color = NA) +
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
  geom_ribbon(alpha = 0.35, color = NA) +
  facet_wrap(~OM, ncol = 1) +
  labs(x = "Years", y = "", color = "Sex", fill = "Sex", title = "Population Sex Ratio") +
  theme_tj() +
  theme(plot.title = element_text(hjust = 0.5))

# Get natural mortality
M_store = data.frame(M = rep(oms$M, (oms$n_years-1)), 
                     year = rep(1:(oms$n_years-1), rep(2, (oms$n_years-1))),
                     sex = c(1, 2), OM = experiment_1_files[i])

M_all = rbind(M_store, M_all)

natmort_plot = ggplot(M_all %>% filter(!str_detect(OM, "No")), 
                      aes(x = year, y = M, color = factor(sex))) +
  geom_line(size = 1.5) +
  facet_wrap(~OM, ncol = 1) +
  labs(x = "Years", y = "", color = "Sex", fill = "Sex",
       title = "Natural Mortality") +
  theme_tj() +
  ylim(0.05, 0.15) +
  theme(plot.title = element_text(hjust = 0.5))

pdf(here("figs", "OM_Exp3.pdf"), width = 17, height = 12)
ggarrange(vonB_plot, waa_plot, fishageselex_plot, natmort_plot, nrow = 1, ncol = 4)
ggarrange(sr_sum_plot, naa_plot, catch_plot, ssb_plot, nrow = 1, ncol = 4)
dev.off()

exp3_sr_plot <- ggplot(sr_sum %>% filter(!str_detect(OM, "No")), 
       aes(x = Years, y = median, color = factor(Sex), fill = factor(Sex), 
           ymin = lwr_95, ymax = upr_95)) +
  geom_line(size = 1.5) +
  scale_color_manual(values = c("#DC3220", "#005AB5")) +
  scale_fill_manual(values = c("#DC3220", "#005AB5")) +
  geom_ribbon(alpha = 0.35, color = NA) +
  facet_wrap(~OM, ncol = 1) +
  labs(x = "Years", y = "Population Sex Ratio", color = "Sex", 
       fill = "Sex", title = "Exp 3 (Sex Ratio Estimation and Misspecification)") +
  theme_tj() +
  theme(plot.title = element_text(hjust = 0.5))

# Combine OM plot ---------------------------------------------------------

om_plots <- ggarrange(prop_plot, exp2_plots, 
                      widths = c(0.5, 0.5), ncol = 1, heights = c(0.3, 0.7),
                      labels = c("A", "B"), font.label = list(size = 25))

pdf(here("figs", "ms_figs", "Fig1_OMScenarios.pdf"), width = 20, height = 13)
ggarrange(om_plots, exp3_sr_plot, widths = c(0.65, 0.35), 
          labels = c("", "C"), font.label = list(size = 25))
dev.off()

