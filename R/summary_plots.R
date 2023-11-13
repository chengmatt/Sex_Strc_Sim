# Purpose: To create summary (relative error) plots for simulations
# Creator: Matthew LH. Cheng
# date: 10/20/23


# Set up ------------------------------------------------------------------

library(here)
library(tidyverse)
library(ggpubr)

# Set up quick theme
theme_tj = function() {
  theme_bw() +
    theme(legend.position = "none",
          axis.text = element_text(size = 18, color = "black"),
          axis.title = element_text(size = 20, color =  "black"),
          strip.text = element_text(size = 13),
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 18),
          plot.title = element_text(size = 20))
}

# Experiment 1 ------------------------------------------------------------

# Read in files
exp1_selex_df = data.table::fread(here("output", "Experiment_1_Selex.csv"))
exp1_growth_df = data.table::fread(here("output", "Experiment_1_Growth.csv"))
exp1_param_df = data.table::fread(here("output", "Experiment_1_Param.csv"))
exp1_ts_df = data.table::fread(here("output", "Experiment_1_TimeSeries.csv"))
exp1_conv_df = data.table::fread(here("output", "Experiment_1_Convergence.csv"))
exp1_naa_df = data.table::fread(here("output", "Experiment_1_NAA.csv"))

### Convergence Summary -----------------------------------------------------

# Convergence summary
conv_df = exp1_conv_df %>% 
  filter(convergence == "Converged") %>% 
  group_by(OM, EM) %>% 
  summarize(sum = n())

# Plot convergence
pdf(here("figs", "Experiment 1", "Convergence.pdf"), width = 15)
ggplot(conv_df, aes(x = OM, y = sum)) +
  geom_point(size = 3) +
  facet_wrap(~EM) +
  theme_tj() +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  labs(x = "Operating Models", y = "Convergence") +
  ylim(150, 200)
dev.off()

### Selectivity Summary -----------------------------------------------------

# summarize relative error in selectivity
exp1_selex_sum = exp1_selex_df %>% 
  filter(Convergence == "Converged") %>% 
  mutate(RE = (Pred - True) / True,
         Magg = case_when( 
           str_detect(EM, "MAgg") ~ "Aggregated M",
           str_detect(EM, "MSex") ~ "Sex-Specific M"),
         Prop = case_when(
           str_detect(EM, "PropAcr") ~ "Proportions Across",
           str_detect(EM, "PropWith_SR_Y") ~ "Proportions Within (SR_Y)",
           str_detect(EM, "PropWith_SR_ALY") ~ "Proportions Within (SR_ALY)"),
         Catch = case_when(
           str_detect(EM, "CatAgg") ~ "Aggregated Catch",
           str_detect(EM, "CatSex") ~ "Sex-Specific Catch")) %>% 
  group_by(OM, EM, Age, Sex, Type, Magg, Prop, Catch) %>% 
  summarize(Median = median(RE),
            lwr_95 = quantile(RE, 0.025),
            upr_95 = quantile(RE, 0.975))

# Fishery Selectivity
pdf(here("figs", "Experiment 1", "RE_Fish_Selex.pdf"), width = 34, height = 15)
fish_selex = print(
  ggplot(exp1_selex_sum %>% filter(Type == "Fishery Selectivity") %>% 
           mutate(Sex = ifelse(Sex == 1, "Female", "Male")), 
         aes(x = Age, y = Median, ymin = lwr_95, ymax = upr_95,
             fill = factor(Sex), color = factor(Sex), group = factor(Sex))) +
    geom_line(size = 3) +
    geom_ribbon(alpha = 0.35) +
    geom_hline(yintercept = 0, lty = 2, size = 1.3) + 
    labs(x = "Age", y = "Relative Error in Fishery Selectivity", 
         color = "Sex", fill = "Sex") +
    facet_grid(OM~EM, scales = "free") +
    theme_tj() +
    theme(legend.position = "top")
)
dev.off()

pdf(here("figs", "Experiment 1", "RE_Fish_Selex_Magg.pdf"), height = 10)

print(
  ggplot(exp1_selex_sum %>% filter(Type == "Fishery Selectivity",
                                   Magg == "Aggregated M",
                                   Prop == 'Proportions Across') %>% 
           mutate(Sex = ifelse(Sex == 1, "Female", "Male")), 
         aes(x = Age, y = Median, ymin = lwr_95, ymax = upr_95,
             fill = factor(Sex), color = factor(Sex), group = factor(Sex))) +
    geom_line(size = 3) +
    geom_ribbon(alpha = 0.35) +
    geom_hline(yintercept = 0, lty = 2, size = 1.3) + 
    labs(x = "Age", y = "Relative Error in Fishery Selectivity", 
         color = "Sex", fill = "Sex", title = "Aggregated M, Proportions Across") +
    facet_grid(OM~EM, scales = "free") +
    theme_tj() +
    theme(legend.position = "top")
)

print(
  ggplot(exp1_selex_sum %>% filter(Type == "Fishery Selectivity",
                                   Magg == "Aggregated M",
                                   Prop == 'Proportions Across') %>% 
           mutate(Sex = ifelse(Sex == 1, "Female", "Male")), 
         aes(x = Age, y = Median, ymin = lwr_95, ymax = upr_95,
             fill = factor(Sex), color = factor(Sex), group = factor(Sex))) +
    geom_line(size = 3) +
    geom_ribbon(alpha = 0.35) +
    geom_hline(yintercept = 0, lty = 2, size = 1.3) + 
    labs(x = "Age", y = "Relative Error in Fishery Selectivity", 
         color = "Sex", fill = "Sex", title = "Aggregated M, Proportions Across, Zoomed") +
    facet_grid(OM~EM, scales = "free") +
    theme_tj() +
    coord_cartesian(ylim = c(-0.5, 0.5)) +
    theme(legend.position = "top")
)

dev.off()


# fishery selectivity for age only EM
pdf(here("figs", "Experiment 1", "RE_Fish_Selex_AgeEM.pdf"), width = 10)
print(
  exp1_selex_df %>% 
    filter(Convergence == "Converged", EM == "Age",
           Type == "Fishery Selectivity") %>% 
    mutate(Sex = ifelse(Sex == 1, "Female", "Male")) %>% 
    ggplot() +
    geom_line(aes(x = Age, y = Pred, group = sim)) +
    geom_line(aes(x = Age, y = True, color = factor(Sex)), size = 1.3) +
    facet_wrap(~OM) +
    labs(x = "Age", y = "Fishery Selectivity", color = "Sex") +
    theme_tj() +
    theme(legend.position = "top")
)
dev.off()

# Survey Selectivity
pdf(here("figs", "Experiment 1", "RE_Srv_Selex.pdf"), width = 34, height = 15)
srv_selex = print(
  ggplot(exp1_selex_sum %>% filter(Type == "Survey Selectivity") %>% 
           mutate(Sex = ifelse(Sex == 1, "Female", "Male")), 
         aes(x = Age, y = Median, ymin = lwr_95, ymax = upr_95,
             fill = factor(Sex), color = factor(Sex), group = factor(Sex))) +
    geom_line(size = 3) +
    geom_ribbon(alpha = 0.35) +
    geom_hline(yintercept = 0, lty = 2, size = 1.3) + 
    labs(x = "Age", y = "Relative Error in Survey Selectivity", color = "Sex", fill = "Sex") +
    facet_grid(OM~EM, scales = "free") +
    theme_tj() +
    theme(legend.position = "top")
)
dev.off()

# survey selectivity for age only EM
pdf(here("figs", "Experiment 1", "RE_Srv_Selex_AgeEM.pdf"), width = 10)
print(
  exp1_selex_df %>% 
    filter(Convergence == "Converged", EM == "Age",
           Type == "Survey Selectivity") %>% 
    mutate(Sex = ifelse(Sex == 1, "Female", "Male")) %>% 
    ggplot() +
    geom_line(aes(x = Age, y = Pred, group = sim)) +
    geom_line(aes(x = Age, y = True, color = factor(Sex)), size = 1.3) +
    facet_wrap(~OM) +
    labs(x = "Age", y = "Survey Selectivity", color = "Sex") +
    theme_tj() +
    theme(legend.position = "top")
)
dev.off()


### Growth Summary ----------------------------------------------------------

# Plot growth curves and estimates
pdf(here("figs", "Experiment 1", "RE_Growth_Age.pdf"))
print(
  exp1_growth_df %>% 
    filter(Convergence == "Converged", EM == "Age") %>% 
    mutate(Sex = ifelse(Sex == 1, "Female", "Male")) %>% 
    ggplot() +
    geom_line(aes(x = Age, y = Pred, group = sim)) +
    geom_line(aes(x = Age, y = True, color = factor(Sex)), size = 1.3) +
    facet_wrap(~OM) +
    labs(x = "Age", y = "WAA", color = "Sex") +
    theme_tj() +
    theme(legend.position = "top")
)
dev.off()
  
# summarize relative error in growth
exp1_growth_sum = exp1_growth_df %>% 
  filter(Convergence == "Converged") %>% 
  mutate(RE = (Pred - True) / True) %>% 
  group_by(OM, EM, Age, Sex) %>% 
  summarize(Median = median(RE),
            lwr_95 = quantile(RE, 0.025),
            upr_95 = quantile(RE, 0.975))


# plot relative error in growth
pdf(here("figs", "Experiment 1", "RE_Growth.pdf"), width = 34, height = 15)
print(
  ggplot(exp1_growth_sum %>%  mutate(Sex = ifelse(Sex == 1, "Female", "Male")), 
         aes(x = Age, y = Median, ymin = lwr_95, ymax = upr_95,
             fill = factor(Sex), color = factor(Sex), group = factor(Sex))) +
    geom_line(size = 3) +
    geom_ribbon(alpha = 0.35) +
    geom_hline(yintercept = 0, lty = 2, size = 1.3) + 
    labs(x = "Age", y = "Relative Error in Growth (WAA)", color = "Sex", fill = "Sex") +
    facet_grid(OM~EM, scales = "free") +
    theme_tj() +
    theme(legend.position = "top")
)
dev.off()

### Parameter Summary -------------------------------------------------------
##### Proportions across vs. within -------------------------------------------
# proportions within vs. across summarize and munge to get em components
prop_param_df = exp1_param_df %>% 
  filter(Convergence == "Converged",
         !str_detect(Type, "Ratio")) %>% 
  mutate(RE = (as.numeric(Pred)-Truth)/Truth,
         Magg = case_when( 
           str_detect(EM, "MAgg") ~ "Aggregated M",
           str_detect(EM, "MSex") ~ "Sex-Specific M"),
         Prop = case_when(
           str_detect(EM, "PropAcr") ~ "Proportions Across",
           str_detect(EM, "PropWith_SR_Y") ~ "Proportions Within (SR_Y)",
           str_detect(EM, "PropWith_SR_ALY") ~ "Proportions Within (SR_ALY)"),
         Catch = case_when(
           str_detect(EM, "CatAgg") ~ "Aggregated Catch",
           str_detect(EM, "CatSex") ~ "Sex-Specific Catch")) %>% 
  drop_na()

# Plot Proportions across vs. within for parameters
pdf(here("figs", "Experiment 1", "RE_Param_Acr_vs_With.pdf"), width = 17, height = 10)
prop_param_df %>% 
  filter(Magg == "Sex-Specific M",
         Catch == "Aggregated Catch",
         Type %in% c("Bmsy", "Fmsy", "Tier 3 HCR Catch"),
         Prop %in% c("Proportions Within (SR_Y)", "Proportions Within (SR_ALY)")) %>% 
  ggplot(aes(x = RE, fill = Prop)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = 0, lty = 2, size = 1.3) + 
  labs(x = "Parameter", y = "Relative Error",
       title = "Sex-Specific M, Aggregated Catch") +
  theme_tj() +
  facet_grid(Type~OM, scales = "free") +
  theme(legend.position = "top")

prop_param_df %>% 
  filter(Magg == "Sex-Specific M",
         Catch == "Aggregated Catch",
         Type %in% c("Bmsy", "Fmsy", "Tier 3 HCR Catch"),
         Prop %in% c("Proportions Across", "Proportions Within (SR_Y)")) %>% 
  ggplot(aes(x = RE, fill = Prop)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = 0, lty = 2, size = 1.3) + 
  labs(x = "Parameter", y = "Relative Error",
       title = "Sex-Specific M, Aggregated Catch") +
  theme_tj() +
  facet_grid(Type~OM, scales = "free") +
  theme(legend.position = "top")

prop_param_df %>% 
  filter(Magg == "Sex-Specific M",
         Catch == "Sex-Specific Catch",
         Type %in% c("Bmsy", "Fmsy", "Tier 3 HCR Catch"),
         Prop %in% c("Proportions Within (SR_Y)", "Proportions Within (SR_ALY)")) %>% 
  ggplot(aes(x = RE, fill = Prop)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = 0, lty = 2, size = 1.3) + 
  labs(x = "Parameter", y = "Relative Error",
       title = "Sex-Specific M, Sex-Specific Catch") +
  theme_tj() +
  facet_grid(Type~OM, scales = "free") +
  theme(legend.position = "top")

prop_param_df %>% 
  filter(Magg == "Sex-Specific M",
         Catch == "Sex-Specific Catch",
         Type %in% c("Bmsy", "Fmsy", "Tier 3 HCR Catch"),
         Prop %in% c("Proportions Within (SR_Y)", "Proportions Across")) %>% 
  ggplot(aes(x = RE, fill = Prop)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = 0, lty = 2, size = 1.3) + 
  labs(x = "Parameter", y = "Relative Error",
       title = "Sex-Specific M, Sex-Specific Catch") +
  theme_tj() +
  facet_grid(Type~OM, scales = "free") +
  theme(legend.position = "top")

prop_param_df %>% 
  filter(Magg == "Aggregated M",
         Catch == "Aggregated Catch",
         Type %in% c("Bmsy", "Fmsy", "Tier 3 HCR Catch"),
         Prop %in% c("Proportions Within (SR_Y)", "Proportions Within (SR_ALY)")) %>% 
  ggplot(aes(x = RE, fill = Prop)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = 0, lty = 2, size = 1.3) + 
  labs(x = "Parameter", y = "Relative Error",
       title = "Aggregated M, Aggregated Catch") +
  theme_tj() +
  facet_grid(Type~OM, scales = "free") +
  theme(legend.position = "top")

prop_param_df %>% 
  filter(Magg == "Sex-Specific M",
         Catch == "Aggregated Catch",
         Type %in% c("Bmsy", "Fmsy", "Tier 3 HCR Catch"),
         Prop %in% c("Proportions Within (SR_Y)", "Proportions Across")) %>% 
  ggplot(aes(x = RE, fill = Prop)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = 0, lty = 2, size = 1.3) + 
  labs(x = "Parameter", y = "Relative Error",
       title = "Sex-Specific M, Aggregated Catch") +
  theme_tj() +
  facet_grid(Type~OM, scales = "free") +
  theme(legend.position = "top")

dev.off()


#### Sex-specific Catch comparison -------------------------------------------

# Plot Proportions across vs. within for parameters
pdf(here("figs", "Experiment 1", "RE_Param_SS_vs_AggCatch.pdf"), width = 17, height = 10)
prop_param_df %>% 
  filter(Magg == "Sex-Specific M",
         Type %in% c("Bmsy", "Fmsy", "Tier 3 HCR Catch"),
         Prop == "Proportions Across") %>% 
  ggplot(aes(x = RE, fill = Catch)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = 0, lty = 2, size = 1.3) + 
  labs(x = "Parameter", y = "Relative Error",
       title = "Sex-Specific M, Sex-Specific Catch (Proportions Across)") +
  theme_tj() +
  facet_grid(Type~OM) +
  theme(legend.position = "top") 

prop_param_df %>% 
  filter(Magg == "Aggregated M",
         Type %in% c("Bmsy", "Fmsy", "Tier 3 HCR Catch"),
         Prop == "Proportions Across") %>% 
  ggplot(aes(x = RE, fill = Catch)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = 0, lty = 2, size = 1.3) + 
  labs(x = "Parameter", y = "Relative Error",
       title = "Aggregated M, Sex-Specific Catch (Proportions Across)") +
  theme_tj() +
  facet_grid(Type~OM, scales = "free") +
  theme(legend.position = "top")

prop_param_df %>% 
  filter(Magg == "Aggregated M",
         Type %in% c("Bmsy", "Fmsy", "Tier 3 HCR Catch"),
         Prop == "Proportions Across",
         OM == "Growth_M (0,10)") %>% 
  ggplot(aes(x = RE, fill = Catch)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = 0, lty = 2, size = 1.3) + 
  labs(x = "Parameter", y = "Relative Error",
       title = "Aggregated M, Sex-Specific Catch (Proportions Across), OM = Growth_M (0,10)") +
  theme_tj() +
  facet_grid(Type~Catch, scales = "free") +
  theme(legend.position = "top") 
dev.off()

##### Other parameters --------------------------------------------------------
# summarize re for parameters
exp1_param_sum = exp1_param_df %>% 
  filter(Convergence == "Converged",
         !str_detect(Type, "Ratio")) %>% 
  mutate(RE = (as.numeric(Pred) - Truth) / Truth,
         Magg = case_when(
           str_detect(EM, "MAgg") ~ "Aggregated M",
           str_detect(EM, "MSex") ~ "Sex-Specific M"),
         Prop = case_when(
           str_detect(EM, "PropAcr") ~ "Proportions Across",
           str_detect(EM, "PropWith_SR_Y") ~ "Proportions Within (SR_Y)",
           str_detect(EM, "PropWith_SR_ALY") ~ "Proportions Within (SR_ALY)"),
         Catch = case_when(
           str_detect(EM, "CatAgg") ~ "Aggregated Catch",
           str_detect(EM, "CatSex") ~ "Sex-Specific Catch")) %>% 
  group_by(OM, EM, Type, Magg, Prop, Catch) %>% 
  summarize(Median = median(RE),
            lwr_95 = quantile(RE, 0.025),
            upr_95 = quantile(RE, 0.975))

# plot (all EMs)
pdf(here("figs", "Experiment 1", "RE_Param_allEMs.pdf"), width = 30, height = 14)
print(ggplot(exp1_param_sum %>% 
               filter(!str_detect(EM, "SR_ALY")), 
             aes(x = Type, y = Median, ymin = lwr_95, ymax = upr_95)) +
        geom_pointrange(position = position_dodge2(width = 0.65), 
                        size = 2, linewidth = 1) +       
        facet_grid(OM~EM, scales = "free") +
        geom_hline(yintercept = 0, lty = 2, size = 1.3) + 
        labs(x = "Parameter", y = "Relative Error") +
        theme_tj() +
        scale_x_discrete(guide = guide_axis(angle = 90)) +
        coord_cartesian(ylim = c(-1.5,1.5)) )

# Zoomed in version
print(ggplot(exp1_param_sum %>% 
               filter(!str_detect(EM, "SR_ALY")), aes(x = Type, y = Median, ymin = lwr_95, ymax = upr_95)) +
        geom_pointrange(position = position_dodge2(width = 0.65), 
                        size = 2, linewidth = 1) +       
        facet_grid(OM~EM, scales = "free") +
        geom_hline(yintercept = 0, lty = 2, size = 1.3) + 
        labs(x = "Parameter", y = "Relative Error", title = "Zoomed") +
        theme_tj() +
        scale_x_discrete(guide = guide_axis(angle = 90)) +
        coord_cartesian(ylim = c(-0.5,0.5)) )
dev.off()  

# plot (Age only EM)
pdf(here("figs", "Experiment 1", "RE_Param_Age.pdf"), width = 13)
print(ggplot(exp1_param_sum %>% filter(EM == "Age"), 
             aes(x = Type, y = Median, ymin = lwr_95, ymax = upr_95)) +
        geom_pointrange(position = position_dodge2(width = 0.65), 
                        size = 1, linewidth = 1) +       
        facet_wrap(~OM, scales = "free_y") +
        geom_hline(yintercept = 0, lty = 2, size = 1.3) + 
        labs(x = "Parameter", y = "Relative Error") +
        theme_tj() +
        scale_x_discrete(guide = guide_axis(angle = 90)) +
        coord_cartesian(ylim = c(-1,1)) )
print(ggplot(exp1_param_sum %>% filter(EM == "Age"), 
             aes(x = Type, y = Median, ymin = lwr_95, ymax = upr_95)) +
        geom_pointrange(position = position_dodge2(width = 0.65), 
                        size = 1, linewidth = 1) +       
        facet_wrap(~OM, scales = "free_y") +
        geom_hline(yintercept = 0, lty = 2, size = 1.3) + 
        labs(x = "Parameter", y = "Relative Error", title = "Zoomed") +
        theme_tj() +
        scale_x_discrete(guide = guide_axis(angle = 90)) +
        coord_cartesian(ylim = c(-0.3,0.3)) )
dev.off()  


# aggregated M EMs
pdf(here("figs", "Experiment 1", "RE_Param_Magg.pdf"), width = 16, height = 13)
print(ggplot(exp1_param_sum %>% 
               filter(!str_detect(EM, "SR_ALY"),
                      Magg == "Aggregated M"),
             aes(x = Type, y = Median, ymin = lwr_95, ymax = upr_95)) +
        geom_pointrange(position = position_dodge2(width = 0.65), 
                        size = 1, linewidth = 1) +       
        facet_grid(OM~EM, scales = "free") +
        geom_hline(yintercept = 0, lty = 2, size = 1.3) + 
        labs(x = "Parameter", y = "Relative Error", title = "Aggregated M EMs") +
        theme_tj() +
        scale_x_discrete(guide = guide_axis(angle = 90)) )  

print(ggplot(exp1_param_sum %>% 
               filter(!str_detect(EM, "SR_ALY"),
                      Magg == "Aggregated M"),
             aes(x = Type, y = Median, ymin = lwr_95, ymax = upr_95)) +
        geom_pointrange(position = position_dodge2(width = 0.65), 
                        size = 1, linewidth = 1) +       
        facet_grid(OM~EM, scales = "free") +
        geom_hline(yintercept = 0, lty = 2, size = 1.3) + 
        labs(x = "Parameter", y = "Relative Error", title = "Aggregated M EMs, Zoomed") +
        theme_tj() +
        scale_x_discrete(guide = guide_axis(angle = 90)) +
        coord_cartesian(ylim = c(-0.85,0.3))) 

print(ggplot(exp1_param_sum %>% 
               filter(!str_detect(EM, "SR_ALY"),
                      Prop == "Proportions Across",
                      Magg == "Aggregated M"),
             aes(x = Type, y = Median, ymin = lwr_95, ymax = upr_95, color = Catch)) +
        geom_pointrange(position = position_dodge2(width = 0.65), 
                        size = 1, linewidth = 1) +       
        facet_wrap(~OM) +
        geom_hline(yintercept = 0, lty = 2, size = 1.3) + 
        labs(x = "Parameter", y = "Relative Error", title = "Aggregated M, Proportions Across, Zoomed") +
        theme_tj() +
        coord_cartesian(ylim = c(-3, 3)) +
        scale_x_discrete(guide = guide_axis(angle = 90)) +
        theme(legend.position = "top")) 
print(ggplot(exp1_param_sum %>% 
               filter(!str_detect(EM, "SR_ALY"),
                      Prop == "Proportions Across",
                      Magg == "Aggregated M",
                      Type %in% c("M_F", "M_M", "R0")),
             aes(x = Type, y = Median, ymin = lwr_95, ymax = upr_95, color = Catch)) +
        geom_pointrange(position = position_dodge2(width = 0.65), 
                        size = 1, linewidth = 1) +       
        facet_wrap(~OM) +
        geom_hline(yintercept = 0, lty = 2, size = 1.3) + 
        labs(x = "Parameter", y = "Relative Error", title = "Aggregated M, Proportions Across, Zoomed") +
        theme_tj() +
        coord_cartesian(ylim = c(-0.8, 0)) +
        scale_x_discrete(guide = guide_axis(angle = 90)) +
        theme(legend.position = "top")) 

dev.off()


### Time Series Summary -----------------------------------------------------
###### Proportions within vs. across -------------------------------------------

prop_exp1 = exp1_ts_df %>% 
  filter(Convergence == "Converged") %>% 
  mutate(RE = (as.numeric(Pred)-Truth)/Truth,
         Magg = case_when(
           str_detect(EM, "MAgg") ~ "Aggregated M",
           str_detect(EM, "MSex") ~ "Sex-Specific M"),
         Prop = case_when(
           str_detect(EM, "PropAcr") ~ "Proportions Across",
           str_detect(EM, "PropWith_SR_Y") ~ "Proportions Within (SR_Y)",
           str_detect(EM, "PropWith_SR_ALY") ~ "Proportions Within (SR_ALY)"),
         Catch = case_when(
           str_detect(EM, "CatAgg") ~ "Aggregated Catch",
           str_detect(EM, "CatSex") ~ "Sex-Specific Catch")) %>% 
  drop_na()


# Plot EM across vs. within
pdf(here("figs", "Experiment 1", "RE_TS_PropAcr_vs_With.pdf"), width = 15, height = 10)
print(
  ggplot(prop_exp1 %>% 
           filter(Magg == "Sex-Specific M",
                  Catch == "Aggregated Catch",
                  Type != "Total Recruitment",
                  Prop %in% c("Proportions Within (SR_Y)", "Proportions Within (SR_ALY)")), 
         aes(x = RE, fill = Prop)) +
    geom_density(alpha = 0.5) +
    geom_vline(xintercept = 0, lty = 2) +
    facet_grid(Type~OM, scales = "free") +
    labs(x = "RE", y = "Relative Error",
         title = "Sex-Specific M, Aggregated Catch") +
    theme_tj() +
    theme(legend.position = "top")
)

print(
  ggplot(prop_exp1 %>% 
           filter(Magg == "Sex-Specific M",
                  Catch == "Aggregated Catch",
                  Type != "Total Recruitment",
                  Prop %in% c("Proportions Within (SR_Y)", "Proportions Across")), 
         aes(x = RE, fill = Prop)) +
    geom_density(alpha = 0.5) +
    geom_vline(xintercept = 0, lty = 2) +
    facet_grid(Type~OM, scales = "free") +
    labs(x = "RE", y = "Relative Error",
         title = "Sex-Specific M, Aggregated Catch") +
    theme_tj() +
    theme(legend.position = "top")
)

print(
  ggplot(prop_exp1 %>% 
           filter(Magg == "Sex-Specific M",
                  Catch == "Sex-Specific Catch",
                  Type != "Total Recruitment",
                  Prop %in% c("Proportions Within (SR_Y)", "Proportions Within (SR_ALY)")), 
         aes(x = RE, fill = Prop)) +
    geom_density(alpha = 0.5) +
    geom_vline(xintercept = 0, lty = 2) +
    facet_grid(Type~OM, scales = "free") +
    labs(x = "RE", y = "Relative Error",
         title = "Sex-Specific M, Sex-Specific Catch") +
    theme_tj() +
    theme(legend.position = "top")
)

print(
  ggplot(prop_exp1 %>% 
           filter(Magg == "Sex-Specific M",
                  Catch == "Sex-Specific Catch",
                  Type != "Total Recruitment",
                  Prop %in% c("Proportions Within (SR_Y)", "Proportions Across")), 
         aes(x = RE, fill = Prop)) +
    geom_density(alpha = 0.5) +
    geom_vline(xintercept = 0, lty = 2) +
    facet_grid(Type~OM, scales = "free") +
    labs(x = "RE", y = "Relative Error",
         title = "Sex-Specific M, Sex-Specific Catch") +
    theme_tj() +
    theme(legend.position = "top")
)

print(
  ggplot(prop_exp1 %>% 
           filter(Magg == "Aggregated M",
                  Catch == "Aggregated Catch",
                  Type != "Total Recruitment",
                  Prop %in% c("Proportions Within (SR_Y)", "Proportions Within (SR_ALY)")), 
         aes(x = RE, fill = Prop)) +
    geom_density(alpha = 0.5) +
    geom_vline(xintercept = 0, lty = 2) +
    facet_grid(Type~OM, scales = "free") +
    labs(x = "RE", y = "Relative Error",
         title = "Aggregated M, Aggregated Catch") +
    theme_tj() +
    theme(legend.position = "top")
)

print(
  ggplot(prop_exp1 %>% 
           filter(Magg == "Aggregated M",
                  Catch == "Aggregated Catch",
                  Type != "Total Recruitment",
                  Prop %in% c("Proportions Within (SR_Y)", "Proportions Across")), 
         aes(x = RE, fill = Prop)) +
    geom_density(alpha = 0.5) +
    geom_vline(xintercept = 0, lty = 2) +
    facet_grid(Type~OM, scales = "free") +
    labs(x = "RE", y = "Relative Error",
         title = "Aggregated M, Aggregated Catch") +
    theme_tj() +
    theme(legend.position = "top")
)

print(
  ggplot(prop_exp1 %>% 
           filter(Magg == "Aggregated M",
                  Catch == "Sex-Specific Catch",
                  Type != "Total Recruitment",
                  Prop %in% c("Proportions Within (SR_Y)", "Proportions Within (SR_ALY)")), 
         aes(x = RE, fill = Prop)) +
    geom_density(alpha = 0.5) +
    geom_vline(xintercept = 0, lty = 2) +
    facet_grid(Type~OM, scales = "free") +
    labs(x = "RE", y = "Relative Error",
         title = "Aggregated M, Sex-Specific Catch") +
    theme_tj() +
    theme(legend.position = "top")
)

print(
  ggplot(prop_exp1 %>% 
           filter(Magg == "Aggregated M",
                  Catch == "Sex-Specific Catch",
                  Type != "Total Recruitment",
                  Prop %in% c("Proportions Within (SR_Y)", "Proportions Across")), 
         aes(x = RE, fill = Prop)) +
    geom_density(alpha = 0.5) +
    geom_vline(xintercept = 0, lty = 2) +
    facet_grid(Type~OM, scales = "free") +
    labs(x = "RE", y = "Relative Error",
         title = "Aggregated M, Sex-Specific Catch") +
    theme_tj() +
    theme(legend.position = "top")
)

dev.off()

###### Relative error all EMs --------------------------------------------------
# time series summary in relative error
ts_exp1_sum = exp1_ts_df %>% 
  filter(Convergence == "Converged") %>%
  mutate(RE = (as.numeric(Pred)-Truth)/Truth,
         Magg = case_when(
           str_detect(EM, "MAgg") ~ "Aggregated M",
           str_detect(EM, "MSex") ~ "Sex-Specific M"),
         Prop = case_when(
           str_detect(EM, "PropAcr") ~ "Proportions Across",
           str_detect(EM, "PropWith_SR_Y") ~ "Proportions Within (SR_Y)",
           str_detect(EM, "PropWith_SR_ALY") ~ "Proportions Within (SR_ALY)"),
         Catch = case_when(
           str_detect(EM, "CatAgg") ~ "Aggregated Catch",
           str_detect(EM, "CatSex") ~ "Sex-Specific Catch")) %>% 
  group_by(Years, Type, OM, EM, Magg, Prop, Catch) %>% 
  summarize(Median = median(RE),
            lwr_95 = quantile(RE, 0.025),
            upr_95 = quantile(RE, 0.95))

# Plot time series EMs
pdf(here("figs", "Experiment 1", "RE_TS_AllEMs.pdf"), width = 34, height = 20)
print(
  ggplot(ts_exp1_sum %>% filter(Type == "Total Biomass",
                                !str_detect(EM, "ALY")), 
         aes(x = Years, y = Median, ymin = lwr_95, ymax = upr_95,
             fill = Catch, color = Catch)) +
    geom_line(size = 2) +
    geom_ribbon(alpha = 0.5) +
    facet_grid(OM~EM, scales = "free") +
    geom_hline(yintercept = 0, lty = 2, size = 1.3) + 
    labs(x = "Year", y = "Relative Error in Total Biomass") +
    theme_tj() +
    coord_cartesian(ylim = c(-0.85,0.85)) 
)

print(
  ggplot(ts_exp1_sum %>% filter(Type == "Spawning Stock Biomass",
                                !str_detect(EM, "ALY")), 
         aes(x = Years, y = Median, ymin = lwr_95, ymax = upr_95,
             fill = Catch, color = Catch)) +
    geom_line(size = 2) +
    geom_ribbon(alpha = 0.5) +
    facet_grid(OM~EM, scales = "free") +
    geom_hline(yintercept = 0, lty = 2, size = 1.3) + 
    labs(x = "Year", y = "Relative Error in Spawning Stock Biomass") +
    theme_tj() +
    coord_cartesian(ylim = c(-0.85,0.85)) 
)

print(
  ggplot(ts_exp1_sum %>% filter(Type == "Total Fishing Mortality",
                                !str_detect(EM, "ALY")), 
         aes(x = Years, y = Median, ymin = lwr_95, ymax = upr_95,
             fill = Catch, color = Catch)) +
    geom_line(size = 2) +
    geom_ribbon(alpha = 0.5) +
    facet_grid(OM~EM, scales = "free") +
    geom_hline(yintercept = 0, lty = 2, size = 1.3) + 
    labs(x = "Year", y = "Relative Error in Total Fishing Mortality") +
    theme_tj() +
    coord_cartesian(ylim = c(-1,1)) 
)

print(
  ggplot(ts_exp1_sum %>% filter(Type == "Total Recruitment",
                                !str_detect(EM, "ALY")), 
         aes(x = Years, y = Median, ymin = lwr_95, ymax = upr_95,
             fill = Catch, color = Catch)) +
    geom_line(size = 2) +
    geom_ribbon(alpha = 0.5) +
    facet_grid(OM~EM, scales = "free") +
    geom_hline(yintercept = 0, lty = 2, size = 1.3) + 
    labs(x = "Year", y = "Relative Error in Total Recruitment") +
    theme_tj() +
    coord_cartesian(ylim = c(-0.85,0.85)) 
)
dev.off()

# Aggregated m ems
pdf(here("figs", "Experiment 1", "RE_TS_MaggEMs.pdf"), width = 20, height = 20)
print(
  ggplot(ts_exp1_sum %>% filter(Magg == "Aggregated M",
                                Prop == "Proportions Across"), 
         aes(x = Years, y = Median, ymin = lwr_95, ymax = upr_95,
             fill = Catch, color = Catch)) +
    geom_line(size = 1.3) +
    geom_ribbon(alpha = 0.5) +
    facet_grid(OM~Type, scales = "free") +
    geom_hline(yintercept = 0, lty = 2, size = 1.3) + 
    labs(x = "Year", y = "Relative Error", title = "Aggregated M, Proportions Across") +
    theme_tj() +
    coord_cartesian(ylim = c(-0.85,0.85)) +
    theme(legend.position = "top")
)

dev.off()

# Age only EM
pdf(here("figs", "Experiment 1", "RE_TS_AgeEM.pdf"), width = 15, height = 13)
print(
  ggplot(ts_exp1_sum %>% filter(EM == "Age"), 
         aes(x = Years, y = Median, ymin = lwr_95, ymax = upr_95)) +
    geom_line(size = 1.3) +
    geom_ribbon(alpha = 0.5) +
    facet_grid(Type~OM, scales = "free") +
    geom_hline(yintercept = 0, lty = 2, size = 1.3) + 
    labs(x = "Year", y = "Relative Error") +
    theme_tj() +
    coord_cartesian(ylim = c(-0.85,0.85)) 
)
dev.off()

pdf(here("figs", "Experiment 1", "RE_TS_AgeEM_ZoomBiom.pdf"), width = 18, height = 6.5)
print(
  ggplot(ts_exp1_sum %>% filter(EM == "Age",Type == "Spawning Stock Biomass"), 
         aes(x = Years, y = Median, ymin = lwr_95, ymax = upr_95)) +
    geom_line(size = 1.3) +
    geom_ribbon(alpha = 0.5) +
    facet_wrap(~OM, scales = "free", nrow = 1) +
    geom_hline(yintercept = 0, lty = 2, size = 1.3) + 
    labs(x = "Year", y = "Relative Error in Spawning Stock Biomass") +
    theme_tj() 
)
print(
  ggplot(ts_exp1_sum %>% filter(EM == "Age",Type == "Total Biomass"), 
         aes(x = Years, y = Median, ymin = lwr_95, ymax = upr_95)) +
    geom_line(size = 1.3) +
    geom_ribbon(alpha = 0.5) +
    facet_wrap(~OM, scales = "free", nrow = 1) +
    geom_hline(yintercept = 0, lty = 2, size = 1.3) + 
    labs(x = "Year", y = "Relative Error in Total Biomass") +
    theme_tj() 
)
dev.off()



# Numbers at age ----------------------------------------------------------

# Get summary of relative error and simulation intervals, etc.
exp1_naa_sum = exp1_naa_df %>% 
  filter(convergence == "Converged") %>% 
  mutate(Pred = as.numeric(Pred),
         Truth = as.numeric(Truth),
         RE = (Pred - Truth) / Truth) %>% 
  group_by(Year, Age, Sex, EM, OM) %>% 
  summarize(Median = median(RE),
            lwr_95 = quantile(RE, 0.025),
            upr_95 = quantile(RE, 0.975))

# unique ems
unique_ems = unique(exp1_naa_sum$EM)
pdf(here("figs", "Experiment 1", "RE_NAA.pdf"), width = 15, height = 13)
for(i in 1:length(unique_ems)) {
  print(
    ggplot(exp1_naa_sum %>% 
             filter(Age %in% c(seq(3, 30, 4)), 
                    EM == unique_ems[i]),
           aes(x = Year, y = Median, ymin = lwr_95, ymax = upr_95, 
               color = factor(Sex), fill = factor(Sex))) +
      geom_line() +
      geom_hline(yintercept = 0, lty = 2) +
      geom_ribbon(alpha = 0.35) +
      facet_grid(OM~Age) +
      coord_cartesian(ylim = c(-1, 1)) +
      theme_tj() +
      labs(x = "Year", y = "RE in NAA", title = unique_ems[i])
  )
} # end i
dev.off()

# Experiment 2 ------------------------------------------------------------

# Read in files
exp2_selex_df = data.table::fread(here("output", "Experiment_2_Selex.csv")) %>% filter(!str_detect(EM, "Magg"))
exp2_growth_df = data.table::fread(here("output", "Experiment_2_Growth.csv")) %>% filter(!str_detect(EM, "Magg"))
exp2_param_df = data.table::fread(here("output", "Experiment_2_Param.csv")) %>% filter(!str_detect(EM, "Magg"))
exp2_ts_df = data.table::fread(here("output", "Experiment_2_TimeSeries.csv")) %>% filter(!str_detect(EM, "Magg"))
exp2_sr_df = data.table::fread(here("output", "Experiment_2_SexRatio.csv")) %>% filter(!str_detect(EM, "Magg"))

##### Convergence summary -----------------------------------------------------
# Convergence summary
conv_df = exp2_growth_df %>% 
  filter(Convergence == "Converged", Age == 1, Sex == 1) %>% 
  group_by(OM, EM) %>% 
  summarize(sum = n())

# Plot convergence
pdf(here("figs", "Experiment 2", "Convergence.pdf"), width = 15)
ggplot(conv_df, aes(x = OM, y = sum)) +
  geom_point(size = 3) +
  facet_wrap(~EM, scales = "free_y") +
  theme_tj() +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  labs(x = "Operating Models", y = "Convergence") 
dev.off()


### Selectivity Summary -----------------------------------------------------
# summarize relative error
exp2_selex_sum = exp2_selex_df %>% 
  filter(Convergence == "Converged") %>% 
  mutate(RE = (Pred - True) / True) %>% 
  group_by(OM, EM, Age, Sex, Type) %>% 
  summarize(Median = median(RE),
            lwr_95 = quantile(RE, 0.025),
            upr_95 = quantile(RE, 0.975))

# Fishery Selectivity with Age EM
pdf(here("figs", "Experiment 2", "RE_Fish_Selex.pdf"), height = 14, width = 30)
fish_selex = print(
  ggplot(exp2_selex_sum %>% filter(Type == "Fishery Selectivity") %>% 
           mutate(Sex = ifelse(Sex == 1, "Female", "Male")), 
         aes(x = Age, y = Median, ymin = lwr_95, ymax = upr_95,
             fill = factor(Sex), color = factor(Sex), group = factor(Sex))) +
    geom_line(size = 3) +
    geom_ribbon(alpha = 0.35) +
    geom_hline(yintercept = 0, lty = 2, size = 1.3) + 
    labs(x = "Age", y = "Relative Error in Fishery Selectivity", color = "Sex", fill = "Sex") +
    facet_grid(EM~OM, scales = "free") +
    theme_tj() +
    theme(legend.position = "top")
)
dev.off()

# Fishery Selectivity with Age EM
pdf(here("figs", "Experiment 2", "RE_Srv_Selex.pdf"),height = 14, width = 30)
srv_selex = print(
  ggplot(exp2_selex_sum %>% filter(Type == "Survey Selectivity") %>% 
           mutate(Sex = ifelse(Sex == 1, "Female", "Male")), 
         aes(x = Age, y = Median, ymin = lwr_95, ymax = upr_95,
             fill = factor(Sex), color = factor(Sex), group = factor(Sex))) +
    geom_line(size = 3) +
    geom_ribbon(alpha = 0.35) +
    geom_hline(yintercept = 0, lty = 2, size = 1.3) + 
    labs(x = "Age", y = "Relative Error in Survey Selectivity", color = "Sex", fill = "Sex") +
    facet_grid(EM~OM, scales = "free") +
    theme_tj() +
    theme(legend.position = "top")
)
dev.off()


### Growth Summary ----------------------------------------------------------

# summarize relative error
exp2_growth_sum = exp2_growth_df %>% 
  filter(Convergence == "Converged") %>% 
  mutate(RE = (Pred - True) / True) %>% 
  group_by(OM, EM, Age, Sex) %>% 
  summarize(Median = median(RE),
            lwr_95 = quantile(RE, 0.025),
            upr_95 = quantile(RE, 0.975))

pdf(here("figs", "Experiment 2", "RE_Growth.pdf"), height = 14, width = 30)
print(
  ggplot(exp2_growth_sum %>%  mutate(Sex = ifelse(Sex == 1, "Female", "Male")), 
         aes(x = Age, y = Median, ymin = lwr_95, ymax = upr_95,
             fill = factor(Sex), color = factor(Sex), group = factor(Sex))) +
    geom_line(size = 3) +
    geom_ribbon(alpha = 0.35) +
    geom_hline(yintercept = 0, lty = 2, size = 1.3) + 
    labs(x = "Age", y = "Relative Error in Growth (WAA)", color = "Sex", fill = "Sex") +
    facet_grid(EM~OM, scales = "free") +
    theme_tj() +
    theme(legend.position = "top")
)
dev.off()

### Parameter Summary -------------------------------------------------------
# summarize re
exp2_param_sum = exp2_param_df %>% 
  filter(Convergence == "Converged") %>% 
  mutate(RE = (as.numeric(Pred) - Truth) / Truth) %>% 
  group_by(OM, EM, Type) %>% 
  summarize(Median = median(RE),
            lwr_95 = quantile(RE, 0.025),
            upr_95 = quantile(RE, 0.975))


# Sex ratio estimability
pdf(here("figs", "Experiment 2", "RE_SexRatio.pdf"), width = 15, height = 10)
exp2_param_df %>% 
  filter(Convergence == "Converged",
         Type == "Female Sex Ratio",
         EM %in% c("Est_PropWith_SR_Y", "Est_PropWith_SR_ALY")) %>% 
       mutate(RE = (as.numeric(Pred) - Truth) / Truth) %>% 
  ggplot(aes(x = Type, y = RE, fill = EM)) +
  geom_violin(alpha = 0.5) +
  geom_boxplot(width = 0.1, outlier.color = NA,
               position = position_dodge(width = 0.9)) +
  geom_hline(yintercept = 0, lty = 2, size = 1.3) + 
  labs(x = "Parameter", y = "Relative Error") +
  theme_tj() +
  facet_wrap(~OM, scales = "free_y") +
  theme(legend.position = "top")

# Comparing across vs. within for a given year
exp2_param_df %>% 
  filter(Convergence == "Converged",
         Type == "Female Sex Ratio",
         str_detect(Type, "Sex"),
         EM %in% c("Est_PropAcr", "Est_PropWith_SR_Y")) %>% 
  mutate(RE = (as.numeric(Pred) - Truth) / Truth) %>% 
  ggplot(aes(x = Type, y = RE, fill = EM)) +
  geom_violin(alpha = 0.5) +
  geom_boxplot(width = 0.1, outlier.color = NA,
               position = position_dodge(width = 0.9)) +
  geom_hline(yintercept = 0, lty = 2, size = 1.3) + 
  labs(x = "Parameter", y = "Relative Error") +
  theme_tj() +
  facet_wrap(~OM, scales = "free_y") +
  theme(legend.position = "top")
dev.off()
  
# plot all other parameters and EMs
pdf(here("figs", "Experiment 2", "RE_ParamAllEMs.pdf"), width = 15, height = 10)
# Comparing proportions within variants
print(ggplot(exp2_param_sum %>% 
               filter(EM %in% c("Est_PropWith_SR_Y", "Est_PropWith_SR_ALY"),
                      Type != "Male Sex Ratio"), 
             aes(x = Type, y = Median, ymin = lwr_95, ymax = upr_95, color = EM)) +
        geom_pointrange(position = position_dodge2(width = 0.65), 
                        size = 1, linewidth = 1) +
        facet_wrap(~OM, scales = "free_y") +
        geom_hline(yintercept = 0, lty = 2, size = 1.3) + 
        labs(x = "Parameter", y = "Relative Error") +
        theme_tj() +
        theme(legend.position = "top") +
        scale_x_discrete(guide = guide_axis(angle = 90))  )

# comparing best proprotions within variant
print(ggplot(exp2_param_sum %>% 
               filter(EM != "Est_PropWith_SR_ALY", Type != "Male Sex Ratio", !str_detect(OM, "No")), 
             aes(x = Type, y = Median, ymin = lwr_95, ymax = upr_95, color = EM)) +
        geom_pointrange(position = position_dodge2(width = 0.65), 
                        size = 1, linewidth = 1) +
        facet_wrap(~OM, scales = "free_y") +
        geom_hline(yintercept = 0, lty = 2, size = 1.3) + 
        labs(x = "Parameter", y = "Relative Error") +
        theme_tj() +
        theme(legend.position = "top") +
        scale_x_discrete(guide = guide_axis(angle = 90)) +
        coord_cartesian(ylim = c(-1, 1)))

# comparing EMs that are fixed for no diff OMs
print(ggplot(exp2_param_sum %>% 
               filter(Type != "Male Sex Ratio", str_detect(OM, "No")), 
             aes(x = Type, y = Median, ymin = lwr_95, ymax = upr_95, color = EM)) +
        geom_pointrange(position = position_dodge2(width = 0.65), 
                        size = 1, linewidth = 1) +
        facet_wrap(~OM, scales = "free_y") +
        geom_hline(yintercept = 0, lty = 2, size = 1.3) + 
        labs(x = "Parameter", y = "Relative Error") +
        theme_tj() +
        theme(legend.position = "top") +
        scale_x_discrete(guide = guide_axis(angle = 90))  )
dev.off()  


### Time Series Summary -----------------------------------------------------
# time series summary
ts_exp2_sum = exp2_ts_df %>% 
  filter(Convergence == "Converged") %>% 
  mutate(RE = (as.numeric(Pred)-Truth)/Truth) %>% 
  group_by(Years, Type, OM, EM) %>% 
  summarize(Median = median(RE),
            lwr_95 = quantile(RE, 0.025),
            upr_95 = quantile(RE, 0.975))

# Density plot of across vs.within
pdf(here("figs", 'Experiment 2', "RE_TSDensity_Acr_vs_With.pdf"), width = 15, height = 12)
# comparing within variants
exp2_ts_df %>% 
  filter(Convergence == "Converged",
         EM %in% c("Est_PropWith_SR_Y", "Est_PropWith_SR_ALY"),
         Years == 1) %>% 
  mutate(RE = (as.numeric(Pred)-Truth)/Truth) %>% 
  ggplot(aes(x = RE, fill = EM)) +
  geom_density(alpha = 0.5) +
  facet_grid(OM~Type, scales = "free") +
  geom_vline(xintercept = 0, lty = 2) +
  labs(x = "Year", y = "Relative Error") +
  theme_tj() +
  theme(legend.position = "top") 

exp2_ts_df %>% 
  filter(Convergence == "Converged",
         EM %in% c("Est_PropWith_SR_Y", "Est_PropAcr"),
         Years == 1) %>% 
  mutate(RE = (as.numeric(Pred)-Truth)/Truth) %>% 
  ggplot(aes(x = RE, fill = EM)) +
  geom_density(alpha = 0.5) +
  facet_grid(OM~Type, scales = "free") +
  geom_vline(xintercept = 0, lty = 2) +
  labs(x = "Year", y = "Relative Error") +
  theme_tj() +
  theme(legend.position = "top") 
dev.off()

# plot RE as a time series
pdf(here("figs", 'Experiment 2', "RE_TS_Acr_vs_With.pdf"), width = 15)
# comparing within variants
print(
  ggplot(ts_exp2_sum %>% filter(Type == "Total Biomass",
         EM %in% c("Est_PropWith_SR_Y", "Est_PropWith_SR_ALY")),
         aes(x = Years, y = Median, ymin = lwr_95, ymax = upr_95, fill = EM)) +
    geom_line(size = 1.3) +
    geom_ribbon(alpha = 0.5) +
    facet_grid(~OM, scales = "free") +
    geom_hline(yintercept = 0, lty = 2, size = 1.3) + 
    labs(x = "Year", y = "Relative Error in Total Biomass") +
    theme_tj() +
    theme(legend.position = "top") 
)

print(
  ggplot(ts_exp2_sum %>% filter(Type == "Total Biomass",
                                EM %in% c("Est_PropWith_SR_Y", "Est_PropAcr")),
         aes(x = Years, y = Median, ymin = lwr_95, ymax = upr_95, fill = EM)) +
    geom_line(size = 1.3) +
    geom_ribbon(alpha = 0.5) +
    facet_grid(~OM, scales = "free") +
    geom_hline(yintercept = 0, lty = 2, size = 1.3) + 
    labs(x = "Year", y = "Relative Error in Total Biomass") +
    theme_tj() +
    theme(legend.position = "top") 
)

print(
  ggplot(ts_exp2_sum %>% filter(Type == "Spawning Stock Biomass", 
                                EM %in% c("Est_PropWith_SR_Y", "Est_PropWith_SR_ALY")),
         aes(x = Years, y = Median, ymin = lwr_95, ymax = upr_95, fill = EM)) +
    geom_line(size = 1.3) +
    geom_ribbon(alpha = 0.5) +
    facet_grid(~OM, scales = "free") +
    geom_hline(yintercept = 0, lty = 2, size = 1.3) + 
    labs(x = "Year", y = "Relative Error in Spawning Stock Biomass") +
    theme_tj() +
    theme(legend.position = "top") 
  
)

print(
  ggplot(ts_exp2_sum %>% filter(Type == "Spawning Stock Biomass", 
                                EM %in% c("Est_PropWith_SR_Y", "Est_PropAcr")),
         aes(x = Years, y = Median, ymin = lwr_95, ymax = upr_95, fill = EM)) +
    geom_line(size = 1.3) +
    geom_ribbon(alpha = 0.5) +
    facet_grid(~OM, scales = "free") +
    geom_hline(yintercept = 0, lty = 2, size = 1.3) + 
    labs(x = "Year", y = "Relative Error in Spawning Stock Biomass") +
    theme_tj() +
    theme(legend.position = "top") 
)

print(ggplot(ts_exp2_sum %>% filter(Type == "Total Fishing Mortality",
              EM %in% c("Est_PropWith_SR_Y", "Est_PropWith_SR_ALY")),
             aes(x = Years, y = Median, ymin = lwr_95, ymax = upr_95, fill = EM)) +
        geom_line(size = 1.3) +
        geom_ribbon(alpha = 0.5) +
        facet_grid(~OM, scales = "free") +
        geom_hline(yintercept = 0, lty = 2, size = 1.3) + 
        labs(x = "Year", y = "Relative Error in Total Fishing Mortality") +
        theme_tj() +
        theme(legend.position = "top") )

print(ggplot(ts_exp2_sum %>% filter(Type == "Total Fishing Mortality",
                                    EM %in% c("Est_PropWith_SR_Y", "Est_PropAcr")),
             aes(x = Years, y = Median, ymin = lwr_95, ymax = upr_95, fill = EM)) +
        geom_line(size = 1.3) +
        geom_ribbon(alpha = 0.5) +
        facet_grid(~OM, scales = "free") +
        geom_hline(yintercept = 0, lty = 2, size = 1.3) + 
        labs(x = "Year", y = "Relative Error in Total Fishing Mortality") +
        theme_tj() +
        theme(legend.position = "top") )

print(
  ggplot(ts_exp2_sum %>% filter(Type == "Total Recruitment",
                                EM %in% c("Est_PropWith_SR_Y", "Est_PropWith_SR_ALY")),
         aes(x = Years, y = Median, ymin = lwr_95, ymax = upr_95, fill = EM)) +
    geom_line(size = 1.3) +
    geom_ribbon(alpha = 0.5) +
    facet_grid(~OM, scales = "free") +
    geom_hline(yintercept = 0, lty = 2, size = 1.3) + 
    labs(x = "Year", y = "Relative Error in Total Recruitment") +
    theme_tj() +
    theme(legend.position = "top") 
)

print(
  ggplot(ts_exp2_sum %>% filter(Type == "Total Recruitment",
                                EM %in% c("Est_PropWith_SR_Y", "Est_PropAcr")),
         aes(x = Years, y = Median, ymin = lwr_95, ymax = upr_95, fill = EM)) +
    geom_line(size = 1.3) +
    geom_ribbon(alpha = 0.5) +
    facet_grid(~OM, scales = "free") +
    geom_hline(yintercept = 0, lty = 2, size = 1.3) + 
    labs(x = "Year", y = "Relative Error in Total Recruitment") +
    theme_tj() +
    theme(legend.position = "top") 
)

dev.off()


# Only plot RE time series for Fixed EM
pdf(here("figs", 'Experiment 2', "RE_TS_FixEM.pdf"), width = 13, height = 15)

print(
  ggplot(ts_exp2_sum %>% filter(EM == "Fix", !str_detect(OM, "No")), 
         aes(x = Years, y = Median, ymin = lwr_95, ymax = upr_95)) +
    geom_line(size = 1.3) +
    geom_ribbon(alpha = 0.5) +
    facet_grid(Type~OM, scales = "free") +
    geom_hline(yintercept = 0, lty = 2, size = 1.3) + 
    labs(x = "Year", y = "Relative Error") +
    theme_tj() +
    theme(legend.position = "top") 
)

ggplot(ts_exp2_sum %>% filter(EM == "Fix", !str_detect(OM, "No")), 
       aes(x = Years, y = Median, ymin = lwr_95, ymax = upr_95)) +
  geom_line(size = 1.3) +
  geom_ribbon(alpha = 0.5) +
  facet_grid(Type~OM, scales = "free") +
  geom_hline(yintercept = 0, lty = 2, size = 1.3) + 
  labs(x = "Year", y = "Relative Error", title = "Zoom") +
  theme_tj() +
  theme(legend.position = "top") +
  coord_cartesian(ylim = c(-1.5, 1.5))

dev.off()

