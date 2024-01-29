# Purpose: To create summary plots for simulations
# Creator: Matthew LH. Cheng
# date: 1/26/23

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
exp_1_om_levs = c("Across_20", "Within_20",
                  "Across_50", "Within_50",
                  "Across_100", "Within_100",
                  "Across_150", "Within_150")
# Read in files
exp1_conv = data.table::fread(here("output", "Experiment_1_Convergence.csv")) %>% 
  mutate(OM = factor(OM, levels = exp_1_om_levs))
exp1_selex_df = data.table::fread(here("output", "Experiment_1_Selex.csv"))  %>% 
  mutate(OM = factor(OM, levels = exp_1_om_levs))
exp1_growth_df = data.table::fread(here("output", "Experiment_1_Growth.csv"))  %>% 
  mutate(OM = factor(OM, levels = exp_1_om_levs))
exp1_param_df = data.table::fread(here("output", "Experiment_1_Param.csv"))  %>% 
  mutate(OM = factor(OM, levels = exp_1_om_levs))
exp1_ts_df = data.table::fread(here("output", "Experiment_1_TimeSeries.csv"))  %>% 
  mutate(OM = factor(OM, levels = exp_1_om_levs))
exp1_cov_df = data.table::fread(here("output", "Experiment_1_Coverage.csv"))  %>% 
  mutate(OM = factor(OM, levels = exp_1_om_levs))

##### Convergence summary -----------------------------------------------------
exp1_conv = exp1_conv %>% 
  mutate(convergence = 
           case_when( # munging gradient stuff
             (pdHess == TRUE & 
                gradient <= 0.001 &
                sdNA == FALSE &
                convergence == "Not Converged") ~ "Converged",
             (pdHess == TRUE & 
                gradient <= 0.001 &
                sdNA == FALSE &
                convergence == "Converged") ~ "Converged",
             (pdHess == FALSE |
                gradient > 0.001 |
                sdNA == TRUE) ~ "Not Converged"
           ))

# Convergence summary
conv_df = exp1_conv %>% 
  filter(convergence == "Converged") %>% 
  group_by(OM, EM) %>% 
  summarize(sum = n())

# Plot convergence
pdf(here("figs", "Experiment 1", "Convergence.pdf"), height = 10, width = 15)
ggplot(conv_df, aes(x = OM, y = sum/750, group = EM, color = EM)) +
  geom_point(size = 5) +
  geom_line(size = 1.3, alpha = 0.5) +
  theme_tj() +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  labs(x = "Operating Models", y = "Convergence Rate") +
  theme(legend.position = "top")
dev.off()


### Selectivity Summary -----------------------------------------------------
# summarize relative error
exp1_selex_sum = exp1_selex_df %>% 
  select(-Convergence) %>% 
  left_join(exp1_conv %>% select(OM, EM, convergence, sim), by = c("OM", "EM", "sim")) %>% 
  filter(convergence == "Converged") %>% 
  mutate(RE = (Pred - True) / True) %>% 
  group_by(OM, EM, Age, Sex, Type) %>% 
  summarize(Median = median(RE),
            lwr_95 = quantile(RE, 0.025),
            upr_95 = quantile(RE, 0.975))

# Fishery Selectivity with Age EM
pdf(here("figs", "Experiment 1", "RE_Fish_Selex.pdf"), height = 8, width = 15)
fish_selex = print(
  ggplot(exp1_selex_sum %>% filter(Type == "Fishery Selectivity",
                                   Sex == 1), 
         aes(x = Age, y = Median, ymin = lwr_95, ymax = upr_95,
             fill = factor(EM), color = factor(EM), group = factor(EM))) +
    geom_line(size = 3) +
    geom_ribbon(alpha = 0.35) +
    geom_hline(yintercept = 0, lty = 2, size = 1.3) + 
    labs(x = "Age", y = "Relative Error in Fishery Selectivity", color = "Sex", fill = "Sex") +
    facet_wrap(~OM, ncol = 4) +
    theme_tj() +
    theme(legend.position = "top")
)
dev.off()


### Parameter Summary -------------------------------------------------------
# summarize re
exp1_param_sum = exp1_param_df %>% 
  select(-Convergence) %>% 
  left_join(exp1_conv %>% select(OM, EM, convergence, sim), by = c("OM", "EM", "sim")) %>% 
  filter(convergence == "Converged") %>% 
  mutate(RE = (as.numeric(Pred) - Truth) / Truth) %>% 
  group_by(OM, EM, Type) %>% 
  summarize(Median = median(RE),
            lwr_75 = quantile(RE, 0.125),
            upr_75 = quantile(RE, 0.875),
            lwr_95 = quantile(RE, 0.025),
            upr_95 = quantile(RE, 0.975))

# plot all other parameters and EMs
pdf(here("figs", "Experiment 1", "RE_ParamAllEMs.pdf"), width = 15, height = 10)
# Comparing proportions within variants
print(ggplot() +
        geom_pointrange(exp1_param_sum,  # 95% quantiles
                        mapping = aes(x = Type, y = Median, ymin = lwr_95, ymax = upr_95, color = EM),
                        position = position_dodge2(width = 0.85), 
                        size = 0, linewidth = 1, alpha = 1) +
        geom_pointrange(exp1_param_sum,  # 75% quantiles
                        mapping = aes(x = Type, y = Median, ymin = lwr_75, ymax = upr_75, color = EM),
                        position = position_dodge2(width = 0.85), 
                        size = 1.3, linewidth = 2, alpha = 0.8) +
        facet_wrap(~OM, ncol = 4) +
        geom_hline(yintercept = 0, lty = 2, size = 1.3) + 
        labs(x = "Parameter", y = "Relative Error") +
        theme_tj() +
        theme(legend.position = "top") +
        scale_x_discrete(guide = guide_axis(angle = 90))  )
dev.off()  


### Time Series Summary -----------------------------------------------------
# time series summary
ts_exp1_sum = exp1_ts_df %>% 
  select(-Convergence) %>%
  left_join(exp1_conv %>% select(OM, EM, convergence, sim), by = c("OM", "EM", "sim")) %>%
  filter(convergence == "Converged") %>% 
  mutate(RE = (as.numeric(Pred)-Truth)/Truth) %>% 
  group_by(Years, Type, OM, EM) %>% 
  summarize(Median = median(RE),
            lwr_75 = quantile(RE, 0.125),
            upr_75 = quantile(RE, 0.875),
            lwr_95 = quantile(RE, 0.025),
            upr_95 = quantile(RE, 0.975))

# plot RE as a time series
pdf(here("figs", 'Experiment 1', "RE_TS_Acr_vs_With.pdf"), width = 15, height = 10)
print(
  ggplot(ts_exp1_sum %>% filter(Type == "Total Biomass"),
         aes(x = Years, y = Median, ymin = lwr_95, ymax = upr_95, fill = EM)) +
    geom_line(size = 1.3) +
    geom_ribbon(alpha = 0.5) +
    facet_wrap(~OM, ncol = 4) +
    geom_hline(yintercept = 0, lty = 2, size = 1.3) + 
    labs(x = "Year", y = "Relative Error in Total Biomass") +
    theme_tj() +
    theme(legend.position = "top") 
)

print(
  ggplot(ts_exp1_sum %>% filter(Type == "Spawning Stock Biomass"),
         aes(x = Years, y = Median, ymin = lwr_95, ymax = upr_95, fill = EM)) +
    geom_line(size = 1.3) +
    geom_ribbon(alpha = 0.5) +
    facet_wrap(~OM, ncol = 4) +
    geom_hline(yintercept = 0, lty = 2, size = 1.3) + 
    labs(x = "Year", y = "Relative Error in Spawning Stock Biomass") +
    theme_tj() +
    theme(legend.position = "top") 
  
)

print(ggplot(ts_exp1_sum %>% filter(Type == "Total Fishing Mortality"),
             aes(x = Years, y = Median, ymin = lwr_95, ymax = upr_95, fill = EM)) +
        geom_line(size = 1.3) +
        geom_ribbon(alpha = 0.5) +
        facet_wrap(~OM, ncol = 4) +
        geom_hline(yintercept = 0, lty = 2, size = 1.3) + 
        labs(x = "Year", y = "Relative Error in Total Fishing Mortality") +
        theme_tj() +
        theme(legend.position = "top") )

print(
  ggplot(ts_exp1_sum %>% filter(Type == "Total Recruitment"),
         aes(x = Years, y = Median, ymin = lwr_95, ymax = upr_95, fill = EM)) +
    geom_line(size = 1.3) +
    geom_ribbon(alpha = 0.5) +
    facet_wrap(~OM, ncol = 4) +
    geom_hline(yintercept = 0, lty = 2, size = 1.3) + 
    labs(x = "Year", y = "Relative Error in Total Recruitment") +
    theme_tj() +
    theme(legend.position = "top") 
)

dev.off()


# Time Series CV plot -----------------------------------------------------

# Set up data for an F-test
exp1_ts_df <- exp1_ts_df %>% 
  select(-Convergence) %>%
  left_join(exp1_conv %>% select(OM, EM, convergence, sim), 
            by = c("OM", "EM", "sim")) %>%
  filter(convergence == "Converged") %>%
  mutate(RE = (as.numeric(Pred)-Truth)/Truth) 

# get statistics for f-test
exp1_f_df <- exp1_ts_df %>% 
  group_by(Type, OM) %>%
  summarize(p_val = var.test(RE ~ factor(EM))$p.value,
            ratio = var.test(RE ~ factor(EM))$estimate)

# plot var
exp1_var_df <- exp1_ts_df %>% 
  group_by(Years, Type, OM, EM) %>%
  summarize(var = var(RE)) %>% 
  ungroup() %>% 
  left_join(exp1_f_df, by = c("OM", "Type"))

pdf(here("figs", 'Experiment 1', "Var_RE_TS_Acr_vs_With.pdf"), width = 13)
print(
  ggplot(exp1_var_df %>% filter(Type == "Spawning Stock Biomass"),
         aes(x = Years, y = var, color = EM)) +
    geom_line(size = 1.3, alpha = 0.85) +
    geom_text(aes(x = 8, y = Inf, label = paste("p = ", round(p_val, 5))), 
              vjust = 4, color = "black", check_overlap = TRUE) +
    geom_text(aes(x = 8, y = Inf, label = paste("F-ratio = ", round(ratio, 5))), 
              vjust = 6, color = "black", check_overlap = TRUE) +
    facet_wrap(~OM, ncol = 4) +
    labs(x = "Year", y = "Variance Relative Error in Spawning Stock Biomass") +
    theme_tj() +
    theme(legend.position = "top") 
)

print(
  ggplot(exp1_var_df %>% filter(Type == "Total Biomass"),
         aes(x = Years, y = var, color = EM)) +
    geom_line(size = 1.3, alpha = 0.85) +
    geom_text(aes(x = 8, y = Inf, label = paste("p = ", round(p_val, 5))), 
              vjust = 4, color = "black", check_overlap = TRUE) +
    geom_text(aes(x = 8, y = Inf, label = paste("F-ratio = ", round(ratio, 5))), 
              vjust = 6, color = "black", check_overlap = TRUE) +
    facet_wrap(~OM, ncol = 4) +
    labs(x = "Year", y = "Variance Relative Error in Total Biomass") +
    theme_tj() +
    theme(legend.position = "top") 
)

print(
  ggplot(exp1_var_df %>% filter(Type == "Total Fishing Mortality"),
         aes(x = Years, y = var, color = EM)) +
    geom_line(size = 1.3, alpha = 0.85) +
    geom_text(aes(x = 8, y = Inf, label = paste("p = ", round(p_val, 5))), 
              vjust = 4, color = "black", check_overlap = TRUE) +
    geom_text(aes(x = 8, y = Inf, label = paste("F-ratio = ", round(ratio, 5))), 
              vjust = 6, color = "black", check_overlap = TRUE) +
    facet_wrap(~OM, ncol = 4) +
    labs(x = "Year", y = "Variance Relative Error in Total Fishing Mortality") +
    theme_tj() +
    theme(legend.position = "top") 
)

print(
  ggplot(exp1_var_df %>% filter(Type == "Total Recruitment"),
         aes(x = Years, y = var, color = EM)) +
    geom_line(size = 1.3, alpha = 0.85) +
    geom_text(aes(x = 8, y = Inf, label = paste("p = ", round(p_val, 5))), 
              vjust = 4, color = "black", check_overlap = TRUE) +
    geom_text(aes(x = 8, y = Inf, label = paste("F-ratio = ", round(ratio, 5))), 
              vjust = 6, color = "black", check_overlap = TRUE) +
    facet_wrap(~OM, ncol = 4) +
    labs(x = "Year", y = "Variance of Relative Error in Total Recruitment") +
    theme_tj() +
    theme(legend.position = "top") 
)
dev.off()


# Coverage ----------------------------------------------------------------

# Compute coverage statistics
coverage_df = exp1_cov_df %>% 
  drop_na() %>% 
  filter(Convergence == "Converged") %>%
  group_by(year, name, EM, OM) %>% 
  summarize(sum = sum(coverage), # get sum of coverage
         unique_rows = length(unique(sim)),# get length of unique simulations
         prop = sum / unique_rows) # get coverage

pdf(here("figs", 'Experiment 1', "Coverage.pdf"), width = 13)
coverage_df %>%
  ggplot(aes(x = year, y = prop, color = EM)) +
  geom_point(size = 2) +
  geom_line(size = 1.3) +
  facet_grid(name~OM) +
  geom_hline(yintercept = 0.95, lty = 2) +
  labs(x = "Year", y = "Coverage") +
  theme_tj() +
  theme(legend.position = "top") 
dev.off()


# Experiment 2 ------------------------------------------------------------
# Read in files
exp2_selex_df = data.table::fread(here("output", "Experiment_2_Selex.csv")) 
exp2_growth_df = data.table::fread(here("output", "Experiment_2_Growth.csv")) 
exp2_param_df = data.table::fread(here("output", "Experiment_2_Param.csv")) 
exp2_ts_df = data.table::fread(here("output", "Experiment_2_TimeSeries.csv"))  
exp2_conv_df = data.table::fread(here("output", "Experiment_2_Convergence.csv")) 
exp2_cov_df = data.table::fread(here("output", "Experiment_2_Coverage.csv")) 

### Convergence Summary -----------------------------------------------------
exp2_conv_df = exp2_conv_df %>%
  mutate(convergence =
           case_when( # munging gradient stuff
             (pdHess == TRUE &
                gradient <= 0.001 &
                sdNA == FALSE &
                convergence == "Not Converged") ~ "Converged",
             (pdHess == TRUE &
                gradient <= 0.001 &
                sdNA == FALSE &
                convergence == "Converged") ~ "Converged",
             (pdHess == FALSE |
                gradient > 0.001 |
                sdNA == TRUE) ~ "Not Converged"
           ))
# Convergence summary
conv_df = exp2_conv_df %>% 
  filter(convergence == "Converged") %>% 
  group_by(OM, EM) %>% 
  summarize(sum = n())

# Plot convergence
pdf(here("figs", "Experiment 2", "Convergence.pdf"), width = 15)
ggplot(conv_df, aes(x = OM, y = sum/750, group = EM, color = EM)) +
  geom_point(size = 3) +
  geom_line() +
  theme_tj() +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  labs(x = "Operating Models", y = "Convergence Rate") +
  theme(legend.position = "top")
dev.off()

### Selectivity Summary -----------------------------------------------------

# summarize relative error in selectivity
exp2_selex_sum = exp2_selex_df %>% 
  select(-Convergence) %>% 
  left_join(exp2_conv_df %>% select(OM, EM, convergence, sim), by = c("OM", "EM", "sim")) %>% 
  filter(convergence == "Converged") %>% 
  mutate(RE = (Pred - True) / True) %>% 
  group_by(OM, EM, Age, Sex, Type) %>% 
  summarize(Median = median(RE),
            lwr_95 = quantile(RE, 0.025),
            upr_95 = quantile(RE, 0.975))

# Fishery Selectivity
pdf(here("figs", "Experiment 2", "RE_Fish_Selex.pdf"), width = 15, height = 12)
fish_selex = print(
  ggplot(exp2_selex_sum %>% filter(Type == "Fishery Selectivity") %>% 
           mutate(Sex = ifelse(Sex == 1, "Female", "Male")), 
         aes(x = Age, y = Median, ymin = lwr_95, ymax = upr_95,
             fill = factor(Sex), color = factor(Sex), group = factor(Sex))) +
    geom_line(size = 3) +
    geom_ribbon(alpha = 0.35, color = NA) +
    geom_hline(yintercept = 0, lty = 2, size = 1.3) + 
    labs(x = "Age", y = "Relative Error in Fishery Selectivity", 
         color = "Sex", fill = "Sex") +
    facet_grid(OM~EM, scales = "free") +
    theme_tj() +
    theme(legend.position = "top")
)
dev.off()

# fishery selectivity for age only EM
pdf(here("figs", "Experiment 2", "RE_Fish_Selex_Lines.pdf"), width = 10, height = 13)
print(
  exp2_selex_df %>% 
    filter(Convergence == "Converged" ,
           Type == "Fishery Selectivity",
           Sex == 1) %>% 
    mutate(Sex = ifelse(Sex == 1, "Female", "Male")) %>% 
    ggplot() +
    geom_line(aes(x = Age, y = Pred, group = sim)) +
    geom_line(aes(x = Age, y = True, color = factor(Sex)), size = 1.3) +
    facet_grid(EM~OM) +
    labs(x = "Age", y = "Fishery Selectivity", color = "Sex") +
    theme_tj() +
    theme(legend.position = "top")
)
print(
  exp2_selex_df %>% 
    filter(Convergence == "Converged" ,
           Type == "Fishery Selectivity",
           Sex == 2) %>% 
    mutate(Sex = ifelse(Sex == 1, "Female", "Male")) %>% 
    ggplot() +
    geom_line(aes(x = Age, y = Pred, group = sim)) +
    geom_line(aes(x = Age, y = True, color = factor(Sex)), size = 1.3) +
    facet_grid(OM~EM) +
    labs(x = "Age", y = "Fishery Selectivity", color = "Sex") +
    theme_tj() +
    theme(legend.position = "top")
)
dev.off()

# Survey Selectivity
pdf(here("figs", "Experiment 2", "RE_Srv_Selex.pdf"), width = 15, height = 12)
srv_selex = print(
  ggplot(exp2_selex_sum %>% filter(Type == "Survey Selectivity") %>% 
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
pdf(here("figs", "Experiment 2", "RE_Srv_Selex_Lines.pdf"), width = 10, height = 13)
print(
  exp2_selex_df %>% 
    filter(Convergence == "Converged", Sex == 1,
           Type == "Survey Selectivity") %>% 
    mutate(Sex = ifelse(Sex == 1, "Female", "Male")) %>% 
    ggplot() +
    geom_line(aes(x = Age, y = Pred, group = sim)) +
    geom_line(aes(x = Age, y = True, color = factor(Sex)), size = 1.3) +
    facet_grid(OM~EM) +
    labs(x = "Age", y = "Survey Selectivity", color = "Sex") +
    theme_tj() +
    theme(legend.position = "top")
)
print(
  exp2_selex_df %>% 
    filter(Convergence == "Converged", Sex == 2,
           Type == "Survey Selectivity") %>% 
    mutate(Sex = ifelse(Sex == 1, "Female", "Male")) %>% 
    ggplot() +
    geom_line(aes(x = Age, y = Pred, group = sim)) +
    geom_line(aes(x = Age, y = True, color = factor(Sex)), size = 1.3) +
    facet_grid(OM~EM) +
    labs(x = "Age", y = "Survey Selectivity", color = "Sex") +
    theme_tj() +
    theme(legend.position = "top")
)
dev.off()

### Parameter Summary -------------------------------------------------------
prop_param_df = exp2_param_df %>% 
  select(-Convergence) %>% 
  left_join(exp2_conv_df %>% select(OM, EM, convergence, sim), by = c("OM", "EM", "sim")) %>% 
  filter(convergence == "Converged") %>% 
  mutate(RE = (as.numeric(Pred) - Truth) / Truth)

##### Other parameters --------------------------------------------------------
# summarize re for parameters
exp2_param_sum = exp2_param_df %>% 
  select(-Convergence) %>% 
  left_join(exp2_conv_df %>% select(OM, EM, convergence, sim), by = c("OM", "EM", "sim")) %>% 
  filter(convergence == "Converged") %>% 
  filter( !str_detect(Type, "Ratio")) %>% 
  mutate(RE = (as.numeric(Pred) - Truth) / Truth) %>% 
  group_by(OM, EM, Type) %>% 
  summarize(Median = median(RE),
            lwr_75 = quantile(RE, 0.125),
            upr_75 = quantile(RE, 0.875),
            lwr_95 = quantile(RE, 0.025),
            upr_95 = quantile(RE, 0.975))

# plot (all EMs)
pdf(here("figs", "Experiment 2", "RE_Param_allEMs.pdf"), width = 15, height = 13)
print(ggplot(exp2_param_sum, 
             aes(x = Type, y = Median, ymin = lwr_95, ymax = upr_95)) +
        geom_pointrange(exp2_param_sum,  # 95% quantiles
                        mapping = aes(x = Type, y = Median, ymin = lwr_95, ymax = upr_95),
                        position = position_dodge2(width = 0.65), 
                        size = 0, linewidth = 1, alpha = 1) +
        geom_pointrange(exp2_param_sum,  # 75% quantiles
                        mapping = aes(x = Type, y = Median, ymin = lwr_75, ymax = upr_75),
                        position = position_dodge2(width = 0.65), 
                        size = 1.5, linewidth = 2, alpha = 0.8) +
        facet_grid(OM~EM, scales = "free") +
        geom_hline(yintercept = 0, lty = 2, size = 0.85) + 
        labs(x = "Parameter", y = "Relative Error") +
        theme_tj() +
        scale_x_discrete(guide = guide_axis(angle = 90)) +
        coord_cartesian(ylim = c(-1.5,1.5)) )
dev.off()  

### Time Series Summary -----------------------------------------------------
###### Relative error all EMs --------------------------------------------------
# time series summary in relative error
ts_exp2_sum = exp2_ts_df %>% 
  select(-Convergence) %>% 
  left_join(exp2_conv_df %>% select(OM, EM, convergence, sim), by = c("OM", "EM", "sim")) %>% 
  filter(convergence == "Converged") %>% 
  mutate(RE = (as.numeric(Pred)-Truth)/Truth) %>% 
  group_by(Years, Type, OM, EM) %>% 
  summarize(Median = median(RE),
            lwr_75 = quantile(RE, 0.125),
            upr_75 = quantile(RE, 0.875),
            lwr_95 = quantile(RE, 0.025),
            upr_95 = quantile(RE, 0.975))

# Plot time series EMs
pdf(here("figs", "Experiment 2", "RE_TS_AllEMs.pdf"), width = 15, height = 13)
print(
  ggplot() +
    geom_ribbon(ts_exp2_sum %>% filter(Type == "Total Biomass"), # 95% simulation intervals
                mapping = aes(x = Years, y = Median, ymin = lwr_95, ymax = upr_95),
                alpha = 1, fill = "grey") +
    geom_ribbon(ts_exp2_sum %>% filter(Type == "Total Biomass"), # 75% simulation intervals
                mapping = aes(x = Years, y = Median, ymin = lwr_75, ymax = upr_75),
                alpha = 0.35, fill = "grey4") +
    geom_point(ts_exp2_sum %>% filter(Type == "Total Biomass"), 
               mapping = aes(x = Years, y = Median), shape = 21, colour = "black", 
               fill = "white", size = 3, stroke = 1, alpha = 1) +
    facet_grid(OM~EM, scales = "free") +
    geom_hline(yintercept = 0, lty = 2, size = 1.3) + 
    labs(x = "Year", y = "Relative Error in Total Biomass") +
    theme_tj() +
    coord_cartesian(ylim = c(-0.85,0.85)) 
)

print(
  ggplot() +
    geom_ribbon(ts_exp2_sum %>% filter(Type == "Spawning Stock Biomass"), # 95% simulation intervals
                mapping = aes(x = Years, y = Median, ymin = lwr_95, ymax = upr_95),
                alpha = 1, fill = "grey") +
    geom_ribbon(ts_exp2_sum %>% filter(Type == "Spawning Stock Biomass"), # 75% simulation intervals
                mapping = aes(x = Years, y = Median, ymin = lwr_75, ymax = upr_75),
                alpha = 0.35, fill = "grey4") +
    geom_point(ts_exp2_sum %>% filter(Type == "Spawning Stock Biomass"), 
               mapping = aes(x = Years, y = Median), shape = 21, colour = "black", 
               fill = "white", size = 3, stroke = 1, alpha = 1) +
    facet_grid(OM~EM, scales = "free") +
    geom_hline(yintercept = 0, lty = 2, size = 1.3) + 
    labs(x = "Year", y = "Relative Error in Spawning Stock Biomass") +
    theme_tj() +
    coord_cartesian(ylim = c(-0.85,0.85)) 
)

print(
  ggplot() +
    geom_ribbon(ts_exp2_sum %>% filter(Type == "Total Fishing Mortality"), # 95% simulation intervals
                mapping = aes(x = Years, y = Median, ymin = lwr_95, ymax = upr_95),
                alpha = 1, fill = "grey") +
    geom_ribbon(ts_exp2_sum %>% filter(Type == "Total Fishing Mortality"), # 75% simulation intervals
                mapping = aes(x = Years, y = Median, ymin = lwr_75, ymax = upr_75),
                alpha = 0.35, fill = "grey4") +
    geom_point(ts_exp2_sum %>% filter(Type == "Total Fishing Mortality"), 
               mapping = aes(x = Years, y = Median), shape = 21, colour = "black", 
               fill = "white", size = 3, stroke = 1, alpha = 1) +
    facet_grid(OM~EM, scales = "free") +
    geom_hline(yintercept = 0, lty = 2, size = 1.3) + 
    labs(x = "Year", y = "Relative Error in Total Fishing Mortality") +
    theme_tj() +
    coord_cartesian(ylim = c(-1,1)) 
)

print(
  ggplot() +
    geom_ribbon(ts_exp2_sum %>% filter(Type == "Total Recruitment"), # 95% simulation intervals
                mapping = aes(x = Years, y = Median, ymin = lwr_95, ymax = upr_95),
                alpha = 1, fill = "grey") +
    geom_ribbon(ts_exp2_sum %>% filter(Type == "Total Recruitment"), # 75% simulation intervals
                mapping = aes(x = Years, y = Median, ymin = lwr_75, ymax = upr_75),
                alpha = 0.35, fill = "grey4") +
    geom_point(ts_exp2_sum %>% filter(Type == "Total Recruitment"), 
               mapping = aes(x = Years, y = Median), shape = 21, colour = "black", 
               fill = "white", size = 3, stroke = 1, alpha = 1) +
    facet_grid(OM~EM, scales = "free") +
    geom_hline(yintercept = 0, lty = 2, size = 1.3) + 
    labs(x = "Year", y = "Relative Error in Total Recruitment") +
    theme_tj() +
    coord_cartesian(ylim = c(-1,1)) 
)
dev.off()

# Age only EM
pdf(here("figs", "Experiment 2", "RE_TS_AgeEM.pdf"), width = 15, height = 13)
print(
  ggplot() +
    geom_ribbon(ts_exp2_sum %>% filter(EM == "Age"), # 95% simulation intervals
                mapping = aes(x = Years, y = Median, ymin = lwr_95, ymax = upr_95),
                alpha = 1, fill = "grey") +
    geom_ribbon(ts_exp2_sum %>% filter(EM == "Age"), # 75% simulation intervals
                mapping = aes(x = Years, y = Median, ymin = lwr_75, ymax = upr_75),
                alpha = 0.35, fill = "grey4") +
    geom_point(ts_exp2_sum %>% filter(EM == "Age"), 
               mapping = aes(x = Years, y = Median), shape = 21, colour = "black", 
               fill = "white", size = 3, stroke = 1, alpha = 1) +
    facet_grid(Type~OM, scales = "free") +
    geom_hline(yintercept = 0, lty = 2, size = 1.3) + 
    labs(x = "Year", y = "Relative Error") +
    theme_tj()) 

dev.off()


# Coverage ----------------------------------------------------------------

# Compute coverage statistics
coverage_df = exp2_cov_df %>% 
  drop_na() %>% 
  filter(Convergence == "Converged") %>%
  group_by(year, name, EM, OM) %>% 
  summarize(sum = sum(coverage), # get sum of coverage
            unique_rows = length(unique(sim)),# get length of unique simulations
            prop = sum / unique_rows) # get coverage

pdf(here("figs", 'Experiment 2', "Coverage.pdf"), width = 15)
coverage_df %>%
  ggplot(aes(x = year, y = prop, color = EM)) +
  geom_point(size = 2) +
  geom_line(size = 1.3) +
  facet_grid(name~OM) +
  geom_hline(yintercept = 0.95, lty = 2) +
  labs(x = "Year", y = "Coverage (SSB)") +
  theme_tj() +
  theme(legend.position = "top") 
dev.off()


# Experiment 3 ------------------------------------------------------------

# Read in files
exp3_conv = data.table::fread(here("output", "Experiment_3_Convergence.csv")) %>% filter(!str_detect(OM, "No"))
exp3_selex_df = data.table::fread(here("output", "Experiment_3_Selex.csv")) %>% filter(!str_detect(OM, "No"))
exp3_growth_df = data.table::fread(here("output", "Experiment_3_Growth.csv"))%>% filter(!str_detect(OM, "No"))
exp3_param_df = data.table::fread(here("output", "Experiment_3_Param.csv")) %>% filter(!str_detect(OM, "No"))
exp3_ts_df = data.table::fread(here("output", "Experiment_3_TimeSeries.csv"))%>% filter(!str_detect(OM, "No"))
exp3_sr_df = data.table::fread(here("output", "Experiment_3_SexRatio.csv")) %>% filter(!str_detect(OM, "No"))
exp3_cov_df = data.table::fread(here("output", "Experiment_3_Coverage.csv")) %>% filter(!str_detect(OM, "No"))

##### Convergence summary -----------------------------------------------------
exp3_conv = exp3_conv %>% 
  mutate(convergence = 
           case_when( # munging gradient stuff
             (pdHess == TRUE & 
                gradient <= 0.001 &
                sdNA == FALSE &
                convergence == "Not Converged") ~ "Converged",
             (pdHess == TRUE & 
                gradient <= 0.001 &
                sdNA == FALSE &
                convergence == "Converged") ~ "Converged",
             (pdHess == FALSE |
                gradient > 0.001 |
                sdNA == TRUE) ~ "Not Converged"
           ))

# Convergence summary
conv_df = exp3_conv %>% 
  filter(convergence == "Converged") %>% 
  group_by(OM, EM) %>% 
  summarize(sum = n())

# Plot convergence
pdf(here("figs", "Experiment 3", "Convergence.pdf"))
ggplot(conv_df, aes(x = OM, y = sum/750, color = EM, group = EM)) +
  geom_point(size = 3) +
  geom_line(size = 1) +
  theme_tj() +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  theme(legend.position = "top") +
  labs(x = "Operating Models", y = "Convergence Rate") 
dev.off()


### Selectivity Summary -----------------------------------------------------
# summarize relative error
exp3_selex_sum = exp3_selex_df %>% 
  select(-Convergence) %>% 
  left_join(exp3_conv %>% select(OM, EM, convergence, sim), by = c("OM", "EM", "sim")) %>% 
  filter(convergence == "Converged") %>% 
  mutate(RE = (Pred - True) / True) %>% 
  group_by(OM, EM, Age, Sex, Type) %>% 
  summarize(Median = median(RE),
            lwr_95 = quantile(RE, 0.025),
            upr_95 = quantile(RE, 0.975))

# Fishery Selectivity with Age EM
pdf(here("figs", "Experiment 3", "RE_Fish_Selex.pdf"), height = 14, width = 20)
fish_selex = print(
  ggplot(exp3_selex_sum %>% filter(Type == "Fishery Selectivity") %>% 
           mutate(Sex = ifelse(Sex == 1, "Female", "Male")), 
         aes(x = Age, y = Median, ymin = lwr_95, ymax = upr_95,
             fill = factor(Sex), color = factor(Sex), group = factor(Sex))) +
    geom_line(size = 3) +
    geom_ribbon(alpha = 0.35, color = NA) +
    geom_hline(yintercept = 0, lty = 2, size = 1.3) + 
    labs(x = "Age", y = "Relative Error in Fishery Selectivity", color = "Sex", fill = "Sex") +
    facet_grid(EM~OM, scales = "free") +
    theme_tj() +
    theme(legend.position = "top")
)
dev.off()

# Fishery Selectivity with Age EM
pdf(here("figs", "Experiment 3", "RE_Srv_Selex.pdf"), height = 14, width = 20)
srv_selex = print(
  ggplot(exp3_selex_sum %>% filter(Type == "Survey Selectivity") %>% 
           mutate(Sex = ifelse(Sex == 1, "Female", "Male")), 
         aes(x = Age, y = Median, ymin = lwr_95, ymax = upr_95,
             fill = factor(Sex), color = factor(Sex), group = factor(Sex))) +
    geom_line(size = 3) +
    geom_ribbon(alpha = 0.35, color = NA) +
    geom_hline(yintercept = 0, lty = 2, size = 1.3) + 
    labs(x = "Age", y = "Relative Error in Survey Selectivity", color = "Sex", fill = "Sex") +
    facet_grid(EM~OM, scales = "free") +
    theme_tj() +
    theme(legend.position = "top")
)
dev.off()

### Parameter Summary -------------------------------------------------------
# summarize re
exp3_param_df = exp3_param_df %>% 
  select(-Convergence) %>% 
  left_join(exp3_conv %>% select(OM, EM, convergence, sim), by = c("OM", "EM", "sim")) %>% 
  filter(convergence == "Converged") 

exp3_param_sum = exp3_param_df%>% 
  mutate(RE = (as.numeric(Pred) - Truth) / Truth) %>% 
  group_by(OM, EM, Type) %>% 
  summarize(Median = median(RE),
            lwr_75 = quantile(RE, 0.125),
            upr_75 = quantile(RE, 0.875),
            lwr_95 = quantile(RE, 0.025),
            upr_95 = quantile(RE, 0.975))

# plot all other parameters and EMs
pdf(here("figs", "Experiment 3", "RE_ParamAllEMs.pdf"), width = 10, height = 10)
# Comparing proportions within variants
print(ggplot() +
        geom_pointrange(exp3_param_sum,  # 95% quantiles
                        mapping = aes(x = Type, y = Median, ymin = lwr_95, ymax = upr_95),
                        position = position_dodge2(width = 0.65), 
                        size = 0, linewidth = 1, alpha = 1) +
        geom_pointrange(exp3_param_sum,  # 75% quantiles
                        mapping = aes(x = Type, y = Median, ymin = lwr_75, ymax = upr_75),
                        position = position_dodge2(width = 0.65), 
                        size = 1.5, linewidth = 2, alpha = 0.8) +
        facet_grid(OM~EM, scales = "free_y") +
        geom_hline(yintercept = 0, lty = 2, size = 1) + 
        labs(x = "Parameter", y = "Relative Error") +
        theme_tj() +
        theme(legend.position = "top") +
        scale_x_discrete(guide = guide_axis(angle = 90))  )
dev.off()  

# Sex ratio estimability
pdf(here("figs", "Experiment 3", "RE_SexRatio.pdf"))
# Comparing proportions within variants
print(ggplot(exp3_param_df %>% filter(Type == "Female Sex Ratio",
                                      convergence == "Converged",
                                      EM == "EstSR"),
             aes(x = Type, y = (Pred - Truth) / Truth)) +
        geom_violin(width = 0.5) +
        geom_boxplot(width = 0.1) +
        facet_wrap(~OM) +
        geom_hline(yintercept = 0, lty = 2, lwd = 1) +
        labs(x = "", y = "Relative Error in Initial Sex Ratio") +
        theme_tj() +
        theme(legend.position = "top",
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank()) +
        scale_x_discrete(guide = guide_axis(angle = 90))  )
dev.off()  



 ### Time Series Summary -----------------------------------------------------
# time series summary
ts_exp3_sum = exp3_ts_df %>% 
  select(-Convergence) %>%
  left_join(exp3_conv %>% select(OM, EM, convergence, sim), by = c("OM", "EM", "sim")) %>%
  filter(convergence == "Converged") %>% 
  mutate(RE = (as.numeric(Pred)-Truth)/Truth) %>% 
  group_by(Years, Type, OM, EM) %>% 
  summarize(Median = median(RE),
            lwr_75 = quantile(RE, 0.125),
            upr_75 = quantile(RE, 0.875),
            lwr_95 = quantile(RE, 0.025),
            upr_95 = quantile(RE, 0.975))

# plot RE as a time series
pdf(here("figs", 'Experiment 3', "RE_TS_All.pdf"), width = 13, height = 13)
# comparing within variants
print(
  ggplot(ts_exp3_sum,
         aes(x = Years, y = Median, ymin = lwr_95, ymax = upr_95, fill = EM, color = EM)) +
    geom_line(size = 1.3) +
    geom_ribbon(alpha = 0.35, color = NA) +
    facet_grid(Type~OM, scales = "free") +
    geom_hline(yintercept = 0, lty = 2, size = 1.3) + 
    labs(x = "Year", y = "Relative Error") +
    theme_tj() +
    theme(legend.position = "top") 
)
dev.off()

# Only plot RE time series for Fixed EM
pdf(here("figs", 'Experiment 3', "RE_TS_FixEM.pdf"), width = 13, height = 11)
print(
  ggplot() +
    geom_ribbon(ts_exp3_sum %>% filter(str_detect(EM, "Fix")), # 95% simulation intervals
                mapping = aes(x = Years, y = Median, ymin = lwr_95, ymax = upr_95),
                alpha = 1, fill = "grey") +
    geom_ribbon(ts_exp3_sum %>% filter(str_detect(EM, "Fix")), # 75% simulation intervals
                mapping = aes(x = Years, y = Median, ymin = lwr_75, ymax = upr_75),
                alpha = 0.35, fill = "grey4") +
    geom_point(ts_exp3_sum %>% filter(str_detect(EM, "Fix")), 
               mapping = aes(x = Years, y = Median), shape = 21, colour = "black", 
               fill = "white", size = 3, stroke = 1, alpha = 1) +
    facet_grid(Type~OM, scales = "free") +
    geom_hline(yintercept = 0, lty = 2, size = 1.3) + 
    labs(x = "Year", y = "Relative Error") +
    theme_tj() +
    theme(legend.position = "top") 
)
dev.off()



 # Coverage ----------------------------------------------------------------

# Compute coverage statistics
coverage_df = exp3_cov_df %>% 
  drop_na() %>% 
  filter(Convergence == "Converged") %>%
  group_by(year, name, EM, OM) %>% 
  summarize(sum = sum(coverage), # get sum of coverage
            unique_rows = length(unique(sim)),# get length of unique simulations
            prop = sum / unique_rows) # get coverage

pdf(here("figs", 'Experiment 3', "Coverage.pdf"), width = 13)
coverage_df %>%
  ggplot(aes(x = year, y = prop, color = EM)) +
  geom_point(size = 3) +
  geom_line(size = 1, alpha = 0.75) +
  geom_hline(yintercept = 0.95, lty = 2) +
  facet_grid(name~OM) +
  labs(x = "Year", y = "Coverage") +
  theme_tj() +
  theme(legend.position = "top") 
dev.off()


