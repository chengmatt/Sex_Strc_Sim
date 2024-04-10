# Purpose: To plot manuscript-ready figures
# Creator: Matthew LH. Cheng (UAF-CFOS)
# Date: 3/28/24

# Set up ------------------------------------------------------------------

library(here)
library(tidyverse)
library(patchwork)

# Set up quick theme
theme_tj <- function() {
  theme_bw() +
    theme(legend.position = "none",
          axis.text = element_text(size = 18, color = "black"),
          axis.title = element_text(size = 20, color =  "black"),
          strip.text = element_text(size = 13),
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 18),
          plot.title = element_text(size = 20))
}

theme_set(theme_tj())
dir.create(here('figs', "ms_figs"))

# Experiment 1 ------------------------------------------------------------
# Read in data for plotting
exp1_ts_df <- data.table::fread(here("output", "Experiment_1_TimeSeries.csv"))  
exp1_cov_df <- data.table::fread(here("output", "Experiment_1_Coverage.csv")) 
exp1_param_df <- data.table::fread(here("output", "Experiment_1_Param.csv")) 
  
# Join time-series with new convergence diagnostics
ts_exp1_sum <- exp1_ts_df %>% 
  filter(Convergence == "Converged") %>% 
  mutate(RE = (as.numeric(Pred)-Truth)/Truth) %>% 
  group_by(Years, Type, OM, EM) %>% 
  summarize(Median = median(RE),
            lwr_75 = quantile(RE, 0.125),
            upr_75 = quantile(RE, 0.875),
            lwr_95 = quantile(RE, 0.025),
            upr_95 = quantile(RE, 0.975)) %>% 
  filter(str_detect(OM, "_50")) %>%  # ISS summed for both sexes = 50
  mutate(OM = str_remove(OM, "_50"))

# summarize re
exp1_param_sum <- exp1_param_df %>% 
  filter(Convergence == "Converged") %>% 
  mutate(RE = (as.numeric(Pred) - Truth) / Truth) %>% 
  group_by(OM, EM, Type) %>% 
  summarize(Median = median(RE),
            lwr_75 = quantile(RE, 0.125),
            upr_75 = quantile(RE, 0.875),
            lwr_95 = quantile(RE, 0.025),
            upr_95 = quantile(RE, 0.975)) %>% 
  filter(str_detect(OM, "_50"),
         Type %in% c("Bmsy", "M_F", "M_M", "R0", "Tier 3 HCR Catch")) %>%  # ISS summed for both sexes = 50
  mutate(OM = str_remove(OM, "_50"),
         Type = case_when(
         Type == "Bmsy" ~ "BMSY",
         Type == "M_F" ~ "F NatMort",
         Type == "M_M" ~ "M NatMort",
         Type == "R0" ~ "R0",
         Type == "Tier 3 HCR Catch" ~ "HCR Catch"
         ), Type = factor(Type, levels = c( "R0", "BMSY", "F NatMort", "M NatMort", "HCR Catch")))

# Figure 2 (Time Series Rel Err + Coverage + Parameters) ------------------------------------------

# Time series relative error
exp1_ts_plot <- ggplot() +
  geom_ribbon(ts_exp1_sum %>% filter(Type == "Spawning Stock Biomass"),
              mapping = aes(x = Years, ymin = lwr_95, ymax = upr_95, fill = EM),
              alpha = 0.3) +
  geom_ribbon(ts_exp1_sum %>% filter(Type == "Spawning Stock Biomass"),
              mapping = aes(x = Years, ymin = lwr_75, ymax = upr_75, fill = EM),
              alpha = 0.4) +
  geom_line(ts_exp1_sum %>% filter(Type == "Spawning Stock Biomass"),
            mapping = aes(x = Years, y = Median, color = EM), size = 1.8) +
  scale_fill_manual(values = c("#e69c4c", "#35577b")) +
  scale_color_manual(values = c("#e69c4c", "#35577b")) +
  facet_wrap(~OM) +
  geom_hline(yintercept = 0, lty = 2, size = 1.3) + 
  labs(x = "Year", y = "Relative Error", color = "Estimation Model", fill = "Estimation Model") +
  theme_tj() +
  theme(legend.position = "top") 

# Comparing proportions within variants
exp1_param_plot <- ggplot() +
  geom_pointrange(exp1_param_sum,  # 95% quantiles
                  mapping = aes(x = Type, y = Median, ymin = lwr_95, ymax = upr_95, color = EM, fill = EM),
                  position = position_dodge2(width = 0.65), 
                  size = 0, linewidth = 1, alpha = 1) +
  geom_pointrange(exp1_param_sum,  # 75% quantiles
                  mapping = aes(x = Type, y = Median, ymin = lwr_75, ymax = upr_75, color = EM, fill = EM),
                  position = position_dodge2(width = 0.65), 
                  size = 1.5, linewidth = 2, alpha = 0.8) +
  facet_wrap(~OM) +
  scale_fill_manual(values = c("#e69c4c", "#35577b")) +
  scale_color_manual(values = c("#e69c4c", "#35577b")) +
  geom_hline(yintercept = 0, lty = 2, size = 1) + 
  labs(x = "Parameter", y = "Relative Error", color = "Estimation Model", fill = "Estimation Model") +
  theme_tj() +
  theme(legend.position = "none") +
  scale_x_discrete(guide = guide_axis(angle = 90))

pdf(here("figs", "ms_figs", "Fig2_Exp1_AcrWith.pdf"), width = 10, height = 10)
exp1_ts_re/exp1_param_plot + plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 23, face = "bold"))
dev.off()



# Experiment 2 ------------------------------------------------------------

# Read in files
exp2_param_df <- data.table::fread(here("output", "Experiment_2_Param.csv")) %>% filter(!str_detect(OM, "Growth_M"))
exp2_ts_df <- data.table::fread(here("output", "Experiment_2_TimeSeries.csv"))  %>% filter(!str_detect(OM, "Growth_M"))
exp2_selex_df <- data.table::fread(here("output", "Experiment_2_Selex.csv")) %>% filter(!str_detect(OM, "Growth_M"))

# Set up EM Names here
exp2_names <- c(
  "AgeStrc", 
  "SexStrc (MAgg, CatAgg)",
  "SexStrc (MAgg, CatSex)",
  "SexStrc (MSex, CatAgg)",
  "SexStrc (MSex, CatSex)"
)

# summarize parameters
exp2_param_df <- exp2_param_df %>% 
  filter(Convergence == "Converged",
         !str_detect(EM, "AgeSel|AgeAgg")) %>% 
  mutate(RE = (as.numeric(Pred) - Truth) / Truth,
         EM = factor(EM, labels = exp2_names)) 

# get relative error distributions and rename stuff
exp2_param_sum <- exp2_param_df %>% 
  group_by(OM, EM, Type) %>% 
  summarize(Median = median(RE),
            lwr_75 = quantile(RE, 0.125),
            upr_75 = quantile(RE, 0.875),
            lwr_95 = quantile(RE, 0.025),
            upr_95 = quantile(RE, 0.975)) %>% 
  filter(Type %in% c("Bmsy", "M_F", "M_M", "R0", "Tier 3 HCR Catch")) %>%  # ISS summed for both sexes = 50
  mutate(Type = case_when(
           Type == "Bmsy" ~ "BMSY",
           Type == "M_F" ~ "F NatMort",
           Type == "M_M" ~ "M NatMort",
           Type == "R0" ~ "R0",
           Type == "Tier 3 HCR Catch" ~ "HCR Catch"
         ), Type = factor(Type, levels = c( "R0", "BMSY", "F NatMort", "M NatMort", "HCR Catch")))

# time series summary in relative error
ts_exp2_sum <- exp2_ts_df %>% 
  filter(Convergence == "Converged",
         !str_detect(EM, "AgeSel|AgeAgg")) %>%
  mutate(RE = (as.numeric(Pred) - Truth) / Truth,
         EM = factor(EM, labels = exp2_names)) %>% 
  group_by(Years, Type, OM, EM) %>% 
  summarize(Median = median(RE),
            lwr_75 = quantile(RE, 0.125),
            upr_75 = quantile(RE, 0.875),
            lwr_95 = quantile(RE, 0.025),
            upr_95 = quantile(RE, 0.975))

# summarize relative error in selectivity
exp2_selex_sum <- exp2_selex_df %>% 
  filter(Convergence == "Converged",
         !str_detect(EM, "AgeSel|AgeAgg")) %>% 
  mutate(EM = factor(EM, labels = exp2_names))

# Figure 3 (Spawning Biomass Trends) --------------------------------------

pdf(here("figs", "ms_figs", "Fig3_Exp2_SSBTrends.pdf"), width = 15, height = 10)
ggplot() +
  geom_ribbon(ts_exp2_sum %>% filter(Type == "Spawning Stock Biomass"), # 95% simulation intervals
              mapping = aes(x = Years, y = Median, ymin = lwr_95, ymax = upr_95),
              alpha = 0.3, fill = "#35577b") +
  geom_ribbon(ts_exp2_sum %>% filter(Type == "Spawning Stock Biomass"), # 75% simulation intervals
              mapping = aes(x = Years, y = Median, ymin = lwr_75, ymax = upr_75),
              alpha = 0.5, fill = "#35577b") +
  geom_line(ts_exp2_sum %>% filter(Type == "Spawning Stock Biomass"), 
            mapping = aes(x = Years, y = Median), 
            colour = "#35577b", size = 1.3, alpha = 1) +
  facet_grid(OM~EM, scales = "free") +
  geom_hline(yintercept = 0, lty = 2, size = 1.3) + 
  labs(x = "Year", y = "Relative Error") +
  theme_tj() +
  coord_cartesian(ylim = c(-0.5, 0.5))
dev.off()


# Figure 4 (Experiment 2 Parameters) ----------------------------------------------------------------

pdf(here("figs", "ms_figs", "Fig4_Exp2_Params.pdf"), width = 15, height = 10)
ggplot(exp2_param_sum, 
             aes(x = Type, y = Median, ymin = lwr_95, ymax = upr_95)) +
        geom_pointrange(exp2_param_sum,  # 95% quantiles
                        mapping = aes(x = Type, y = Median, ymin = lwr_95, ymax = upr_95),
                        position = position_dodge2(width = 0.65), 
                        size = 0, linewidth = 1, alpha = 1, 
                        fill = "#35577b", color = "#35577b") +
        geom_pointrange(exp2_param_sum,  # 75% quantiles
                        mapping = aes(x = Type, y = Median, ymin = lwr_75, ymax = upr_75),
                        position = position_dodge2(width = 0.65), 
                        size = 1.5, linewidth = 2, alpha = 0.8,
                        fill = "#35577b", color = "#35577b") +
        facet_grid(OM~EM, scales = "free") +
        geom_hline(yintercept = 0, lty = 2, size = 0.85) + 
        labs(x = "Parameter", y = "Relative Error") +
        theme_tj() +
        scale_x_discrete(guide = guide_axis(angle = 90)) +
        coord_cartesian(ylim = c(-1.5,1.5))
dev.off()


# Experiment 3 ------------------------------------------------------------
exp3_param_df <- data.table::fread(here("output", "Experiment_3_Param.csv")) %>% filter(!str_detect(OM, "No"))
exp3_ts_df <- data.table::fread(here("output", "Experiment_3_TimeSeries.csv"))%>% filter(!str_detect(OM, "No"))

# summarize re
exp3_param_df <- exp3_param_df %>% 
  filter(Convergence == "Converged") 

exp3_param_sum <- exp3_param_df %>% 
  mutate(RE = (as.numeric(Pred) - Truth) / Truth) %>% 
  group_by(OM, EM, Type) %>% 
  summarize(Median = median(RE),
            lwr_75 = quantile(RE, 0.125),
            upr_75 = quantile(RE, 0.875),
            lwr_95 = quantile(RE, 0.025),
            upr_95 = quantile(RE, 0.975)) %>% 
  filter(Type %in% c("Female Sex Ratio", "Bmsy", "M_F", "M_M", "R0", "Tier 3 HCR Catch")) %>%  # ISS summed for both sexes = 50
  mutate(Type = case_when(
    Type == "Female Sex Ratio" ~ "F SexRatio",
    Type == "Bmsy" ~ "BMSY",
    Type == "M_F" ~ "F NatMort",
    Type == "M_M" ~ "M NatMort",
    Type == "R0" ~ "R0",
    Type == "Tier 3 HCR Catch" ~ "HCR Catch"
  ), Type = factor(Type, levels = c("R0", "BMSY", "F SexRatio", "F NatMort", "M NatMort", "HCR Catch")))

# time series summary
ts_exp3_sum <- exp3_ts_df %>% 
  filter(Convergence == "Converged") %>% 
  mutate(RE = (as.numeric(Pred)-Truth)/Truth) %>% 
  group_by(Years, Type, OM, EM) %>% 
  summarize(Median = median(RE),
            lwr_75 = quantile(RE, 0.125),
            upr_75 = quantile(RE, 0.875),
            lwr_95 = quantile(RE, 0.025),
            upr_95 = quantile(RE, 0.975))


# Figure 5 ----------------------------------------------------------------

# Time series relative error
exp3_ts_plot <- ggplot() +
    geom_ribbon(ts_exp3_sum %>% filter(Type == "Spawning Stock Biomass"),
                mapping = aes(x = Years, ymin = lwr_95, ymax = upr_95, fill = EM),
                alpha = 0.3) +
    geom_ribbon(ts_exp3_sum %>% filter(Type == "Spawning Stock Biomass"),
                mapping = aes(x = Years, ymin = lwr_75, ymax = upr_75, fill = EM),
                alpha = 0.4) +
    geom_line(ts_exp3_sum %>% filter(Type == "Spawning Stock Biomass"),
              mapping = aes(x = Years, y = Median, color = EM), size = 1.8) +
    scale_fill_manual(values = c("#e69c4c", "#35577b")) +
    scale_color_manual(values = c("#e69c4c", "#35577b")) +
    facet_wrap(~OM) +
    geom_hline(yintercept = 0, lty = 2, size = 1.3) + 
    labs(x = "Year", y = "Relative Error", color = "Estimation Model", fill = "Estimation Model") +
    theme_tj() +
    theme(legend.position = "top") 

# Comparing proportions within variants
exp3_param_plot <- ggplot() +
        geom_pointrange(exp3_param_sum,  # 95% quantiles
                        mapping = aes(x = Type, y = Median, ymin = lwr_95, ymax = upr_95, color = EM, fill = EM),
                        position = position_dodge2(width = 0.65), 
                        size = 0, linewidth = 1, alpha = 1) +
        geom_pointrange(exp3_param_sum,  # 75% quantiles
                        mapping = aes(x = Type, y = Median, ymin = lwr_75, ymax = upr_75, color = EM, fill = EM),
                        position = position_dodge2(width = 0.65), 
                        size = 1.5, linewidth = 2, alpha = 0.8) +
        facet_wrap(~OM) +
        scale_fill_manual(values = c("#e69c4c", "#35577b")) +
        scale_color_manual(values = c("#e69c4c", "#35577b")) +
        geom_hline(yintercept = 0, lty = 2, size = 1) + 
        labs(x = "Parameter", y = "Relative Error", color = "Estimation Model", fill = "Estimation Model") +
        theme_tj() +
        theme(legend.position = "none") +
        scale_x_discrete(guide = guide_axis(angle = 90))

pdf(here("figs", "ms_figs", "Fig5_Exp3_SR.pdf"), width = 13, height = 12)
exp3_ts_plot/exp3_param_plot + plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 23, face = "bold"))
dev.off()


# Figure 6 (Experiment 2 Age Structure Bias) ----------------------------------------------------------------

# fishing mortality plot
exp2_fplot <- ggplot() +
  geom_ribbon(ts_exp2_sum %>% filter(Type == "Total Fishing Mortality", EM == "AgeStrc"), # 95% simulation intervals
              mapping = aes(x = Years, y = Median, ymin = lwr_95, ymax = upr_95),
              alpha = 0.3, fill = "#35577b") +
  geom_ribbon(ts_exp2_sum %>% filter(Type == "Total Fishing Mortality", EM == "AgeStrc"), # 75% simulation intervals
              mapping = aes(x = Years, y = Median, ymin = lwr_75, ymax = upr_75),
              alpha = 0.5, fill = "#35577b") +
  geom_line(ts_exp2_sum %>% filter(Type == "Total Fishing Mortality", EM == "AgeStrc"), 
            mapping = aes(x = Years, y = Median), 
            colour = "#35577b", size = 1.3, alpha = 1) +
  facet_grid(OM~., scales = "free") +
  geom_hline(yintercept = 0, lty = 2, size = 1.3) + 
  labs(x = "Year", y = "Relative Error (Total Fishing Mortality)") +
  theme_tj() +
  theme(strip.background = element_blank(),
        strip.text = element_blank())

# selex plot
exp2_selex_plot <- exp2_selex_sum %>% 
    filter(Type == "Fishery Selectivity", EM == "AgeStrc") %>% 
    mutate(Sex = ifelse(Sex == 1, "Female", "Male")) %>% 
    ggplot() +
    geom_line(aes(x = Age, y = Pred, group = sim), color = "grey") +
    geom_line(aes(x = Age, y = True, color = Sex), size = 1) +
    scale_color_manual(values = c("#DC3220", "#005AB5")) +
    facet_grid(OM~Sex) +
    labs(x = "Age", y = "Fishery Selectivity") +
    theme_tj() +
    theme(legend.position = c(0.8, 0.88),
          strip.background.x = element_blank(),
          strip.text.x = element_blank())

pdf(here("figs", "ms_figs", "Fig6_Exp2_FSel_Age.pdf"), width = 13, height = 10)
(exp2_fplot & exp2_selex_plot) + plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 23, face = "bold"))
dev.off()
