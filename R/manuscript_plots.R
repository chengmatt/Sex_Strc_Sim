# Purpose: To plot manuscript-ready figures
# Creator: Matthew LH. Cheng (UAF-CFOS)
# Date: 3/28/24

# Set up ------------------------------------------------------------------

library(here)
library(tidyverse)
library(patchwork)
library(ggtext)

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
n_sims <- 500 # number of simulations run

# Experiment 1 ------------------------------------------------------------
# Read in data for plotting
exp1_ts_df <- data.table::fread(here("output", "Experiment_1_TimeSeries.csv"))  
exp1_cov_df <- data.table::fread(here("output", "Experiment_1_Coverage.csv")) 
exp1_param_df <- data.table::fread(here("output", "Experiment_1_Param.csv")) 
exp1_conv_df <- data.table::fread(here("output", "Experiment_1_Convergence.csv")) 

# Figure out convergence rates
exp1_convergence_diag <- exp1_conv_df %>% 
  filter(str_detect(OM, "_100"),
         convergence == "Converged") %>% 
  group_by(OM, EM) %>% 
  summarize(conv = n() / n_sims) %>% 
  mutate(OM = str_remove(OM, "_100"))

exp1_convergence_diag %>% group_by(EM) %>% summarize(mean(conv))
  
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
  filter(str_detect(OM, "_100")) %>%  # ISS summed for both sexes = 50
  mutate(OM = str_remove(OM, "_100"))

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
  filter(str_detect(OM, "_100"),
         Type %in% c("Bmsy", "M_F", "M_M", "R0", "Tier 3 HCR Catch")) %>%  # ISS summed for both sexes = 50
  mutate(OM = str_remove(OM, "_100"),
         Type = case_when(
         Type == "Bmsy" ~ "BMSY",
         Type == "M_F" ~ "F NatMort",
         Type == "M_M" ~ "M NatMort",
         Type == "R0" ~ "R0",
         Type == "Tier 3 HCR Catch" ~ "HCR Catch"
         ), Type = factor(Type, levels = c( "R0", "BMSY", "F NatMort", "M NatMort", "HCR Catch")))

# Figure 3 (Time Series Rel Err + Coverage + Parameters) ------------------------------------------

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
  scale_fill_manual(values = c("#e69c4c", "#35577b"), labels = c("*Across*", "*Within*")) +
  scale_color_manual(values = c("#e69c4c", "#35577b"), labels = c("*Across*", "*Within*")) +
  coord_cartesian(ylim = c(-0.3, 0.3)) +
  facet_wrap(~OM) +
  geom_hline(yintercept = 0, lty = 2, size = 1.3) + 
  labs(x = "Year", y = "Relative Error (Spawning Stock Biomass)", color = "Integrated Population Model", fill = "Integrated Population Model") +
  theme_tj() +
  theme(legend.position = "top",
        legend.text = element_markdown()) 

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
  scale_fill_manual(values = c("#e69c4c", "#35577b"), labels = c("*Across*", "*Within*")) +
  scale_color_manual(values = c("#e69c4c", "#35577b"), labels = c("*Across*", "*Within*")) +
  geom_hline(yintercept = 0, lty = 2, size = 1) + 
  coord_cartesian(ylim = c(-0.85, 0.85)) +
  labs(x = "", y = "Relative Error (Est & Deriv Pars)", color = "Integrated Population Model", fill = "Integrated Population Model") +
  theme_tj() +
  theme(legend.position = "none",
        legend.text = element_markdown()) +
  scale_x_discrete(guide = guide_axis(angle = 90))

# combine plots
exp1_plot <- exp1_ts_plot/exp1_param_plot + plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 23, face = "bold"))

ggsave(plot = exp1_plot, filename = here("figs", "ms_figs", "Fig3_Exp1_AcrWith.png"), width = 12, height = 13)

# Experiment 2 ------------------------------------------------------------

# Read in files
exp2_param_df <- data.table::fread(here("output", "Experiment_2_Param.csv")) %>% filter(!str_detect(OM, "Growth_M|LargErr"))
exp2_ts_df <- data.table::fread(here("output", "Experiment_2_TimeSeries.csv"))  %>% filter(!str_detect(OM, "Growth_M|LargErr"))
exp2_selex_df <- data.table::fread(here("output", "Experiment_2_Selex.csv")) %>% filter(!str_detect(OM, "Growth_M|LargErr"))
exp2_conv_df <- data.table::fread(here("output", "Experiment_2_Convergence.csv")) %>% filter(!str_detect(OM, "Growth_M|LargErr"))

# Figure out convergence rates
exp2_convergence_diag <- exp2_conv_df %>% 
  filter(!str_detect(EM, "AgeSel|AgeAgg"), convergence == "Converged") %>% 
  group_by(OM, EM) %>% 
  summarize(conv = n() / n_sims) 

exp2_convergence_diag %>% group_by(EM) %>% summarize(mean(conv))

# Set up EM Names here
exp2_labels <- c(
  "SglSx", 
  "MltSx_AggM_AggC",
  "MltSx_AggM_SxC",
  "MltSx_SxM_AggC",
  "MltSx_SxM_SxC"
)

exp2_levels <- c(
  "SglSx_LenSel", 
  "MltSx_AggM_AggC",
  "MltSx_AggM_SxC",
  "MltSx_SxM_AggC",
  "MltSx_SxM_SxC"
)


# summarize parameters
exp2_param_df <- exp2_param_df %>% 
  filter(Convergence == "Converged",
         !str_detect(EM, "AgeSel")) %>% 
  mutate(RE = (as.numeric(Pred) - Truth) / Truth,
         EM = factor(EM, labels = exp2_labels, levels = exp2_levels)) 

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
         !str_detect(EM, "AgeSel")) %>%
  mutate(RE = (as.numeric(Pred) - Truth) / Truth,
         EM = factor(EM, labels = exp2_labels, levels = exp2_levels)) %>% 
  group_by(Years, Type, OM, EM) %>% 
  summarize(Median = median(RE),
            lwr_75 = quantile(RE, 0.125),
            upr_75 = quantile(RE, 0.875),
            lwr_95 = quantile(RE, 0.025),
            upr_95 = quantile(RE, 0.975))

# summarize relative error in selectivity
exp2_selex_sum <- exp2_selex_df %>% 
  filter(Convergence == "Converged",
         !str_detect(EM, "AgeSel")) %>% 
  mutate(EM = factor(EM, labels = exp2_labels, levels = exp2_levels))

# Figure 4 (Spawning Biomass Trends) --------------------------------------

exp2_ssb_plot <- ggplot() +
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
  labs(x = "Year", y = "Relative Error (Spawning Stock Biomass)") +
  theme_tj() +
  theme(strip.text.x = element_text(face = "italic")) +
  coord_cartesian(ylim = c(-0.5, 0.5))

ggsave(plot = exp2_ssb_plot, filename = here("figs", "ms_figs", "Fig4_Exp2_SSBTrends.png"), width = 12.5, height = 11)

# Figure 5 (Experiment 2 Parameters) ----------------------------------------------------------------

exp2_param_plot <- ggplot(exp2_param_sum, 
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
        labs(x = "", y = "Relative Error (Est & Deriv Pars)") +
        theme_tj() +
        scale_x_discrete(guide = guide_axis(angle = 90)) +
        coord_cartesian(ylim = c(-1, 2)) +
        theme(strip.text.x = element_text(face = "italic")) 
  
ggsave(plot = exp2_param_plot, filename = here("figs", "ms_figs", "Fig5_Exp2_Params.png"), width = 12.5, height = 11)

# Experiment 3 ------------------------------------------------------------
exp3_param_df <- data.table::fread(here("output", "Experiment_3_Param.csv")) %>% filter(!str_detect(OM, "No"), EM != 'AgeStrc')
exp3_ts_df <- data.table::fread(here("output", "Experiment_3_TimeSeries.csv"))%>% filter(!str_detect(OM, "No"), EM != 'AgeStrc')
exp3_conv_df <- data.table::fread(here("output", "Experiment_3_Convergence.csv"))

# Figure out convergence rates
exp3_convergence_diag <- exp3_conv_df %>% 
  filter(convergence == "Converged") %>% 
  group_by(OM, EM) %>% 
  summarize(conv = n() / n_sims) 

exp3_convergence_diag %>% group_by(EM) %>% summarize(mean(conv))

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


# Figure 6 ----------------------------------------------------------------

# Time series relative error
exp3_ts_plot <- ggplot() +
    geom_ribbon(ts_exp3_sum %>% filter(Type == "Spawning Stock Biomass",
                                       OM %in% c("Fem40_Mal60", "Fem50_Mal50", "Fem60_Mal40")),
                mapping = aes(x = Years, ymin = lwr_95, ymax = upr_95, fill = EM),
                alpha = 0.3) +
    geom_ribbon(ts_exp3_sum %>% filter(Type == "Spawning Stock Biomass",
                                       OM %in% c("Fem40_Mal60", "Fem50_Mal50", "Fem60_Mal40")),
                mapping = aes(x = Years, ymin = lwr_75, ymax = upr_75, fill = EM),
                alpha = 0.4) +
    geom_line(ts_exp3_sum %>% filter(Type == "Spawning Stock Biomass",
                                     OM %in% c("Fem40_Mal60", "Fem50_Mal50", "Fem60_Mal40")),
              mapping = aes(x = Years, y = Median, color = EM), size = 1.8) +
    scale_fill_manual(values = c("#e69c4c", "#35577b"), labels = c("*EstSxRat*", "*FixSxRat*")) +
    scale_color_manual(values = c("#e69c4c", "#35577b"), labels = c("*EstSxRat*", "*FixSxRat*")) +
    facet_wrap(~OM) +
    coord_cartesian(ylim = c(-0.85, 0.85)) +
    geom_hline(yintercept = 0, lty = 2, size = 1.3) + 
    labs(x = "Year", y = "Relative Error (Spawning Stock Biomass)", color = "Integrated Population Model", fill = "Integrated Population Model") +
    theme_tj() +
    theme(legend.position = "top",
          legend.text = element_markdown())  


# Comparing proportions within variants
exp3_param_plot <- ggplot() +
        geom_pointrange(exp3_param_sum %>% filter(OM %in% c("Fem40_Mal60", "Fem50_Mal50", "Fem60_Mal40")),  # 95% quantiles
                        mapping = aes(x = Type, y = Median, ymin = lwr_95, ymax = upr_95, color = EM, fill = EM),
                        position = position_dodge2(width = 0.65), 
                        size = 0, linewidth = 1, alpha = 1) +
        geom_pointrange(exp3_param_sum %>% filter(OM %in% c("Fem40_Mal60", "Fem50_Mal50", "Fem60_Mal40")),  # 75% quantiles
                        mapping = aes(x = Type, y = Median, ymin = lwr_75, ymax = upr_75, color = EM, fill = EM),
                        position = position_dodge2(width = 0.65), 
                        size = 1.5, linewidth = 2, alpha = 0.8) +
        facet_wrap(~OM) +
        scale_fill_manual(values = c("#e69c4c", "#35577b"), labels = c("*EstSxRat*", "*FixSxRat*")) +
        scale_color_manual(values = c("#e69c4c", "#35577b"), labels = c("*EstSxRat*", "*FixSxRat*")) +
        geom_hline(yintercept = 0, lty = 2, size = 1) + 
        coord_cartesian(ylim = c(-1.5, 1.5)) +
        labs(x = "", y = "Relative Error (Est & Deriv Pars)", color = "Integrated Population Model", fill = "Integrated Population Model") +
        theme_tj() +
        theme(legend.position = "none",
              legend.text = element_markdown()) +
        scale_x_discrete(guide = guide_axis(angle = 90))

exp3_comb_plot <- exp3_ts_plot/exp3_param_plot + plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 23, face = "bold"))

ggsave(plot = exp3_comb_plot, filename = here("figs", "ms_figs", "Fig6_Exp3_SR.png"), width = 12, height = 13)

# Figure S1 (Exp 1 Total Biomass) ------------------------------------------------------

# Time series relative error
exp1_ts_total_biom_plot <- ggplot() +
  geom_ribbon(ts_exp1_sum %>% filter(Type == "Total Biomass"),
              mapping = aes(x = Years, ymin = lwr_95, ymax = upr_95, fill = EM),
              alpha = 0.3) +
  geom_ribbon(ts_exp1_sum %>% filter(Type == "Total Biomass"),
              mapping = aes(x = Years, ymin = lwr_75, ymax = upr_75, fill = EM),
              alpha = 0.4) +
  geom_line(ts_exp1_sum %>% filter(Type == "Spawning Stock Biomass"),
            mapping = aes(x = Years, y = Median, color = EM), size = 1.8) +
  scale_fill_manual(values = c("#e69c4c", "#35577b"), labels = c("*Across*", "*Within*")) +
  scale_color_manual(values = c("#e69c4c", "#35577b"), labels = c("*Across*", "*Within*")) +
  coord_cartesian(ylim = c(-0.3, 0.3)) +
  facet_wrap(~OM) +
  geom_hline(yintercept = 0, lty = 2, size = 1.3) + 
  labs(x = "Year", y = "Relative Error (Total Biomass)", color = "Integrated Population Model", fill = "Integrated Population Model") +
  theme_tj() +
  theme(legend.position = "top",
        legend.text = element_markdown()) 

ggsave(plot = exp1_ts_total_biom_plot, width = 8, height = 5,
       filename = here("figs", "ms_figs", "FigS1_Exp1_TotalBiom.png"))


# Figure S2 (Exp 2 Total Biomass) -----------------------------------------

exp2_total_biom_plot <- ggplot() +
  geom_ribbon(ts_exp2_sum %>% filter(Type == "Total Biomass"), # 95% simulation intervals
              mapping = aes(x = Years, y = Median, ymin = lwr_95, ymax = upr_95),
              alpha = 0.3, fill = "#35577b") +
  geom_ribbon(ts_exp2_sum %>% filter(Type == "Total Biomass"), # 75% simulation intervals
              mapping = aes(x = Years, y = Median, ymin = lwr_75, ymax = upr_75),
              alpha = 0.5, fill = "#35577b") +
  geom_line(ts_exp2_sum %>% filter(Type == "Total Biomass"), 
            mapping = aes(x = Years, y = Median), 
            colour = "#35577b", size = 1.3, alpha = 1) +
  facet_grid(OM~EM, scales = "free") +
  geom_hline(yintercept = 0, lty = 2, size = 1.3) + 
  labs(x = "Year", y = "Relative Error (Total Biomass)") +
  theme_tj() +
  theme(strip.text.x = element_text(face = "italic")) +
  coord_cartesian(ylim = c(-0.5, 0.5)) 
  
ggsave(plot = exp2_total_biom_plot,  width = 12.5, height = 11,
       filename = here("figs", "ms_figs", "FigS2_Exp2_TotalBiomTrends.png"))

# Figure S3 (Exp2, AgeStrc, Avg Growth and Selex) -------------------------------
# Growth
exp2_growth_df <- data.table::fread(here("output", "Experiment_2_Growth.csv")) %>% 
  filter(Convergence == "Converged") %>% 
  filter(!str_detect(EM, "AgeSel|AgeAgg")) %>% 
  mutate(RE = (as.numeric(Pred) - True) / True,
         EM = factor(EM, labels = exp2_labels, levels = exp2_levels)) %>% 
  group_by(OM, EM, Sex, Age, Type) %>% 
  summarize(Median_Pred = median(Pred),
            Median_True = median(True),
            lwr_95 = quantile(Pred, 0.025),
            upr_95 = quantile(Pred, 0.975))

# Selex
exp2_selex_df <- data.table::fread(here("output", "Experiment_2_Selex.csv")) %>% 
  filter(Convergence == "Converged") %>% 
  filter(!str_detect(EM, "AgeSel|AgeAgg")) %>% 
  mutate(RE = (as.numeric(Pred) - True) / True,
         EM = factor(EM, labels = exp2_labels, levels = exp2_levels)) %>% 
  group_by(OM, EM, Sex, Age, Type) %>% 
  summarize(Median_Pred = median(Pred),
            Median_True = median(True),
            lwr_95 = quantile(Pred, 0.025),
            upr_95 = quantile(Pred, 0.975))

# Growth plots
exp2_laa_age_plot <- exp2_growth_df %>% 
  filter(EM == "SglSx", str_detect(Type, "Length")) %>% 
  mutate(Sex = ifelse(Sex == 1, "Female", "Male")) %>% 
  ggplot() +
  geom_line(aes(x = Age, y = Median_Pred), color = 'black', alpha = 1, lwd = 1) +
  geom_ribbon(aes(x = Age, y = Median_Pred, ymin = lwr_95, ymax = upr_95), alpha = 0.4) +
  geom_line(aes(x = Age, y = Median_True, color = factor(Sex)), size = 1.3, alpha = 0.85) +
  scale_color_manual(values = c("#DC3220", "#005AB5")) +
  facet_grid(OM~.) +
  labs(x = "Age", y = "Length-at-age (cm)", color = "Sex") +
  theme_tj() +
  theme(strip.background = element_blank(),
        strip.text = element_blank())

exp2_waa_age_plot <- exp2_growth_df %>% 
  filter(EM == "SglSx",
         Type == "Weight at age") %>% 
  mutate(Sex = ifelse(Sex == 1, "Female", "Male")) %>% 
  ggplot() +
  geom_line(aes(x = Age, y = Median_Pred), color = 'black', alpha = 1, lwd = 1) +
  geom_ribbon(aes(x = Age, y = Median_Pred, ymin = lwr_95, ymax = upr_95), alpha = 0.4) +
  geom_line(aes(x = Age, y = Median_True, color = factor(Sex)), size = 1.3, alpha = 0.85) +
  scale_color_manual(values = c("#DC3220", "#005AB5")) +
  facet_grid(OM~.) +
  labs(x = "Age", y = "Weight-at-age (g)", color = "Sex") +
  theme_tj() +
  theme(strip.background = element_blank(),
        strip.text = element_blank())

# Survey selectivity plot
exp2_srv_selex_age_plot <- exp2_selex_df %>% 
  filter(Type == "Survey Selectivity",
         EM == "SglSx") %>% 
  mutate(Sex = ifelse(Sex == 1, "Female", "Male")) %>% 
  ggplot() +
  geom_line(aes(x = Age, y = Median_Pred), color = 'black', alpha = 1, lwd = 1) +
  geom_ribbon(aes(x = Age, y = Median_Pred, ymin = lwr_95, ymax = upr_95), alpha = 0.4) +
  geom_line(aes(x = Age, y = Median_True, color = factor(Sex)), size = 1.3, alpha = 0.85) +
  scale_color_manual(values = c("#DC3220", "#005AB5")) +
  facet_grid(OM~.) +
  labs(x = "Age", y = "Survey Selectivity", color = "Sex") +
  ylim(0, NA) +
  theme_tj() +
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        legend.position = "top")

# fishery selectivity plot
exp2_fish_selex_age_plot <- exp2_selex_df %>% 
  filter(Type == "Fishery Selectivity",
         EM == "SglSx") %>% 
  mutate(Sex = ifelse(Sex == 1, "Female", "Male")) %>% 
  ggplot() +
  geom_line(aes(x = Age, y = Median_Pred), color = 'black', alpha = 1, lwd = 1) +
  geom_ribbon(aes(x = Age, y = Median_Pred, ymin = lwr_95, ymax = upr_95), alpha = 0.4) +
  geom_line(aes(x = Age, y = Median_True, color = factor(Sex)), size = 1.3, alpha = 0.85) +
  scale_color_manual(values = c("#DC3220", "#005AB5")) +
  facet_grid(OM~.) +
  ylim(0, NA) +
  labs(x = "Age", y = "Fishery Selectivity", color = "Sex") +
  theme_tj()

comb_exp2_age_plot <- ggpubr::ggarrange(exp2_laa_age_plot,
                                        exp2_waa_age_plot, 
                                        exp2_srv_selex_age_plot, 
                                        exp2_fish_selex_age_plot, common.legend = TRUE, ncol = 4)

ggsave(plot = comb_exp2_age_plot, width = 17,
       filename = here("figs", "ms_figs", "FigS3_Exp2.png"))

# Figure S4 (Exp 3 Total Biomass) -----------------------------------------

# Time series relative error
exp3_ts_total_biom_plot <- ggplot() +
  geom_ribbon(ts_exp3_sum %>% filter(Type == "Total Biomass", OM %in% c("Fem40_Mal60", "Fem50_Mal50", "Fem60_Mal40")),
              mapping = aes(x = Years, ymin = lwr_95, ymax = upr_95, fill = EM),
              alpha = 0.3) +
  geom_ribbon(ts_exp3_sum %>% filter(Type == "Total Biomass", OM %in% c("Fem40_Mal60", "Fem50_Mal50", "Fem60_Mal40")),
              mapping = aes(x = Years, ymin = lwr_75, ymax = upr_75, fill = EM),
              alpha = 0.4) +
  geom_line(ts_exp3_sum %>% filter(Type == "Total Biomass", OM %in% c("Fem40_Mal60", "Fem50_Mal50", "Fem60_Mal40")),
            mapping = aes(x = Years, y = Median, color = EM), size = 1.8) +
  scale_fill_manual(values = c("#e69c4c", "#35577b"), labels = c("*EstSxRat*", "*FixSxRat*")) +
  scale_color_manual(values = c("#e69c4c", "#35577b"), labels = c("*EstSxRat*", "*FixSxRat*")) +
  facet_wrap(~OM) +
  coord_cartesian(ylim = c(-0.85, 0.85)) +
  geom_hline(yintercept = 0, lty = 2, size = 1.3) + 
  labs(x = "Year", y = "Relative Error (Total Biomass)", color = "Integrated Population Model", fill = "Integrated Population Model") +
  theme_tj() +
  theme(legend.position = "top",
        legend.text = element_markdown()) 

ggsave(plot = exp3_ts_total_biom_plot, width = 10, height = 5,
       filename = here("figs", "ms_figs", "FigS4_Exp3_TotalBiom.png"))

# Figure S5 (Exp 5 Total Fishing Mortality) -----------------------------------------

# Time series relative error
exp3_ts_total_f_plot <- ggplot() +
  geom_ribbon(ts_exp3_sum %>% filter(Type == "Total Fishing Mortality", OM %in% c("Fem40_Mal60", "Fem50_Mal50", "Fem60_Mal40")),
              mapping = aes(x = Years, ymin = lwr_95, ymax = upr_95, fill = EM),
              alpha = 0.3) +
  geom_ribbon(ts_exp3_sum %>% filter(Type == "Total Fishing Mortality", OM %in% c("Fem40_Mal60", "Fem50_Mal50", "Fem60_Mal40")),
              mapping = aes(x = Years, ymin = lwr_75, ymax = upr_75, fill = EM),
              alpha = 0.4) +
  geom_line(ts_exp3_sum %>% filter(Type == "Total Fishing Mortality", OM %in% c("Fem40_Mal60", "Fem50_Mal50", "Fem60_Mal40")),
            mapping = aes(x = Years, y = Median, color = EM), size = 1.8) +
  scale_fill_manual(values = c("#e69c4c", "#35577b"), labels = c("*EstSxRat*", "*FixSxRat*")) +
  scale_color_manual(values = c("#e69c4c", "#35577b"), labels = c("*EstSxRat*", "*FixSxRat*")) +
  facet_wrap(~OM) +
  coord_cartesian(ylim = c(-0.85, 0.85)) +
  geom_hline(yintercept = 0, lty = 2, size = 1.3) + 
  labs(x = "Year", y = "Relative Error (Total Fishing Mortality)", color = "Integrated Population Model", fill = "Integrated Population Model") +
  theme_tj() +
  theme(legend.position = "top",
        legend.text = element_markdown()) 

ggsave(plot = exp3_ts_total_f_plot, width = 10, height = 7,
       filename = here("figs", "ms_figs", "FigS5_Exp3_TotalF.png"))

 
# Figure S6 (Exp 3 Senstivity Eexploration) -----------------------------------------

# Time series relative error
exp3_senstivity_ts_plot <- ggplot() +
  geom_ribbon(ts_exp3_sum %>% filter(Type == "Spawning Stock Biomass", 
                                     !OM %in% c("Fem40_Mal60", "Fem50_Mal50", "Fem60_Mal40")),
              mapping = aes(x = Years, ymin = lwr_95, ymax = upr_95, fill = EM),
              alpha = 0.3) +
  geom_ribbon(ts_exp3_sum %>% filter(Type == "Spawning Stock Biomass", 
                                     !OM %in% c("Fem40_Mal60", "Fem50_Mal50", "Fem60_Mal40")),
              mapping = aes(x = Years, ymin = lwr_75, ymax = upr_75, fill = EM),
              alpha = 0.4) +
  geom_line(ts_exp3_sum %>% filter(Type == "Spawning Stock Biomass", 
                                   !OM %in% c("Fem40_Mal60", "Fem50_Mal50", "Fem60_Mal40")),
            mapping = aes(x = Years, y = Median, color = EM), size = 1.8) +
  scale_fill_manual(values = c("#e69c4c", "#35577b"), labels = c("*EstSxRat*", "*FixSxRat*")) +
  scale_color_manual(values = c("#e69c4c", "#35577b"), labels = c("*EstSxRat*", "*FixSxRat*")) +
  facet_wrap(~OM, nrow = 1) +
  coord_cartesian(ylim = c(-0.85, 0.85)) +
  geom_hline(yintercept = 0, lty = 2, size = 1.3) + 
  labs(x = "Year", y = "Relative Error", color = "Integrated Population Model", fill = "Integrated Population Model") +
  theme_tj() +
  theme(legend.position = "top",
        legend.text = element_markdown()) 

# Comparing proportions within variants
exp3_sensitivity_param_plot <- ggplot() +
  geom_pointrange(exp3_param_sum %>% filter(!OM %in% c("Fem40_Mal60", "Fem50_Mal50", "Fem60_Mal40")),  # 95% quantiles
                  mapping = aes(x = Type, y = Median, ymin = lwr_95, ymax = upr_95, color = EM, fill = EM),
                  position = position_dodge2(width = 0.65), 
                  size = 0, linewidth = 1, alpha = 1) +
  geom_pointrange(exp3_param_sum %>% filter(!OM %in% c("Fem40_Mal60", "Fem50_Mal50", "Fem60_Mal40")),  # 75% quantiles
                  mapping = aes(x = Type, y = Median, ymin = lwr_75, ymax = upr_75, color = EM, fill = EM),
                  position = position_dodge2(width = 0.65), 
                  size = 1.5, linewidth = 2, alpha = 0.8) +
  facet_wrap(~OM, nrow = 1) +
  scale_fill_manual(values = c("#e69c4c", "#35577b"), labels = c("*EstSxRat*", "*FixSxRat*")) +
  scale_color_manual(values = c("#e69c4c", "#35577b"), labels = c("*EstSxRat*", "*FixSxRat*")) +
  geom_hline(yintercept = 0, lty = 2, size = 1) + 
  coord_cartesian(ylim = c(-1.5, 1.5)) +
  labs(x = "", y = "Relative Error (Est & Deriv Pars)", color = "Integrated Population Model", fill = "Integrated Population Model") +
  theme_tj() +
  theme(legend.position = "none",
        legend.text = element_markdown()) +
  scale_x_discrete(guide = guide_axis(angle = 90))

exp3_sens_comb_plot <- exp3_senstivity_ts_plot/exp3_sensitivity_param_plot + plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 23, face = "bold"))

ggsave(plot = exp3_sens_comb_plot, width = 15, height = 13, 
       filename = here("figs", "ms_figs", "FigS6_Exp3_Senstivity.png"))
