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
  mutate(OM = str_remove(OM, "_100"),
         EM = ifelse(EM == 'Across', 'Joint', 'Split'),
         OM = ifelse(OM == 'Across', 'Joint', 'Split'))

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
         Type %in% c("Bmsy", "M_F", "M_M", "R0", "Tier 3 HCR Catch", "SSBcur / BMSY")) %>%  # ISS summed for both sexes = 50
  mutate(OM = str_remove(OM, "_100"),
         Type = case_when(
         Type == "Bmsy" ~ "BMSY",
         Type == "M_F" ~ "F NatMort",
         Type == "M_M" ~ "M NatMort",
         Type == "R0" ~ "R0",
         Type == "Tier 3 HCR Catch" ~ "HCR Catch",
         Type == "SSBcur / BMSY" ~ "SSBcur / BMSY"
         ), Type = factor(Type, levels = c( "R0", "BMSY", "SSBcur / BMSY", "F NatMort", "M NatMort", "HCR Catch")),
         EM = ifelse(EM == 'Across', 'Joint', 'Split'),
         OM = ifelse(OM == 'Across', 'Joint', 'Split'))

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
  scale_fill_manual(values = c("#e69c4c", "#35577b"), labels = c("*Joint*", "*Split*")) +
  scale_color_manual(values = c("#e69c4c", "#35577b"), labels = c("*Joint*", "*Split*")) +
  coord_cartesian(ylim = c(-0.3, 0.3)) +
  facet_wrap(~OM) +
  geom_hline(yintercept = 0, lty = 2, size = 1.3) + 
  labs(x = "Year", y = "Relative Error (Spawning Stock Biomass)", color = "Integrated Population Model", fill = "Integrated Population Model") +
  theme_tj() +
  theme(legend.position = "top",
        legend.text = element_markdown()) 

ggsave(here("figs", "Presentation_Figures", "Exp1_SSB.png"),
       exp1_ts_plot + labs(color = "EM", fill = "EM") + coord_cartesian(ylim = c(-0.5, 0.5)) + theme(strip.text = element_text(size = 17)))

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
  filter(Type %in% c("Bmsy", "M_F", "M_M", "R0", "Tier 3 HCR Catch", "SSBcur / BMSY")) %>%  # ISS summed for both sexes = 50
  mutate(OM = str_remove(OM, "_100"),
         Type = case_when(
           Type == "Bmsy" ~ "BMSY",
           Type == "M_F" ~ "F NatMort",
           Type == "M_M" ~ "M NatMort",
           Type == "R0" ~ "R0",
           Type == "Tier 3 HCR Catch" ~ "HCR Catch",
           Type == "SSBcur / BMSY" ~ "SSBcur / BMSY"
         ), Type = factor(Type, levels = c( "R0", "BMSY", "SSBcur / BMSY", "F NatMort", "M NatMort", "HCR Catch")))

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

# Save for presentation
ggsave(here("figs", "Presentation_Figures", "Exp2_SSB.png"), 
       ggplot(ts_exp2_sum %>% filter(Type == "Spawning Stock Biomass",!str_detect(EM, "SxC")) %>% 
                mutate(EM = str_remove(EM, "_AggC"),
                       EM = factor(EM, levels = c("SglSx", "MltSx_AggM", "MltSx_SxM")))) +
         geom_ribbon(mapping = aes(x = Years, y = Median, ymin = lwr_95, ymax = upr_95),
                     alpha = 0.3, fill = "#35577b") +
         geom_ribbon(mapping = aes(x = Years, y = Median, ymin = lwr_75, ymax = upr_75),
                     alpha = 0.5, fill = "#35577b") +
         geom_line(mapping = aes(x = Years, y = Median), 
                   colour = "#35577b", size = 1.3, alpha = 1) +
         facet_grid(EM~OM, scales = "free") +
         geom_hline(yintercept = 0, lty = 2, size = 1.3) + 
         labs(x = "Year", y = "Relative Error (Spawning Stock Biomass)", color = "EM", fill = "EM") +
         theme_tj() +
         theme(strip.text.y = element_text(face = "italic"),
               strip.text = element_text(size = 17)) +
         coord_cartesian(ylim = c(-0.5, 0.5))
)

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
  filter(Type %in% c("Female Sex Ratio", "Bmsy", "M_F", "M_M", "R0", "Tier 3 HCR Catch", "SSBcur / BMSY")) %>%  # ISS summed for both sexes = 50
  mutate(Type = case_when(
    Type == "Female Sex Ratio" ~ "F SexRatio",
    Type == "Bmsy" ~ "BMSY",
    Type == "M_F" ~ "F NatMort",
    Type == "M_M" ~ "M NatMort",
    Type == "R0" ~ "R0",
    Type == "Tier 3 HCR Catch" ~ "HCR Catch",
    Type == "SSBcur / BMSY" ~ "SSBcur / BMSY"
  ), Type = factor(Type, levels = c("R0", "BMSY", "SSBcur / BMSY", "F SexRatio", "F NatMort", "M NatMort", "HCR Catch")))

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

ggsave(here("figs", "Presentation_Figures", "Exp3_SSB.png"),
       exp3_ts_plot + labs(color = "EM", fill = "EM") + theme(strip.text = element_text(size = 17)))

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


# Figure 7 (Sex-Ratio Demonstration) --------------------------------------
load(here("output", 'Experiment 3', "Fem60_Mal40", "Fem60_Mal40.RData")) # load in OM

# Get composition data
fish_agecomps = reshape2::melt(oms$Fish_AgeComps) %>% mutate(Type = 'Fishery Ages')
names(fish_agecomps) = c("Years", "Age", "Sex", "Fleet", 'Sim', "Numbers", 'Type')
srv_agecomps = reshape2::melt(oms$Srv_AgeComps) %>% mutate(Type = 'Survey Ages')
names(srv_agecomps) = c("Years", "Age", "Sex", "Fleet", 'Sim', "Numbers", 'Type')
comp_store = rbind(fish_agecomps, srv_agecomps) # bind together

# Get Sex ratio information over time
sr_store = comp_store %>% 
  group_by(Years, Sim, Type) %>% 
  mutate(Total_N = sum(Numbers)) %>% # get total numbers
  ungroup() %>% 
  filter(Total_N != 0) %>% 
  group_by(Years, Sim, Sex, Type) %>% 
  summarize(Total_Sex = sum(Numbers) / Total_N) %>% 
  unique() %>% 
  group_by(Years, Sex, Type) %>% 
  summarize(median = median(Total_Sex),
            lwr_95 = quantile(Total_Sex, 0.025),
            upr_95 = quantile(Total_Sex, 0.975))

# Get true population sex ratio
naa = reshape2::melt(oms$NAA) %>% mutate(Type = 'Population Sex Ratio')
names(naa) = c("Years", "Age", "Sex", 'Sim', "Numbers", 'Type')

# Get Sex ratio information over time
naa_store = naa %>% 
  group_by(Years, Sim, Type) %>% 
  mutate(Total_N = sum(Numbers)) %>% # get total numbers
  ungroup() %>% 
  filter(Total_N != 0) %>% 
  group_by(Years, Sim, Sex, Type) %>% 
  summarize(Total_Sex = sum(Numbers) / Total_N) %>% 
  unique() %>% 
  group_by(Years, Sex, Type) %>% 
  summarize(median = median(Total_Sex),
            lwr_95 = quantile(Total_Sex, 0.025),
            upr_95 = quantile(Total_Sex, 0.975))

sr_all = rbind(naa_store, sr_store)

# sex ratio across time summary
sr_plot = ggplot(sr_all %>% mutate(Type = factor(Type, levels = c("Population Sex Ratio", "Fishery Ages", "Survey Ages"))) %>% 
                   filter(Type != 'Population Sex Ratio'), 
                 aes(x = Years, y = median, color = factor(Sex), fill = factor(Sex), ymin = lwr_95, ymax = upr_95)) +
  geom_line(size = 1.5) +
  geom_ribbon(alpha = 0.35, color = NA) +
  scale_color_manual(values = c("#DC3220", "#005AB5")) +
  scale_fill_manual(values = c("#DC3220", "#005AB5")) +
  facet_wrap(~Type) +
  labs(x = "Year", y = "Sex Ratio", color = "Sex", fill = "Sex") +
  theme_tj() +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(plot = sr_plot, filename = here("figs", "ms_figs", "Fig7_SexRatio.png"), width = 10, height = 5)

# Weight at age
waa_df = reshape2::melt(oms$waa) %>% mutate(Type = 'Weight-at-age')
names(waa_df) = c("Age", "Sex", "Value", 'Type')
# Age-based selectivity - fishery
fishageselex_df = reshape2::melt(oms$FishAge_Selex) %>% select(-Var3) %>% mutate(Type = 'Fishery Selectivity')
names(fishageselex_df) = c("Age", "Sex", "Value", 'Type')
grw_sel = rbind(waa_df, fishageselex_df)

# fishery selex age plot
grw_sel_plot = ggplot(grw_sel, aes(x = Age, y = Value, color = factor(Sex))) +
  geom_line(size = 2, alpha = 0.8) +
  facet_wrap(~Type, scales = 'free') +
  scale_color_manual(values = c("#DC3220", "#005AB5")) +
  labs(x = "Age", y = "Value", color = "Sex") +
  theme_tj() +
  theme(legend.position = "none",
        plot.background = element_rect(fill = "transparent", colour = NA),
        panel.background = element_rect(fill = "transparent", colour = NA),
        legend.text = element_markdown())  

ggsave(plot = grw_sel_plot, filename = here("figs", "ms_figs", "Fig7_GrwSel.png"), width = 10, height = 5)

# Time series relative error
exp3_biom_plot <- ggplot() +
  geom_ribbon(ts_exp3_sum %>% filter(Type == "Spawning Stock Biomass", EM == 'FixSR', OM %in% c("Fem60_Mal40")),
              mapping = aes(x = Years, ymin = lwr_95, ymax = upr_95),
              alpha = 0.3) +
  geom_ribbon(ts_exp3_sum %>% filter(Type == "Spawning Stock Biomass", EM == 'FixSR', OM %in% c("Fem60_Mal40")),
              mapping = aes(x = Years, ymin = lwr_75, ymax = upr_75),
              alpha = 0.4) +
  geom_line(ts_exp3_sum %>% filter(Type == "Spawning Stock Biomass", EM == 'FixSR', OM %in% c("Fem60_Mal40")),
            mapping = aes(x = Years, y = Median), size = 1.8) +
  coord_cartesian(ylim = c(-0.85, 0.85)) +
  geom_hline(yintercept = 0, lty = 2, size = 1.3) + 
  labs(x = "Year", y = "Relative Error (SSB)") +
  theme_tj() +
  theme(legend.position = "top",
        plot.background = element_rect(fill = "transparent", colour = NA),
        panel.background = element_rect(fill = "transparent", colour = NA),
        legend.text = element_markdown())  

ggsave(plot = exp3_biom_plot, filename = here("figs", "ms_figs", "Fig7_RE.png"), width = 10, height = 5)


# Table of "Lost Yield or OverHarvest" ------------------------------------

# Compute potential yield difference compared to "best" model
yield_diff_table <- rbind(
  exp1_param_sum %>% filter(Type == "HCR Catch") %>% select(OM, EM, Median) %>% mutate(Exp = 1), 
  exp2_param_sum %>% filter(Type == "HCR Catch") %>% select(OM, EM, Median) %>% mutate(Exp = 2), 
  exp3_param_sum %>% filter(Type == "HCR Catch") %>% select(OM, EM, Median) %>% mutate(Exp = 3)
) %>% 
  filter(!str_detect(OM, "Same")) 

write.csv(yield_diff_table, file = here("output", "yield_bias.csv"))

# Figure S1 (Exp 1 Total Biomass) ------------------------------------------------------

# Time series relative error
exp1_ts_total_biom_plot <- ggplot() +
  geom_ribbon(ts_exp1_sum %>% filter(Type == "Total Biomass"),
              mapping = aes(x = Years, ymin = lwr_95, ymax = upr_95, fill = EM),
              alpha = 0.3) +
  geom_ribbon(ts_exp1_sum %>% filter(Type == "Total Biomass"),
              mapping = aes(x = Years, ymin = lwr_75, ymax = upr_75, fill = EM),
              alpha = 0.4) +
  geom_line(ts_exp1_sum %>% filter(Type == "Total Biomass"),
            mapping = aes(x = Years, y = Median, color = EM), size = 1.8) +
  scale_fill_manual(values = c("#e69c4c", "#35577b"), labels = c("*Joint*", "*Split*")) +
  scale_color_manual(values = c("#e69c4c", "#35577b"), labels = c("*Joint*", "*Split*")) +
  coord_cartesian(ylim = c(-0.3, 0.3)) +
  facet_wrap(~OM) +
  geom_hline(yintercept = 0, lty = 2, size = 1.3) + 
  labs(x = "Year", y = "Relative Error (Total Biomass)", color = "Integrated Population Model", fill = "Integrated Population Model") +
  theme_tj() +
  theme(legend.position = "top",
        legend.text = element_markdown()) 

ggsave(plot = exp1_ts_total_biom_plot, width = 8, height = 5,
       filename = here("figs", "ms_figs", "FigS1_Exp1_TotalBiom.png"))

# Figure S2 (Exp 1 Depletion) ------------------------------------------------------

# Time series relative error
exp1_depl_plot <- ggplot() +
  geom_ribbon(ts_exp1_sum %>% filter(Type == "Depletion"),
              mapping = aes(x = Years, ymin = lwr_95, ymax = upr_95, fill = EM),
              alpha = 0.3) +
  geom_ribbon(ts_exp1_sum %>% filter(Type == "Depletion"),
              mapping = aes(x = Years, ymin = lwr_75, ymax = upr_75, fill = EM),
              alpha = 0.4) +
  geom_line(ts_exp1_sum %>% filter(Type == "Depletion"),
            mapping = aes(x = Years, y = Median, color = EM), size = 1.8) +
  scale_fill_manual(values = c("#e69c4c", "#35577b"), labels = c("*Joint*", "*Split*")) +
  scale_color_manual(values = c("#e69c4c", "#35577b"), labels = c("*Joint*", "*Split*")) +
  coord_cartesian(ylim = c(-0.3, 0.3)) +
  facet_wrap(~OM) +
  geom_hline(yintercept = 0, lty = 2, size = 1.3) + 
  labs(x = "Year", y = "Relative Error (Depletion)", color = "Integrated Population Model", fill = "Integrated Population Model") +
  theme_tj() +
  theme(legend.position = "top",
        legend.text = element_markdown()) 

ggsave(plot = exp1_depl_plot, width = 8, height = 5,
       filename = here("figs", "ms_figs", "FigS2_Exp1_Depletion.png"))


# Figure S3 (Exp 2 Total Biomass) -----------------------------------------

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
       filename = here("figs", "ms_figs", "FigS3_Exp2_TotalBiomTrends.png"))

# Figure S4 (Exp 2 Depletion) -----------------------------------------

exp2_depl_plot <- ggplot() +
  geom_ribbon(ts_exp2_sum %>% filter(Type == "Depletion"), # 95% simulation intervals
              mapping = aes(x = Years, y = Median, ymin = lwr_95, ymax = upr_95),
              alpha = 0.3, fill = "#35577b") +
  geom_ribbon(ts_exp2_sum %>% filter(Type == "Depletion"), # 75% simulation intervals
              mapping = aes(x = Years, y = Median, ymin = lwr_75, ymax = upr_75),
              alpha = 0.5, fill = "#35577b") +
  geom_line(ts_exp2_sum %>% filter(Type == "Depletion"), 
            mapping = aes(x = Years, y = Median), 
            colour = "#35577b", size = 1.3, alpha = 1) +
  facet_grid(OM~EM, scales = "free") +
  geom_hline(yintercept = 0, lty = 2, size = 1.3) + 
  labs(x = "Year", y = "Relative Error (Depletion)") +
  theme_tj() +
  theme(strip.text.x = element_text(face = "italic")) +
  coord_cartesian(ylim = c(-0.5, 0.5)) 

ggsave(plot = exp2_depl_plot,  width = 12.5, height = 11,
       filename = here("figs", "ms_figs", "FigS4_Exp2_Depletion.png"))

# Figure S5 (Exp2, AgeStrc, Avg Growth and Selex) -------------------------------
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
       filename = here("figs", "ms_figs", "FigS5_Exp2.png"))

# Figure S6 (Exp 3 Total Biomass) -----------------------------------------

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
       filename = here("figs", "ms_figs", "FigS6_Exp3_TotalBiom.png"))

# Figure S7 (Exp 3 Depletion) -----------------------------------------

# Time series relative error
exp3_ts_total_biom_plot <- ggplot() +
  geom_ribbon(ts_exp3_sum %>% filter(Type == "Depletion", OM %in% c("Fem40_Mal60", "Fem50_Mal50", "Fem60_Mal40")),
              mapping = aes(x = Years, ymin = lwr_95, ymax = upr_95, fill = EM),
              alpha = 0.3) +
  geom_ribbon(ts_exp3_sum %>% filter(Type == "Depletion", OM %in% c("Fem40_Mal60", "Fem50_Mal50", "Fem60_Mal40")),
              mapping = aes(x = Years, ymin = lwr_75, ymax = upr_75, fill = EM),
              alpha = 0.4) +
  geom_line(ts_exp3_sum %>% filter(Type == "Depletion", OM %in% c("Fem40_Mal60", "Fem50_Mal50", "Fem60_Mal40")),
            mapping = aes(x = Years, y = Median, color = EM), size = 1.8) +
  scale_fill_manual(values = c("#e69c4c", "#35577b"), labels = c("*EstSxRat*", "*FixSxRat*")) +
  scale_color_manual(values = c("#e69c4c", "#35577b"), labels = c("*EstSxRat*", "*FixSxRat*")) +
  facet_wrap(~OM) +
  # coord_cartesian(ylim = c(-0.85, 0.85)) +
  geom_hline(yintercept = 0, lty = 2, size = 1.3) + 
  labs(x = "Year", y = "Relative Error (Depletion)", color = "Integrated Population Model", fill = "Integrated Population Model") +
  theme_tj() +
  theme(legend.position = "top",
        legend.text = element_markdown()) 

ggsave(plot = exp3_ts_total_biom_plot, width = 10, height = 5,
       filename = here("figs", "ms_figs", "FigS7_Exp3_TotalBiom.png"))

# Figure S8 (Exp 5 Total Fishing Mortality) -----------------------------------------

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
       filename = here("figs", "ms_figs", "FigS8_Exp3_TotalF.png"))

 
# Figure S9 (Exp 3 Senstivity Eexploration) -----------------------------------------

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
       filename = here("figs", "ms_figs", "FigS9_Exp3_Senstivity.png"))
