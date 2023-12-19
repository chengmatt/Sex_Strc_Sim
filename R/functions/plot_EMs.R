#' Title To plot performance of individual EMs
#'
#' @param time_series Dataframe for time seris
#' @param NAA_sexratio Df for sex ratios
#' @param growth_df Df for growth
#' @param selex_df Df for selectivifty
#' @param params Df for parameters
#'
#' @return
#' @export
#'
#' @examples
plot_EMs = function(time_series, 
                    NAA_sexratio,
                    growth_df,
                    selex_df,
                    params,
                    path
                    ) {
  
  pdf(here(path, "EM_Summary.pdf"))

# Time Series -------------------------------------------------------------

  time_series_sum = time_series %>% 
    mutate(RE = (Pred - Truth) / Truth) %>% 
    group_by(Years, Type) %>% 
    summarize(median = median(RE),
              lwr_95 = quantile(RE, 0.025),
              upr_95 = quantile(RE, 0.975))
  
  print(
    ggplot(time_series_sum, aes(x = Years, y = median,
                                ymin = lwr_95, ymax = upr_95)) +
      geom_line(size = 1.5) +
      geom_ribbon(alpha = 0.5) +
      geom_hline(yintercept = 0, lty = 2) +
      facet_wrap(~Type) 
      # coord_cartesian(ylim = c(-0.75, 0.75))
  )
  

# Sex Ratio ---------------------------------------------------------------

  NAA_sexratio_sum = NAA_sexratio %>% 
    mutate(RE = (Pred - Truth) / Truth) %>% 
    group_by(Years, Type, Age) %>% 
    summarize(median = median(RE),
              lwr_95 = quantile(RE, 0.025),
              upr_95 = quantile(RE, 0.975))  

  print(
    ggplot(NAA_sexratio_sum, aes(x = Age, y = median,
                                ymin = lwr_95, ymax = upr_95)) +
      geom_line(size = 1.5) +
      geom_ribbon(alpha = 0.5) +
      geom_hline(yintercept = 0, lty = 2) +
      facet_wrap(~Years) +
      coord_cartesian(ylim = c(-1, 1))
  )
  

# Growth ------------------------------------------------------------------
  growth_sum = growth_df %>% 
    mutate(RE = (Pred - True) / True) %>% 
    group_by(Type, Age, Sex) %>% 
    summarize(median = median(RE),
              lwr_95 = quantile(RE, 0.025),
              upr_95 = quantile(RE, 0.975))  
  
  print(
    ggplot(growth_sum, aes(x = Age, y = median,
                                 ymin = lwr_95, ymax = upr_95)) +
      geom_line(size = 1.5) +
      geom_ribbon(alpha = 0.5) +
      geom_hline(yintercept = 0, lty = 2) +
      facet_wrap(~Sex) +
      coord_cartesian(ylim = c(-0.5, 0.5))
  )
  
  print(
    ggplot() +
      geom_line(growth_df, mapping = aes(x = Age, 
                                         y = Pred, group = sim), alpha = 0.5) +
      geom_line(growth_df, mapping = aes(x = Age, 
                                         y = True), color = "red",lty = 2, size = 1.3) +
      facet_wrap(~Sex) 
  )
  

# Selectivity -------------------------------------------------------------
  selex_sum = selex_df %>% 
    mutate(RE = (Pred - True) / True) %>% 
    group_by(Type, Age, Sex) %>% 
    summarize(median = median(RE),
              lwr_95 = quantile(RE, 0.025),
              upr_95 = quantile(RE, 0.975))  
  
  print(
    ggplot(selex_sum, aes(x = Age, y = median,
                           ymin = lwr_95, ymax = upr_95)) +
      geom_line(size = 1.3) +
      geom_ribbon(alpha = 0.5) +
      geom_hline(yintercept = 0, lty = 2) +
      facet_grid(Type~Sex) +
      coord_cartesian(ylim = c(-0.5, 0.5))
  )
  
  print(
    ggplot() +
      geom_line(selex_df, mapping = aes(x = Age, 
                                         y = Pred, group = sim), alpha = 0.5) +
      geom_line(selex_df, mapping = aes(x = Age, 
                                         y = True), color = "red",lty = 2, size = 1.5) +
      facet_grid(Type~Sex) 
  )

  
# Parameters --------------------------------------------------------------
  # a = params %>%
  #   filter(Convergence == "Converged") %>%
  #   mutate(RE = (Pred - Truth) / Truth) %>%
  #   filter(Type == "R0")
  # median(a$RE)
  
  print(
    params %>% 
      mutate(RE = (Pred - Truth) / Truth) %>% 
      ggplot(aes(x = Type, y = RE)) +
      geom_violin() +
      geom_boxplot(width = 0.1, outlier.colour = NA) +
      facet_wrap(~Type, scales = "free") +
      geom_hline(yintercept = 0, lty = 2)
  )
  
  dev.off()
}