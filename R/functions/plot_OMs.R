#' Title Plot operating models
#'
#' @param oms OM list
#' @param path Path to OM folder
#'
#' @return
#' @export
#'
#' @examples
plot_OMs = function(oms, path) {

  pdf(here(path, "OM_Plots.pdf"))
  # Plot NAA ----------------------------------------------------------------
  NAA = reshape2::melt(oms$NAA)
  print(
    ggplot(NAA %>% filter(Var4 == 1), aes(x = Var1, y = value, fill = factor(Var2))) +
      geom_col() +
      facet_wrap(~Var3) +
      labs(x = "Years", y = "Numbers", fill = "Ages")
  )
  
  # Plot CAA ----------------------------------------------------------------
  CAA = reshape2::melt(oms$CAA)
  print(
    ggplot(CAA %>% filter(Var5 == 1), aes(x = Var1, y = value, fill = factor(Var2))) +
      geom_col() +
      facet_wrap(~Var3) +
      labs(x = "Years", y = "Numbers", fill = "Ages")
  )
  
  # Plot Catch By Sex -------------------------------------------------------
  Catch_Sex = reshape2::melt(oms$Total_Catch_Sex)
  print(
    ggplot(Catch_Sex %>% filter(value != 0, Var3 == 1), aes(x = Var1, y = value, 
                                                 color = factor(Var2), group = Var4)) +
      geom_line(size = 1.8, alpha = 0.5) +
      facet_wrap(~Var2) +
      labs(x = "Year", y = "Catch", color = "Sex")
  )
  
  # Plot Fmort --------------------------------------------------------------
  fmort = reshape2::melt(oms$Fmort)
  print(
    ggplot(fmort, aes(x = Var1, y = value)) +
      geom_line() +
      labs(x = "Year", y = "F")
  )
  
  # Plot SSB ----------------------------------------------------------------
  SSB = reshape2::melt(oms$SSB)
  print(
    ggplot(SSB, aes(x = Var1, y = value, group = factor(Var2))) +
      geom_line() +
      labs(x = "Year", y = "SSB")
  )
  
  # Plot Total Biomass ------------------------------------------------------
  Total_Biom = reshape2::melt(oms$Total_Biom)
  print(
    ggplot(Total_Biom, aes(x = Var1, y = value, group = factor(Var2))) +
      geom_line() +
      labs(x = "Year", y = "Total Biomass")
  )
  
  # Plot WAA and LAA --------------------------------------------------------
  laa = reshape2::melt(oms$vonB)
  print(
    ggplot(laa, aes(x = Var1, y = value, color = factor(Var2))) +
      geom_line(size = 1.8) +
      labs(x = "Age", y = "Length", color = "Sex")
  )
  
  waa = reshape2::melt(oms$waa)
  print(
    ggplot(waa, aes(x = Var1, y = value, color = factor(Var2))) +
      geom_line(size = 1.8) +
      labs(x = "Age", y = "Weight", color = "Sex")
  )
  
  print(
    ggplot(oms$Srv_LAA %>% filter(sim == 1), aes(x = ages, y = lens, color = factor(sex))) +
      geom_point() +
      facet_wrap(~sex, scales = "free_x")
  )
  
  # Natural Mortality -------------------------------------------------------
  M = data.frame(M = rep(oms$M, dim(oms$NAA)[1]), year = 1:dim(oms$NAA)[1],
                 sex = c(1, 2))
  print(
    ggplot(M, aes(x = year, y = M, color = factor(sex))) +
      geom_line(size = 1.3) +
      labs(x = "Year", y = "Natural Mortality", color = "Sex")
  )
  
  # Sex Ratio -------------------------------------------------------
  sr = data.frame(sr = rep(oms$sexRatio, dim(oms$NAA)[1]), 
                  year = 1:dim(oms$NAA)[1], sex = c(1, 2))
  print(
    ggplot(sr, aes(x = year, y = sr, color = factor(sex))) +
      geom_line(size = 1.3) +
      labs(x = "Year", y = "Sex Ratio", color = "Sex")
  )
  
  # Plot selex --------------------------------------------------------------
  fish_selex = reshape2::melt(oms$FishAge_Selex)
  print(
    ggplot(fish_selex, aes(x = Var1, y = value, color = factor(Var2))) +
      geom_line(size = 1.8) +
      labs(x = "Age", y = "Fishery Selex", color = "Sex")
  )
  
  srv_selex = reshape2::melt(oms$SrvAge_Selex)
  print(
    ggplot(srv_selex, aes(x = Var1, y = value, color = factor(Var2))) +
      geom_line(size = 1.8) +
      labs(x = "Age", y = "Survey Selex", color = "Sex")
  )
  
  # Plot comps --------------------------------------------------------------
  fish_acomps = reshape2::melt(oms$Fish_AgeComps)
  print(
    ggplot(fish_acomps %>% filter(Var4 == 1, Var5 == 1), 
           aes(x = Var1, y = value, fill = factor(Var2))) +
      geom_col() +
      facet_wrap(~Var3) +
      labs(x = "Year", y = "Numbers - Fishery Comps (Ages)", fill = "Age")
  )
  
  srv_acomps = reshape2::melt(oms$Srv_AgeComps)
  print(
    ggplot(srv_acomps %>% filter(Var4 == 1, Var5 == 1), 
           aes(x = Var1, y = value, fill = factor(Var2))) +
      geom_col() +
      facet_wrap(~Var3) +
      labs(x = "Year", y = "Numbers - Survey Comps (Ages)", fill = "Age")
  )
  
  fish_lcomps = reshape2::melt(oms$Fish_LenComps)
  print(
    ggplot(fish_lcomps %>% filter(Var4 == 1, Var5 == 1), 
           aes(x = Var1, y = value, fill = factor(Var2))) +
      geom_col() +
      facet_wrap(~Var3) +
      labs(x = "Year", y = "Numbers - Fishery Comps (Lengths)", fill = "Lengths")
  )
  
  srv_lcomps = reshape2::melt(oms$Srv_LenComps)
  print(
    ggplot(srv_lcomps %>% filter(Var4 == 1, Var5 == 1), 
           aes(x = Var1, y = value, fill = factor(Var2))) +
      geom_col() +
      facet_wrap(~Var3) +
      labs(x = "Year", y = "Numbers - Survey Comps (Lengths)", fill = "Lengths")
  )
  
  Srv_Index = reshape2::melt(oms$Srv_Index)
  print(
    ggplot(Srv_Index %>% filter(value != 0), aes(x = Var1, y = value, group = factor(Var3))) +
      geom_line() +
      labs(x = "Year", y = "Survey Index")
  )
  
  
  dev.off()
    
}

