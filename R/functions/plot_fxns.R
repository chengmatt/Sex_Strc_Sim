# Purpose: Utility functions but for plotting purposes
# Creator: Matthew LH. Cheng (UAF-CFOS)
# Date: 9/11//23

#' Title Plot relative error as a time series
#'
#' @param df Dataframe
#' @param ylab label for y axis 
#'
#' @return
#' @export
#'
#' @examples
plot_re_ts = function(df, ylab) {
  with(df,  plot(1, type = "n", xlab = "Years",
                     ylab = ylab,  ylim = c(-0.8, 0.8), xlim = c(1, max(years)), xaxs="i", cex.lab = 1.75,
                     cex.axis = 2)) # xaxs i gets rid of bordering
  with(df, polygon(c(years, rev(years)), c(lwr95, rev(upr95)), border = NA, col = adjustcolor("#9fcae1", 1))) # 95% quantiles
  with(df, polygon(c(years, rev(years)), c(lwr50, rev(upr50)), border = NA, col = adjustcolor("#1a59a1", 0.55))) # 50% quantiles
  with(df,  lines(years, median, type = "l", lwd = 2.5, col = adjustcolor("#044391", 0.75))) # plot median
  abline(h = 0, lty = 2, lwd = 2) # get horizontal line (no bias)
}

#' Title Plot relative error as a bean plot
#'
#' @param df Dataframe 
#'
#' @return
#' @export
#'
#' @examples
  plot_re_par = function(df) {
  vioplot::vioplot(df, rectCol = "#1a59a1", lineCol = "#1a59a1", 
                   colMed = "white", border = NA, col = alpha("#9fcae1", 0.5), 
                   cex.axis = 2)
  mtext(text = "Relative Error in Parameters", side = 2, padj = -2.75, cex = 1.2)
  abline(h = 0, lty = 2, lwd = 2) # get horizontal line (no bias)
}

