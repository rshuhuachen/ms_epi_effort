plot_posterior <- function(posteriors, response, predictor, name, label, minx=-1, maxx=1){
  pacman::p_load(brms, bayesplot, ggridges)
  source("scripts/plotting_theme.R")
  
  # output complete interval
  write.csv(as.data.frame(mcmc_intervals_data(posteriors, prob =0.8, prob_outer = 0.95)), 
            file = paste0("results/modeloutput/pheno/intervals_", name, ".csv"), quote=F, row.names=F)
  
  # get interval of betas
  interval <- mcmc_intervals_data(posteriors, prob =0.8, prob_outer = 0.95) %>% 
    subset(grepl(paste0("b_", response, "_", predictor), parameter))

  interval$response <- response
  
  # get areas
  
  area <- mcmc_areas_data(posteriors) %>%
    subset(grepl(paste0("b_", response, "_", predictor), parameter))
  
  area$response <- response
  
  # split by interval
  brms <- split(area, area$interval)
  
  brms$bottom <- brms$outer %>%
    summarise(
      ll = min(.data$x),
      hh = max(.data$x),
      .groups = "drop_last"
    ) %>%
    ungroup()
  
  clr = clr_high
  
  if (interval$m > 0 & interval$ll > 0 & interval$hh > 0){
    clr = clr_froh}

  if (interval$m < 0 & interval$ll < 0 & interval$hh < 0){
    clr = clr_gerp}

  if(label == "top"){
    label_y = 1.3
  }

   if(label == "bottom"){
    label_y = 0.7
  }

  # plot
  plot <- ggplot(data = brms$outer) +  
    aes(x = .data$x, y = .data$response) + 
    geom_ridgeline(aes(scale = 0.4, height = scaled_density, fill = response, col = response))+
    geom_segment(data=interval, aes(x = l, xend = h, yend = response), col = "black", linewidth=3)+
    geom_segment(data=interval, aes(x = ll, xend = hh, yend = response), col = "black")+
    geom_point(data=interval, aes(x = m, y = response), fill="white",  col = "black", shape=21, size = 8) + 
    geom_vline(xintercept = 0, col = "#ca562c", linetype="longdash", linewidth=0.8)+
    scale_fill_manual(values =alpha(c(clr), 0.7)) +
    xlim(minx, maxx)+
    scale_color_manual(values =c(clr)) +
    annotate("label", x = (minx + maxx)/2, y = label_y, label = paste0("Posterior median = ", format(round(interval$m, digits=2), nsmall=2), ", 
    95% CI = [", format(round(interval$ll, digits=2), nsmall=2), ", ", format(round(interval$hh, digits=2), nsmall=2), "]"),
                                                        size = 10, label.padding = unit(1, "lines"))+
    theme(panel.border = element_blank(),
          panel.grid = element_blank(),
          strip.background = element_blank(),
          legend.position = "none",
          axis.line = element_blank(),
          axis.text.y = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_blank(),
          panel.background = element_rect(fill='transparent'), #transparent panel bg
          plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
          panel.grid.major = element_blank(), #remove major gridlines
          panel.grid.minor = element_blank(), #remove minor gridlines
          legend.background = element_rect(fill='transparent'), #transparent legend bg
          legend.box.background = element_rect(fill='transparent')) #transparent legend panel
  
  ggsave(plot, file = paste0("plots/pheno/posteriors_", name, ".png"), width=10, height=10)
  return(plot)
}


