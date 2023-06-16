
#--------- Melding -------------

get_melded <- function(alpha_mean = 0.9,
                       alpha_sd = 0.04,
                       alpha_bounds = NA,
                       beta_mean = .15,
                       beta_sd =.09,
                       beta_bounds = NA,
                       s_untested_mean = .025,
                       s_untested_sd = .0225,
                       s_untested_bounds = NA,
                       p_s0_pos_mean = .4,
                       p_s0_pos_sd = .1225,
                       p_s0_pos_bounds = NA,
                       pre_nsamp = 1e4,
                       post_nsamp = 1e3,
                       include_corrected = TRUE,
                       bde = FALSE,
                       quiet=TRUE) {
  
  
  given_args <- as.list(environment())
  
  if (!quiet) {
    cat("Arguments to get_melded:\n")
    print(given_args)
  } 
  
  
  theta <- tibble(alpha = sample_gamma_density(pre_nsamp,
                                               mean = alpha_mean,
                                               sd = alpha_sd,
                                               bounds = alpha_bounds),
                  beta= sample_beta_density(pre_nsamp,
                                            mean = beta_mean,
                                            sd = beta_sd,
                                            bounds = beta_bounds),
                  P_S_untested = sample_beta_density(pre_nsamp,
                                                     mean = s_untested_mean,
                                                     sd = s_untested_sd,
                                                     bounds = s_untested_bounds)) %>%
    mutate(phi_induced = est_P_A_testpos(P_S_untested = P_S_untested,
                                         alpha = alpha,
                                         beta=beta))
  
  # theta contains values sampled from alpha, beta, P_S_untested, and M(theta) = phi_induced
  # induced phi
  phi <- theta$phi_induced
  
  # approximate induced distribution via a density approximation
  if(bde) {
    message("Beta Kernel Density Estimation")
    phi_induced_density <- kdensity::kdensity(x = phi,  bw="SJ", adjust = 1, start = "beta", kernel = "beta", support = c(0,1))
    phi_sampled_density <- phi_induced_density(phi) 
  }
  else {
    phi_induced_density <- density(x = phi, n = pre_nsamp, adjust = 2, kernel = "gaussian")
    indexes <- findInterval(phi, phi_induced_density$x)
    phi_sampled_density <- phi_induced_density$y[indexes]
  }
  
  dp_s0_pos <- function(x) {
    
    beta_density(x,
                 mean=p_s0_pos_mean,
                 sd = p_s0_pos_sd,
                 bounds=p_s0_pos_bounds)
  }
  
  
  
  weights <- (phi_sampled_density/ dp_s0_pos(phi))^(.5)
  
  
  post_samp_ind <-sample.int(n=pre_nsamp,
                             size=post_nsamp,
                             prob=1/weights,
                             replace=TRUE)
  
  
  pi_samp <- cbind(theta[post_samp_ind,],
                   P_A_testpos =  phi[post_samp_ind]) %>%
    select(-phi_induced)
  
  
  return(list(post_melding = pi_samp, pre_melding = theta))
  
  
}



get_corrected_counts <- function(county_df, 
                                 melded_df, 
                                 post_nsamp,
                                 num_reps_counts = 1e3) {
  
  # melded_df <- melded_df %>%
  #   rename(Z_S = alpha,
  #          Z_A = beta)
  
  corrected <- pmap_df(county_df, function(posrate,
                                           population,
                                           total, 
                                           positive,
                                           biweek,
                                           fips,
                                           ...) {
    process_priors_per_county(
      priors = melded_df,
      county_df = tibble(posrate, 
                         population,
                         total, 
                         positive,
                         biweek,
                         fips),
      nsamp = post_nsamp) %>%
      generate_corrected_sample(., num_reps = 1e3) %>%
      summarize_corrected_sample() })
  
  corrected %>%
    left_join(dates)
}


#--------- Plotting -------------



#' reformat for plot generation
reformat_melded <- function(melded_df,
                            theta_df,
                            pre_nsamp,
                            p_s0_pos_mean,
                            p_s0_pos_sd,
                            p_s0_pos_bounds) {
  
  melded_df_long <- melded_df %>%
    pivot_longer(cols=everything()) %>%
    mutate(type = "After Melding")
  
  
  melded <- theta_df %>%
    mutate(P_A_testpos = sample_beta_density(pre_nsamp,
                                             mean = p_s0_pos_mean,
                                             sd = p_s0_pos_sd,
                                             bounds = p_s0_pos_bounds)) %>%
    pivot_longer(cols=everything()) %>%
    mutate(type = ifelse(
      name == "phi_induced",
      "Induced", "Before Melding")) %>%
    mutate(name = ifelse(name == "phi_induced",
                         "P_A_testpos",
                         name)) %>%
    bind_rows(melded_df_long) %>%
    mutate(name = case_when(
      name == "alpha" ~"$\\alpha$",
      name == "beta" ~"$\\beta$",
      name == "P_A_testpos" ~ "$P(S_0|test+,untested)$",
      name == "P_S_untested" ~ "$P(S_1|untested)$")
    ) %>%
    mutate(name = factor(name,
                         levels = c(
                           "$\\alpha$",
                           "$\\beta$",
                           "$P(S_1|untested)$",
                           "$P(S_0|test+,untested)$")))
  
}


plot_melded <- function(melded, custom_title="", nsamp) {
  
  
  p1 <- melded %>%
    filter(name != "$P(S_0|test+,untested)$") %>%
    ggplot(aes(x = value, fill = type)) +
    geom_density(alpha = .5, show.legend=FALSE) +
    facet_wrap(~name,
               labeller = as_labeller(
                 TeX,   default = label_parsed),
               ncol = 3,
               scales = "fixed") +
    theme_bw() +
    theme(
      # axis.text.y = element_blank(),
      # axis.ticks.y = element_blank(),
      axis.title = element_text(size = 11),
      axis.text.x = element_text(size = 10),
      plot.title =element_text(size = 11,
                               margin =margin(0,0, .5,0, 'cm')),
      strip.text = element_text(size = 11,
                                color="white"),
      legend.text = element_text(size = 16)) +
    labs(title = TeX(custom_title,bold=TRUE),
         subtitle =paste0("Number of Samples: ", nsamp),
         fill = "",
         y = "Density") +
    scale_fill_manual(values = c("#5670BF", "#418F6A","#B28542")) +
    guides(fill = guide_legend(keyheight = 2,  keywidth = 2))
  
  p2 <- melded %>%
    filter(name == "$P(S_0|test+,untested)$") %>%
    ggplot(aes(x = value, fill = type)) +
    geom_density(alpha = .5) +
    facet_wrap(~name,
               labeller = as_labeller(
                 TeX,   default = label_parsed),
               ncol = 3,
               scales = "fixed") +
    theme_bw() +
    theme(
      # axis.text.y = element_blank(),
      # axis.ticks.y = element_blank(),
      axis.title = element_text(size = 13),
      axis.text.x = element_text(size = 10),
      plot.title =element_text(size = 13,
                               margin =margin(0,0, .5,0, 'cm')),
      strip.text = element_text(size = 11, color="white"),
      legend.text = element_text(size = 13)) +
    labs(
      #title = paste0("Number of Samples: ", nsamp),
      fill = "",
      y = "Density") +
    scale_fill_manual(values = c("#5670BF", "#418F6A","#B28542")) +
    guides(fill = guide_legend(keyheight = 2,  keywidth = 2)) +
    xlim(0,1)
  
  
  p1 / p2 +  plot_layout(nrow =2, widths = c(4,1))
  
}








plot_melded <- function(melded, custom_title="", pre_nsamp, post_nsamp) {
  
  
  p1 <- melded %>%
    filter(name != "$P(S_0|test+,untested)$") %>%
    ggplot(aes(x = value, fill = type)) +
    geom_density(alpha = .5, show.legend=FALSE) +
    facet_wrap(~name,
               labeller = as_labeller(
                 TeX,   default = label_parsed),
               ncol = 3,
               scales = "fixed") +
    theme_bw() +
    theme(
      # axis.text.y = element_blank(),
      # axis.ticks.y = element_blank(),
      axis.title = element_text(size = 13),
      axis.text.x = element_text(size = 10),
      plot.title =element_text(size = 13,
                               margin =margin(0,0, .5,0, 'cm')),
      strip.text = element_text(size = 11,color="white"),
      strip.background = element_rect(fill = "#3E3D3D"),
      legend.text = element_text(size = 16)) +
    labs(title = TeX(custom_title,bold=TRUE),
         subtitle =paste0("Sample Size: ", pre_nsamp, "\nResample Size: ", post_nsamp),
         fill = "",
         y = "Density") +
    scale_fill_manual(values = c("#5670BF", "#418F6A","#B28542")) +
    guides(fill = guide_legend(keyheight = 2,  keywidth = 2))
  
  p2 <- melded %>%
    filter(name == "$P(S_0|test+,untested)$") %>%
    ggplot(aes(x = value, fill = type)) +
    geom_density(alpha = .5) +
    facet_wrap(~name,
               labeller = as_labeller(
                 TeX,   default = label_parsed),
               ncol = 3,
               scales = "fixed") +
    theme_bw() +
    theme(
      strip.text = element_text(size = 11, color = "white"),
      strip.background = element_rect(fill = "#3E3D3D"),
      # axis.text.y = element_blank(),
      # axis.ticks.y = element_blank(),
      axis.title = element_text(size = 13),
      axis.text.x = element_text(size = 10),
      plot.title =element_text(size = 13,
                               margin =margin(0,0, .5,0, 'cm')),
      legend.text = element_text(size = 13)) +
    labs(
      #title = paste0("Number of Samples: ", nsamp),
      fill = "",
      y = "Density") +
    scale_fill_manual(values = c("#5670BF", "#418F6A","#B28542")) +
    guides(fill = guide_legend(keyheight = 2,  keywidth = 2)) +
    xlim(0,1)
  
  
  p1 / p2 +  plot_layout(nrow =2, widths = c(4,1))
  
}
