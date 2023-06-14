
# Packages ---------------------------------------------------------------
library(shiny)
library(tidyverse)
library(latex2exp)
library(future)
library(truncdist)
library(shinyscreenshot)
library(shinyjs)
library(patchwork)


# Read Data -----------
state_name <- "ma"
covid_county <- "./example_county_data/ma_county_biweekly.RDS" %>%
  readRDS()  %>%
  select(-date) %>%
  distinct() %>%
  filter(fips == "25025")

# covid_county <- covid_county[1:15,]

corrected_original <- readRDS(
               "./example_county_data/example_county.RDS")


dates <- readRDS( "./example_county_data/date_to_biweek.RDS")



background <- "<b></i><span style='font-size:120%;' >Background:</b></i> </span> In the implementation of Bayesian melding,
we sample values of $\\theta = \\{ \\alpha, \\beta, P(S_1|\\text{untested})\\}$, where the sample size here is selected via the 
<i>Sample Size</i> option.
Then, we use the function $M: \\theta \\to \\phi$ defined by
$M(\\theta) = \\dfrac{\\beta(1- P(S_1|\\text{untested}))}{\\beta(1-P(S_1|\\text{untested})) + \\alpha P(S_1|\\text{untested})}$,
yielding the induced distribution $f_\\phi^{induced}$. Because we do not have an analytical form for $f_\\phi^{induced}$,
we estimate it using kernel density estimation.

Then, we compute the sampling weights $\\left(\\dfrac{f_\\phi^{direct}(\\phi_i)}{f_{\\phi}^{induced}(\\phi_i)}\\right)^{0.5}$.
Here, $f_\\phi^{direct}(\\phi_i)$ is calculated using the density function $f_\\phi^{direct}$, which is defined based on
information from meta-analyses on the asymptomatic rate, and $f_{\\phi}^{induced}(\\phi_i)$ is from the kernel density estimation.

Sampling with these weights from our sample of $M(\\theta) = \\phi$ and our sample of $\\theta$ yields the post-melding distributions,
where the post-melding density of $\\phi$ is the the logarithmic pool of $f_\\phi^{direct}$ and $f_\\phi^{induced}$, that is,
$f_\\phi^{post} (\\phi_i) = \\left({f_\\phi^{direct}(\\phi_i)}\\right)^{0.5} \\; \\left( f_{\\phi}^{induced}(\\phi_i)\\right)^{0.5}.$
The number of samples we draw from the pooled distribution is defined by the <i>Resample Size</i> option. The resample size should always
be smaller than the sample size, and a higher ratio of sample size to resample size will lead to less irregularities in the estimated melded 
distribution." 



source("./base_functions.R")




#' INDUCED PRIOR ON ASYMPTOMATIC RATE  P(S_0|test+,untested)
#' input sampled values of theta and compute M(\theta)
est_P_A_testpos = function(P_S_untested, alpha, beta){
  beta * (1 - P_S_untested) / (( beta * (1 - P_S_untested)) + (alpha * P_S_untested))
}


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
                       bde = FALSE) {
  
  progress_denominator = ifelse(include_corrected, 3, 1)
  
  given_args <- as.list(environment())
  cat("Arguments to get_melded:\n")
  print(given_args)
  
  
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
  
  incProgress(.2/progress_denominator, detail = paste("Computing induced density of phi..."))
  
  # approximate induced distribution via a density approximation
  if(bde) {
    phi_induced_density <- bde::bde(dataPoints = phi, dataPointsCache = phi, estimator = "betakernel")
    phi_sampled_density <- phi_induced_densityy@densityCache
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
  

  incProgress(.4/progress_denominator, detail = paste("Calculating weights..."))
  
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
  
  corrected <- pmap_df(county_df, ~ {
    incProgress(.66/nrow(county_df),
                detail = paste("Biweek ",
                               unique(list(...)$biweek)))
    process_priors_per_county(
      priors = melded_df,
      county_df = list(...),
      nsamp = post_nsamp) %>%
      generate_corrected_sample(., num_reps = num_reps_counts) %>%
      summarize_corrected_sample() })
  
  corrected %>%
    left_join(dates)
}





#' reformat for plot generation
reformat_melded <- function(melded_df,
                            theta_df,
                            pre_nsamp,
                            p_s0_pos_mean,
                            p_s0_pos_sd,
                            p_s0_pos_bounds,
                            include_corrected = TRUE) {
  
  progress_denominator= ifelse(include_corrected, 3, 1)
  incProgress(.4/progress_denominator, detail = paste("Generating plot..."))
  
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







# Shiny App Functions  ------

# Define UI for application that draws a histogram
ui <- fluidPage(

    tags$head(
      tags$div(HTML("<script type='text/x-mathjax-config' >
            MathJax.Hub.Config({
            tex2jax: {inlineMath: [['$','$'], ['\\(','\\)']]}
            });
            </script >
            ")),
        # Note the wrapping of the string in HTML()
        tags$style(HTML("
        #beta_plot, #alpha_plot, #s_untested_plot, #asymp_plot { overflow:hidden;  }
        #submit_melding, #reset_input {color: #fff;
                        background-color: #337ab7;
                        border-color: #2e6da4;
                        font-weight:bold;
                        size:2vw;}
        #submit_melding:hover, #reset_input:hover {
                        background-color:#78CFDE;
                        border-color:gold;
                        border-width:2px;}
        #screenshot {
                        float:right;
                        color: #fff;
                        background-color: #4A706D;
                        border-color: #2e6da4;
                        font-weight:bold;
                        size:2vw;}
        #screenshot:hover {
                        background-color:#659D99;
                        border-color:gold;
                        border-width: 2px;}
        input[type=checkbox] { transform: scale(1.8);}
        #loadmessage {
           //  position: fixed;
             float:center;
             top: 50%;
             left: 0px;
             width: 100%;
             padding: 10px 0px 10px 0px;
             text-align: center;
             font-weight: bold;
             font-size: 180%;
             color: #000000;
             background-color: #659D99;
             z-index: 105;
             border-radius: 15px;
           }",
        ))),


    # Application title
    titlePanel("Bayesian Melding"),
    withMathJax(),
    useShinyjs(),

    # Sidebar with a slider inputs for bias parameters
    sidebarLayout(
    # USER INPUTS -----------
        sidebarPanel(


          selectInput("pre_nsamp",
                      label =  "$$\\textbf{Sample Size for } \\\\ \\textbf{Bayesian Melding:}$$",
                      choices = list(1e3,
                                     1e4,
                                     1e5,
                                     1e6),
                      selected = 1e5),
          selectInput("post_nsamp",
                      label =  "$$\\textbf{Resample Size for } \\\\ \\textbf{Bayesian Melding:}$$",
                      choices = list(1e3,
                                     1e4,
                                     1e5,
                                     1e6),
                      selected = 1e4),
          h6("$$\\textbf{Include Corrected Estimates} \\\\ \\textbf{Using Melded Distributions:}$$"),
          checkboxInput(label =  "",
                        inputId = "include_corrected", value = FALSE),

            h4("$$\\underline{\\;\\textbf{Distributions for }\\theta\\;\\;}$$"),
            h5("$$\\textbf{Distribution of } \\alpha\\\\(\\textit{Gamma Distribution})$$"),

            sliderInput("alpha_mean",
                        "$$\\text{ Mean for } \\alpha$$",
                        min = 0,
                        max = 5,
                        value = 0.9,
                        step = 0.1),
            sliderInput("alpha_sd",
                        "$$\\text{ Standard deviation for } \\alpha$$",
                        min = 0,
                        max = 2,
                        value = 0.04,
                        step = 0.01),
            p("$$\\text{ Truncation Bounds for } \\alpha$$"),

            checkboxInput(inputId = "trunc_alpha", label = "Truncate", value = TRUE),
            sliderInput(inputId = "alpha_bounds",
                        label = "",
                        min = 0,
                        max = 4,
                        value=c(.8, 1),
                        step = .1),
            br(),

            h5("$$\\textbf{Distribution of } \\beta\\\\(\\textit{Beta Distribution})$$"),
            sliderInput(inputId = "beta_mean",
                        "$$\\text{ Mean for } \\beta$$",
                        min = 0,
                        max = 1,
                        value = 0.15),
            sliderInput("beta_sd",
                        "$$\\text{ Standard deviation for } \\beta$$",
                        min = 0,
                        max = 0.5,
                        value = 0.09),
            p("$$\\text{ Truncation Bounds for } \\beta$$"),

            checkboxInput(inputId = "trunc_beta", label = "Truncate", value = TRUE),
            sliderInput(inputId = "beta_bounds",
                        label = "",
                        min =0,
                        max = 1,
                        value=c( 0.002, 0.4),
                        step = 0.001),
            br(),

            h5("$$\\textbf{Distribution of } P(S_1|untested)\\\\(\\textit{Beta Distribution})$$"),
            sliderInput("s_untested_mean",
                        "$$\\text{ Mean for } P(S_1|untested)$$",
                        min = 0,
                        max = 1,
                        value = 0.025,
                        step = .001),
            sliderInput("s_untested_sd",
                        "$$\\text{ Standard deviation for }\\\\ P(S_1|untested)$$",
                        min = 0,
                        max = 0.5,
                        value =  0.0225,
                        step = .0001),
            p("$$\\text{ Truncation Bounds for } \\\\P(S_1|untested)$$"),

            checkboxInput(inputId = "trunc_s_untested", label = "Truncate", value = TRUE),
            sliderInput(inputId = "s_untested_bounds",
                        label = "",
                        min = 0,
                        max = 1,
                        value=c(0, .15)),

            h4("$$\\underline{\\;\\textbf{Distribution for }\\phi\\;\\;}\\\\(\\textit{Beta Distribution})$$"),
            sliderInput("p_s0_pos_mean",
                        "$$\\text{ Mean for }\\\\ P(S_0|untested, test +)$$",
                        min = 0,
                        max = 1,
                        value =  0.4),
            sliderInput("p_s0_pos_sd",
                        "$$\\text{ Standard deviation for }\\\\ P(S_0|untested, test +)$$",
                        min = 0,
                        max = 0.5,
                        value =  0.1225,
                        step = .0001),

            p("$$\\text{ Truncation Bounds for }\\\\ P(S_0|untested, test +)$$"),
            checkboxInput(inputId = "trunc_p_s0_pos",
                          label = "Truncate",
                          value = TRUE),
            sliderInput(inputId = "p_s0_pos_bounds",
                        label = "",
                        min = 0,
                        max = 1,
                        value=c(.25, .7),
                        step = .01),
            actionButton(inputId ="reset_input",
                         label ="Reset inputs"),
            width = 3
        ),

         # PLOT OUTPUTS -----------

        # Show a plot of the generated distribution
        mainPanel(
         # LOADING
            conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                             tags$div("Loading...",id="loadmessage") ),
            br(),

            # run melding upon click of button
            actionButton(
                inputId = "submit_melding",
                label = "Run Melding"
            ),
            actionButton(inputId = "screenshot",
                         label ="Take a screenshot"),
            br(), br(),
            p(HTML(paste0("<b></i><span style='font-size:120%;' >Instructions:</b></i> </span>",
            "Change input parameters and press <i>Submit Melding</i> to generate ",
            "the melded distributions. Melding is only performed once <i>Submit Melding</i> is pressed.",
            "<br>Use <i>Include Correct Estimates</i> to see the effect on the final corrected estimates for Suffolk County in Massachusetts. ",
            "Note that run time will be substantially longer when including corrected estimates."))),
            br(),
            p(HTML((background))),
            h4(paste0("$$\\textbf{Distributions for }\\mathbf{\\boldsymbol{\\theta}",
            "= \\{ \\boldsymbol{\\alpha, \\beta,} P(S_1|untested)\\}}$$")),
            br(),
            fluidRow(
                splitLayout(cellWidths = c("33%", "33%", "33%"),
                            plotOutput("alpha_plot",  height = "200px"),
                            plotOutput("beta_plot",   height = "200px"),
                            plotOutput("s_untested_plot",   height = "200px"))),
            fluidRow(
                splitLayout(
                    cellWidths = c("50%"),
                    plotOutput("induced_plot", height = "300px"),
                    plotOutput("asymp_plot", height = "300px"))),
            br(),
            h4("$$\\textbf{Distributions After Bayesian Melding:}$$"),
            fluidRow(
                plotOutput("melded_plot", height = "700px")),
            width = 9

        )
    )
)

server <- function(input, output, session) {

    # reset button --------
    observeEvent(input$screenshot, {
        time <- gsub(":", ".", Sys.time())
        time <- gsub(" ", "_", time)
        file_name=paste0("shiny_", time)
        screenshot(filename=file_name)
    })

    observeEvent( eventExpr = input$reset_input,
                  priority = 1,
                  handlerExpr = {

        updateNumericInput(session, "alpha_mean", value = 0.9)
        updateNumericInput(session, "alpha_sd", value = 0.04)
        updateNumericInput(session, "alpha_bounds", value=c(.8, 1))
        updateCheckboxInput(session, "trunc_alpha", value = TRUE)


        updateNumericInput(session, "beta_mean", value = 0.15)
        updateNumericInput(session, "beta_sd", value = 0.09)
        updateNumericInput(session, "beta_bounds", value = c(0.002,.4))
        updateCheckboxInput(session, "trunc_beta", value = TRUE)


        updateNumericInput(session, "s_untested_mean", value = .025)
        updateNumericInput(session, "s_untested_sd", value = 0.0225)
        updateNumericInput(session, "s_untested_bounds", value = c(0,.15))
        updateCheckboxInput(session, "trunc_s_untested", value = TRUE)

        updateNumericInput(session, "p_s0_pos_mean", value = 0.4)
        updateNumericInput(session, "p_s0_pos_sd", value = 0.1225)
        updateNumericInput(session, "p_s0_pos_bounds", value = c(0.25,.7))
        updateCheckboxInput(session, "trunc_p_s0_pos", value = TRUE)

        output$melded_plot <<- NULL

    })

    # alpha plot --------

    output$alpha_plot <- renderPlot({


        # print(input$alpha_mean)
        # print(input$alpha_sd)


        inp_bounds <- if(input$trunc_alpha) input$alpha_bounds else NA


        # message(input$alpha_bounds)
        # message(input$trunc_alpha)
        # message("alpha bounds", inp_bounds)
        # message(is.na(inp_bounds))

        tibble(Value = sample_gamma_density(1e5,
                                            mean = input$alpha_mean,
                                            sd = input$alpha_sd,
                                            bounds = inp_bounds)) %>%
            ggplot(aes(x = Value)) +
            geom_density(fill = "#4F4B68",
                         color = "#4F4B68",
                         linewidth = .2) +
            theme_bw() +
            theme(plot.title = element_text(hjust = .5, size = 16),
                  plot.subtitle = element_text(hjust = .5, size = 11),
                  axis.title = element_text(size = 16),
                  axis.text.x = element_text(size = 12)) +
            labs(title = latex2exp::TeX("Prior for $\\alpha$"),
                 subtitle = ifelse(length(inp_bounds) ==2,
                                   paste0("Mean = ",
                                   input$alpha_mean,
                                   ", Standard Deviation = ",
                                   input$alpha_sd,
                                   "\nTruncation Bounds: ",
                                   "(", inp_bounds[1], ", ", inp_bounds[2], ")"),
                                   paste0("Mean = ",
                                          input$alpha_mean,
                                          ", Standard Deviation = ",
                                          input$alpha_sd)),

                 y = "Density",
                 x = "Value")  +
            # add point just to ensure axis starts at 0
            geom_point(aes(x=0, y = 0), size=0) +
            scale_x_continuous(n.breaks = 6,
                               limits = c(0,
                                          input$alpha_mean + 6*input$alpha_sd))


    })

    # beta plot --------
    output$beta_plot <- renderPlot({

        inp_bounds <- if(input$trunc_beta) input$beta_bounds else NA


        ggplot() +
            stat_function(fun = beta_density,
                          geom="area",
                          args = list("mean" = input$beta_mean,
                                      "sd" = input$beta_sd,
                                      "bounds" = inp_bounds),
                          fill = "#4F4B68") +
            theme_bw() +
            theme(plot.title = element_text(hjust = .5, size = 16),
                  plot.subtitle = element_text(hjust = .5, size = 11),
                  axis.title = element_text(size = 16),
                  axis.text.x = element_text(size = 12)) +
            labs(title = latex2exp::TeX("Prior for $\\beta$"),
                 # subtitle = paste0("Mean = ",
                 #                   input$beta_mean,
                 #                   ", Standard Deviation = ",
                 #                   input$beta_sd ),
                 subtitle = ifelse(length(inp_bounds) ==2,
                                   paste0("Mean = ",
                                          input$beta_mean,
                                          ", Standard Deviation = ",
                                          input$beta_sd,
                                          "\nTruncation Bounds: ",
                                          "(", inp_bounds[1], ", ", inp_bounds[2], ")"),
                                   paste0("Mean = ",
                                          input$beta_mean,
                                          ", Standard Deviation = ",
                                          input$beta_sd)),
                 y = "Density",
                 x = "Value") +
            xlim(0,1)
    })


    # s_untested plot --------
    output$s_untested_plot <- renderPlot({
        inp_bounds <- if(input$trunc_s_untested) input$s_untested_bounds else NA

        ggplot() +
            stat_function(fun = beta_density,
                          geom="area",
                          args = list("mean" = input$s_untested_mean,
                                      "sd" = input$s_untested_sd,
                                      "bounds" = inp_bounds),
                          fill = "#4F4B68") +
            theme_bw() +
            theme(plot.title = element_text(hjust = .5, size = 16),
                  plot.subtitle = element_text(hjust = .5, size = 11),
                  axis.title = element_text(size = 16),
                  axis.text.x = element_text(size = 12)) +
            labs(title = latex2exp::TeX("Prior for $\\P(S_1|untested)$"),
                 # subtitle = paste0("Mean = ",
                 #                   input$s_untested_mean,
                 #                   ", Standard Deviation = ",
                 #                   input$s_untested_sd ),
                subtitle = ifelse(length(inp_bounds) ==2,
                        paste0("Mean = ",
                               input$s_untested_mean,
                               ", Standard Deviation = ",
                               input$s_untested_sd,
                               "\nTruncation Bounds: ",
                               "(", inp_bounds[1], ", ", inp_bounds[2], ")"),
                        paste0("Mean = ",
                               input$s_untested_mean,
                               ", Standard Deviation = ",
                               input$s_untested_sd)),
                 y = "Density",
                 x = "Value") +
            xlim(0,1)
    })

    # induced plot --------
    output$induced_plot <- renderPlot({

        nsamp <- 1e5
        theta <- tibble(alpha = sample_gamma_density(nsamp,
                                                    mean = input$alpha_mean,
                                                    sd = input$alpha_sd,
                                                    bounds = input$alpha_bounds),
                        beta = sample_beta_density(nsamp,
                                                   mean = input$beta_mean,
                                                   sd = input$beta_sd,
                                                   bounds = input$beta_bounds),
                        s_untested = sample_beta_density(nsamp,
                                                         mean = input$s_untested_mean,
                                                         sd = input$s_untested_sd,
                                                         bounds = input$s_untested_bounds),
                        induced = est_P_A_testpos(P_S_untested = s_untested,
                                                  alpha = alpha,
                                                  beta=beta))
        theta %>%
            ggplot(aes(x = induced)) +
            geom_density(fill ="#DEB578", color = "#DEB578", alpha = .7) +
            theme_bw() +
            theme(plot.title = element_text(hjust = .5, size = 16),
                  plot.subtitle = element_text(size = 14, hjust = .5),
                  axis.title = element_text(size = 16),
                  axis.text.x = element_text(size = 12)) +
            labs( title = latex2exp::TeX("Induced Prior for $P(S_0|untested,test+)$"),
                  subtitle = latex2exp::TeX("Density of $M(\\theta)$"),
                  y = "Density",
                  x = "Value")
    })


    output$asymp_plot <- renderPlot({

     # message("checked?",input$trunc_p_s0_pos)

      inp_bounds <- if(input$trunc_p_s0_pos) input$p_s0_pos_bounds else NA


        ggplot() +
            stat_function(fun = beta_density,
                          geom="area",
                          args = list("mean" = input$p_s0_pos_mean,
                                      "sd" = input$p_s0_pos_sd,
                                      "bounds" = inp_bounds),
                          fill = "#166C89",
                          alpha = .8) +
            theme_bw() +
            theme(plot.title = element_text(hjust = .5, size = 16),
                  axis.title = element_text(size = 16),
                  plot.subtitle = element_text(size = 11, hjust = .5),
                  axis.text.x = element_text(size = 12)) +
            labs(
                title = latex2exp::TeX("Prior for $\\P(S_0|untested, test +)$"),
                # subtitle = paste0("Based on Meta-analyses on Asymptomatic Rate",
                #                   "\n", "Mean = ",  input$p_s0_pos_mean,
                #                   ", Standard Deviation = ", input$p_s0_pos_sd),
               subtitle = ifelse(length(inp_bounds) ==2,
                       paste0("Based on Meta-analyses on Asymptomatic Rate\n",
                              "Mean = ",
                              input$p_s0_pos_mean,
                              ", Standard Deviation = ",
                              input$p_s0_pos_sd,
                              "\nTruncation Bounds: ",
                              "(", inp_bounds[1], ", ", inp_bounds[2], ")"),
                       paste0("Mean = ",
                              input$p_s0_pos_mean,
                              ", Standard Deviation = ",
                              input$p_s0_pos_sd)),
                y = "Density",
                x = "Value")
    })


    observeEvent(

      eventExpr = input$submit_melding,
      handlerExpr = {
      message("reset input: ", input$reset_input)
      message("submit melding: ", input$submit_melding)

      # melded plot -------
      output$melded_plot <- renderPlot({

        isolate( {

          alpha_bounds <- if(input$trunc_alpha) input$alpha_bounds else NA
          beta_bounds <- if(input$trunc_beta) input$beta_bounds else NA
          p_s0_pos_bounds <- if(input$trunc_p_s0_pos) input$p_s0_pos_bounds else NA
          s_untested_bounds <- if(input$trunc_s_untested) input$s_untested_bounds else NA

          cat("----| Bounds: |----\n",
              "alpha:", alpha_bounds, "\n",
              "beta:", beta_bounds, "\n",
              "P(S1|untested)", s_untested_bounds, "\n",
              "P(S0|untested,+)",p_s0_pos_bounds, "\n",
              "----------------")

            cat(paste0("-----\n",
                       "Sample size (pre):", input$pre_nsamp, "\n",
                       "Sample size (post):", input$post_nsamp, "\n",
                       "Alpha mean: ", input$alpha_mean, "\n",
                       "Alpha sd: ", input$alpha_sd, "\n",
                       "Beta mean: ", input$beta_mean, "\n",
                       "Beta sd: ", input$beta_sd, "\n",
                       "P(S1|untested) mean: ", input$s_untested_mean, "\n",
                       "P(S1|untested) sd: ", input$s_untested_sd, "\n",
                       "P(S0|untested,+) mean: ", input$p_s0_pos_mean, "\n",
                       "P(S0|untested,+) sd: ", input$p_s0_pos_sd, "\n",
                       "-----\n"))



            withProgress(message = "Performing melding computation",
                         value = 0,
                         {

                             melded <- get_melded(alpha_mean =input$alpha_mean,
                                                  alpha_sd = input$alpha_sd,
                                                  alpha_bounds = alpha_bounds,
                                                  beta_mean = input$beta_mean,
                                                  beta_sd = input$beta_sd,
                                                  beta_bounds = beta_bounds,
                                                  s_untested_mean = input$s_untested_mean,
                                                  s_untested_sd = input$s_untested_sd,
                                                  s_untested_bounds = s_untested_bounds,
                                                  p_s0_pos_mean = input$p_s0_pos_mean,
                                                  p_s0_pos_sd = input$p_s0_pos_sd,
                                                  p_s0_pos_bounds = p_s0_pos_bounds,
                                                  pre_nsamp = as.numeric(input$pre_nsamp),
                                                  post_nsamp = as.numeric(input$post_nsamp),
                                                  include_corrected = input$include_corrected
                                                  )

                             melded_long <- reformat_melded(melded_df = melded$post_melding,
                                                            theta_df = melded$pre_melding,
                                                            pre_nsamp = input$pre_nsamp,
                                                            p_s0_pos_mean = input$p_s0_pos_mean,
                                                            p_s0_pos_sd = input$p_s0_pos_sd,
                                                            p_s0_pos_bounds = p_s0_pos_bounds,
                                                            include_corrected = input$include_corrected)


                            p1 <- melded_long %>%
                              filter(name != "$P(S_0|test+,untested)$") %>%
                                 ggplot(aes(x = value, fill = type)) +
                                 geom_density(alpha = .5, show.legend=FALSE) +
                                 facet_wrap(~name,
                                            labeller = as_labeller(
                                                TeX,
                                                default = label_parsed),
                                            ncol = 3,
                                            scales = "free_y") +
                                 theme_bw() +
                                 theme(
                                       axis.title = element_text(size = 18),
                                       axis.text.x = element_text(size = 11),
                                       plot.title =element_text(size = 18,
                                                                margin =margin(0,0, .5,0, 'cm')),
                                       axis.text.y = element_text(size = 11),
                                       strip.text = element_text(size = 16),
                                       legend.text = element_text(size = 16)) +
                                 labs(title = paste0("Sample Size: ", input$pre_nsamp,
                                                    ", Resample Size: ", input$post_nsamp),
                                      fill = "",
                                      y = "Density") +
                                 scale_fill_manual(values = c("#5670BF",
                                                              "#418F6A",
                                                              "#B28542")) +
                                 guides(fill = guide_legend(keyheight = 2,
                                                            keywidth = 2))  +
                                scale_x_continuous(n.breaks = 5,
                                                   limits = c(0,input$alpha_mean + 6*input$alpha_sd))

                          p2 <- melded_long %>%
                            filter(name == "$P(S_0|test+,untested)$") %>%
                            ggplot(aes(x = value, fill = type)) +
                            geom_density(alpha = .5) +
                            facet_wrap(~name,
                                       labeller = as_labeller(
                                         TeX,
                                         default = label_parsed),
                                       ncol = 3,
                                       scales = "fixed") +
                            theme_bw() +
                            theme(
                              axis.title = element_text(size = 18),
                              axis.text.x = element_text(size = 11),
                              axis.text.y = element_text(size = 11),
                              plot.title =element_text(size = 18,
                                                       margin =margin(0,0, .5,0, 'cm')),
                              strip.text = element_text(size = 16),
                              legend.text = element_text(size = 19)) +
                            labs(
                                 fill = "",
                                 y = "Density") +
                            scale_fill_manual(values = c("#5670BF",
                                                         "#418F6A",
                                                         "#B28542")) +
                            guides(fill = guide_legend(keyheight = 3,
                                                       keywidth = 3))  +
                            scale_x_continuous(n.breaks = 6, limits = c(0, 1))


                          if(input$include_corrected) {

                            corrected <- get_corrected_counts(
                              county_df = covid_county,
                              melded_df = melded$post_melding,
                              post_nsamp = input$post_nsamp) %>%
                              mutate(version = "User-specified Priors")


                            p3 <- corrected_original  %>%
                              bind_rows(corrected) %>%
                              ggplot(aes(x = date,
                                         ymin = exp_cases_lb,
                                         ymax= exp_cases_ub,
                                         fill = version)) +
                              geom_ribbon(alpha = .7) +
                              labs(x = "Date",
                                   y = "Estimated Total Cases\nfor 2-week Interval",
                                   fill = "",
                                   title = "Corrected Estimates for Suffolk County") +
                              theme_bw() +
                              theme(axis.text.x = element_text(angle = 30,
                                                               vjust = .5,
                                                               size = 12),
                                    axis.title = element_text(size = 16,
                                                              face = "bold"),
                                    legend.text = element_text(size = 16),
                                    legend.position = "right",
                                    plot.title = element_text(face = "bold",
                                                              hjust = .5,
                                                              size =22)) +
                              scale_y_continuous(labels = scales::comma) +
                              scale_x_date(date_breaks = "2 months", date_labels = "%b %Y") +
                              scale_fill_manual(values = c("#C94136", "#4297F8")) +
                              guides(fill = guide_legend(keyheight = 2, keywidth = 2))




                            cowplot::plot_grid(p1, p2,p3,
                                               nrow=3,
                                               rel_widths = c(4,1,1),
                                               rel_heights = c(3,3,4))
                          }else{
                            cowplot::plot_grid(p1, p2,
                                               nrow=2,
                                               rel_widths = c(4,1))
                          }

                         # p1 / p2 +  plot_layout(nrow =2,widths = c(4,1))
                         }  # withProgress bracket
            ) # close withProgress
          }) # close isolate
       # } # else
    }) # close renderPlot

     # }

    })

#    observeEvent(input$reset_input, {hide("melded_plot")})
 #   observeEvent(input$submit_melded, {message('here');show("melded_plot")})

} # end server function


# Run the application
shinyApp(ui = ui, server = server)



# to update
# library(rsconnect); rsconnect::deployApp(appTitle = "bayesian_melding_priors" )




# RESOURCES ---------
# https://stackoverflow.com/questions/70002107/how-to-build-a-shiny-app-that-shows-the-output-only-when-the-user-clicks-a-click
# strange things happen with rounding
# https://stackoverflow.com/questions/59923693/avoid-sliderinput-rounding
