#Authors: Upendo & Deus
#Date: from 03/05/2024
#Title: SARS-CoV-2 transmission dynamics and implications for non-pharmaceutical interventions before COVID-19 vaccination in Malawi

#====================================================================

#considerations
#compute effective reproduction number, Rt
#compute epidemic doubling time
#relate transmission dynamics to changes of Malawi Government COVID-19 policies

#====================================================================

#load a package "pacman" used for installing and loading other packages
if(!require(pacman)) install.packages("pacman")

#load packages for analysis
pacman::p_load(char = c("rio", "here", "tidyverse", "epicontacts", "EpiNow2", "EpiEstim", "projections", "incidence2", 
                        "scales", "boot", "magrittr", "MASS", "patchwork", "tsibble", "readxl", "epitrix", "distcrete"))

#set seed using a task call for entire session to ensure reproducibility
addTaskCallback(function(...) {set.seed(1988); TRUE})

#turn off the task call to reset seed if needed
#removeTaskCallback(1)

#====================================================================

#load datasets
covid_cases <- rio::import(here("data", "covid_cases.xlsx"))
covid_admit <- rio::import(here("data", "covid_admit.xlsx"))
linelist <- rio::import(here("data", "linelist.rds"))

#====================================================================
#EPINOW2
#====================================================================

#epinow2 requires delays on log scale
#assumed/literature incubation period (time delay between infection and symptom onset)
incub_per1 <- 
  list(
    mean = log(9.1),
    mean_sd = log(0.1),
    sd = log(7.3),
    sd_sd = log(0.1),
    max = 30)

#alternatively estimate incubation period from the dataset if available e.g., linelist
## estimate incubation period
incub_per2 <- 
  EpiNow2::bootstrapped_dist_fit(
    linelist$date_onset - linelist$date_infection,
    dist = "lognormal",
    max_value = 100,
    bootstraps = 1)

#generation time requires data on infection times (infector-infectee pair) and transmission links
#generate contacts
contacts <- 
  linelist %>%
  dplyr::transmute(from = infector, to = case_id) %>%
  tidyr::drop_na()

#generate epicontacts object
epic <- 
  epicontacts::make_epicontacts(
    linelist = linelist,
    contacts = contacts, 
    directed = TRUE)

#estimate gamma generation time
gen_time <- 
  EpiNow2::bootstrapped_dist_fit(
    epicontacts::get_pairwise(epic, "date_infection"),
    dist = "gamma",
    max_value = 20,
    bootstraps = 1)

#calculate daily incidence from onset dates
cases <- 
  linelist %>%
  dplyr::group_by(date = date_onset) %>%
  dplyr::summarise(confirm = n())

#run epinow inference machine
epinow_res <- 
  EpiNow2::epinow(
  reported_cases = cases,
  generation_time = gen_time,
  delays = delay_opts(incub_per2),
  return_output = TRUE,
  verbose = TRUE,
  horizon = 21,
  stan = stan_opts(samples = 750, chains = 4))

#plot summary figure
plot(epinow_res)

#summary table
epinow_res$summary

#extract summary and convert to tibble
estimates <- as_tibble(epinow_res$estimates$summarised)
estimates

#make wide df for median plotting
df_wide <- estimates %>%
  filter(
    variable %in% c("growth_rate", "R"),
    date < as.Date("2014-09-01")) %>%
  
#convert growth rates to doubling times
  mutate(
    across(
      c(median, lower_90:upper_90),
      ~ case_when(
        variable == "growth_rate" ~ log(2)/.x,
        TRUE ~ .x
      )
    ),
#rename variable to reflect transformation
    variable = replace(variable, variable == "growth_rate", "doubling_time")
  )

#make long df for quantile plotting
df_long <- df_wide %>%
  ## here we match matching quantiles (e.g. lower_90 to upper_90)
  pivot_longer(
    lower_90:upper_90,
    names_to = c(".value", "quantile"),
    names_pattern = "(.+)_(.+)"
  )

#make plot
ggplot() +
  geom_ribbon(
    data = df_long,
    aes(x = date, ymin = lower, ymax = upper, alpha = quantile),
    color = NA
  ) +
  geom_line(
    data = df_wide,
    aes(x = date, y = median)
  ) +
#use label_parsed to allow subscript label
  facet_wrap(
    ~ variable,
    ncol = 1,
    scales = "free_y",
    labeller = as_labeller(c(R = "R[t]", doubling_time = "Doubling~time"), label_parsed),
    strip.position = 'left'
  ) +
#manually define quantile transparency
  scale_alpha_manual(
    values = c(`20` = 0.7, `50` = 0.4, `90` = 0.2),
    labels = function(x) paste0(x, "%")
  ) +
  labs(
    x = NULL,
    y = NULL,
    alpha = "Credibel\ninterval"
  ) +
  scale_x_date(
    date_breaks = "1 month",
    date_labels = "%b %d\n%Y"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    strip.background = element_blank(),
    strip.placement = 'outside'
  )


#====================================================================
#EPIESTIM
#====================================================================



