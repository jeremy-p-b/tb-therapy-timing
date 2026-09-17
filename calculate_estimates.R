# Set globals ---------------------------------------------------------------------------------

# Define data dictionary
DATA_DIR <- "/Users/k2587030/Library/CloudStorage/OneDrive-King\'sCollegeLondon/Data/TB"

# Load libraries and helper functions
load_globals_and_helpers <- function() {
  source("scripts/load_globals.R")
  source("R/analysis_helpers.R")
  blas_set_num_threads(1)
}
load_globals_and_helpers()
# Load formulas
source("scripts/formulas.R")

# Load and inspect data -----------------------------------------------------------------------

# Read in data
tb_cohort <- read_csv(glue("{DATA_DIR}/patient2.csv"))

# Inspect code
tb_cohort %>% select(!pid) %>% skim() %>% print()

# Apply exclusions
tb_cohort %<>% filter(!baseline_on_rx)
tb_cohort %<>% filter(!(futime==1 & event==1))

# Clean data ----------------------------------------------------------------------------------

# Fix coding of init_day
tb_cohort %<>% mutate(init_day = init_day -1)

# Add new identifier (regenerated when resampling)
tb_cohort %<>% mutate(newid=row_number())

# Add treatment categorization 
tb_cohort %<>% 
  mutate(init_immediately=as.integer(!is.na(init_day) & init_day==1)) %>%
  mutate(pretty_init_immediately=if_else(init_immediately==1, "Immediate initiation", "Not immediate initiation")) %>%
  mutate(female = as.integer(sex == "F")) %>%
  mutate(hb = exp(hb_log),
         cd4 = exp(cd4_log),
         creat = exp(creat_log))

# Add labels
var_label(tb_cohort) <- list(
  pid = "Patient ID",
  futime = "Follow-up time",
  event = "Outcome",
  baseline_on_rx = "Baseline on TB RX",
  init_day="Day of treatment initiation",
  bac_load_grm = "Bacterial load GRM",
  creat_log = "Log creatinine",
  hb_log = "Log haemoglobin",
  cd4_log = "Log CD4 count",
  impaired_conscious = "Impaired consciousness",
  hypoxia = "Hypoxia",
  age = "Age",
  sex = "Sex",
  hiv_vl_suppressed = "HIV VL supressed",
  rifresist = "Rifampicin resistance",
  female = "Female",
  hb = "Hemoglobin g/dL",
  cd4 = "CD4 count, cells/\u03BCL",
  creat = "Creatinine, \u03BCmol/L"
)

# Preliminary descriptive statistics ----------------------------------------------------------

descriptive_outputs <- list()

# Graph day of initiation
descriptive_outputs$graph_initiation <- tb_cohort %>% 
  filter(!is.na(init_day)) %>%
  ggplot(aes(init_day)) +
  geom_histogram(binwidth = 1) +
  theme_minimal()

# Count how many are lTFU
tb_cohort %<>% mutate(ltfu = as.integer(futime < 84 & event!=1))
descriptive_outputs$number_ltfu <- tb_cohort %>% filter(ltfu==1) %>% nrow()

# Plot survival
descriptive_outputs$survival_plot <- survfit2(Surv(futime, event) ~ 1, data = tb_cohort) %>%
  ggsurvfit() +
  labs(
    x = "Days",
    y = "Overall survival probability"
  )

# Characteristics of those who do and do not initate immediate
descriptive_outputs$baseline_characteristics <- tb_cohort %>%
  tbl_summary(by=pretty_init_immediately, include=c(female, age, hiv_vl_suppressed, cd4, hb, creat, rifresist, bac_load_grm,   impaired_conscious, hypoxia)) %>% 
  add_difference(everything() ~ "smd") %>% 
  remove_abbreviation() %>%
  modify_column_hide(conf.low) 

# Crude event counts among those not censored
descriptive_outputs$crude_counts <- tb_cohort %>% 
  filter(ltfu == 0) %>%
  tabyl(init_immediately, event) %>%
  adorn_percentages() %>% 
  adorn_pct_formatting() %>% 
  adorn_ns(position="front") %>%
  adorn_title()

# Compare immediate versus two to six ---------------------------------------------------------

# set globals 
immediate_vs_two_to_six_outputs <- list() 
TRUNCATE_PERCENTILE <- NULL

# Estimate risks
tb_cohort_long <- create_long_data(tb_cohort) 
treat_model_day1 <- fit_model_day1(tb_cohort_long, formula_model_day1)
treat_model_day2to6 <- fit_model_day2to6(tb_cohort_long, formula_model_day2to6)
tb_cohort_long %<>% calculate_weights_immediate_vs_delayed(treat_model_day1, treat_model_day2to6, treat_window_end=6) %>% truncate_weights(TRUNCATE_PERCENTILE)
tb_cohort_long %>%  filter(pweight > 0) %>%
  calculate_weighted_risks_combined() %>% 
  calculate_risk_comparisons() %>% 
  pretty_print_res()  %>% 
  filter(day == 83)

# Plot weights
immediate_vs_two_to_six_outputs$plot_weights <- tb_cohort_long %>% plot_final_weights() 

# Create weighted table
immediate_vs_two_to_six_outputs$weighted_table <- tb_cohort_long %>% create_weighted_table()

# Compare predicted versus observed treatment hazards to assess model fit
tb_cohort_long %>% 
  filter(day <= 6 & (is.na(init_day) | day <= init_day) & prior_outcome==0 & censor==0 & prior_censor==0) %>%
  mutate(ptreat=if_else(day==1, ptreat_day1, ptreat_day2to6)) %>%
  group_by(day) %>%
  summarise(observed_ptreat=mean(treatment_initiation), estimated_ptreat=mean(ptreat)) %>%
  mutate(across(c(observed_ptreat, estimated_ptreat), ~scales::number(., accuracy=0.001))) %>%
  rename(`Day`=day, `Observed hazard of treatment`=observed_ptreat, `Estimated hazard of treatment`=estimated_ptreat)

# Bootstrap CIs
clust <- makeCluster(5)
clusterSetRNGStream(clust, 123)
clusterCall(clust, load_globals_and_helpers)
res <- pblapply(1:1000, sample_compute_weighted_risks, cohort_data=tb_cohort, weight_function=calculate_weights_immediate_vs_delayed, formula_day1=formula_model_day1, 
                formula_day2to6=formula_model_day2to6, truncate_percentile=TRUNCATE_PERCENTILE, cl=clust)
stopCluster(clust)
immediate_vs_two_to_six_processed_res <- res %>% summarise_boot_weighted(cohort_data=tb_cohort, weight_function=calculate_weights_immediate_vs_delayed, formula_day1=formula_model_day1,
                                formula_day2to6 = formula_model_day2to6, truncate_percentile=TRUNCATE_PERCENTILE) 
immediate_vs_two_to_six_processed_res %>% filter(day==83) %>% pretty_print_res() %>% select(analysis, starts_with("risk_ratio"))

# Plot maximum weights in bootstrap resample
res %>% plot_max_weights()

# Plot weighted KM curve
immediate_vs_two_to_six_processed_res %>% filter(analysis=="Fully weighted") %>% create_km_plot()

# Compare bootstrap SEs with robust SE for risk ratio (as a sense check)
tb_cohort_long %>% calculate_gee_risk_ratio()


# Secondary analysis - by bacterial load GRM --------------------------------------------------

TRUNCATE_PERCENTILE <- 0.99

# Estimate risks
tb_cohort_long <- create_long_data(tb_cohort) %>% filter(bac_load_grm>0)
treat_model_day1 <- fit_model_day1(tb_cohort_long, formula_model_day1)
treat_model_day2to6 <- fit_model_day2to6(tb_cohort_long, formula_model_day2to6)
tb_cohort_long %<>% calculate_weights_immediate_vs_delayed(treat_model_day1, treat_model_day2to6, treat_window_end=6) %>% truncate_weights(TRUNCATE_PERCENTILE)
tb_cohort_long %>%  filter(pweight > 0) %>%
  calculate_weighted_risks_combined() %>% 
  calculate_risk_comparisons() %>% 
  pretty_print_res()  %>% 
  filter(day == 83)
tb_cohort_long %>% create_weighted_table()

# Estimate risks
clust <- makeCluster(5)
clusterSetRNGStream(clust, 123)
clusterCall(clust, load_globals_and_helpers)
res <- pblapply(1:1000, sample_compute_weighted_risks, cohort_data=tb_cohort %>% filter(bac_load_grm>0), weight_function=calculate_weights_immediate_vs_delayed, formula_day1=formula_model_day1, 
                formula_day2to6=formula_model_day2to6, truncate_percentile=TRUNCATE_PERCENTILE, cl=clust)
stopCluster(clust)
grmabove0_processed_res <- res %>% summarise_boot_weighted(cohort_data=tb_cohort %>% filter(bac_load_grm > 0), weight_function=calculate_weights_immediate_vs_delayed, formula_day1=formula_model_day1,
                                                                         formula_day2to6 = formula_model_day2to6, truncate_percentile=TRUNCATE_PERCENTILE) 
grmabove0_processed_res %>% filter(day==83) %>% pretty_print_res() %>% select(analysis, starts_with("risk"))
res %>% plot_max_weights()

clust <- makeCluster(5)
clusterSetRNGStream(clust, 123)
clusterCall(clust, load_globals_and_helpers)
res <- pblapply(1:1000, sample_compute_weighted_risks, cohort_data=tb_cohort %>% filter(bac_load_grm<=0), weight_function=calculate_weights_immediate_vs_delayed, formula_day1=formula_model_day1, 
                formula_day2to6=formula_model_day2to6, truncate_percentile=TRUNCATE_PERCENTILE, cl=clust)
stopCluster(clust)
grm0orless_processed_res <- res %>% summarise_boot_weighted(cohort_data=tb_cohort %>% filter(bac_load_grm <= 0), weight_function=calculate_weights_immediate_vs_delayed, formula_day1=formula_model_day1,
                                                           formula_day2to6 = formula_model_day2to6, truncate_percentile=TRUNCATE_PERCENTILE) 
grm0orless_processed_res %>% filter(day==83) %>% pretty_print_res() %>% select(analysis, starts_with("risk_ratio"))
res %>% plot_max_weights()

# Plot weighted KM curve
grmabove0_processed_res %>% filter(analysis=="Fully weighted") %>% create_km_plot(50)



# Sensitivity analysis - without resistance  ---------------------------------------------------------

TRUNCATE_PERCENTILE <- NULL

without_resistance <- list()

clust <- makeCluster(5)
clusterSetRNGStream(clust, 389)
clusterCall(clust, load_globals_and_helpers)
res <- pblapply(1:1000, sample_compute_weighted_risks, cohort_data=tb_cohort %>% filter(rifresist==0), weight_function=calculate_weights_immediate_vs_delayed, formula_day1=formula_model_day1, 
                formula_day2to6=formula_model_day2to6, truncate_percentile=TRUNCATE_PERCENTILE, cl=clust)
stopCluster(clust)
without_resistance_processed_res <- res %>% summarise_boot_weighted(cohort_data=tb_cohort %>% filter(rifresist==0), weight_function=calculate_weights_immediate_vs_delayed, formula_day1=formula_model_day1,
                                                 formula_day2to6 = formula_model_day2to6, truncate_percentile=TRUNCATE_PERCENTILE) 
without_resistance_processed_res %>% filter(day==83) %>% pretty_print_res() %>% select(analysis, starts_with("risk_ratio"))


# Sensitivity analysis - With truncation -----------------------------------------------------------------------------

immediate_vs_two_to_six_processed_res_truncated <- list()

for (TRUNCATE_PERCENTILE in c(0.99, 0.95, 0.90)) {
  clust <- makeCluster(5)
  clusterSetRNGStream(clust, 123)
  clusterCall(clust, load_globals_and_helpers)
  res <- pblapply(1:1000, sample_compute_weighted_risks, cohort_data=tb_cohort, weight_function=calculate_weights_immediate_vs_delayed, formula_day1=formula_model_day1, 
                  formula_day2to6=formula_model_day2to6, truncate_percentile=TRUNCATE_PERCENTILE, cl=clust)
  stopCluster(clust)
  immediate_vs_two_to_six_processed_res_truncated[[glue(TRUNCATE_PERCENTILE)]] <- res %>% summarise_boot_weighted(cohort_data=tb_cohort, weight_function=calculate_weights_immediate_vs_delayed, formula_day1=formula_model_day1,
                                                                           formula_day2to6 = formula_model_day2to6, truncate_percentile=TRUNCATE_PERCENTILE) 
}

bind_rows(immediate_vs_two_to_six_processed_res %>% mutate(analysis_type = "Immediate vs. 2-6 day delay - no truncation"), 
          immediate_vs_two_to_six_processed_res_truncated[[glue(0.99)]] %>% mutate(analysis_type="Immediate vs. 2-6 day delay - weights truncated at 99th percentile"),
          immediate_vs_two_to_six_processed_res_truncated[[glue(0.95)]] %>% mutate(analysis_type="Immediate vs. 2-6 day delay - weights truncated at 95th percentile"),
          immediate_vs_two_to_six_processed_res_truncated[[glue(0.90)]] %>% mutate(analysis_type="Immediate vs. 2-6 day delay - weights truncated at 90th percentile")) %>% 
  pretty_print_res() %>% filter(day == 83 & analysis == "Fully weighted") %>% 
  select(analysis_type, analysis, risk_immediate, risk_notimmediate, starts_with("risk_ratio")) %>%
  flextable()




# Plot maximum weights in bootstrap resample
res %>% plot_max_weights()

bind_rows(immediate_vs_two_to_six_processed_res %>% 
            mutate(analysis_type = "Immediate vs. 2-6 day delay"), 
          without_resistance_processed_res %>% mutate(analysis_type="Immediate vs. 2-6 day delay - no rifampicin resistance")) %>% 
            pretty_print_res() %>% filter(day == 83 & analysis %in% c("Unweighted", "Fully weighted")) %>% 
            select(analysis_type, analysis, risk_immediate, risk_notimmediate, starts_with("risk_ratio")) %>%
  rename(`Analysis`=analysis_type, `Weighting`=analysis, `Risk immediate`=risk_immediate, `Risk delayed`=risk_notimmediate,
         `Risk ratio (RR)`=risk_ratio, `SE log(RR)`=risk_ratio_se_log) %>%
  flextable()


          
          
          
