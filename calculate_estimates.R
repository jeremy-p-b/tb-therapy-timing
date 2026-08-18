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

descriptive_outputs <- list()

# Read in data
tb_cohort <- read_csv(glue("{DATA_DIR}/patient2.csv"))

# Inspect code
tb_cohort %>% select(!pid) %>% skim() %>% print()

# 
descriptive_outputs$inclusions <- tribble(
  ~`Inclusion criteria`, ~`No. patients remaining`,
  "None", nrow(tb_cohort),
  "No baseline RX", nrow(tb_cohort %>% filter(!baseline_on_rx)),
  "No death on enrollment day", nrow(tb_cohort %>% filter(!baseline_on_rx) %>% filter(!(futime==1 & event==1)))
) %>% mutate(`No. exclusions` = lag(`No. patients remaining`) - `No. patients remaining` )

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
  mutate(pretty_init_immediately=if_else(init_immediately==1, "Immediate initiation", "Delayed initiation")) %>%
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
  mutate(pretty_init_immediately= if_else(init_immediately==1, "Immediate initiation", "Not immediate initiation")) %>%
  tbl_summary(by=pretty_init_immediately, include=c(female, age, hiv_vl_suppressed, cd4_log, hb_log, creat_log, rifresist, bac_load_grm,   impaired_conscious, hypoxia)) %>% 
  add_difference(everything() ~ "smd") %>% 
  remove_abbreviation() %>%
  modify_column_hide(conf.low) 

# Crude event counts among those not censored
descriptive_outputs$crude_counts <- tb_cohort %>% 
  mutate(pretty_init_immediately= if_else(init_immediately==1, "Immediate initiation", "Not immediate initiation")) %>%
  filter(ltfu == 0) %>%
  group_by(pretty_init_immediately) %>%
  summarise(N = n(), Deaths=sum(event)) %>%
  rename(`Treatment initiation`=pretty_init_immediately) %>%
  mutate(Deaths = glue("{Deaths} ({scales::number(100*Deaths/N, accuracy=0.1)}%)"))

saveRDS(descriptive_outputs, "output/descriptive_outputs.rds")

# Compare immediate versus two to six ---------------------------------------------------------

# set globals 
immediate_vs_two_to_six_outputs <- list() 
TRUNCATE_PERCENTILE <- 0.99

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
immediate_vs_two_to_six_outputs$treatment_hazards <- tb_cohort_long %>% 
  filter(day <= 6 & (is.na(init_day) | day <= init_day) & prior_outcome==0 & censor==0 & prior_censor==0) %>%
  mutate(ptreat=if_else(day==1, ptreat_day1, ptreat_day2to6)) %>%
  group_by(day) %>%
  summarise(observed_ptreat=mean(treatment_initiation), estimated_ptreat=mean(ptreat)) %>%
  mutate(across(c(observed_ptreat, estimated_ptreat), ~scales::number(., accuracy=0.001))) %>%
  rename(`Day`=day, `Observed hazard of treatment`=observed_ptreat, `Estimated hazard of treatment`=estimated_ptreat)

# Bootstrap CIs
clust <- makeCluster(5)
clusterSetRNGStream(clust, 743)
clusterCall(clust, load_globals_and_helpers)
res <- pblapply(1:1000, sample_compute_weighted_risks, cohort_data=tb_cohort, weight_function=calculate_weights_immediate_vs_delayed, formula_day1=formula_model_day1, 
                formula_day2to6=formula_model_day2to6, truncate_percentile=TRUNCATE_PERCENTILE, cl=clust)
stopCluster(clust)
immediate_vs_two_to_six_processed_res <- res %>% summarise_boot_weighted(cohort_data=tb_cohort, weight_function=calculate_weights_immediate_vs_delayed, formula_day1=formula_model_day1,
                                formula_day2to6 = formula_model_day2to6, truncate_percentile=TRUNCATE_PERCENTILE) 


immediate_vs_two_to_six_outputs$results <- immediate_vs_two_to_six_processed_res %>% filter(day==83) %>% pretty_print_res() %>% select(analysis, risk_delayed, risk_immediate, risk_difference, risk_ratio) %>% create_risk_forest_plot()

immediate_vs_two_to_six_outputs$results_day28 <- immediate_vs_two_to_six_processed_res %>% filter(day==27) %>% pretty_print_res() %>% select(analysis, risk_delayed, risk_immediate, risk_difference, risk_ratio) %>% create_risk_forest_plot()


# Plot maximum weights in bootstrap resample
immediate_vs_two_to_six_outputs$max_weight_plot  <- res %>% plot_max_weights()

# Plot weighted KM curve
immediate_vs_two_to_six_outputs$survival_curve <- immediate_vs_two_to_six_processed_res %>% filter(analysis=="Fully weighted") %>% create_km_plot()

# Compare bootstrap SEs with robust SE for risk ratio (as a sense check)
tb_cohort_long %>% calculate_gee_risk_ratio()

saveRDS(immediate_vs_two_to_six_outputs, "output/immediate_vs_two_to_six_outputs.rds")

# Secondary analysis - by bacterial load GRM --------------------------------------------------

secondary_analysis_outputs <- list()

TRUNCATE_PERCENTILE <- 0.99

# Estimate risks
tb_cohort_long <- create_long_data(tb_cohort) %>% filter(bac_load_grm>0)
treat_model_day1 <- fit_model_day1(tb_cohort_long, formula_model_day1)
treat_model_day2to6 <- fit_model_day2to6(tb_cohort_long, formula_model_day2to6)
tb_cohort_long %<>% calculate_weights_immediate_vs_delayed(treat_model_day1, treat_model_day2to6, treat_window_end=6) %>% 
  truncate_weights(TRUNCATE_PERCENTILE)
tb_cohort_long %>%  filter(pweight > 0) %>%
  calculate_weighted_risks_combined() %>% 
  calculate_risk_comparisons() %>% 
  pretty_print_res()  %>% 
  filter(day == 83)
tb_cohort_long %>% create_weighted_table()

# Compare bootstrap SEs with robust SE for risk ratio (as a sense check)
tb_cohort_long %>% calculate_gee_risk_ratio()

tb_cohort_long <- create_long_data(tb_cohort) %>% filter(bac_load_grm<=0)
treat_model_day1 <- fit_model_day1(tb_cohort_long, formula_model_day1)
treat_model_day2to6 <- fit_model_day2to6(tb_cohort_long, formula_model_day2to6)
tb_cohort_long %<>% calculate_weights_immediate_vs_delayed(treat_model_day1, treat_model_day2to6, treat_window_end=6) %>% 
  truncate_weights(TRUNCATE_PERCENTILE)
tb_cohort_long %>%  filter(pweight > 0) %>%
  calculate_weighted_risks_combined() %>% 
  calculate_risk_comparisons() %>% 
  pretty_print_res()  %>% 
  filter(day == 83) %>% select(analysis, starts_with("risk_ratio"))
tb_cohort_long %>% create_weighted_table()

# Compare bootstrap SEs with robust SE for risk ratio (as a sense check)
tb_cohort_long %>% calculate_gee_risk_ratio()

# Estimate risks and bootstrap CI - GRM > 0
clust <- makeCluster(5)
clusterSetRNGStream(clust, 357)
clusterCall(clust, load_globals_and_helpers)
res <- pblapply(1:1000, sample_compute_weighted_risks, cohort_data=tb_cohort %>% filter(bac_load_grm>0), weight_function=calculate_weights_immediate_vs_delayed, formula_day1=formula_model_day1, 
                formula_day2to6=formula_model_day2to6, truncate_percentile=TRUNCATE_PERCENTILE, cl=clust)
stopCluster(clust)
grmabove0_processed_res <- res %>% summarise_boot_weighted(cohort_data=tb_cohort %>% filter(bac_load_grm > 0), weight_function=calculate_weights_immediate_vs_delayed, formula_day1=formula_model_day1,
                                                                         formula_day2to6 = formula_model_day2to6, truncate_percentile=TRUNCATE_PERCENTILE) 
# Plot maximum weights 
res %>% plot_max_weights()

# Estimate risks and bootstrap CI - GRM <= 0
clust <- makeCluster(5)
clusterSetRNGStream(clust, 364)
clusterCall(clust, load_globals_and_helpers)
res <- pblapply(1:1000, sample_compute_weighted_risks, cohort_data=tb_cohort %>% filter(bac_load_grm<=0), weight_function=calculate_weights_immediate_vs_delayed, formula_day1=formula_model_day1, 
                formula_day2to6=formula_model_day2to6, truncate_percentile=TRUNCATE_PERCENTILE, cl=clust)
stopCluster(clust)
grm0orless_processed_res <- res %>% summarise_boot_weighted(cohort_data=tb_cohort %>% filter(bac_load_grm <= 0), weight_function=calculate_weights_immediate_vs_delayed, formula_day1=formula_model_day1,
                                                           formula_day2to6 = formula_model_day2to6, truncate_percentile=TRUNCATE_PERCENTILE) 
# Plot maximum weights
res %>% plot_max_weights()

# create combined table
secondary_analysis_outputs$grm_estimates <- bind_rows(grmabove0_processed_res %>% mutate(grm_category="Greater than 0"), grm0orless_processed_res %>% mutate(grm_category="Less than or equal to 0")) %>% 
  filter(day==83 & analysis != "Censor weighted") %>% pretty_print_res() %>% select(grm_category, analysis, risk_delayed, risk_immediate, risk_difference, risk_ratio)


secondary_analysis_outputs$grm_survival_curve <- plot_grid(
 grmabove0_processed_res %>% 
  filter(analysis=="Fully weighted") %>% create_km_plot(50)
,
grm0orless_processed_res %>% 
  filter(analysis=="Fully weighted") %>% create_km_plot(50)
, ncol=2, labels=c("A", "B")) 


# Secondary analysis - naive approach ---------------------------------------------------------

TRUNCATE_PERCENTILE <- 0.99

tb_cohort_survivors <- tb_cohort %>% filter(!is.na(init_day) & init_day <= 6)

# Estimate risks
tb_cohort_survivors_long <- create_long_data(tb_cohort_survivors)
treat_model_day1 <- fit_model_day1(tb_cohort_survivors_long, formula_model_day1)
tb_cohort_survivors_long %<>% calculate_weights_immediate_vs_no_immediate(treat_model_day1) %>% truncate_weights(TRUNCATE_PERCENTILE)
tb_cohort_survivors_long %>%  filter(pweight > 0) %>%
  calculate_weighted_risks_combined() %>% 
  calculate_risk_comparisons() %>% 
  pretty_print_res()  %>% 
  filter(day == 83)
tb_cohort_survivors_long %>% create_weighted_table()
tb_cohort_survivors_long %>% calculate_gee_risk_ratio()

clust <- makeCluster(5)
clusterSetRNGStream(clust, 242)
clusterCall(clust, load_globals_and_helpers)
res <- pblapply(1:1000, sample_compute_weighted_risks, cohort_data=tb_cohort_survivors, weight_function=calculate_weights_immediate_vs_no_immediate, formula_day1=formula_model_day1, 
                truncate_percentile=TRUNCATE_PERCENTILE, cl=clust)
stopCluster(clust)
survivors_processed_res <- res %>% summarise_boot_weighted(cohort_data=tb_cohort_survivors, weight_function=calculate_weights_immediate_vs_no_immediate, formula_day1=formula_model_day1,
                                                           truncate_percentile=TRUNCATE_PERCENTILE) 
secondary_analysis_outputs$survivors_estimates <- survivors_processed_res %>% filter(day==83 & analysis != "Censor weighted") %>% pretty_print_res() %>% select(analysis, risk_delayed, risk_immediate, risk_difference, risk_ratio) %>% create_risk_forest_plot()

# Plot weighted KM curve
secondary_analysis_outputs$survivors_survival_curve <- survivors_processed_res %>% 
  filter(analysis=="Fully weighted") %>% create_km_plot()

# Secondary analysis - initiate immediately versus not---------------------------------------------------------

TRUNCATE_PERCENTILE <- 0.99

clust <- makeCluster(5)
clusterSetRNGStream(clust, 927)
clusterCall(clust, load_globals_and_helpers)
res <- pblapply(1:1000, sample_compute_weighted_risks, cohort_data=tb_cohort, weight_function=calculate_weights_immediate_vs_no_immediate, formula_day1=formula_model_day1, 
                truncate_percentile=TRUNCATE_PERCENTILE, cl=clust)
stopCluster(clust)
notimmediate_processed_res <- res %>% summarise_boot_weighted(cohort_data=tb_cohort, weight_function=calculate_weights_immediate_vs_no_immediate, formula_day1=formula_model_day1,
                                                              truncate_percentile=TRUNCATE_PERCENTILE) 
secondary_analysis_outputs$not_immediate_estimates <- notimmediate_processed_res %>% filter(day==83 & analysis != "Censor weighted") %>% pretty_print_res() %>% select(analysis, risk_delayed, risk_immediate, risk_difference, risk_ratio) %>% create_risk_forest_plot()


saveRDS(secondary_analysis_outputs, "output/secondary_analysis_outputs.rds")


# Sensitivity analysis - With truncation -----------------------------------------------------------------------------

sensitivity_analysis_outputs <- list()

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

sensitivity_analysis_outputs$truncation_results <- bind_rows(immediate_vs_two_to_six_processed_res %>% mutate(analysis_type = "Immediate vs. 2-6 day delay - no truncation"), 
          immediate_vs_two_to_six_processed_res_truncated[[glue(0.99)]] %>% mutate(analysis_type="Immediate vs. 2-6 day delay - weights truncated at 99th percentile"),
          immediate_vs_two_to_six_processed_res_truncated[[glue(0.95)]] %>% mutate(analysis_type="Immediate vs. 2-6 day delay - weights truncated at 95th percentile"),
          immediate_vs_two_to_six_processed_res_truncated[[glue(0.90)]] %>% mutate(analysis_type="Immediate vs. 2-6 day delay - weights truncated at 90th percentile")) %>% 
  pretty_print_res() %>% filter(day == 83 & (analysis == "Fully weighted" | (analysis=="Unweighted" & str_detect(analysis_type, "no truncation")))) %>% 
  select(analysis_type, analysis, risk_delayed, risk_immediate, starts_with("risk_ratio"))

# Sensitivity analysis - without resistance  ---------------------------------------------------------

TRUNCATE_PERCENTILE <- 0.99

tb_cohort_without_resistance <- tb_cohort %>% filter(rifresist==0)

clust <- makeCluster(5)
clusterSetRNGStream(clust, 389)
clusterCall(clust, load_globals_and_helpers)
res <- pblapply(1:1000, sample_compute_weighted_risks, cohort_data=tb_cohort_without_resistance, weight_function=calculate_weights_immediate_vs_delayed, formula_day1=formula_model_day1_norifresist, 
                formula_day2to6=formula_model_day2to6_norifresist, truncate_percentile=TRUNCATE_PERCENTILE, cl=clust)
stopCluster(clust)
without_resistance_processed_res <- res %>% summarise_boot_weighted(cohort_data=tb_cohort_without_resistance, weight_function=calculate_weights_immediate_vs_delayed, formula_day1=formula_model_day1_norifresist,
                                                 formula_day2to6 = formula_model_day2to6_norifresist, truncate_percentile=TRUNCATE_PERCENTILE) 
sensitivity_analysis_outputs$no_resistance_results <- without_resistance_processed_res %>% filter(day==83 & analysis != "Censor weighted") %>% pretty_print_res() %>% select(analysis, risk_delayed, risk_immediate, risk_difference, risk_ratio) %>% create_risk_forest_plot()


# Sensitivity analysis - uniform  ---------------------------------------------------------

TRUNCATE_PERCENTILE <- 0.99

clust <- makeCluster(5)
clusterSetRNGStream(clust, 389)
clusterCall(clust, load_globals_and_helpers)
res <- pblapply(1:1000, sample_compute_weighted_risks, cohort_data=tb_cohort, weight_function=calculate_weights_immediate_vs_delayed_uniform, formula_day1=formula_model_day1, 
                formula_day2to6=formula_model_day2to6, truncate_percentile=TRUNCATE_PERCENTILE, cl=clust)
stopCluster(clust)
uniform_processed_res <- res %>% summarise_boot_weighted(cohort_data=tb_cohort, weight_function=calculate_weights_immediate_vs_delayed_uniform, formula_day1=formula_model_day1,
                                                                    formula_day2to6 = formula_model_day2to6_norifresist, truncate_percentile=TRUNCATE_PERCENTILE) 
sensitivity_analysis_outputs$uniform_results <- uniform_processed_res %>% filter(day==83 & analysis != "Censor weighted") %>% pretty_print_res() %>% select(analysis, risk_delayed, risk_immediate, risk_difference, risk_ratio) %>% create_risk_forest_plot()


saveRDS(sensitivity_analysis_outputs, "output/sensitivity_analysis_outputs.rds")

# Render qmd ----------------------------------------------------------------------------------

quarto::quarto_render("early_vs_delayed.qmd")




          
          
          
