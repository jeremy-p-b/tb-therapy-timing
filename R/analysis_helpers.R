# outcome is yk+1, censor is ck+1, treatment_initiation is ak. Assumed temporal ordering is ck, yk, ak

# Create long version of data
create_long_data <- function(cohort_data, max_day=84) {
  cohort_data <- cohort_data %>% 
    mutate(day=list(1:(max_day-1))) %>% 
    unnest(day) %>%
    mutate(outcome = as.integer(day == futime-1 & event==1)) %>%
    mutate(censor = as.integer(day == futime & event==0)) %>%
    mutate(treatment_initiation = as.integer(!is.na(init_day) & day==init_day)) %>%
    group_by(newid) %>%
    arrange(newid, day) %>%
    mutate(prior_outcome = cummax(lag(outcome, default=0L)),
           prior_censor = cummax(lag(censor, default=0L)),
           prior_treatment = cummax(lag(treatment_initiation, default=0L))) %>%
    ungroup()
  return(cohort_data)
}

# Fit a treatment model to day 1 data
fit_model_day1 <- function(cohort_data_long, chr_formula) {
  day1_data <- cohort_data_long %>% filter(day==1)
  treat_model <- glm(as.formula(chr_formula), data=day1_data, family="binomial")
  return(treat_model)
}

# Fit a treatment model to day 2 to 6 data
fit_model_day2to6 <- function(cohort_data_long, chr_formula) {
  day2to6_data <- cohort_data_long %>% filter(day > 1 & day <= 6 & prior_treatment==0 & prior_censor==0 & prior_outcome==0)
  treat_model <- glm(as.formula(chr_formula),
                     data=day2to6_data, family="binomial")
  return(treat_model)
}

add_nonparametric_censor_weight <- function(cohort_data_long) {
    cohort_data_long %<>% 
      mutate(current_or_prior_treatment = pmax(treatment_initiation, prior_treatment)) %>%
      group_by(day, current_or_prior_treatment*pmax(init_day,0, na.rm = TRUE)) %>%
      mutate(pcensor = weighted_mean(censor, 1-pmax(prior_censor, prior_outcome))) %>%
      mutate(censor_weight = case_when(prior_outcome == 1 ~ 1,
                                       censor == 0 & prior_censor == 0 ~ 1/(1-pcensor), 
                                       .default=0))  %>%
      ungroup()
  return(cohort_data_long)
}

# Calculate weights for treatment strategy of initiate immediately vs. 
calculate_weights_immediate_vs_no_immediate <- function(cohort_data_long, treatment_model) {
  cohort_data_long %<>% mutate(ptreat = predict(treatment_model, newdata=., type="response"))
  cohort_data_long %<>% mutate(week_weight = case_when(init_immediately==1 & day == 1 ~ 1/ptreat, 
                                                        init_immediately==0 & day == 1 ~ 1/(1-ptreat),
                                                        .default=1))
  cohort_data_long %<>% add_nonparametric_censor_weight()
  cohort_data_long %<>% 
    group_by(newid) %>% 
    arrange(newid, day) %>% 
    mutate(censor_weight=cumprod(censor_weight),
           treat_weight=cumprod(week_weight)) %>%
    ungroup() %>%
    mutate(pweight=censor_weight*treat_weight)
  return(cohort_data_long)
}

calculate_weights_immediate_vs_delayed_uniform <- function(cohort_data_long, treat_model_day1, treat_model_day2to6, treat_window_end=6) {
  cohort_data_long %<>% mutate(ptreat_day1 = predict(treat_model_day1, newdata=., type="response"))
  cohort_data_long %<>% mutate(ptreat_day2to6 = predict(treat_model_day2to6, newdata=., type="response"))
  cohort_data_long %<>% 
    mutate(week_weight = case_when(prior_outcome == 1 ~ 1,
                                   censor == 1  | prior_censor == 1 ~ 0,
                                   day == 1 & treatment_initiation==1 ~ 1/ptreat_day1,
                                   day >= 1 & init_immediately==1 ~ 1,     
                                   day == 1 & treatment_initiation == 0 ~ 1/(1-ptreat_day1),
                                   day >= 2 & day <= treat_window_end & treatment_initiation == 1 ~ (1/((treat_window_end+1)-day))/ptreat_day2to6,
                                   day >= 2 & day <= treat_window_end & (is.na(init_day) | init_day > day) ~ (1-(1/(treat_window_end+1-day)))/(1-ptreat_day2to6),
                                   day >= 2 & day <= treat_window_end & !is.na(init_day) & init_day < day ~ 1,
                                   day > treat_window_end ~ 1)) 
  cohort_data_long %<>% add_nonparametric_censor_weight()
  cohort_data_long %<>% 
    group_by(newid) %>% 
    arrange(newid, day) %>% 
    mutate(censor_weight=cumprod(censor_weight),
           treat_weight=cumprod(week_weight)) %>%
    ungroup() %>%
    mutate(pweight=censor_weight*treat_weight)
  return(cohort_data_long)
}

calculate_weights_immediate_vs_delayed <- function(cohort_data_long, treat_model_day1, treat_model_day2to6, treat_window_end=6) {
  cohort_data_long %<>% mutate(ptreat_day1 = predict(treat_model_day1, newdata=., type="response"))
  cohort_data_long %<>% mutate(ptreat_day2to6 = predict(treat_model_day2to6, newdata=., type="response"))
  cohort_data_long %<>% 
    mutate(week_weight = case_when(prior_outcome == 1 ~ 1,
                                   censor == 1  | prior_censor == 1 ~ 0,
                                   day == 1 & treatment_initiation == 1 ~ 1/ptreat_day1,
                                   day == 1 & treatment_initiation == 0 ~ 1/(1-ptreat_day1),
                                   day >= 2 & init_immediately==1 ~ 1,     
                                   day >= 2 & day < treat_window_end ~ 1,
                                   day == treat_window_end & !is.na(init_day) & init_day == treat_window_end ~ 1/ptreat_day2to6,
                                   day == treat_window_end & (is.na(init_day) | init_day > treat_window_end) ~ 0, 
                                   day == treat_window_end & !is.na(init_day) & init_day < treat_window_end ~ 1, 
                                   day > treat_window_end ~ 1)) 
  cohort_data_long %<>% add_nonparametric_censor_weight()
  cohort_data_long %<>% 
    group_by(newid) %>% 
    arrange(newid, day) %>% 
    mutate(censor_weight=cumprod(censor_weight),
           treat_weight=cumprod(week_weight)) %>%
    ungroup() %>%
    mutate(pweight=censor_weight*treat_weight)
  return(cohort_data_long)
}

weighted_mean <- function(x, w) {
  if (is.null(w)) {
    return(mean(x))
  } else {
    return(weighted.mean(x,w))
  }
}

calculate_weighted_risks <- function(cohort_data_long, weight) {
  res <- cohort_data_long %>%
    filter(censor == 0 & prior_censor == 0 & prior_outcome == 0) %>%
    group_by(init_immediately, day) %>% 
    summarise(hazard = weighted_mean(outcome, {{weight}})) %>% 
    mutate(one_minus_hazard = 1 - hazard) %>%
    group_by(init_immediately) %>%
    arrange(day,.by_group=TRUE) %>% 
    mutate(survival = cumprod(one_minus_hazard))  %>%
    mutate(risk = 1 - survival) %>%
    ungroup()
  return(res)
}

calculate_weighted_risks_combined <- function(cohort_data_long) {
  res_censor_weighted <- calculate_weighted_risks(cohort_data_long, censor_weight) %>% mutate(analysis="Censor weighted")
  res_weighted <- calculate_weighted_risks(cohort_data_long, pweight) %>% mutate(analysis="Fully weighted")
  res_unweighted <- calculate_weighted_risks(cohort_data_long, NULL) %>% mutate(analysis="Unweighted")
  res <- bind_rows(res_unweighted, res_censor_weighted, res_weighted)
  return(res)
}

calculate_risk_comparisons <- function(cohort_data_long) {
  cohort_data_long %<>% mutate(exposure=if_else(init_immediately==1, "immediate", "delayed")) %>%
    pivot_wider(id_cols=c(analysis, day), names_from=c(exposure), names_sep="_", values_from=c(survival, risk)) %>%
    mutate(risk_difference = risk_immediate-risk_delayed, risk_ratio=risk_immediate/risk_delayed) 
  return(cohort_data_long)
}

# truncate overall probability weight
truncate_weights <- function(cohort_data_long, percentile) {
  if (!is.null(percentile)) {
    quantile_weight <- quantile(cohort_data_long %>% filter(pweight > 0 & day == 83 & censor==0 & prior_censor==0) %>% {.[["pweight"]]}, percentile)
    cohort_data_long %<>% mutate(pweight=case_when(pweight > quantile_weight[1] ~ quantile_weight[1],
                                                   .default=pweight))
  }
  return(cohort_data_long)
}

compute_weighted_risks <- function(cohort_data, weight_function, formula_day1, formula_day2to6=NULL, restrict_survival=FALSE, truncate_percentile=NULL) {
  cohort_data_long <- create_long_data(cohort_data)
  if (restrict_survival) {
    cohort_data_long %<>% filter(init_day==1 | init_day > 1 & init_day <=6) 
  }
  treatment_model_day1 <- fit_model_day1(cohort_data_long, formula_day1)
  if (!is.null(formula_day2to6)) {
    treatment_model_day2to6 <- fit_model_day2to6(cohort_data_long, formula_day2to6)
    cohort_data_long %<>% weight_function(treatment_model_day1, treatment_model_day2to6)
  } else {
    cohort_data_long %<>% weight_function(treatment_model_day1)
  }
  cohort_data_long %<>% filter(pweight > 0)
  cohort_data_long %<>% truncate_weights(truncate_percentile)
  res <- cohort_data_long  %>% calculate_weighted_risks_combined() %>%
    calculate_risk_comparisons() %>% 
    mutate(min_pweight = min(cohort_data_long$pweight), max_pweight = max(cohort_data_long$pweight))
  return(res)
}

sample_compute_weighted_risks <- function(iter_no, cohort_data, weight_function, formula_day1, formula_day2to6=NULL, restrict_survival=FALSE, truncate_percentile=NULL) {
  sample_cohort <- cohort_data[sample(1:nrow(cohort_data), nrow(cohort_data), replace=TRUE),]
  sample_cohort %<>% mutate(newid=row_number())
  res <- sample_cohort %>% 
    compute_weighted_risks(weight_function=weight_function, formula_day1=formula_day1, formula_day2to6=formula_day2to6, restrict_survival=restrict_survival, truncate_percentile=truncate_percentile) %>%
    mutate(iter = iter_no)
return(res)
}

create_km_plot <- function(res, ylimit=40) {
  day0 <- tibble_row(day=0, risk_delayed=0, risk_delayed_lb=0, risk_delayed_ub=0,
                     risk_immediate=0, risk_immediate_lb=0, risk_immediate_ub=0)
  kmplot <- res %>%
    bind_rows(day0) %>%
    select(day, starts_with("risk_immediate"), starts_with("risk_delayed")) %>%
    rename(risk_immediate_estimate=risk_immediate, risk_delayed_estimate=risk_delayed) %>%
    pivot_longer(starts_with("risk"), names_prefix="risk_", names_to="category", values_to="risk") %>%
    mutate(Group = str_extract(category, "\\w+(?=_)"), quantity=str_extract(category, "(?<=_)\\w+")) %>%
    pivot_wider(id_cols=c(day, Group), names_from=quantity, values_from=risk, names_glue="risk_{quantity}") %>%
    mutate(Group = if_else(Group == "immediate", "Immediate", "Delayed")) %>%
    mutate(across(starts_with("risk"), ~.x*100)) %>%
    rename(`Treatment\nstrategy`=Group) %>%
    ggplot(aes(x=day, y=risk_estimate)) +
    geom_step(aes(color=`Treatment\nstrategy`)) +
    geom_stepribbon(aes(ymin = risk_lb, ymax=risk_ub, fill=`Treatment\nstrategy`), alpha = 0.2) +
    scale_y_continuous(limits=c(0,ylimit), expand = expansion(mult = c(0, 0.05))) +
    scale_x_continuous(limits=c(0,83), expand = expansion(mult = c(0, 0.05))) +
    xlab("Day") +
    ylab("Cumulative risk (%)") +
    theme_classic() +
    theme(
      axis.title = element_text(size=12),
      axis.text = element_text(size=12),
      legend.title = element_text(size=12),
      legend.text = element_text(size=12)
    )
  return(kmplot)
}

# utility function to calculate 95% CI lower bound from vector of bootstrap estimates
lb <- function(vec) {
  return(if_else(all(is.na(vec)), NA_real_, quantile(vec, 0.025, na.rm=TRUE, type=2)))
}

# utility function to calculate 95% CI upper bound from vector of bootstrap estimates
ub <- function(vec) {
  return(if_else(all(is.na(vec)), NA_real_, quantile(vec, 0.975, na.rm=TRUE, type=2)))
}

# Summarise bootstrap estimates
summarise_boot_weighted<- function(res, cohort_data, weight_function, formula_day1, formula_day2to6=NULL, restrict_survival=FALSE, truncate_percentile=NULL) {
  res <- res %>%
    bind_rows() %>%
    group_by(analysis, day) %>%
    summarise(across(c(starts_with("risk_ratio")), list(lb=lb, ub=ub, se_log=~sqrt(var(log(.)))), .names="{.col}_{.fn}"),
              across(c(starts_with("survival"), starts_with("risk_difference"), matches("^risk_immediate$"), matches("^risk_delayed$")),
                     list(lb = lb, ub = ub, se=~sqrt(var(.))), .names = "{.col}_{.fn}"),
              .groups="drop") %>%
    ungroup()
  # get nonbootstrap estimates
  res <- compute_weighted_risks(cohort_data=cohort_data, weight_function=weight_function, formula_day1=formula_day1, formula_day2to6=formula_day2to6, restrict_survival=restrict_survival, truncate_percentile=truncate_percentile) %>%
    left_join(res, by=c("day", "analysis")) %>%
    {select(., order(colnames(.)))} %>%
    relocate(c(day, analysis))
  gc() 
  return(res)
}

pretty_print_res <- function(res) {

  interval_vars <- names(res)[
    paste0(names(res), "_lb") %in% names(res) &
      paste0(names(res), "_ub") %in% names(res)
  ]
  
  if (length(interval_vars) > 0) {
    res <- purrr::reduce(
      interval_vars,
      \(data, variable) {
        digits <- case_when(str_detect(variable, "ratio") ~ 2,
                            TRUE ~ 1)
        multiplier <- case_when(str_detect(variable, "ratio|_se") ~ 1,
                                TRUE ~ 100)
        
        number_format <- paste0(
          "%.", digits, "f (%.", digits, "f, %.", digits, "f)"
        )
        
        data %>%
          mutate(
            "{variable}" := sprintf(
              number_format,
              .data[[variable]] * multiplier,
              .data[[paste0(variable, "_lb")]] * multiplier,
              .data[[paste0(variable, "_ub")]] * multiplier
            )
          )
      },
      .init = res
    ) %>%
      select(-ends_with("_lb"), -ends_with("_ub"))
  } else {
    res %<>% 
      mutate(across(matches("^(risk_immediate|risk_delayed|risk_difference|survival_immediate|survival_delayed)"), ~ scales::number(100*., accuracy=0.1))) %>%
      mutate(across(matches("^risk_ratio"), ~scales::number(., accuracy=0.01)))
  }
  res %<>% mutate(across(matches("se_log$"), ~scales::number(., accuracy=0.001))) 
  return(res)
}


# Plot probability weights at end of follow-up
plot_final_weights <- function(cohort_data_long) {
  final_weight_plot <- cohort_data_long %>% 
    filter(pweight > 0) %>%
    filter(day == 83 & censor==0 & prior_censor==0) %>%
    ggplot(aes(pweight)) +
    geom_histogram(binwidth = 0.1) +
    theme_minimal() + 
    xlab("Weight") + 
    ylab("Count")
  return(final_weight_plot)
}

plot_max_weights <- function(res) {
  max_weight_plot <- res %>% bind_rows() %>% filter(day==1 & analysis=="Fully weighted") %>% 
    ggplot(aes(max_pweight)) + geom_histogram(binwidth=1) + 
    theme_minimal() + 
    xlab("Weight") + 
    ylab("Count") 
  return(max_weight_plot)
}

create_weighted_table <- function(cohort_data_long) {
  svy_obj <-
    survey::svydesign(
      id = ~pid,
      weights = ~pweight,
      data = cohort_data_long %>% filter(pweight > 0 & censor==0 & prior_censor==0 & day == 83),
    )
  weighted_table <- svy_obj %>%
    tbl_svysummary(by=pretty_init_immediately, include=c(female, age, days_from_adm, hiv_vl_suppressed, cd4_log, hb_log, creat_log, rifresist, bac_load_grm, impaired_conscious, hypoxia)) %>% 
    add_difference(everything() ~ "smd") %>% 
    remove_abbreviation() %>%
    modify_column_hide(conf.low) 
  return(weighted_table)
}

calculate_gee_risk_ratio <- function(cohort_data_long) {
  res <- cohort_data_long %>% filter(pweight >0) %>% 
    filter(day == 83 & prior_censor == 0 & censor==0) %>% 
    mutate(outcome = pmax(outcome, prior_outcome)) %>% 
    {geeglm(outcome ~ init_immediately, id=pid, family=poisson, data=., weights=.[["pweight"]])} %>% 
    tidy(conf.int=TRUE, exponentiate=TRUE)
  return(res)
}


create_risk_forest_plot <- function(res) {
  res %<>%
    mutate(mean=str_extract(risk_ratio, "([0-9.]+)"), lower=str_extract(risk_ratio, "(?<=\\()[0-9.]+"),
           upper=str_extract(risk_ratio, "(?<=, )[0-9.]+")) %>%
    mutate(across(starts_with("risk"), ~str_replace_all(.," \\(", "\\\n\\("))) %>%
    mutate(across(c(mean, lower, upper), as.numeric)) %>%
    rename(`Analysis`=analysis, `Risk delayed (%)`=risk_delayed, `Risk immediate (%)`=risk_immediate, `Risk difference (%)`=risk_difference,
           `Risk ratio`=risk_ratio) 
    
    res %<>% select(`Analysis`, `Risk delayed (%)`, `Risk immediate (%)`, `Risk difference (%)`, `Risk ratio`, mean, lower, upper)
    forest_rr <- res %>% forestplot(labeltext=c(`Analysis`, `Risk delayed (%)`, `Risk immediate (%)`, `Risk difference (%)`, `Risk ratio`), align=c("l", "l","l","l"),colgap=unit(3, "mm"), boxsize=0.1, xlog=TRUE, xticks=c(0.2, 1, 2), vertices=TRUE, clip=c(0.2, 2)) %>%
      fp_add_header(`Analysis`="Analysis", `Risk delayed (%)` = "Risk\ndelayed - %\n(95% CI)",
                    `Risk immediate (%)` = "Risk\nimmediate - %\n(95% CI)", `Risk difference (%)` = "Risk\ndifference - %\n(95% CI)", `Risk ratio` = "Risk\nratio")

  forest_rr %<>%
    fp_set_style(box="black", line="black", txt_gp=fpTxtGp(summary = gpar(fontfamily = "", cex = 0.7),
                                                           label = gpar(fontfamily = "", cex = 0.7),ticks = gpar(fontfamily = "", cex = 0.8))) %>%
    fp_set_zebra_style("#F5F9F9")
  return(forest_rr)
}
