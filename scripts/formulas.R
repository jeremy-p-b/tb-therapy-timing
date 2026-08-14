formula_model_day1 <- "treatment_initiation ~ female + rcs(age, 3) + bac_load_grm + creat_log + hb_log + cd4_log + rifresist + impaired_conscious + hypoxia"
formula_model_day2to6 <- "treatment_initiation ~ poly(day, 3) + female + rcs(age, 3) + bac_load_grm + creat_log + hb_log + cd4_log + rifresist + impaired_conscious + hypoxia"
formula_model_day2to6_simple <- "treatment_initiation ~ day + female + rcs(age, 3) + bac_load_grm + creat_log + hb_log + cd4_log + rifresist + impaired_conscious + hypoxia"
formula_model_day2to6_rcs <- "treatment_initiation ~ rcs(day, 4) + female + rcs(age, 3) + bac_load_grm + creat_log + hb_log + cd4_log + rifresist + impaired_conscious + hypoxia"
formula_model_day2to6_factor <- "treatment_initiation ~ factor(day, levels=c(2:6)) + female + rcs(age, 3) + bac_load_grm + creat_log + hb_log + cd4_log + rifresist + impaired_conscious + hypoxia"

formula_model_day1_norifresist <- "treatment_initiation ~ female + rcs(age, 3) + bac_load_grm + creat_log + hb_log + cd4_log + impaired_conscious + hypoxia"
formula_model_day2to6_norifresist <- "treatment_initiation ~ poly(day, 3) + female + rcs(age, 3) + bac_load_grm + creat_log + hb_log + cd4_log + impaired_conscious + hypoxia"