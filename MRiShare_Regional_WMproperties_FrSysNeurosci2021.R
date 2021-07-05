library(stats)
library(tidyverse)
library(broom)
library(effectsize)
library(afex)
library(emmeans)
library(modelr)
library(glue)
library(grid)
library(ggpubr)
library(scales)
library(ggtext)
library(gt)
library(gtsummary)
library(GGally)
options(na.action = na.warn)

setwd("/Users/amitsuchida/Dropbox (GIN)/Ami-work/MRiShare/NODDI_paper/20210201_FrSystNSci_Tsuchida_MRI-Share_NODDI/Rcode_data_for_Frontiers_submission/")

dir.create("Results")

### 1) Load data 

# Primary data
u26_dat <- read_csv("data/MRiShare_WM_IDPs_FrontSysNeurosci2021.csv",
                    col_types = cols(Sex=readr::col_factor(levels = c("F", "M"))))

# Add centered data for the main covariates and classify columns for easier manipulation
#####
# Add centered Age, eTIV, WM volume and squared age column
u26_dat <- u26_dat %>%                                          
  mutate(Age_c = as.numeric(scale(Age, scale=FALSE)),
         eTIV_c = as.numeric(scale(eTIV, scale=FALSE)),
         wmV_c = as.numeric(scale(SPM_WM_Volume, scale=FALSE)),
         Sq_Age_c = as.numeric(Age_c^2)) %>%                                      
  select(ID, Sex, Age, eTIV, SPM_WM_Volume,
         Age_c, eTIV_c, wmV_c, Sq_Age_c, everything())  

# Classify columns to make it easy to manipulate
covar_cols <- colnames(u26_dat)[2:9]
dep_cols <- colnames(u26_dat)[10:225]
qc_cols <- colnames(u26_dat)[226:263]
#####

# Useful information for JHU ROIs (used for plotting purposes)
mori_ordered_27 <-read_csv("data/JHU_27ROI_info.csv",
                           col_type = cols(abb_tract_name = readr::col_factor(NULL),
                                           mori_group = readr::col_factor(
                                             levels = c("Brainstem", "Projection", "Association", "Commissural")),
                                           size_group = readr::col_factor(
                                             levels = c("XS", "S", "M", "L", "XL"))))
# Some useful variables and ROI abbreviation table
#####
# For various plots
size_ordered_p <- mori_ordered_27 %>%
  arrange(size_group, desc(atlas_vol_hemi_averaged))

size_group_count <- mori_ordered_27 %>%
  group_by(size_group) %>%
  summarize(n = n())

mori_ordered_p <- mori_ordered_27 %>%
  arrange(mori_group, desc(abb_tract_name))

mori_group_count <- mori_ordered_27 %>%
  group_by(mori_group) %>%
  summarize(n = n())

# Save the tract name/abbreviations as a table
roi_tbl <- mori_ordered_27 %>%
  select(tract_name, abb_tract_name, mori_group) %>%
  mutate(tract_name = stringr::str_replace_all(tract_name, "_", " "),
         tract_name = stringr::str_replace_all(tract_name, "ant", "anterior"),
         tract_name = stringr::str_replace_all(tract_name, "post", "posterior"),
         tract_name = stringr::str_replace_all(tract_name, "sup", "superior"),
         tract_name = stringr::str_replace_all(tract_name, "inf", "inferior"),
         tract_name = stringr::str_replace_all(tract_name, " cc", " corpus callosum"),
         tract_name = stringr::str_replace_all(tract_name, "internal ", "of the internal "),
         tract_name = stringr::str_replace_all(tract_name, "fronto occipital", "fronto-occipital")) %>%
  group_by(mori_group) %>%
  gt(rowname_col = "tract_name") %>%
  cols_label(abb_tract_name = "abbreviation") %>%
  gtsave(filename = "JHU_ROI_abbreviations.html",
         path = glue("{getwd()}/Results/"))
#####

# Create a long form DF for the primary data
eight_metrics <- c("volume", "FA", "MD", "AD", "RD",
                   "NDI", "ODI", "ISOVF")

long_u26 <- u26_dat %>%
  pivot_longer(
    cols = all_of(dep_cols),
    names_to = c("abb_tract_name", "metric"),
    names_sep = "_",
    values_to = "value"
  ) %>%
  mutate(abb_tract_name = factor(abb_tract_name, levels = mori_ordered_27$abb_tract_name),
         metric = factor(metric, levels = eight_metrics)) %>%
  select(ID, all_of(covar_cols), metric, abb_tract_name, value, everything())

  
### 2) Descriptive stats and QC ################################################
# Before the actual analyses, summarize demographic data and descriptive stats
# for both volume and DTI/NODDI data.
# 
# We check for outliers based on IQR method for 
# volume as well as DTI/NODDI data, and create columns w/o outliers.
# In the case of volume, we replace volume=0 to NA before computing outliers.
#
# Some useful functions used in this section
# -Fxn to return outlier upper and lower values 
getIQROutlierLims <- function(data, iqr_thr=3.0) {
  lowerq = unname(quantile(data, na.rm = TRUE)[2])
  upperq = unname(quantile(data, na.rm = TRUE)[4])
  iqr = upperq - lowerq #Or use IQR(data)
  threshold.upper = (iqr * iqr_thr) + upperq
  threshold.lower = lowerq - (iqr * iqr_thr)
  return(list("upper"=threshold.upper, "lower"=threshold.lower))
}
#
# -Fxn to get the number of outliers for a given data 
numIQROutliers <- function(data, iqr_thr=3.0) {
  thr <- getIQROutlierLims(data, iqr_thr)
  result <- length(which(data > thr$upper | data < thr$lower))
}
#
# -Fxn to make the plot of distributions, ordered/grouped by size or Mori
plotTractDistSex <- function(long_df, val, val_label, title,
                             OL_iqr = NULL, fix_scale=FALSE,
                             color="grey", grouping="mori",
                             save_fname = NULL, save_path = NULL) {
  if (grouping == "mori") {
    groupby_col = "mori_group"
    group_count = mori_group_count
    # reverse abb_tract_name within group
    ordered_info = mori_ordered_p
  } else if (grouping == "size") {
    groupby_col = "size_group"
    group_count = size_group_count
    ordered_info = size_ordered_p
  }
  
  if (!is.null(OL_iqr)) {
    dat <- long_df %>%
      group_by(abb_tract_name) %>%
      mutate(outliers = case_when({{val}} > getIQROutlierLims({{val}}, iqr_thr = OL_iqr)$upper ~ {{val}},
                                  {{val}} < getIQROutlierLims({{val}}, iqr_thr = OL_iqr)$lower ~ {{val}},
                                  TRUE ~ as.numeric(NA)))
  }
  
  facet_scale = if (fix_scale) "free_y" else "free"
  p <- dat %>%
    left_join(ordered_info, by = "abb_tract_name") %>%
    group_by(.data[[groupby_col]], abb_tract_name) %>%
    ggplot(mapping = aes(abb_tract_name, {{val}}, fill = Sex)) +
    labs(title = title, x = "JHU ROI", y = val_label) +
    geom_boxplot(outlier.colour = color, outlier.alpha = 0.5, na.rm = TRUE) +
    geom_point(mapping = aes(abb_tract_name, outliers, color = Sex),
               position = position_dodge(width = 0.9), 
               color = "firebrick3", na.rm = TRUE) +
    facet_wrap(~.data[[groupby_col]], scales = facet_scale, ncol = 1,strip.position = "right") +
    coord_flip() +
    theme_linedraw() +
    theme(strip.background = element_rect(fill = "white", color="black",size=1),
          strip.text = element_text(size=10, color="black"),
          axis.text = element_text(size = 10),
          axis.line = element_line(colour = "black"),
          panel.border = element_rect(colour = "black", fill=NA, size=1))
  
  gt = ggplotGrob(p)
  gt$heights[7] = unit(group_count$n[[1]], "null")
  gt$heights[11] = unit(group_count$n[[2]], "null")
  gt$heights[15] = unit(group_count$n[[3]], "null")
  gt$heights[19] = unit(group_count$n[[4]], "null")
  if (grouping == "size") {
    gt$heights[23] = unit(group_count$n[[5]], "null")
  }
  grid.newpage()
  
  if (!is.null(save_fname)) {
    ggsave(filename = save_fname,
           path = save_path,
           gt, width = 6, height = 10, units = "in", dpi = 300)
  }
  return(gt)
}
##############################################################################

### Dir to save Descriptive stats/QC info
desc_dir <- "Results/Descriptive_stats_and_QC"
dir.create(desc_dir)

### 2-1) Make a table of age/sex in this sample
#####
demog_tbl <- u26_dat %>%
  select(Age, Sex) %>%
  tbl_summary(by = Sex,
              type = all_continuous() ~ "continuous2",
              statistic = all_continuous() ~ c("{mean} ({sd})", 
                                               "{min}, {max}")) %>%
  add_p(test = everything() ~ "t.test") %>%
  add_overall() %>%
  as_gt() %>%
  gtsave(filename = "demographic_summary.html",
         path = glue("{getwd()}/{desc_dir}"))
#####

### 2-2) Summarize descriptive stats and create a column without outliers
# For volume data, remove those who have volume of 0, then for each
# metric and abb_tract_name, create a column without outliers
long_u26 <- long_u26 %>%
  mutate(value = if_else((metric == "volume" & value == 0.0), as.numeric(NA), value)) %>%
  group_by(metric, abb_tract_name) %>%
  mutate(noIQR3OL_val = case_when(value > getIQROutlierLims(value, iqr_thr = 3.0)$upper ~ as.numeric(NA),
                                  value < getIQROutlierLims(value, iqr_thr = 3.0)$lower ~ as.numeric(NA),
                                  TRUE ~ value))  %>%
  select(ID, all_of(covar_cols), metric, abb_tract_name, value, noIQR3OL_val, everything())

# Save descriptive stats for non NA vals
#####
summary_u26 <- long_u26 %>%
  group_by(metric, abb_tract_name) %>%
  summarize(n = sum(!is.na(value)),
            min = min(value, na.rm = TRUE),
            max = max(value, na.rm = TRUE),
            mean = mean(value, na.rm = TRUE),
            sd = sd(value, na.rm = TRUE),
            median = median(value, na.rm = TRUE),
            iqr = IQR(value, na.rm = TRUE),
            iqr_outliers = numIQROutliers(value, iqr_thr = 3.0),
            shapiro_stat = shapiro.test(value)$statistic,
            shapiro_p = shapiro.test(value)$p.value) %>%
  write_csv(glue("{desc_dir}/descriptive_stats.csv"))

# Create and save gt
desc_tbl <- summary_u26 %>%
  left_join(mori_ordered_27, by = "abb_tract_name") %>%
  mutate(
    "scaled_mean" = case_when(
      metric == "volume" ~ mean*10^-3,
      metric %in% c("MD", "RD", "AD") ~ mean*10^4,
      TRUE ~ mean),
    "fmt_mean" = if_else(
      metric %in% c("FA", "NDI", "ODI", "ISOVF"),
      format(scaled_mean, digit = 0, nsmall = 2, scientific = FALSE),
      format(scaled_mean, digit = 0, nsmall = 1, scientific = FALSE)
    ),
    "scaled_sd" = case_when(
      metric == "volume" ~ sd*10^-3,
      metric %in% c("MD", "RD", "AD") ~ sd*10^4,
      TRUE ~ sd),
    "fmt_sd" = if_else(
      metric %in% c("FA", "NDI", "ODI", "ISOVF"),
      format(scaled_sd, digit = 0, nsmall = 2, scientific = FALSE),
      format(scaled_sd, digit = 0, nsmall = 1, scientific = FALSE)
    ),
    "scaled_min" = case_when(
      metric == "volume" ~ min*10^-3,
      metric %in% c("MD", "RD", "AD") ~ min*10^4,
      TRUE ~ min),
    "fmt_min" = if_else(
      metric %in% c("FA", "NDI", "ODI", "ISOVF"),
      format(scaled_min, digit = 0, nsmall = 2, scientific = FALSE),
      format(scaled_min, digit = 0, nsmall = 1, scientific = FALSE)
    ),
    "scaled_max" = case_when(
      metric == "volume" ~ max*10^-3,
      metric %in% c("MD", "RD", "AD") ~ max*10^4,
      TRUE ~ max),
    "fmt_max" = if_else(
      metric %in% c("FA", "NDI", "ODI", "ISOVF"),
      format(scaled_max, digit = 0, nsmall = 2, scientific = FALSE),
      format(scaled_max, digit = 0, nsmall = 1, scientific = FALSE)
    )
  ) %>%
  pivot_wider(id_cols = c(mori_group, abb_tract_name),
              names_from = metric,
              values_from = starts_with("fmt_"),
              names_sep = "_") %>%
  arrange(factor(abb_tract_name, levels = mori_ordered_27$abb_tract_name)) %>%
  group_by(mori_group) %>%
  gt(rowname_col = "abb_tract_name") %>%
  tab_spanner(
    label = html("Mean (SD)<br>[min, max]"),
    columns = starts_with("fmt_")
  ) %>%
  cols_merge(
    columns = vars(fmt_mean_volume, fmt_sd_volume, fmt_min_volume, fmt_max_volume),
    hide_columns = vars(fmt_sd_volume, fmt_min_volume, fmt_max_volume),
    pattern = "{1} ({2})<br>[{3}, {4}]" 
  ) %>%
  cols_merge(
    columns = vars(fmt_mean_FA, fmt_sd_FA, fmt_min_FA, fmt_max_FA),
    hide_columns = vars(fmt_sd_FA, fmt_min_FA, fmt_max_FA),
    pattern = "{1} ({2})<br>[{3}, {4}]" 
  ) %>%
  cols_merge(
    columns = vars(fmt_mean_MD, fmt_sd_MD, fmt_min_MD, fmt_max_MD),
    hide_columns = vars(fmt_sd_MD, fmt_min_MD, fmt_max_MD),
    pattern = "{1} ({2})<br>[{3}, {4}]" 
  ) %>%
  cols_merge(
    columns = vars(fmt_mean_AD, fmt_sd_AD, fmt_min_AD, fmt_max_AD),
    hide_columns = vars(fmt_sd_AD, fmt_min_AD, fmt_max_AD),
    pattern = "{1} ({2})<br>[{3}, {4}]" 
  ) %>%
  cols_merge(
    columns = vars(fmt_mean_RD, fmt_sd_RD, fmt_min_RD, fmt_max_RD),
    hide_columns = vars(fmt_sd_RD, fmt_min_RD, fmt_max_RD),
    pattern = "{1} ({2})<br>[{3}, {4}]" 
  ) %>%
  cols_merge(
    columns = vars(fmt_mean_NDI, fmt_sd_NDI, fmt_min_NDI, fmt_max_NDI),
    hide_columns = vars(fmt_sd_NDI, fmt_min_NDI, fmt_max_NDI),
    pattern = "{1} ({2})<br>[{3}, {4}]" 
  ) %>%
  cols_merge(
    columns = vars(fmt_mean_ISOVF, fmt_sd_ISOVF, fmt_min_ISOVF, fmt_max_ISOVF),
    hide_columns = vars(fmt_sd_ISOVF, fmt_min_ISOVF, fmt_max_ISOVF),
    pattern = "{1} ({2})<br>[{3}, {4}]" 
  ) %>%
  cols_merge(
    columns = vars(fmt_mean_ODI, fmt_sd_ODI, fmt_min_ODI, fmt_max_ODI),
    hide_columns = vars(fmt_sd_ODI, fmt_min_ODI, fmt_max_ODI),
    pattern = "{1} ({2})<br>[{3}, {4}]" 
  ) %>%
  cols_label(
    abb_tract_name = "ROI",
    fmt_mean_volume = html("Volume <br>(x10<sup>3</sup> mm<sup>3</sup>)"),
    fmt_mean_FA = "FA",
    fmt_mean_MD = html("MD <br>(x10<sup>-4</sup> mm<sup>2</sup>/sec)"),
    fmt_mean_AD = html("AD <br>(x10<sup>-4</sup> mm<sup>2</sup>/sec)"),
    fmt_mean_RD = html("RD <br>(x10<sup>-4</sup> mm<sup>2</sup>/sec)"),
    fmt_mean_NDI = "NDI",
    fmt_mean_ODI = "ODI",
    fmt_mean_ISOVF = "IsoVF"
  ) %>%
  gtsave(filename = "descriptive_stats_summary.html",
         path = glue("{getwd()}/{desc_dir}"))
#####

# Plot the distributions for volume and WM diff metrics, volume ordered by size
# and the rest ordered by Mori group.
for (metric_name in eight_metrics) {
  metric_dat <- long_u26 %>%
    filter(metric == metric_name)
  
  grouping = if (metric_name == "volume") "size" else "mori"
  fix_scale = if (metric_name == "volume") FALSE else TRUE
  plotTractDistSex(metric_dat,
                   value,
                   grouping = grouping,
                   OL_iqr = 3.0,
                   title = glue("Summary of JHU tract ROI {metric_name} distributions"),
                   val_label=metric_name,
                   fix_scale = fix_scale,
                   save_fname = glue("JHU_{metric_name}_distributions_bySex.png"),
                   save_path = desc_dir)
}

### 3) Primary analysis #######################################################
# For the primary analyses, we will describe the age and sex effects for JHU ROI
# volume as well as each of DTI/NODDI metric (7 metrics), for each tract (27 ROIs). 
#   
#  All models will contain Age*Sex interaction term.. 
#
#  Thus, the model will be;
#  y ~ Age_c + Sex + Age_c:Sex 
#
#  Perform the analysis without removing the outliers, but later perform 
# comparisons for how the removal impacts the outcome.
#
# Functions used in this section
# -Fxn to get tidied results for a specific model {model_name}_model with 
# standardized estimate (based on refit, robust and two-sd) and generalized eta 
# squared (but requires {model_name}_af in the nested data)
getTidied <- function(model_name, nested_dat, save_fname = NULL, save_path = "") {
  mod_name <- glue("{model_name}_model")
  # generalized eta use afex package
  af_name <- glue("{model_name}_af")
  ges_dat <- nested_dat %>%
    mutate(af_nice = map(!!sym(af_name), afex::nice)) %>%
    select(-data, -ends_with("model"), -ends_with("af")) %>%
    unnest(cols = af_nice) %>%
    select(Effect, ges) %>%
    rename(term = Effect) %>%
    mutate(term = stringr::str_replace(term, "Sex", "Sex1"))
  
  # standardized using effectsize package
  stdbeta_dat <- nested_dat %>%
    mutate(std_est = map(!!sym(mod_name), effectsize::standardize_parameters,
                         robust = TRUE, two_sd = TRUE)) %>%
    select(-data, -ends_with("model"), -ends_with("af")) %>%
    unnest(cols = std_est) %>%
    select(-CI) %>%
    rename(term = Parameter,
           std_estimate = Std_Coefficient,
           std_est_CI_low = CI_low,
           std_est_CI_high = CI_high)
  
  tidied <- nested_dat %>%
    arrange(metric, abb_tract_name) %>%
    mutate(tidied = map(!!sym(mod_name), broom::tidy, conf.int=TRUE)) %>%
    select(-data, -ends_with("model"), -ends_with("af")) %>%
    unnest(cols = tidied) %>%
    rename(est_CI_low = conf.low,
           est_CI_high = conf.high) %>%
    left_join(stdbeta_dat, by = c("metric", "abb_tract_name", "term")) %>%
    left_join(ges_dat, by = c("metric", "abb_tract_name", "term")) %>%
    mutate(model = model_name) %>%
    select(metric, abb_tract_name, model, term, starts_with("est"), everything()) 
  
  if (!is.null(save_fname)) {
    write_csv(tidied, glue("{save_path}/{save_fname}"))
  }
  
  return(tidied)
}
#
# -Fxn to create a table for a given term
makeTable4Term <- function(term_name, save_fname = NULL) {
  mori_group_info <- mori_ordered_27 %>%
    select(abb_tract_name, mori_group)
  
  tbl_dat <- tidied_dat %>%
    filter(term == {{term_name}}) %>%
    select(abb_tract_name, estimate, p.value, metric) %>%
    drop_na() %>%
    mutate(
      "scaled_est" = case_when(
        metric %in% c("MD", "RD", "AD") ~ estimate*10^6,
        metric %in% c("FA", "NDI", "ODI", "ISOVF") ~ estimate*10^3,
        TRUE ~ estimate),
      "fmt_estimate" = 
        format(scaled_est, digit = 0, nsmall = 1, scientific = FALSE)
    ) %>%
    mutate("fmt_estimate" = case_when(
      p.value < 0.0001 ~ glue("{fmt_estimate}***"),
      p.value < 0.001 & p.value >= 0.0001 ~ glue("{fmt_estimate}**"),
      p.value < 0.05 & p.value >= 0.001 ~ glue("{fmt_estimate}*"),
      TRUE ~ glue("{fmt_estimate}"))) %>%
    pivot_wider(id_cols = abb_tract_name,
                names_from = metric,
                values_from = c(fmt_estimate, p.value),
                names_sep = "_") %>%
    left_join(mori_group_info, by = "abb_tract_name")
  
  tbl <- tbl_dat %>%
    arrange(factor(abb_tract_name, levels = mori_group_info$abb_tract_name)) %>%
    group_by(mori_group) %>%
    gt(rowname_col = "abb_tract_name") %>%
    tab_spanner(
      label = html("<i>&beta;</i>"),
      columns = starts_with("fmt_")
    ) %>%
    cols_label(
      abb_tract_name = "ROI",
      fmt_estimate_volume = "volume",
      fmt_estimate_FA = html("FA (x10<sup>-3</sup>)"),
      fmt_estimate_MD = html("MD (x10<sup>-6</sup>)"),
      fmt_estimate_AD = html("AD (x10<sup>-6</sup>)"),
      fmt_estimate_RD = html("RD (x10<sup>-6</sup>)"),
      fmt_estimate_NDI = html("NDI (x10<sup>-3</sup>)"),
      fmt_estimate_ODI = html("ODI (x10<sup>-3</sup>)"),
      fmt_estimate_ISOVF = html("IsoVF (x10<sup>-3</sup>)")
    ) %>%
    tab_style(
      style = list(cell_text(weight = "bold")),
      locations = cells_body(
        columns = vars(fmt_estimate_volume),
        rows = p.value_volume < (0.05/216))
    ) %>%
    tab_style(
      style = list(cell_text(weight = "bold")),
      locations = cells_body(
        columns = vars(fmt_estimate_FA),
        rows = p.value_FA < (0.05/216))
    ) %>%
    tab_style(
      style = list(cell_text(weight = "bold")),
      locations = cells_body(
        columns = vars(fmt_estimate_MD),
        rows = p.value_MD < (0.05/216))
    ) %>%
    tab_style(
      style = list(cell_text(weight = "bold")),
      locations = cells_body(
        columns = vars(fmt_estimate_AD),
        rows = p.value_AD < (0.05/216))
    ) %>%
    tab_style(
      style = list(cell_text(weight = "bold")),
      locations = cells_body(
        columns = vars(fmt_estimate_RD),
        rows = p.value_RD < (0.05/216))
    ) %>%
    tab_style(
      style = list(cell_text(weight = "bold")),
      locations = cells_body(
        columns = vars(fmt_estimate_NDI),
        rows = p.value_NDI < (0.05/216))
    ) %>%
    tab_style(
      style = list(cell_text(weight = "bold")),
      locations = cells_body(
        columns = vars(fmt_estimate_ODI),
        rows = p.value_ODI < (0.05/216))
    ) %>%
    tab_style(
      style = list(cell_text(weight = "bold")),
      locations = cells_body(
        columns = vars(fmt_estimate_ISOVF),
        rows = p.value_ISOVF < (0.05/216))
    ) %>%
    cols_hide(
      columns = starts_with("p.value")
    ) %>%
    cols_move_to_end(
      columns = vars(fmt_estimate_ISOVF)
    )
  
  
  if (!is.null(save_fname)) {
    gtsave(tbl,
           filename = save_fname,
           path = glue("{getwd()}/{main_dir}"))
  }
  
  return(tbl)
}
#
# -Fxn to create a table for a given metric
makeTable4Metric <- function(metric_name, save_fname = NULL) {
  mori_group_info <- mori_ordered_27 %>%
    select(abb_tract_name, mori_group)
  
  sq_r <- glance_dat %>%
    filter(metric == {{metric_name}}) %>%
    select(abb_tract_name, adj.r.squared) %>%
    rename(sqR = adj.r.squared) %>%
    ungroup() %>%
    select(-metric)
  
  tbl_dat <- tidied_dat %>%
    filter(metric == {{metric_name}},
           term %in% c("Age_c", "Sex1", "Age_c:Sex1")) %>%
    select(abb_tract_name, metric, term, starts_with("est"), p.value, ges) %>%
    mutate(term = stringr::str_replace(term, "Age_c", "Age"),
           term = stringr::str_replace(term, "Sex1", "Sex"),
           term = stringr::str_replace(term, ":", "X")) %>%
    mutate(
      "scaled_est" = case_when(
        metric %in% c("MD", "RD", "AD") ~ estimate*10^6,
        metric %in% c("FA", "NDI", "ODI", "ISOVF") ~ estimate*10^3,
        TRUE ~ estimate),
      "fmt_estimate" = 
        format(scaled_est, digit = 0, nsmall = 1, scientific = FALSE),
      "scaled_est_CI_low" = case_when(
        metric %in% c("MD", "RD", "AD") ~ est_CI_low*10^6,
        metric %in% c("FA", "NDI", "ODI", "ISOVF") ~ est_CI_low*10^3,
        TRUE ~ est_CI_low),
      "fmt_est_CI_low" = 
        format(scaled_est_CI_low, digit = 0, nsmall = 1, scientific = FALSE),
      "scaled_est_CI_high" = case_when(
        metric %in% c("MD", "RD", "AD") ~ est_CI_high*10^6,
        metric %in% c("FA", "NDI", "ODI", "ISOVF") ~ est_CI_high*10^3,
        TRUE ~ est_CI_high),
      "fmt_est_CI_high" = 
        format(scaled_est_CI_high, digit = 0, nsmall = 1, scientific = FALSE)
    ) %>%
    mutate("fmt_estimate" = case_when(
      p.value < 0.0001 ~ glue("{fmt_estimate}***"),
      p.value < 0.001 & p.value >= 0.0001 ~ glue("{fmt_estimate}**"),
      p.value < 0.05 & p.value >= 0.001 ~ glue("{fmt_estimate}*"),
      TRUE ~ glue("{fmt_estimate}"))) %>%
    ungroup() %>%
    select(-metric, -starts_with("est"), -starts_with("scaled")) %>%
    pivot_wider(id_cols = abb_tract_name,
                names_from = term,
                values_from = c(starts_with("fmt"), p.value, ges),
                names_sep = "_") %>%
    left_join(mori_group_info, by = "abb_tract_name") %>%
    left_join(sq_r, by = "abb_tract_name")
  
  tbl <- tbl_dat %>%
    arrange(factor(abb_tract_name, levels = mori_group_info$abb_tract_name)) %>%
    group_by(mori_group) %>%
    gt(rowname_col = "abb_tract_name") %>%
    cols_merge(
      columns = vars(fmt_estimate_Age, fmt_est_CI_low_Age, fmt_est_CI_high_Age),
      hide_columns = vars(fmt_est_CI_low_Age, fmt_est_CI_high_Age),
      pattern = "{1}<br>[{2}, {3}]"
    ) %>%
    cols_merge(
      columns = vars(fmt_estimate_Sex, fmt_est_CI_low_Sex, fmt_est_CI_high_Sex),
      hide_columns = vars(fmt_est_CI_low_Sex, fmt_est_CI_high_Sex),
      pattern = "{1}<br>[{2}, {3}]"
    ) %>%
    cols_merge(
      columns = vars(fmt_estimate_AgeXSex, fmt_est_CI_low_AgeXSex, fmt_est_CI_high_AgeXSex),
      hide_columns = vars(fmt_est_CI_low_AgeXSex, fmt_est_CI_high_AgeXSex),
      pattern = "{1}<br>[{2}, {3}]"
    ) %>%
    tab_spanner(
      label = "Age",
      columns = ends_with("_Age")
    ) %>%
    tab_spanner(
      label = "Sex",
      columns = ends_with("_Sex")
    ) %>%
    tab_spanner(
      label = "Age X Sex",
      columns = ends_with("_AgeXSex")
    ) %>%
    cols_label(
      abb_tract_name = "ROI",
      fmt_estimate_Age = html("<i>&beta;</i> [95%CI]"),
      fmt_estimate_Sex = html("<i>&beta;</i> [95%CI]"),
      fmt_estimate_AgeXSex = html("<i>&beta;</i> [95%CI]"),
      p.value_Age = html("<i>p value"),
      p.value_Sex = html("<i>p value"),
      p.value_AgeXSex = html("<i>p value"),
      ges_Age = html("<i>&eta;<sup>2</sup><sub>G</sub>"),
      ges_Sex = html("<i>&eta;<sup>2</sup><sub>G</sub>"),
      ges_AgeXSex = html("<i>&eta;<sup>2</sup><sub>G</sub>"),
      sqR = html("<i>adj. R<sup>2</sup>")
    ) %>%
    fmt_number(
      columns = starts_with("p.value"),
      decimals = 4
    ) %>%
    fmt_number(
      columns = starts_with("sq"),
      decimals = 3
    ) %>%
    text_transform(
      locations = cells_body(
        columns = starts_with("p.value")),
      fn = function(x) {
        if_else(
          x > 0.1, ">0.1", x)
      }
    ) %>%
    tab_style(
      style = list(cell_text(weight = "bold")),
      locations = cells_body(
        columns = ends_with("_Age"),
        rows = p.value_Age < (0.05/216))
    ) %>%
    tab_style(
      style = list(cell_text(weight = "bold")),
      locations = cells_body(
        columns = ends_with("_Sex"),
        rows = p.value_Sex < (0.05/216))
    ) %>%
    tab_style(
      style = list(cell_text(weight = "bold")),
      locations = cells_body(
        columns = ends_with("_AgeXSex"),
        rows = p.value_AgeXSex < (0.05/216))
    )
  
  if (!is.null(save_fname)) {
    gtsave(tbl,
           filename = glue("{save_fname}_age_sex.html"),
           path = glue("{getwd()}/{main_dir}"))
  }
  
  return(tbl)
}
#
# -Fxn to plot any of the selected variable in the tidied data (std_estimate,
# t-stat, ges) for a given beta term
plotTidied <- function(tidied_dat, beta_term, var_to_plot,
                       title = "", var_label = "", var_max = NULL,
                       save_fname = NULL, save_path = NULL) {
  div_vars <- c("std_estimate", "statistic")
  
  n.breaks = if (!is.null(var_max)) 5 else NULL
  if (!is.null(var_max)) {
    if (var_to_plot %in% div_vars) {
      limits = c(-var_max, var_max)
      if (var_max > 10) {
        labels = c(sprintf("< %.0f", -var_max),
                   sprintf("%.0f", -var_max/2),
                   "0",
                   sprintf("%.0f", var_max/2),
                   sprintf("> %.0f", var_max))
      } else {
        labels = c(sprintf("< %.1f", -var_max),
                   sprintf("%.1f", -var_max/2),
                   "0",
                   sprintf("%.1f", var_max/2),
                   sprintf("> %.1f", var_max))
      }
    } else {
      limits = c(0, var_max)
      labels = c("0.00",
                 sprintf("%.2f", var_max/4),
                 sprintf("%.2f", var_max/2),
                 sprintf("%.2f", 3*(var_max/4)),
                 sprintf("> %.2f", var_max))
    }
  }
  else {
    limits = NULL
    labels = waiver()
  }
  
  beta_term_dat <- tidied_dat %>%
    filter(term == beta_term) %>%
    mutate(ges = stringr::str_replace(ges, "<.001", "0.000"),
           ges = as.numeric(ges))
  
  p <- beta_term_dat %>%
    mutate(filtered_var = ifelse(p.value <= (0.05/216), !!sym({{var_to_plot}}), 0)) %>%
    left_join(mori_ordered_p, by="abb_tract_name") %>%
    mutate(abb_tract_name = factor(
      abb_tract_name,
      levels = c(as.character(mori_ordered_p$abb_tract_name)))) %>%
    group_by(mori_group, abb_tract_name) %>%
    ggplot(mapping = aes(x = metric, y = abb_tract_name)) +
    geom_tile(mapping = aes(fill = filtered_var)) +
    facet_wrap(~mori_group, ncol=1, scales="free_y", strip.position = "right") +
    labs(title = title, x = "", y = "", fill = var_label) +
    theme(plot.title = element_markdown(),
          legend.title = element_markdown(),
          axis.text.x = element_text(size=12, angle = 45, hjust = 1),
          axis.text.y = element_text(size=12))
  
  if (var_to_plot %in% div_vars) { 
    p <- p + scale_fill_gradient2(high = "#d01c8b", mid = "white", low = "#4dac26",
                                  midpoint = 0, na.value = "grey", limits = limits,
                                  n.breaks = n.breaks, labels = labels, oob = scales::squish)
  } else {
    p <- p + scale_fill_gradient(low = "white", high = muted("blue"), na.value = "grey",
                                 limits = limits, n.breaks = n.breaks,
                                 labels = labels, oob = scales::squish) 
  }
  
  g_count <- mori_ordered_27 %>%
    group_by(mori_group) %>%
    summarize(n = n())
  
  gt = ggplotGrob(p)
  gt$heights[7] = unit(g_count$n[[1]], "null")
  gt$heights[11] = unit(g_count$n[[2]], "null")
  gt$heights[15] = unit(g_count$n[[3]], "null")
  gt$heights[19] = unit(g_count$n[[4]], "null")
  
  grid.newpage()
  
  if (!is.null(save_fname)) {
    ggsave(filename = save_fname, path = save_path,
           gt, width = 4.2, height = 6, units = "in", dpi = 300)
  }
  return(gt)
  
}
#
# -Fxn to plot age effect for a given tract (abb_tract_name), together with 
# predicted trajectory for males and females separately
plotAgeScatter <- function(nested_dat, metric, abb_tract_name, yval_label,
                           scaling_factor = 1, yval = "value",
                           title = NULL, save_fname = NULL) {
  dat <- nested_dat %>%
    filter(metric == {{metric}},
           abb_tract_name == {{abb_tract_name}}) %>%
    select(abb_tract_name, data) %>%
    unnest(data)
  
  # Extract model as a list of length 1
  mod_list <- nested_dat %>%
    filter(metric == {{metric}},
           abb_tract_name == {{abb_tract_name}}) %>%
    '[['('main_model')
  
  # data grid 
  d_grid <- ref_grid(
    mod_list[[1]],
    at = list(Age_c = seq(-3.65, 4.2, length.out = 20))
  )
  
  p_dat <- d_grid %>%
    emmip(Age_c~Sex,
          CIs=TRUE, plotit=FALSE) %>%
    mutate(Age = mean(dat$Age) + Age_c)
  
  p <- p_dat %>%
    ggplot(aes(x=Age, y=yvar*scaling_factor)) +
    geom_point(aes(x=Age, y=get(yval)*scaling_factor, fill=Sex), 
               shape = 21, stroke = 0,
               data = dat, alpha = 0.2, inherit.aes = FALSE) +
    geom_line(aes(color=Sex)) +
    geom_ribbon(aes(ymin=LCL*scaling_factor, ymax=UCL*scaling_factor, color=Sex), alpha=0.2, size=0.1) +
    scale_fill_manual(values = c("F" = "coral1", "M" = "cyan3")) +
    scale_color_manual(values = c("F" = "firebrick4", "M" = "darkslategrey"))
  
  
  # if scaling factor is used, control the number formatting in y axis
  if (scaling_factor != 1) {
    p <- p + scale_y_continuous(
      labels = scales::number_format(accuracy = 0.1))
  }
  
  # add title
  title = if (is.null(title)) glue("{metric} in {abb_tract_name}") else title
  p <- p + labs(title = title, y = yval_label, x = "Age (in years)") +
    theme_linedraw() +
    theme(axis.title.y = ggtext::element_markdown(),
          axis.text=element_text(size=10))
  
  if (!is.null(save_fname)) {
    ggsave(plot = p,
           filename = save_fname,
           path = save_path,
           device="png",
           width=4,
           height=3,
           units="in",
           dpi=300)
  }
  return(p)
}
###############################################################################
main_dir <- "Results/Primary_analyses"
dir.create(main_dir)

### 3-1) Gather data and apply models
#####
# Apply the main model
nested_dat <- long_u26 %>%
  group_by(metric, abb_tract_name) %>%
  nest() %>%
  mutate(
    main_model = map(data, ~lm(value ~ Age_c + Sex + Age_c:Sex,
                               contrasts = list(Sex="contr.sum"),
                               data = .))
  )
######

### 3-2) Get model info and save
#####
glance_dat <- nested_dat %>%
  arrange(metric, abb_tract_name) %>%
  mutate(glance = map(main_model, broom::glance)) %>%
  select(-data, -main_model) %>%
  unnest(cols = glance) %>%
  write_csv(glue("{main_dir}/model_info_summary.csv"))
#####

### 3-3) Get tidied results and save
# First add afex anova equivalent of the model for getting generalized etasq with
# getTidied fxn. Since adding afex anova is memory expensive, delete after use if
# necessary.
#####
nested_af <- nested_dat %>%
  group_by(metric, abb_tract_name) %>%
  mutate(
    main_af = map(data, ~aov_car(value ~ Age_c + Sex + Age_c:Sex
                                 + Error(ID),
                                 data = .,
                                 factorize = FALSE,
                                 observed = c("Age_c", "Sex")))
  )

tidied_dat <- getTidied(model_name = "main",
                        nested_dat = nested_af,
                        save_fname = "model_results_summary.csv",
                        save_path = main_dir)
rm(nested_af)

# check for significance of interaction terms
AxS_sig <- tidied_dat %>%
  filter(term == "Age_c:Sex1") %>%
  arrange(p.value)

# not significant. Minimum uncorrected p is 0.0008 for SCP ODI.
#####

## For age, compute percent change from age 18 and save a table
#####
pred18_dat <- nested_dat %>%
  mutate(d_grid = map(main_model, emmeans::ref_grid,
                      at = list(Age_c = -3.65)),
         pred = map(d_grid, emmeans::emmip,
                    formula = Age_c~Sex,
                    plotit = FALSE)) %>%
  select(metric, abb_tract_name, pred) %>%
  unnest(pred) %>%
  select(metric, abb_tract_name, Sex, yvar) 

pred18_acrossSex <- pred18_dat %>%
  group_by(metric, abb_tract_name) %>%
  summarize(pred18 = mean(yvar))

age_perc_change <- tidied_dat %>%
  filter(term == "Age_c") %>%
  select(metric, abb_tract_name, estimate, p.value) %>%
  left_join(pred18_acrossSex, by=c("metric", "abb_tract_name")) %>%
  mutate(perc_change = (estimate/pred18)*100)

apc_by_mori_group <- age_perc_change %>%
  left_join(mori_ordered_27, by=c("abb_tract_name")) %>%
  group_by(metric, mori_group) %>%
  summarize(mean_18val = mean(pred18),
            min_18val = min(pred18),
            max_18val = max(pred18),
            mean_est = mean(estimate),
            min_est = min(estimate),
            max_est = max(estimate),
            mean_apc = mean(perc_change),
            min_apc = min(perc_change),
            max_apc = max(perc_change))

# Age perc change table
age_perc_table <- age_perc_change %>%
  mutate(
    "fmt_perc_change" = 
      format(perc_change, digit = 0, nsmall = 2, scientific = FALSE)
  ) %>%
  mutate("fmt_perc_change" = case_when(
    p.value < 0.0001 ~ glue("{fmt_perc_change}***"),
    p.value < 0.001 & p.value >= 0.0001 ~ glue("{fmt_perc_change}**"),
    p.value < 0.05 & p.value >= 0.001 ~ glue("{fmt_perc_change}*"),
    TRUE ~ glue("{fmt_perc_change}"))) %>%
  pivot_wider(id_cols = abb_tract_name,
              names_from = metric,
              values_from = c(fmt_perc_change, p.value),
              names_sep = "_") %>%
  left_join(mori_ordered_27, by = "abb_tract_name") %>%
  select(mori_group, abb_tract_name, starts_with("fmt"), starts_with("p.val")) %>%
  arrange(factor(abb_tract_name, levels = mori_ordered_27$abb_tract_name)) %>%
  group_by(mori_group) %>%
  gt(rowname_col = "abb_tract_name") %>%
  cols_label(
    abb_tract_name = "ROI",
    fmt_perc_change_volume = "volume",
    fmt_perc_change_FA = html("FA"),
    fmt_perc_change_MD = html("MD"),
    fmt_perc_change_AD = html("AD"),
    fmt_perc_change_RD = html("RD"),
    fmt_perc_change_NDI = html("NDI"),
    fmt_perc_change_ODI = html("ODI"),
    fmt_perc_change_ISOVF = html("IsoVF")
  ) %>%
  tab_style(
    style = list(cell_text(weight = "bold")),
    locations = cells_body(
      columns = vars(fmt_perc_change_volume),
      rows = p.value_volume < (0.05/216))
  ) %>%
  tab_style(
    style = list(cell_text(weight = "bold")),
    locations = cells_body(
      columns = vars(fmt_perc_change_FA),
      rows = p.value_FA < (0.05/216))
  ) %>%
  tab_style(
    style = list(cell_text(weight = "bold")),
    locations = cells_body(
      columns = vars(fmt_perc_change_MD),
      rows = p.value_MD < (0.05/216))
  ) %>%
  tab_style(
    style = list(cell_text(weight = "bold")),
    locations = cells_body(
      columns = vars(fmt_perc_change_AD),
      rows = p.value_AD < (0.05/216))
  ) %>%
  tab_style(
    style = list(cell_text(weight = "bold")),
    locations = cells_body(
      columns = vars(fmt_perc_change_RD),
      rows = p.value_RD < (0.05/216))
  ) %>%
  tab_style(
    style = list(cell_text(weight = "bold")),
    locations = cells_body(
      columns = vars(fmt_perc_change_NDI),
      rows = p.value_NDI < (0.05/216))
  ) %>%
  tab_style(
    style = list(cell_text(weight = "bold")),
    locations = cells_body(
      columns = vars(fmt_perc_change_ODI),
      rows = p.value_ODI < (0.05/216))
  ) %>%
  tab_style(
    style = list(cell_text(weight = "bold")),
    locations = cells_body(
      columns = vars(fmt_perc_change_ISOVF),
      rows = p.value_ISOVF < (0.05/216))
  ) %>%
  cols_hide(
    columns = starts_with("p.value")
  ) %>%
  cols_move_to_end(
    columns = vars(fmt_perc_change_ISOVF)
  ) %>%
  gtsave(filename = "age_perc_change.html",
         path = glue("{getwd()}/{main_dir}"))


#####

### 3-4) Summary of main analysis
# Gather all the results for the main model and save the following;
# 1) Tables of PE/CI/pval for each comparison
# 2) t stats for each tract/metric as a heatmap for Age and Sex
# 3) Heatmaps for effect size (ges)
# 4) Individual plots for age and sex effects for selected comparisons

##### 3-4-1) Tables
# For main table, put PE for 8 phenotypes for a given term (i.e. 5 tables)
age_tbl <- makeTable4Term("Age_c", "Age_effects_table.html")
sex_tbl <- makeTable4Term("Sex1", "Sex_effects_table.html")
ageXsextbl <- makeTable4Term("Age_c:Sex1", "AgeXSex_effects_table.html")

# Supplemental tables
for (metric_name in eight_metrics) {
  makeTable4Metric(metric_name, glue("{metric_name}_table"))
}

##### 3-4-2) t stat summary
#####
# Age
age_t <- plotTidied(tidied_dat = tidied_dat,
                    beta_term = "Age_c",
                    var_to_plot = "statistic",
                    var_label = "*t*",
                    title = "(A) *t* statistics for age",
                    save_fname = "Age_t_stat.png",
                    save_path = main_dir)

# Sex
sex_t <- plotTidied(tidied_dat = tidied_dat,
                    beta_term = "Sex1",
                    var_to_plot = "statistic",
                    var_label = "*t*",
                    var_max = 20.0,
                    title = "(A) *t* statistics for sex",
                    save_fname = "Sex_t_stat.png",
                    save_path = main_dir)
#####

##### 3-4-3) Effect size summary
#####
# Age
age_eta <- plotTidied(tidied_dat = tidied_dat,
                      beta_term = "Age_c",
                      var_to_plot = "ges",
                      var_label = "*&eta;<sup>2</sup><sub>G</sub>*",
                      title = "(B) *&eta;<sup>2</sup><sub>G</sub>* for age",
                      save_fname = "Age_ges.png",
                      save_path = main_dir)

# Sex
sex_eta <- plotTidied(tidied_dat = tidied_dat,
                      beta_term = "Sex1",
                      var_to_plot = "ges",
                      var_label = "*&eta;<sup>2</sup><sub>G</sub>*",
                      var_max = 0.2,
                      title = "(B) *&eta;<sup>2</sup><sub>G</sub>* for sex",
                      save_fname = "Sex_ges.png",
                      save_path = main_dir)

#####

##### 3-4-4) Individual plots of age/sex effects
#####
age_scatter_dir = glue("{main_dir}/Age_scatter_plots")
dir.create(age_scatter_dir)

# We will show some example age scatter plots for the following ROI and metrics;
# 1) CgH: volume, MD, NDI
# 2) UNC: volume, MD, NDI
# 3) PLIC: volume, AD, ODI 
# 4) SCP: volume, AD, ODI -- legend
# create individual plots, and arrange them to save a single image per ROI.

## 1) CgH
cgh_vol_age <- plotAgeScatter(nested_dat,
                              metric = "volume",
                              abb_tract_name = "CgH",
                              yval_label = glue("Volume (mm<sup>3</sup>)"),
                              scaling_factor = 1,
                              yval = "value",
                              title = "Volume")

cgh_md_age <- plotAgeScatter(nested_dat,
                             metric = "MD",
                             abb_tract_name = "CgH",
                             yval_label = glue("MD (x10<sup>-4</sup> mm<sup>2</sup>/sec)"),
                             scaling_factor = 10^4,
                             yval = "value",
                             title = "MD")

cgh_ndi_age <- plotAgeScatter(nested_dat,
                             metric = "NDI",
                             abb_tract_name = "CgH",
                             yval_label = glue("NDI"),
                             scaling_factor = 1,
                             yval = "value",
                             title = "NDI")

cgh_fig <- ggarrange(cgh_vol_age, cgh_md_age, cgh_ndi_age,
                     nrow = 1, legend = "none")

ggsave(glue("{age_scatter_dir}/CgH_fig_poster.png"),
       cgh_fig,
       device="png",
       width=6,
       height=2,
       units="in",
       dpi=300)

## 2) UNC
unc_vol_age <- plotAgeScatter(nested_dat,
                             metric = "volume",
                             abb_tract_name = "UNC",
                             yval_label = glue("Volume (mm<sup>3</sup>)"),
                             scaling_factor = 1,
                             yval = "value",
                             title = "Volume")

unc_md_age <- plotAgeScatter(nested_dat,
                            metric = "MD",
                            abb_tract_name = "UNC",
                            yval_label = glue("MD (x10<sup>-4</sup> mm<sup>2</sup>/sec)"),
                            scaling_factor = 10^4,
                            yval = "value",
                            title = "MD")

unc_ndi_age <- plotAgeScatter(nested_dat,
                             metric = "NDI",
                             abb_tract_name = "UNC",
                             yval_label = glue("NDI"),
                             scaling_factor = 1,
                             yval = "value",
                             title = "NDI")

unc_fig <- ggarrange(unc_vol_age, unc_md_age, unc_ndi_age,
                     nrow = 1, legend = "none")

ggsave(glue("{age_scatter_dir}/UNC_fig.png"),
       unc_fig,
       device="png",
       width=6.6,
       height=2,
       units="in",
       dpi=300)

## 3) PLIC
plic_vol_age <- plotAgeScatter(nested_dat,
                              metric = "volume",
                              abb_tract_name = "PLIC",
                              yval_label = glue("Volume (x10<sup>3</sup> mm<sup>3</sup>)"),
                              scaling_factor = 10^-3,
                              yval = "value",
                              title = "Volume")

plic_ad_age <- plotAgeScatter(nested_dat,
                             metric = "AD",
                             abb_tract_name = "PLIC",
                             yval_label = glue("AD (x10<sup>-4</sup> mm<sup>2</sup>/sec)"),
                             scaling_factor = 10^4,
                             yval = "value",
                             title = "AD")

plic_odi_age <- plotAgeScatter(nested_dat,
                              metric = "ODI",
                              abb_tract_name = "PLIC",
                              yval_label = glue("ODI"),
                              scaling_factor = 1,
                              yval = "value", 
                              title = "ODI")

plic_fig <- ggarrange(plic_vol_age, plic_ad_age, plic_odi_age,
                     nrow = 1, legend = "none")

ggsave(glue("{age_scatter_dir}/PLIC_fig.png"),
       plic_fig,
       device="png",
       width=6.6,
       height=2,
       units="in",
       dpi=300)

## 4) SCP
scp_vol_age <- plotAgeScatter(nested_dat,
                              metric = "volume",
                              abb_tract_name = "SCP",
                              yval_label = glue("Volume (mm<sup>3</sup>)"),
                              scaling_factor = 1,
                              yval = "value",
                              title = "Volume")

scp_ad_age <- plotAgeScatter(nested_dat,
                             metric = "AD",
                             abb_tract_name = "SCP",
                             yval_label = glue("AD (x10<sup>-4</sup> mm<sup>2</sup>/sec)"),
                             scaling_factor = 10^4,
                             yval = "value",
                             title = "AD")

scp_odi_age <- plotAgeScatter(nested_dat,
                              metric = "ODI",
                              abb_tract_name = "SCP",
                              yval_label = glue("ODI"),
                              scaling_factor = 1,
                              yval = "value", 
                              title = "ODI")

scp_fig <- ggarrange(scp_vol_age, scp_ad_age, scp_odi_age,
                     nrow = 1, common.legend = TRUE, legend = "bottom")

ggsave(glue("{age_scatter_dir}/SCP_fig.png"),
       scp_fig,
       device="png",
       width=6.6,
       height=2.5,
       units="in",
       dpi=300)

# 5) Plot all age effects for each ROI per metric
for (metric_name in eight_metrics){
  metric_nested_dat <- nested_dat %>%
    filter(metric == metric_name) 

  sf <- if (metric_name %in% c("MD", "AD", "RD")) 10^4 else 1
  
  metric_plots <- c()
  print(metric_name)
  for (t_name in mori_ordered_27$abb_tract_name) {
    print(t_name)
    
    metric_plots[[t_name]] <- plotAgeScatter(metric_nested_dat,
                                           metric = metric_name,
                                           abb_tract_name = t_name,
                                           yval_label = "",
                                           scaling_factor = sf,
                                           yval = "value",
                                           title = t_name)

  }
  
  metric_plot <- ggarrange(plotlist = metric_plots,
                           ncol=5, nrow = 6,
                           common.legend = TRUE, legend = "bottom")
  
  ggsave(glue("{age_scatter_dir}/{metric_name}_age_scatter.png"),
         metric_plot,
         device="png",
         width=8.5,
         height=11,
         units="in",
         dpi=300)
}
#####

### 3-5) Inter-relationships between the age effects ###########################
# Explore the relationships between patterns of age effects in the ROIs 
# Additional functions and variables used in this section
# -color groups used in the plots
group_colors =  c("Brainstem" = "dodgerblue",
                  "Projection" = "magenta2",
                  "Association" = "green3",
                  "Commissural" = "gold")
## -Fxns to plot DTINODDI std beta corplot for a given term
clean_breaks <- function(data, mapping, breaks, lims, ...){ 
  ggplot(data = data, mapping = mapping, ...) + 
    geom_point(...) + 
    geom_vline(xintercept = 0, linetype = 2, color = "grey") +
    geom_hline(yintercept = 0, linetype = 2, color = "grey") +
    scale_y_continuous(breaks = breaks,
                       limits = lims) +
    scale_x_continuous(breaks = breaks,
                       limits = lims) + 
    theme(panel.background = element_blank())
}

clean_breaks_diag <- function(data, mapping, lims, ...){ 
  ggplot(data = data, mapping = mapping, ...) + 
    geom_density() +
    geom_vline(xintercept = 0, linetype = 2, color = "grey") +
    scale_x_continuous(limits = lims)  + 
    theme(panel.background = element_blank())
}

cor_func <- function(data, mapping, method, symbol, ...){
  x <- eval_data_col(data, mapping$x)
  y <- eval_data_col(data, mapping$y)
  
  corr <- cor.test(x, y, method=method)
  est <- corr$estimate
  stars <- c("***", "**", "*", "")[findInterval(corr$p.value, c(0, 0.0001, 0.001, 0.05, 1))]
  lbl <- paste0(round(est, 2), stars)
  colFn <- colorRampPalette(c("dodgerblue", "white", "brown1"), 
                            interpolate ='spline')
  fill <- colFn(100)[findInterval(est, seq(-1, 1, length = 100))]
  ff <- if (corr$p.value < 0.05/28) 2 else 1
  
  ggally_text(
    label = lbl, 
    mapping = ggplot2::aes(size=32),
    xP = 0.5, yP = 0.5,
    color = 'black',
    fontface = ff,
    ...
  ) + #removed theme_void()
    theme(panel.background = element_rect(fill = fill))
}

make_corplot <- function(wide_dat,
                         save_fname,
                         breaks = c(-0.4, 0, 0.4),
                         lims = c(-0.5, 0.5),
                         save_path=int_dir) {

  corplot <- ggpairs(wide_dat,
                     columns = eight_metrics,
                     lower = list(mapping = aes(color = mori_group),
                                  continuous = wrap(clean_breaks, breaks = breaks, lims = lims)),
                     diag = list(continuous = wrap(clean_breaks_diag, lims = lims)),
                     upper = list(continuous = wrap(cor_func,
                                                    method = 'pearson'))) +
    scale_color_manual(values = group_colors) +
    theme(strip.background = element_rect(fill = "white"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text = element_text(size = 10),
          axis.line = element_line(colour = "black"),
          panel.border = element_rect(colour = "black", fill=NA, size=1),
          legend.position = "bottom")
  
  corplot[1, 1] <- corplot[1, 1] +
    theme(axis.text.y = element_blank(), 
          axis.ticks = element_blank())
  
  corplot[8, 8] <- corplot[8, 8] +
    theme(axis.text.x = element_blank(), 
          axis.ticks = element_blank())
  
  ggsave(corplot,
         filename = save_fname,
         path = save_path,
         device = "png",
         width = 8, height = 8, units = "in", dpi = 300)
  
}

stdBeta_corplot <- function(term,
                            tidied_dat,
                            save_fname = "stdBeta_corplot.png",
                            save_path = int_dir) {
  wide_dat <- tidied_dat %>%
      filter(term == {{term}}) %>%
      drop_na() %>%
      select(metric, abb_tract_name, std_estimate) %>%
      pivot_wider(names_from = metric, values_from = std_estimate) %>%
      left_join(mori_ordered_27, by = c("abb_tract_name"))

  #plot with default break and lims
  make_corplot(wide_dat,
               save_fname = save_fname,
               save_path = save_path)

}
################################################################################

stdBeta_corplot(term = "Age_c", tidied_dat, save_fname = "Age_stdBeta_corplot.jpeg")
stdBeta_corplot(term = "Sex1", tidied_dat, save_fname = "Sex_stdBeta_corplot.jpeg")

# Corr plot of standardized mean
#####
std_mean_wide <- long_u26 %>%
  filter(!is.na(value)) %>%
  group_by(metric) %>%
  mutate(std_value = as.numeric(scale(value, scale=T))) %>%
  group_by(metric, abb_tract_name) %>%
  summarize(mean = mean(std_value, na.rm = TRUE)) %>%
  ungroup() %>%
  left_join(mori_ordered_27, by="abb_tract_name") %>%
  select(metric, abb_tract_name, mori_group, mean) %>%
  pivot_wider(names_from = metric, values_from = mean)

make_corplot(std_mean_wide,
             breaks = c(-2, 0, 2), 
             lims = c(-3.1, 3.1),
             save_fname = "ROI_std_mean_val_corplot.jpeg",
             save_path = int_dir)
#####

# To control for sex effect, create a similar corr plot using the intercept 
# from standardized (across ROIs) value
#####
std_nested_dat_af <- long_u26 %>%
  group_by(metric) %>%
  mutate(std_value = as.numeric(scale(value, scale=T))) %>%
  group_by(metric, abb_tract_name) %>%
  nest() %>%
  mutate(
    main_model = map(data, ~lm(std_value ~ Age_c + Sex + Age_c:Sex,
                               contrasts = list(Sex="contr.sum"),
                               data = .)),
    main_af = map(data, ~aov_car(std_value ~ Age_c + Sex + Age_c:Sex
                                  + Error(ID),
                                  data = .,
                                  factorize = FALSE,
                                observed = c("Age_c", "Sex")))
    )

std_tidied_dat <- getTidied(model_name = "main",
                            nested_dat = std_nested_dat_af)

std_intercept_beta <- std_tidied_dat %>%
  filter(str_detect(term, "Intercept")) %>%
  select(metric, abb_tract_name, estimate) %>%
  pivot_wider(names_from = metric, values_from = estimate) %>%
  left_join(mori_ordered_27, by = c("abb_tract_name"))

make_corplot(wide_dat = std_intercept_beta,
             breaks = c(-2, 0, 2), 
             lims = c(-3.1, 3.1),
             save_fname = "stdIntercept_Beta_corplot.png",
             save_path = int_dir)
#####

### 4) Check for non-linear trends using quadratic model #######################
# Here I will compare linear vs quadratic age effect models to check for 
# non-linearity of age effects in both WM volume and DWI metrics
#
# -models to be compared;
#   1) lin: linear age effect (main model)
#   2) quad: additional quadratic age effect
#
# More functions used in this and the following sections
# -Fxn to compare multiple models in tidied data for a given beta term, and for
# any of the selected variable in the tidied data (estimate, std_estimate, t-stat)
plotTidiedModelComp <- function(tidied_dat, beta_term, var_to_plot,
                                title = "",
                                title_hjust= 0,
                                var_label = "",
                                var_max = NULL,
                                model_labels = NULL,
                                multi_metric = FALSE,
                                save_fname = NULL,
                                save_path = "") {
  div_vars <- c("std_estimate", "statistic")
  limits = if (!is.null(var_max)) c(-var_max, var_max) else NULL
  n.breaks = if (!is.null(var_max)) 5 else NULL
  if (!is.null(var_max)) {
    if (var_max > 10) {
      labels = c(sprintf("< %.0f", -var_max),
                 sprintf("%.0f", -var_max/2),
                 "0",
                 sprintf("%.0f", var_max/2),
                 sprintf("> %.0f", var_max))
    } else {
      labels = c(sprintf("< %.1f", -var_max),
                 sprintf("%.1f", -var_max/2),
                 "0",
                 sprintf("%.1f", var_max/2),
                 sprintf("> %.1f", var_max))
    }
  } else {
    labels = waiver()
  }
  
  beta_term_dat <- tidied_dat %>%
    filter(term == beta_term) %>%
    mutate(ges = stringr::str_replace(ges, "<.001", "0.000"),
           ges = as.numeric(ges))
  
  p <- beta_term_dat %>%
    mutate(filtered_var = ifelse(p.value <= (0.05/216), !!sym({{var_to_plot}}), 0)) %>%
    left_join(mori_ordered_p, by="abb_tract_name") %>%
    mutate(abb_tract_name = factor(
      abb_tract_name,
      levels = c(as.character(mori_ordered_p$abb_tract_name)))) %>%
    group_by(mori_group, abb_tract_name) %>%
    ggplot(mapping = aes(x = model, y = abb_tract_name)) +
    geom_tile(mapping = aes(fill = filtered_var)) +
    labs(title = title, x = "", y = "", fill = var_label) +
    theme(plot.title = element_markdown(hjust = title_hjust),
          legend.title = element_markdown(),
          axis.text.x = element_text(size=12, angle = 45, hjust = 1),
          axis.text.y = element_text(size=12))
  
  if (var_to_plot %in% div_vars) { 
    p <- p + scale_fill_gradient2(high = "#d01c8b", mid = "white", low = "#4dac26",
                                  midpoint = 0, na.value = "grey", limits = limits,
                                  n.breaks = n.breaks, labels = labels, oob = scales::squish)
  } else {
    p <- p + scale_fill_gradient(low = "white", high = muted("blue"), na.value = "grey",
                                 limits = limits, n.breaks = n.breaks,
                                 labels = labels, oob = scales::squish) 
  }
  
  if (!is.null(model_labels)) { 
    p <- p + scale_x_discrete(labels = model_labels)
  }
  
  if (multi_metric) {
    p <- p + facet_grid(mori_group~metric, scales = "free", space="free", drop = TRUE)
    save_w = if (length(unique(beta_term_dat$model)) < 4) 6 else 8
  } else {
    p <- p + facet_grid(mori_group~., scales = "free", space="free", drop = TRUE)
    save_w = if (length(unique(beta_term_dat$model)) < 4) 3 else 3.5
  }
  
  
  if (!is.null(save_fname)) {
    ggsave(glue("{save_path}/{save_fname}"),
           p, width = save_w, height = 6, units = "in", dpi = 300)
  }
  return(p)
  
}

# -Fxn to plot rel bic or adj sq R from respective data to compare multiple models
plotModelQualityComp <- function(quality_dat, quality_metric = "relative_BIC", 
                                 title = "", title_hjust =0, q_max = NULL, model_labels = NULL,
                                 multi_metric = FALSE, save_fname = NULL, save_path = "") {
  
  scale_x_labels = if (!is.null(model_labels)) model_labels else waiver()
  
  if (quality_metric == "relative_BIC") {
    fill_grad_low = "white"
    fill_grad_high = "forestgreen"
    fill_lab = "relative BIC"
  } else if (quality_metric == "sqR") {
    fill_grad_low = "white"
    fill_grad_high = "darkmagenta"
    fill_lab = "adj *R<sup>2</sup>*"
  }
  
  limits = if (!is.null(q_max)) c(0, q_max) else NULL
  n.breaks = if (!is.null(q_max)) 5 else NULL
  if (!is.null(q_max)) {
    labels = c("0.0",
               sprintf("%.1f", q_max/4),
               sprintf("%.1f", (q_max/4)*2),
               sprintf("%.1f", (q_max/4)*3),
               sprintf("> %.1f", q_max)
    ) 
  } else {
    labels = waiver()
  }
  
  p <- quality_dat %>%
    left_join(mori_ordered_p, by="abb_tract_name") %>%
    mutate(
      abb_tract_name = factor(
        abb_tract_name,
        levels = c(as.character(mori_ordered_p$abb_tract_name))),
    ) %>%
    ggplot(mapping = aes(x = model, y = abb_tract_name)) +
    geom_tile(mapping = aes_string(fill = quality_metric)) +
    scale_x_discrete(labels=scale_x_labels) +
    scale_fill_gradient(low = fill_grad_low, high = fill_grad_high,
                        limits = limits, n.breaks = n.breaks,
                        labels = labels,
                        oob = scales::squish) +
    labs(title = title, x = "", y = "", fill = fill_lab) +
    theme(plot.title = element_markdown(hjust = title_hjust),
          legend.title = element_markdown(),
          axis.text.x = element_text(size=12, angle = 45, hjust = 1),
          axis.text.y = element_text(size=12))
  
  if (multi_metric) {
    p <- p + 
      facet_grid(mori_group ~ metric, scales = "free", space="free", drop = TRUE)
    save_w = if (length(unique(quality_dat$model)) < 4) 6 else 8
  } else {
    p <- p + 
      facet_grid(mori_group ~., scales = "free", space="free", drop = TRUE)
    save_w = if (length(unique(quality_dat$model)) < 4) 3 else 3.5
  }
  
  if (!is.null(save_fname)) {
    ggsave(glue("{save_path}/{save_fname}"),
           p, width = save_w, height = 6, units = "in", dpi = 300)
  }
  return(p)
}
################################################################################
quad_dir <- "Results/Lin_Quad"
dir.create(quad_dir)

models_quad <- c("lin", "quad")

### 4-1) Gather data and apply models
#####
# Apply the quad model to nonzero_vol, rename the main model to lin
nested_quad <- nested_dat %>%
  mutate(
    quad_model = map(data, ~lm(value ~ Age_c + Sex + Sq_Age_c +
                                 Age_c:Sex,
                               contrasts = list(Sex="contr.sum"),
                               data = .))
  ) %>%
  rename(lin_model = main_model)
#####

### 4-2) Get model info for quad model and save
#####
glance_quad <- nested_quad %>%
  arrange(metric, abb_tract_name) %>%
  mutate(glance = map(quad_model, broom::glance)) %>%
  select(-data, -ends_with("model")) %>%
  unnest(cols = glance) %>%
  write_csv(glue("{quad_dir}/quad_model_info_summary.csv"))

# Combine glance data from main model 
linQuad_glance_dat <- bind_rows(
  "lin" = glance_dat,
  "quad" = glance_quad,
  .id = "model"
) %>%
  rename("sqR" = adj.r.squared) %>%
  mutate(model = factor(model, levels = models_quad),
         metric = factor(metric, levels = eight_metrics))
#####

### 4-3) Get tidied results and save
#####
nested_quad_af <- nested_quad %>%
  mutate(
    quad_af = map(data, ~aov_car(value ~ Age_c + Sex + Sq_Age_c +
                                   Age_c:Sex + Error(ID),
                                 data = .,
                                 factorize = FALSE,
                                 observed = c("Sq_Age_c", "Age_c", "Sex")))
  )

tidied_quad <- getTidied(model_name = "quad",
                         nested_dat = nested_quad_af,
                         save_fname = 'quad_results_summary.csv',
                         save_path = quad_dir)

### Combine tidied data
linQuad_tidied <- bind_rows(
  "lin" = tidied_dat,
  "quad" = tidied_quad,
  .id = "model"
) %>%
  mutate(model = factor(model, levels = models_quad))
#####

### 4-4) Plot model quality comparisons
# 1) relative BIC
#####
linQuad_rel_BIC <- linQuad_glance_dat %>%
  select(model, metric, abb_tract_name, BIC) %>%
  pivot_wider(names_from = model,
              values_from = BIC) %>%
  rowwise() %>%
  mutate(
    min_BIC = min(lin, quad),
    lin_rel_BIC = (lin - min_BIC),
    quad_rel_BIC = (quad - min_BIC)
  ) %>%
  select(metric, abb_tract_name, ends_with("rel_BIC")) %>%
  pivot_longer(
    cols = contains("_rel_"),
    names_to = "model",
    values_to = "relative_BIC"
  ) %>%
  mutate(model = stringr::str_replace(model, "_rel_BIC", ""),
         model = factor(model, levels = models_quad)) 

linQuad_rel_bic_p <- plotModelQualityComp(linQuad_rel_BIC,
                                          title = "relative BIC values",
                                          title_hjust = 0,
                                          multi_metric = TRUE,
                                          save_fname = "rel_BIC.png",
                                          save_path = quad_dir)
#####

# 2) adjusted sq R -- plot separately for volume and dti/noddi
#####
linQuad_glance_vol <- linQuad_glance_dat %>%
  filter(metric == "volume")

linQuad_glance_dn <- linQuad_glance_dat %>%
  filter(metric != "volume")

linQuad_sqR_vol_p <- plotModelQualityComp(linQuad_glance_vol,
                                          title = html("adjusted R<sup>2</sup> for volume"),
                                          title_hjust = 0,
                                          quality_metric = "sqR",
                                          multi_metric = FALSE,
                                          save_fname = "volume_sqR.png",
                                          save_path = quad_dir)
linQuad_sqR_dn_p <- plotModelQualityComp(linQuad_glance_dn,
                                         title = html("adjusted R<sup>2</sup> for DTI/NODDI"),
                                         title_hjust = 0,
                                         quality_metric = "sqR",
                                         multi_metric = TRUE,
                                         save_fname = "dtinoddi_sqR.png",
                                         save_path = quad_dir)
#####

#### 5) Effects of outliers and QC metric ######################################
# Here I will compare the impact of QC-related info on the analysis.
# For WM vol and DTI/NODDI, we compare the effect of
# 1) QC and phenotypic outliers
# 2) OL removal plus adding a generic quality measure as a covariate
#
# For WM volumetry, Euler number is used as the generic covariate.
#
# For DTI/NODDI, RMS is used.
#
################################################################################

# first split the volume and DTI/NODDI data in long_u26 since we will perform 
# different analyses on these data.
long_vol <- long_u26 %>%
  filter(metric == "volume")

long_dn <- long_u26 %>%
  filter(metric != "volume")

### 5-1) DTI/NODDI
dnQCdir <- "Results/DTINODDI_QC"
dir.create(dnQCdir)

models_dnQC <- c("main", "noOL", "noOL_RMS")
dnQC_labels <- c("main", "noOL", "noOL + RMS")

#### 5-1-1) Gather data and apply models
#####
# First we need to define QC outliers for DTI/NODDI data
dwiqc_cols <- qc_cols[13:38]
dwiqc <- u26_dat %>%
  select(ID, all_of(dwiqc_cols))

long_dwiqc <- dwiqc%>%
  pivot_longer(
    cols = contains("QC"),
    names_to = c("QCtype", "QCmetric"),
    names_sep = "QC_",
    values_to = "value"
  )

# Outliers only on the upper side
dwiqc_OL <- long_dwiqc %>%
  group_by(QCmetric) %>%
  mutate(upper = getIQROutlierLims(value, iqr_thr = 3.0)$upper) %>%
  filter(value > getIQROutlierLims(value, iqr_thr = 3.0)$upper) 

dwiqc_OL_bySub <- dwiqc_OL %>%
  group_by(ID, QCtype) %>%
  summarise(n=n())

dwiqc_OLsubs <- unique(dwiqc_OL_bySub$ID)

# Now create a new nested data that removed all dwiqc_OLsubs
# also select DTI/NODDI metrics, and add models with and w/o RMS
dnQC_nested_dat <- long_dn %>%
  filter(!ID %in% dwiqc_OLsubs) %>%
  group_by(metric, abb_tract_name) %>%
  mutate(rms_c = as.numeric(scale(DWIQC_shellNS_relative_rms, scale = FALSE))) %>%
  nest() %>%
  mutate(
    noOL_model = map(data, ~lm(noIQR3OL_val ~ Age_c + Sex + Age_c:Sex,
                               contrasts = list(Sex="contr.sum"),
                               data = .)),
    noOL_RMS_model = map(data, ~lm(noIQR3OL_val ~ Age_c + Sex + rms_c + Age_c:Sex,
                                   contrasts = list(Sex="contr.sum"),
                                   data = .))
  )
#####    

#### 5-1-2) Get model info for noOL/noOL_RMS and save
#####
glance_OL <- dnQC_nested_dat %>%
  arrange(metric, abb_tract_name) %>%
  mutate(glance = map(noOL_model, broom::glance)) %>%
  select(-data, -noOL_model, -noOL_RMS_model) %>%
  unnest(cols = glance) %>%
  write_csv(glue("{dnQCdir}/noOL_model_info_summary.csv"))

glance_OL_RMS <- dnQC_nested_dat %>%
  arrange(metric, abb_tract_name) %>%
  mutate(glance = map(noOL_RMS_model, broom::glance)) %>%
  select(-data, -noOL_model, -noOL_RMS_model) %>%
  unnest(cols = glance) %>%
  write_csv(glue("{dnQCdir}/noOL_RMS_model_info_summary.csv"))

# Combine glance data from main model 
dnQC_glance_dat <- bind_rows(
  "main" = glance_dat,
  "noOL" = glance_OL,
  "noOL_RMS" = glance_OL_RMS,
  .id = "model"
) %>%
  rename("sqR" = adj.r.squared) %>%
  filter(metric != "volume")
#####

#### 5-1-3) Get tidied results and save
#####
dnQC_nested_af <- dnQC_nested_dat %>%
  mutate(
    noOL_af = map(data, ~aov_car(noIQR3OL_val ~ Age_c + Sex + Age_c:Sex
                                 + Error(ID),
                                 data = .,
                                 factorize = FALSE,
                                 observed = c("Age_c", "Sex"))),
    noOL_RMS_af = map(data, ~aov_car(noIQR3OL_val ~ Age_c + Sex + rms_c + Age_c:Sex
                                     + Error(ID),
                                     data = .,
                                     factorize = FALSE,
                                     observed = c("Age_c", "Sex", "rms_c")))
  )

tidied_OL <- getTidied(model_name = "noOL",
                       nested_dat = dnQC_nested_af,
                       save_fname = "noOL_results_summary.csv",
                       save_path = dnQCdir)

tidied_OL_RMS <- getTidied(model_name = "noOL_RMS",
                           nested_dat = dnQC_nested_af,
                           save_fname = "noOL_RMS_results_summary.csv",
                           save_path = dnQCdir)

rm(dnQC_nested_af)

# Combine tidied data from main model 
dnQC_tidied_dat <- bind_rows(
  "main" = tidied_dat,
  "noOL" = tidied_OL,
  "noOL_RMS" = tidied_OL_RMS,
  .id = "model"
) %>%
  mutate(model = factor(model, levels = models_dnQC)) %>%
  filter(metric!="volume")
#####

#### 5-1-4) Plot comparisons
#####
## adj sq R
dnQC_sqR_comp <- plotModelQualityComp(dnQC_glance_dat,
                                      title = html("(A) adjusted R<sup>2</sup> for DTI/NODDI"),
                                      title_hjust = 0,
                                      quality_metric = "sqR",
                                      multi_metric = TRUE,
                                      model_labels = dnQC_labels,
                                      save_fname = "dtinoddi_sqR.png",
                                      save_path = dnQCdir)

## age effect
dnQC_age_t_comp <- plotTidiedModelComp(dnQC_tidied_dat,
                                    beta_term = "Age_c",
                                    var_to_plot = "statistic",
                                    title = "(A) *t* statistics for age",
                                    var_label = "*t*",
                                    model_labels = dnQC_labels,
                                    multi_metric = T,
                                    save_fname = "age_t_stat.png",
                                    save_path = dnQCdir)

dnQC_age_ges_comp <- plotTidiedModelComp(dnQC_tidied_dat,
                                      beta_term = "Age_c",
                                      var_to_plot = "ges",
                                      title = "Age effect sizes",
                                      var_label = "*&eta;<sup>2</sup><sub>G</sub>*",
                                      model_labels = dnQC_labels,
                                      multi_metric = T,
                                      save_fname = "age_ges.png",
                                      save_path = dnQCdir)
## sex effect
dnQC_sex_t_comp <- plotTidiedModelComp(dnQC_tidied_dat,
                                       beta_term = "Sex1",
                                       var_to_plot = "statistic",
                                       title = "(B) *t* statistics for sex",
                                       var_label = "*t*",
                                       model_labels = dnQC_labels,
                                       multi_metric = T,
                                       save_fname = "sex_t_stat.png",
                                       save_path = dnQCdir)

dnQC_sex_ges_comp <- plotTidiedModelComp(dnQC_tidied_dat,
                                      beta_term = "Sex1",
                                      var_to_plot = "ges",
                                      title = "Sex effect sizes",
                                      var_label = "*&eta;<sup>2</sup><sub>G</sub>*",
                                      model_labels = dnQC_labels,
                                      multi_metric = T,
                                      save_fname = "sex_ges.png",
                                      save_path = dnQCdir)

######

#### 5-2) anat QC
volQCdir <- "Results/WMvol_QC"
dir.create(volQCdir)

models_volQC <- c("all", "noOL", "noOL_Euler")
volQC_labels <- c("all", "noOL", "noOL + Euler")

#### 5-2-1) Gather data and apply models
#####
# First we need to define QC outliers for volumetric data
structqc_cols <- qc_cols[c(1:12, 31)] 
structqc <- u26_dat %>%
  select(ID, all_of(structqc_cols))

# Get list of extreme IQR outliers in ANY QC metrics
long_structqc <- structqc %>%
  pivot_longer(
    cols = contains("QC"),
    names_to = c("QCtype", "QCmetric"),
    names_sep = "QC_",
    values_to = "value"
  )

structqc_OL <- long_structqc %>%
  group_by(QCmetric) %>%
  mutate(upper = getIQROutlierLims(value, iqr_thr = 3.0)$upper) %>%
  filter(value > getIQROutlierLims(value, iqr_thr = 3.0)$upper) 

structqc_OL_bySub <- structqc_OL %>%
  group_by(ID, QCtype) %>%
  summarise(n=n())

structqc_OLsubs <- unique(structqc_OL_bySub$ID) 

# Now create a new nested data that removed all structqc_OLsubs
# and add models with and w/o Euler as a covariate
volQC_nested_dat <- long_vol %>%
  filter(!ID %in% structqc_OLsubs) %>%
  group_by(metric, abb_tract_name) %>%
  mutate(Eul_c = as.numeric(scale(structQC_FS6_invEuler, scale=FALSE))) %>%
  nest() %>%
  mutate(
    noOL_model = map(data, ~lm(noIQR3OL_val ~ Age_c + Sex + Age_c:Sex,
                               contrasts = list(Sex="contr.sum"),
                               data = .)),
    noOL_Euler_model = map(data, ~lm(noIQR3OL_val ~ Age_c + Sex + Eul_c + Age_c:Sex,
                                   contrasts = list(Sex="contr.sum"),
                                   data = .))
  )
#####

#### 5-2-2) Get model info for noOL/noOL_Euler and save
#####
vol_glance_OL <- volQC_nested_dat %>%
  arrange(metric, abb_tract_name) %>%
  mutate(glance = map(noOL_model, broom::glance)) %>%
  select(-data, -noOL_model, -noOL_Euler_model) %>%
  unnest(cols = glance) %>%
  write_csv(glue("{volQCdir}/noOL_model_info_summary.csv"))

vol_glance_OL_Euler <- volQC_nested_dat %>%
  arrange(metric, abb_tract_name) %>%
  mutate(glance = map(noOL_Euler_model, broom::glance)) %>%
  select(-data, -noOL_model, -noOL_Euler_model) %>%
  unnest(cols = glance) %>%
  write_csv(glue("{volQCdir}/noOL_Euler_model_info_summary.csv"))

# Combine glance data from main model 
volQC_glance_dat <- bind_rows(
  "main" = glance_dat,
  "noOL" = vol_glance_OL,
  "noOL_Euler" = vol_glance_OL_Euler,
  .id = "model"
) %>%
  rename("sqR" = adj.r.squared) %>%
  filter(metric == "volume")
#####

#### 5-2-3) Get tidied results and save
#####
volQC_nested_af <- volQC_nested_dat %>%
  mutate(
    noOL_af = map(data, ~aov_car(noIQR3OL_val ~ Age_c + Sex + Age_c:Sex
                                 + Error(ID),
                                 data = .,
                                 factorize = FALSE,
                                 observed = c("Age_c", "Sex"))),
    noOL_Euler_af = map(data, ~aov_car(noIQR3OL_val ~ Age_c + Sex + Eul_c + Age_c:Sex
                                     + Error(ID),
                                     data = .,
                                     factorize = FALSE,
                                     observed = c("Age_c", "Sex", "Eul_c")))
  )

vol_tidied_OL <- getTidied(model_name = "noOL",
                       nested_dat = volQC_nested_af,
                       save_fname = "noOL_results_summary.csv",
                       save_path = volQCdir)

vol_tidied_OL_Euler <- getTidied(model_name = "noOL_Euler",
                           nested_dat = volQC_nested_af,
                           save_fname = "noOL_Euler_results_summary.csv",
                           save_path = volQCdir)

rm(volQC_nested_af)

# Combine tidied data from main model 
volQC_tidied_dat <- bind_rows(
  "main" = tidied_dat,
  "noOL" = vol_tidied_OL,
  "noOL_Euler" = vol_tidied_OL_Euler,
  .id = "model"
) %>%
  mutate(model = factor(model, levels = models_volQC)) %>%
  filter(metric=="volume")
#####

#### 5-2-4) Plot comparisons
#####
## adj sq R
volQC_sqR_comp <- plotModelQualityComp(volQC_glance_dat,
                                      title = html("adjusted R<sup>2</sup> for volume"),
                                      title_hjust = 0,
                                      quality_metric = "sqR",
                                      multi_metric = FALSE,
                                      model_labels = volQC_labels,
                                      save_fname = "WMvol_sqR.png",
                                      save_path = volQCdir)

## age effect
volQC_age_t_comp <- plotTidiedModelComp(volQC_tidied_dat,
                                       beta_term = "Age_c",
                                       var_to_plot = "statistic",
                                       title = "(A) *t* statistics for age",
                                       var_label = "*t*",
                                       model_labels = volQC_labels,
                                       multi_metric = F,
                                       save_fname = "age_t_stat.png",
                                       save_path = volQCdir)

volQC_age_ges_comp <- plotTidiedModelComp(volQC_tidied_dat,
                                         beta_term = "Age_c",
                                         var_to_plot = "ges",
                                         title = "Age effect sizes",
                                         var_label = "*&eta;<sup>2</sup><sub>G</sub>*",
                                         model_labels = volQC_labels,
                                         multi_metric = F,
                                         save_fname = "age_ges.png",
                                         save_path = volQCdir)
## sex effect
volQC_sex_t_comp <- plotTidiedModelComp(volQC_tidied_dat,
                                       beta_term = "Sex1",
                                       var_to_plot = "statistic",
                                       title = "(B) *t* statistics for sex",
                                       var_label = "*t*",
                                       model_labels = volQC_labels,
                                       multi_metric = F,
                                       save_fname = "sex_t_stat.png",
                                       save_path = volQCdir)

volQC_sex_ges_comp <- plotTidiedModelComp(volQC_tidied_dat,
                                         beta_term = "Sex1",
                                         var_to_plot = "ges",
                                         title = "Sex effect sizes",
                                         var_label = "*&eta;<sup>2</sup><sub>G</sub>*",
                                         model_labels = volQC_labels,
                                         multi_metric = F,
                                         save_fname = "sex_ges.png",
                                         save_path = volQCdir)

#####

#### 6) Effects of volume corrections ##########################################
# Here I will investigate the effects of global volume correction on volumetry
# and global/local volume correction on DTI/NODDI metrics.
#
# 1) Effects of global volume on regional WM volumes
# -models to be compared;
#   1) noVol: no TIV or other global volume (main model)
#   2) TIV: addition of TIV
#   3) TWMV: TWMV instead of TIV
#   4) TIV_TWMV: TIV + WMvol
#
# 2) Effects of global and regional volume on DTI/NODDI
# -models to be compared;
#   1) noVol: no TIV or other volumes (main model)
#   2) TIV: addition of TIV
#   3) ROIV: ROI vol instead of TIV
#   4) TIV_ROIV: TIV + ROI vol
#
################################################################################

### 6-1) Effects of global volume on WM volumetry
volTIV_dir <- "Results/volTIV"
dir.create(volTIV_dir)

models_volTIV <- c("noVol", "TIV", "TWMV", "TIV_TWMV")
volTIV_labels <- c("noVol", "TIV", "TWMV", "TIV + TWMV")

#### 6-1-1) Gather data and apply models
#####
# Add 3 additional models (TIV, TWMV, combined) and rename main_model to noVol
volTIV_nested_dat <- long_vol %>%
  group_by(metric, abb_tract_name) %>%
  nest() %>%
  mutate(
    TIV_model = map(data, ~lm(value ~ Age_c + Sex + eTIV_c + Age_c:Sex,
                                contrasts = list(Sex="contr.sum"),
                                data = .)),
    TWMV_model = map(data, ~lm(value ~ Age_c + Sex + wmV_c + Age_c:Sex,
                               contrasts = list(Sex="contr.sum"),
                               data = .)),
    TIV_TWMV_model = map(data, ~lm(value ~ Age_c + Sex + eTIV_c + wmV_c + 
                                   Age_c:Sex,
                                   contrasts = list(Sex="contr.sum"),
                                   data = .))
  )  
#####

#### 6-1-2) Get model info for each model and save
#####
volTIV_glance_TIV <- volTIV_nested_dat %>%
  arrange(metric, abb_tract_name) %>%
  mutate(glance = map(TIV_model, broom::glance)) %>%
  select(-data, -ends_with("model")) %>%
  unnest(cols = glance) %>%
  write_csv(glue("{volTIV_dir}/TIV_model_info_summary.csv"))

volTIV_glance_TWMV <- volTIV_nested_dat %>%
  arrange(metric, abb_tract_name) %>%
  mutate(glance = map(TWMV_model, broom::glance)) %>%
  select(-data, -ends_with("model")) %>%
  unnest(cols = glance) %>%
  write_csv(glue("{volTIV_dir}/TWMV_model_info_summary.csv"))

volTIV_glance_combined <- volTIV_nested_dat %>%
  arrange(metric, abb_tract_name) %>%
  mutate(glance = map(TIV_TWMV_model, broom::glance)) %>%
  select(-data, -ends_with("model")) %>%
  unnest(cols = glance) %>%
  write_csv(glue("{volTIV_dir}/TIV_TWMV_model_info_summary.csv"))

# Combine glance data from main model 
volTIV_glance_dat <- bind_rows(
  "noVol" = glance_dat,
  "TIV" = volTIV_glance_TIV,
  "TWMV" = volTIV_glance_TWMV,
  "TIV_TWMV" = volTIV_glance_combined,
  .id = "model"
) %>%
  rename("sqR" = adj.r.squared) %>%
  filter(metric == "volume")
#####

#### 6-1-3) Get tidied results and save
#####
volTIV_nested_af <- volTIV_nested_dat %>%
  mutate(
    TIV_af = map(data, ~aov_car(value ~ Age_c + Sex + eTIV_c + Age_c:Sex
                                 + Error(ID),
                                 data = .,
                                 factorize = FALSE,
                                 observed = c("Age_c", "Sex", "eTIV_c"))),
    TWMV_af = map(data, ~aov_car(value ~ Age_c + Sex + wmV_c + Age_c:Sex
                                  + Error(ID),
                                  data = .,
                                  factorize = FALSE,
                                  observed = c("Age_c", "Sex", "wmV_c"))),
    TIV_TWMV_af = map(data, ~aov_car(value ~ Age_c + Sex + eTIV_c + wmV_c
                                     + Age_c:Sex + Error(ID),
                                    data = .,
                                    factorize = FALSE,
                                    observed = c("Age_c", "Sex", "eTIV_c", "wmV_c")))
  )

volTIV_tidied_TIV <- getTidied(model_name = "TIV",
                           nested_dat = volTIV_nested_af,
                           save_fname = "TIV_results_sumamry.csv",
                           save_path = volTIV_dir)

volTIV_tidied_TWMV <- getTidied(model_name = "TWMV",
                               nested_dat = volTIV_nested_af,
                               save_fname = "TWMV_results_sumamry.csv",
                               save_path = volTIV_dir)

volTIV_tidied_combined <- getTidied(model_name = "TIV_TWMV",
                               nested_dat = volTIV_nested_af,
                               save_fname = "TIV_TWMV_results_sumamry.csv",
                               save_path = volTIV_dir)

rm(volTIV_nested_af)

# Combine tidied data from main model 
volTIV_tidied_dat <- bind_rows(
  "noVol" = tidied_dat,
  "TIV" = volTIV_tidied_TIV,
  "TWMV" = volTIV_tidied_TWMV,
  "TIV_TWMV" = volTIV_tidied_combined,
  .id = "model"
) %>%
  mutate(model = factor(model, levels = models_volTIV)) %>%
  filter(metric=="volume")
#####

#### 6-1-4) Plot comparisons
#####
## adj sq R
volTIV_sqR_comp <- plotModelQualityComp(volTIV_glance_dat,
                                       title = html("(A) adjusted R<sup>2</sup> for volume"),
                                       title_hjust = 0,
                                       quality_metric = "sqR",
                                       multi_metric = FALSE,
                                       model_labels = volTIV_labels,
                                       save_fname = "WMvol_sqR.png",
                                       save_path = volTIV_dir)

## age effect
volTIV_age_t_comp <- plotTidiedModelComp(volTIV_tidied_dat,
                                        beta_term = "Age_c",
                                        var_to_plot = "statistic",
                                        title = "(B) *t* statistics for age",
                                        var_label = "*t*",
                                        model_labels = volTIV_labels,
                                        multi_metric = F,
                                        save_fname = "age_t_stat.png",
                                        save_path = volTIV_dir)

volTIV_age_ges_comp <- plotTidiedModelComp(volTIV_tidied_dat,
                                          beta_term = "Age_c",
                                          var_to_plot = "ges",
                                          title = "Age effect sizes",
                                          var_label = "*&eta;<sup>2</sup><sub>G</sub>*",
                                          model_labels = volTIV_labels,
                                          multi_metric = F,
                                          save_fname = "age_ges.png",
                                          save_path = volTIV_dir)
## sex effect
volTIV_sex_t_comp <- plotTidiedModelComp(volTIV_tidied_dat,
                                        beta_term = "Sex1",
                                        var_to_plot = "statistic",
                                        title = "(C) *t* statistics for sex",
                                        var_label = "*t*",
                                        model_labels = volTIV_labels,
                                        multi_metric = F,
                                        save_fname = "sex_t_stat.png",
                                        save_path = volTIV_dir)

volTIV_sex_ges_comp <- plotTidiedModelComp(volTIV_tidied_dat,
                                          beta_term = "Sex1",
                                          var_to_plot = "ges",
                                          title = "Sex effect sizes",
                                          var_label = "*&eta;<sup>2</sup><sub>G</sub>*",
                                          model_labels = volTIV_labels,
                                          multi_metric = F,
                                          save_fname = "sex_ges.png",
                                          save_path = volTIV_dir)

#####

### 6-2) Effects of global/local volume on regional DTI/NODDI metrics
dnVol_dir <- "Results/dnVol"
dir.create(dnVol_dir)

models_dnVol <- c("noVol", "TIV", "ROIV", "TIV_ROIV")
dnVol_labels <- c("noVol", "TIV", "ROIV", "TIV + ROIV")

#### 6-2-1) Gather data and apply models
#####
# Create ROI vol dat to add to DWI metric DF
ROIvol_dat <- long_vol %>%
  group_by(metric, abb_tract_name) %>%
  mutate(ROIvol_c = as.numeric(scale(value, scale = FALSE))) %>%
  ungroup() %>%
  select(ID, abb_tract_name, value, ROIvol_c) %>%
  rename(ROIvol = value)

# Add 3 additional models (TIV, ROIV, combined) 
dnVol_nested_dat <- long_dn %>%
  left_join(ROIvol_dat, by=c("ID", "abb_tract_name")) %>%
  group_by(metric, abb_tract_name) %>%
  nest() %>%
  mutate(
    TIV_model = map(data, ~lm(value ~ Age_c + Sex + eTIV_c + Age_c:Sex,
                              contrasts = list(Sex="contr.sum"),
                              data = .)),
    ROIV_model = map(data, ~lm(value ~ Age_c + Sex + ROIvol_c + Age_c:Sex,
                               contrasts = list(Sex="contr.sum"),
                               data = .)),
    TIV_ROIV_model = map(data, ~lm(value ~ Age_c + Sex + eTIV_c + ROIvol_c + 
                                     Age_c:Sex,
                                   contrasts = list(Sex="contr.sum"),
                                   data = .))
  )  
#####

#### 6-2-2) Get model info for each model and save
#####
dnVol_glance_TIV <- dnVol_nested_dat %>%
  arrange(metric, abb_tract_name) %>%
  mutate(glance = map(TIV_model, broom::glance)) %>%
  select(-data, -ends_with("model")) %>%
  unnest(cols = glance) %>%
  write_csv(glue("{dnVol_dir}/TIV_model_info_summary.csv"))

dnVol_glance_ROIV <- dnVol_nested_dat %>%
  arrange(metric, abb_tract_name) %>%
  mutate(glance = map(ROIV_model, broom::glance)) %>%
  select(-data, -ends_with("model")) %>%
  unnest(cols = glance) %>%
  write_csv(glue("{dnVol_dir}/ROIV_model_info_summary.csv"))

dnVol_glance_combined <- dnVol_nested_dat %>%
  arrange(metric, abb_tract_name) %>%
  mutate(glance = map(TIV_ROIV_model, broom::glance)) %>%
  select(-data, -ends_with("model")) %>%
  unnest(cols = glance) %>%
  write_csv(glue("{dnVol_dir}/TIV_ROIV_model_info_summary.csv"))

# Combine glance data from main model 
dnVol_glance_dat <- bind_rows(
  "noVol" = glance_dat,
  "TIV" = dnVol_glance_TIV,
  "ROIV" = dnVol_glance_ROIV,
  "TIV_ROIV" = dnVol_glance_combined,
  .id = "model"
) %>%
  rename("sqR" = adj.r.squared) %>%
  filter(metric != "volume") %>%
  mutate(model = factor(model, levels = models_dnVol))
#####

#### 6-2-3) Get tidied results and save
#####
dnVol_nested_af <- dnVol_nested_dat %>%
  mutate(
    TIV_af = map(data, ~aov_car(value ~ Age_c + Sex + eTIV_c + Age_c:Sex
                                + Error(ID),
                                data = .,
                                factorize = FALSE,
                                observed = c("Age_c", "Sex", "eTIV_c"))),
    ROIV_af = map(data, ~aov_car(value ~ Age_c + Sex + ROIvol_c + Age_c:Sex
                                 + Error(ID),
                                 data = .,
                                 factorize = FALSE,
                                 observed = c("Age_c", "Sex", "ROIvol_c"))),
    TIV_ROIV_af = map(data, ~aov_car(value ~ Age_c + Sex + eTIV_c + ROIvol_c
                                     + Age_c:Sex + Error(ID),
                                     data = .,
                                     factorize = FALSE,
                                     observed = c("Age_c", "Sex", "eTIV_c", "ROIvol_c")))
  )

dnVol_tidied_TIV <- getTidied(model_name = "TIV",
                               nested_dat = dnVol_nested_af,
                               save_fname = "TIV_results_sumamry.csv",
                               save_path = dnVol_dir)

dnVol_tidied_ROIV <- getTidied(model_name = "ROIV",
                                nested_dat = dnVol_nested_af,
                                save_fname = "ROIV_results_sumamry.csv",
                                save_path = dnVol_dir)

dnVol_tidied_combined <- getTidied(model_name = "TIV_ROIV",
                                    nested_dat = dnVol_nested_af,
                                    save_fname = "TIV_ROIV_results_sumamry.csv",
                                    save_path = dnVol_dir)

rm(dnVol_nested_af)

# Combine tidied data from main model 
dnVol_tidied_dat <- bind_rows(
  "noVol" = tidied_dat,
  "TIV" = dnVol_tidied_TIV,
  "ROIV" = dnVol_tidied_ROIV,
  "TIV_ROIV" = dnVol_tidied_combined,
  .id = "model"
) %>%
  filter(metric !="volume") %>%
  mutate(model = factor(model, levels = models_dnVol))
#####

#### 6-2-4) Plot comparisons
#####
## relative BIC
dnVol_rel_BIC <- dnVol_glance_dat %>%
  select(model, metric, abb_tract_name, BIC) %>%
  pivot_wider(names_from = model,
              values_from = BIC) %>%
  rowwise() %>%
  mutate(
    min_BIC = min(noVol, TIV, ROIV, TIV_ROIV),
    noVol_rel_BIC = (noVol - min_BIC),
    TIV_rel_BIC = (TIV - min_BIC),
    ROIV_rel_BIC = (ROIV - min_BIC),
    TIV_ROIV_rel_BIC = (TIV_ROIV - min_BIC),
  ) %>%
  select(metric, abb_tract_name, ends_with("rel_BIC")) %>%
  pivot_longer(
    cols = contains("_rel_"),
    names_to = "model",
    values_to = "relative_BIC"
  ) %>%
  mutate(model = stringr::str_replace(model, "_rel_BIC", ""),
         model = factor(model, levels = models_dnVol)) 

dnVol_rel_bic_comp <- plotModelQualityComp(dnVol_rel_BIC,
                                    title = "(A) relative BIC values",
                                    title_hjust = 0,
                                    q_max = 100,
                                    quality_metric = "relative_BIC",
                                    model_labels = dnVol_labels,
                                    multi_metric = TRUE,
                                    save_fname = "rel_BIC.png",
                                    save_path = dnVol_dir)

## adj sq R
dnVol_sqR_comp <- plotModelQualityComp(dnVol_glance_dat,
                                       title = html("(B) adjusted R<sup>2</sup>"),
                                       title_hjust = 0,
                                       quality_metric = "sqR",
                                       multi_metric = TRUE,
                                       model_labels = dnVol_labels,
                                       save_fname = "DTINODDI_sqR.png",
                                       save_path = dnVol_dir)

## effect of eTIV/ROIV
eTIV_OR_ROIV_dat <- dnVol_tidied_dat %>%
  filter(term == "eTIV_c" | term == "ROIvol_c",
         model != "TIV_ROIV")

limits = c(0, 10) 
n.breaks = if (!is.null(q_max)) 5 else NULL

eTIV_ROIV_t_stats <- eTIV_OR_ROIV_dat %>%
  mutate(filtered_var = ifelse(p.value <= (0.05/216), statistic, 0)) %>%
  left_join(mori_ordered_p, by="abb_tract_name") %>%
  mutate(abb_tract_name = factor(
                              abb_tract_name,
                              levels = c(as.character(mori_ordered_p$abb_tract_name)))) %>%
  group_by(mori_group, abb_tract_name) %>%
  ggplot(mapping = aes(x = model, y = abb_tract_name)) +
  geom_tile(mapping = aes(fill = filtered_var)) +
  scale_fill_gradient2(high = "#d01c8b", mid = "white", low = "#4dac26",
                       midpoint = 0, na.value = "grey", limits = c(-20, 20),
                       n.breaks = 5, labels = c("-20", "-10", "0", "10", "20"), oob = scales::squish) +
  facet_grid(mori_group~metric, scales = "free", space="free", drop = TRUE) +
  scale_x_discrete(labels = c("TIV", "ROIV")) +
  labs(x = "", y = "", fill = "*t*", title = "(A) Effects of TIV and ROIV") +
  theme(plot.title = element_markdown(),
        legend.title = element_markdown(),
        axis.text.x = element_text(size=12, angle = 45, hjust = 1),
        axis.text.y = element_text(size=12))

ggsave(glue("{dnVol_dir}/eTIV_ROIV_t_stat.png"),
       eTIV_ROIV_t_stats, width = 6, height = 6, units = "in", dpi = 300)

## age effect
dnVol_age_t_comp <- plotTidiedModelComp(dnVol_tidied_dat,
                                        beta_term = "Age_c",
                                        var_to_plot = "statistic",
                                        title = "(B) *t* statistics for age",
                                         var_label = "*t*",
                                         model_labels = dnVol_labels,
                                         multi_metric = T,
                                         save_fname = "age_t_stat.png",
                                         save_path = dnVol_dir)

dnVol_age_ges_comp <- plotTidiedModelComp(dnVol_tidied_dat,
                                           beta_term = "Age_c",
                                           var_to_plot = "ges",
                                           title = "Age effect sizes",
                                           var_label = "*&eta;<sup>2</sup><sub>G</sub>*",
                                           model_labels = dnVol_labels,
                                           multi_metric = T,
                                           save_fname = "age_ges.png",
                                           save_path = dnVol_dir)

## sex effect
dnVol_sex_t_comp <- plotTidiedModelComp(dnVol_tidied_dat,
                                         beta_term = "Sex1",
                                         var_to_plot = "statistic",
                                         title = "(C) *t* statistics for sex",
                                         var_label = "*t*",
                                         model_labels = dnVol_labels,
                                         multi_metric = T,
                                         save_fname = "sex_t_stat.png",
                                         save_path = dnVol_dir)

dnVol_sex_ges_comp <- plotTidiedModelComp(dnVol_tidied_dat,
                                           beta_term = "Sex1",
                                           var_to_plot = "ges",
                                           title = "Sex effect sizes",
                                           var_label = "*&eta;<sup>2</sup><sub>G</sub>*",
                                           model_labels = dnVol_labels,
                                           multi_metric = T,
                                           save_fname = "sex_ges.png",
                                           save_path = dnVol_dir)
#####

### THE END!!! #################################################################

