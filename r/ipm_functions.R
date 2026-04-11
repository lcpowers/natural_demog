## ipm_functions.R
## Site × aspect × year vital rate models (climate-first then microtopo).
## Source this file at the top of IPM analysis scripts.
## Assumes working directory is the project root (natural_demog/).

library(tidyverse)
library(lme4)

##### 1. Data loading #####
sites <- c("guanella","niwot","jones")
aspects <- c("north","south","top")
years <- c(2023:2025)

mod_data   <- readRDS("data/demog/vr_modData/mod_data.rds")
demog_data <- readRDS("data/demog/vr_modData/mod_demog_data.rds")

clim <- read_csv("data/climate/for_analyses/clim_summary.csv", show_col_types = FALSE)
clim_vars <- names(clim)[!names(clim) %in% c("site","pop","year")]

microtopo <- read_csv("data/output/model_microtopo_25cm.csv", show_col_types = FALSE)
mt_vars <- names(microtopo)[!names(microtopo) %in% c("site","pop","plantID","elv_pt","prcv_25")]

scale_params <- read_csv("data/demog/vr_modData/scale_params.csv", show_col_types = FALSE)

# Germination rates by site × aspect.cat; replace zeros with 0.0005
germ_rates <- demog_data$germ_rates %>%
  mutate(aspect.cat = case_when(
    str_starts(pop, "North") ~ "north",
    str_starts(pop, "South") ~ "south",
    str_starts(pop, "Top")   ~ "top"
  )) %>%
  group_by(site, aspect.cat) %>%
  summarise(germ_rate = mean(germ_rate), .groups = "drop") %>%
  mutate(germ_rate = ifelse(germ_rate == 0, 0.0005, germ_rate))

# Default covariate values: climate by site × year, microtopo by site × aspect.cat
# Both in raw (unscaled) units; scaled inside kernel functions using scale_params.
clim_means <- mod_data$growth_mod_data %>%
  group_by(site, year) %>%
  summarise(across(all_of(clim_vars), ~mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  mutate(year = as.character(year))

mt_means <- mod_data$growth_mod_data %>%
  group_by(site, aspect.cat) %>%
  summarise(across(all_of(mt_vars), ~mean(.x, na.rm = TRUE)), .groups = "drop")

# Pop-level microtopo means for pop_kernel() / mega_kernel()
mt_means_pop <- mod_data$growth_mod_data %>%
  mutate(aspect.cat = case_when(
    str_starts(pop, "North") ~ "north",
    str_starts(pop, "South") ~ "south",
    str_starts(pop, "Top")   ~ "top"
  )) %>%
  group_by(site, pop, aspect.cat) %>%
  summarise(across(all_of(mt_vars), ~mean(.x, na.rm = TRUE)), .groups = "drop")
##### End 1 #####

##### 2. Load models #####
models <- readRDS("r/modelFitResults/best_clim_then_mt_mods_Apr9.RDS")
list2env(models, envir = .GlobalEnv)
##### End 2 #####

##### 3. IPM bin setup #####
log_sizes_all <- c(mod_data$growth_mod_data$log_size_t0,
                   mod_data$growth_mod_data$log_size_t0 + mod_data$growth_mod_data$growth)

minsize <- min(log_sizes_all, na.rm = TRUE)
maxsize <- max(log_sizes_all, na.rm = TRUE) + 0.1

n.bin <- 100
vec.bin <- seq(minsize, maxsize, length.out = n.bin + 1)
h <- vec.bin[2] - vec.bin[1]
binmids <- vec.bin[1:n.bin] + h / 2

recruit.mean.log <- demog_data$recruit_log_size
recruit.sd.log <- 0.25
sdlgszcdf <- pnorm(q = vec.bin, mean = recruit.mean.log, sd = recruit.sd.log)
sdlgszprobs <- diff(sdlgszcdf)
sdlgszprobs <- sdlgszprobs / sum(sdlgszprobs)

nutlets_per_flower <- 1.2
##### End 3 #####

##### 4. Kernel functions #####

## nat_kernel(site, aspect.cat, year)
## Returns K = P + F using observed covariate means for the given site × aspect × year.
nat_kernel <- function(site, aspect.cat, year) {
  log_size_center <- scale_params$center[scale_params$var == "log_size_t0"]
  log_size_scale  <- scale_params$scale[scale_params$var  == "log_size_t0"]

  clim_row <- clim_means %>% filter(site == !!site, year == as.character(year))
  mt_row   <- mt_means   %>% filter(site == !!site, aspect.cat == !!aspect.cat)

  pred_data <- data.frame(log_size_t0 = (binmids - log_size_center) / log_size_scale,
                          site = site, aspect.cat = aspect.cat,
                          year = factor(as.character(year)))
  for (v in clim_vars) {
    pred_data[[v]] <- (clim_row[[v]] - scale_params$center[scale_params$var == v]) /
                       scale_params$scale[scale_params$var  == v]
  }
  for (v in mt_vars) {
    pred_data[[v]] <- (mt_row[[v]] - scale_params$center[scale_params$var == v]) /
                       scale_params$scale[scale_params$var  == v]
  }

  surv_prob    <- predict(surv_clim2mt_best_mod,       newdata = pred_data, type = "response", re.form = NA)
  growth_mean  <- binmids + predict(growth_clim2mt_best_mod,    newdata = pred_data, re.form = NA)
  growth_var   <- pmax(predict(growth_var_clim2mt_best_mod, newdata = pred_data, re.form = NA), 0)
  flower_prob  <- predict(repro_clim2mt_best_mod,      newdata = pred_data, type = "response", re.form = NA)
  flower_count <- predict(nflrs_clim2mt_best_mod,      newdata = pred_data, type = "response", re.form = NA)

  site_germ_rate <- germ_rates$germ_rate[germ_rates$site == site & germ_rates$aspect.cat == aspect.cat]

  growth_mx <- matrix(0, n.bin, n.bin)
  for (ss in 1:n.bin) {
    growcdf <- pnorm(vec.bin, growth_mean[ss], sqrt(growth_var[ss]))
    grows <- diff(growcdf)
    if (sum(grows) > 0) growth_mx[, ss] <- grows / sum(grows)
  }

  P_kernel <- t(t(growth_mx) * surv_prob)
  F_kernel <- outer(sdlgszprobs, flower_prob * flower_count * nutlets_per_flower * site_germ_rate)
  P_kernel + F_kernel
}

## predict_lambda(df)
## df must have columns: site, aspect.cat, year.
## Any additional numeric predictor columns (raw/unscaled) override the default
## site × aspect × year means. Returns df with lambda column appended.
predict_lambda <- function(df) {
  log_size_center <- scale_params$center[scale_params$var == "log_size_t0"]
  log_size_scale  <- scale_params$scale[scale_params$var  == "log_size_t0"]

  df$lambda <- NA_real_

  for (i in seq_len(nrow(df))) {
    row <- df[i,]
    site_i <- row$site; aspect_i <- row$aspect.cat; year_i <- as.character(row$year)

    clim_row <- clim_means %>% filter(site == site_i, year == year_i)
    mt_row   <- mt_means   %>% filter(site == site_i, aspect.cat == aspect_i)

    pred_data <- data.frame(log_size_t0 = (binmids - log_size_center) / log_size_scale,
                            site = site_i, aspect.cat = aspect_i,
                            year = factor(year_i))
    for (v in clim_vars) {
      raw_v <- if (v %in% names(row) && !is.na(row[[v]])) row[[v]] else clim_row[[v]]
      pred_data[[v]] <- (raw_v - scale_params$center[scale_params$var == v]) /
                         scale_params$scale[scale_params$var  == v]
    }
    for (v in mt_vars) {
      raw_v <- if (v %in% names(row) && !is.na(row[[v]])) row[[v]] else mt_row[[v]]
      pred_data[[v]] <- (raw_v - scale_params$center[scale_params$var == v]) /
                         scale_params$scale[scale_params$var  == v]
    }

    surv_prob    <- predict(surv_clim2mt_best_mod,       newdata = pred_data, type = "response", re.form = NA)
    growth_mean  <- binmids + predict(growth_clim2mt_best_mod,    newdata = pred_data, re.form = NA)
    growth_var   <- pmax(predict(growth_var_clim2mt_best_mod, newdata = pred_data, re.form = NA), 0)
    flower_prob  <- predict(repro_clim2mt_best_mod,      newdata = pred_data, type = "response", re.form = NA)
    flower_count <- predict(nflrs_clim2mt_best_mod,      newdata = pred_data, type = "response", re.form = NA)

    site_germ_rate <- germ_rates$germ_rate[germ_rates$site == site_i & germ_rates$aspect.cat == aspect_i]

    growth_mx <- matrix(0, n.bin, n.bin)
    for (ss in 1:n.bin) {
      growcdf <- pnorm(vec.bin, growth_mean[ss], sqrt(growth_var[ss]))
      grows <- diff(growcdf)
      if (sum(grows) > 0) growth_mx[, ss] <- grows / sum(grows)
    }

    P_kernel <- t(t(growth_mx) * surv_prob)
    F_kernel <- outer(sdlgszprobs, flower_prob * flower_count * nutlets_per_flower * site_germ_rate)
    K <- P_kernel + F_kernel
    df$lambda[i] <- Re(eigen(K, only.values = TRUE)$values[1])
  }
  df
}
## pop_kernel(site, pop, aspect.cat, year)
## Like nat_kernel but uses pop-level microtopo means and pop-specific germ rate.
pop_kernel <- function(site, pop, aspect.cat, year) {
  log_size_center <- scale_params$center[scale_params$var == "log_size_t0"]
  log_size_scale  <- scale_params$scale[scale_params$var  == "log_size_t0"]

  clim_row <- clim_means   %>% filter(site == !!site, year == as.character(year))
  mt_row   <- mt_means_pop %>% filter(site == !!site, pop  == !!pop)

  pred_data <- data.frame(log_size_t0 = (binmids - log_size_center) / log_size_scale,
                          site = site, aspect.cat = aspect.cat,
                          year = factor(as.character(year)))
  for (v in clim_vars) {
    pred_data[[v]] <- (clim_row[[v]] - scale_params$center[scale_params$var == v]) /
                       scale_params$scale[scale_params$var  == v]
  }
  for (v in mt_vars) {
    pred_data[[v]] <- (mt_row[[v]] - scale_params$center[scale_params$var == v]) /
                       scale_params$scale[scale_params$var  == v]
  }

  surv_prob    <- predict(surv_clim2mt_best_mod,       newdata = pred_data, type = "response", re.form = NA)
  growth_mean  <- binmids + predict(growth_clim2mt_best_mod,    newdata = pred_data, re.form = NA)
  growth_var   <- pmax(predict(growth_var_clim2mt_best_mod, newdata = pred_data, re.form = NA), 0)
  flower_prob  <- predict(repro_clim2mt_best_mod,      newdata = pred_data, type = "response", re.form = NA)
  flower_count <- predict(nflrs_clim2mt_best_mod,      newdata = pred_data, type = "response", re.form = NA)

  pop_germ <- demog_data$germ_rates$germ_rate[demog_data$germ_rates$pop == pop &
                                               demog_data$germ_rates$site == site]
  if (length(pop_germ) == 0 || pop_germ == 0) pop_germ <- 0.0005

  growth_mx <- matrix(0, n.bin, n.bin)
  for (ss in 1:n.bin) {
    growcdf <- pnorm(vec.bin, growth_mean[ss], sqrt(growth_var[ss]))
    grows <- diff(growcdf)
    if (sum(grows) > 0) growth_mx[, ss] <- grows / sum(grows)
  }

  P_kernel <- t(t(growth_mx) * surv_prob)
  F_kernel <- outer(sdlgszprobs, flower_prob * flower_count * nutlets_per_flower * pop_germ)
  P_kernel + F_kernel
}

## mega_kernel(site_group, year, m)
## site_group: e.g. "guanella upper", "niwot lower", "jones east"
## Builds a (3*n.bin) × (3*n.bin) megamatrix for the 3 pops in that group
## with symmetric movement rate m between every pair of pops.
mega_kernel <- function(site_group, year, m) {
  parts   <- str_split(site_group, " ")[[1]]
  site_i  <- parts[1]
  suffix  <- str_to_title(parts[2])  # "upper"->"Upper", "east"->"East", etc.

  pop_df  <- mt_means_pop %>% filter(site == site_i, str_ends(pop, suffix))
  pops    <- pop_df$pop
  asps    <- pop_df$aspect.cat
  n_pops  <- length(pops)  # should be 3

  K_list  <- lapply(seq_len(n_pops), function(i) pop_kernel(site_i, pops[i], asps[i], year))

  N <- n.bin * n_pops
  M <- matrix(0, N, N)
  for (from in seq_len(n_pops)) {
    for (to in seq_len(n_pops)) {
      disp <- if (from == to) 1 - (n_pops - 1) * m else m
      row_idx <- ((to   - 1) * n.bin + 1):(to   * n.bin)
      col_idx <- ((from - 1) * n.bin + 1):(from * n.bin)
      M[row_idx, col_idx] <- disp * K_list[[from]]
    }
  }

  lambda <- Re(eigen(M, only.values = TRUE)$values[1])
  list(M = M, lambda = lambda, K_list = setNames(K_list, pops))
}
##### End 4 #####


test_mx <- matrix(data = 0,nrow=6,ncol=6)
test_mx[1,2] <- 0.018
test_mx[2,2] <- 0.5
test_mx[3,2] <- 0.177
test_mx[3,3] <- 0.950
test_mx[1,3] <- 7.8

test_mx[5,4] <- 0.018
test_mx[5,5] <- 0.48
test_mx[6,5] <- 0.15
test_mx[6,6] <- 0.96
test_mx[4,6] <- 7.5
test_mx
Re(eigen(test_mx)$values[1])

# movement cells
test_mx[4,3] <- 1.5
test_mx[1,6] <- 3
test_mx
Re(eigen(test_mx)$values[1])

