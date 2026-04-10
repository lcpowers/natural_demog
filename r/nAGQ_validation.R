## nAGQ_validation.R
## Validates that dredge with nAGQ=0 selects the same top models as nAGQ=1.
## Uses the clim-phase global model structure from file 5 (all 6 clim vars,
## same subset constraints) on real data.
## Reference: Phillip Alday (lme4 developer), R-sig-mixed-models 2022:
##   https://stat.ethz.ch/pipermail/r-sig-mixed-models/2022q1/029941.html

library(tidyverse)
library(lme4)
library(MuMIn)

options(na.action = "na.fail")

readRDS("data/demog/vr_modData/mod_data.rds") %>% list2env(.GlobalEnv)
clim <- read_csv("data/climate/for_analyses/clim_summary.csv", show_col_types = FALSE)
clim_vars <- names(clim)[!names(clim) %in% c("site", "pop", "year")]

surv_mod_data_scaled$aspect.cat  <- relevel(as.factor(surv_mod_data_scaled$aspect.cat),  ref = "top")
surv_mod_data_scaled$site <- relevel(as.factor(surv_mod_data_scaled$site),        ref = "niwot")
growth_mod_data_scaled$aspect.cat <- relevel(as.factor(growth_mod_data_scaled$aspect.cat), ref = "top")
growth_mod_data_scaled$site <- relevel(as.factor(growth_mod_data_scaled$site),       ref = "niwot")
repro_mod_data_scaled$aspect.cat  <- relevel(as.factor(repro_mod_data_scaled$aspect.cat),  ref = "top")
repro_mod_data_scaled$site <- relevel(as.factor(repro_mod_data_scaled$site),        ref = "niwot")
nflrs_mod_data_scaled$aspect.cat  <- relevel(as.factor(nflrs_mod_data_scaled$aspect.cat),  ref = "top")
nflrs_mod_data_scaled$site <- relevel(as.factor(nflrs_mod_data_scaled$site),        ref = "niwot")

# Global formula mirrors file 5 clim phase: log_size_t0 * (year + site + aspect.cat + clim_vars)
# Simplified: drop year to reduce model space and runtime
clim_formula <- function(response) {
  as.formula(paste(c(response, "~ log_size_t0*(site + aspect.cat +",
    paste(clim_vars, collapse = " + "), ") + (1|site.pop)"), collapse = " "))
}

surv_0 <- glmer(clim_formula("surv"), family = binomial, data = surv_mod_data_scaled,  nAGQ = 0)
surv_1 <- glmer(clim_formula("surv"), family = binomial, data = surv_mod_data_scaled,  nAGQ = 1)
growth <- lmer(clim_formula("growth"),      data = growth_mod_data_scaled, REML = FALSE)
repro_0 <- glmer(clim_formula("flowered"),   family = binomial, data = repro_mod_data_scaled, nAGQ = 0)
repro_1 <- glmer(clim_formula("flowered"),   family = binomial, data = repro_mod_data_scaled, nAGQ = 1)
nflrs_0 <- glmer(clim_formula("total_flrs"), family = poisson,  data = nflrs_mod_data_scaled, nAGQ = 0)
nflrs_1  <- glmer(clim_formula("total_flrs"), family = poisson,  data = nflrs_mod_data_scaled, nAGQ = 1)

# Same subset constraints as file 5
dredge_subset <- quote(
  (ann_snow_days + win_meanT_mean + win_meanT_range + gs_snow_days + gs_meanT_mean + gs_meanT_range) <= 2 &
  !(ann_snow_days & win_meanT_mean) &
  !(ann_snow_days & gs_snow_days)
)

clust <- makeCluster(8, type = "FORK")

t_surv_0 <- system.time(d_surv_0  <- dredge(surv_0,  cluster = clust, fixed = "log_size_t0", subset = dredge_subset, trace = 2))
t_surv_1 <- system.time(d_surv_1  <- dredge(surv_1,  cluster = clust, fixed = "log_size_t0", subset = dredge_subset, trace = 2))
t_growth <- system.time(d_growth  <- dredge(growth,  cluster = clust, fixed = "log_size_t0", subset = dredge_subset, trace = 2))
t_repro_0 <- system.time(d_repro_0 <- dredge(repro_0, cluster = clust, fixed = "log_size_t0", subset = dredge_subset, trace = 2))
t_repro_1 <- system.time(d_repro_1 <- dredge(repro_1, cluster = clust, fixed = "log_size_t0", subset = dredge_subset, trace = 2))
t_nflrs_0 <- system.time(d_nflrs_0 <- dredge(nflrs_0, cluster = clust, fixed = "log_size_t0", subset = dredge_subset, trace = 2))
t_nflrs_1 <- system.time(d_nflrs_1 <- dredge(nflrs_1, cluster = clust, fixed = "log_size_t0", subset = dredge_subset, trace = 2))

stopCluster(clust)

runtime_table <- data.frame(
  vital_rate = c("survival", "survival", "growth", "repro", "repro", "nflrs", "nflrs"),
  nAGQ       = c(0, 1, NA, 0, 1, 0, 1),
  seconds    = round(c(t_surv_0["elapsed"], t_surv_1["elapsed"], t_growth["elapsed"],
                       t_repro_0["elapsed"], t_repro_1["elapsed"],
                       t_nflrs_0["elapsed"], t_nflrs_1["elapsed"]), 1)
)
print(runtime_table)
write_csv(runtime_table, "r/modelFitResults/nAGQ_validation_runtime.csv")

# Top 10 comparison
top10_tables <- function(d0, d1, vr_name) {
  t0 <- as.data.frame(d0) %>% slice_head(n = 10) %>%
    mutate(vital_rate = vr_name, nAGQ = 0, rank = row_number())
  t1 <- as.data.frame(d1) %>% slice_head(n = 10) %>%
    mutate(vital_rate = vr_name, nAGQ = 1, rank = row_number())
  bind_rows(t0, t1)
}

top10_surv  <- top10_tables(d_surv_0,  d_surv_1,  "survival")
top10_repro <- top10_tables(d_repro_0, d_repro_1, "repro")
top10_nflrs <- top10_tables(d_nflrs_0, d_nflrs_1, "nflrs")
top10_growth <- as.data.frame(d_growth) %>% slice_head(n = 10) %>%
  mutate(vital_rate = "growth", nAGQ = NA, rank = row_number())

top10_all <- bind_rows(top10_surv, top10_repro, top10_nflrs, top10_growth)
write_csv(top10_all, "r/modelFitResults/nAGQ_validation_top10.csv")

# Summary overlap table
model_overlap <- function(d0, d1, vr_name) {
  rows0 <- rownames(as.data.frame(d0))[1:10]
  rows1 <- rownames(as.data.frame(d1))[1:10]
  data.frame(vital_rate     = vr_name,
             top_model_match = rows0[1] == rows1[1],
             overlap_top10   = sum(rows0 %in% rows1))
}

overlap_summary <- bind_rows(
  model_overlap(d_surv_0,  d_surv_1,  "survival"),
  model_overlap(d_repro_0, d_repro_1, "repro"),
  model_overlap(d_nflrs_0, d_nflrs_1, "nflrs")
)
print(overlap_summary)
write_csv(overlap_summary, "r/modelFitResults/nAGQ_validation_overlap.csv")
