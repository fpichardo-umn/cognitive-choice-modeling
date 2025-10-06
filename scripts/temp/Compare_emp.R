main_dir = "/Users/icd/Library/CloudStorage/Dropbox/Projects/cognitive-choice-modeling/Outputs/igt_mod/fits/fit/ahrb/ses-00"

indiv_only = "task-igt_mod_cohort-ahrb_ses-00_group-batch_001_model-ev_type-fit_desc-output.rds"
hier_only = "task-igt_mod_cohort-ahrb_ses-00_group-hier_model-ev_type-fit_desc-output.rds"
emp = "task-igt_mod_cohort-ahrb_ses-00_group-emp_model-ev_type-fit_desc-output.rds"

fit_indiv = readRDS(file.path(main_dir, indiv_only))
fit_hier = readRDS(file.path(main_dir, hier_only))
fit_emp = readRDS(file.path(main_dir, emp))

# Get subs
subs_indiv = unname(sapply(names(fit_indiv), function(index) fit_indiv[[index]]$subid))
subs_hier = fit_hier$subject_list
subs_emp = unname(sapply(names(fit_emp), function(index) fit_emp[[index]]$subid))

# Same subs
all_subs = unique(c(subs_indiv, subs_hier, subs_emp))
shared_subs = c()
for (sub in all_subs){
  bool_indiv = sub %in% subs_indiv
  bool_hier = sub %in% subs_hier
  bool_emp = sub %in% subs_emp
  if (bool_indiv & bool_hier & bool_emp){
    shared_subs[[length(shared_subs) + 1]] = sub
  }
}
shared_subs = unlist(shared_subs)

# Single sub param meds
sub = shared_subs[1]
model_params = fit_hier$model_params

## Indiv
sub_fit = fit_indiv[[which(subs_indiv == sub)]]
sub_indiv_stats = sub_fit$summary_stats %>% filter(variable %in% model_params)

## Hier
sub_hier_idx = which(subs_hier == sub)
sub_hier_stats = fit_hier$summary_stats %>% 
  filter(variable %in% paste0(model_params, "[", sub_hier_idx, "]"))  %>%
  mutate(variable = sub("\\[.*\\]", "", variable))

## Emp
sub_fit = fit_emp[[which(subs_emp == sub)]]
sub_emp_stats = sub_fit$summary_stats %>% filter(variable %in% model_params)


