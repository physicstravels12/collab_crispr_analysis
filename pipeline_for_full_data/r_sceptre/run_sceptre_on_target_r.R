suppressPackageStartupMessages({
  library(Matrix)
  library(sceptre)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) {
  root <- args[[1]]
} else {
  file_arg <- grep("^--file=", commandArgs(), value = TRUE)
  script_path <- normalizePath(sub("^--file=", "", file_arg[1]))
  root <- dirname(script_path)
}
inp <- file.path(root, "inputs")
out <- file.path(root, "outputs")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

response_matrix <- readMM(file.path(inp, "response_matrix.mtx"))
grna_matrix <- readMM(file.path(inp, "grna_matrix.mtx"))
response_names <- read.csv(file.path(inp, "response_names.csv"), stringsAsFactors = FALSE)$response_id
grna_names <- read.csv(file.path(inp, "grna_names.csv"), stringsAsFactors = FALSE)$grna_id

rownames(response_matrix) <- response_names
rownames(grna_matrix) <- grna_names

grna_target_data_frame <- read.csv(file.path(inp, "grna_target_data_frame.csv"), stringsAsFactors = FALSE)
discovery_pairs <- read.csv(file.path(inp, "discovery_pairs_on_target_expressed.csv"), stringsAsFactors = FALSE)

obj <- import_data(
  response_matrix = response_matrix,
  grna_matrix = grna_matrix,
  grna_target_data_frame = grna_target_data_frame,
  moi = "low"
)

obj <- set_analysis_parameters(
  sceptre_object = obj,
  discovery_pairs = discovery_pairs,
  side = "both",
  control_group = "nt_cells",
  resampling_mechanism = "permutations",
  multiple_testing_method = "BH",
  multiple_testing_alpha = 0.1
)

obj <- assign_grnas(
  sceptre_object = obj,
  method = "thresholding",
  threshold = 1
)

obj <- run_qc(
  sceptre_object = obj,
  response_n_umis_range = c(0.01, 0.99),
  response_n_nonzero_range = c(0.01, 0.99),
  n_nonzero_trt_thresh = 7,
  n_nonzero_cntrl_thresh = 7
)

obj <- run_discovery_analysis(
  sceptre_object = obj,
  parallel = FALSE
)

res <- get_result(obj, analysis = "run_discovery_analysis")
write.csv(res, file.path(out, "discovery_r_on_target.csv"), row.names = FALSE)

summary_tbl <- data.frame(
  n_rows = nrow(res),
  n_pass_qc = sum(res$pass_qc, na.rm = TRUE),
  n_sig = sum(res$significant, na.rm = TRUE),
  frac_neg_lfc = mean(res$log_2_fold_change < 0, na.rm = TRUE),
  mean_lfc = mean(res$log_2_fold_change, na.rm = TRUE)
)
write.csv(summary_tbl, file.path(out, "summary_r_on_target.csv"), row.names = FALSE)

cat("R run complete\n")
print(summary_tbl)
