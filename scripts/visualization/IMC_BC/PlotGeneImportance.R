#!/usr/bin/env Rscript

# PlotGeneImportance.R
#
# Identify hub markers from BSNMani subnetworks and evaluate the survival
# association of selected markers using patient-level mean IMC expression.
#
# Usage:
#   Rscript PlotGeneImportance.R /path/to/Github_BSNMani
#
# Alternatively, set the repository root with:
#   export BSNMANI_ROOT=/path/to/Github_BSNMani
#   Rscript PlotGeneImportance.R

# -----------------------------------------------------------------------------
# 1. Configuration
# -----------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1L) {
  args[[1L]]
} else {
  Sys.getenv("BSNMANI_ROOT", unset = ".")
}
project_root <- normalizePath(project_root, mustWork = FALSE)

q <- 2L
analysis_output <- "output_22"
run_id <- "s123"

config <- list(
  project_root = project_root,
  q = q,
  n_roi = 44L,
  top_n = 10L,
  # NULL = automatically use the union of the top hub genes across subnetworks.
  # Replace NULL with a character vector to analyze a custom gene set.
  target_genes = NULL,
  model_run_dir = file.path(
    project_root, "IMC", "bc", analysis_output, paste0("q", q), run_id
  ),
  clinical_path = file.path(
    project_root, "IMC", "bc", analysis_output,
    paste0("clinical_df_coxph_q", q, "_stacked.csv")
  ),
  marker_path = file.path(
    project_root, "IMC", "Preprocessed", "Filtered_and_Renamed_Markers.csv"
  ),
  patient_data_dir = file.path(
    project_root, "IMC", "SpaceX", "BC", "expression_matrix"
  ),
  cpp_path = file.path(
    project_root,
    "SEA-AD", "BSNMani", "SEA-AD", "BSNMani-dev",
    "hybrid_M0_MALA_LR_FAST_v2.cpp"
  ),
  utils_path = file.path(
    project_root, "IMC", "Preprocessed", "Utils.R"
  ),
  output_dir = file.path(
    project_root, "IMC", "bc", analysis_output, "gene_importance"
  )
)

# Ordered metal tags used by the 29-marker IMC analysis.
selected_metals <- c(
  "Lu175", "Er167", "Eu151", "Nd143", "Dy163", "Nd144", "Gd158",
  "Gd156", "In113", "Er166", "La139", "Sm152", "Dy164", "Er168",
  "Yb172", "Sm147", "Gd155", "Pr141", "Tm169", "Yb174", "Nd145",
  "Nd150", "Yb176", "Dy162", "Gd160", "Nd148", "Nd146", "Sm149",
  "Nd142"
)

# All 44 metal tags represented in the posterior loading matrix.
metals_all <- c(
  "Ar80", "Dy162", "Dy163", "Dy164", "Er166", "Er167", "Er168",
  "Eu151", "Gd155", "Gd156", "Gd158", "Gd160", "Hg202", "In113",
  "In115", "Ir191", "Ir193", "La139", "Lu175", "Nd142", "Nd143",
  "Nd144", "Nd145", "Nd146", "Nd148", "Nd150", "Pb204", "Pb206",
  "Pr141", "Ru100", "Ru101", "Ru102", "Ru104", "Ru96", "Ru98",
  "Ru99", "Sm147", "Sm149", "Sm152", "Tm169", "Xe134", "Yb172",
  "Yb174", "Yb176"
)

# -----------------------------------------------------------------------------
# 2. Dependencies and input validation
# -----------------------------------------------------------------------------

required_packages <- c(
  "Rcpp", "ggplot2", "patchwork", "survival", "survminer"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0L) {
  stop(
    "Missing required R package(s): ",
    paste(missing_packages, collapse = ", "),
    ".\nInstall them before running this script.",
    call. = FALSE
  )
}

required_files <- c(
  config$cpp_path,
  config$utils_path,
  config$marker_path,
  config$clinical_path,
  file.path(
    config$model_run_dir, "train", "g1", "diagnostics", "g1_res_GQN_train.RDS"
  )
)
missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0L) {
  stop(
    "Required input file(s) not found:\n- ",
    paste(missing_files, collapse = "\n- "),
    "\n\nCheck the repository root supplied to this script.",
    call. = FALSE
  )
}

if (!dir.exists(config$patient_data_dir)) {
  stop(
    "Patient expression directory not found: ", config$patient_data_dir,
    call. = FALSE
  )
}

fig_dir <- file.path(config$output_dir, "fig")
km_dir <- file.path(config$output_dir, "km_curves")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(km_dir, recursive = TRUE, showWarnings = FALSE)

message("Repository root: ", config$project_root)
message("Output directory: ", config$output_dir)

# Compile/source project functions only after validating the paths.
Rcpp::sourceCpp(config$cpp_path)
source(config$utils_path)

if (!exists("polar_expansion", mode = "function")) {
  stop(
    "Function 'polar_expansion' was not found after sourcing project code.",
    call. = FALSE
  )
}

# -----------------------------------------------------------------------------
# 3. Helper functions
# -----------------------------------------------------------------------------

load_marker_map <- function(marker_path) {
  marker_data <- read.csv(marker_path, stringsAsFactors = FALSE)

  required_cols <- c("Original.Metal", "Gene")
  missing_cols <- setdiff(required_cols, names(marker_data))
  if (length(missing_cols) > 0L) {
    stop(
      "Marker mapping file is missing column(s): ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  marker_data
}

calculate_u_mean <- function(
    model_run_dir,
    marker_data,
    q,
    n_roi,
    metals_all,
    selected_metals) {

  posterior_path <- file.path(
    model_run_dir, "train", "g1", "diagnostics", "g1_res_GQN_train.RDS"
  )
  posterior <- readRDS(posterior_path)

  if (is.null(posterior$X)) {
    stop("Posterior object does not contain an 'X' element.", call. = FALSE)
  }

  expected_length <- n_roi * q
  if (length(posterior$X) != expected_length) {
    stop(
      "Unexpected posterior X length: ", length(posterior$X),
      ". Expected ", expected_length, " (= n_roi * q).",
      call. = FALSE
    )
  }

  x_mean <- matrix(posterior$X, nrow = n_roi, ncol = q)
  u_mean <- polar_expansion(x_mean)

  if (nrow(u_mean) != length(metals_all)) {
    stop(
      "Number of rows in U (", nrow(u_mean),
      ") does not match the number of metal tags (", length(metals_all), ").",
      call. = FALSE
    )
  }

  rownames(u_mean) <- metals_all

  missing_metals <- setdiff(selected_metals, rownames(u_mean))
  if (length(missing_metals) > 0L) {
    stop(
      "Selected metal tag(s) absent from U: ",
      paste(missing_metals, collapse = ", "),
      call. = FALSE
    )
  }

  u_mean <- u_mean[selected_metals, , drop = FALSE]

  gene_names <- marker_data$Gene[
    match(rownames(u_mean), marker_data$Original.Metal)
  ]

  if (anyNA(gene_names) || any(gene_names == "")) {
    unmapped <- rownames(u_mean)[is.na(gene_names) | gene_names == ""]
    stop(
      "Could not map the following metal tag(s) to genes: ",
      paste(unmapped, collapse = ", "),
      call. = FALSE
    )
  }

  rownames(u_mean) <- gene_names
  u_mean
}

calculate_hub_tables <- function(u_mean, top_n = 10L) {
  hub_tables <- vector("list", ncol(u_mean))

  for (i in seq_len(ncol(u_mean))) {
    u <- u_mean[, i, drop = FALSE]
    adjacency <- u %*% t(u)
    diag(adjacency) <- 0

    hub_df <- data.frame(
      Subnetwork = i,
      Gene = rownames(adjacency),
      HubScore = rowSums(abs(adjacency)),
      stringsAsFactors = FALSE
    )
    hub_df <- hub_df[order(hub_df$HubScore, decreasing = TRUE), , drop = FALSE]
    hub_tables[[i]] <- head(hub_df, top_n)
  }

  hub_tables
}

plot_hub_markers <- function(hub_tables, output_path) {
  plots <- lapply(seq_along(hub_tables), function(i) {
    hub_df <- hub_tables[[i]]

    ggplot2::ggplot(
      hub_df,
      ggplot2::aes(
        x = stats::reorder(Gene, HubScore),
        y = HubScore,
        fill = HubScore
      )
    ) +
      ggplot2::geom_col(color = "black", alpha = 0.8) +
      ggplot2::coord_flip() +
      ggplot2::scale_fill_gradient(low = "#377EB8", high = "#E41A1C") +
      ggplot2::theme_minimal(base_size = 16) +
      ggplot2::labs(
        title = paste0(
          "Top ", nrow(hub_df), " Hub Markers (Subnetwork ", i, ")"
        ),
        x = "Marker / Gene",
        y = "Total Connectivity (Hub Score)"
      ) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
        legend.position = "none"
      )
  })

  combined_plot <- patchwork::wrap_plots(plots, nrow = 1)

  ggplot2::ggsave(
    filename = output_path,
    plot = combined_plot,
    width = max(8, 7 * length(plots)),
    height = 6,
    dpi = 300
  )

  invisible(combined_plot)
}

resolve_target_markers <- function(marker_data, target_genes) {
  gene_idx <- match(target_genes, marker_data$Gene)

  if (anyNA(gene_idx)) {
    missing_genes <- target_genes[is.na(gene_idx)]
    warning(
      "No metal mapping found for target gene(s): ",
      paste(missing_genes, collapse = ", "),
      call. = FALSE
    )
  }

  mapped <- marker_data[gene_idx[!is.na(gene_idx)], c("Original.Metal", "Gene"), drop = FALSE]
  rownames(mapped) <- NULL

  if (anyDuplicated(mapped$Gene)) {
    stop("Target gene mapping contains duplicated gene names.", call. = FALSE)
  }

  mapped
}

aggregate_patient_expression <- function(
    patient_ids,
    patient_data_dir,
    target_markers) {

  patient_ids <- unique(as.character(patient_ids))
  genes <- target_markers$Gene
  metals <- target_markers$Original.Metal

  expression_matrix <- matrix(
    NA_real_,
    nrow = length(patient_ids),
    ncol = length(genes),
    dimnames = list(patient_ids, genes)
  )

  missing_files <- character(0)
  missing_columns <- list()

  for (i in seq_along(patient_ids)) {
    core_id <- patient_ids[[i]]
    file_path <- file.path(patient_data_dir, paste0(core_id, ".csv"))

    if (!file.exists(file_path)) {
      missing_files <- c(missing_files, core_id)
      next
    }

    single_cell <- read.csv(
      file_path,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )

    available <- metals %in% names(single_cell)

    if (any(available)) {
      means <- vapply(
        metals[available],
        function(metal) mean(single_cell[[metal]], na.rm = TRUE),
        numeric(1)
      )
      means[is.nan(means)] <- NA_real_
      expression_matrix[i, genes[available]] <- means
    }

    if (any(!available)) {
      missing_columns[[core_id]] <- metals[!available]
    }
  }

  if (length(missing_files) > 0L) {
    warning(
      length(missing_files),
      " patient expression file(s) were not found. Example(s): ",
      paste(utils::head(missing_files, 5L), collapse = ", "),
      call. = FALSE
    )
  }

  if (length(missing_columns) > 0L) {
    n_affected <- length(missing_columns)
    warning(
      "One or more target metal columns were absent for ", n_affected,
      " patient file(s). Missing values were recorded as NA.",
      call. = FALSE
    )
  }

  data.frame(
    core = patient_ids,
    expression_matrix,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

save_survival_plot <- function(surv_plot, output_path) {
  grDevices::png(
    filename = output_path,
    width = 2100,
    height = 2100,
    res = 300
  )
  on.exit(grDevices::dev.off(), add = TRUE)
  print(surv_plot)
  invisible(NULL)
}

plot_km_for_gene <- function(data, gene, output_dir) {
  required_cols <- c("OSmonth", "Patientstatus", gene)
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0L) {
    warning(
      "Skipping ", gene, ": missing column(s): ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
    return(NULL)
  }

  complete <- stats::complete.cases(
    data[, c("OSmonth", "Patientstatus", gene), drop = FALSE]
  )
  plot_df <- data[complete, , drop = FALSE]

  if (nrow(plot_df) < 2L) {
    warning("Skipping ", gene, ": insufficient complete observations.", call. = FALSE)
    return(NULL)
  }

  median_expr <- stats::median(plot_df[[gene]], na.rm = TRUE)
  plot_df$Marker_Group <- ifelse(plot_df[[gene]] > median_expr, "High", "Low")
  plot_df$Marker_Group <- factor(plot_df$Marker_Group, levels = c("Low", "High"))

  if (length(unique(stats::na.omit(plot_df$Marker_Group))) < 2L) {
    warning(
      "Skipping ", gene,
      ": median split did not produce both High and Low groups.",
      call. = FALSE
    )
    return(NULL)
  }

  fit <- survival::survfit(
    survival::Surv(OSmonth, Patientstatus) ~ Marker_Group,
    data = plot_df
  )

  km_plot <- survminer::ggsurvplot(
    fit,
    data = plot_df,
    pval = TRUE,
    conf.int = TRUE,
    risk.table = TRUE,
    palette = c("#377EB8", "#E41A1C"),
    title = paste("Patient Survival by Mean", gene, "Expression"),
    xlab = "Time (Months)",
    ylab = "Survival Probability",
    legend.title = paste(gene, "Level"),
    ggtheme = ggplot2::theme_minimal(base_size = 14)
  )

  output_path <- file.path(output_dir, paste0("KM_Curve_", gene, ".png"))
  save_survival_plot(km_plot, output_path)

  invisible(km_plot)
}

# -----------------------------------------------------------------------------
# 4. Hub-marker analysis
# -----------------------------------------------------------------------------

marker_data <- load_marker_map(config$marker_path)

u_mean <- calculate_u_mean(
  model_run_dir = config$model_run_dir,
  marker_data = marker_data,
  q = config$q,
  n_roi = config$n_roi,
  metals_all = metals_all,
  selected_metals = selected_metals
)

hub_tables <- calculate_hub_tables(u_mean, top_n = config$top_n)
hub_table_all <- do.call(rbind, hub_tables)
rownames(hub_table_all) <- NULL

utils::write.csv(
  hub_table_all,
  file.path(config$output_dir, "Top_Hub_Markers.csv"),
  row.names = FALSE
)

plot_hub_markers(
  hub_tables,
  file.path(fig_dir, "Top_Hub_Markers.png")
)

message("Saved hub-marker results.")

target_genes <- if (is.null(config$target_genes)) {
  unique(hub_table_all$Gene)
} else {
  unique(config$target_genes)
}
message("Genes selected for KM analysis: ", paste(target_genes, collapse = ", "))

# -----------------------------------------------------------------------------
# 5. Patient-level marker expression and Kaplan-Meier analysis
# -----------------------------------------------------------------------------

clinical_df <- read.csv(
  config$clinical_path,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_clinical_cols <- c("core", "OSmonth", "Patientstatus")
missing_clinical_cols <- setdiff(required_clinical_cols, names(clinical_df))
if (length(missing_clinical_cols) > 0L) {
  stop(
    "Clinical file is missing required column(s): ",
    paste(missing_clinical_cols, collapse = ", "),
    call. = FALSE
  )
}

mapped_targets <- resolve_target_markers(marker_data, target_genes)

if (nrow(mapped_targets) == 0L) {
  stop("None of the requested target genes could be mapped to metal tags.", call. = FALSE)
}

message(
  "Mapped ", nrow(mapped_targets), " of ", length(target_genes),
  " requested target genes."
)

patient_expression <- aggregate_patient_expression(
  patient_ids = clinical_df$core,
  patient_data_dir = config$patient_data_dir,
  target_markers = mapped_targets
)

utils::write.csv(
  patient_expression,
  file.path(config$output_dir, "Patient_Level_Marker_Expression.csv"),
  row.names = FALSE
)

merged_data <- merge(
  clinical_df,
  patient_expression,
  by = "core",
  all = FALSE,
  sort = FALSE
)

message("Generating Kaplan-Meier curves for mapped target genes...")
for (gene in mapped_targets$Gene) {
  message("  - ", gene)
  plot_km_for_gene(merged_data, gene, km_dir)
}

message("Analysis complete.")
