#' @importFrom withr with_envvar
#' @importFrom withr with_options
#'
julia_spca_open_pipe <- function(
    mat, output_file, K, D, search_cap, eigen_gap, uint64_seed = NULL) {
  project_dir <- find.package("spcaFeatures")
  julia_depot <- file.path(project_dir, "julia_depot")
  mat_output <- tempfile("cov", fileext = ".csv.gz")
  csv_gz_conn <- file(mat_output, "wb")
  csv_conn <- gzcon(csv_gz_conn)
  write.csv(mat, csv_conn, row.names = F)
  close(csv_conn)

  with_options(
    list(scipen = 100),
    with_envvar(
      c(
        JULIA_DEPOT_PATH = file.path(project_dir, "julia_depot"),
        LD_LIBRARY_PATH = NA
      ),
      p <- pipe(
        paste0(
          "julia --project=\"",
          project_dir,
          "\" \"",
          file.path(project_dir, "julia_project", "SPCA.jl"),
          "\" \"",
          mat_output,
          "\" \"",
          output_file,
          "\" ",
          K,
          " ",
          D,
          " ",
          search_cap,
          " ",
          eigen_gap,
          ifelse(
            is.null(uint64_seed),
            "",
            paste0(" ", uint64_seed)
          )
        ),
        "rb",
        blocking = TRUE
      )
    )
  )
}

#' @importFrom magrittr %>%
#' @importFrom progress progress_bar
#'
run_optimal_spca <- function(
    mat, K, D, search_cap = 500000, eigen_gap = 0.001, uint64_seed = NULL) {
  output_file <- tempfile("SPCA", fileext = ".csv")
  p <- julia_spca_open_pipe(
    mat, output_file, K, D, search_cap, eigen_gap, uint64_seed
  )
  num_progress_dots <- 20

  # Assume that the total number of ticks is D * num_progress_dots. If one PC
  # is optimized quickly (newline), then remove ticks from total_ticks, but
  # don't jump ahead a full percentage point (value of D).
  num_ticks <- D * num_progress_dots

  message("Searching the covariance for optimal SPCA")
  prog <- progress_bar$new("SPCA [:bar] :percent :eta", force = T)

  total_progress <- 0
  pc_progress <- 0
  num_pcs <- 0
  while (length(chr <- readChar(p, 1)) == 1) {
    if (chr == "\n") {
      # We have already seen between 0 and num_progress_dots progress updates.
      # We don't want to show no progress at all (progress happens when we have
      # a dot indicating a certain number of search iterations). So if the
      # current PC exited extremely easily, then we will increment pc_progress.
      if (pc_progress == 0) {
        total_progress <- total_progress + 1
        pc_progress <- 1
      }
      # Remove the remaining dots yet to be seen from num_ticks.
      num_ticks <- num_ticks - (num_progress_dots - pc_progress)
      pc_progress <- 0
      # Tracked so that we can set prog$finished and avoid a progress bar error.
      num_pcs <- num_pcs + 1
    } else {
      # Update with a progress dot (searches performed within the current PC).
      # This updates both the total progress (work done) and pc progress
      # (work done on the current PC).
      total_progress <- total_progress + 1
      pc_progress <- pc_progress + 1
    }
    if (num_pcs == D) {
      stopifnot(total_progress == num_ticks)
      if (!prog$finished) prog$update(1)
    } else {
      prog$update(total_progress / num_ticks)
    }
  }
  close(p)
  feature.loadings <- read.csv(output_file, check.names = F) %>% as.matrix()
}
