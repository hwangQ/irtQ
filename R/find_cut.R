# find_cut.R
# ---------------------------------------------------------------------------
# find_cut() — Identify TIF-crossing cut scores for MST routing
#
# For each pair of adjacent modules (sorted by mean item difficulty) at each
# stage transition, this function finds the theta value where the harder
# module's TIF curve crosses the easier module's TIF curve with a positive
# slope (a "proper" crossing).  The resulting cut scores can be passed
# directly to run_mst() as the cut_score argument.
# ---------------------------------------------------------------------------


#' Find TIF-Crossing Cut Scores for MST Routing
#'
#' @description
#' Computes routing cut scores for a multistage test (MST) by identifying the
#' theta values at which adjacent modules' Test Information Functions (TIFs)
#' cross in the correct direction. The resulting cut scores can be passed
#' directly to \code{\link{run_mst}} as its \code{cut_score} argument
#' (with \code{route_method = NULL}), providing a principled and psychometrically
#' grounded alternative to heuristic cut-score selection.
#'
#' @param x A data frame of item metadata in the standard \pkg{irtQ} format
#'   (columns: \code{id}, \code{cats}, \code{model}, \code{par.1},
#'   \code{par.2}, \ldots). See \code{\link{shape_df}} for details.
#' @param module A binary integer matrix of dimensions \emph{J} \eqn{\times}
#'   \emph{M}, where \emph{J} is the number of items and \emph{M} is the
#'   number of modules. A value of 1 at row \emph{j}, column \emph{m}
#'   indicates that item \emph{j} belongs to module \emph{m}. This is the
#'   same \code{module} argument used in \code{\link{run_mst}}.
#' @param route_map A binary square matrix of dimensions \emph{M} \eqn{\times}
#'   \emph{M} defining the MST transition structure. A value of 1 at row
#'   \emph{i}, column \emph{j} indicates that examinees can be routed from
#'   module \emph{i} to module \emph{j}. This is the same \code{route_map}
#'   argument used in \code{\link{run_mst}}.
#' @param D A numeric scaling constant for the IRT model. Default is
#'   \code{1.702}.
#' @param theta_range A numeric vector of length 2 specifying the theta range
#'   over which to search for TIF crossings. Default is \code{c(-6, 6)}.
#' @param n_grid An integer specifying the number of equally spaced theta
#'   points used for the initial sign-change scan. A finer grid reduces the
#'   chance of missing a crossing but increases computation time. Default is
#'   \code{2001}.
#' @param ref_theta A single numeric value used to break ties when multiple
#'   proper crossings are found for an adjacent module pair. The crossing
#'   closest to \code{ref_theta} is selected. Default is \code{0}.
#'
#' @return An object of class \code{"find_cut"}, which is a named list with
#'   three elements:
#' \describe{
#'   \item{\code{cut_score}}{A list of length \emph{S} - 1 (where \emph{S}
#'     is the number of stages), compatible with the \code{cut_score} argument
#'     of \code{\link{run_mst}}. Each element is a numeric vector of cut
#'     scores for the transition into that stage, sorted in ascending order.
#'     Stages with only one module have an empty numeric vector
#'     (\code{numeric(0)}).}
#'   \item{\code{details}}{A named list of per-stage diagnostic information,
#'     including module indices, mean item locations, difficulty-sorted module
#'     order, and per-pair crossing details (proper crossings, anomalous
#'     crossings, and the selected cut score).}
#'   \item{\code{tif_data}}{A \code{tibble} with columns \code{theta},
#'     \code{module}, \code{stage}, and \code{tif}, providing TIF curves for
#'     all modules at stages 2 and above. Useful for visualizing crossings with
#'     \pkg{ggplot2}.}
#' }
#'
#' @details
#' \subsection{Motivation: the MFI routing problem}{
#' When \code{\link{run_mst}} uses Maximum Fisher Information (MFI) routing
#' (\code{route_method = "mfi"}), it selects the next-stage module that
#' maximises TIF at the examinee's current ability estimate. In principle, a
#' harder module should always have a higher TIF for high-ability examinees,
#' and a lower TIF for low-ability examinees, making MFI routing work as
#' intended. However, in practice, a harder module's TIF curve can peak at
#' a theta value that is unexpectedly low, so that the harder module
#' \emph{also} has the highest TIF near the low end of the ability scale.
#' In this case, a low-ability examinee is incorrectly routed into the harder
#' module — the opposite of what the MST is designed to do. This is called an
#' \emph{anomalous routing} or \emph{path reversal}.
#' }
#'
#' \subsection{Solution: TIF-crossing cut scores}{
#' \code{find_cut()} resolves this problem by computing cut scores from the
#' TIF curves themselves, rather than relying on pure MFI at the time of
#' routing. Specifically, for each pair of adjacent modules (ordered by mean
#' item difficulty), the function finds the theta value at which the harder
#' module's TIF first overtakes the easier module's TIF \emph{from below}
#' as theta increases. This crossing point is the natural boundary between
#' the two modules: below it, the easier module is more informative; above it,
#' the harder module is more informative. Using this point as a fixed cut
#' score pre-empts the path reversal problem entirely.
#' }
#'
#' \subsection{Proper vs. anomalous crossings}{
#' Two TIF curves may cross more than once. \code{find_cut()} distinguishes
#' two types of crossing:
#' \itemize{
#'   \item \strong{Proper crossing}: the difference
#'     \eqn{\text{TIF}_{\text{harder}}(\theta) - \text{TIF}_{\text{easier}}(\theta)}
#'     changes sign from negative to positive (positive slope at the crossing).
#'     This is the correct, expected crossing in the middle of the theta scale.
#'   \item \strong{Anomalous crossing}: the difference changes sign from
#'     positive to negative (negative slope at the crossing). This corresponds
#'     to the pathological situation in which the harder module temporarily
#'     dominates at the low end of the theta scale — exactly the pattern that
#'     causes path reversals under MFI routing.
#' }
#' Only proper crossings are used as cut scores. Anomalous crossings are
#' reported in \code{$details} for diagnostic purposes but excluded from the
#' returned \code{cut_score}. If no proper crossing is found for a module pair,
#' an informative error is thrown. If multiple proper crossings are found
#' (unusual for well-designed modules), a warning is issued and the one
#' closest to \code{ref_theta} is selected.
#' }
#'
#' \subsection{Algorithm}{
#' The function proceeds stage by stage (from stage 2 onward). At each stage:
#' \enumerate{
#'   \item Modules are sorted by their mean item location parameter (ascending
#'     difficulty), using the internal \code{mean_loc()} helper also used by
#'     the \code{"bmat"} routing method in \code{\link{run_mst}}.
#'   \item For each adjacent pair (easier, harder) in difficulty order, the
#'     difference \eqn{\text{TIF}_{\text{harder}}(\theta) -
#'     \text{TIF}_{\text{easier}}(\theta)} is evaluated on a fine theta grid
#'     of \code{n_grid} points. Sign changes on the grid are detected and
#'     then refined with \code{\link[stats]{uniroot}}.
#'   \item The slope at each refined root is estimated by a finite difference
#'     (step size 0.01) to classify it as proper or anomalous.
#'   \item The selected proper crossing is stored as the cut score for that pair.
#' }
#' }
#'
#' \subsection{Module index vs. difficulty order}{
#' \code{\link{run_mst}} internally orders candidate modules by ascending
#' module index and maps routing rank 1 to the lowest-index module. For
#' cut-score routing to assign low-theta examinees to the easiest module, the
#' modules at each stage must be numbered in ascending difficulty order
#' (easiest module = lowest index). \code{find_cut()} checks this assumption
#' and issues a warning if it is violated. The returned cut scores are always
#' sorted in ascending order, as required by \code{\link[base]{cut}}.
#' }
#'
#' @seealso \code{\link{run_mst}}, \code{\link{panel_info}},
#'   \code{\link{reval_mst}}, \code{\link{info}}
#'
#' @examples
#' ## ── Setup: use the built-in simMST 1-3-3 panel ──────────────────────────
#' ## simMST is a 7-module, 3-stage MST panel (8 dichotomous 3PLM items each).
#' ## Modules 1 (routing), 2-4 (stage 2: easy/medium/hard),
#' ##         5-7 (stage 3: easy/medium/hard).
#' x         <- simMST$item_bank
#' module    <- simMST$module
#' route_map <- simMST$route_map
#'
#' ## ── Find TIF-crossing cut scores ─────────────────────────────────────────
#' ## For each adjacent module pair at stages 2 and 3, find_cut() identifies
#' ## the theta at which the harder module's TIF first exceeds the easier
#' ## module's TIF (proper crossing), and returns it as a cut score.
#' cut_result <- find_cut(x = x, module = module, route_map = route_map)
#'
#' ## Print a summary: crossing points found, anomalous crossings excluded,
#' ## and the final selected cut scores per stage transition.
#' print(cut_result)
#'
#' ## Inspect the cut_score element — a list directly compatible with run_mst()
#' ## cut_score[[1]]: two cut scores for the stage-1 -> stage-2 transition
#' ## cut_score[[2]]: cut scores for the stage-2 -> stage-3 transition
#' cut_result$cut_score
#'
#' ## Compare with the manually specified cut scores stored in simMST
#' simMST$cut_score
#'
#' ## ── Use the cut scores in run_mst() ──────────────────────────────────────
#' ## Pass cut_result$cut_score directly to run_mst() with route_method = NULL.
#' ## This replaces pure MFI routing with TIF-crossing-based fixed cut scores,
#' ## preventing anomalous path reversals at the extremes of the theta scale.
#' \dontrun{
#' set.seed(1)
#' theta_true <- rnorm(500)
#' result <- run_mst(
#'   x            = x,
#'   route_map    = route_map,
#'   module       = module,
#'   theta        = theta_true,
#'   route_method = NULL,
#'   cut_score    = cut_result$cut_score
#' )
#' }
#'
#' @export
find_cut <- function(x,
                     module,
                     route_map,
                     D           = 1.702,
                     theta_range = c(-6, 6),
                     n_grid      = 2001L,
                     ref_theta   = 0) {

  # ── 0. Input validation ────────────────────────────────────────────────────
  if (!is.data.frame(x))
    stop("'x' must be a data frame of item metadata in irtQ format.")
  if (!is.matrix(module) || !is.numeric(module))
    stop("'module' must be a numeric binary matrix (J x M).")
  if (nrow(x) != nrow(module))
    stop("Number of rows in 'x' must equal number of rows in 'module' (J items).")
  # route_map may be a matrix or a data.frame (panel_info() accepts both)
  if (!(is.matrix(route_map) || is.data.frame(route_map)) ||
      !all(vapply(as.data.frame(route_map), is.numeric, logical(1L))))
    stop("'route_map' must be a numeric binary square matrix (M x M).")
  if (ncol(module) != nrow(route_map))
    stop("Number of columns in 'module' (M modules) must equal dimensions of 'route_map'.")
  if (!is.numeric(D) || length(D) != 1L || !is.finite(D))
    stop("'D' must be a single finite numeric scalar.")
  if (!is.numeric(theta_range) || length(theta_range) != 2L ||
      theta_range[1L] >= theta_range[2L])
    stop("'theta_range' must be a numeric vector of length 2 with theta_range[1] < theta_range[2].")
  n_grid <- as.integer(n_grid)                    # coerce to integer
  if (n_grid < 100L)
    stop("'n_grid' must be an integer >= 100.")
  if (!is.numeric(ref_theta) || length(ref_theta) != 1L || !is.finite(ref_theta))
    stop("'ref_theta' must be a single finite numeric scalar.")

  # ── 1. Panel structure ─────────────────────────────────────────────────────
  # Use panel_info() to identify which modules belong to each stage
  pinfo  <- panel_info(route_map)
  config <- pinfo$config    # named list: config[[s]] = integer vector of module indices
  n_stg  <- pinfo$n.stage   # total number of stages (integer)
  n_mods <- ncol(module)    # total number of modules M

  # Build a list of item metadata subsets, one per module
  # item_mod[[m]]: rows of x assigned to module m (module[,m] == 1)
  item_mod <- lapply(seq_len(n_mods), function(m) {
    x[module[, m] == 1L, , drop = FALSE]
  })

  # Fine theta grid for sign-change scan (root detection)
  theta_grid <- seq(theta_range[1L], theta_range[2L], length.out = n_grid)

  # ── 2. Allocate output containers ──────────────────────────────────────────
  # cut_list[[k]]: cut score vector for transition from stage k to stage k+1
  # Indexed 1 … n_stg-1, matching the cut_score argument of run_mst()
  cut_list     <- vector("list", n_stg - 1L)
  names(cut_list) <- paste0("stage.", seq(2L, n_stg))

  # details[[k]]: per-stage diagnostic information
  details      <- vector("list", n_stg - 1L)
  names(details) <- names(cut_list)

  # tif_rows: accumulator for tif_data tibble rows
  tif_rows     <- list()

  # ── 3. Main loop: process each stage from stage 2 onward ──────────────────
  for (s in seq(2L, n_stg)) {

    mod_ids <- config[[s]]    # module indices at stage s (integer vector)
    n_mod_s <- length(mod_ids) # number of modules at this stage

    # ── 3a. Single module: no routing decision at this stage ─────────────────
    if (n_mod_s == 1L) {
      cut_list[[s - 1L]] <- numeric(0L)   # empty: no cut score needed
      details[[s - 1L]]  <- list(
        modules      = mod_ids,
        mean_locs    = NULL,
        sorted_order = NULL,
        pairs        = NULL,
        message      = "Single module at this stage; no cut score required."
      )
      # Still collect TIF data for the single module (for tif_data tibble)
      tif_vals <- info(x = item_mod[[mod_ids]], theta = theta_grid,
                       D = D, tif = TRUE)$tif
      tif_rows[[length(tif_rows) + 1L]] <- tibble::tibble(
        theta  = theta_grid,
        module = mod_ids,
        stage  = s,
        tif    = tif_vals
      )
      next   # move to next stage
    }

    # ── 3b. Precompute TIF matrix for all modules at stage s ─────────────────
    # tif_mat: n_grid × n_mod_s matrix; column k = TIF of mod_ids[k]
    tif_mat <- vapply(
      X         = mod_ids,
      FUN       = function(m) info(x = item_mod[[m]], theta = theta_grid,
                                   D = D, tif = TRUE)$tif,
      FUN.VALUE = numeric(n_grid)
    )

    # Accumulate TIF rows for visualization tibble
    for (ki in seq_len(n_mod_s)) {
      tif_rows[[length(tif_rows) + 1L]] <- tibble::tibble(
        theta  = theta_grid,
        module = mod_ids[ki],   # module index
        stage  = s,
        tif    = tif_mat[, ki]  # TIF values for this module
      )
    }

    # ── 3c. Sort modules by ascending mean item location (difficulty) ─────────
    # mean_loc() returns a vector of per-item mean location parameters
    mean_locs <- vapply(
      X         = mod_ids,
      FUN       = function(m) mean(mean_loc(item_mod[[m]])),
      FUN.VALUE = numeric(1L)
    )

    diff_ord    <- order(mean_locs)       # ascending difficulty → index in mod_ids
    sorted_mods <- mod_ids[diff_ord]      # module ids sorted by ascending difficulty
    sorted_locs <- mean_locs[diff_ord]    # corresponding mean locations

    # ── 3d. Check: index order should match difficulty order ──────────────────
    # run_mst() maps routing rank 1 → next_possible[1] (lowest module index).
    # For cut-score routing to assign low theta → easy module, the modules at
    # each stage must have ascending indices in ascending difficulty order.
    if (!identical(sorted_mods, sort(mod_ids))) {
      warning(sprintf(paste0(
        "Stage %d: difficulty order of modules does not match their index order.\n",
        "  Ascending index order    : modules %s\n",
        "  Ascending difficulty order: modules %s (mean locations: %s)\n",
        "  run_mst() routes by ascending module index (rank 1 -> lowest index).\n",
        "  If index order differs from difficulty order, cut-score routing may\n",
        "  not produce the intended module assignments.\n",
        "  Consider redesigning the panel so that easier modules have lower indices."
      ),
      s,
      paste(sort(mod_ids),         collapse = ", "),
      paste(sorted_mods,           collapse = ", "),
      paste(round(sorted_locs, 3L), collapse = ", ")
      ))
    }

    # ── 3e. Find cut score for each adjacent difficulty-sorted pair ───────────
    n_pairs   <- n_mod_s - 1L                  # number of adjacent pairs
    pair_cuts <- numeric(n_pairs)              # one cut score per pair
    pair_info <- vector("list", n_pairs)       # diagnostic info per pair

    for (k in seq_len(n_pairs)) {

      easier_m <- sorted_mods[k]          # easier module index
      harder_m <- sorted_mods[k + 1L]    # harder module index

      # Column positions of easier and harder modules within tif_mat
      col_e    <- which(mod_ids == easier_m)
      col_h    <- which(mod_ids == harder_m)

      # Difference of TIF values on the grid: positive means harder > easier
      diff_grid <- tif_mat[, col_h] - tif_mat[, col_e]

      # Inline evaluation function for uniroot() and derivative estimation
      # Captures easier_m and harder_m from the current loop iteration
      e_m <- easier_m   # local copy to avoid closure pitfalls
      h_m <- harder_m
      diff_tif_fn <- function(theta_val) {
        # Evaluate TIF difference at an arbitrary theta (scalar or vector)
        tif_h <- info(x = item_mod[[h_m]], theta = theta_val, D = D, tif = TRUE)$tif
        tif_e <- info(x = item_mod[[e_m]], theta = theta_val, D = D, tif = TRUE)$tif
        as.numeric(tif_h - tif_e)
      }

      # Detect sign-change positions on the fine grid
      # diff() of sign() gives +2, -2, or 0; nonzero means a sign change occurred
      chng_idx <- which(diff(sign(diff_grid)) != 0L)

      # Error: no sign changes → no crossing in the search range
      if (length(chng_idx) == 0L) {
        stop(sprintf(paste0(
          "Stage %d, pair (module %d [easier] vs module %d [harder]): ",
          "no TIF crossing detected in theta range [%.2f, %.2f].\n",
          "The TIF curves do not intersect within this range.\n",
          "Consider widening 'theta_range' or reviewing module item composition."
        ), s, easier_m, harder_m, theta_range[1L], theta_range[2L]))
      }

      # Refine each sign-change interval with uniroot() and classify the root
      all_roots  <- numeric(0L)      # all refined roots (proper + anomalous)
      root_types <- character(0L)    # "proper" or "anomalous" for each root

      for (ci in chng_idx) {
        lo <- theta_grid[ci]        # left endpoint of sign-change interval
        hi <- theta_grid[ci + 1L]  # right endpoint

        # Re-evaluate at endpoints (guard against numerical ties on grid)
        val_lo <- diff_tif_fn(lo)
        val_hi <- diff_tif_fn(hi)
        if (sign(val_lo) == sign(val_hi)) next   # not a genuine crossing

        # Refine the root position with uniroot()
        root_i <- tryCatch(
          stats::uniroot(
            f        = diff_tif_fn,
            interval = c(lo, hi),
            tol      = .Machine$double.eps^0.5   # tight tolerance
          )$root,
          error = function(e) NA_real_
        )
        if (is.na(root_i)) next    # uniroot failed — skip this candidate

        # Classify the crossing by the sign of the slope (finite difference, h = 0.01)
        slope_i <- (diff_tif_fn(root_i + 0.01) - diff_tif_fn(root_i - 0.01)) / 0.02
        type_i  <- if (slope_i > 0) "proper" else "anomalous"

        all_roots  <- c(all_roots,  root_i)
        root_types <- c(root_types, type_i)
      }

      proper_roots    <- all_roots[root_types == "proper"]    # valid cut candidates
      anomalous_roots <- all_roots[root_types == "anomalous"] # pathological crossings

      # ── Error: no proper crossings ──────────────────────────────────────────
      if (length(proper_roots) == 0L) {
        if (length(anomalous_roots) > 0L) {
          # Crossings exist but ALL are anomalous (harder module dominates at low theta)
          stop(sprintf(paste0(
            "Stage %d, pair (module %d [easier] vs module %d [harder]): ",
            "no valid cut score found.\n",
            "All %d crossing(s) detected have a negative slope (anomalous pattern):\n",
            "  theta = [%s]\n",
            "This indicates that the harder module unexpectedly has higher TIF at\n",
            "low theta values, which would cause misrouting of low-ability examinees.\n",
            "Review module item composition or consult a test design specialist."
          ),
          s, easier_m, harder_m,
          length(anomalous_roots),
          paste(round(anomalous_roots, 4L), collapse = ", ")))
        } else {
          # Sign changes existed on grid but uniroot() found nothing — numerical issue
          stop(sprintf(paste0(
            "Stage %d, pair (module %d [easier] vs module %d [harder]): ",
            "sign changes were detected on the grid but no root could be refined.\n",
            "Try increasing 'n_grid' or widening 'theta_range'."
          ), s, easier_m, harder_m))
        }
      }

      # ── Warning: multiple proper crossings → select closest to ref_theta ───
      if (length(proper_roots) > 1L) {
        selected_root <- proper_roots[which.min(abs(proper_roots - ref_theta))]
        warning(sprintf(paste0(
          "Stage %d, pair (module %d [easier] vs module %d [harder]): ",
          "%d proper crossings found at theta = [%s].\n",
          "Selecting the crossing closest to ref_theta = %.4f: theta = %.4f.\n",
          "To override the selection, specify a different 'ref_theta'."
        ),
        s, easier_m, harder_m,
        length(proper_roots),
        paste(round(proper_roots, 4L), collapse = ", "),
        ref_theta, selected_root))
      } else {
        selected_root <- proper_roots[1L]    # exactly one proper crossing
      }

      pair_cuts[k] <- selected_root    # store selected cut score for this pair

      # Record per-pair diagnostic information
      pair_info[[k]] <- list(
        easier_module       = easier_m,
        harder_module       = harder_m,
        proper_crossings    = proper_roots,      # all valid crossing points
        anomalous_crossings = anomalous_roots,   # excluded pathological crossings
        selected_cut        = selected_root      # the cut score used
      )
      names(pair_info)[k] <- sprintf("mod%d_vs_mod%d", easier_m, harder_m)

    }  # end pair loop (k)

    # ── 3f. Store stage-level results ─────────────────────────────────────────
    # Cut scores must be sorted ascending: give_path() uses them as breaks in
    # cut(x, breaks = c(-Inf, cut_sc, Inf)), so ascending order is required.
    cut_list[[s - 1L]] <- sort(pair_cuts)

    # Collect diagnostic information for this stage
    details[[s - 1L]] <- list(
      modules      = mod_ids,                              # module indices at this stage
      mean_locs    = setNames(mean_locs, paste0("mod", mod_ids)),  # named mean locations
      sorted_order = sorted_mods,                         # difficulty-sorted module order
      pairs        = pair_info                            # per-pair crossing details
    )

  }  # end stage loop (s)

  # ── 4. Build tif_data tibble ───────────────────────────────────────────────
  # Tidy data frame with one row per (theta, module) combination
  # Suitable for ggplot2 visualisation of TIF curves and crossing points
  tif_data <- dplyr::bind_rows(tif_rows)

  # ── 5. Assemble and return the result object ───────────────────────────────
  rst <- list(
    cut_score = cut_list,    # list: pass directly to run_mst(cut_score = ...)
    details   = details,     # list: per-stage diagnostic info
    tif_data  = tif_data     # tibble: TIF curves for all stages >= 2
  )
  class(rst) <- "find_cut"   # S3 class for print dispatch
  rst
}
