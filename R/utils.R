.integer_scalar <- function(x, name, lower = 1, upper = Inf) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      x != floor(x) || x < lower || x > upper)
    stop(name, " must be an integer in [", lower, ", ", upper, "].", call. = FALSE)
  as.integer(x)
}

.probability <- function(x, name, closed_right = FALSE) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      x <= 0 || x > 1 || (!closed_right && x == 1))
    stop(name, " must be in (0, ", if (closed_right) "1]." else "1).", call. = FALSE)
  x
}

.check_data <- function(z, y, allow_infinite = FALSE) {
  if (!(is.numeric(z) || is.logical(z)) || anyNA(z) ||
      !all(z %in% c(0, 1)) || length(unique(z)) != 2L)
    stop("z must contain both 0 and 1, with no missing values.", call. = FALSE)
  if (!is.numeric(y) || length(y) != length(z) || anyNA(y) ||
      (!allow_infinite && any(!is.finite(y))))
    stop("y must be a numeric vector of the same length as z, without missing",
         if (allow_infinite) " values." else " or infinite values.", call. = FALSE)
  invisible(NULL)
}

.tie_order <- function(tie_order, n) {
  # A fresh random priority is safe even when rows were sorted by assignment.
  # Resolve once per analysis and retain it throughout imputation/inversion.
  if (is.null(tie_order)) return(sample.int(n))
  if (!is.numeric(tie_order) || length(tie_order) != n ||
      any(!is.finite(tie_order)) || anyDuplicated(tie_order))
    stop("tie_order must have one distinct finite number per unit.", call. = FALSE)
  tie_order
}

.adjust_y <- function(y, adjustment) {
  if (is.null(adjustment)) return(y)
  if (!is.numeric(adjustment) || length(adjustment) != length(y) ||
      anyNA(adjustment) || any(adjustment == -Inf))
    stop("adjustment must have one finite or +Inf value per unit.", call. = FALSE)
  y - adjustment
}

.with_seed <- function(seed, expr) {
  if (is.null(seed)) return(eval.parent(substitute(expr)))
  seed <- .integer_scalar(seed, "seed", 0, .Machine$integer.max)
  existed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (existed) old <- get(".Random.seed", envir = .GlobalEnv)
  on.exit(if (existed) assign(".Random.seed", old, envir = .GlobalEnv)
          else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
            rm(".Random.seed", envir = .GlobalEnv))
  set.seed(seed)
  eval.parent(substitute(expr))
}

.statistic_name <- function(statistic) {
  match.arg(statistic, c("rks", "wilcoxon", "stephenson"))
}

.scores <- function(n, statistic, s) {
  if (statistic == "wilcoxon") return(seq_len(n))
  s <- .integer_scalar(s, "s", 2, n)
  # A positive common scale leaves all tests unchanged and prevents overflow.
  exp(lchoose(seq_len(n) - 1, s - 1) - lchoose(n - 1, s - 1))
}

.snap_integer <- function(x) {
  near <- is.finite(x) & abs(x - round(x)) <= 32 * .Machine$double.eps * pmax(1, abs(x))
  x[near] <- round(x[near])
  x
}

.midpoint <- function(lo, hi) {
  width <- hi-lo
  mid <- lo+width/2
  overflow <- !is.finite(width)
  mid[overflow] <- lo[overflow]/2+hi[overflow]/2
  mid
}
