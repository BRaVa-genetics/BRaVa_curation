library(data.table)
library(dplyr)

reconstruct_saige_counts <- function(beta, se, n_case, n_ctrl) {
  Cc <- 2*n_case; Ct <- 2*n_ctrl; OR <- exp(beta)
  f <- function(ac) {
    at <- ac * Ct / (OR*Cc + ac*(1-OR))
    if (ac<=0 || ac>=Cc || at<=0 || at>=Ct) return(NA_real_)
    (1/ac + 1/(Cc-ac) + 1/at + 1/(Ct-at)) - se^2
  }
  # bracket and root-find for ac
  xs <- seq(1e-6, Cc-1e-6, length.out=400)
  vals <- sapply(xs, f)
  idx <- which(diff(sign(vals))!=0)[1]
  if (is.na(idx)) stop("Could not bracket a root; check inputs.")
  ac <- uniroot(f, lower=xs[idx], upper=xs[idx+1])$root
  at <- ac * Ct / (OR*Cc + ac*(1-OR))
  list(
    AC_case = ac, AC_ctrl = at,
    AF_case = ac/Cc, AF_ctrl = at/Ct,
    AC = ac+at, AF = (ac+at)/(Cc+Ct)
  )
}

react_analytic_both <- function(beta, se, n_case, n_ctrl) {
  x <- 2 * n_case; y <- 2 * n_ctrl; z <- exp(beta); w0 <- se^2

  # gentle inflation like ReACt (non-cumulative)
  D <- NA_real_; w <- w0
  for (i in 0:49) {
    w <- w0 * (1.001^i)
    D <- (2*z*y*(1 - z) - w*x*y*z)^2 - 4*y*z*(y*z + x)*((z - 1)^2 + w*x*z)
    if (D >= 0) break
  }
  if (!is.finite(D) || D < 0) return(list(root1 = NULL, root2 = NULL))

  tmp1 <- w*x*y*z - 2*y*z*(1 - z)
  tmp2 <- sqrt(D)
  tmp3 <- 2 * ((z - 1)^2 + w*x*z)

  build <- function(d) {
    c <- y - d
    b <- x * d / (z*y - z*d + d)
    a <- x - b
    ok <- all(is.finite(c(a,b,c,d))) && all(c(a,b,c,d) > 0)
    if (!ok) return(NULL)
    list(
      AC_case = a, AC_ctrl = c,
      AF_case = a/(a+b), AF_ctrl = c/(c+d),
      AC = a+c, AF = (a+c)/(x+y)
    )
  }

  list(
    root1 = build((tmp1 - tmp2)/tmp3),
    root2 = build((tmp1 + tmp2)/tmp3)
  )
}

se_prior_overall_AF <- function(se, n_case, n_ctrl) {
  Neff <- 4*n_case*n_ctrl / (n_case + n_ctrl)
  a <- 2 / (Neff * se^2)
  a <- pmin(a, 0.25 - 1e-12)   # clamp; max of f(1-f) is 0.25
  (1 - sqrt(1 - 4*a)) / 2      # choose <= 0.5 root
}

# prior_type: "control" to match AF_ctrl, "overall" to match AF
# prior_allele: "allele2" if prior is for SAIGE Allele2 (effect allele),
#               "other"   if prior refers to the *other* allele (we'll flip 1-prior)
pick_root <- function(roots, prior, prior_type = c("control","overall"), prior_allele = c("allele2","other")) {
  prior_type <- match.arg(prior_type)
  prior_allele <- match.arg(prior_allele)
  if (is.null(roots$root1) && is.null(roots$root2)) return(NULL)
  if (prior_allele == "other") prior <- 1 - prior

  dist <- function(r) {
    if (is.null(r)) return(Inf)
    target <- if (prior_type == "control") r$AF_ctrl else r$AF
    abs(target - prior)
  }
  d1 <- dist(roots$root1); d2 <- dist(roots$root2)
  if (d1 <= d2) roots$root1 else roots$root2
}

reconstruct_analytic_prefer_correct <- function(beta, se, n_case, n_ctrl,
                                                prior_AF = NA_real_,
                                                prior_type = c("control","overall"),
                                                prior_allele = c("allele2","other")) {
  roots <- react_analytic_both(beta, se, n_case, n_ctrl)
  if (is.na(prior_AF)) {
    prior_AF <- se_prior_overall_AF(se, n_case, n_ctrl)
    prior_type <- "overall"
    prior_allele <- "allele2"
  } else {
    prior_type <- match.arg(prior_type)
    prior_allele <- match.arg(prior_allele)
  }
  pick_root(roots, prior_AF, prior_type, prior_allele)
}

reconstruct_df <- function(df,
                           beta="BETA", se="SE", ncase="N_case", nctrl="N_ctrl",
                           prior_col = NULL,         # e.g. "AF_ctrl_Allele2" or "AF_alt"
                           prior_type = "control",   # "control" or "overall"
                           prior_is_for = "allele2"  # "allele2" or "other"
                           ) {
  B <- df[[beta]]; S <- df[[se]]; Nc <- df[[ncase]]; Nt <- df[[nctrl]]
  if (!is.null(prior_col) && prior_col %in% names(df)) {
    P <- df[[prior_col]]
  } else {
    P <- rep(NA_real_, nrow(df))   # triggers SE-based overall prior
  }
  res <- lapply(seq_along(B), function(i)
    reconstruct_analytic_prefer_correct(B[i], S[i], Nc[i], Nt[i],
                                        prior_AF = P[i],
                                        prior_type = prior_type,
                                        prior_allele = prior_is_for))
  # bind
  to_row <- function(x) if (is.null(x)) rep(NA_real_, 6) else unlist(x)[c("AC_case","AC_ctrl","AF_case","AF_ctrl","AC","AF")]
  mat <- t(vapply(res, to_row, numeric(6)))
  colnames(mat) <- c("AC_case","AC_ctrl","AF_case","AF_ctrl","AC","AF")
  cbind(df, as.data.frame(mat))
}
