# irtfit() computes observed proportions, residuals, and fit statistics for each item

prm_file <- system.file("extdata", "flexmirt_sample-prm.txt", package = "irtQ")
x_mix <- bring.flexmirt(file = prm_file, "par")$Group1$full_df[c(1:15, 51:55), ]
set.seed(2026)
theta_fit <- rnorm(1000)
resp_fit <- simdat(x = x_mix, theta = theta_fit, D = 1)

test_that("irtfit() observed proportions are the category frequencies divided by the group size", {
  for (gm in c("equal.width", "equal.freq")) {
    fit <- irtfit(
      x = x_mix, score = theta_fit, data = resp_fit, group.method = gm,
      n.width = 8, loc.theta = "average", range.score = c(-4, 4), D = 1
    )
    for (i in seq_along(fit$contingency.plot)) {
      tb <- fit$contingency.plot[[i]]
      freq <- tb[, grep("^obs\\.freq\\.[0-9]+$", names(tb)), drop = FALSE]
      prop <- tb[, grep("^obs\\.prop\\.[0-9]+$", names(tb)), drop = FALSE]
      # proportions in each interval sum to one
      expect_equal(unname(rowSums(prop)), rep(1, nrow(tb)))
      # proportions equal the frequencies over the interval totals
      expect_equal(unname(as.matrix(prop)), unname(as.matrix(freq)) / tb$total)
    }
  }
})

test_that("irtfit() raw residuals equal observed proportions minus expected probabilities", {
  fit <- irtfit(
    x = x_mix, score = theta_fit, data = resp_fit, group.method = "equal.freq",
    n.width = 8, loc.theta = "average", range.score = c(-4, 4), D = 1
  )
  tb <- fit$contingency.plot[[1]]
  prop <- tb[, grep("^obs\\.prop\\.[0-9]+$", names(tb))]
  expp <- tb[, grep("^exp\\.prob\\.[0-9]+$", names(tb))]
  rsd <- tb[, grep("^raw\\.rsd\\.[0-9]+$", names(tb))]
  expect_equal(unname(as.matrix(rsd)), unname(as.matrix(prop - expp)))
})

test_that("irtfit() drops empty score groups and keeps proportions that sum to one", {
  # two extreme scores leave empty intervals in the equal-width grouping
  sc <- theta_fit
  sc[1] <- -6
  sc[2] <- 6
  fit <- irtfit(
    x = x_mix, score = sc, data = resp_fit, group.method = "equal.width",
    n.width = 10, loc.theta = "average", range.score = c(-7, 7), D = 1
  )
  n_group <- vapply(fit$contingency.plot, nrow, integer(1))
  expect_true(any(n_group < 10L))
  for (tb in fit$contingency.plot) {
    prop <- tb[, grep("^obs\\.prop\\.[0-9]+$", names(tb))]
    expect_equal(unname(rowSums(prop)), rep(1, nrow(tb)))
  }
})

test_that("plot.irtfit() Wald intervals use the two-sided critical value", {
  # the last plot is read back, which needs get_last_plot() in ggplot2
  skip_if_not("get_last_plot" %in% getNamespaceExports("ggplot2"))
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  fit <- irtfit(
    x = x_mix, score = theta_fit, data = resp_fit, group.method = "equal.freq",
    n.width = 8, loc.theta = "average", range.score = c(-4, 4), D = 1, alpha = 0.05
  )
  plot(x = fit, item.loc = 1, type = "icc", ci.method = "wald", show.table = FALSE)
  built <- ggplot2::ggplot_build(ggplot2::get_last_plot())
  # the segment layer holds the interval (y = upper limit, yend = lower limit)
  seg <- Filter(function(d) all(c("y", "yend", "xend") %in% names(d)), built$data)[[1]]
  unclipped <- seg$yend > 0 & seg$y < 1
  half_width <- (seg$y - seg$yend)[unclipped] / 2
  tb <- fit$contingency.plot[[1]]
  se_all <- c(tb$se.0, tb$se.1)
  # each half width equals the 97.5th percentile of the normal times a standard error
  expect_true(length(half_width) > 0)
  expect_true(all(vapply(half_width, function(h) any(abs(h / stats::qnorm(0.975) - se_all) < 1e-6), logical(1))))
})

# ---- Input checks ---------------------------------------------------------------

x_ft <- x_mix[c(1:3, 20), ]
resp_ft <- resp_fit[, c(1:3, 20)]

test_that("irtfit() stops for data that do not match the items or the score categories", {
  # a column is missing
  expect_error(irtfit(x_ft, score = theta_fit, data = resp_ft[, 1:3], D = 1), "number of columns")

  # responses above the highest category and non-integer responses
  r <- resp_ft
  r[1, 1] <- 2
  expect_error(irtfit(x_ft, score = theta_fit, data = r, D = 1), "outside the score categories")
  r <- resp_ft
  r[2, 4] <- 0.5
  expect_error(irtfit(x_ft, score = theta_fit, data = r, D = 1), "outside the score categories")

  # a code declared with 'missing' is accepted and gives the result for NA
  r <- resp_ft
  r[1, 1] <- -9
  expect_error(irtfit(x_ft, score = theta_fit, data = r, D = 1), "outside the score categories")
  expect_identical(
    irtfit(x_ft, score = theta_fit, data = r, D = 1, missing = -9)$fit_stat,
    irtfit(x_ft, score = theta_fit, data = replace(r, r == -9, NA), D = 1)$fit_stat
  )
})

test_that("irtfit() reads character and factor responses like numeric responses", {
  r_na <- resp_ft
  r_na[1, 2] <- NA
  r_chr <- matrix(as.character(r_na), nrow(r_na))
  r_chr[is.na(r_na)] <- "."
  df_chr <- as.data.frame(r_chr, stringsAsFactors = FALSE)
  df_fac <- as.data.frame(lapply(df_chr, factor))
  ref <- irtfit(x_ft, score = theta_fit, data = r_na, D = 1)
  for (dat in list(r_chr, df_chr, df_fac)) {
    fit <- irtfit(x_ft, score = theta_fit, data = dat, D = 1, missing = ".")
    expect_identical(fit[names(fit) != "call"], ref[names(ref) != "call"])
  }

  # a response that cannot be read as a number stops
  r_chr[1, 1] <- "a"
  expect_error(irtfit(x_ft, score = theta_fit, data = r_chr, D = 1, missing = "."), "numeric scores")
})

test_that("irtfit() drops examinees with missing ability estimates", {
  sc <- theta_fit
  sc[c(5, 10)] <- NA
  fit <- irtfit(x_ft, score = sc, data = resp_ft, D = 1)
  expect_equal(fit$fit_stat$N, rep(998, 4))
  ref <- irtfit(x_ft, score = theta_fit[-c(5, 10)], data = resp_ft[-c(5, 10), ], D = 1)
  expect_equal(fit$fit_stat, ref$fit_stat)
})

test_that("irtfit() keeps working when only one item remains after excluding items", {
  resp2 <- resp_ft[, 1:2]
  resp2[, 2] <- NA
  expect_warning(
    fit <- irtfit(x_ft[1:2, ], score = theta_fit, data = resp2, D = 1),
    "fewer than two responses.*V2|fewer than two responses.*column 2"
  )
  expect_equal(fit$fit_stat$id, x_ft$id[1])
  expect_identical(fit$fit_stat, irtfit(x_ft[1, ], score = theta_fit, data = resp_ft[, 1], D = 1)$fit_stat)

  # the warning names the id and the column of the excluded item
  expect_warning(irtfit(x_ft[1:2, ], score = theta_fit, data = resp2, D = 1), paste0(x_ft$id[2], " [(]column 2[)]"))

  # every item without responses stops
  resp3 <- resp_ft[, 1:2]
  resp3[] <- NA
  expect_error(irtfit(x_ft[1:2, ], score = theta_fit, data = resp3, D = 1), "fewer than two responses")
})

test_that("irtfit() stops for invalid grouping and score arguments", {
  # a misspelled location of the theta point
  expect_error(irtfit(x_ft, score = theta_fit, data = resp_ft, D = 1, loc.theta = "midle"), "should be one of")
  expect_identical(
    irtfit(x_ft, score = theta_fit, data = resp_ft, D = 1, loc.theta = "MIDDLE")$fit_stat,
    irtfit(x_ft, score = theta_fit, data = resp_ft, D = 1, loc.theta = "middle")$fit_stat
  )

  # scores that do not match the rows of the data
  expect_error(irtfit(x_ft, score = theta_fit[1:400], data = resp_ft, D = 1), "length of 'score'")

  # scores without any variation
  expect_error(irtfit(x_ft, score = rep(0.5, 1000), data = resp_ft, D = 1), "no variation")

  # a reversed or incomplete score range
  expect_error(irtfit(x_ft, score = theta_fit, data = resp_ft, D = 1, range.score = c(3, -3)), "range.score")
  expect_error(irtfit(x_ft, score = theta_fit, data = resp_ft, D = 1, range.score = 3), "range.score")
})

test_that("plot.irtfit() stops for an unknown type or an item location outside the evaluated items", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  fit <- irtfit(x_ft, score = theta_fit, data = resp_ft, D = 1)
  expect_error(plot(fit, item.loc = 1, type = "bothh"), "should be one of")
  expect_error(plot(fit, item.loc = 9), "item.loc")
  expect_error(plot(fit, item.loc = 0), "item.loc")
  expect_error(plot(fit, item.loc = 1.5), "item.loc")
  expect_no_error(plot(fit, item.loc = 4, type = "SR", show.table = FALSE))
})
