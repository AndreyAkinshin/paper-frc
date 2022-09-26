source("utils.R")
source("mad-factors.R")

knitr::opts_chunk$set(echo = FALSE, fig.align = "center", fig.pos = "ht!", fig.height = 3)
options(knitr.kable.NA = "-")
theme_set(theme_bw())
options(scipen = 999)

# Settings ---------------------------------------------------------------------
settings <- list(
  consistency = list(
    rebuild = FALSE,
    filename = "data-consistency.csv",
    repetitions1 = 25 * 1000 * 1000,
    repetitions2 = 2 * 1000 * 1000,
    ns = c(
      2:100,
      c(109, 110) + rep(seq(0, 90, by = 10), each = 2),
      c(249, 250) + rep(seq(0, 1050, by = 50), each = 2),
      c(1499, 1500) + rep(seq(0, 8500, by = 500), each = 2)
    )
  ),
  efficiency = list(
    rebuild = FALSE,
    filename = "data-efficiency.csv",
    repetitions1 = 10 * 1000 * 1000,
    repetitions2 = 2 * 1000 * 1000,
    ns = c(
      2:100,
      c(109, 110) + rep(seq(0, 90, by = 10), each = 2),
      c(249, 250) + rep(seq(0, 250, by = 50), each = 2),
      600, 700, 800, 900, 1000, 1500, seq(2000, 10000, by = 1000),
      25000, 50000, 100000
    )
  )
)

# Functions --------------------------------------------------------------------

format_int_latex <- function(n) {
  s <- format(n, scientific = FALSE)
  res <- ""
  while (nchar(s) > 0) {
    if (nchar(res) > 0) {
      res <- paste0("\\,", res)
    }
    res <- paste0(substr(s, nchar(s) - 2, nchar(s)), res)
    s <- substr(s, 1, nchar(s) - 3)
  }
  res
}

format_double_latex <- function(x, digits = 3, print_unary = TRUE) {
  s <- format(round(x, digits), scientific = FALSE)
  if (x < 0) s else paste0("+", s)
}

to_tex <- function(estimator) {
  unname(c("sn" = "$S_n$", "qn" = "$Q_n$")[estimator])
}

c4 <- function(n) ifelse(n < 300, sqrt(2 / (n - 1)) * gamma(n / 2) / gamma((n - 1) / 2), 1)
sd.unbiased <- function(x) sd(x) / c4(length(x))

sn_bias_rc <- function(n) {
  sapply(n, function(n) {
    if (n <= 9) {
      c(0.743, 1.851, 0.954, 1.351, 0.993, 1.198, 1.005, 1.131)[n - 1]
    } else if (n %% 2) {
      n / (n - 0.9)
    } else {
      1
    }
  })
}

qn_bias_rc <- function(n) {
  sapply(n, function(n) {
    Qn.finite.c <- function(n) {
      (if (n %% 2) {
        n / (n + 1.4)
      } # n is odd
      else {
        n / (n + 3.8)
      } # n is even
      )
    }
    if (n <= 9) {
      c(0.399, 0.994, 0.512, 0.844, 0.611, 0.857, 0.669, 0.872)[n - 1L]
    } else {
      Qn.finite.c(n)
    }
  })
}

qn_bias_rb <- function(n) {
  sapply(n, function(n) {
    Qn.finite.c <- function(n) {
      (if (n %% 2) {
        1.60188 + (-2.1284 - 5.172 / n) / n
      } # n is odd
      else {
        3.67561 + (1.9654 + (6.987 - 77 / n) / n) / n
      } # n is even
      ) / n + 1
    }
    if (n <= 12) {
      c(
        0.399356, 0.99365, 0.51321, 0.84401, 0.6122,
        0.85877, 0.66993, 0.87344, 0.72014, 0.88906,
        0.75743
      )[n - 1L]
    } else {
      1 / Qn.finite.c(n)
    }
  })
}


asymptotic_sn_constant <- uniroot(
  function(c) pnorm(qnorm(3 / 4) + 1 / c) - pnorm(qnorm(3 / 4) - 1 / c) - 1 / 2,
  c(1.19, 1.20),
  tol = 1e-15
)$root
asymptotic_qn_constant <- 1 / (sqrt(2) * qnorm(5 / 8))
get_asymptotic_constant <- function(estimator) {
  unname(c("sn" = asymptotic_sn_constant, "qn" = asymptotic_qn_constant)[estimator])
}

get_bias_coefficients <- function(estimator, parity) {
  df <- consistency_full
  df <- df[df$n > 100 & df$n <= 1000 & df$n %% 2 == parity, ]
  df$bias <- df[, paste0(estimator, "_factor")] - 1
  fit <- lm(bias ~ 0 + I(n^(-1)) + I(n^(-2)), data = df)
  fit$coefficients
}

get_factor <- function(estimator, ns, always_predict = FALSE) {
  sapply(ns, function(n) {
    df <- consistency_full
    if (n %in% df$n && !always_predict) {
      return(df[df$n == n, paste0(estimator, "_factor")])
    }

    coefficients <- get_bias_coefficients(estimator, n %% 2)
    1 + coefficients[1] / n + coefficients[2] / n^2
  })
}

Sn_new <- function(x) {
  factor <- get_factor("sn", length(x))
  robustbase::Sn(x, asymptotic_sn_constant * factor)
}

Qn_new <- function(x) {
  factor <- get_factor("qn", length(x))
  robustbase::Qn(x, asymptotic_qn_constant * factor)
}

# Data -------------------------------------------------------------------------

simulation_consistency <- function(rebuild = NULL, filename = NULL, repetitions = NULL, ns = NULL) {
  apply_settings(settings$consistency)

  process <- function(n) {
    repetitions <- if (n <= 100) repetitions1 else repetitions2
    sn_constant <- 1 / mean(future_replicate(repetitions, robustbase::Sn(rnorm(n), 1)))
    qn_constant <- 1 / mean(future_replicate(repetitions, robustbase::Qn(rnorm(n), 1)))
    data.frame(
      n = n,
      constant.sn = round(sn_constant, 6),
      constant.qn = round(qn_constant, 6),
      factor.sn = round(sn_constant / asymptotic_sn_constant, 6),
      factor.qn = round(qn_constant / asymptotic_qn_constant, 6)
    )
  }

  df <- multi_estimate(rebuild, filename, ns, process)
  consistency_full <<- data_full()
  df
}

simulation_efficiency <- function(rebuild = NULL, filename = NULL, repetitions = NULL, ns = NULL) {
  apply_settings(settings$efficiency)

  estimate <- function(x) {
    c(
      n = length(x),
      sdn = sd.unbiased(x),
      madn = mad.sm(x),
      sn = Sn_new(x),
      qn = Qn_new(x)
    )
  }

  process <- function(n) {
    repetitions <- if (n <= 100) repetitions1 else repetitions2
    df_n <- data.frame(t(future_replicate(repetitions, estimate(rnorm(n)))))

    df <- data.frame(
      n = n,
      bias.sdn = mean(df_n$sdn),
      bias.madn = mean(df_n$madn),
      bias.sn = mean(df_n$sn),
      bias.qn = mean(df_n$qn),
      svar.sdn = n * var(df_n$sdn) / mean(df_n$sdn)^2,
      svar.madn = n * var(df_n$madn) / mean(df_n$madn)^2,
      svar.sn = n * var(df_n$sn) / mean(df_n$sn)^2,
      svar.qn = n * var(df_n$qn) / mean(df_n$qn)^2
    )
    df$eff.madn <- df$svar.sdn / df$svar.madn
    df$eff.sn <- df$svar.sdn / df$svar.sn
    df$eff.qn <- df$svar.sdn / df$svar.qn

    round(df, 6)
  }

  multi_estimate(rebuild, filename, ns, process)
}

data_full <- function() {
  if (!file.exists(settings$consistency$filename)) {
    return(data.frame())
  }
  df_constants <- read.csv(settings$consistency$filename) %>% extract_columns("constant")
  data.frame(
    n = df_constants$n,
    parity = ifelse(df_constants$n %% 2 == 0, "Even", "Odd"),
    sn_constant = df_constants$sn,
    sn_factor = df_constants$sn / asymptotic_sn_constant,
    sn_bias = df_constants$sn / asymptotic_sn_constant - 1,
    qn_constant = df_constants$qn,
    qn_factor = df_constants$qn / asymptotic_qn_constant,
    qn_bias = df_constants$qn / asymptotic_qn_constant - 1
  )
}
consistency_full <- data_full()

data_factor_compare_reference <- function(estimator) {
  df_full <- consistency_full
  df <- data.frame(n = df_full$n, new = df_full[, paste0(estimator, "_factor")])
  if (estimator == "sn") {
    df$rc <- sn_bias_rc(df$n)
    df$diff_rc <- abs(df$new - df$rc)
  }
  if (estimator == "qn") {
    df$rc <- qn_bias_rc(df$n)
    df$diff_rc <- abs(df$new - df$rc)
    df$rb <- qn_bias_rb(df$n)
    df$diff_rb <- abs(df$new - df$rb)
  }
  df
}

data_factor_compare_predicted <- function(estimator) {
  df_full <- consistency_full
  df <- data.frame(n = df_full$n, actual = df_full[, paste0(estimator, "_factor")])
  df <- df[df$n > 100, ]
  df$predicted <- get_factor(estimator, df$n, TRUE)
  df$diff <- abs(df$actual - df$predicted)
  df
}

# Inlines ----------------------------------------------------------------------

inline_factor_predicted_diff <- function(estimator) {
  df <- data_factor_compare_predicted(estimator)
  format(round(max(df$diff), 6), scientific = FALSE)
}

inline_factor_rc_diff <- function(estimator, type = "d") {
  df <- data_factor_compare_reference(estimator)
  max_diff <- max(df$diff_rc)
  if (type == "d") {
    format(round(max_diff, 6), scientific = FALSE)
  } else if (type == "n") {
    df[df$diff_rc == max_diff, ]$n
  } else {
    stop(paste("Unrecognized type:", type))
  }
}

inline_factor_rb_diff <- function(type = "d", minn = 13) {
  df <- data_factor_compare_reference("qn")
  df <- df[df$n >= minn, ]
  max_diff <- max(df$diff_rb)
  if (type == "d") {
    format(round(max_diff, 6), scientific = FALSE)
  } else if (type == "n") {
    df[df$diff_rb == max_diff, ]$n
  } else {
    stop(paste("Unrecognized type:", type))
  }
}

# Tables -----------------------------------------------------------------------

table_factors <- function() {
  df <- consistency_full
  df <- data.frame(
    n = df$n,
    cn = df$sn_factor,
    dn = df$qn_factor
  )

  size <- 50
  caption <- "Finite-sample bias-correction factors for $\\Sn$ and $\\Qn$."

  if (nrow(df[df$n == 1, ]) == 0) {
    df <- rbind(data.frame(n = 1, cn = NA, dn = NA), df)
  }

  slice <- function(index) df[((index - 1) * size + 1):(index * size), ]

  slice_count <- ceiling(nrow(df) / size)
  df2 <- slice(1)
  for (i in 2:slice_count) {
    df2 <- cbind(df2, slice(i))
  }

  header <- rep(c("n", "$c_n$", "$d_n$"), slice_count)
  knitr::kable(df2, caption = caption, col.names = header, escape = FALSE, digits = 4) %>%
    kable_styling(latex_options = "hold_position")
}

table_efficiency <- function() {
  df <- read.csv("data-efficiency.csv") %>% extract_columns("eff")
  caption <- "Finite-sample Gaussian efficiency of $\\MAD_n$, $\\Sn$, $\\Qn$."

  if (nrow(df[df$n == 2, ]) == 0) {
    df <- rbind(data.frame(n = 2, madn = NA, sn = NA, qn = NA), df)
  }
  if (nrow(df[df$n == 1, ]) == 0) {
    df <- rbind(data.frame(n = 1, madn = NA, sn = NA, qn = NA), df)
  }

  size <- 50
  slice <- function(index) df[((index - 1) * size + 1):(index * size), ]

  slice_count <- ceiling(nrow(df) / size)
  df2 <- slice(1)
  for (i in 2:slice_count) {
    df2 <- cbind(df2, slice(i))
  }

  header <- rep(c("n", "$\\MAD_n$", "$\\Sn$", "$\\Qn$"), slice_count)

  knitr::kable(df2, caption = caption, col.names = header, escape = FALSE, digits = 4) %>%
    kable_styling(latex_options = "hold_position")
}

# Figures ----------------------------------------------------------------------

figure_factors1 <- function(estimator) {
  df <- consistency_full
  df <- df[df$n <= 100, ]
  df$value <- df[, paste0(estimator, "_factor")]
  color <- if (estimator == "sn") cbp$blue else cbp$green

  ggplot(df, aes(n, value, color = parity)) +
    geom_hline(yintercept = 1, linetype = "dashed", col = cbp$grey) +
    geom_point(aes(shape = parity), col = color) +
    labs(
      title = TeX(paste0("(a) Bias-correction factors for ", to_tex(estimator), " ($n \\leq 100$)")),
      x = "Sample size (n)",
      y = "Bias-correction factor",
      col = "Parity of n",
      shape = "Parity of n"
    ) +
    scale_color_manual(values = cbp$values)
}

figure_factors2 <- function(estimator) {
  df <- consistency_full
  df <- df[df$n > 100, ]
  df$value <- df[, paste0(estimator, "_factor")]
  df$predicted <- get_factor(estimator, df$n, TRUE)
  color <- if (estimator == "sn") cbp$blue else cbp$green

  ggplot(df, aes(n, value, col = parity, shape = parity)) +
    geom_hline(yintercept = 1, linetype = "dashed", col = cbp$grey) +
    geom_line(aes(n, predicted), alpha = 0.5, col = color) +
    geom_point(col = color) +
    labs(
      title = TeX(paste0("(b) Bias-correction factors for ", to_tex(estimator), ": actual and predicted (n > 100)")),
      x = "Sample size (n)",
      y = "Bias-correction factor",
      col = "Parity of n",
      shape = "Parity of n"
    )
}

figure_factors <- function(estimator) {
  p1 <- figure_factors1(estimator)
  p2 <- figure_factors2(estimator)
  ggarrange(p1, p2, ncol = 1)
}

figure_efficiency1 <- function() {
  df <- simulation_efficiency() %>%
    extract_columns("eff") %>%
    gather("type", "value", -n)
  df$type <- factor(df$type, levels = c("madn", "sn", "qn"))
  df$parity <- factor(ifelse(df$n %% 2 == 0, "Even", "Odd"), levels = c("Even", "Odd"))
  df <- df[3 <= df$n & df$n <= 100, ]
  y_breaks <- c(0.3675, 0.4, 0.5, 0.5823, 0.6, 0.7, 0.8, 0.8227)
  ggplot(df, aes(n, value, col = type, shape = parity, linetype = parity)) +
    geom_hline(yintercept = 0.3675, col = cbp$red, linetype = "dotted", alpha = 0.5) +
    geom_hline(yintercept = 0.5823, col = cbp$blue, linetype = "dotted", alpha = 0.5) +
    geom_hline(yintercept = 0.8227, col = cbp$green, linetype = "dotted", alpha = 0.5) +
    geom_line(alpha = 0.3) +
    geom_point() +
    labs(
      title = TeX("(a) Gaussian efficiency of $\\MAD_n, \\S_n, \\Q_n\\; (3 \\leq n \\leq 100)$"),
      x = "Sample size (n)",
      y = "Gaussian efficiency",
      col = "Estimator",
      shape = "Parity of n",
      linetype = "Parity of n"
    ) +
    scale_color_manual(values = cbp$values, labels = c(TeX("$\\MAD_n$"), TeX("$\\S_n$"), TeX("$\\Q_n$"))) +
    scale_y_continuous(breaks = y_breaks)
}

figure_efficiency2 <- function() {
  df <- simulation_efficiency() %>%
    extract_columns("eff") %>%
    gather("type", "value", -n)
  df$type <- factor(df$type, levels = c("madn", "sn", "qn"))
  df$parity <- factor(ifelse(df$n %% 2 == 0, "Even", "Odd"), levels = c("Even", "Odd"))
  df <- df[100 <= df$n & df$n <= 1000, ]
  x_breaks <- seq(100, 1000, by = 100)
  y_breaks <- c(0.3675, 0.4, 0.5, 0.5823, 0.6, 0.7, 0.8, 0.8227)
  ggplot(df, aes(n, value, col = type, shape = parity, linetype = parity)) +
    geom_hline(yintercept = 0.3675, col = cbp$red, linetype = "dotted", alpha = 0.5) +
    geom_hline(yintercept = 0.5823, col = cbp$blue, linetype = "dotted", alpha = 0.5) +
    geom_hline(yintercept = 0.8227, col = cbp$green, linetype = "dotted", alpha = 0.5) +
    geom_point(size = 1.5) +
    labs(
      title = TeX("(b) Gaussian efficiency of $\\MAD_n, \\S_n, \\Q_n\\; (100 \\leq n \\leq 1000)$"),
      x = "Sample size (n)",
      y = "Gaussian efficiency",
      col = "Estimator",
      shape = "Parity of n",
      linetype = "Parity of n"
    ) +
    scale_color_manual(values = cbp$values, labels = c(TeX("$\\MAD_n$"), TeX("$\\S_n$"), TeX("$\\Q_n$"))) +
    scale_x_continuous(breaks = x_breaks) +
    scale_y_continuous(breaks = y_breaks)
}

figure_efficiency <- function() {
  p1 <- figure_efficiency1()
  p2 <- figure_efficiency2()
  ggarrange(p1, p2, ncol = 1)
}

# Simulations ------------------------------------------------------------------
df_sc <- simulation_consistency()
df_se <- simulation_efficiency()
