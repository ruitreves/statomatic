# Simulations for the "Approach and limitations" vignette.
#
# They measure, with statomatic's own functions:
#   1. how often test_norm() sorts rows into the right set (ne, nu, nn)
#   2. the false-positive rate of ANOVA, Welch and Kruskal-Wallis when the
#      group means are equal but the variances are not
#   3. the power of the three tests when the variances are equal
#
# Run from the package root to update the results used by the vignette:
#     Rscript inst/simulations/sorting_simulations.R
# The script is installed with the package:
#     system.file("simulations", "sorting_simulations.R", package = "statomatic")

library(statomatic)
set.seed(20261019)

out_file <- "vignettes/simulation_results.rds"
n_rows <- 500      # simulated rows per scenario for the sorting simulation
n_reps <- 5000     # simulated rows per scenario for error rates and power

groups_of <- function(sizes) rep(paste0("g", seq_along(sizes)), sizes)

# rows x samples matrix, each row drawn by draw(n_total, group labels)
simulate <- function(n, sizes, draw) {
    g <- groups_of(sizes)
    x <- t(replicate(n, draw(length(g), g)))
    dimnames(x) <- list(paste0("row", seq_len(n)), paste0("s", seq_along(g)))
    list(x = x, g = g)
}

# 1. Sorting accuracy ---------------------------------------------------------
# Three groups of equal size. One group's SD is `sd_ratio` times the others'.
truths <- list(
    ne = function(N, g) stats::rnorm(N, 10, 1),
    nu = function(N, g, sd_ratio) stats::rnorm(N, 10, ifelse(g == "g3", sd_ratio, 1)),
    nn = function(N, g) stats::rlnorm(N, 0, 0.8)                 # skewed
)
sorting <- do.call(rbind, lapply(c(4, 5, 6, 8, 12, 20), function(n_per) {
    do.call(rbind, lapply(c(3, 5), function(sd_ratio) {
        do.call(rbind, lapply(names(truths), function(truth) {
            draw <- if (truth == "nu") function(N, g) truths$nu(N, g, sd_ratio) else truths[[truth]]
            sim <- simulate(n_rows, rep(n_per, 3), draw)
            set <- test_norm(sim$x, sim$g)$tests$set
            pct <- 100 * as.numeric(table(set)) / n_rows
            data.frame(n_per_group = n_per, sd_ratio = sd_ratio, truth = truth,
                       ne = pct[1], nu = pct[2], nn = pct[3])
        }))
    }))
}))
# the ne and nn scenarios don't depend on sd_ratio; keep one copy of each
sorting <- sorting[sorting$truth == "nu" | sorting$sd_ratio == 3, ]
sorting$sd_ratio[sorting$truth != "nu"] <- NA
rownames(sorting) <- NULL

# 2. False-positive rates -----------------------------------------------------
# All group means are equal, so every p-value below 0.05 is a false positive.
p_values <- function(sim) {
    cbind(anova = run_anova(sim$x, sim$g, padj = FALSE)$pvalue,
          welch = run_welch(sim$x, sim$g, padj = FALSE)$pvalue,
          kruskal_wallis = run_kruskal(sim$x, sim$g, padj = FALSE)$pvalue)
}
designs <- list(c(4, 4, 4), c(6, 6, 6), c(6, 5, 4), c(12, 8, 5))
false_positives <- do.call(rbind, lapply(designs, function(sizes) {
    do.call(rbind, lapply(c(1, 3, 5), function(sd_ratio) {
        sim <- simulate(n_reps, sizes, function(N, g) stats::rnorm(N, 10, ifelse(g == "g3", sd_ratio, 1)))
        rate <- 100 * colMeans(p_values(sim) < 0.05)
        data.frame(group_sizes = paste(sizes, collapse = "/"),
                   sds = paste(c(1, 1, sd_ratio), collapse = "/"), t(rate))
    }))
}))

# 3. Power with equal variances -----------------------------------------------
# Group 3's mean is 2 SD higher.
power <- do.call(rbind, lapply(c(4, 5, 6, 8), function(n_per) {
    sim <- simulate(n_reps, rep(n_per, 3), function(N, g) stats::rnorm(N, ifelse(g == "g3", 2, 0), 1))
    rate <- 100 * colMeans(p_values(sim) < 0.05)
    data.frame(n_per_group = n_per, t(rate))
}))

results <- list(sorting = sorting, false_positives = false_positives, power = power,
                settings = list(n_rows = n_rows, n_reps = n_reps, seed = 20261019,
                                statomatic = as.character(utils::packageVersion("statomatic")),
                                r_version = R.version.string, date = Sys.Date()))
saveRDS(results, out_file)
cat("Saved simulation results to", out_file, "\n")
