# --- test_numerical_equivalence.R ------------------------------------------
# Asserts the vectorised covariate regression in 1.2_setup_expression.R gives the same
# residuals as the original per-gene loop.
#
#   in : nothing -- generates its own random matrices
#   out: stops on mismatch
#   run: Rscript scripts/tests/test_numerical_equivalence.R
# --------------------------------------------------------------------------
set.seed(42)

expression <- matrix(rnorm(240), nrow = 20, ncol = 12)
covariates <- matrix(rnorm(36), nrow = 12, ncol = 3)

loop_residuals <- matrix(nrow = nrow(expression), ncol = ncol(expression))
for (row in seq_len(nrow(expression))) {
  fit <- summary(lm(expression[row, ] ~ covariates))
  loop_residuals[row, ] <- scale(fit$residuals)
}

design <- cbind(1, covariates)
design_qr <- qr(design)
residuals_by_gene <- t(expression) - design %*% qr.coef(design_qr, t(expression))
qr_residuals <- t(scale(residuals_by_gene))

stopifnot(isTRUE(all.equal(loop_residuals, qr_residuals, tolerance = 1e-12,
                           check.attributes = FALSE)))

annotation_ids <- c("ENSG000001", "ENSG0000010")
stopifnot(identical(match("ENSG000001", annotation_ids), 1L))

suppressPackageStartupMessages(library(data.table))
parts <- list(data.table(gene = "g1", value = 1), data.table(gene = "g2", value = 2))
stopifnot(isTRUE(all.equal(do.call(rbind, parts), rbindlist(parts))))

cat("Round 4 numerical equivalence tests passed\n")
