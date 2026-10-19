# statomatic rebuild — handoff

Last updated 2026-10-19. Branch `rebuild`, last commit `5673d8e`.

**To resume with Claude Code:** open this repo, run `git pull` on `rebuild`, and tell Claude:
> Read dev/HANDOFF.md and pick up where it leaves off. Save the "Working rules" and "Direction" sections to your memory.

---

## Working rules (for Claude)

- **Suggest only.** Never create, edit or delete project files without my explicit OK for that change.
- **Commit messages go in the chat** as copyable text, with the `git add` command. Never write them to a file, never use `git commit -F`. I do the add and commit myself.
- **The original statomatic must keep working exactly as it does now.** Don't touch existing `R/` files, exports or old tests until I say we're integrating. New functions live next to the old ones.
- **Don't run `devtools::document()`.** roxygen2 8.1.0 rewrites `DESCRIPTION` (RoxygenNote) and `man/p.adjustment.methods.Rd` no matter what. Use `devtools::load_all()` and `devtools::test()`.
- If the chat drifts, steer back to the next step below.
- I work on two computers: `git pull` at the start of a session, `git push` at the end.

## Direction (decided 2026-10-15)

**A faithful rebuild of statomatic** that keeps its pre-test sorting approach, with updated documentation of the approach and its limitations. statomatic is my first first-author publication; don't re-argue replacing it with DESeq2/edgeR/limma. "Faithful" still allows improvements: Type II ANOVA, name-based sample matching, the extra `new_test_norm` output, `log_input` for fold changes.

## What statomatic does

For each gene (row), test whether the model residuals are normal (Shapiro-Wilk) and the group variances are equal (Levene). That sorts genes into three sets:

| Set | Meaning | 2 groups | >2 groups: omnibus | >2 groups: post-hoc |
|---|---|---|---|---|
| ne | normal, equal variance | t-test | ANOVA | Tukey |
| nu | normal, unequal variance | Welch | Welch | Dunnett T3 |
| nn | not normal | Wilcoxon | Kruskal-Wallis | Dunn |

In the multigroup workflow, each omnibus result is combined with its post-hoc result and the fold changes.

## Conventions in the new code

- New functions are **internal**: roxygen header ending in `@noRd`, no `@export`. Change to `@export` only at integration.
- Names: `new_ttest`, `new_run_<test>`, `new_*_main`. Files: `R/new_*.R`. Tests: `tests/testthat/test-new_*.R`. Never reuse an existing file name in `R/` (`dunnett_test.R`, `tukey_test.R`, `welch_test.R`, `wilcox_test.R` already exist).
- All p-value results go through `make_res()` (`R/new_utils.R`):
  - named vector (one p-value per gene) → columns `pvalue`, `padj`
  - matrix (one column per comparison/term) → `pvalue_<name>`, `padj_<name>`
  - `padj` = Benjamini-Hochberg across genes; optional, default `TRUE`
- **Factors and groups.** The main functions take the sample table and the names of its factor columns: `new_multigroup_main(x, sample_info, c("condition", "sex"))` (1–3 factors), `new_twogroup_main(x, sample_info, "condition")`. `select_factors()` finds the sample IDs (`sample_col` if given, else the row names, else the one column whose values contain every sample of `x`), drops rows for samples not in `x`, and gives clear errors for typos, ambiguous or missing IDs, and for passing a column (`sample_info$group`) instead of its name. `check_factors()` then validates the design and builds the combined group (`condition_sex`, e.g. `Control_F`) with **radix-sorted levels**. ANOVA and `new_test_norm` use the factors; Welch, Kruskal-Wallis, all post-hoc tests and fold changes use the combined group.
- **Comparison names** are `A_vs_B`, with A before B in the group's level order, in every function (Tukey, Dunnett, Dunn, fold changes). `new_run_dunn` re-orients `Adunn.test`'s comparisons to achieve this. Fold changes are attached to p-values **by comparison name** (`add_fold_changes()`), never by position.
- `include_test_name` defaults to **FALSE for omnibus tests** (ANOVA, Welch, KW) because the post-hoc `test` column already shows the path, and TRUE for post-hoc tests. `new_twogroup_main` defaults it to TRUE (no post-hoc there).
- Post-hoc p-values (Tukey, Dunnett T3, Dunn) are already multiple-comparison corrected within each gene; that is always on. `padj` adds across-gene BH on top, **within each set** (ne, nu, nn) as in the original.
- Old functions use `Pval`/`Padj` and default `padj = FALSE`; new ones use `pvalue`/`padj` and default `padj = TRUE`.

## Done (all tested)

| Function | File | Notes |
|---|---|---|
| `make_res`, `select_factors`, `check_factors`, `add_fold_changes` | `R/new_utils.R` | shared helpers. `check_factors(..., two_group = TRUE)` for the 2-group workflow |
| `new_ttest`, `new_run_welch`, `new_run_wilcox` | `R/new_two_group.R` | Welch also works for >2 groups; Wilcoxon suppresses the ties warning like `run_wilcox` |
| `new_run_tukey`, `new_run_dunnett`, `new_run_dunn` | `R/new_multi_group.R` | Dunn requires >2 groups |
| `dunn_bonferroni`, `dunnett_t3` | `R/new_multi_group.R` | own implementations (~20 lines each) used by `new_run_dunn` / `new_run_dunnett`; identical results to the vendored `Adunn.test` / `AdunnettT3Test` (same seed for Dunnett), incl. ties and constant groups (p = 1) |
| `new_run_anova`, `new_run_kw` | `R/new_anovas.R` | ANOVA is **always Type II** via `car::Anova` |
| `new_test_norm` | `R/new_test_norm.R` | returns `ne`/`nu`/`nn` + `tests` (per-gene Shapiro/Levene p, set); separate thresholds; NAs dropped per gene; constant genes not handled (filter beforehand) |
| `new_fold_change` | `R/new_fold_change.R` | log2 of ratio of means; `log_input = TRUE` → difference of means; NAs ignored; no `FUN` argument |
| `new_multigroup_main` | `R/new_multigroup_main.R` | one-way results identical to `multigroup_main`; `verbose`; 0-row tables for empty sets; no `write_files` |
| `new_twogroup_main` | `R/new_twogroup_main.R` | results identical to `twogroup_main`; same conventions |

Result table layout (both main functions): omnibus `pvalue`/`padj` (per ANOVA term for multi-way), then for each comparison `pvalue`, `padj`, `log2FoldChange`, then `<group>_mean` columns, then `test`.

Tests: one `test-new_*.R` per R file, helpers in `helper-new.R` (`with_seed()`, `seeded_dummy_data()`, `quietly()`). Last run: 73 test blocks, 0 failures. Old tests all pass.

Dev folder: `load_data.R`, `new_stuff.R` (scratch), `baseline_tests.txt` (old test results before the rebuild), `counts.csv` / `sample_info.csv` (anonymized; sample_info is made up). The main functions find the `sample` column automatically.

## Gotchas found so far

- Old tests use `dummy_data()` with no seed, so they're random. `run_dunn` crashes whenever the random data has exactly 2 groups — pre-existing, left alone.
- Dunnett T3 p-values come from Monte Carlo integration (`mvtnorm::pmvt`): repeated runs differ by up to ~0.0003. Compare old vs new under the same seed.
- **R 4.6.0 changed `wilcox.test`:** with ties it now computes exact conditional p-values; before 4.6 it fell back to the normal approximation (with a warning). Wilcoxon p-values for tied data therefore differ between R versions.
- **dunn.test 1.4.1 gives different p-values than the 1.3.6 copied into statomatic**, so calling the package directly would change results; `dunn_bonferroni` reproduces 1.3.6 exactly. Unmodified `PMCMRplus::dunnettT3Test` still errors on two identical constant groups; `dunnett_t3` keeps the fix.
- `apply()` on one row, or `as.data.frame()` on a vector, changes shapes and names — the cause of earlier `pvalue_res` and one-row `colnames` bugs.
- For 2+ factor ANOVA, Type II intentionally differs from the old `run_anova` (Type I), and the post-hoc tests now use the combined group (the old default was the first factor), so only one-way designs can be compared to `multigroup_main`.
- `dev/new_stuff.R` defines old drafts of `add_fold_changes` and `new_multigroup_main` (with a `group` argument). Sourcing it after `load_all()` replaces the package versions — delete those definitions.

## Next steps

1. **Integration** (started 2026-10-19). Checklist:
   - Export the new functions under their final names (`@noRd` → `@export`); internal helpers (`make_res`, `select_factors`, `check_factors`, `add_fold_changes`, `dunn_bonferroni`, `dunnett_t3`) stay internal.
   - Remove the old implementations they replace, and their old tests.
   - Remove the S4 layer (decided 2026-10-14): `class_setup.R`, `methods.R`, `sds_analyze.R` (`make_sds`, `sds_*`, `get_*`), and `import(SummarizedExperiment)` from `NAMESPACE` (never declared in `DESCRIPTION`). Add one line in the docs showing how to use a SummarizedExperiment anyway.
   - Remove the vendored copies `R/altered_dunn.test.R` and `R/altered_dunnettT3test.R` (`Adunn.test`, `p.adjustment.methods`, `AdunnettT3Test`; ~900 lines of GPL code in an MIT package).
   - Remove `combine_fc_pvals` (replaced by `add_fold_changes`), and `chop` and `takeN` (deprecated, decided 2026-10-19). Keep `cf` and `dummy_data`.
   - `DESCRIPTION`: add `mvtnorm` to Imports (used by Dunnett T3, never declared); drop `pheatmap` (never used); drop `PMCMRplus` if nothing uses it any more (it needs the GMP/MPFR system libraries); remove the Bioconductor `Additional_repositories` line; fix the placeholder maintainer email and ORCID.
   - Breaking change: bump the version and add a NEWS entry (new names, defaults, `sample_info` + `factors` interface, Type II ANOVA, combined-group post-hoc, `log_input`, removed `write_files`, S4, `chop`, `takeN`, `combine_fc_pvals`, `Adunn.test`, `AdunnettT3Test`).
   - Run `document()` and `R CMD check`.
2. **Documentation**: rewrite the README and both vignettes around the `*_main` functions, including:
   - an **"Approach and limitations"** section: the pre-test reasoning; low power of Shapiro/Levene at 4–6 samples per group (most genes sort to ne); unequal-variance genes often sort to nn (Shapiro on pooled residuals); misrouted heteroscedastic genes can get inflated false positives (ANOVA up to ~17% in our simulations); interaction p-values only for ne genes; post-hoc p-values are within-gene corrected and padj adds BH; Wilcoxon results with ties depend on R version; Dunnett T3 p-values vary slightly between runs unless a seed is set; for RNA-seq interaction analysis, DESeq2/edgeR/limma are the standard tools.
   - a vignette built from the simulations ("How well does statomatic sort genes?"): sorting accuracy and false-positive rates by group size.
3. Open questions (decide or document): BH pooled vs per set; optional `scale_residuals` (helps at ≥6 per group, not at 4–5); optional "always Welch for normal genes"; optional exact (deterministic) Dunnett T3 p-values instead of Monte Carlo.

## Setting up on another computer

```
git clone https://github.com/ruitreves/statomatic.git
cd statomatic
git switch rebuild
```
In R: `install.packages(c("devtools", "testthat", "car", "onewaytests", "PMCMRplus"))`, then `devtools::load_all()` and `devtools::test()`.
