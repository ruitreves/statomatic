# statomatic rebuild — handoff

Last updated 2026-10-19. Branch `rebuild`, last commit `84bbee6`. **Integration is done**: statomatic 1.0.0 passes `R CMD check` with 0 errors, 0 warnings, 0 notes. Next: documentation, then merge `rebuild` into `main`.

**To resume with Claude Code:** open this repo, run `git pull` on `rebuild`, and tell Claude:
> Read dev/HANDOFF.md and pick up where it leaves off. Save the "Working rules" and "Direction" sections to your memory.

---

## Working rules (for Claude)

- **Suggest only.** Never create, edit or delete project files without my explicit OK for that change.
- **Commit messages go in the chat** as copyable text, with the `git add` command. Never write them to a file, never use `git commit -F`. I do the add and commit myself.
- **Vignettes and user-facing docs get extra review by me.** Show the rendered result; call `head()` on whole data.frames (never crop columns), for inputs and outputs.
- If the chat drifts, steer back to the next step below.
- I work on two computers: `git pull` at the start of a session, `git push` at the end.

## Direction (decided 2026-10-15)

**A faithful rebuild of statomatic** that keeps its pre-test sorting approach, with updated documentation of the approach and its limitations. statomatic is my first first-author publication (Treves, Gripshover & Hardesty 2025, *Stats* 8(1), doi:10.3390/stats8010018); don't re-argue replacing it with DESeq2/edgeR/limma. "Faithful" still allows improvements: Type II ANOVA, name-based sample matching, the `tests` table, `log_input`, design checks.

## What statomatic does

For each row (gene, taxon, ...), test whether the model residuals are normal (Shapiro-Wilk) and the group variances are equal (Levene). That sorts rows into three sets:

| Set | Meaning | 2 groups | >2 groups: omnibus | >2 groups: post-hoc |
|---|---|---|---|---|
| ne | normal, equal variance | t-test | ANOVA | Tukey |
| nu | normal, unequal variance | Welch | Welch | Dunnett T3 |
| nn | not normal | Wilcoxon | Kruskal-Wallis | Dunn |

In the multigroup workflow, each omnibus result is combined with its post-hoc result and the fold changes.

## The package (1.0.0)

| Exported function | File | Notes |
|---|---|---|
| `multigroup_main(x, sample_info, factors, ...)` | `R/multigroup_main.R` | 1–3 factors; one-way results identical to the original |
| `twogroup_main(x, sample_info, factors, ...)` | `R/twogroup_main.R` | results identical to the original; `include_test_name = TRUE` by default |
| `test_norm(x, factors)` | `R/test_norm.R` | returns `ne`/`nu`/`nn` + `tests`; separate thresholds; NAs dropped per row; constant rows not handled (filter beforehand) |
| `run_ttest`, `run_welch`, `run_wilcox` | `R/two_group_tests.R` | Welch also works for >2 groups; Wilcoxon suppresses the ties warning |
| `run_tukey`, `run_dunnett`, `run_dunn` | `R/multi_group_tests.R` | Dunn needs >2 groups |
| `run_anova(x, factors)`, `run_kruskal` | `R/anovas.R` | ANOVA is **always Type II** via `car::Anova` |
| `fold_change(x, var1, log_input)` | `R/fold_change.R` | log2 of ratio of means; `log_input = TRUE` → difference of means |
| `cf`, `dummy_data` | `R/odds.R` | kept from the original |

Internal helpers: `make_res`, `select_factors`, `check_factors`, `add_fold_changes`, `n_rows` (`R/utils.R`); `dunn_bonferroni`, `dunnett_t3` (`R/multi_group_tests.R`) — statomatic's own implementations, identical to the modified copies of dunn.test 1.3.6 and PMCMRplus::dunnettT3Test that 1.0.0 removed.

Imports: `car`, `mvtnorm`, `onewaytests`, `stats`, `utils`. Maintainer: rui.treves@louisville.edu (no ORCID listed). `NEWS.md` has the full 1.0.0 entry.

## Conventions

- Results go through `make_res()`: named vector → `pvalue`, `padj`; matrix → `pvalue_<name>`, `padj_<name>`. `padj` = Benjamini-Hochberg across rows, **within each set**, default `TRUE`.
- **Main functions** take the sample table and factor column names. `select_factors()` finds the sample IDs (`sample_col`, else row names, else the one column containing every sample of `x`), ignores extra rows, and errors clearly on typos, ambiguous/missing IDs, or a column passed instead of its name. `check_factors()` validates the design (1–3 factors, fully crossed, ≥2 per group, no NAs; 2-group mode needs one factor with 2 levels) and builds the combined group (`condition_sex`) with **radix-sorted levels**. ANOVA and `test_norm` use the factors; Welch, KW, post-hoc tests and fold changes use the combined group.
- **Individual test functions** take group labels as a plain vector in the column order of `x` — no name matching (documented in the additional workflows vignette).
- **Comparison names** are `A_vs_B` in the group's level order everywhere; fold changes are attached by name (`add_fold_changes()`).
- Result layout: omnibus `pvalue`/`padj` (per ANOVA term for multi-way), then per comparison `pvalue`, `padj`, `log2FoldChange`, then `<group>_mean`, then `test`.
- `include_test_name` defaults to FALSE for omnibus tests (the post-hoc `test` column shows the path), TRUE for post-hoc and two-group tests.
- Post-hoc p-values are within-row multiple-comparison corrected; `padj` adds across-row BH.
- Progress messages say "rows" (`verbose = TRUE`).

## Tests

One `test-*.R` per R file, helpers in `tests/testthat/helper.R` (`with_seed()`, `seeded_dummy_data()`, `quietly()`, `original()`). The faithfulness tests compare against **saved results of the original functions** in `tests/testthat/fixtures/original_results.rds` (338 results from commit `5673d8e`, R 4.6.1). `original("key", expr)` returns the saved result; `expr` is only run when recording (`dev/record_original_results.R`, which needs a commit that still has the original code, e.g. `5673d8e`). If a test's input data changes, its saved results must be re-recorded from that commit. Last run: 74 test blocks, 0 failures.

Dev folder (excluded from the build): `load_data.R`, `new_stuff.R` (scratch; still uses the old `new_*` names), `record_original_results.R`, `baseline_tests.txt`, `counts.csv` / `sample_info.csv` (anonymized; sample_info is made up).

## Gotchas

- Dunnett T3 p-values come from Monte Carlo integration (`mvtnorm::pmvt`): repeated runs differ by up to ~0.0003 unless a seed is set.
- **R 4.6.0 changed `wilcox.test`:** with ties it computes exact conditional p-values; before 4.6 it fell back to the normal approximation. Wilcoxon p-values for tied data differ between R versions.
- **dunn.test 1.4.1 gives different p-values than the 1.3.6 copy statomatic used**; `dunn_bonferroni` reproduces 1.3.6. Unmodified `PMCMRplus::dunnettT3Test` errors on two identical constant groups; `dunnett_t3` keeps the fix (p = 1).
- `apply()` on one row, or `as.data.frame()` on a vector, changes shapes and names.
- For 2+ factors, Type II ANOVA and combined-group post-hoc tests intentionally differ from the original, so only one-way designs are compared to it.
- `dummy_data()` is unseeded; tests use `seeded_dummy_data()`.

## Next steps

1. **Documentation**:
   - Rewrite `README.md` around the main functions (it still shows the S4 interface; it uses the top-level `normalized_counts.csv`, excluded from the build).
   - An **"Approach and limitations"** section (README and/or its own vignette): the pre-test reasoning; low power of Shapiro/Levene at 4–6 samples per group (most rows sort to ne); unequal-variance rows often sort to nn (Shapiro on pooled residuals); misrouted heteroscedastic rows can get inflated false positives (ANOVA up to ~17% in our simulations); interaction p-values only for ne rows; post-hoc p-values are within-row corrected and padj adds BH within each set; Wilcoxon results with ties depend on R version; Dunnett T3 needs a seed for reproducibility; for RNA-seq interaction analysis, DESeq2/edgeR/limma are the standard tools.
   - A vignette built from the simulations ("How well does statomatic sort genes?"): sorting accuracy and false-positive rates by group size.
   - Rebuild the pkgdown site (`docs/`, `_pkgdown.yml`): it still documents the old functions.
   - Optional: `@examples` in the help pages.
2. **Merge** `rebuild` into `main` with a pull request; tag `v1.0.0`.
3. Open questions (decide or document): BH pooled vs per set; a `combine_results()` helper (deferred); optional `scale_residuals` (helps at ≥6 per group, not at 4–5); optional "always Welch for normal rows"; optional exact (deterministic) Dunnett T3 p-values.

Done vignettes: `basic_workflow.Rmd` ("Basic workflow") and `additional_workflows.Rmd` ("Additional usage of statomatic"), both rewritten for 1.0.0.

## Setting up on another computer

```
git clone https://github.com/ruitreves/statomatic.git
cd statomatic
git switch rebuild
```
In R: `install.packages(c("devtools", "testthat", "car", "mvtnorm", "onewaytests", "rmarkdown", "knitr"))`, then `devtools::load_all()` and `devtools::test()`.
