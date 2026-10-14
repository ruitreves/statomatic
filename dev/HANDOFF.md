# statomatic rebuild — handoff

Last updated 2026-10-14. Branch `rebuild`, last commit `2bff938`.

**To resume with Claude Code:** open this repo, run `git switch rebuild`, and tell Claude:
> Read dev/HANDOFF.md and pick up where it leaves off. Save the "Working rules" section to your memory.

---

## Working rules (for Claude)

- **Suggest only.** Never create, edit or delete project files without my explicit OK for that change.
- **Commit messages go in the chat** as copyable text, with the `git add` command. Never write them to a file, never use `git commit -F`. I do the add and commit myself.
- **The original statomatic must keep working exactly as it does now.** Don't touch existing `R/` files, exports or old tests until I say we're integrating. New functions live next to the old ones.
- **Don't run `devtools::document()`.** roxygen2 8.1.0 rewrites `DESCRIPTION` (RoxygenNote) and `man/p.adjustment.methods.Rd` no matter what. Use `devtools::load_all()` and `devtools::test()`.
- If the chat drifts, steer back to the next step below.

## What statomatic does

For each gene (row), test whether the model residuals are normal (Shapiro-Wilk) and the group variances are equal (Levene). That sorts genes into three sets:

| Set | Meaning | 2 groups | >2 groups: omnibus | >2 groups: post-hoc |
|---|---|---|---|---|
| ne | normal, equal variance | t-test | ANOVA | Tukey |
| nu | normal, unequal variance | Welch | Welch | Dunnett T3 |
| nn | not normal | Wilcoxon | Kruskal-Wallis | Dunn |

In the multigroup workflow, each omnibus result is `cbind`-ed with its post-hoc result.

## Conventions in the new code

- New functions are **internal**: roxygen header ending in `@noRd`, no `@export`. Change to `@export` only at integration.
- Names: `new_ttest`, `new_run_<test>`. Files: `R/new_*.R`. Tests: `tests/testthat/test-new_*.R`. Never reuse an existing file name in `R/` (`dunnett_test.R`, `tukey_test.R`, `welch_test.R`, `wilcox_test.R` already exist).
- All results go through `make_res()` (`R/new_utils.R`):
  - named vector (one p-value per gene) → columns `pvalue`, `padj`
  - matrix (one column per comparison/term) → `pvalue_<name>`, `padj_<name>`
  - `padj` = Benjamini-Hochberg across genes; optional, default `TRUE`
  - `sort_columns` pairs each pvalue/padj column (only when `padj = TRUE`)
- Comparison names are `A_vs_B`. Sort them with `order(..., method = "radix")` so the order is the same in every locale.
- `include_test_name` defaults to **FALSE for omnibus tests** (ANOVA, Welch, KW) on purpose: their post-hoc test's `test` column already shows the path. Post-hoc and two-group tests default to TRUE. Exception to remember: in the 2-group workflow Welch has no post-hoc, so `twogroup_main` must call `new_run_welch(..., include_test_name = TRUE)`.
- Post-hoc p-values (Tukey, Dunnett T3, Dunn) are already multiple-comparison corrected within each gene; that is always on. `padj` adds across-gene BH on top. Say this clearly in the docs.
- Old functions use `Pval`/`Padj` and default `padj = FALSE`; new ones use `pvalue`/`padj` and default `padj = TRUE`. Code that reads results by name must be updated at integration.

## Done (all tested against the originals)

| Function | File | Notes |
|---|---|---|
| `make_res` | `R/new_utils.R` | shared result builder |
| `new_ttest`, `new_run_welch`, `new_run_wilcox` | `R/new_two_group.R` | Welch also works for >2 groups |
| `new_run_tukey`, `new_run_dunnett`, `new_run_dunn` | `R/new_multi_group.R` | `sort_columns` defaults TRUE; Dunn requires >2 groups |
| `new_run_anova`, `new_run_kw` | `R/new_anovas.R` | ANOVA is **always Type II** via `car::Anova`; `factors` is a vector (one-way) or a data.frame (any number of factors, names become columns) |
| `new_test_norm` | `R/new_test_norm.R` | returns `ne`/`nu`/`nn` + `tests` (per-gene Shapiro/Levene p, set); `shapiro_threshold`, `levene_threshold`; NAs dropped per gene; constant genes not handled (filter beforehand) |

Tests: `test-new_two_group.R`, `test-new_multi_group.R`, `test-new_anovas.R`, helpers in `helper-new.R` (`with_seed()`, `seeded_dummy_data()`, `quietly()`). Last run: 21 test blocks, 0 failures. Old tests all pass.

Dev folder: `load_data.R` (reads `dev/counts.csv` + `dev/sample_info.csv`, matches samples by name, filters low counts), `new_stuff.R` (scratch), `baseline_tests.txt` (old test results before the rebuild). The data is anonymized and `sample_info.csv` is made up.

## Gotchas found so far

- Old tests use `dummy_data()` with no seed, so they're random. `run_dunn` crashes whenever the random data has exactly 2 groups — pre-existing, left alone.
- Dunnett T3 p-values come from Monte Carlo integration (`mvtnorm::pmvt`): repeated runs differ by up to ~0.0003. Compare old vs new under the same seed.
- `apply()` on one row, or `as.data.frame()` on a vector, changes shapes and names — the cause of earlier `pvalue_res` and one-row `colnames` bugs.
- For 2+ factor ANOVA, Type II intentionally differs from the old `run_anova` (Type I), so only one-way can be compared to the old function.
- `dev/new_stuff.R` may still define old drafts of `new_run_dunn` / `new_run_kw`. Sourcing it after `load_all()` replaces the package versions — delete those definitions.

## Next steps

1. ~~Rebuild `test_norm`~~ **Done** → `R/new_test_norm.R`, tests in `test-new_test_norm.R` (see Done table).
2. Rebuild `fold_change` and `combine_fc_pvals` (`R/fold_change.R`, `R/odds.R`).
3. Rebuild `twogroup_main` and `multigroup_main`. Decide first:
   - **BH pooled or per subset?** Today padj is computed separately within ne, nu and nn. Pooling all omnibus p-values and correcting once is the more common choice.
   - **Two-way designs:** only ne genes get a two-way ANOVA; nu/nn genes get one-way tests on `var3`. Accept and document, or change?
   - **Shapiro power:** with ~8 samples per group many non-normal genes pass as normal. Document it; maybe expose the threshold.
4. ~~Rebuild the S4 layer~~ **Remove the S4 layer at integration** (decided 2026-10-14). Nobody uses it, and it's
   the only reason for the SummarizedExperiment dependency, which `NAMESPACE` imports but `DESCRIPTION`
   never declares. Don't build new versions of it. At integration: delete `class_setup.R`, `methods.R`
   and `sds_analyze.R` (`make_sds`, `sds_*`, `get_*`), and remove `import(SummarizedExperiment)`. Rewrite
   the README and both vignettes around the `*_main` functions. Add one line showing how to use a
   SummarizedExperiment anyway: `new_test_norm(assay(se), as.data.frame(colData(se))[, c(...)])`.
5. Decide whether to keep the modified copies `Adunn.test` / `AdunnettT3Test` or call dunn.test / PMCMRplus directly.
6. Utilities (`chop`, `cf`, `takeN`, `dummy_data`) as needed.
7. Integration: swap `@noRd` → `@export`, replace old functions, update names/defaults, run `document()`.
   Also: drop `pheatmap` from `DESCRIPTION` Imports (listed but never used) and remove the
   Bioconductor `Additional_repositories` line. This is a breaking change: bump the version and add a NEWS entry.

## Setting up on another computer

```
git clone https://github.com/ruitreves/statomatic.git
cd statomatic
git switch rebuild
```
In R: `install.packages(c("devtools", "testthat", "car", "onewaytests", "PMCMRplus"))`, then `devtools::load_all()` and `devtools::test()`.
