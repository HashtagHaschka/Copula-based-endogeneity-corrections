# Copula-based endogeneity corrections in R

Instrument-free corrections for endogenous regressors in linear models, using
the Gaussian copula to model the dependence between the regressors and the
structural error. Eight estimators, one interface.

Estimation is by augmented least squares (PG, 2sCOPE, 2sCOPE-np, IMA, BMW,
JAMS), by maximum likelihood (PANEL), or by MCMC (BAYES).

| Function | Estimator | Reference |
|---|---|---|
| `CopRegPG()` | PG | Park & Gupta (2012) |
| `CopReg2sCOPE()` | 2sCOPE | Yang, Qian & Xie (2025) |
| `CopReg2sCOPEnp()` | 2sCOPE-np | Hu, Qian & Xie (2025) |
| `CopRegIMA()` | IMA | Haschka (2025a) |
| `CopRegBMW()` | BMW | Breitung, Mayer & Wied (2024) |
| `CopRegJAMS()` | JAMS | Liengaard et al. (2025) |
| `CopRegPANEL()` | PANEL | Haschka (2022) |
| `CopRegBAYES()` | BAYES | Haschka (2025b) |

---

## Getting started

Download **`1EXAMPLES.R`** and run it. That is the only file you need: it reads
the estimators straight from this repository.

```r
gh <- paste0("https://raw.githubusercontent.com/HashtagHaschka/",
             "Copula-based-endogeneity-corrections/main/")

for (f in c("Copreg_core.R",        # the core has to come first
            "Copreg_pg.R", "Copreg_2scope.R", "Copreg_ima.R",
            "Copreg_jams.R", "Copreg_bmw.R", "Copreg_2scope_np.R",
            "Copreg_panel.R", "Copreg_bayes.R"))
  source(paste0(gh, f))
```

`Copreg_core.R` holds the shared machinery — formula handling, the CDF
estimators, the bootstrap, the diagnostics and the S3 methods — and must be
sourced before the others.

**Requirements.** `Formula` for every estimator, `np` in addition for
2sCOPE-np. Nothing else: the Bayesian sampler is written in base R. The
examples also use `bayesm`, `AER` and `ISLR` for their data.

---

## The formula

All estimators share the same two-part syntax:

```
y ~ endogenous_1 + endogenous_2 + ... | exogenous_1 + exogenous_2 + ...
```

**Position decides.** A term written before the `|` is treated as endogenous
and receives a copula correction; written after it, the same term is exogenous
and receives none. Everything `lm()` accepts works: transformations
(`log(price)`), interactions (`price:feat`), polynomials (`I(x^2)`), factors
(`as.factor(store)` or a factor column), and `-1` to drop the intercept.
Missing values are handled by `na.omit`.

Writing an interaction or a power of an endogenous regressor before the `|`
gives it its own copula term and produces a warning: such terms are strongly
collinear with the base term, and Qian, Koschmann & Xie (2025) find the added
variance usually outweighs the gain. It is a warning, not a veto.

---

## Arguments shared by the cross-sectional estimators

`CopRegPG()`, `CopReg2sCOPE()`, `CopRegIMA()`, `CopRegBMW()` and
`CopRegJAMS()` take the same arguments.

| Argument | Default | Meaning |
|---|---|---|
| `formula` | — | two-part formula, see above |
| `data` | — | a `data.frame` |
| `cdf` | per estimator | how the marginal CDF is estimated, see below |
| `ties` | `"max"` | `"max"` for the counting function, `"average"` for midranks |
| `nboots` | `199` | bootstrap replicates for the standard errors |
| `subset` | `NULL` | logical vector, as in `lm()` |
| `contrasts` | `NULL` | passed to `model.matrix()` |
| `parallel` | `FALSE` | `TRUE`, `"multicore"` or `"snow"` to spread the bootstrap |
| `ncores` | `NULL` | workers; `NULL` uses one less than detected |
| `verbose` | `interactive()` | progress bar during the bootstrap |

No seed is set anywhere. Set one yourself beforehand if you want reproducible
standard errors.

### `cdf` — estimating the marginal distribution

The literature disagrees on how to estimate the CDF that enters the copula
transformation, so all proposals are implemented and the choice is free
(except in BMW). Seven options:

| Value | What it is |
|---|---|
| `"kde.silverman"` | integral of a Gaussian kernel density, Silverman bandwidth; Park & Gupta (2012) |
| `"kde.cv"` | same, cross-validated bandwidth; Li, Li & Racine (2017) |
| `"kde.plugin"` | same, plug-in bandwidth; Polansky & Baker (2000) |
| `"ecdf.fixed"` | empirical CDF with replaced boundary; Becker, Proksch & Ringle (2022) |
| `"ecdf.adj"` | adjusted ECDF; Liengaard et al. (2025) |
| `"rank.n"` | rescaled ECDF, Rank/n with a correction; Qian, Koschmann & Xie (2025) |
| `"rank.n1"` | Rank/(n+1); Breitung, Mayer & Wied (2024) |

Defaults follow each estimator's own paper: `"kde.silverman"` for PG,
`"rank.n"` for 2sCOPE and IMA, `"ecdf.adj"` for JAMS, `"rank.n1"` for BMW,
`"kde.plugin"` for PANEL. BMW is the one place where the choice is not free:
its Proposition 3.1 is derived for Rank/(n+1), so anything else warns.

---

## The estimators

### `CopRegPG()` — Park & Gupta (2012)

The original. Transforms each endogenous regressor to a normal score,
C = Φ⁻¹(F̂(P)), and adds it to the regression. Assumes the endogenous and
exogenous regressors are uncorrelated; `summary()` and `validity()` test that
assumption, because violating it is what motivated the later estimators.

### `CopReg2sCOPE()` — Yang, Qian & Xie (2025)

Relaxes that assumption. Transforms the exogenous regressors as well, runs a
first-stage regression of P\* on W\* and uses its residual as the copula term,
so the correction is orthogonal to W by construction.

### `CopRegIMA()` — Haschka (2025a)

2sCOPE without an intercept in the first stage, which is what the paper's
derivation of ρ requires. The two are usually very close; the gap widens the
further the transformed exogenous regressors sit from mean zero.

### `CopRegBMW()` — Breitung, Mayer & Wied (2024)

Inverts the order: the first stage runs on the raw variables and the rank
transform is applied to its **residuals**. Two consequences show in the
output. The identification requirement falls on the first-stage residuals
rather than on P, so `validity()` tests those; and their Corollary 3.2 makes
the textbook *t* statistic valid for testing ρ = 0, so `summary()` reports a
Durbin–Hausman–Wu test next to the bootstrap one.

### `CopRegJAMS()` — Liengaard et al. (2025)

Lets the copula structure differ across the categories of the discrete
exogenous regressors.

| Argument | Default | Meaning |
|---|---|---|
| `conditional` | `TRUE` | `TRUE` uses the joint categories of all factors in the exogenous part; `FALSE` imposes one common structure; a character vector names the variables to condition on |

`summary()` adds a bootstrap Wald test of the copula terms jointly, and a
second asking whether they differ across the cells.

### `CopReg2sCOPEnp()` — Hu, Qian & Xie (2025)

Replaces the whole first stage by a nonparametric estimate of the
**conditional** CDF of P given W, C = Φ⁻¹(F̂(P|W)), computed with the kernel
estimator of Li & Racine (2008) as implemented in `np::npcdist()`. It assumes
nothing about the marginals or the functional form linking the regressors, and
it is the only estimator here that tolerates **discrete** endogenous regressors:
the kernel smooths P as well, so the transformation stays well defined for
counts and even binary treatments.

Has neither `cdf` nor `ties` — there is no marginal CDF to choose — and
supplying either is an error rather than a silently ignored argument.

| Argument | Default | Meaning |
|---|---|---|
| `condition` | `NULL` | which exogenous variables the conditional CDF conditions on; `NULL` uses all of them except polynomial and interaction terms, which the kernel handles by itself |
| `groups` | `NULL` | estimate the conditional CDF separately within each joint category of the named variables instead of smoothing them |
| `heterogeneous` | `NULL` | let the endogeneity itself differ across groups; the copula terms are parameterised against the first group, so the remaining coefficients test the differences |
| `demean` | `NULL` | treat the named categories as fixed effects and sweep them out of P and the continuous conditioning variables first |
| `bwmethod` | `"cv.ls"` | cross-validated bandwidths; `"normal-reference"` is the cheap alternative |
| `bw.refit` | `TRUE` | reselect the bandwidths in every bootstrap replicate, which is what the paper prescribes; `FALSE` freezes them |
| `...` | — | passed to `np::npcdistbw()`, e.g. `ckertype`, `ukertype` |

**This estimator is slow.** Cross-validation is O(n²) per evaluation and is
repeated in every bootstrap replicate. With `verbose = TRUE` the estimator
prints the expected runtime before it starts.

### `CopRegPANEL()` — Haschka (2022)

Linear panel model with fixed effects, estimated by maximum likelihood after a
forward orthogonal deviations transformation. A separate model class, the way
`plm` is separate from `lm`; it is deliberately not reachable through
`copreg()`.

| Argument | Default | Meaning |
|---|---|---|
| `index` | — | `c("panelvariable", "timevariable")`, id first |
| `intercept` | `FALSE` | a free constant in the transformed regression; the structural intercept is absorbed by the individual effects |
| `method` | `"BFGS"` | optimiser; falls back to Nelder–Mead by itself |
| `start` | `NULL` | starting values |
| `maxit` | `5000` | iteration cap |

Data go in untransformed; the function transforms them. `lag()`, `lead()` and
`diff()` work inside the formula and never reach across a panel boundary. Time
dummies stay in the model but are dropped from the copula, since after the
transformation a full set of them is one dimension short. The bootstrap
resamples whole panels, not rows. `fixef()`, `logLik()`, `AIC()` and `BIC()`
exist for this class.

### `CopRegBAYES()` — Haschka (2025b)

Everything the others do in stages happens here in one. The marginal CDFs are
not plugged in but are unknown parameters — probability masses on the uniquely
observed values of each regressor — drawn along with the coefficients, the
error variance and the full copula correlation matrix. Nothing is estimated a
priori, so there is no plug-in uncertainty to bootstrap back in; the posterior
carries it. Ties are the normal case rather than a problem, and a binary
regressor simply has two masses.

Has no `cdf`, no `ties` and no `nboots`; supplying any of them is an error.

| Argument | Default | Meaning |
|---|---|---|
| `iterations` | `102000` | MCMC iterations |
| `burnin` | `2000` | discarded at the start |
| `thin` | `100` | keep every *n*-th draw; thinning happens while sampling, so only the retained draws are held in memory |
| `horseshoe` | `TRUE` | horseshoe prior on every coefficient including the intercept; `FALSE` uses an independent normal with an inverse Gamma variance |
| `prior.args` | `list()` | overrides any hyperparameter, all defaulting to the paper |
| `start` | `NULL` | starting values; the default is OLS, the identity, and one Dir(1) draw per regressor |

Exogeneity is imposed by holding the corresponding entries of the copula
covariance matrix at zero. That is done by a block-recursive reparameterisation
rather than by an inverse Wishart draw, so the zeros hold exactly rather than
approximately.

`plot()` shows the chain — `type = "trace"`, `"acf"`, `"pacf"`, `"density"` —
with `which =` selecting parameters by name or position. `type = "cdf"` shows
the estimated marginal distribution of a regressor with a pointwise credible
band, an object that does not exist in the frequentist estimators.

### `copreg()`

Generic entry point: `copreg(formula, data, method = "pg", ...)` with `method`
one of `"pg"`, `"2scope"`, `"ima"`, `"bmw"`, `"jams"`, `"np"`. It only rewrites
the call and hands it to the estimator's own function, so both routes behave
identically down to the warnings. PANEL and BAYES are deliberately outside:
the first needs an index and a different likelihood, the second replaces
bootstrap inference by a posterior.

---

## Working with a fitted model

```r
coef(mod); vcov(mod); nobs(mod); formula(mod)
confint(mod, level = 0.95, type = "normal")      # or "percentile"
residuals(mod, type = "structural")              # or "augmented"
fitted(mod); predict(mod, newdata = ...)         # structural model only
summary(mod)
update(mod, . ~ . | . - x3)                      # or method = "jams"
```

Copula terms never enter `predict()` or `fitted()`: they are endogeneity
controls, not part of the causal model.

Each object also carries the uncorrected benchmark from the same bootstrap
resamples (`$ols.coefficients`, `$ols.std.error`, or `$within.coefficients`
for PANEL), the raw bootstrap draws (`$boot`), the diagnostics
(`$diagnostics`), and the standard error inflation `$icon`.

For BAYES, `coef()` is the posterior mean, `vcov()` the posterior covariance
and `confint()` a credible interval taken from the draws themselves. The draws
are in `$draws`, `$Sigma.draws` and `$lambda.draws`.

### `validity()`

Walks the identification requirements and reports what the data say about
each: the decision tree of Yang, Qian & Xie (2025, Fig. 2), the sample-size
and nonnormality boundaries of Becker, Proksch & Ringle (2022, Fig. 8), and
the standard error inflation relative to uncorrected OLS.

```r
validity(mod, level = 0.05, power = 0.8)
```

What it checks depends on the estimator. For BMW the nonnormality requirement
falls on the first-stage residuals; for 2sCOPE-np it is Assumption 3 rather
than nonnormality of P; for BAYES the convergence of the chain comes with it —
Geweke, effective sample size and autocorrelation, plus Gelman–Rubin on
request via `chains = TRUE`, which costs one further run of the sampler per
chain.

---

## References

Arellano, M. (1993). On the testing of correlated effects with panel data.
*Journal of Econometrics* 59, 87–97.

Becker, J.-M., D. Proksch, and C. M. Ringle (2022). Revisiting Gaussian
copulas to handle endogenous regressors. *Journal of the Academy of Marketing
Science* 50, 46–66.

Breitung, J., A. Mayer, and D. Wied (2024). Asymptotic properties of
endogeneity corrections using nonlinear transformations. *The Econometrics
Journal* 27(3), 362–383.

Carvalho, C. M., N. G. Polson, and J. G. Scott (2010). The horseshoe estimator
for sparse signals. *Biometrika* 97(2), 465–480.

Gamerman, D. (1997). Sampling from the posterior distribution in generalized
linear mixed models. *Statistics and Computing* 7(1), 57–68.

Gonçalves, S. and L. Kilian (2004). Bootstrapping autoregressions with
conditional heteroskedasticity of unknown form. *Journal of Econometrics* 123,
89–120.

Haschka, R. E. (2022). Handling endogenous regressors using copulas: A
generalization to linear panel models with fixed effects and correlated
regressors. *Journal of Marketing Research* 59(4), 861–880.

Haschka, R. E. (2025a). Robustness of copula-correction models in causal
analysis: Exploiting between-regressor correlation. *IMA Journal of Management
Mathematics* 36(1), 161–180.

Haschka, R. E. (2025b). Bayesian inference for joint estimation models using
copulas to handle endogenous regressors. *Oxford Bulletin of Economics and
Statistics*.

Hayfield, T. and J. S. Racine (2008). Nonparametric econometrics: The np
package. *Journal of Statistical Software* 27(5), 1–32.

Hu, X., Y. Qian, and H. Xie (2025). Correcting endogeneity via nonparametric
copula control functions. NBER Working Paper 33607.
<http://www.nber.org/papers/w33607>

Li, Q. and J. S. Racine (2008). Nonparametric estimation of conditional CDF
and quantile functions with mixed categorical and continuous data. *Journal of
Business & Economic Statistics* 26(4), 423–434.

Li, Q., J. Li, and J. S. Racine (2017). Cross-validated mixed-datatype
bandwidth selection for nonparametric cumulative distribution/survivor
functions. *Econometric Reviews* 36, 970–987.

Li, Q., J. Lin, and J. S. Racine (2013). Optimal bandwidth selection for
nonparametric conditional distribution and quantile functions. *Journal of
Business & Economic Statistics* 31(1), 57–65.

Liengaard, B. D., J.-M. Becker, M. Bennedsen, P. Heiler, L. N. Taylor, and
C. M. Ringle (2025). Dealing with regression models' endogeneity by means of
an adjusted estimator for the Gaussian copula approach. *Journal of the
Academy of Marketing Science* 53, 279–299.

Ng, K. W., G. Tian, and M. Tang (2011). *Dirichlet and Related Distributions*.
Chichester: Wiley.

Park, S. and S. Gupta (2012). Handling endogenous regressors by joint
estimation using copulas. *Marketing Science* 31(4), 567–586.

Polansky, A. M. and E. R. Baker (2000). Multistage plug-in bandwidth selection
for kernel distribution function estimates. *Journal of Statistical
Computation and Simulation* 65, 63–80.

Qian, Y., A. Koschmann, and H. Xie (2025). A practical guide to endogeneity
correction using copulas. *Journal of Marketing*.

Trivedi, P. K. and D. M. Zimmer (2017). A note on identification of bivariate
copulas for discrete count data. *Econometrics* 5(1), 10.

Yang, F., Y. Qian, and H. Xie (2025). Addressing endogeneity using a two-stage
copula generated regressor approach. *Journal of Marketing Research* 62(4),
601–623.

Zhang, X., W. J. Boscardin, and T. R. Belin (2006). Sampling correlation
matrices in Bayesian models with correlated latent variables. *Journal of
Computational and Graphical Statistics* 15(4), 880–896.
