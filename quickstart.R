## =============================================================================
##  QUICKSTART.R
##
##  The short tour: one dataset, every estimator, nothing else.  For the full
##  interface with all arguments spelled out, see 1EXAMPLES.R.
## =============================================================================

## -----------------------------------------------------------------------------
##  Copula-based endogeneity corrections in R
##  https://github.com/HashtagHaschka/Copula-based-endogeneity-corrections
##
##  Copyright (C) 2026 Rouven E. Haschka
##  ORCID: https://orcid.org/0000-0002-2916-9745
##
##  If this code contributes to work you publish, please cite the software
##
##    Haschka, R. E. (2026). Copula-based endogeneity corrections in R.
##    https://github.com/HashtagHaschka/Copula-based-endogeneity-corrections
##
##  and, for the estimators it implements, the papers named in each file.
##
##  This program is free software: you can redistribute it and/or modify it
##  under the terms of the GNU General Public License as published by the Free
##  Software Foundation, either version 3 of the License, or (at your option)
##  any later version.
##
##  This program is distributed in the hope that it will be useful, but WITHOUT
##  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
##  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
##  more details.
##
##  You should have received a copy of the GNU General Public License along
##  with this program.  If not, see <https://www.gnu.org/licenses/>.
## -----------------------------------------------------------------------------


gh <- paste0("https://raw.githubusercontent.com/HashtagHaschka/",
             "Copula-based-endogeneity-corrections/functions/")

for (f in c("Copreg_core.R",        # the core has to come first
            "Copreg_pg.R", "Copreg_2scope.R", "Copreg_ima.R",
            "Copreg_jams.R", "Copreg_bmw.R", "Copreg_2scope_np.R",
            "Copreg_panel.R", "Copreg_bayes.R"))
  source(paste0(gh, f))

need <- c("Formula", "np", "bayesm")
miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss) > 0L) install.packages(miss)


## =============================================================================
##  Data: Dominick's orange juice, Tropicana
## =============================================================================

library(bayesm)
data("orangeJuice")

dat <- orangeJuice[[1]]
dat <- dat[dat$brand == 4, ]                       # Tropicana
pc  <- grep("^price", names(dat), value = TRUE)
dat[pc] <- lapply(dat[pc], log)                    # log prices

one <- dat[dat$store == 54, ]                      # one store, 121 weeks

## Tropicana's own price is endogenous: it responds to demand shocks that also
## sit in the error.  Competing prices and the promotion variables are the
## controls.  Position around the "|" is what decides which is which.
f <- logmove ~ price4 | price1 + price2 + price3 + feat


## =============================================================================
##  Cross-section: one store
## =============================================================================

ols <- lm(logmove ~ price4 + price1 + price2 + price3 + feat, data = one)

## Each estimator defaults to the CDF its own paper uses; cdf = overrides that.
pg   <- CopRegPG(f, one)
sc   <- CopReg2sCOPE(f, one)
ima  <- CopRegIMA(f, one)
bmw  <- CopRegBMW(f, one)
jams <- CopRegJAMS(f, one)

## 2sCOPE-np uses the conditional CDF of price4 given the controls, estimated
## by kernel, so it has no cdf argument at all.  It is much slower than the
## others: the bandwidths are cross-validated in every bootstrap replicate.
np <- CopReg2sCOPEnp(f, one, verbose = TRUE)

## BAYES estimates the marginal CDFs jointly with everything else and returns a
## posterior rather than a bootstrap distribution.  iterations, burnin and thin
## take the place of nboots.
bay <- CopRegBAYES(f, one, verbose = TRUE)

summary(pg)
summary(sc)
summary(ima)
summary(bmw)
summary(jams)
summary(np)
summary(bay)

## validity() walks the identification requirements and says what the data
## show for each.  What it checks differs: BMW needs its first-stage residuals
## to be nonnormal rather than price4, 2sCOPE-np needs its Assumption 3, and
## for BAYES the convergence of the chain comes with it.
validity(pg)
validity(sc)
validity(ima)
validity(bmw)
validity(jams)
validity(np)
validity(bay)

## The correction should move the price coefficient away from OLS, and the
## standard errors should not explode while doing it.  pg$icon reports that
## inflation; a factor beyond five or six means the model is barely identified.
round(c(OLS   = coef(ols)["price4"],
        PG    = coef(pg)["price4"],
        s2COPE = coef(sc)["price4"],
        IMA   = coef(ima)["price4"],
        BMW   = coef(bmw)["price4"],
        JAMS  = coef(jams)["price4"],
        np    = coef(np)["price4"],
        BAYES = coef(bay)["price4"]), 3)

## the uncorrected benchmark from the same bootstrap resamples, so the standard
## errors are comparable rather than merely printed side by side
pg$ols.coefficients["price4"]
pg$ols.std.error["price4"]
pg$icon

## the chain behind the Bayesian fit
plot(bay, type = "trace")
plot(bay, type = "cdf")                            # the estimated margin of price4


## =============================================================================
##  Panel: all 83 stores
## =============================================================================
##
##  Same formula, all stores, week as the time index.  Data go in untransformed;
##  the function does the fixed-effects transformation itself and bootstraps
##  whole panels rather than rows.

pan <- CopRegPANEL(f, dat, index = c("store", "week"))

summary(pan)
validity(pan)

## structural = FALSE gives the residuals of the transformed model, which are
## the ones the likelihood assumes to be normal
qqnorm(residuals(pan, structural = FALSE))
qqline(residuals(pan, structural = FALSE))

fixef(pan)                                         # the store effects
round(c(within = pan$within.coefficients["price4"],
        PANEL  = coef(pan)["price4"]), 3)
