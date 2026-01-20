The functions implement copula-based endogeneity corrections using the least-squares-based correction function approach, maximum likelihood (PANEL), or MCMC sampling (BAYES):
- PG is the estimator by Park & Gupta (2012)
- 2sCOPE is the estimator by Yang et al. (2025)
- 2sCOPE-np is the estimator by Hu et al. (2025)
- IMA is the estimator by Haschka (2024)
- BWM is the estimator by Breitung et al. (2024)
- JAMS is the estimator by Liengaard et al. (2025)
- PANEL is the estimator by Haschka (2022)
- BAYES is the estimator by Haschka (2025)

The functions for cross-sectional estimators are similar in terms of required arguments. In light of the current discussion in the literature on the estimation of the cumulative distribution function, 
all approaches considered in the literature are implemented and can be selected via cdf = c("kde", "ecdf", "resc.ecdf", "adj.ecdf")
- kde is the integral of a density estimator used in Park & Gupta (2012) and Haschka (2022)
- ecdf is the empirical cumulative distribution function (ecdf) with replaced boundary proposed by Becker et al. (2022)
- resc.ecdf is a rescaled ecdf proposed by Qian et al. (2025), used in Haschka (2024) and Yang et al. (2024)
- adj.ecdf is an adjusted ecdf proposed by Liengaard et al. (2024)

The required arguments should be specified as follows:
- formula should be depvar ~ endog_var1 + endog_var2 + ... | exog_var1 + exog_var2 + ...
- depvar is the dependent variable, endog_var1, etc., are the continuous endogenous variables, exog_var1, etc., are the exogenous variables
- formula accepts -1 to remove intercept
- dummy variables can be modelled using as.factor(exog_var1), etc.
- data should be a data.frame
- Example: CopRegJAMS(formula = Y ~ endog1 + endog2 + endog3 | exog1 + exog2 + as.factor(exog3), cdf = "adj.ecdf", data = data1)

For the copula panel regression model:
- dummy variables must be time-varying because it is a fixed-effects model
- index represents the panel structure and is specified by index = c("panelvariable", "timevariable")
- the panel id variable should come first in the argument, as the model transformation and bootstrapping is done within each panel
- panelvariable and timevariable must be numeric
- data should be untransformed panel data, as the function will perform model transformation automatically
- function will return error messages if the identifying assumptions that can be tested are not met. This includes normality of endogenous regressors, symmetry of residual distribution, and small differences between the residual distribution and the distribution of endogenous regressor(s).
- bootstrapping is done in parallel
- Example: CopRegML_par(formula = Y ~ endog1 + endog2 + endog3 | exog1 + exog2 + as.factor(year), index = c("id", "year"), data = data1)


REFERENCES
- Becker, J.-M., D. Proksch, and C. M. Ringle (2021). Revisiting Gaussian copulas to handle endogenous regressors. Journal of the Academy of Marketing Science 50, 46–66.
- Breitung, J., A. Mayer, and D. Wied (2024). Asymptotic properties of endogeneity corrections using nonlinear transformations. The Econometrics Journal 27 (3), 362–383.
- Haschka, R. E (2022). Handling endogenous regressors using copulas: A generalisation to linear panel models with fixed effects and correlated regressors. Journal of Marketing Research 59(4), 860–881.
- Haschka, R. E. (2024). Robustness of copula-correction models in causal analysis: Exploiting between-regressor correlation. IMA Journal of Management Mathematics 36 (1), 161–180.
- Haschka, R. E. (2025). Bayesian Inference for Joint Estimation Models Using Copulas to Handle Endogenous Regressors. Oxford Bulletin of Economics and Statistics 0, 1–16.
- Hu, X., Qian, Y., A., and H. Xie (2025). Correcting endogeneity via instrument-free two-stage nonparametric copula control functions. http://www.nber.org/papers/w33607
- Liengaard, B. D., J.-M. Becker, M. Bennedsen, P. Heiler, L. N. Taylor, and C. M. Ringle (2025). Dealing with regression models’ endogeneity by means of an adjusted estimator for the Gaussian copula approach. Journal of the Academy of Marketing Science 53, 279–299.
- Park, S. and S. Gupta (2012). Handling endogenous regressors by joint estimation using copulas. Marketing Science 31 (4), 567–586.
- Qian, Y., A. Koschmann, and H. Xie (2025). A practical guide to endogeneity correction using copulas. Journal of Marketing.
- Yang, F., Y. Qian, and H. Xie (2025). Addressing endogeneity using a two-stage copula generated regressor approach. Journal of Marketing Research 62(4), 601-623.
