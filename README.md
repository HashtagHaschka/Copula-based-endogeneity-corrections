The functions implement copula-based endogeneity corrections using the least-squares-based correction function approach.

PG is the estimator by Park & Gupta (2012)
2sCOPE is the estimator by Yang et al. (2024)
IMA is the estimator by Haschka (2024)
BWM is the estimator by Breitung et al. (2024)

All functions are similar in terms of required arguments. In light of the current discussion in the literature on the estimation of the cumulative distribution function, 
all approaches considered in the literature are implemented and can be selected via cdf = c("kde", "ecdf", "resc.ecdf", "adj.ecdf")
- kde is the integral of a density estimator used in Park & Gupta (2012) and Haschka (2022)
- ecdf is the empirical cumulative distribution function (ecdf) with replaced boundary proposed by Becker et al. (2022)
- resc.ecdf is a rescaled ecdf proposed by Qian et al. (2024)
- adj.ecdf is an adjusted ecdf proposed by Liengaard (2024)


REFERENCES

- Becker, J.-M., D. Proksch, and C. M. Ringle (2021). Revisiting Gaussian copulas to handle endogenous regressors. Journal of the Academy of Marketing Science 50, 46–66.
- Breitung, J., A. Mayer, and D. Wied (2024). Asymptotic properties of endogeneity corrections using nonlinear transformations. The Econometrics Journal 27 (3), 362–383.
- Haschka, R. E (2022). Handling endogenous regressors using copulas: A generalisation to linear panel models with fixed effects and correlated regressors. Journal of Marketing Research 59(4), 860–881.
- Haschka, R. E. (2024). Robustness of copula-correction models in causal analysis: Exploiting between-regressor correlation. IMA Journal of Management Mathematics 36 (1), 161–180.
- Liengaard, B. D., J.-M. Becker, M. Bennedsen, P. Heiler, L. N. Taylor, and C. M. Ringle (2024). Dealing with regression models’ endogeneity by means of an adjusted estimator for the Gaussian copula approach. Journal of the Academy of Marketing Science, 1–21.
- Park, S. and S. Gupta (2012). Handling endogenous regressors by joint estimation using copulas. Marketing Science 31 (4), 567–586.
- Qian, Y., A. Koschmann, and H. Xie (2024). A practical guide to endogeneity correction using copulas. NBER Working Paper. https://www.nber.org/system/files/workingpapers/w32231/w32231.pdf. 
- Yang, F., Y. Qian, and H. Xie (2024). Addressing endogeneity using a two-stage copula generated regressor approach. EXPRESS: Journal of Marketing Research
