# Chapter 3 — Maximum Likelihood Estimation and Statistical Inference

*Discrete Choice Modelling: Theory, Mathematics, and Behavioural Interpretation*

## Learning goals

By the end of this chapter, the reader should be able to:

1. Distinguish a probability model, a likelihood function, and a log-likelihood function.
2. Derive the multinomial logit (MNL) log-likelihood from observed choices.
3. Explain how maximum likelihood estimates are found and what identification requires.
4. Interpret the score, Hessian, information matrix, and standard errors.
5. Construct confidence intervals and perform Wald and likelihood-ratio tests.
6. Understand why repeated observations call for cluster-robust inference.
7. Evaluate model fit without confusing statistical significance with behavioural importance.

---

## 3.1 From choice probabilities to parameter estimation

A discrete choice model assigns a probability to each alternative in a decision-maker's choice set. For individual or choice occasion \(n\), let \(\mathcal C_n\) be the set of available alternatives and \(y_n\) the alternative actually selected.

The random utility specification is

$$
U_{ni}=V_{ni}(\boldsymbol\theta)+\varepsilon_{ni}, \qquad i\in\mathcal C_n,
$$

where \(\boldsymbol\theta\) is the vector of unknown parameters. Under the usual MNL assumptions, the choice probability is

$$
P_{ni}(\boldsymbol\theta)=
\frac{\exp(V_{ni}(\boldsymbol\theta))}
{\displaystyle\sum_{j\in\mathcal C_n}\exp(V_{nj}(\boldsymbol\theta))}.
$$

**Prediction** computes these probabilities for specified parameter values. **Estimation** uses observed choices to learn which parameter values are most consistent with the data.

For example, if time enters utility as \(\beta_T T_{ni}\), the estimation problem determines \(\beta_T\) rather than choosing it arbitrarily.

## 3.2 Likelihood: how plausible are the observed choices?

Define a choice indicator

$$
y_{ni}=\begin{cases}
1,&\text{if alternative }i\text{ was chosen on occasion }n,\\
0,&\text{otherwise.}
\end{cases}
$$

Exactly one indicator equals one in each choice occasion. Conditional on the explanatory data and under **independence across choice occasions**, the joint probability of the observed choices is

$$
L(\boldsymbol\theta)=
\prod_{n=1}^N\prod_{i\in\mathcal C_n}
\big[P_{ni}(\boldsymbol\theta)\big]^{y_{ni}}.
$$

This expression is the **likelihood function**. The observed choices are treated as fixed; \(\boldsymbol\theta\) varies. The maximum likelihood estimate is

$$
\widehat{\boldsymbol\theta}_{\mathrm{ML}}=
\arg\max_{\boldsymbol\theta}L(\boldsymbol\theta).
$$

### Worked example 1: comparing two candidate models

Suppose three travelers choose Car, Public Transport, and Walk, respectively. Two candidate parameter vectors give the following probabilities **to the alternatives actually chosen**:

| Choice occasion | Candidate A | Candidate B |
|---|---:|---:|
| 1 | 0.60 | 0.80 |
| 2 | 0.70 | 0.80 |
| 3 | 0.50 | 0.70 |

Their likelihoods are

$$
L_A=0.60(0.70)(0.50)=0.210,
$$

$$
L_B=0.80(0.80)(0.70)=0.448.
$$

Candidate B gives the observed sequence a higher likelihood. This does **not** mean that every individual prediction made by B is correct or that B is necessarily a better model for unseen data.

## 3.3 Why use the log-likelihood?

Multiplying many probabilities quickly produces extremely small numbers. Because the natural logarithm is increasing, maximizing likelihood is equivalent to maximizing log-likelihood:

$$
\ell(\boldsymbol\theta)=\log L(\boldsymbol\theta)
=\sum_{n=1}^N\sum_{i\in\mathcal C_n}
 y_{ni}\log P_{ni}(\boldsymbol\theta).
$$

For the two candidate models in Example 1:

$$
\ell_A=\log(0.210)\approx-1.561,
\qquad
\ell_B=\log(0.448)\approx-0.803.
$$

The larger value, \(-0.803\), is preferred. For probabilities strictly between zero and one, their logarithms are negative; **less negative means larger**.

Substituting the MNL probability into the log-likelihood gives

$$
\boxed{
\ell(\boldsymbol\theta)=
\sum_{n=1}^N\left[
V_{n,y_n}(\boldsymbol\theta)
-\log\left(\sum_{j\in\mathcal C_n}
\exp(V_{nj}(\boldsymbol\theta))\right)
\right].}
$$

This is the central estimation equation for MNL.

### Numerical stability

The term \(\log\sum_j\exp(V_{nj})\) should be evaluated using the **log-sum-exp identity**:

$$
\log\sum_j e^{V_{nj}}=
 m_n+\log\sum_j e^{V_{nj}-m_n},
\qquad m_n=\max_j V_{nj}.
$$

Subtracting the maximum leaves the probabilities unchanged while avoiding numerical overflow and many underflow problems.

## 3.4 The score and first-order conditions

The **score** is the gradient (vector of first derivatives) of log-likelihood:

$$
\mathbf s(\boldsymbol\theta)
=\frac{\partial\ell(\boldsymbol\theta)}{\partial\boldsymbol\theta}.
$$

Assume utility is linear in the parameters:

$$
V_{ni}=\mathbf x_{ni}^{\top}\boldsymbol\theta,
$$

where \(\mathbf x_{ni}\) contains attributes such as time, cost, and alternative-specific indicators. The score contributed by choice occasion \(n\) is

$$
\boxed{
\mathbf s_n(\boldsymbol\theta)=
\mathbf x_{n,y_n}-
\sum_{j\in\mathcal C_n}P_{nj}(\boldsymbol\theta)\mathbf x_{nj}.}
$$

This has a useful interpretation:

> Score = attributes of the chosen alternative − probability-weighted expected attributes.

The sample score is \(\sum_n\mathbf s_n\). At an unconstrained interior maximum, the first-order condition is

$$
\mathbf s(\widehat{\boldsymbol\theta})=\mathbf0.
$$

An optimizer searches for parameters that satisfy this condition and maximize log-likelihood. Convergence should be checked numerically; a success message alone is not a proof of global optimality for an arbitrary model.

### Worked example 2: the score for a time coefficient

Suppose three alternatives have travel times \(10\), \(20\), and \(30\) minutes and probabilities \(0.2\), \(0.5\), and \(0.3\). The traveler chooses the 20-minute alternative.

The probability-weighted expected time is

$$
0.2(10)+0.5(20)+0.3(30)=21\text{ minutes}.
$$

Therefore, the contribution to the score for a time coefficient measured **per minute** is

$$
\frac{\partial\ell_n}{\partial\beta_T}=20-21=-1.
$$

This is the slope of this observation's log-likelihood at the current parameter vector, not the estimated coefficient itself. With time in hours, the same numerical slope would be \(-1/60\).

## 3.5 Curvature, the Hessian, and information

The **Hessian** of log-likelihood is the matrix of second derivatives:

$$
\mathcal H(\boldsymbol\theta)=
\frac{\partial^2\ell(\boldsymbol\theta)}
{\partial\boldsymbol\theta\partial\boldsymbol\theta^\top}.
$$

For linear-in-parameters MNL, the negative Hessian is

$$
\boxed{
\mathbf H(\boldsymbol\theta)=
-\mathcal H(\boldsymbol\theta)=
\sum_n\left[
\sum_j P_{nj}\mathbf x_{nj}\mathbf x_{nj}^{\top}
-\bar{\mathbf x}_n\bar{\mathbf x}_n^\top
\right],}
$$

where

$$
\bar{\mathbf x}_n=\sum_jP_{nj}\mathbf x_{nj}.
$$

Thus each observation contributes a probability-weighted covariance matrix of alternative attributes. The negative log-likelihood is convex for a standard linear-utility MNL, and its Hessian is positive semidefinite. A **positive definite** Hessian near the optimum supports local identification of the estimated parameters.

Under standard regularity conditions and correct model specification, the inverse information matrix approximates the covariance of the estimator:

$$
\widehat{\operatorname{Var}}(\widehat{\boldsymbol\theta})
\approx\mathbf H(\widehat{\boldsymbol\theta})^{-1}.
$$

For numerical optimization of \(-\ell\), \(\mathbf H\) is the Hessian of the minimized objective; no extra minus sign is needed.

## 3.6 Parameter identification and normalizations

Identification means that different permissible parameter vectors yield distinguishable model probabilities.

### Utility location

Adding the same quantity \(a_n\) to all available utilities for an individual leaves MNL probabilities unchanged:

$$
\frac{e^{V_{ni}+a_n}}{\sum_j e^{V_{nj}+a_n}}
=\frac{e^{V_{ni}}}{\sum_j e^{V_{nj}}}.
$$

Consequently, only **differences in utility** are identified. With alternative-specific constants (ASCs), one ASC must be fixed, often at zero, as the reference.

### Utility scale

In random utility theory, the common scale of utility coefficients and the unobserved error scale cannot both be identified from choices alone. Standard MNL fixes the error scale by normalizing the Type-I extreme-value errors. Estimated coefficients are consequently expressed in that normalized utility scale.

### Individual-level characteristics

A characteristic such as age that appears with **the same coefficient in every alternative** adds the same quantity to every utility and cancels. Individual and household characteristics are therefore commonly interacted with alternative indicators or alternative-specific coefficients.

### Other identification concerns

Perfect collinearity, an alternative that is never chosen, or complete separation may prevent finite or unique estimates. Choice-set availability, attribute variation, and specification should be checked before interpreting parameters.

## 3.7 Statistical inference: standard errors and confidence intervals

A point estimate reports the fitted coefficient; a standard error quantifies the estimated uncertainty around that coefficient under specified assumptions.

For parameter \(k\), the conventional standard error is

$$
SE(\widehat\theta_k)=
\sqrt{\left[\mathbf H^{-1}\right]_{kk}}.
$$

A large-sample approximate 95% confidence interval is

$$
\widehat\theta_k\pm1.96\,SE(\widehat\theta_k).
$$

### Worked example 3: a confidence interval

Suppose a time coefficient is \(-5.31\) per hour, with standard error \(0.063\). Its approximate 95% interval is

$$
-5.31\pm1.96(0.063)
=\boxed{[-5.4335,\,-5.1865]}.
$$

The interval excludes zero. This supports a negative coefficient under the estimation and variance assumptions; it does not establish that time is causally responsible for choices.

## 3.8 Wald tests and statistical versus practical significance

For testing \(H_0:\theta_k=0\), the large-sample Wald statistic is

$$
z_k=\frac{\widehat\theta_k}{SE(\widehat\theta_k)}.
$$

A two-sided approximate \(p\)-value under a standard normal reference distribution is

$$
p_k=2\big[1-\Phi(|z_k|)\big].
$$

A common 5% two-sided cutoff is \(|z_k|>1.96\). Extremely small \(p\)-values should be reported as, for example, \(p<0.001\), not interpreted as a literal zero probability.

**Statistical significance is not the same as behavioural importance.** A small effect can be highly significant in a large sample. To judge practical importance, examine changes in utilities, odds, predicted probabilities, willingness to pay, or policy outcomes.

For a more general null hypothesis \(\mathbf R\boldsymbol\theta=\mathbf r\), a Wald statistic is

$$
W=(\mathbf R\widehat{\boldsymbol\theta}-\mathbf r)^\top
\left[\mathbf R\widehat{\operatorname{Var}}(\widehat{\boldsymbol\theta})\mathbf R^\top\right]^{-1}
(\mathbf R\widehat{\boldsymbol\theta}-\mathbf r),
$$

asymptotically compared with a chi-squared distribution with degrees of freedom equal to the number of independent restrictions.

## 3.9 Repeated choices and cluster-robust inference

Travel surveys can contain several trips from the same person. Those observations may be correlated because they share unobserved tastes, circumstances, or constraints.

A conventional independent-observation variance calculation may understate uncertainty. One common alternative is a **cluster-robust sandwich estimator** at the person level.

Let \(g=1,\ldots,G\) index independent clusters (people), and let the cluster score be

$$
\mathbf S_g(\widehat{\boldsymbol\theta})=
\sum_{n\in g}\mathbf s_n(\widehat{\boldsymbol\theta}).
$$

Define the observed-information matrix \(\mathbf H\) from Section 3.5. The cluster-robust covariance is

$$
\boxed{
\widehat{\operatorname{Var}}_{\mathrm{cluster}}(\widehat{\boldsymbol\theta})
=\mathbf H^{-1}
\left(\sum_{g=1}^G\mathbf S_g\mathbf S_g^\top\right)
\mathbf H^{-1}.}
$$

A common optional finite-sample adjustment, when appropriate, is

$$
\frac{G}{G-1}\frac{N-1}{N-K},
$$

where \(N\) is the number of choice occasions and \(K\) is the number of free parameters. The validity of cluster-robust inference relies on having sufficiently many independent clusters and suitable regularity assumptions.

**Important distinction:** clustering standard errors changes uncertainty estimates but does not change the fitted choice probabilities or explicitly model within-person taste dependence. A panel mixed logit, for example, addresses the dependence structurally.

### Worked example 4: why cluster scores matter

Suppose one person reports two trips and the score for a single coefficient is \(2\) on each trip. Treating them as separate contributions to the sandwich *meat* gives

$$
2^2+2^2=8.
$$

Summing scores within the person first gives

$$
(2+2)^2=16.
$$

The cluster calculation preserves the within-person cross-product. The example illustrates the calculation, not a universal rule that clustered standard errors always double or always increase.

## 3.10 Likelihood-ratio tests and nested models

Suppose a **restricted model** has \(K_R\) parameters and an **unrestricted model** adds \(q=K_U-K_R\) identifiable parameters. Both must be estimated on **the same observations** with comparable likelihood definitions.

The likelihood-ratio statistic is

$$
LR=2\left[\ell_U(\widehat{\boldsymbol\theta}_U)
-\ell_R(\widehat{\boldsymbol\theta}_R)\right].
$$

Under suitable null regularity conditions, \(LR\) is asymptotically compared with \(\chi_q^2\).

### Worked example 5: adding household characteristics

Assume two nested models estimated on the same data have log-likelihoods

$$
\ell_R=-6000, \qquad \ell_U=-5800,
$$

and the unrestricted model adds three parameters. Then

$$
LR=2[-5800-(-6000)]=400,
$$

which is compared with a chi-squared distribution with three degrees of freedom. The improvement is statistically strong under the usual test assumptions.

**Caution:** The standard chi-squared reference may not hold for boundary parameters or certain mixture specifications, and this test does not by itself measure out-of-sample predictive performance.

## 3.11 Model fit and predictive evaluation

### Null models and McFadden's pseudo-\(R^2\)

A null or reference model must be clearly defined. Two common choices are:

- **Equal-share model:** every available alternative receives equal probability (when choice sets are equal).
- **Constants-only model:** only ASCs are fitted, reproducing sample-wide alternative shares under the standard interior optimum conditions.

These are **not the same baseline**.

Given a reference log-likelihood \(\ell_0\), McFadden's pseudo-\(R^2\) is

$$
\rho^2=1-\frac{\ell(\widehat{\boldsymbol\theta})}{\ell_0}.
$$

### Worked example 6: pseudo-\(R^2\)

Suppose \(\ell_0=-1000\) and the fitted model gives \(\ell=-800\). Then

$$
\rho^2=1-\frac{-800}{-1000}=0.20.
$$

This indicates an improvement relative to the defined reference likelihood. It does **not** mean that 20% of individual travel behaviour has been explained or that classification accuracy is 20%.

### AIC and BIC

For \(K\) fitted parameters and \(N\) observations, common selection criteria are

$$
AIC=-2\ell(\widehat{\boldsymbol\theta})+2K,
$$

$$
BIC=-2\ell(\widehat{\boldsymbol\theta})+K\ln N.
$$

Smaller values are preferred when applied to comparable models and the same data; the usual \(N\) convention and model assumptions should be stated, especially for panel data.

### Out-of-sample predictions

A model should also be evaluated on independent held-out decision-makers or another defensible validation split. Two common metrics are:

**Classification accuracy**

$$
Accuracy=\frac1{N_{test}}\sum_n
\mathbb 1\left\{y_n=\arg\max_i\widehat P_{ni}\right\}.
$$

**Multiclass log loss**

$$
LogLoss=-\frac1{N_{test}}
\sum_n\log\widehat P_{n,y_n}.
$$

Accuracy evaluates the most probable alternative; log loss evaluates the quality of probabilities assigned to actual choices. For behavioural policy analysis, it is also important to study predicted responses to time, fare, and service changes. Higher held-out accuracy does not automatically establish that counterfactual predictions are correct.

## 3.12 Interpretation of an extended utility function

Consider four alternatives and use Walking as the ASC reference. A possible car utility is

$$
V_{n,car}=\alpha_{car}+\beta_T T_{n,car}
+\beta_C C_{n,car}
+\gamma_L License_n
+\gamma_1 Own1_n+\gamma_2 Own2_n.
$$

Here \(Own1_n\) and \(Own2_n\) are indicators for two ownership categories, with the omitted category as the reference.

For a change from \(License=0\) to \(License=1\), holding other covariates fixed, the odds of Car versus any *unchanged* other alternative are multiplied by

$$
\exp(\gamma_L).
$$

For example, if \(\gamma_L=0.7\), the factor is \(e^{0.7}\approx2.01\). This is an **odds multiplier**, not a doubling of Car choice probability.

The ASC for Car now refers to the baseline socioeconomic category; one should not compare the ASC numerically across specifications without considering the changed reference interpretation and other coefficients.

## 3.13 Common distinctions worth remembering

| Concept | What it answers |
|---|---|
| Probability | Given parameters and attributes, how likely is each alternative? |
| Likelihood | Given observed choices, how plausible are particular parameters? |
| MLE | Which parameters maximize observed-data likelihood? |
| Score | In what direction does log-likelihood change locally? |
| Hessian | How sharply does the objective curve near an estimate? |
| Standard error | How uncertain is a coefficient estimate under a variance model? |
| Cluster-robust SE | How does uncertainty change when scores may be correlated within clusters? |
| Wald / LR test | Does the data support particular parameter restrictions? |
| Log loss | How well does the model assign probabilities to unseen choices? |
| Elasticity | How strongly do predicted choices respond to attribute changes? |

## 3.14 Review questions and worked answers

**Q1. Why maximize log-likelihood rather than the product of probabilities?**

**Answer:** The logarithm preserves the ordering of positive likelihood values, turns products into sums, and avoids severe numerical underflow.

**Q2. If the log-likelihood improves from \(-900\) to \(-750\), by how much has the negative log-likelihood decreased?**

**Answer:** By \(150\). The fitted log-likelihood has increased by \(150\).

**Q3. A coefficient is \(-0.25\) with an estimated standard error of \(0.10\). Calculate its Wald z-statistic and approximate 95% interval.**

**Answer:** \(z=-0.25/0.10=-2.5\). The interval is \(-0.25\pm1.96(0.10)=[-0.446,-0.054]\).

**Q4. Why can individual-level age not be inserted with one identical coefficient into all four MNL utilities?**

**Answer:** It adds the same amount to all four utilities for that individual, which cancels from their probability ratios. Alternative-specific interactions are needed to identify a relative effect.

**Q5. Does applying cluster-robust standard errors change MNL coefficients?**

**Answer:** No. Given the same parameter estimates and likelihood, it changes the estimated covariance matrix and associated inference, not the fitted coefficients.

**Q6. Can a more accurate classifier necessarily predict fare-change responses more accurately?**

**Answer:** No. Classification accuracy measures the most probable observed choice, whereas a fare-change forecast depends on the behavioural structure and coefficient response to counterfactual costs.

**Q7. What is the difference between a constants-only null model and a four-way equal-share null model?**

**Answer:** Constants-only probabilities are fitted to the observed mode frequencies; an equal-share model assigns \(25\%\) to each of four available alternatives regardless of their observed frequency.

---

## Chapter summary

Maximum likelihood estimation connects individual observed choices to parameters in a probabilistic utility model. MNL estimation maximizes the sum of logged probabilities assigned to the alternatives actually chosen. The score describes the local direction of improvement, the Hessian describes curvature, and its inverse supports conventional large-sample uncertainty estimates under suitable assumptions. Repeated observations from the same person motivate cluster-robust inference. Likelihood-based tests assess restrictions and fit, while held-out log loss and behavioural sensitivity analysis evaluate additional dimensions of model usefulness.

The key principle is that **a precisely estimated or highly predictive model is not automatically a behaviourally valid counterfactual model**. Reliable travel choice analysis requires careful utility specification, valid inference, and defensible interpretation.
