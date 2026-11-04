# Chapter 4 — Nested Logit Models

**Discrete Choice Modelling: Theory, Mathematics, and Behavioural Interpretation**

> **Purpose.** This chapter develops the nested logit (NL) model from the random utility perspective, derives its choice probabilities, explains substitution patterns and the dissimilarity parameter, and introduces estimation and model assessment. It is a self-contained theoretical chapter; no programming background is required.

## Learning objectives

After reading this chapter, you should be able to:

1. Explain why the independence of irrelevant alternatives (IIA) property of multinomial logit can be restrictive.
2. Specify a behavioural hypothesis about which alternatives should belong to the same nest.
3. Derive and calculate the conditional, nest-level, and unconditional NL choice probabilities.
4. Interpret the inclusive value and the dissimilarity parameter under an explicitly stated parameterization.
5. Distinguish within-nest from across-nest substitution.
6. Explain how NL is estimated, compared with MNL, and evaluated on new observations.
7. Recognize what a statistically supported nesting structure does—and does not—tell us about actual behaviour.

---

## 4.1 The problem that nested logit addresses

Discrete choice models describe decisions between mutually exclusive alternatives. Consider a commuter choosing among walking, cycling, public transport (PT), and a car. The observed part of each alternative's utility can include travel time, money cost, reliability, and characteristics of the decision-maker. Some influences will remain unobserved.

A conventional **multinomial logit (MNL)** model assumes that the unobserved utility disturbances are independent and identically distributed type-I extreme value (Gumbel) random variables, with a common scale. Its choice probability is

$$
P_{ni}^{\mathrm{MNL}}=\frac{\exp(V_{ni})}{\sum_{j\in\mathcal C_n}\exp(V_{nj})},
$$

where $n$ indexes choice observations, $i$ indexes alternatives, $\mathcal C_n$ is the choice set, and $V_{ni}$ is systematic utility. The normalized common error scale is absorbed into the coefficients.

For two alternatives $i$ and $j$, MNL implies

$$
\frac{P_{ni}^{\mathrm{MNL}}}{P_{nj}^{\mathrm{MNL}}}
=\exp(V_{ni}-V_{nj}).
$$

This odds ratio is independent of the attributes of *other* alternatives. This property is known as **independence of irrelevant alternatives (IIA)**.

### Why can IIA matter?

Imagine one bus route and one substantially different rail route. A second bus route appears that is nearly a copy of the first. If all three have equal systematic utilities, elementary MNL assigns one third of demand to each, even if the new bus option adds little independent attractiveness. This is an illustration of the **red-bus/blue-bus problem**. It is not a claim that a genuinely new bus service can never attract additional passengers: more vehicles may reduce waiting, improve reliability, or expand accessibility. The problem arises when the *choice representation* treats highly dependent options as independent draws.

Nested logit allows some alternatives to share a structured component of unobserved utility, relaxing the global IIA restriction. It remains a probabilistic model: it does not require passengers to explicitly think in two consecutive stages.

---

## 4.2 Behavioural nesting hypotheses

A **nest** is a group of alternatives hypothesized to have a greater degree of unobserved similarity with one another than with alternatives outside the group.

For example, two candidate structures for a four-mode problem are:

**Structure A — Active versus motorized**

$$
\mathcal N_A=\{\text{Walk, Cycle}\},\qquad
\mathcal N_M=\{\text{PT, Car}\}.
$$

**Structure B — Non-car versus car**

$$
\mathcal N_{NC}=\{\text{Walk, Cycle, PT}\},\qquad
\mathcal N_C=\{\text{Car}\}.
$$

Neither classification is automatically correct. Physical similarity, vehicle technology, and the label “sustainable transport” do not establish the required *unobserved utility correlation*. A proposed nest should follow a substantive behavioural argument, then be confronted with data.

A nest with a single alternative is called a **singleton nest**. Its within-nest choice probability is one. Its dissimilarity coefficient cannot be identified separately and is conventionally fixed at one.

### Two kinds of similarity

- **Observed similarity:** alternatives have similar travel times, costs, routes, or service attributes. These can enter systematic utility directly.
- **Unobserved similarity:** alternatives share omitted influences, producing dependence in their utility disturbances. Nested logit represents a particular pattern of this dependence.

Observed similarity does not, by itself, prove unobserved correlation. Likewise, strong estimated dependence does not identify a particular psychological trait unless the model and data justify that interpretation.

---

## 4.3 Random utility foundation and notation

The random utility formulation is

$$
U_{ni}=V_{ni}+\varepsilon_{ni}.
$$

The individual chooses alternative $i$ if its total utility exceeds those of all other available alternatives:

$$
P_{ni}=\Pr\big(U_{ni}\ge U_{nj}\quad \forall j\in\mathcal C_n\big).
$$

Nested logit can be derived from a **generalized extreme value (GEV)** joint distribution for the errors. Within a nest, the errors are allowed to be dependent; across different nests, the structure is more restricted. The resulting choice probabilities have a convenient decomposition.

We use the following **two-level dissimilarity-parameter formulation** throughout:

- Each alternative $i$ belongs to one and only one nest $m$.
- The upper-level (between-nest) scale is normalized to one.
- $\lambda_m$ denotes the dissimilarity parameter for nest $m$.
- For the standard random-utility-consistent formulation, $0<\lambda_m\le 1$.
- $V_{ni}$ is systematic utility on the upper-level normalized scale.

This convention matters: other textbooks and software packages can use the reciprocal of $\lambda_m$, or differently scaled “inclusive value coefficients.” Never interpret an estimated number without checking the parameterization.

---

## 4.4 Derivation of the nested logit probabilities

Let $\mathcal C_{nm}$ denote the alternatives available to individual $n$ in nest $m$. Define the **within-nest logsum**:

$$
I_{nm}=\ln\left(\sum_{j\in\mathcal C_{nm}}
\exp\left(\frac{V_{nj}}{\lambda_m}\right)\right).
\tag{4.1}
$$

It is sometimes called the *inclusive value* (IV), although authors may use “inclusive value” for $\lambda_m I_{nm}$ instead. We distinguish the unscaled logsum $I_{nm}$ from the *composite nest utility* below.

### Step 1 — Conditional choice within a nest

Once nest $m$ is selected, the conditional probability of choosing $i$ is

$$
\boxed{P_n(i\mid m)=
\frac{\exp(V_{ni}/\lambda_m)}{
\sum_{j\in\mathcal C_{nm}}\exp(V_{nj}/\lambda_m)}.}
\tag{4.2}
$$

This is a logit formula *within* the nest, with utilities divided by $\lambda_m$. The conditional probabilities add to one within every nonempty nest.

### Step 2 — Nest-level utility

Define

$$
W_{nm}=\lambda_m I_{nm}
=\lambda_m\ln\left(\sum_{j\in\mathcal C_{nm}}
\exp(V_{nj}/\lambda_m)\right).
\tag{4.3}
$$

$W_{nm}$ summarizes the attractiveness of all available alternatives in nest $m$. If the model has separately specified upper-level explanatory variables $Z_{nm}$ with coefficients $\gamma$, one may instead write $W_{nm}=\gamma^\top Z_{nm}+\lambda_m I_{nm}$, provided the overall specification is properly identified. For simplicity, we use (4.3).

The probability of choosing the nest is

$$
\boxed{P_n(m)=\frac{\exp(W_{nm})}{\sum_{k\in\mathcal M_n}\exp(W_{nk})}.}
\tag{4.4}
$$

### Step 3 — Combine the probabilities

The unconditional probability of alternative $i\in m$ is

$$
\boxed{P_{ni}^{\mathrm{NL}}=P_n(m)P_n(i\mid m).}
\tag{4.5}
$$

In an equivalent single expression,

$$
P_{ni}^{\mathrm{NL}}=
\underbrace{\frac{\exp(V_{ni}/\lambda_m)}{
\sum_{j\in\mathcal C_{nm}}\exp(V_{nj}/\lambda_m)}}_{\text{within-nest choice}}
\underbrace{\frac{\left[\sum_{j\in\mathcal C_{nm}}\exp(V_{nj}/\lambda_m)\right]^{\lambda_m}}{
\sum_k\left[\sum_{j\in\mathcal C_{nk}}\exp(V_{nj}/\lambda_k)\right]^{\lambda_k}}}_{\text{nest choice}}.
\tag{4.6}
$$

Because $\sum_{i\in m}P_n(i\mid m)=1$ and $\sum_m P_n(m)=1$, summing (4.5) over all alternatives gives exactly one.

---

## 4.5 Fully worked numerical example

Suppose a traveller has these hypothetical systematic utilities:

| Alternative | Nest | Systematic utility |
|:--|:--|--:|
| Walk | Non-car | $-2.0$ |
| Cycle | Non-car | $-3.0$ |
| Public transport | Non-car | $-1.0$ |
| Car | Car (singleton) | $-2.5$ |

Let $\lambda_{NC}=0.60$ and $\lambda_C=1$.

### A. Conditional probabilities inside Non-car

First divide utilities by $0.60$:

$$
\left(\frac{-2}{0.6},\frac{-3}{0.6},\frac{-1}{0.6}\right)
=(-3.3333,-5.0000,-1.6667).
$$

The common denominator for the three exponential terms is approximately

$$
e^{-3.3333}+e^{-5}+e^{-1.6667}\approx0.23127.
$$

Thus, to four decimals,

$$
\begin{aligned}
P(\text{Walk}\mid NC)&\approx0.1542,\\
P(\text{Cycle}\mid NC)&\approx0.0291,\\
P(\text{PT}\mid NC)&\approx0.8166.
\end{aligned}
$$

### B. Inclusive value and probability of the nest

The unscaled logsum is

$$
I_{NC}=\ln(0.23127)\approx-1.4641.
$$

The composite utility is

$$
W_{NC}=0.60(-1.4641)\approx-0.8785.
$$

For the singleton car nest, $W_C=V_C=-2.5$. Therefore,

$$
P(NC)=\frac{e^{-0.8785}}{e^{-0.8785}+e^{-2.5}}
\approx0.8350,
\qquad P(C)\approx0.1650.
$$

### C. Final mode probabilities

Multiply within-nest and nest probabilities:

| Alternative | $P(i\mid NC)$ where applicable | $P(i)$ |
|:--|--:|--:|
| Walk | 0.1542 | **0.1288** |
| Cycle | 0.0291 | **0.0243** |
| Public transport | 0.8166 | **0.6819** |
| Car | — | **0.1650** |
| **Total** | | **1.0000** |

For comparison, MNL ($\lambda_{NC}=1$) produces approximately $(0.2131,0.0784,0.5793,0.1293)$ in the same mode order. The two models use the *same systematic utilities* in this example; only their dependence structures differ.

The result is not a universal claim that NL always raises PT demand. Direction and magnitude depend on the utilities, nest definitions, and parameter values.

---

## 4.6 Interpreting the dissimilarity parameter

The parameter $\lambda_m$ affects the relative weight of alternatives **within** nest $m$ compared with the upper-level decision.

- **$\lambda_m=1$:** the nest has the standard MNL scaling. If all nests have $\lambda_m=1$, the entire NL model is algebraically identical to MNL.
- **$0<\lambda_m<1$:** relative utility differences inside the nest are amplified through the division by $\lambda_m$, leading to a different substitution structure and a model consistent with positive within-nest error dependence under the standard GEV interpretation.
- **$\lambda_m\to0^+$:** within-nest choices become increasingly concentrated on the alternative(s) with the highest systematic utility. The precise limiting model requires care, especially for ties and identification.

Under the familiar homoscedastic two-level GEV normalization, an often-quoted expression for the within-nest error correlation is

$$
\operatorname{Corr}(\varepsilon_{ni},\varepsilon_{nj})=1-\lambda_m^2
\qquad(i\ne j,\;i,j\in m).
\tag{4.7}
$$

This is a **model-implied latent-error correlation under that normalization**, not a measured correlation between observed choices, travel times, or individual responses. It is not automatically portable to other NL parameterizations.

For illustration, $\lambda=0.60$ gives $1-0.60^2=0.64$. A value $\lambda=0.4742$ gives approximately $0.775$. Those numerical translations are meaningful only under the assumptions stated above.

### Does NL imply a real two-stage mental process?

Not necessarily. The decomposition $P(i)=P(i\mid m)P(m)$ is an exact probability identity under the model. It does **not** show that people consciously choose “non-car” first and “walking versus cycling versus PT” second. The nesting structure is best interpreted as a testable statistical and behavioural representation.

---

## 4.7 The IIA property: relaxed, not abolished

MNL has the global IIA property: for any $i,j$, the odds ratio depends only on $V_i-V_j$.

NL modifies this in a specific way.

### Alternatives in the same nest

For $i,j\in m$, the common nest probability cancels:

$$
\frac{P_{ni}}{P_{nj}}
=\frac{P_n(i\mid m)}{P_n(j\mid m)}
=\exp\left(\frac{V_{ni}-V_{nj}}{\lambda_m}\right).
\tag{4.8}
$$

Thus **IIA still holds *within* a nest**, conditional on the same nest membership and the specified utility. Changes to a third option inside the nest do not change this pairwise odds ratio.

### Alternatives in different nests

For $i\in m$ and $j\in k$, the ratio generally includes the two different inclusive values and therefore depends on other alternatives in their nests. Global IIA no longer holds.

**Key principle:** NL replaces unrestricted MNL substitution with a *structured* pattern, not with completely arbitrary substitution. More flexible structures, such as cross-nested logit or mixed logit, may be needed when an alternative belongs naturally to multiple overlapping groups.

---

## 4.8 How substitution changes: a small analytic derivation

Let $i$ and $j$ be two distinct alternatives in the same nest $m$. Suppress the individual index, and write $q_i=P(i\mid m)$ and $s_m=P(m)$, so $P_i=s_mq_i$.

For a change in the systematic utility of alternative $i$, the conditional derivative is

$$
\frac{\partial q_i}{\partial V_i}=\frac{1}{\lambda_m}q_i(1-q_i),
\qquad
\frac{\partial q_j}{\partial V_i}=-\frac{1}{\lambda_m}q_iq_j.
$$

From (4.3), $\partial W_m/\partial V_i=q_i$. The upper-level derivative is

$$
\frac{\partial s_m}{\partial V_i}=s_m(1-s_m)q_i.
$$

By the product rule, for $j\ne i$ **within the same nest**,

$$
\boxed{\frac{\partial P_j}{\partial V_i}
=-P_j q_i\left(\frac{1}{\lambda_m}-1+s_m\right).}
\tag{4.9}
$$

For an alternative $k$ in a different nest,

$$
\boxed{\frac{\partial P_k}{\partial V_i}=-P_kP_i.}
\tag{4.10}
$$

For the chosen alternative $i$ itself,

$$
\boxed{\frac{\partial P_i}{\partial V_i}
=P_i\left[\frac{1-q_i}{\lambda_m}+(1-s_m)q_i\right].}
\tag{4.11}
$$

These derivatives show the structure explicitly: as $\lambda_m$ falls below one, the within-nest cross-response has an extra term that is absent in MNL. However, **net aggregate changes also depend on the initial probabilities, utility coefficients, and trip characteristics**; a universal numerical substitution share cannot be inferred from $\lambda_m$ alone.

### A policy interpretation

If PT becomes less attractive, a non-car nest may predict relatively stronger increases in walking and cycling than a comparable MNL. That is a model-based pattern, not an observed person-by-person transfer. Predicted changes in aggregate mode shares describe *net changes in expected probabilities*.

---

## 4.9 Observable attributes, socioeconomic characteristics, and nests

The systematic utility may include alternative attributes and person-specific variables:

$$
V_{ni}=\alpha_i+\beta_T T_{ni}+\beta_C C_{ni}
+\boldsymbol\gamma_i^{\top}\mathbf z_n,
\tag{4.12}
$$

where $\mathbf z_n$ might include age, driving licence status, household car ownership, or other observed characteristics. These variables can influence particular alternatives through alternative-specific coefficients.

For example, a licence indicator can enter car utility:

$$
V_{n,Car}=\alpha_{Car}+\beta_T T_{n,Car}
+\beta_C C_{n,Car}+\gamma_L License_n.
$$

If the identical term $\gamma z_n$ is added to **all** alternatives, it cancels from every logit probability and is not identified. The same identification issue applies in NL, which depends ultimately on utility differences and scale normalizations.

### Why hold explanatory variables fixed across model families?

To compare MNL and NL meaningfully, it is useful to estimate both with the same observed attributes, choice-set rules, and train/test partition. Their **coefficient estimates are allowed to differ**, because changing the stochastic structure changes the likelihood function. This design isolates the contribution of model structure more clearly than changing both covariates and error dependence at once.

Do not compare time coefficients across differently scaled models as if they had an invariant cardinal interpretation.

---

## 4.10 Maximum likelihood estimation

For $N$ observed choices, let $y_n$ be the alternative selected by observation $n$, and let $\boldsymbol\theta$ contain utility coefficients and dissimilarity parameters. The log-likelihood is

$$
\ell(\boldsymbol\theta)=\sum_{n=1}^{N}
\ln P_{n,y_n}^{\mathrm{NL}}(\boldsymbol\theta).
\tag{4.13}
$$

The maximum likelihood estimator is

$$
\hat{\boldsymbol\theta}=\arg\max_{\boldsymbol\theta\in\Theta}
\ell(\boldsymbol\theta),
\tag{4.14}
$$

with $0<\lambda_m\le1$ for each non-singleton nest under our convention.

NL's probability expression is nonlinear in both utility coefficients and $\lambda_m$, so numerical optimization is normally needed. Good practice includes reasonable starting values, numerical-stability calculations (log-sum-exp), checking the gradient and objective near the solution, and comparing multiple initializations when appropriate.

### MNL as an embedded special case

If $\lambda_m=1$ for every nest, then

$$
W_{nm}=\ln\sum_{j\in\mathcal C_{nm}}e^{V_{nj}},
$$

and the NL formula simplifies to

$$
P_{ni}=\frac{e^{V_{ni}}}{\sum_{j\in\mathcal C_n}e^{V_{nj}}}.
$$

Thus NL **nests MNL mathematically**. If both are estimated on the same observations and with the same systematic utility specification, the *global maximum* NL log-likelihood cannot be below the MNL maximum. A tiny reversal in a numerical run is a sign to check convergence tolerances, initial values, and optimization accuracy; it is not a theoretical NL disadvantage.

### How much evidence does $\hat\lambda<1$ provide?

An interior estimate below one suggests that the specified dependence structure may fit the data better. But an estimate alone does not establish statistical significance; uncertainty also matters. For example, a conventional approximate confidence interval is

$$
\hat\lambda\pm1.96\,SE(\hat\lambda).
$$

A large improvement in log-likelihood also deserves attention, but the usual likelihood-ratio reference distribution must be considered carefully: **$\lambda=1$ is a boundary of the allowed parameter space**. Standard interior-point $\chi^2$ approximations may not apply directly; a suitable boundary-aware test or bootstrap can be used. With several boundary parameters, the null distribution can be more complicated.

---

## 4.11 Conventional and cluster-robust inference

Let $H$ be the observed Hessian of the **negative** log-likelihood at the optimum:

$$
H=\left.\frac{\partial^2[-\ell(\boldsymbol\theta)]}
{\partial\boldsymbol\theta\,\partial\boldsymbol\theta^\top}
\right|_{\boldsymbol\theta=\hat{\boldsymbol\theta}}.
$$

Under standard regularity assumptions, the conventional covariance estimate is approximately

$$
\widehat{\operatorname{Var}}_{conv}(\hat{\boldsymbol\theta})=H^{-1}.
\tag{4.15}
$$

If one individual supplies several trip observations, the associated scores may be dependent. Define the score for observation $n$ as $s_n=\partial\ln P_{n,y_n}/\partial\boldsymbol\theta$, and aggregate them for individual $g$:

$$
S_g=\sum_{n\in g}s_n.
$$

Then the basic cluster-robust sandwich covariance is

$$
\boxed{\widehat{\operatorname{Var}}_{cl}(
\hat{\boldsymbol\theta})=H^{-1}
\left(\sum_{g=1}^{G}S_gS_g^{\top}\right)H^{-1}.}
\tag{4.16}
$$

Finite-sample adjustments may be applied. Its square-root diagonal elements give cluster-robust standard errors.

**Important distinction:** cluster-robust inference changes the estimated uncertainty, not the fitted coefficients or the underlying model's substitution structure. Panel mixed logit and other longitudinal choice models can explicitly model repeated-choice dependence; clustered standard errors do not replace such modelling.

### Example: interpreting an estimated lambda

Suppose $\hat\lambda=0.4742$, with conventional $SE=0.0078$ and cluster-robust $SE=0.0101$. Then the respective approximate 95% intervals are

$$
[0.4589,0.4895]\quad\text{and}\quad[0.4544,0.4940].
$$

The latter is wider because the estimated uncertainty is larger after accounting for within-person clustering. These are *illustrative inputs for interpretation*: they do not, by themselves, validate the nest's behavioural meaning or the accuracy of policy responses.

---

## 4.12 Comparing NL with MNL responsibly

Different evaluation criteria answer different questions:

| Criterion | Main question | Better direction |
|:--|:--|:--|
| Training log-likelihood | How well does the fitted model describe estimation observations? | Higher |
| AIC / BIC | Is extra fit worth added parameter complexity under the criterion's assumptions? | Lower |
| Held-out log loss | How much probability is assigned to choices not seen during estimation? | Lower |
| Accuracy | How often is the highest-probability option the actual choice? | Higher |
| Aggregate share error | How close are predicted shares to observed shares? | Smaller absolute error |
| Calibration | Do predicted probabilities correspond to observed frequencies? | Better agreement |
| Policy sensitivity | What behavioural changes does the model imply when attributes change? | Plausible and independently supported |

For $K$ estimated parameters and training sample size $N$,

$$
AIC=2K-2\ell(\hat\theta),\qquad
BIC=K\ln N-2\ell(\hat\theta).
$$

When repeated choices are present, the usual BIC effective-sample-size interpretation deserves care. For predictive validation, **split by decision-maker**, not randomly by trip, if the goal is prediction for unseen individuals. Keep the same holdout set across candidate models.

An NL model may improve held-out log loss but not accuracy: it can give higher probability to the observed choices without changing—or even while slightly worsening—the identity of the most-probable option. For travel demand, probabilities, aggregate shares, and substitution patterns often matter more than a single classification score.

### Model fit does not equal policy validity

Suppose NL predicts more substitution from PT to walking after a fare increase than MNL does. If there are no observed before-and-after policy data, **neither forecast is verified by model fit alone**. Scenario analysis is useful, but it is a prediction conditional on assumptions about the model and the intervention.

---

## 4.13 Elasticities and policy sensitivity

Suppose utility of alternative $i$ depends on attribute $x_{ni}$ linearly:

$$
V_{ni}=\cdots+\beta_x x_{ni}.
$$

For a small change, the individual-level probability derivative is

$$
\frac{\partial P_{ni}}{\partial x_{ni}}
=\beta_x\frac{\partial P_{ni}}{\partial V_{ni}}.
$$

A point elasticity is

$$
E_{ni,x}=\frac{\partial P_{ni}}{\partial x_{ni}}
\frac{x_{ni}}{P_{ni}}.
\tag{4.17}
$$

NL derivatives differ depending on whether a changed alternative and an affected alternative share a nest (see Section 4.8). Therefore, a single MNL-style cross-elasticity formula should **not** automatically be applied to all NL alternatives.

In practical scenario work, one can calculate probabilities with original attributes ($P_{ni}^0$) and counterfactual attributes ($P_{ni}^1$), then average:

$$
S_i^0=\frac{1}{N}\sum_n P_{ni}^0,
\qquad S_i^1=\frac{1}{N}\sum_n P_{ni}^1.
$$

The change in share in **percentage points** is $100(S_i^1-S_i^0)$. If an attribute is increased by a common 10%, the baseline-referenced *finite-change* elasticity is

$$
E_{i,x}^{(10\%)}\approx
\frac{(S_i^1-S_i^0)/S_i^0}{0.10}.
\tag{4.18}
$$

For a 10% decrease the denominator is $-0.10$. This is not a point derivative and not necessarily an arc-elasticity using midpoint denominators. If some initial prices are zero, a proportional increase leaves those prices unchanged; the scenario must be described accordingly.

---

## 4.14 Identification, specification, and practical cautions

1. **Normalization:** one alternative-specific constant typically serves as a reference. Utility scale is normalized, so coefficient magnitudes are not automatically comparable across differently scaled specifications.
2. **Availability:** probabilities must be formed over alternatives actually available in the specified choice set. Household ownership or possession of a driving licence may help explain choices but need not determine availability perfectly.
3. **Nest assignment:** standard two-level NL assigns each alternative to only one nest. Overlapping memberships need another model class, such as cross-nested logit.
4. **Data support:** rare alternatives and weak variation in observed attributes can make some nest parameters difficult to estimate reliably.
5. **Boundary estimates:** a fitted $\lambda=1$ can indicate that the added nest dependence provides no fit improvement under the current specification; it does not prove that all unobserved utility errors are independent in reality.
6. **Causal language:** an estimated licence or ownership coefficient is a conditional association unless the design supports a causal interpretation. Likewise, changes predicted under a fare scenario are model-based counterfactuals.
7. **Likelihood optimization:** convergence flags and reported parameter values should be checked against the objective at the nested MNL solution and, where warranted, alternative initializations.
8. **Welfare and value of time:** the standard top-level expected-maximum utility is related to $\ln\sum_m e^{W_m}$ under the stated scale, but utility-to-money conversion must respect the full structure and the monetary coefficient used. A coefficient ratio can be informative for a common linear time and cost term; interpretation requires care when costs, time components, or scales differ across levels and alternatives.

---

## 4.15 Conceptual and numerical exercises

Try these problems before reading the answers.

**Exercise 1 — Nest membership.** Give one plausible reason to group walking and cycling and one reason that grouping might fail after observed time, distance, and other characteristics have been accounted for.

**Exercise 2 — Probability decomposition.** If $P(\mathrm{NonCar})=0.70$ and $P(\mathrm{Cycle}\mid\mathrm{NonCar})=0.15$, calculate $P(\mathrm{Cycle})$.

**Exercise 3 — MNL equivalence.** Why must NL collapse to MNL when all $\lambda_m=1$? Show the cancellation algebraically.

**Exercise 4 — Correlation interpretation.** In the parameterization of this chapter, an estimated $\lambda=0.8$ is reported. Calculate $1-\lambda^2$ and explain why this is not a correlation between people's observed choices.

**Exercise 5 — Substitution.** For alternatives $i$ and $j$ in the same nest, derive $P_i/P_j$ and state whether this odds ratio depends on the utility of a third alternative in that nest.

**Exercise 6 — Likelihood comparison.** A matched-sample MNL fit has $\ell=-510$ and a one-parameter-higher NL model has $\ell=-490$. Calculate the LR statistic and the AIC for the two models if MNL has eight parameters. State one caveat about formal inference.

**Exercise 7 — Uncertainty.** An NL model estimates $\lambda=0.55$ with conventional $SE=0.04$ and clustered $SE=0.06$. Form approximate 95% intervals and explain the change.

**Exercise 8 — Model evaluation.** NL has slightly lower accuracy than MNL but lower held-out log loss. Is this mathematically possible? Which model is necessarily better for fare-policy forecasting?

### Answers

**1.** Walking and cycling are both active, usually nonmotorized modes and might share unobserved preferences; the assumption may fail if riders distinguish them sharply on weather exposure, physical exertion, safety, or bicycle access, especially when these are already represented by observed covariates.

**2.** $0.70\times0.15=0.105$, or **10.5%**.

**3.** With $\lambda_m=1$, $e^{W_m}=\sum_{j\in m}e^{V_j}$. Multiplying conditional probability $e^{V_i}/\sum_{j\in m}e^{V_j}$ by nest probability $(\sum_{j\in m}e^{V_j})/(\sum_k\sum_{j\in k}e^{V_j})$ cancels the within-nest sum, giving MNL.

**4.** $1-0.8^2=0.36$; the value concerns the latent errors in a specific GEV normalization, not empirical correlations in realized binary or categorical choices.

**5.** $P_i/P_j=\exp[(V_i-V_j)/\lambda_m]$; it does not depend on the third alternative. NL retains IIA **within** a nest.

**6.** $LR=2(-490+510)=40$. $AIC_{MNL}=2(8)-2(-510)=1036$; $AIC_{NL}=2(9)-2(-490)=998$. The null $\lambda=1$ is a boundary point, so the usual interior $\chi^2_1$ reference must not be assumed without justification.

**7.** Conventional: $0.55\pm1.96(0.04)=[0.4716,0.6284]$. Clustered: $0.55\pm1.96(0.06)=[0.4324,0.6676]$. The second is wider because it reflects greater estimated uncertainty when within-cluster dependence is permitted.

**8.** Yes. Accuracy tests the identity of the highest-probability option; log loss evaluates assigned probability to the observed option. Neither result alone demonstrates which model's *policy counterfactuals* are more accurate.

---

## 4.16 Chapter summary

Nested logit extends MNL by providing a structured form of dependence among the random components of utility. Its core probability identity is

$$
\boxed{P(i)=P(i\mid m)P(m).}
$$

Under the stated dissimilarity normalization, $0<\lambda_m\le1$; MNL is recovered when all $\lambda_m=1$. The inclusive value summarizes within-nest attractiveness, while the nest probability allocates probability across groups. NL retains IIA within a nest but relaxes it across nests. Good empirical practice combines a defensible nesting hypothesis, careful maximum likelihood estimation, appropriate clustered inference for repeated observations, held-out predictive evaluation, and cautious interpretation of policy sensitivities.

### Selected references for further reading

- McFadden, D. (1978). *Modelling the choice of residential location.* In A. Karlqvist et al. (Eds.), **Spatial Interaction Theory and Planning Models**. North-Holland.
- Williams, H. C. W. L. (1977). *On the formation of travel demand models and economic evaluation measures of user benefit.* **Environment and Planning A**, 9(3), 285–344.
- Ben-Akiva, M., & Lerman, S. R. (1985). **Discrete Choice Analysis: Theory and Application to Travel Demand**. MIT Press.
- Train, K. E. (2009). **Discrete Choice Methods with Simulation**, 2nd ed. Cambridge University Press. See especially the treatment of GEV and nested logit models.

*End of Chapter 4.*
