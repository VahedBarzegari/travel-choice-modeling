# Chapter 2 — The Multinomial Logit Model

*Discrete Choice Modelling: Theory, Mathematics, and Behavioural Interpretation*

## Learning objectives

By the end of this chapter, you should be able to:

1. Explain how the multinomial logit (MNL) model follows from random utility maximization.
2. Calculate choice probabilities and interpret changes in utility.
3. Explain why utility differences—not absolute utility levels—determine choices.
4. Distinguish alternative-specific constants from observed attributes.
5. Understand the independence of irrelevant alternatives (IIA) property and its practical consequences.
6. Calculate and interpret marginal effects and simple choice elasticities.
7. Recognize what MNL can and cannot establish about traveller behaviour.

---

## 2.1 From utility to probability

Consider traveller $n$ choosing one alternative from a set $\mathcal C_n$, such as walking, cycling, public transport, or driving. Under **random utility maximization (RUM)**, the traveller chooses the alternative with the highest utility:

$$
y_n=\arg\max_{i\in\mathcal C_n} U_{ni}.
$$

The analyst does not observe all determinants of utility. We therefore write

$$
U_{ni}=V_{ni}+\varepsilon_{ni},
$$

where $V_{ni}$ is **systematic utility** (constructed from observed variables and model parameters), and $\varepsilon_{ni}$ represents unobserved influences. “Random” describes the analyst's incomplete information; it need not imply that the traveller is making a random or irrational decision.

For a mode-choice model, a possible specification is

$$
V_{ni}=\alpha_i+\beta_T T_{ni}+\beta_C C_{ni},
$$

where $T_{ni}$ is travel time, $C_{ni}$ is monetary cost, $\alpha_i$ is an alternative-specific constant (ASC), and the coefficients measure changes in systematic utility per unit of the relevant variable.

Choice probabilities arise because the analyst does not know the realized errors:

$$
P_{ni}=\Pr\!\left(U_{ni}\ge U_{nj}\ \text{for all }j\in\mathcal C_n\right).
$$

## 2.2 The distributional assumption behind MNL

The standard MNL model assumes that the $\varepsilon_{ni}$ are **independent and identically distributed (i.i.d.) Type-I extreme-value (Gumbel)** errors across alternatives, conditional on the observed explanatory variables. The Gumbel cumulative distribution function with location zero and scale $\mu>0$ is

$$
F(\varepsilon)=\exp\left[-\exp\left(-\frac{\varepsilon}{\mu}\right)\right].
$$

Its variance is $\pi^2\mu^2/6$. Under these assumptions, integrating over the unobserved errors gives

$$
\boxed{P_{ni}=\frac{\exp(V_{ni}/\mu)}{\displaystyle\sum_{j\in\mathcal C_n}\exp(V_{nj}/\mu)}}.
$$

Only the ratios of utility coefficients to the scale $\mu$ are identified. We therefore **normalize $\mu=1$**, giving the familiar formula

$$
\boxed{P_{ni}=\frac{e^{V_{ni}}}{\sum_{j\in\mathcal C_n}e^{V_{nj}}}}.
$$

The common scale normalization matters when comparing coefficients across samples or separately estimated models: coefficients can differ partly because the scale of unobserved utility differs.

### A short derivation

The Type-I extreme-value density is

$$
f(\varepsilon)=\frac{1}{\mu}\exp(-\varepsilon/\mu)\exp[-\exp(-\varepsilon/\mu)].
$$

Condition on the error in alternative $i$. For $i$ to be chosen, all other errors must satisfy $\varepsilon_{nj}\le V_{ni}+\varepsilon_{ni}-V_{nj}$. By independence,

$$
P_{ni}=\int_{-\infty}^{\infty} f(\varepsilon)\prod_{j\ne i}F(V_{ni}+\varepsilon-V_{nj})\,d\varepsilon.
$$

Let $A=\sum_j\exp[(V_{nj}-V_{ni})/\mu]$ and $z=\exp(-\varepsilon/\mu)$. Substitution reduces the integral to

$$
P_{ni}=\int_0^\infty e^{-Az}\,dz=\frac{1}{A}
=\frac{e^{V_{ni}/\mu}}{\sum_j e^{V_{nj}/\mu}}.
$$

This is why the logit formula is a consequence of the specified **error distribution**, rather than an arbitrary rule for converting utility into probabilities.

## 2.3 A worked numerical example

Suppose three travel alternatives have the following **hypothetical systematic utilities**:

| Alternative | Utility $V_i$ |
|---|---:|
| Walk | $-2$ |
| Public transport | $-1$ |
| Car | $-2.5$ |

Their exponentials are approximately $0.1353$, $0.3679$, and $0.0821$. The sum is $0.5853$, so

$$
P_{\mathrm{walk}}=\frac{e^{-2}}{e^{-2}+e^{-1}+e^{-2.5}}\approx0.2312,
$$

$$
P_{\mathrm{PT}}\approx0.6286,\qquad P_{\mathrm{car}}\approx0.1402.
$$

Thus public transport has the highest predicted choice probability, but its selection is not certain. These are **probabilities conditional on the assumed attributes, coefficients, and choice set**, not an observed distribution of trips.

### Utility is not probability

A utility change of $+0.5$ does **not** mean a 50-percentage-point increase in probability. The effect depends on the other utilities. If only the PT utility rises from $-1$ to $-0.5$, then

$$
P'_{\mathrm{PT}}=\frac{e^{-0.5}}{e^{-2}+e^{-0.5}+e^{-2.5}}\approx0.7361,
$$

an increase of roughly **10.8 percentage points**, not 50 points.

## 2.4 Utility normalization and the role of ASCs

Adding the same constant $k$ to every alternative's utility leaves all probabilities unchanged:

$$
\frac{e^{V_{ni}+k}}{\sum_j e^{V_{nj}+k}}
=\frac{e^{V_{ni}}}{\sum_j e^{V_{nj}}}.
$$

Therefore, only **utility differences** can be identified. If an ASC is included for each alternative, one must be fixed as a reference, typically $\alpha_{\mathrm{walk}}=0$.

An ASC captures average differences in systematic utility between an alternative and the reference **after accounting for variables included in the specification**. It may absorb omitted attributes or persistent differences in the alternatives. A negative car ASC does not, by itself, imply that people prefer walking to driving: time, cost, availability, and other attributes also contribute to utility.

Personal characteristics (such as age or car ownership) are constant across alternatives within a single choice situation. Adding the **same** term $\gamma Z_n$ to every $V_{ni}$ cancels from the probability formula. To estimate its effect in MNL, it must vary by alternative—for example,

$$
V_{n,\mathrm{car}}=\alpha_{\mathrm{car}}+\beta_TT_{n,\mathrm{car}}+\beta_CC_{n,\mathrm{car}}+\gamma\,\mathrm{CarOwnership}_n.
$$

## 2.5 Availability and the choice set

MNL probabilities are calculated only over the alternatives assumed available to a traveller:

$$
P_{ni}=\frac{e^{V_{ni}}}{\sum_{j\in\mathcal C_n}e^{V_{nj}}}
\quad\text{if }i\in\mathcal C_n;
\qquad P_{ni}=0\quad\text{otherwise}.
$$

The choice set is not necessarily identical for all people. For example, an itinerary may be infeasible at a particular departure time. Nevertheless, **having no driving licence does not necessarily make every car trip unavailable**, since travelling as a passenger may remain possible. Availability rules should reflect the definitions of alternatives and the observed travel situation, not merely a convenient shortcut.

## 2.6 Independence of irrelevant alternatives (IIA)

For any two available alternatives $i$ and $j$,

$$
\frac{P_{ni}}{P_{nj}}=e^{V_{ni}-V_{nj}}.
$$

The ratio depends only on their own utilities, not on attributes or the number of other alternatives. This is the **IIA property** of standard MNL.

### A simple example: adding a similar option

Suppose two alternatives, route A and route B, have equal systematic utility. MNL assigns each probability $0.5$. Now a third alternative C is introduced, with the same systematic utility as A and B. MNL assigns each $1/3$.

If C is merely a **duplicate representation** of B in the analyst's choice list, this changes the predicted demand for A without changing the traveller's actual opportunities. That is the classic red-bus/blue-bus concern.

If C is instead a **real new service**, the choice opportunities may genuinely improve—through frequency, waiting time, reliability, or access. A change in the overall share of B and C is then not automatically wrong. The modelling question is whether the alternatives' similarities and service effects have been represented appropriately.

IIA is a statement about the *model's substitution structure*, not a claim that real travellers lack habits or that similar routes must have a fixed combined market share.

## 2.7 Marginal effects: how probabilities respond to utility

Because $P_{ni}=e^{V_{ni}}/\sum_j e^{V_{nj}}$, differentiation yields

$$
\frac{\partial P_{ni}}{\partial V_{ni}}=P_{ni}(1-P_{ni})
$$

and, for $j\ne i$,

$$
\frac{\partial P_{nj}}{\partial V_{ni}}=-P_{nj}P_{ni}.
$$

If an attribute $x_{ni}$ enters only the utility of alternative $i$ linearly, with coefficient $\beta_x$, then

$$
\boxed{\frac{\partial P_{ni}}{\partial x_{ni}}=\beta_xP_{ni}(1-P_{ni})}
$$

and

$$
\boxed{\frac{\partial P_{nj}}{\partial x_{ni}}=-\beta_xP_{nj}P_{ni},\qquad j\ne i.}
$$

If $\beta_x<0$ for travel time or cost, increasing that attribute lowers the own-choice probability and raises the other probabilities. The formula is a local derivative: a large policy change should be evaluated by recalculating probabilities at the new attributes.

### The proportional substitution implication

For any change affecting only alternative $i$, the MNL odds among the **unchanged** alternatives remain constant. In the infinitesimal case, the changes in their probabilities are proportional to their existing shares. This is one reason MNL can misrepresent substitution between similar services.

## 2.8 Direct and cross elasticities

An elasticity expresses the percentage change in an outcome associated with a **1% change** in an explanatory variable, locally:

$$
E_{P_i,x_i}=\frac{\partial P_i}{\partial x_i}\frac{x_i}{P_i}.
$$

For an attribute appearing only in the utility of alternative $i$ with coefficient $\beta_x$,

$$
\boxed{E_{P_i,x_i}=\beta_x x_i(1-P_i).}
$$

For an alternative $j\ne i$,

$$
\boxed{E_{P_j,x_i}=-\beta_x x_iP_i.}
$$

**Example.** Let PT have probability $P_{\mathrm{PT}}=0.40$, fare $x_{\mathrm{PT}}=£2$, and cost coefficient $\beta_C=-0.20$ per pound. Then its direct fare elasticity is

$$
E_{P_{\mathrm{PT}},C_{\mathrm{PT}}}=(-0.20)(2)(1-0.40)=-0.24.
$$

Thus, *for a sufficiently small percentage increase in fare*, the model predicts an approximately proportional reduction in PT choice probability, with an elasticity of $-0.24$.

For a finite change—for example, a 10% fare increase—a simple baseline-referenced measure is

$$
E^{\mathrm{finite}}\approx
\frac{(P_i'-P_i)/P_i}{(x_i'-x_i)/x_i}.
$$

This is a **finite-change percentage ratio**, not the exact point elasticity or the midpoint arc elasticity. If fares differ across travellers, aggregate sensitivity should be calculated from the model-predicted aggregate shares after applying a well-defined scenario. It should not be obtained by blindly substituting average values into individual elasticity formulas.

## 2.9 Value of travel time savings (VTTS)

If travel time is measured in hours and cost in pounds, with both entering utility linearly using common coefficients, then the monetary value of a one-hour travel-time saving is

$$
\boxed{VTTS=\frac{\beta_T}{\beta_C}\quad\text{GBP/hour}.}
$$

Both coefficients are usually negative, so the ratio is positive. For example, if $\beta_T=-6$ per hour and $\beta_C=-0.2$ per pound,

$$
VTTS=\frac{-6}{-0.2}=£30/\text{hour}.
$$

This value represents a marginal willingness to pay for travel-time savings **within the assumed model specification**. It is not a wage rate, and it may vary across people or modes if the model allows heterogeneous coefficients.

## 2.10 Numerical stability

When utilities are very large, direct computation of $e^{V_i}$ can overflow; when very negative, it can underflow. Subtracting the maximum utility before exponentiation solves the common problem without changing the probabilities:

$$
P_{ni}=\frac{e^{V_{ni}-m_n}}{\sum_j e^{V_{nj}-m_n}},
\quad m_n=\max_j V_{nj}.
$$

For likelihood calculations, the same principle leads to a stable **log-sum-exp** expression:

$$
\ln P_{n,y_n}=V_{n,y_n}-\ln\!\left(\sum_j e^{V_{nj}}\right).
$$

An equivalent stabilized denominator is

$$
\ln\!\left(\sum_j e^{V_{nj}}\right)
=m_n+\ln\!\left(\sum_j e^{V_{nj}-m_n}\right).
$$

This is a computational technique, not a behavioural assumption.

## 2.11 What MNL tells us—and what it does not

MNL offers a coherent probabilistic description of discrete decisions and interpretable trade-offs between observed attributes. It is often a useful starting point and benchmark. Its standard form, however, imposes the i.i.d. error assumption and IIA substitution structure. It also represents all travellers with common taste coefficients unless heterogeneity is explicitly introduced through interactions or segmentation.

A statistically significant coefficient indicates evidence of association under a model and inference assumptions; it does not by itself establish causality. Likewise, matching observed aggregate shares does not establish accurate individual choice predictions or trustworthy responses to future policies.

## 2.12 Concept checks (with answers)

**Q1. Why must one ASC be normalized?**  
Because adding the same constant to all utilities leaves the probabilities unchanged. Without normalization, the absolute level of all ASCs is unidentified.

**Q2. If the car-time coefficient is negative, does a 10-minute increase in car time always reduce car choice probability by the same number of percentage points?**  
No. The utility change can be the same, but the probability change depends on baseline probabilities and competing alternatives.

**Q3. Can a traveller without a driving licence choose a car trip?**  
Yes, depending on how “car” is defined; they might travel as a passenger. Licence status alone need not determine availability.

**Q4. If a near-duplicate route is added, does MNL necessarily predict an unrealistic change?**  
Not always. If the route is a real new service, it can improve access or waiting times. The concern arises when modelled alternatives are nearly duplicate representations without corresponding independent improvements in actual opportunities.

**Q5. Is a direct cost elasticity of $-0.3$ the same as a 30-percentage-point decline in mode share after a 1% fare increase?**  
No. It is a local proportional relationship: roughly a 0.3% reduction in choice probability for a 1% fare increase, holding the other conditions fixed.

**Q6. What changes in a standard MNL when the errors' common scale changes?**  
The strength of systematic utility relative to unobserved variation changes. Only coefficients divided by the scale can be separately identified from choice probabilities.

---

### Summary

MNL converts systematic utilities into choice probabilities through an i.i.d. Gumbel error assumption. Its major strengths are clarity, tractability, and interpretable trade-offs; its central structural restriction is IIA. Sound behavioural interpretation requires attention to availability, units, the choice set, and the distinction between utility changes, probability changes, elasticities, and causal effects.

**Next chapter:** Maximum Likelihood Estimation — how observed choices are used to estimate the utility parameters.
