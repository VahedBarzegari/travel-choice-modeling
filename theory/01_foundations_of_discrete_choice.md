# Chapter 1 — Foundations of Discrete Choice Modelling

*Discrete Choice Modelling: Theory, Mathematics, and Behavioural Interpretation*

## Learning objectives

By the end of this chapter, you should be able to:

1. Explain what a discrete choice model attempts to describe and distinguish a choice from a choice probability.
2. Define the decision maker, alternatives, choice set, attributes, and observed choice.
3. State the random utility model and explain why the modeller treats part of utility as unobserved.
4. Derive the probability of choosing an alternative as a comparison of utilities.
5. Interpret systematic utility, alternative-specific constants, coefficient signs, and trade-offs.
6. Recognize the roles of utility normalization, scale, availability, and repeated choices.

---

## 1.1 What is a discrete choice?

A **discrete choice** is a selection from a finite set of distinct alternatives. Examples include choosing a mode of travel (walking, cycling, public transport, or car), a departure time (morning, afternoon, or evening), or a route between two locations.

Discrete choice modelling seeks to explain and predict **which alternative is chosen**, typically by estimating the probability of each available option as a function of its attributes and the characteristics of the decision maker.

This is different from predicting a continuous quantity such as travel time or the number of boardings on a bus route. A mode-choice model starts with alternatives and a decision, rather than directly predicting a numerical outcome.

### Essential notation

Let:

- $n=1,\ldots,N$ index decision-making observations (often trips, made by individuals).
- $i$ index a possible alternative.
- $\mathcal C_n$ be the **choice set** containing alternatives available to observation $n$.
- $y_n\in\mathcal C_n$ be the alternative actually chosen.
- $x_{ni}$ denote measured attributes of alternative $i$ for observation $n$.
- $z_n$ denote characteristics of the decision maker, such as age or household vehicle ownership.

For a commuter, $x_{ni}$ might include time, fare, and transfers for each mode. The same commuter's $z_n$ remains unchanged across alternatives for that decision.

**Important distinction:** An option can exist in the transportation system without being available to a particular person. For example, cycling may not be feasible for some trips. Availability and preference are not the same thing.

---

## 1.2 Utility and the decision rule

**Utility** is a numerical representation of the desirability of an alternative to a decision maker. In the standard random utility framework, the decision maker is assumed to select the available alternative with the highest utility:

$$
y_n=i \quad\text{if}\quad U_{ni}\ge U_{nj}\quad \text{for all } j\in\mathcal C_n.
$$

Here, $U_{ni}$ is the utility of alternative $i$ for observation $n$.

Two cautions matter:

1. Utility is a **model construct**; it is not usually a directly observed psychological score.
2. A larger utility indicates a more preferred alternative **relative to competing alternatives**, not necessarily an objectively good or pleasant journey.

A traveller may rationally choose a slower service because it is cheaper, more comfortable, easier to access, or available when the traveller needs it.

---

## 1.3 Random utility theory

The modeller cannot observe every factor affecting a decision. We therefore write:

$$
\boxed{U_{ni}=V_{ni}+\varepsilon_{ni}}
$$

where:

- $V_{ni}$ is **systematic (observed or explainable) utility**, computed from variables and parameters in the model.
- $\varepsilon_{ni}$ is the **unobserved component** of utility, treated as random by the modeller.

The random component can reflect omitted characteristics, imperfect measurements, variation in unobserved preferences, and other influences that are not specified in $V_{ni}$. Calling it random **does not mean travellers choose randomly**. It means the modeller has incomplete information about their utilities.

A simple specification is:

$$
V_{ni}=\alpha_i+\beta_T T_{ni}+\beta_C C_{ni},
$$

with:

- $\alpha_i$: alternative-specific constant (ASC),
- $T_{ni}$: travel time,
- $C_{ni}$: monetary cost,
- $\beta_T$: coefficient on time,
- $\beta_C$: coefficient on cost.

Typically, $\beta_T<0$ and $\beta_C<0$, other things equal. The coefficient magnitudes depend on the units of measurement and the model's utility scale.

### Why not simply pick the smallest travel time?

Because time is only one component of utility. For two alternatives $a$ and $b$:

$$
V_{na}-V_{nb}
=(\alpha_a-\alpha_b)
+\beta_T(T_{na}-T_{nb})
+\beta_C(C_{na}-C_{nb}).
$$

Even if $T_{na}>T_{nb}$, option $a$ can have greater utility if its cost or other attributes compensate for the time difference.

---

## 1.4 From utility to choice probability

The choice rule involves utilities that are partly unobserved to the modeller. Consequently, the modeller cannot generally know the chosen alternative with certainty.

The probability that alternative $i$ is chosen is:

$$
\boxed{P_{ni}=\Pr\!\left(U_{ni}\ge U_{nj}\ \forall j\in\mathcal C_n\right).}
$$

Substituting the random utility decomposition gives:

$$
P_{ni}
=\Pr\!\left(
\varepsilon_{nj}-\varepsilon_{ni}
\le V_{ni}-V_{nj}
\quad\forall j\in\mathcal C_n
\right).
$$

This expression contains the central insight of discrete choice modelling:

> Choice probabilities are determined by **differences in systematic utility** and by assumptions about the **joint distribution of unobserved utilities**.

The probability formula therefore depends on the assumed distribution of $\boldsymbol\varepsilon_n$.

For example, if the errors are independent and identically distributed Type-I extreme value with a common scale, we obtain the familiar **multinomial logit (MNL)** formula (derived in Chapter 2):

$$
P_{ni}=\frac{\exp(V_{ni})}{\sum_{j\in\mathcal C_n}\exp(V_{nj})},
$$

where the error scale has been normalized to one. Other error assumptions lead to other choice models.

---

## 1.5 Worked numerical example

Suppose one traveller considers three modes. After accounting for the measured attributes, their systematic utilities (in normalized utility units) are:

| Alternative | Systematic utility $V_i$ |
|---|---:|
| Walk | $-1.5$ |
| Public transport | $-1.0$ |
| Car | $-2.0$ |

Public transport has the highest systematic utility, but that does **not** imply it will be chosen with probability 1.

For illustration, assume the MNL error distribution. The probabilities are:

$$
P_{\mathrm{walk}}=
\frac{e^{-1.5}}{e^{-1.5}+e^{-1.0}+e^{-2.0}}
\approx 0.3072,
$$

$$
P_{\mathrm{PT}}=
\frac{e^{-1.0}}{e^{-1.5}+e^{-1.0}+e^{-2.0}}
\approx 0.5065,
$$

$$
P_{\mathrm{car}}=
\frac{e^{-2.0}}{e^{-1.5}+e^{-1.0}+e^{-2.0}}
\approx 0.1863.
$$

| Alternative | Probability |
|---|---:|
| Walk | 30.72% |
| Public transport | 50.65% |
| Car | 18.63% |
| **Total** | **100.00%** |

**Interpretation:** Public transport is the most probable choice, not a certain choice. These percentages are the model's probabilistic predictions, not an observed division of travellers.

### A small improvement in public transport

Suppose public transport becomes more attractive and its systematic utility increases by $0.5$, from $-1.0$ to $-0.5$, while the others remain unchanged. Its probability becomes:

$$
P'_{\mathrm{PT}}
=\frac{e^{-0.5}}{e^{-1.5}+e^{-0.5}+e^{-2.0}}
\approx 0.6285.
$$

The improved option becomes more likely to be chosen. *How* the shares of the other alternatives change depends on the substitution structure of the selected model.

---

## 1.6 Alternative-specific constants

An **alternative-specific constant**, $\alpha_i$, represents the baseline level of an alternative's systematic utility relative to a reference alternative, conditional on the included attributes.

With three modes, we might write:

$$
\begin{aligned}
V_{n,\mathrm{walk}} &= \beta_T T_{n,\mathrm{walk}},\\
V_{n,\mathrm{PT}} &= \alpha_{\mathrm{PT}}+\beta_T T_{n,\mathrm{PT}}+\beta_C C_{n,\mathrm{PT}},\\
V_{n,\mathrm{car}} &= \alpha_{\mathrm{car}}+\beta_T T_{n,\mathrm{car}}+\beta_C C_{n,\mathrm{car}}.
\end{aligned}
$$

Here the walking ASC is normalized to zero. This is a **reference convention**, not an assumption that walking has zero utility.

A negative car ASC does **not** mean that all travellers dislike cars. Car may still have the highest total systematic utility for many journeys.

### Individual characteristics enter through alternatives

A characteristic such as a driving licence is constant across the alternatives considered by the *same* traveller. Adding $\gamma\,\mathrm{Licence}_n$ to **every** utility has no effect on that traveller's choice probabilities: the term cancels from all utility differences.

To study whether licence status is associated with car choice, a possible specification is:

$$
V_{n,\mathrm{car}}
=\alpha_{\mathrm{car}}+\beta_T T_{n,\mathrm{car}}
+\beta_C C_{n,\mathrm{car}}
+\gamma\,\mathrm{Licence}_n.
$$

This is an **interaction with a specific alternative**, even though the variable itself is person-specific.

---

## 1.7 Identification and normalization

Choice data reveal **relative**, not absolute, utility.

### Additive normalization

If the same constant $k$ is added to every utility for one observation, the ranking does not change:

$$
(U_{ni}+k)-(U_{nj}+k)=U_{ni}-U_{nj}.
$$

In MNL, the probabilities are likewise unchanged:

$$
\frac{e^{V_{ni}+k}}{\sum_j e^{V_{nj}+k}}
=\frac{e^{V_{ni}}}{\sum_j e^{V_{nj}}}.
$$

Therefore, one ASC must ordinarily be fixed (often to zero); otherwise, all ASCs cannot be identified separately.

### Utility scale

The magnitude of utility is not inherently measured in pounds or minutes. In a random utility model, multiplying both systematic utilities and the random errors by the same positive constant does not change the choices. Thus, the overall utility scale must be normalized by an assumption on the error distribution.

For MNL, this is usually handled by fixing the Type-I extreme-value error scale. **Raw coefficient magnitudes from models with different scale normalizations are not automatically comparable.** Ratios of coefficients sharing the same scale, such as a time-to-cost trade-off, can be more interpretable when well defined.

---

## 1.8 Choice sets and availability

The probability of choosing an option is defined over its available choice set $\mathcal C_n$, not automatically over every imaginable alternative.

For a four-mode application:

$$
\mathcal C_n\subseteq
\{\mathrm{walk},\mathrm{cycle},\mathrm{PT},\mathrm{car}\}.
$$

A distinction is useful:

- **Availability:** Is the option feasible for the decision maker and this decision?
- **Attributes:** What are its time, cost, transfers, and other measurable characteristics?
- **Preference:** How strongly does the decision maker value the alternative, given those attributes?

For example, household car ownership and possession of a driving licence may affect access to driving, but they do not necessarily tell us whether an individual can travel **as a car passenger**. Availability rules need to match the definition of the alternative.

---

## 1.9 Repeated observations and the unit of analysis

A dataset can contain multiple trips by the same person. Each trip may represent one decision, but the decisions need not be statistically independent: individual tastes, habitual choices, and shared circumstances can persist across trips.

It is useful to distinguish:

- **Observation:** One recorded choice (such as a trip).
- **Decision maker:** The person or household making that choice.
- **Panel:** Several observed choices from the same decision maker.

Later chapters will distinguish methods that **adjust uncertainty** for repeated observations (cluster-robust standard errors) from methods that **explicitly model** dependence and preference heterogeneity (panel mixed logit).

---

## 1.10 Key takeaways

1. Discrete choice models describe choices among a set of available alternatives.
2. Random utility theory separates utility into systematic and unobserved components.
3. The decision maker is assumed to choose the option with the highest total utility; the modeller predicts probabilities because total utility is only partly observed.
4. Choice probabilities depend on systematic **utility differences** and assumptions about unobserved utility.
5. Alternative-specific constants and coefficients must be interpreted within a clearly defined utility specification and reference convention.
6. Availability, individual characteristics, and repeated observations matter when translating a behavioural theory into an empirical model.

## Concept-check questions

1. Why can a slower travel mode still be chosen under utility maximization?
2. Does the random utility error imply irrational or random behaviour by travellers?
3. Why is one alternative-specific constant normally fixed to zero?
4. Why does adding the same person-specific characteristic to all alternative utilities fail to influence MNL probabilities?
5. What is the difference between an available alternative and an attractive alternative?
6. If two models have different utility scales, can their estimated time coefficients be compared directly? Why or why not?

### Short answers

1. Other attributes and unobserved influences can make its total utility greater.
2. No. Randomness represents uncertainty from the modeller's perspective.
3. Choice data identify only utility differences, so a reference ASC is required.
4. It adds the same value to every utility and cancels from the differences.
5. Availability concerns feasibility; attractiveness concerns relative utility.
6. Not necessarily. The scale normalization affects the coefficient magnitudes.

---

**Next chapter:** *Multinomial Logit: assumptions, choice probabilities, identification, and the independence of irrelevant alternatives (IIA).* 
