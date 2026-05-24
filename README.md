# Dependent Pitman-Yor Process Mixture

This project studies a Bayesian nonparametric mixture model for two related
populations. The goal is to estimate two possibly different densities while still
allowing them to share latent structure. This is useful when groups are not
identical, but are expected to have related clustering behavior.

The model is based on a dependent Pitman-Yor process with a beta-product
dependent stick-breaking construction. The implementation focuses on the
symmetric H1 specification described in Bassetti, Casarin, and Leisen (2014),
combined with an ANOVA-style hierarchical structure for the mixture atoms.

## Motivation

Standard Dirichlet and Pitman-Yor process mixtures are powerful priors for
unknown densities, but in their basic form they model one random distribution at
a time. In grouped or repeated-measure settings, this is often too limited:
each group may have its own distribution, but the groups should not be treated
as completely unrelated.

The dependent construction used here addresses this by defining two random
distributions, `G1` and `G2`, whose mixture weights and atom locations are
correlated. This induces borrowing of strength between groups while preserving
group-specific heterogeneity.

## Dirichlet Process vs Pitman-Yor Process

The Pitman-Yor process generalizes the Dirichlet process by adding a discount
parameter. In stick-breaking form, the random weights are generated from beta
variables whose parameters depend on this discount.

When the discount parameter is zero, the Pitman-Yor process reduces to the
Dirichlet process. When the discount is positive, the model produces heavier
tails in the cluster-size distribution. This is important when the data contain
many small latent components, rare subpopulations, or long-tail behavior.

In practical terms:

- the Dirichlet process tends to produce a smaller number of dominant clusters;
- the Pitman-Yor process can support more small clusters;
- the Pitman-Yor process is better suited to power-law or heavy-tail structure.

## Dependent Stick-Breaking Construction

For each group, the random distribution is represented as an infinite mixture:

```text
Gi = sum_k Wik delta_{psi_ik}
```

The dependence between groups is introduced through the weights. In the H1
construction, the stick-breaking terms are:

```text
S1k = V0k V1k
S2k = V0k V2k
```

The variable `V0k` is shared across groups and acts as a common latent factor.
The variables `V1k` and `V2k` are group-specific. This creates dependence
between the weights `W1k` and `W2k`: if a component is important in one group,
the model encourages the corresponding component to also be relevant in the
other group.

This is the first level of dependence in the model: shared complexity.

## ANOVA Structure for Atoms

The second level of dependence concerns the atoms of the mixture. Each Gaussian
kernel has group-specific parameters, but those parameters are built from common
and group-specific components:

```text
mu_ik = mu_0k + mu_specific,ik
sigma^2_ik = sigma^2_0k sigma^2_specific,ik
```

The common component represents the baseline position and scale of the cluster.
The group-specific component allows each population to deviate from that
baseline.

This ANOVA-style hierarchy makes the model interpretable:

- common atoms capture structure shared by both groups;
- specific atoms capture group-level shifts;
- the prior controls how strongly clusters are aligned across groups.

This is the second level of dependence: shared location.

## Inference

The theoretical model has infinitely many mixture components. For computation,
the process is approximated with a finite truncation level `K_max`. The final
stick is closed so that the truncated weights sum to one.

Posterior inference is performed by MCMC. The sampler alternates between:

- updating cluster allocations for observations;
- updating common and group-specific atom parameters;
- updating dependent stick-breaking variables;
- updating the hyperparameters controlling dependence and Pitman-Yor behavior.

The atom updates use conjugate Normal-Inverse-Gamma steps. The stick-breaking
variables and hyperparameters are updated with Metropolis-Hastings steps.

## Interpretation

The model separates two related ideas that are often conflated:

- whether the groups use the same number and importance of latent components;
- whether those components are located in the same places.

The dependent weights answer the first question. If posterior weights for the
same component index are highly correlated across groups, this suggests that the
groups share latent complexity.

The ANOVA atom structure answers the second question. If posterior component
means and variances are similar across groups, this suggests that the groups
share not only complexity, but also physical cluster location.

Together, these two mechanisms allow the model to represent:

- nearly identical group distributions;
- distributions with shared clusters but shifted locations;
- distributions with similar global complexity but different local behavior;
- heavy-tail settings with many small components.

## Simulation Study

The project reproduces four synthetic mixture scenarios inspired by the
reference paper. These scenarios test the ability of the model to recover:

- balanced mixtures with common components;
- mixtures with unequal weights across groups;
- mixtures with shifted component locations;
- mixtures with different numbers of active components across groups.

Posterior mean density plots are used to assess density recovery. Scatter plots
of posterior weights are used to inspect dependence between groups.

## Real-Data Application

The empirical example uses IBM credit-card transaction data. The model is applied
to transaction amounts from two consecutive years, 2008 and 2009. This setting is
natural for dependent random measures: spending behavior is expected to be
correlated across adjacent years, but the distributions may still differ because
of economic changes and individual spending variation.

The model is especially relevant here because transaction amounts are
multimodal, include refunds, and may contain rare high-value transactions. These
features make a flexible nonparametric mixture more appropriate than a simple
unimodal parametric model.

## Practical Use

The main report is `NP_Favagrossa.pdf`. The R implementation is
`BNP_Favagrossa.R`.

To load the model functions without running the full analysis:

```r
source("BNP_Favagrossa.R")
```

To run the complete set of simulations and applications:

```sh
RUN_FULL_ANALYSIS=true Rscript BNP_Favagrossa.R
```

Required R packages:

```r
install.packages(c("ggplot2", "gridExtra", "tidyverse"))
```

The real-data section expects the file:

```text
IBM Credit Data/User0_credit_card_transactions.csv
```

inside the project folder.

