# Dependent Pitman-Yor Process Mixture

This project studies a Bayesian nonparametric mixture for two related populations. The objective is to estimate two unknown distributions while allowing them to share latent clustering structure. Dependence is introduced both in the mixture weights and in the mixture atoms.

## Model

For groups $i=1,2$, observations are modeled as:

$$y_{ij} \mid z_{ij}=k \sim \mathcal{N}(\mu_{ik}, \sigma_{ik}^{2})$$

$$\Pr(z_{ij}=k) = \pi_{ik}$$

The random measure for group $i$ is:

$$G_i = \sum_{k=1}^{\infty} \pi_{ik} \delta_{\psi_{ik}}.$$

Here, $\pi_{ik}$ are random weights and $\psi_{ik}$ are the atoms of the mixture.

## Dependent Stick-Breaking

The model uses a beta-product dependent stick-breaking construction. In the symmetric $H_1$ case:

$$S_{1k} = V_{0k} V_{1k}, \qquad S_{2k} = V_{0k} V_{2k}.$$

The variable $V_{0k}$ is shared and induces dependence between the two groups. The weights are:

$$\pi_{ik} = S_{ik} \prod_{\ell < k}(1 - S_{i\ell}).$$

Equivalently:

$$\pi_{1k} = V_{0k} V_{1k} \prod_{\ell < k}(1 - V_{0\ell} V_{1\ell}),$$

$$\pi_{2k} = V_{0k} V_{2k} \prod_{\ell < k}(1 - V_{0\ell} V_{2\ell}).$$

The beta variables are:

$$V_{0k} \sim \mathrm{Beta}(1-c,\, \theta_1)$$

$$V_{1k} \sim \mathrm{Beta}(1-c+\theta_1,\, \theta_2 + c(k-1))$$

$$V_{2k} \sim \mathrm{Beta}(1-c+\theta_1,\, \theta_2 + c(k-1)).$$

The parameter $c$ is the Pitman-Yor discount:

$$c = 0 \quad\Rightarrow\quad \text{dependent Dirichlet-process behavior}$$

$$0 < c < 1 \quad\Rightarrow\quad \text{Pitman-Yor behavior with heavier tails}.$$

Thus, positive $c$ allows more small active clusters and is better suited to long-tail data.

## ANOVA Atom Structure

The atoms are decomposed into a shared component and a group-specific deviation:

$$\mu_{ik} = \mu_{0k} + \mu_{ik}^{*}$$

$$\sigma_{ik}^{2} = \sigma_{0k}^{2} \, \sigma_{ik}^{2*}.$$

The priors are:

$$\mu_{0k} \sim \mathcal{N}(0, s_0^2), \qquad \mu_{ik}^{*} \sim \mathcal{N}(0, s_i^2)$$

$$\sigma_{0k}^{2} \sim \mathrm{InvGamma}\left(\frac{\varepsilon}{2}, \frac{\varepsilon}{2}\right), \qquad \sigma_{ik}^{2*} \sim \mathrm{InvGamma}\left(\frac{\lambda}{2}, \frac{\lambda}{2}\right).$$

This hierarchy separates common structure from group-specific variation. The shared atom $(\mu_{0k}, \sigma_{0k}^2)$ aligns clusters across groups, while $(\mu_{ik}^{*}, \sigma_{ik}^{2*})$ allows local deviations.

## Posterior Allocation

Given weights and atoms, the allocation probability is:

$$\Pr(z_{ij} = k \mid {-}) \propto \pi_{ik} \, \mathcal{N}(y_{ij} \mid \mu_{ik}, \sigma_{ik}^{2}).$$

After normalization:

$$\Pr(z_{ij} = k \mid {-}) = \frac{\pi_{ik} \, \mathcal{N}(y_{ij} \mid \mu_{ik}, \sigma_{ik}^{2})}{\displaystyle\sum_{h=1}^{K} \pi_{ih} \, \mathcal{N}(y_{ij} \mid \mu_{ih}, \sigma_{ih}^{2})}.$$

The infinite mixture is approximated with a finite truncation level $K$.

## Posterior Updates

Conditional on allocations, atom parameters have conjugate Normal-Inverse-Gamma updates. For the group-specific mean:

$$\mu_{ik}^{*} \mid {-} \sim \mathcal{N}(m_{ik}, V_{ik})$$

with:

$$V_{ik}^{-1} = \frac{1}{s_i^2} + \frac{n_{ik}}{\sigma_{0k}^{2} \, \sigma_{ik}^{2*}},$$

$$m_{ik} = V_{ik} \frac{\displaystyle\sum_{j:\, z_{ij}=k}(y_{ij} - \mu_{0k})}{\sigma_{0k}^{2} \, \sigma_{ik}^{2*}}.$$

The common mean uses observations from both groups:

$$\mu_{0k} \mid {-} \sim \mathcal{N}(m_{0k}, V_{0k})$$

with:

$$V_{0k}^{-1} = \frac{1}{s_0^2} + \frac{n_{1k}}{\sigma_{0k}^{2} \, \sigma_{1k}^{2*}} + \frac{n_{2k}}{\sigma_{0k}^{2} \, \sigma_{2k}^{2*}}.$$

This is the main borrowing-of-strength mechanism for the atoms.

## Stick and Hyperparameter Updates

The stick variables are updated with Metropolis-Hastings proposals on the logit scale:

$$\mathrm{logit}(V_k^{\text{new}}) = \mathrm{logit}(V_k^{\text{old}}) + \eta, \qquad \eta \sim \mathcal{N}(0, \tau^2).$$

The acceptance probability is:

$$
\alpha =
\min\left(
1,\;
\frac{p(V^{\mathrm{new}} \mid z, \theta_1, \theta_2, c)}
     {p(V^{\mathrm{old}} \mid z, \theta_1, \theta_2, c)}
\frac{J(V^{\mathrm{new}})}
     {J(V^{\mathrm{old}})}
\right).
$$

The hyperparameters $\theta_1$ and $\theta_2$ control concentration of the stick-breaking weights. The discount $c$ controls the departure from Dirichlet-process behavior and determines how strongly the model supports small clusters.

## Interpretation

The model captures two kinds of dependence.

**Dependence in cluster importance:** $\pi_{1k} \approx \pi_{2k}$ means that component $k$ has similar relevance in both groups.

**Dependence in cluster location:** $\mu_{1k} \approx \mu_{2k}$ and $\sigma_{1k}^{2} \approx \sigma_{2k}^{2}$ mean that the same latent component appears in similar regions of the sample space.

The dependent Pitman-Yor process is therefore able to represent shared clusters, shifted clusters, unequal weights, and heavy-tailed behavior with many small components.
