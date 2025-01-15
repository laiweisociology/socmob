# socmob
An R package that implements the proposed methods in the research article "Social Mobility as Causal Intervention" by Lai Wei and Yu Xie.

The package is at an early developmental phase so comments, feedbacks, and reporting of BUGs would be most appreciated!

The estimation procedures of the package is based on the "sl3" environment. To install the package, run:

```
remotes::install_github("tlverse/sl3")
remotes::install_github("laiweisociology/socmob")
```

A list of functions is below:

calculate_odds_ratio(): calculating the odds ratios of a given table
stepstable(): generating a counterfactual mobility table by defining the strength of an intervention that drags OD association toward independence
contrasttable(): generating a counterfactual mobility table by incorporating the odds-ratios from a comparison table while maintaining the observed marginal distributions unchanged
pips(): estimating the post-intervention propensity scrore of ending up in specific destinations, given a counterfactual table and the conditional odds ratio preservation assumption
conditionaleffect(): estimating the conditional destination effect by origin
socmob(): estimating social mobility effect, given the outputs from pips() and conditionaleffect()
