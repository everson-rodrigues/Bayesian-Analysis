# Bayesian-Analysis
This code uses the Bayesian statistics to look for a set of parameters for the dark matter.
Dark matter is believed to be one of the main constituents of Universe, with its total mass beeing greater than the massa of the visible matter (all that you see is just 5% of the universe).
In this code I use the package pymultinest to evaluate a Bayesian analysis looking for the most probable physical parameters that describes fermionic dark model. 

How that can be applied to business? Suppose you have a handfull of variables and that generates models (e.g., churn prediction, pricing, credit risk). Bayesian analysis provides not only the most likely parameter values but also the uncertainty around them, enabling risk-aware decisions. Examples include predicting customer churn with confidence intervals, testing pricing strategies under market uncertainty, or assessing credit risk with probabilistic outputs.

Besides estimating the most probable parameters, Bayesian analysis also provides the Bayesian Evidence, which measures how well a model explains the data. This makes it possible to define which model is better given our prior knowledge. One key advantage is that Bayesian methods naturally penalize overly complex models with too many variables, rewarding those that achieve a balance between explanatory power and simplicity.

In short, the same statistical tools used here to understand the Universe can be applied to optimize strategies, manage risks, and make data-driven decisions in business.


What the code takes: It import your EOS and your Likelihood function and gives back the parameters. 
