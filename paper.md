---
title: 'skipTrack: An R package for Identifying Skips in Self-Tracked Mobile Menstrual Cycle Data'
tags:
  - menstrual cycle
  - Bayesian hierarchical model
  - MCMC
authors:
  - name: Luke Duttweiler
    orcid: 0000-0002-0467-995X
    affiliation: 1 # (Multiple affiliations must be quoted)
    corresponding: true
  - name: Shruthi Mahalingaiah
    orcid: 0000-0002-5527-5787
    affiliation: 2
  - name: Brent Coull
    orcid: 0000-0002-1808-4156
    affiliation: 1
affiliations:
 - name: Department of Biostatistics, Harvard T.H. Chan School of Public Health, Boston, MA, United States
   index: 1
 - name: Department of Epidemiology, Harvard T.H. Chan School of Public Health, Boston, MA, United States
   index: 2
date: 21 May 2024
bibliography: bibliography.bib
---

# Summary

Mobile apps that allow users to self-track menstrual cycle lengths and symptoms are now widely available and frequently used [@fox2010mobile]. Multiple studies (consider [@bull2019real; @li2020characterizing;@mahalingaiah2022design]) have taken advantage of these uniquely large data sets to gain insight on characteristics of the menstrual cycle, which is an important vital sign [@diaz2006menstruation]. Due to the self-reported nature of the gathered data, recorded cycle lengths may be inflated if users skip tracking any cycle related bleeding events in the app. A non-trivial number of incorrectly inflated cycle lengths in a data set will be damaging to the reliability and reproducibility of analysis results. 

Current solutions to this problem of non-adherence (skipped tracking) in cycle length reporting include removing implausibly long cycles that exhibit no user-app interaction [@li2020characterizing], identifying possibly inaccurate cycles based on user-specific average cycle lengths [@li2022predictive], or *ad hoc* removal of cycles based on well-established menstrual cycle characteristics such as average cycle length or cycle length difference. The `skipTrack` package implements a Bayesian hierarchical model that is the first to explicitly use information on both an individual's cycle length **and** regularity to identify errors in recorded cycle lengths that arise from user non-adherence in logging one or more bleeding events.

# Statement of Need

Analyses involving large amount of user-tracked menstrual cycle data sets are becoming more prevalent. Identifying recorded cycle lengths that result from skips in tracking one or more period bleeding events (hereafter referred to as 'skipped cycles') is crucially important for maintaining the validity of these studies. The `skipTrack` package provides easy to use software in R that can identify skipped cycles in menstrual cycle data based on a pre-specified Bayesian hierarchical model. The resulting inference on possible skipped cycles may then be included by a researcher *a priori* in an analysis, or may be used to develop a multiple-imputation scheme. 

Additionally, while based on the Bayesian hierarchical model from [@li2022predictive], the model used by `skipTrack` includes components for both cycle length mean and regularity. This allows the model to correctly adjust for individuals with irregular cycles who are often excluded from menstrual cycle analyses, despite the important information their data contains. 
Finally, the `skipTrack` model and software lead to many possible useful extensions including the addition of regression models for both cycle length mean and regularity, an auto-regressive modeling structure for sequential cycle lengths from the same individual, and a method for the inclusion of user-app interaction or other external data to help with skip identification. These updates, along with open availability and ease-of-use, will provide researchers easy access to high level modeling techniques for mobile menstrual cycle data.

# The SkipTrack Model

We present a short overview of the SkipTrack model and notation here. 

Let $y_{ij}$ be the $j$th recorded cycle length provided by participant $i$. We assume that 

$$
y_{ij} \sim \text{LogNormal}\Big(\mu_i + \text{log}(c_{ij}), \tau_i\Big)
$$

where $\mu_i$ is the natural log of individual $i$'s cycle length median, $\tau_i$ is the precision of the distribution (providing a measure of regularity), and $c_{ij}$ is an integer-valued parameter that represents the number of **true** cycles occurring in recorded cycle $y_{ij}$. For example, if $c_{ij} = 1$, then $y_{ij}$ is a true cycle length, if $c_{ij} = 2$ then $y_{ij}$ is the length of two true cycles added together, and so on. 

Then we assume,

$$
\mu_i \sim \text{Normal}(\mu, \rho) \mspace{100mu} \tau_i \sim \text{Gamma}(\theta, \phi)
$$

where the natural log of $\mu$ gives the overall population cycle length median, $\rho$ is a precision parameter, $\theta$ is the mean of the Gamma distribution and $\phi$ is the rate. 

Finally, 

$$
c_{ij} \sim \text{Categorical}(\pi_1, \pi_2, \dots, \pi_{K})
$$

where $\pi_k = \text{Pr}(c_{ij} = k)$ and $K$ is the maximum number of skips allowed in the model. 

# Package Description

The `skipTrack` package contains tools for fitting the SkipTrack model, visualizing model results, diagnosing model convergence, and simulating example data. 

 - **Model Fitting**: In order to fit the SkipTrack model, the code employs a Markov Chain Monte Carlo (MCMC) algorithm composed of Gibbs sampling steps. Model fitting may be accessed through the `skipTrack.fit()` function, and is accomplished through this easy-to-use interface that allows users to select the number of MCMC chains to run, the number of iterations to run per chain, and the parameters used to initialize each chain. 
 
 - **Visualizing Results**: Model results may be visualized or retrieved through standard reporting and visualization functions (`summary()`, `plot()`, etc.).

 - **Diagnosing Convergence**: MCMC convergence diagnostics (traceplots, effective sample size, and the Gelman-Rubin potential scale reduction factor) are multivariate and multi-chain and are provided using the R package `genMCMCDiag` [@genMCMCPackage], accessible through `skipTrack.diagnostics()`.

 - **Simulating Data**: Data simulation options using `skipTrack.simulate()` are included which allow a user to simulate example data from the SkipTrack model, the generative model provided in @li2022predictive, or a provided mixture model. 

# Availability 

A stable version of `skipTrack` is available on CRAN, and a development version is publicly available on GitHub (https://github.com/LukeDuttweiler/skipTrack).

# Acknowledgements

Research reported in this publication was supported by the National Institute of Environmental Health Sciences (NIEHS) Grants T32ES007142, P30 ES000002, and R01 ES035106.  The content is solely the responsibility of the authors and does not necessarily represent the official views of the NIH.

# References
