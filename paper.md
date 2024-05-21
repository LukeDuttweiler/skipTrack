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

Mobile apps that allow users to self-track menstrual cycle lengths and symptoms are now widely available and frequently used [@fox2010mobile]. Multiple studies (consider [@bull2019real; @li2020characterizing;@mahalingaiah2022design]) have taken advantage of these uniquely large data sets to gain insight on characteristics of the menstrual cycle, which is an important vital sign [@diaz2006menstruation]. Unfortunately, due to the self-tracking nature of the gathered data, recorded cycle lengths may be inflated if users do not accurately document all period dates in the app. A non-trivial number of incorrectly inflated cycle lengths in a data set will be damaging to the reliability and reproducibility of analysis results. 

Current solutions to this problem of non-adherence in cycle tracking include removing cycles that exhibit no user-app interaction [@li2020characterizing], identifying possibly inaccurate cycles based on user-specific average cycle lengths [@li2022predictive], or *ad hoc* removal of cycles based on well-established menstrual cycle characteristics. The `skipTrack` package adapts and advances the Bayesian approach of @li2022predictive by identifying possible skips in cycle tracking based on user-specific average cycle lengths **and** user-specific cycle regularity. 

# Statement of need

Analyses involving large amount of user-tracked menstrual cycle data sets are becoming more prevalent. Identifying skips in cycle tracking is crucially important for maintaining the validity of these studies. The `skipTrack` package provides easy to use software in R that can identify skips in menstrual cycle data based on a pre-specified Bayesian hierarchical model. The resulting inference on possible skipped cycles may then be included by a researcher *a priori* in an analysis, or may be used to develop a multiple-imputation scheme. 

Additionally, while based on the Bayesian hierarchical model from [@li2022predicitve], the model used by `skipTrack` includes components for both cycle length mean and regularity. This allows the model to correctly adjust for individuals with irregular cycles who are often excluded from menstrual cycle analyses, despite the important information their data contains. 
Finally, several extensions to the current `skipTrack` model and software are planned. These include the addition of regression models for both cycle length mean and regularity, an auto-regressive modeling structure for sequential cycle lengths from the same individual, and a method for the inclusion of user-app interaction or other external data to help with skip identification. These updates, along with open availability and ease-of-use, will provide researchers easy access to high level modeling techniques for mobile menstrual cycle data.

# The SkipTrack Model

We present a short overview of the SkipTrack model and notation here. 

Let $y_{ij}$ be the $j$th recorded cycle length provided by participant $i$. We assume that 

\[
y_{ij} \sim \text{LogNormal}\Big(\mu_i + \log(c_{ij}), \tau_i\Big)
\]

where $\mu_i$ is the natural log of individual $i$'s cycle length median, $\tau_i$ is the precision of the distribution (providing a measure of regularity), and $c_{ij}$ is an integer-valued parameter that represents the number of **true** cycles occurring in recorded cycle $y_{ij}$. For example, if $c_{ij} = 1$, then $y_{ij}$ is a true cycle length, if $c_{ij} = 2$ then $y_{ij}$ is the length of two true cycles added together, and so on. 

Then we assume,

\[
\mu_i \sim \text{Normal}(\mu, \rho) \mspace{100mu} \tau_i \sim \text{Gamma}(\theta, \phi)
\]

where the natural log of $\mu$ gives the overall population cycle length median, $\rho$ is a precision parameter, $\theta$ is the mean of the Gamma distribution and $\phi$ is the rate. 

Finally, 

\[
c_{ij} \sim \text{Categorical}(\pi_1, \pi_2, \dots, \pi_{NS})
\]

where $pi_k = \text{Pr}(c_{ij} = k)$ and $NS$ is the maximum number of skips allowed in the model. 

# Package Description

The `skipTrack` package contains tools for fitting the SkipTrack model, visualizing model results, diagnosing model convergence, and simulating example data. 

The model fit is accomplished using an MCMC algorithm composed mainly of Gibbs sampling steps with a small number of Metropolis-Hastings steps. Model fitting is accomplished through an easy-to-use interface that allows users to select the number of MCMC chains to run, the number of iterations to run per chain, and the parameters used to initialize each chain. Model results may be visualized or retrived through standard interaction functions (`summary()`, `plot()`, etc.).

MCMC convergence diagnostics are multivariate and multi-chain and are provided using the R package @genMCMCPackage. 

Example data may be simulated from the SkipTrack model, the generative model provided in @li2022predictive, or a provided mixture model. 

# Availability 

A stable version of `skipTrack` is available on CRAN, and a development version is publicly available on GitHub (https://github.com/LukeDuttweiler/skipTrack).

# Acknowledgements

Research reported in this publication was supported by the National Institute of Environmental Health Sciences of the National Institutes of Health (NIH) under award number T32ES007142.  The content is solely the responsibility of the authors and does not necessarily represent the official views of the NIH.

(Shruthi, Brent, grant info?)

# References
