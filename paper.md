---
title: 'SkipTrack: An R package for Identifying Skips in Self-Tracked Mobile Menstrual Cycle Data'
tags:
  - menstrual cycle
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

Mobile apps that allow users to self-track menstrual cycle lengths and symptoms are now widely available and frequently used [@fox2010mobile]. Multiple studies (consider [@bull2019real; @li2020characterizing;@mahalingaiah2022design]) have taken advantage of these uniquely large data sets to gain insight on characteristics of the menstrual cycle, which is an important vital sign [@diaz2006menstruation]. Unfortunately, due to the self-tracking nature of the gathered data, recorded cycle lengths may be over-inflated as users fail to identify period dates in the app. A non-trivial number of incorrectly inflated cycle lengths in a data set will be damaging to the reliability and reproducibility of analysis results. 

Current solutions to this problem of non-adherence in cycle tracking include removing cycles that contain no user-app interaction [@li2020characterizing], identifying possibly inaccurate cycles based on user-specific average cycle lengths [@li2022predictive], or ad-hoc removal of cycles based on well-established menstrual cycle characteristics. The `SkipTrack` package adapts and advances the Bayesian approach of @li2022predictive by identifying possible skips in cycle tracking based on user-specific average cycle lengths *and* user-specific cycle regularity. 

# Statement of need

`Gala` is an Astropy-affiliated Python package for galactic dynamics. Python
enables wrapping low-level languages (e.g., C) for speed without losing
flexibility or ease-of-use in the user-interface. The API for `Gala` was
designed to provide a class-based and user-friendly interface to fast (C or
Cython-optimized) implementations of common operations such as gravitational
potential and force evaluation, orbit integration, dynamical transformations,
and chaos indicators for nonlinear dynamics. `Gala` also relies heavily on and
interfaces well with the implementations of physical units and astronomical
coordinate systems in the `Astropy` package [@astropy] (`astropy.units` and
`astropy.coordinates`).

`Gala` was designed to be used by both astronomical researchers and by
students in courses on gravitational dynamics or astronomy. It has already been
used in a number of scientific publications [@Pearson:2017] and has also been
used in graduate courses on Galactic dynamics to, e.g., provide interactive
visualizations of textbook material [@Binney:2008]. The combination of speed,
design, and support for Astropy functionality in `Gala` will enable exciting
scientific explorations of forthcoming data releases from the *Gaia* mission
[@gaia] by students and experts alike.

# Mathematics

Single dollars ($) are required for inline mathematics e.g. $f(x) = e^{\pi/x}$

Double dollars make self-standing equations:

$$\Theta(x) = \left\{\begin{array}{l}
0\textrm{ if } x < 0\cr
1\textrm{ else}
\end{array}\right.$$

You can also use plain \LaTeX for equations
\begin{equation}\label{eq:fourier}
\hat f(\omega) = \int_{-\infty}^{\infty} f(x) e^{i\omega x} dx
\end{equation}
and refer to \autoref{eq:fourier} from text.

# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

If you want to cite a software repository URL (e.g. something on GitHub without a preferred
citation) then you can do it with the example BibTeX entry below for @fidgit.

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

# Figures

Figures can be included like this:
![Caption for example figure.\label{fig:example}](figure.png)
and referenced from text using \autoref{fig:example}.

Figure sizes can be customized by adding an optional second parameter:
![Caption for example figure.](figure.png){ width=20% }

# Acknowledgements

We acknowledge contributions from Brigitta Sipocz, Syrtis Major, and Semyeong
Oh, and support from Kathryn Johnston during the genesis of this project.

# References
