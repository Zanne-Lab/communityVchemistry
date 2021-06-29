Community v chemistry
================
10/23/2017

## Overview

This repository contains code associated with the manuscript “Initial
wood trait variation overwhelms endophyte community effects for
explaining decay trajectories” (in review)

Author names, institutions and addresses: Marissa Lee 1 Jeff R. Powell 2
Brad Oberle 3 Faride Unda 4 Shawn D. Mansfield 4 Kylie Brice 2 Jessica
Rigg 5 *Will K. Cornwell 6 *Amy E. Zanne 1

\*co-senior authorship 1 Department of Biological Sciences, The George
Washington University, Washington, DC, United States 2 Hawkesbury
Institute for the Environment, University of Western Sydney, Hawkesbury,
Australia 3 Division of Natural Sciences, New College of Florida,
Sarasota, FL, United States 4 Department of Wood Science, University of
British Columbia, Vancouver, Canada 5 Elizabeth Macarthur Agricultural
Institute, NSW Department of Primary Industries, Meanagle, Australia 6
Evolution & Ecology Research Centre, School of Biological Earth and
Environmental Sciences, University of New South Wales, Sydney, Australia

## Data availability

Raw DNA sequencing data are available under NCBI BioProject ID
PRJNA690172

## Repository organization

Code to address three key analyses aims are in R markdown files. A brief
description of analyses and associated script are below.

1.  How do decay trajectories vary between wood species and sizes?
    **1\_decayPatterns.Rmd**

2.  How do wood traits covary? **2a\_woodtraits\_correlated.Rmd** Do
    wood traits explain decay? **2b\_woodtraits\_explainDecay.Rmd**

3.  Summarize endophyte communities across wood species and sizes
    **3a\_endoComp\_summary.Rmd** Does endophyte community composition
    explain decay? **3b\_endoComp\_explainDecay.Rmd** Does endophyte
    richness explain decay? **3c\_endoDiv\_explainDecay.Rmd**

-   How does including t=0 points affect the decay model fit? Including
    it affects the liklihood and the model selection criteria, but the
    curve fits are identical. Excluding t=0 leads to prefering simpler
    models, which is the same effect as increasing the penalty for model
    complexity **code/testing\_time\_zero.Rmd**

-   Do endophyte community dissimilarities correlate strongly with decay
    trajectory parameters?
    **code/initialDist\_vs\_decayDist\_btwCode.Rmd**
