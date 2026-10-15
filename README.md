# Ecosphere_Submission_2026
This repository contains R scripts to conduct the Bayesian spatial capture recapture modeling described in "Density and habitat use of the southeastern beach mouse Peromyscus polionotus niveiventris" by Eric D. Stolen, Stephanie A. Legare, Timothy J. Kozusko, Stephanie Weiss, Dave Breininger, Tammy Foster, and Karen Holloway-Adkins. This manuscript is for consideration as a research article in Ecosphere. 

A brief description of the files:
1) "dAC_inhomog_hab.r": Nimble function for the probability of activity center placement as a function of habitat (needed to sample from the distribution according to an inhomogeneous point process).
2) "dMultiCapture.r": Nimble function for the multinomial likelihood calculation (used to speed up processing).
3) "Nimble model.R": R code to create the Nimble model.     
4) "R prep.r": R code to prepare the data.
5) "Run model in nimble.r": R code to set up and run the Nimble model in R.
6) "beach_mouse_SCR.RData": RData file with data needed to run the Bayesian analysis.

To run the model in R:
1. load(beach_mouse_SCR.RData)  # load the RData file: 
2. source("R prep.r) # run data prep script
3. source("dAC_inhomog_hab.r")  # load Nimble function
4. source("dMultiCapture.r")  # load Nimble function
5. source("Nimble model.R")  # load Nimble model code

