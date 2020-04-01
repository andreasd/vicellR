[![DOI](https://zenodo.org/badge/252146142.svg)](https://zenodo.org/badge/latestdoi/252146142)

# Introduction

ViCellR is a library for R, which reads ViCell-txt-files into R, and can generate batch figures, calculate growth rates and specific productivities.

# Materials
* R (v3.2 or higher)
* Optional: An easy to use program for R(eg. RStudio)
* Dependencies: Execute install.packages(c("lubridate", "qdap", "xlsx", "ggplot2")) to install dependencies
* Txt-Files produced by ViCell-measurements used following file naming scheme

# Naming scheme
Files have to be named following a specific format, so that the script can extract relevant information like samples belonging to the same timepoint or replicates. Use this format to name your samples when entering them in the ViCell. Alternatively, the names can be changed later during the script execution, but this requires additional code/manual editing.

[YYYYMMDD]-[INITIALS]-[SAMPLENAME]-TP[##]-R[##]

Text in [] should be replaced with your own values (without the []).

TP = Timepoint of this measurement, to group together measurements that were taken at the same time point. TP followed by a number. Time points should be named with increment

R = Replicate, to group together replicates of the same measurement. R followed by a number. Those can be combined for plotting and growth calculations. For this to work, sample names in the file name have to be identical

Example: 20170530-AD-CHO_K1-TP03-R03

# Usage
For examples, please see testRuns added to the library, as they should be always up to date.

# Authors
Andreas B. Diendorfer and Gerald Klanert
