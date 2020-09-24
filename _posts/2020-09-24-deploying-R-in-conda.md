---
layout: post
title: "Deploying R in conda"
author: "RRBIII"
categories: misc
tags: [misc,R,conda]
image: condaR-1.png
---


## a solution to the problem of portability

Working with R, you may have come across the frustration of figuring out where are put stuff, or the more frustrating 'trying to install R as user'. Conda has spent a fair but of time setting it up so that conda virtual environments can store all of the necessary files, and access local libraries very efficiently.  

```sh
conda create -n r4.0 -y python=3 r-base r-essentials r-devtools
```

#### So now we are done right?

Not exactly, as a plus, you can specify your R version with `r-base==3.6.3` and if you are wondering about bioconductor, then simply `bioconductor-package` will install any package in the right place. However, you will eventually run into the problem that there is an R package that isn't available in conda. In this case, you can do as below, where I build a custom r-tigris and r-gg.gap package for a shiny app with maps. You can actually push these to your private conda channel in the cloud if you register for one, but more often, you will build it where you want it, and specify `--use-local` in your install to have conda check the locally built things. Also, you can clear out the temo build files with `conda build purge`

```sh
#!/usr/bin/env bash
exec 1>> dashboard_build.log 2>&1
set -ex

# construct tigris
conda skeleton cran tigris
conda-build r-tigris

# construct gg.gap
conda skeleton cran gg.gap
conda-build r-gg.gap

# cleanup after builds
conda-build purge
rm -rf r-tigris/
rm -rf r-gg.gap/

# build initial
conda create -n shiny -y --use-local python=3 r-tigris r-gg.gap r-shiny \
  r-readr r-tidyr r-dplyr r-ggplot2 r-data.table r-shinyalert \
  r-rcolorbrewer r-gridextra r-stringr r-shinydashboard r-ggpubr r-ggplotify \
  r-leaflet r-units r-sf r-rgdal r-leafsync r-reshape2

# # Can be launched from Rscript (run in background screen)
# screen -S covid
# source $(conda info --base)/etc/profile.d/conda.sh
# conda activate shiny || { echo "shiny Conda environment not activated"; exit; }
# Rscript covid_19_data_viz_v6.R
```

#### Next time
One of the things this doesn't address is installing github packages. If they are simple, then hopefully after you are done installing R (include `r-devtools`), you can run R and `install_github('', dependencies=FALSE)` from inside. However, beware, because if it has dependencies, R will try to install them where it knows to: `/usr/local` and this will create some havoc. Instead, you may have to research those dependencies and install them prior with conda, as I did here with r-units, r-sf, r-rgdal. Then, go back into R and try again. Repeat until done. Ugh, never an easy answer is there.  


