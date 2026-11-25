# Sex differences in human music perception are negligible 

> [!WARNING]  
> This repository is currently under construction 🚧

This repository contains all data, code, and text for the manuscript **Sex differences in human musicality are negligible**. The preprint for this manuscript is publicly available at https://www.biorxiv.org/content/10.1101/2023.05.23.541970v3.

For assistance, please contact the corresponding authors: Mila Bertolo (mila.bertolo@mail.mcgill.ca) and Samuel Mehr (mehr@hey.com).

## Anatomy of the repo
To render the paper, `knit` the .Rmd file `./writing/manuscript-rev1.Rmd`. Running this file will import the pre-processed data stored in `./data` in order to run all analyses and build all visualizations reported in the manuscript. This will produce the manuscript `./writing/manuscript-rev1.pdf`.

> [!IMPORTANT]
> This repository requires `git-lfs` as it contains some large data files. By default, `git-lfs` is tracking `*.rds` files.

#### System Requirements 

This code was last tested in R version 4.6.0 (2026-04-24). 

#### Installation Guide 

The `load_libraries_and_functions` chunk in `manuscript-rev1.Rmd` handles package installation. This chunk executes as part of knitting the file. Depending on your pre-existing R / RStudio environment, this installation may take a couple minutes. Knitting `manuscript-rev1.Rmd` to its output .pdf may take <5 minutes.
