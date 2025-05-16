# Lag Profile Inversion (LPI)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6405720.svg)](https://doi.org/10.5281/zenodo.6405720)

LPI is an R package for deconvolving lag profiles from voltage level incoherent scatter radar data. This version is the old version with socket cluster. See the master branch for an MPI version that works in HPC environment.

Data I/O routines need to be provided by other packages. See the packages [LPI.gdf](https://github.com/ilkkavir/LPI.gdf) and [LPI.KAIRA](https://github.com/ilkkavir/LPI.KAIRA) for examples. 


Ilkka Virtanen (ilkka.i.virtanen@oulu.fi) 2025

## INSTALLATION

1. Clone the github repository with
```
git clone https://github.com/ilkkavir/LPI
```

2. Install from *nix command line
```
R CMD INSTALL <path-to-LPI-package>
```

3. OR install from R console
```
install.packages(pkgs=<path-to-LPI-package>,repos=NULL)
```

## User instructions
See the attached manuals, LPI-manual.pdf and LPI-tutorial.pdf, for instructions.

