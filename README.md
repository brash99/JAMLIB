![jamlogo](gallery/jam.jpg)

## About
 
The repository contains a collection of codes/scripts (in fortran, python,
mathematica) along with interpolation tables for the collinear parton 
distributions in the nucleon as well as the collinear parton to hadron 
fragmentation functions. In addition we include: 

* polarized twist 3 up and down distributions.
* polarized structure functions g1 and g2 for proton, neutron, deuterium and helium.


## Quick start
The codes can be downloaded in two ways

* [Download the latest release](https://github.com/JeffersonLab/JAMLIB/archive/master.zip).
*  Clone the repo:  `$ git https://github.com/JeffersonLab/JAMLIB.git`.

To get the latest update, just pull from your local repo, e.g. `$ git pull`.

## Documentation
JAMLIB documentation is available in the [wiki](https://github.com/JeffersonLab/JAMLIB/wiki). 

## Theory
* The JAM analysis uses collinear factorization at NLO in perturbative QCD.
* For DIS data we include a treatment of higher twist as well as target mass corrections.

## Roadmap
Tha main gloal of the JAM analysis is to perform a universal fit to extract
collinear (un)polarized parton densities as well as parton to hadron
fragmentation functions from all existing available data. At present the
observables with checkmarks has been included in the JAM analysis.
<img src="gallery/roadmap.jpg" width="500">


## Iterative Monte Carlo method 
The JAM analysis uses a novel fitting procedure based on Monte Carlo techniques
to have robust estimates of expectation values and variances of the extracted 
distributions.  

<img src="gallery/workflow.jpg" width="400">





## Authors:
* Nobuo Sato  (nsato@jlab.org)
* Jake Ethier 
* Wally Melnitchouk  
* Alberto Accardi


