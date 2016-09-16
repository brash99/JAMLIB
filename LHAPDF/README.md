![jamlogo](../gallery/jam.jpg)

# LHAPDF grids 

## Overview

Grids can be downlaoded in two ways:

* Choose your grids and download the corresponding link in the first column (or click on "View raw" if you get a new page).
* Clone the JAMLIB repo, and get all the available grids:  `$ git https://github.com/JeffersonLab/JAMLIB.git`. With this option you can get the latest version by issuing `$ git pull`

Either way, you need to (unzip and) copy the grids to your local LHAPDF grid folder. Teh default location can be found, e.g., in the output of the `$ lhapdf` command. Further info on the LHAPDF interface installation and usge can be found on the [LHAPDF project web page](https://lhapdf.hepforge.org/).

A test ipython notebook is provided for convenience and to demonstrate usage of JAMLIB python code and LHAPDF grids. In the future, we'd like to add a python, C++, and fortran script, as well - help is welcome! - but for the moment you may have a look at the code xamples in the [LHAPDF project web page](https://lhapdf.hepforge.org/).

## Grids

**Note:** Grids with "preliminary" status or with version "v.-1" are in development and should not be downloaded.

| Name                                         | Info                                            | Reference                                                      | Status | Notes                                       |
| :--:                                         | :--:                                            | :--:                                                           | :--:   | :--:                                        |
| [JAM16FF_pi_Ceven](zip/JAM16FF_pi_Ceven.zip) | [.info](JAM16FF_pi_Ceven/JAM16FF_pi_Ceven.info) | [arXiv:1609.00899](http://inspirehep.net/record/1485196?ln=en) | v1     | Only glue and q^+ = q+qbar at NLO  |
| [JAM16FF_K_Ceven](zip/JAM16FF_K_Ceven.zip)   | [.info](JAM16FF_K_Ceven/JAM16FF_K_Ceven.info)   | [arXiv:1609.00899](http://inspirehep.net/record/1485196?ln=en) | v1     | Only glue and q^+ = q+qbar at NLO  |


## Questions/bugs

Please send us feedback or questions using the github 
[issue tracker system](https://github.com/JeffersonLab/JAMLIB/issues).


## Acknowledgments

Maintainer:
* Alberto Accardi

Many thanks to the following testers:
* Juan Guerrero (Hampton U. and Jefferson Lab)
* ... and all the other JAMLIB authors! 
