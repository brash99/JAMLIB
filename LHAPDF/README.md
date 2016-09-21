[![jamlogo](../gallery/jam.jpg)](http://www.jlab.org/jam)

# LHAPDF grids 

## Overview

Grids can be downlaoded in two ways:

* Choose your grids and download the corresponding link in the first column (or click on "View raw" if you get a new page).
* Clone the JAMLIB repo, and get all the available grids:  `$ git https://github.com/JeffersonLab/JAMLIB.git`. With this option you can get the latest version by issuing `$ git pull`

Either way, you need to (unzip and) copy the grids to your local LHAPDF grid folder. The default location can be found, e.g., in the output of the `$ lhapdf` command. Further info on the LHAPDF interface installation and usage can be found on the [LHAPDF project web page](https://lhapdf.hepforge.org/).

A test ipython notebook is provided for convenience and to demonstrate usage of JAMLIB python code and LHAPDF grids. In the future python, C++ and fortran scripts will be added; for the moment you may have a look at the  [code examples](https://lhapdf.hepforge.org/codeexamples.html) in the LHAPDF project web page.

## Grids

**Note:** Grids with "preliminary" status or with version "v.-1" are in development and should be used with caution.

| Name                                         | Info                                            | Reference                                                      | Status | Notes                                       |
| :--:                                         | :--:                                            | :--:                                                           | :--:   | :--:                                        |
| [JAM15_PPDF_Ceven](zip/JAM15_PPDF_Ceven.zip) | [.info](GRIDS/JAM15_PPDF_Ceven/JAM15_PPDF_Ceven.info) | [Phys.Rev. D93 (2016) 074005](http://inspirehep.net/record/1418180?ln=en) | v.-1  PRELIMINARY   | q+ (=q+qbar) & glue polarized PDFs at NLO     |
| [JAM15_T3PPDF_Ceven](zip/JAM15_PPDF_Ceven.zip) | [.info](GRIDS/JAM15_T3PPDF_Ceven/JAM15_T3PPDF_Ceven.info) | [Phys.Rev. D93 (2016) 074005](http://inspirehep.net/record/1418180?ln=en) | v.-1  PRELIMINARY   | q+ (=q+qbar) & glue polarized twist-3 PDFs  at NLO     |
| [JAM16FF_pi_Ceven](zip/JAM16FF_pi_Ceven.zip) | [.info](GRIDS/JAM16FF_pi_Ceven/JAM16FF_pi_Ceven.info) | [arXiv:1609.00899](http://inspirehep.net/record/1485196) | v.-1  PRELIMINARY   | q+ (=q+qbar) & glue to pions at NLO     |
| [JAM16FF_K_Ceven](zip/JAM16FF_K_Ceven.zip)   | [.info](GRIDS/JAM16FF_K_Ceven/JAM16FF_K_Ceven.info)   | [arXiv:1609.00899](http://inspirehep.net/record/1485196) | v.-1  PRELIMINARY   | q+ (=q+qbar) and glue to kaons at NLO             |


## Questions/bugs

Send feedback or questions using the github 
[issue tracker system](https://github.com/JeffersonLab/JAMLIB/issues).


## Acknowledgments

Maintainer:
* Alberto Accardi

Many thanks to the following testers:
* Juan Guerrero (Hampton U. and Jefferson Lab)
* ... and all the other JAMLIB authors! 
