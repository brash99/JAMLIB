[![jamlogo](gallery/jam.jpg)](http://www.jlab.org/jam)

## About
 
The repository contains interpolation tables for the collinear parton
distributions in the nucleon, as well as the collinear parton to hadron
fragmentation functions. The interpolation codes and grids are available in:

* [python](https://github.com/JeffersonLab/JAMLIB/tree/master/python)

* [LHAPDF format](https://github.com/JeffersonLab/JAMLIB/tree/master/LHAPDF) (grids only)

* [fortran](https://github.com/JeffersonLab/JAMLIB/tree/master/fortran)

The usage is described in the README files within each subfolder. 

## Available sets
| Table         | Reference         | Notes                                       |
| :--           | :--:              | :--                                         |
| JAM15 PPDF    | [inspirehep][jam15] | q^+ (=q+qbar) and glue at NLO                |
| JAM15 T3PPDF  | [inspirehep][jam15] | u and d type twist-3 distributions          |
| JAM15 T4g1    | [inspirehep][jam15] | twist-4 part of g1 for proton and neutron   |
| JAM16 FFpion  | [inspirehep][jam16] | q^+ and glue at NLO                |
| JAM16 FFkaon  | [inspirehep][jam16] | q^+ and glue at NLO                |
[jam15]:https://inspirehep.net/record/1418180
[jam16]:http://inspirehep.net/record/1485196?ln=en

## Quick start
The JAMLIB can be downloaded in two ways

* [Download the latest release](https://github.com/JeffersonLab/JAMLIB/archive/master.zip).
*  Clone the repo:  `$ git https://github.com/JeffersonLab/JAMLIB.git`.

To get the latest update, pull from your local repo, e.g. `$ git pull`.




## Questions/bugs
Please send us feedback or questions using the github 
[issue tracker system](https://github.com/JeffersonLab/JAMLIB/issues).


## Authors
* Nobuo Sato (Jefferson Lab)
* Alberto Accardi (Hampton U. and Jefferson Lab)
* Jacob Ethier (College of William and Mary)
* Wally Melnitchouk (Jefferson Lab)

