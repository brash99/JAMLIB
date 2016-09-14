![jamlogo](../gallery/jam.jpg)

# Fortran Interface

In this directory you will find two fortran scripts:

   - JAMLIB.f: contains FUNCTION get_xF which returns
     	       x times the PDF or FF distribution as well
	       as the interpolation routine
   - example.f: contains an example on how to use JAMLIB.f

## Comments

- The tables are located inside JAM libraries (JAM15,JAM16,etc.)
- To run example:
  - $gfortran -o example example.f JAMLIB.f
  - $./example
- NOTE: The available fortran code was only tested with the
  	gfortran compiler!

## Questions/Comments

- Please contact [Jacob Ethier](mailto:jethier@jlab.org) with
any questions or comments regarding the fortran codes/tables
 
## Authors
* Nobuo Sato (Jefferson Lab)
* Alberto Accardi (Hampton University and Jefferson Lab)
* Jacob Ethier (College of William and Mary)
* Wally Melnitchouk (Jefferson Lab)

