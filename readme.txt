### Quantile Regression with Constrained Splines 
### R code branch 
### 04/2026
[![Licence: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.6+](https://img.shields.io/badge/python-3.6+-blue.svg)](https://www.python.org/downloads/)

Set of programmes for quantile regression with constrained splines, based on
the article  ‘Quantile regression with cubic splines under shape constraints’ by Alexandree Abbes,  https://zenodo.org/records/16999785. 

### Prerequisites : install first the cvxR library

library(CVXR)
library(pracma)
library(spline) #only for the comarison test test_bspline() in the file PP_spline2.R
⚠️ Known issues: Several functions in the  libraries pracma have inconsistencies. Refer to Bug_polynomial.txt for documentation. The necessary corrections and reimplementations are included in this package.
The bs() function of the R library spline , calculating the B-spline basis values at given x, gives the same results our function bs_direct(). see the test_bspline() function.

⚠️ THE CONVENTION FOR POLYNOMIALS IS NOT THE STANDARD ONE:
#FOR ALL THE FUNCTIONS HERE p0+p1x+p2x^2  <-> c(p0,P1,p2)
#THUS ALL POLYNOMIAL OPERATIONS FOLLOW THIS CONVENTION AND HAVE BEEN RE-WRITTEN



### First use:
Copy the files PPèspline2.R and the main file SplineCubicQuantBspkn3_Kar_V10.2.R
in a directory.
First execute the code PP-spline2.R (for the calculation of the bspline basis coefficients)
then run the main code to see the little demo.
You shoul see the same as : Contained_quantile_reg_Bspline_test_R.png

Feel free to modify the constrains, the function, the noisy part, and the tau. 


### Licence
MIT Licence - see the LICENCE file for more details
Citation

### Citation
# If you use this code in your research, please cite:
ABBES, A. (2025). Quantile regression with cubic polynomial splines under shape constraints with applications. https://doi.org/10.5281/zenodo.16999785

### Contributions
#Contributions are welcome! Feel free to:

    Open an issue

    Submit a pull request

    Suggest improvements
