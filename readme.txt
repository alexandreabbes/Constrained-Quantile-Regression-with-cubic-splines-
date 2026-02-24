### Quantile Regression with Constrained Splines

[![Licence: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.6+](https://img.shields.io/badge/python-3.6+-blue.svg)](https://www.python.org/downloads/)

Set of programmes for quantile regression with constrained splines, based on
the article  ‘Quantile regression with cubic splines under shape constraints’ by Alexandree Abbes,  https://zenodo.org/records/16999785. 

### Prerequisites (the cvxpy library is more difficult to install)

```bash
pip install pandas numpy scipy warnings matplotlib cvxpy

###Installation
#Copy all files to the same directory. 


### File structure
#Graphical interface
    Quant_reg_tk.py - Tkinter graphical interface for interactive testing of parameters (spline degree, nodes, derivative constraints)

#Main algorithms

#Each file is self-contained and includes tests:
#File        Description
SplineCubicQuantBspkn3__V2_der3.py    Cubic splines with constraints on 1st, 2nd, and 3rd derivatives
SplineQuarticQuantBspkn4_V3.py    Quartic splines with exact constraints over the entire interval
SplineQuadQuantBspkn2_V1.py	Quadratic splines with constraints on first and second derivatives
SplineLinQuantBspkn1_V1.py    Linear splines with constraints on first derivative

#Examples and data
    Test_temp.py - Replication of the test on global temperatures (data in temp.xls)



## Quick start
# Example

Translated with DeepL.com (free version)
## Quick use
# Example of using a file
python3 SplineCubicQuantBspkn3__V2_der3.py
# Example of global use
python3 Quantregtk.py


 Features

    ✅ Quantile regression with splines of degree 1 to 4

    ✅ Constraints on derivatives (theoretically implemented over the entire interval)

    ✅ Graphical interface for interactive testing

    ✅ Standalone code with built-in tests


### Licence
MIT Licence - see the LICENCE file for more details

### Contributions
#Contributions are welcome! Feel free to:

    Open an issue

    Submit a pull request

    Suggest improvements

Translated with DeepL.com (free version)
