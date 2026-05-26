### Quantile Regression with Constrained Splines

[![Licence: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.6+](https://img.shields.io/badge/python-3.6+-blue.svg)](https://www.python.org/downloads/)

Set of programmes for quantile regression with constrained splines, based on
the article  ‘Quantile regression with cubic splines under shape constraints’ by Alexandree Abbes,  https://doi.org/10.5281/zenodo.17427913 

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

# Example of using a file
python3 SplineCubicQuantBspkn3_V2_der3.py
# Example of global use
python3 Quantregtk.py


 Features

    ✅ Quantile regression with splines of degree 1 to 4

    ✅ Constraints on derivatives (valid over the entire interval or selected regions, not just at nodes)

    ✅ Graphical interface for interactive testing

    ✅ Standalone code with built-in tests


### Licence
MIT Licence - see the LICENCE file for more details
Citation

### Citation
# If you use this code in your research, please cite:
ABBES, A. (2025). Quantile regression with cubic polynomial splines under shape constraints with applications. https://doi.org/10.5281/zenodo.17427913

### Contributions
#Contributions are welcome! Feel free to:

    Open an issue

    Submit a pull request

    Suggest improvements


###  AI assistance

#This project was carried out with the assistance of DeepSeek (https://deepseek.com/), which helped the developer with:
- Implementing constrained spline algorithms from Matlab to Python
- Translating mathematical concepts into Python code
- Structuring programmes and documentation

DeepSeek is a language model developed by 深度求索 [1].

[1] DeepSeek-AI. (2024). DeepSeek-V3 Technical Report. arXiv:2412.19437.
