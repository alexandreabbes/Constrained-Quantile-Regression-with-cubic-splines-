Content 
File                                Description
---------------------------------------------------------------------------------------------------------------
SplineCubicQuantBspkn3_V2_OK.py		===========	python version of the quantile regression with cubic splines under constraints (monotone and convex). A simple test is provided.

SplineCubicQuantBspkn3__V2_der3.py ==========		Same as SplineCubicQuantBspkn3_V2_OK.py ,also implementing constraints on the	3rd  derivative

Test_temp.py	============				test on temperature identical to the test in the article, but implemented with python

temp.xls	===========				data set of global temperatures

SplineQuarticQuantBspkn4_V2_OK.py	============	Quantile Regression under constraints with quartic splines.

These codes use some ordinary python libraries, (numpy and scipy).
It also needs **cvxpy**, that has to be specificaly downloaded.
