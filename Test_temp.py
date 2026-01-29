import numpy as np
import cvxpy as cp
import cvxpy_gurobi
from scipy.interpolate import BSpline, splev, PPoly
from scipy.stats import percentileofscore
from matplotlib import pyplot as plt
import warnings
import mosek

import pandas as ps

from SplineCubicQuantBspkn3_V2_OK import *

temp=ps.read_excel('temp.xls')
temp_val=temp.values
year=temp_val[:,0]
ytab=temp_val[:,1]



xtab=(year-year[0])/(year[-1]-year[0])


year_kn=np.array([ 1880, 1889, 1900, 1910, 1930,1940, 1965, 1992])
knots=(year_kn-1880)/(1992-1880)
kn=7


cv=0#[1,1,-1,-1]

list_result=[]

fig, axes = plt.subplots(1,3, figsize=(15, 4))

axes[0].set_xlabel('year')
axes[0].set_ylabel('Diff. Temp.')
axes[0].set_title('temp with full constraints')

axes[0].grid(True, alpha=0.3)

axes[1].set_xlabel('year')
axes[1].set_ylabel('Diff. Temp.')
axes[1].set_title('temp without constraints')

axes[1].grid(True, alpha=0.3)

axes[2].set_xlabel('year')
axes[2].set_ylabel('Diff. Temp.')
axes[2].set_title('temp with partial constraints')

axes[2].grid(True, alpha=0.3)


x_eval = np.linspace(0, 1, 200)
def yr(x):
    return 1880+x*(1992-1880) 

for tau in [0.01,0.5,0.99]:
     
      spline_result= SplineCubicQuantBspkn3(xtab, ytab, knots, tau, 1, cv)
      y_eval = spline_result(x_eval)
      axes[0].plot(yr(x_eval), y_eval, linewidth=1, label='tau'+str(tau))
      spline_result= SplineCubicQuantBspkn3(xtab, ytab, knots, tau, 0, cv)
      y_eval = spline_result(x_eval)
      axes[1].plot(yr(x_eval), y_eval, linewidth=1, label='tau'+str(tau))
      monot=[0,0,0,0,0,-1,0]
      
      spline_result= SplineCubicQuantBspkn3(xtab, ytab, knots, tau, monot, cv)
      y_eval = spline_result(x_eval)
      axes[2].plot(yr(x_eval), y_eval, linewidth=1, label='tau'+str(tau))

axes[0].scatter(yr(xtab), ytab, alpha=0.5, s=20, label='data')
axes[0].legend()
axes[1].scatter(yr(xtab), ytab, alpha=0.5, s=20, label='data')
axes[1].legend()
axes[2].scatter(yr(xtab), ytab, alpha=0.5, s=20, label='data')
axes[2].legend()
plt.show()
