# -*- coding: utf-8 -*-
"""
Created on Sat Nov 18 00:00:58 2017

@author: AlexMialdun
"""

# %reset -f
import numpy as np
import numpy.polynomial.polynomial as poly
import matplotlib.pyplot as plt
import libPolyRK as rk


#_load original data_:
fname_data = 'rho_TEG-Wat_Guo2012.txt'
data = np.genfromtxt(fname_data, delimiter='\t', 
                     dtype=['f8', 'f8', 'f8'], names=True)
#_extraction of the data columns from the 'data' object :
x_1 = data['x_TEG'] # molar fraction of triethylene glycol
# rho_1 = data['rho_g_cm3'] # density of the mixture
VE_1 = data['V_E_cm3_mol'] # excess molar volume of the mixture

#_plot original data_:
fig, ax1 = plt.subplots() # figsize=(11, 9)
ax1.set_xlabel('$x$ (TEG) / mole frac.', fontsize=14)
ax1.set_ylabel(r'$V^E$ / cm$^{3}$ mol$^{-1}$', fontsize=14)
line1, = ax1.plot(x_1, VE_1, 'ro', label='Guo 2012')
ax1.legend(handles=[line1], loc=0, fontsize=13)
ax1.tick_params(axis='both', which='major', labelsize=12)
fig.tight_layout()
# plt.savefig('VE-x_TEG-Wat.png')
plt.show()

#_fit to Redlich-Kister polynomial_:
A_VE = rk.polyfitRK(x_1, VE_1, 3)
print(A_VE)
# x_extr = rk.findExtremaRK(A_VE)

#_convert to power series polynomial_:
P_VE = rk.polyConvRK2power(A_VE)
print(P_VE)

#_estimate deviation_:
VE_dev = rk.polyvalRK(A_VE, x_1)
delta = VE_1 - VE_dev

#_generate fitted curve_:
xi = np.linspace(0.0, 1.0, 101)
VE_calc_RK = rk.polyvalRK(A_VE, xi)
VE_calc_Pow = poly.polyval(xi, P_VE)


#_plot data with fit curve and deviations_:
fig = plt.figure(figsize=(6, 6)) # figsize=(11, 9)
ax1 = fig.add_subplot(3,1,1)
ax1.set_ylabel('$\delta$ / cm$^{3}$ mol$^{-1}$', fontsize=14)
ax1.plot(x_1, delta, 'bo-', ms=3)
Xlim = [-0.02, 1.02]
ax1.set_xlim(Xlim)
# ax1.set_ylim(-0.026, 0.026)
ax1.tick_params(axis='both', which='major', labelsize=12)
# 
ax2 = fig.add_subplot(3,1,(2,3))
ax2.set_xlabel('$x$ (TEG) / mole frac.', fontsize=14)
ax2.set_ylabel('$V^E$ / cm$^{3}$ mol$^{-1}$', fontsize=14)
line1, = ax2.plot(x_1, VE_1, 'bo', ms=5, label='data points')
line2, = ax2.plot(xi, VE_calc_RK, 'g--', label='fit by R-K polynomial')
line3, = ax2.plot(xi, VE_calc_Pow, 'r:', label='power series polynomial')
ax2.legend(handles=[line1, line2, line3], loc=0, fontsize=13)
ax2.set_xlim(Xlim)
ax2.set_ylim(-0.89, 0.04)
ax2.tick_params(axis='both', which='major', labelsize=12)
fig.tight_layout()
# plt.savefig('VE-x_TEG-Wat_fit.png')
plt.show()
