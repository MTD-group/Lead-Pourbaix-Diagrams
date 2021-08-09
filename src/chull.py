#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 22 11:45:30 2020

@author: alex
"""

import matplotlib.pyplot as plt

x = [0, 0.33333, 0.43, 0.5, 1]

PBESolnospin = [0,-1.198046802365,-1.28165956625446,-1.2980779290629,0]

PBESolspin = [0,-0.781454585698333,-0.902264773397317,-0.9488537040629,0]

PBE = [0,-0.689102119031667,-0.861826944825888,-0.9212704790629,0]

LDA = [0,-1.059241252365,-1.12240284482589,-1.1483685790629,0]

Exp = [0, -0.75106433333,  -0.890143, -0.9737235, 0]

MP = [0, -1.31, -1.456, -1.477, 0]

SCAN = [0,-0.773937335698331,-0.915539451968747,-0.978315204062902,0]

HSE = [0, -0.704489645698334, -0.855641074825888, -0.9231525615629, 0]


# print('difference btwn PBESolspin and Exp')

# print("PbO2 ", (PBESolspin[1]-Exp[1]), "eV", (PBESolspin[1]-Exp[1])/Exp[1]*100, "%")
# print("Pb3O4", (PBESolspin[2]-Exp[2]), "eV", (PBESolspin[2]-Exp[2])/Exp[2]*100, "%")
# print("PbO  ", (PBESolspin[3]-Exp[3]), " eV", (PBESolspin[3]-Exp[3])/Exp[3]*100, "%")



fig, ax = plt.subplots()
ax.plot(x, SCAN, label='SCAN, SOC', marker = "D", color='darkorange')
ax.plot(x, PBE, label = 'PBE, SOC', marker = "o")
ax.plot(x, HSE, label = 'HSE, no SOC', marker = "^", color='purple')
ax.plot(x, PBESolspin, label = 'PBEsol, SOC', marker = ">", color='limegreen')

for i in range(5):
    ax.plot(x[i], Exp[i], 'x', color='k', markersize=10)

ax.plot(x, LDA, label = 'LDA, SOC', marker = "<", color='r')
ax.plot(x,PBESolnospin, label = 'PBEsol, no SOC', marker = "s",linestyle='dashed', color='darkgreen')

#ax.plot(x, MP, label = "Materials Project", marker = ".")

ax.legend()
ax.set(xlim=(0, 1), ylim=(-1.5, 0))
plt.xlabel('Pb$_{1-x}$O$_{x}$ Composition', fontsize=13)
plt.ylabel('eV atom$^{-1}$', fontsize=13)
plt.savefig('Energies', dpi=300)