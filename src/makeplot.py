import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

colordict = {
    "Pb": "#85807b",
    "Pb++": "#db653d",
    "PbO": "#e04a16",
    "PbO2": "#3ec760",
    "Pb3O4": "#2a8744",
    "Pb++++": "#c92840",
    "HPbO2-": "#db7d42", 
    "PbO3--": "#bf3939",
    "Hc": "#191cf7",
    "PbCO3": "#007BA7", 
    "Pb3(PO4)2": "#d142f5"
    }

pHexp = [7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	10,	10,	10,	10,	10,	10,	10,	10,	10,	10,	10,	10,	10,	10,	10]
Vexp = [-0.62,	-0.62,	-0.16,	-0.12,	0.06,	0.48,	0.58,	0.67,	0.74,	1.08,	1.08,	-0.8,	-0.42,	-0.28,	-0.26,	-0.03,	0.02,	0.24,	0.34,	0.46,	0.75,	0.75,	0.84,	0.94,	0.94,	1.24]

def makeplot(pH_, V_, mesh, title, exp = False):
    
    print("making plot")
    
    #identify unique species, turn grid into numbers and make colormap
    nmesh = np.empty(np.shape(mesh), dtype=float)
    present = pd.unique(mesh.flatten())
    print("Species Present: ", present)
    colors_ = []
    for i in range(len(present)):
        nmesh[mesh == present[i]] = i
        colors_.append(colordict[present[i]])
    print(nmesh)
    print(colors_)
    levels = np.arange(len(present)+1)
    levels = levels - .5
    print(levels)
    
    fig, ax = plt.subplots()
    
    #fill regions
    CS = ax.contourf(pH_, V_, nmesh, levels, colors = colors_)
    ax.contour(pH_, V_, nmesh, colors= 'k', linewidths=0.25, antialiased=True)
    
    #experimental data
    if exp: 
        for i in list(range(3)) + list(range(11,12)):
            plt.plot(pHexp[i], Vexp[i], 'o', color = 'k')
        for i in list(range(3,11))+list(range(13,26)):
            plt.plot(pHexp[i], Vexp[i], '^', color = 'red', markeredgewidth=0.5, markeredgecolor='k')
    
    #water stability lines
    a = 1.299-0.0592*pH_[0]
    b = -0.059*pH_[0]
    plt.plot(pH_[0], a, '--', color='b')
    plt.plot(pH_[0], b, '--', color='b')
    plt.text(0,1.1,'Water Oxidation', fontsize=11, rotation=-7, color='b')
    plt.text(0,-0.2,'Water Reduction', fontsize=11, rotation=-7, color='b')
    
    #axes    
    plt.xticks(np.arange(-2,17,2))
    plt.yticks(np.arange(-2,4.5,1))
    plt.xlabel('pH')
    plt.ylabel('$\mathrm{V_{SHE}}$')
    
    plt.savefig(title, dpi = 300)
    
    