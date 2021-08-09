from generate_mesh import generate
from makeplot import makeplot
import numpy as np

pH_, V_, mesh = generate(conc = 1.5e-3, carbon = None, phosphate = 1e-3)
print(mesh)
makeplot(pH_, V_, mesh, 'plot', exp = True)