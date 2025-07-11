from uetools import *
from uetools.UeVacuum import VacuumRegion
from UeVacuum.VacuumRegion import Surface

import matplotlib.pyplot as plt

c = Case()
c.vacuum.tokamakPlot()


'''Triangle Geometry'''
#c.vacuum.trianglePlot()

'''Square Geometry'''
# c.vacuum.squarePlot()

'''Shaded Geometry'''
# startingSurface = Surface((3, 1), (1, 1), 0)
# c.vacuum.geometries(startingSurface.lineOfSightVertices(), 1)

'''Full Tokamak Geometry'''
# c.vacuum.tokamakPlot()

input("Press 'Enter' to close plots.")
plt.close('all')