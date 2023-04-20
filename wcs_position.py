#convert the column data in the form of circle(x, y, 2)
#from astropysics.coords import ICRSCoordinates,GalacticCoordinates
#gcoords = ICRSCoordinates('2h34m12.32s',12.3).convert(GalacticCoordinates)


import numpy as np

#load the coordinates

data = np.loadtxt('data_xy.dat')

reg_file = open("region_xy.reg","w")
reg_file.write("# Region file format: DS9 version 4.1\n")
reg_file.write("global color=green dashlist=8 3 width=1 font=\"helvetica 10 normal roman\" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n")
reg_file.write("fk5\n")

# save in a new file in a new format
#np.savetxt(reg_file, data, fmt = 'circle(%.7f , %.7f, 1\")', newline = '\n')
