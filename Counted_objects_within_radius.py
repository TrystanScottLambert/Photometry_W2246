import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u

pixel_scale=0.16*u.arcsec
#rings = [
#           (1566, 1042, 250)
#]
ring_center = [1566, 1042]

PATH= '/Users/dejenewoldeyes/Documents/SExtractor/Main_catalog_data/The_final_plots/'

position=np.loadtxt(PATH + 'LBG_pos.dat')
block_LBG=np.zeros(position.shape[0], dtype=np.bool)


def count_LBG(pos, ring_center, r_in, r_out):
  dx = pos[:,0]-ring_center[0]
  dy = pos[:,1]-ring_center[1]
  r = (dx**2+dy**2)**0.5
  return len(r[(r>r_in) & (r<=r_out)])
  #return np.where(r<ring[2], True, block_LBG)
  
r_ins = np.arange(10, 1260, 125)
for r_in in r_ins:
  print(r_in, r_in+125, count_LBG(position, ring_center, r_in, r_in+125))
  print(r_in*pixel_scale)
  