""""This activity is intended to correct intergalactic medium (IGM) absorption for our selection criteria.
 We are using the Lyman break galaxy (LBG) composite spectrum from Shaply et al. (2003); 
 in their study, they used the IGM absorption for redshift 3, 
 to apply the Madau (1995) Lyman alpha forest at targeted redshift, 
 we removed the IGM absorption from the LBG composite Spectrum at z~3 and we applied 
the IGM absorption at the Hot DOG redshift, by taking into account the IGM absorption from Madau (1995)."""


import numpy as np
import matplotlib.pyplot as plt
from synphot import SourceSpectrum, Observation, units
from synphot.spectrum import SpectralElement
from synphot.models import Empirical1D
from synphot import SourceSpectrum, etau_madau

#Load the spectrum.
#data = np.loadtxt("LBG_Composite_Spectrum_Without_IGM.dat")
data= np.loadtxt("Composite_LBG.dat")
LBG = SourceSpectrum(Empirical1D, points=data[:,0], lookup_table=data[:,1]*units.FNU, keep_neg=True)





imacs_g = np.loadtxt("LCO_IMACS.sloan_g.dat")
imacs_r = np.loadtxt("LCO_IMACS.sloan_r.dat")
imacs_i = np.loadtxt("LCO_IMACS.sloan_i.dat")
imacs_z = np.loadtxt("LCO_IMACS.sloan_z.dat")



z=np.arange(3.10, 4.0, 0.1)
z =4.0

LBG.z =z
wave_shi = data[:,0] * (1 + z)
print(wave_shi)
extcurve = etau_madau(wave_shi, z)
tau_eff = -np.log(extcurve(data[:,0]*(1+z)))
kappa = 1
tau_use = kappa*tau_eff
lbg_fnu_with_igm = data[:,1] * np.exp(-tau_use) * units.FNU
LBG_with_IGM = SourceSpectrum(Empirical1D, points=wave_shi, lookup_table=lbg_fnu_with_igm, keep_neg=True)


bp_imacs_g = SpectralElement(Empirical1D, points=imacs_g[:,0], lookup_table=imacs_g[:,1])
bp_imacs_r = SpectralElement(Empirical1D, points=imacs_r[:,0], lookup_table=imacs_r[:,1])
bp_imacs_i = SpectralElement(Empirical1D, points=imacs_i[:,0], lookup_table=imacs_i[:,1])
bp_imacs_z = SpectralElement(Empirical1D, points=imacs_z[:,0], lookup_table=imacs_z[:,1])

obs_g = Observation(LBG_with_IGM, bp_imacs_g, force = 'extrap')#, binset='array-like')
obs_r = Observation(LBG_with_IGM, bp_imacs_r, force = 'extrap')#, binset='array-like')
obs_i = Observation(LBG_with_IGM, bp_imacs_i, force = 'extrap')#,  binset='array-like')
obs_z = Observation(LBG_with_IGM, bp_imacs_z, force = 'extrap')#,  binset='array-like')




#plots

plt.plot(LBG.waveset, LBG(LBG.waveset, flux_unit='fnu')*10**(29)/(5.6),label=r' LBG Composite at z=4.601',color='xkcd:grey',  linestyle='--')

plt.plot(imacs_g[:,0], imacs_g[:,1],color='xkcd:blue',label=r'IMACS-g, -r, -i, z-band$')
plt.plot(imacs_r[:,0], imacs_r[:,1], color='xkcd:red')
plt.plot(imacs_i[:,0], imacs_i[:,1], color='xkcd:brown')
#plt.plot(imacs_z[:,0], imacs_z[:,1], color='xkcd:magenta')

#plt.xlim(3800, 11300)
plt.ylim(0.001, 1.11015)
plt.xlabel(r'Obseverd wavelength $(\AA)$ ', fontsize=16)
plt.ylabel(r'$f\nu$ ($\mu Jy$)', fontsize=16)

print('redshift', z)
print('g-band = ', obs_g.effstim('abmag'))
print('r-band = ',obs_r.effstim('abmag'))
print('i-band = ',obs_i.effstim('abmag'))
print('z-band = ',obs_z.effstim('abmag'))
plt.savefig('LBG_4.0.png')
plt.show()