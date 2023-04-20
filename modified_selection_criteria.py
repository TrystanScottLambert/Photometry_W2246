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


PATH = '/Users/dejenewoldeyes/Documents/IMACS_LBG_Hot_DOGs/IMACS_proposal/my_project/'
#PATH = 'data/'
#Load the spectrum.
#data = np.loadtxt("LBG_Composite_Spectrum_Without_IGM.dat")
data= np.loadtxt(PATH +"Composite_LBG.dat")
LBG = SourceSpectrum(Empirical1D, points=data[:,0], lookup_table=data[:,1]*units.FNU, keep_neg=True)


data = np.loadtxt("LBG_Composite_Spectrum_Without_IGM.dat")
LBG = SourceSpectrum(Empirical1D, points=data[:,0], lookup_table=data[:,1]*units.FNU, keep_neg=True)
LBG.plot(flux_unit='fnu')

Subaru_r = np.loadtxt("Subaru_Suprime_Rc.dat")
Subaru_i = np.loadtxt("Subaru_Suprime_i.dat")
Subaru_z = np.loadtxt("Subaru_Suprime_z.dat")


Gemini_GMOS_r = np.loadtxt("Gemini_GMOS_r_filter.dat")
Gemini_GMOS_i = np.loadtxt("Gemini_GMOS_i_filter.dat")
Gemini_GMOS_z = np.loadtxt("Gemini_GMOS_z_filter.dat")

#z=np.arange(4.50, 5.0, 0.1)
z =4.601

LBG.z =z
wave_shi = data[:,0] * (1 + z)
print(wave_shi)
extcurve = etau_madau(wave_shi, z)
tau_eff = -np.log(extcurve(data[:,0]*(1+z)))
kappa = 1
tau_use = kappa*tau_eff
lbg_fnu_with_igm = data[:,1] * np.exp(-tau_use) * units.FNU
LBG_with_IGM = SourceSpectrum(Empirical1D, points=wave_shi, lookup_table=lbg_fnu_with_igm, keep_neg=True)
#np.savetxt('LBG_Composite_Spectrum_With_IGM.dat', lbg_fnu_with_igm, fmt='%1.4e') 
bp_subaru_r = SpectralElement(Empirical1D, points=Subaru_r[:,0], lookup_table=Subaru_r[:,1])
bp_subaru_i = SpectralElement(Empirical1D, points=Subaru_i[:,0], lookup_table=Subaru_i[:,1])
bp_subaru_z = SpectralElement(Empirical1D, points=Subaru_z[:,0], lookup_table=Subaru_z[:,1])

obs_sr = Observation(LBG_with_IGM, bp_subaru_r, force = 'extrap')#, binset=binset)
obs_si = Observation(LBG_with_IGM, bp_subaru_i, force = 'extrap')#, binset=binset)
obs_sz = Observation(LBG_with_IGM, bp_subaru_z, force = 'extrap')#, binset=binset)

ris=obs_sr.effstim('abmag')- obs_si.effstim('abmag')
izs=obs_si.effstim('abmag')- obs_sz.effstim('abmag')
print(z)
print(izs)
print(ris)
print(obs_sr.effstim('abmag'))
print(obs_si.effstim('abmag'))
print(obs_sz.effstim('abmag'))

bp_gmos_r = SpectralElement(Empirical1D, points=Gemini_GMOS_r[:,0], lookup_table=Gemini_GMOS_r[:,1])
bp_gmos_i = SpectralElement(Empirical1D, points=Gemini_GMOS_i[:,0], lookup_table=Gemini_GMOS_i[:,1])
bp_gmos_z = SpectralElement(Empirical1D, points=Gemini_GMOS_z[:,0], lookup_table=Gemini_GMOS_z[:,1])

obs_gr = Observation(LBG_with_IGM, bp_gmos_r, force = 'extrap')#, binset='array-like')
obs_gi = Observation(LBG_with_IGM, bp_gmos_i, force = 'extrap')#, binset='array-like')
obs_gz = Observation(LBG_with_IGM, bp_gmos_z, force = 'extrap')#,  binset='array-like')

rig=obs_gr.effstim('abmag')- obs_gi.effstim('abmag')
izg=obs_gi.effstim('abmag')- obs_gz.effstim('abmag')
print('redshift',z)
print(izg)
print(rig)
print(obs_gr.effstim('abmag'))
print(obs_gi.effstim('abmag'))
print(obs_gz.effstim('abmag'))
plt.show()

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