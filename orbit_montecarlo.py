#This calculates n random orbits from the given data and errors via a Monte Carlo method for a single star and plots them in the same diagram
#Now with correlation!

# Importing packages
import numpy as np
import matplotlib.pyplot as plt
from astropy import units
from astropy.coordinates import SkyCoord
from astropy.coordinates import ICRS, Galactic, FK4, FK5
from astropy.coordinates import Angle, Latitude, Longitude
from astropy.io import ascii
from galpy.orbit import Orbit
from galpy.potential import MiyamotoNagaiPotential
from galpy.potential import PowerSphericalPotentialwCutoff
from galpy.potential import NFWPotential
from galpy.potential import KeplerPotential
from galpy.util import bovy_conversion

#Default solar values
#Schoenrich (2012): v_c = 238 +/- 9 km/s, R_0 = 8.27+/-0.29 kpc
voS= 238.
roS= 8.27

radcon = 3.1415926/180  # radian conversion factor

#Source file
t = ascii.read("data_BHB2")

#define potential
bp= PowerSphericalPotentialwCutoff(alpha=-1.8,rc=1.9*units.kpc, amp=0.5*10**10./bovy_conversion.mass_in_msol(220.,8.))
mp= MiyamotoNagaiPotential(a=3.*units.kpc,b=0.28*units.kpc, amp=6.8*10**10./bovy_conversion.mass_in_msol(220.,8.))
nfwp= NFWPotential(a=16*units.kpc, amp=9*10**11./bovy_conversion.mass_in_msol(220.,8.))
MWPotential= [bp,mp,nfwp]
MWPotentialwBH= [MWPotential,KeplerPotential(amp=4*10**6./bovy_conversion.mass_in_msol(220.,8.))]

#Choose row and number of random draws
print("Choose row from 1 to " + str(len(t['ra'])))
i = int(input())-1

n = 25

#Read out table and produce random values from gaussian distributions
source_id = t['distance_type'][i]

mean = [t['ra'][i], t['dec'][i], t['parallax'][i], t['pmra'][i], t['pmdec'][i]]
cov = [[t['ra_error'][i]**2, t['ra_dec_corr'][i]*t['ra_error'][i]*t['dec_error'][i], t['ra_parallax_corr'][i]*t['ra_error'][i]*t['parallax_error'][i] ,t['ra_pmra_corr'][i]*t['ra_error'][i]*t['pmra_error'][i], t['ra_pmdec_corr'][i]*t['ra_error'][i]*t['pmdec_error'][i]], \
	[t['ra_dec_corr'][i]*(t['ra_error'][i]*t['dec_error'][i]), t['dec_error'][i]**2, t['dec_parallax_corr'][i]*(t['dec_error'][i]*t['parallax_error'][i]), t['dec_pmra_corr'][i]*(t['dec_error'][i]*t['pmra_error'][i]), t['dec_pmdec_corr'][i]*(t['dec_error'][i]*t['pmdec_error'][i])], \
	[t['ra_parallax_corr'][i]*(t['ra_error'][i]*t['parallax_error'][i]), t['dec_parallax_corr'][i]*(t['dec_error'][i]*t['parallax_error'][i]), t['parallax_error'][i]**2, t['parallax_pmra_corr'][i]*t['parallax_error'][i]*t['pmra_error'][i], t['parallax_pmdec_corr'][i]*t['parallax_error'][i]*t['pmdec_error'][i]], \
	[t['ra_pmra_corr'][i]*(t['ra_error'][i]*t['pmra_error'][i]), t['dec_pmra_corr'][i]*(t['dec_error'][i]*t['pmra_error'][i]), t['parallax_pmra_corr'][i]*t['parallax_error'][i]*t['pmra_error'][i], t['pmra_error'][i]**2, t['pmra_pmdec_corr'][i]*t['pmra_error'][i]*t['pmdec_error'][i]], \
	[t['ra_pmdec_corr'][i]*(t['ra_error'][i]*t['pmdec_error'][i]), t['dec_pmdec_corr'][i]*(t['dec_error'][i]*t['pmdec_error'][i]), t['parallax_pmdec_corr'][i]*t['parallax_error'][i]*t['pmdec_error'][i], t['pmra_pmdec_corr'][i]*t['pmra_error'][i]*t['pmdec_error'][i], t['pmdec_error'][i]**2]]
RA, Dec, par, pmra, pmdec = np.random.multivariate_normal(mean, cov, size=n).T
RV = np.random.normal(t['RV'][i], t['RV_err'][i], n)

#transform coordinates from radec to galactic
c = SkyCoord(ra=RA, dec=Dec, unit='deg')
c_gal = c.transform_to('galactic')
l = c_gal.l.deg * np.pi/180
b = c_gal.b.deg * np.pi/180

#calculate Vlos and distance
Vlos = RV + 10 * np.cos(l * radcon) * np.cos(b * radcon) + 7.2 * np.sin(b * radcon) + (voS + 5.2) * (np.sin(l * radcon) * np.cos(b * radcon))
distance = np.random.normal(t['d'][i], t['d_err'][i], n)*units.kpc

#Function:
#o= Orbit(vxvv=[RA,Dec,distance,pmRA,pmDec,Vlos],radec=True)
# RA and Dec are expressed in degrees
# the distance is expressed in kpc,
# proper motions are expressed in mas/yr (pmra = pmra' * cos[Dec] )
# Vlos is the heliocentric line-of-sight velocity given in km/s

#initialize orbit array
o = np.zeros(n, dtype=Orbit)

#initialize diagram
fig, ax = plt.subplots(1, sharex=True, squeeze=True, figsize=(7, 7))
ax.set_ylim(-20, 20)
ax.set_xlim(-5, 100)

ts= np.linspace(0.,2000.,20000)*units.Myr #Define time and steps

#Plot n orbits, taking the j-th entry for each variable and adding the result to the orbit array
for j in range(n):
	o[j]=Orbit(vxvv=[RA[j],Dec[j],distance[j],pmra[j],pmdec[j],Vlos[j]],radec=True,ro=roS,vo=voS)
	#o[j]=o[j].flip()

#Integrate orbits over time
	o[j].integrate(ts,MWPotentialwBH)

#Plot the orbits
	o[j].plot(overplot=True)
	ax.plot([o[j].R()],[o[j].z()],'ro')

#Save the plot
plt.savefig(str(source_id) + '_err.pdf')
plt.close()
