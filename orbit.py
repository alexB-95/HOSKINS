#Calculates a single orbit for each of the stars in the source file

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
#from galpy.potential import MWPotential2014
from galpy.potential import KeplerPotential
from galpy.util import bovy_conversion

#from galpy.potential import plotPotentials

#Default solar values
#Schoenrich (2012): v_c = 238 +/- 9 km/s, R_0 = 8.27+/-0.29 kpc
voS= 238.
roS= 8.27

radcon = 3.1415926/180

#msol = 1 * bovy_conversion.mass_in_msol(220.,8.)

#Source file
t = ascii.read("data/data.dat")

#Initialize orbit
#Define potential
bp= PowerSphericalPotentialwCutoff(alpha=-1.8,rc=1.9*units.kpc, amp=0.5*10**10./bovy_conversion.mass_in_msol(220.,8.))
mp= MiyamotoNagaiPotential(a=3.*units.kpc,b=0.28*units.kpc, amp=6.8*10**10./bovy_conversion.mass_in_msol(220.,8.))
nfwp= NFWPotential(a=16*units.kpc, amp=1.5*10**12./bovy_conversion.mass_in_msol(220.,8.))
MWPotential= [bp,mp,nfwp]
MWPotentialwBH= [MWPotential,KeplerPotential(amp=4*10**6./bovy_conversion.mass_in_msol(220.,8.))]

#plotPotentials(MWPotentialwBH,rmin=0.01)
#plt.show()

#Original MWPotential2014:
#bp= PowerSphericalPotentialwCutoff(alpha=1.8,rc=1.9/8.,normalize=0.05)
#mp= MiyamotoNagaiPotential(a=3./8.,b=0.28/8.,normalize=.6)
#nfwp= NFWPotential(a=16/8.,normalize=.35)
#MWPotential2014= [bp,mp,nfwp]
#MWPotential2014wBH=[MWPotential2014,KeplerPotential(amp=4*10**6./bovy_conversion.mass_in_msol(220.,8.))]

#Define time and steps
ts= np.linspace(0.,5000.,20000)*units.Myr

#Choose i
i = 0
for line in t:

#Read out table
	source_id = t['Name'][i]
	RA = t['ra'][i]
	Dec = t['dec'][i]
	par = t['parallax'][i]
	RV = t['RV'][i]
	pmra = t['pmRA'][i]
	pmdec = t['pmDE'][i]

	c = SkyCoord(ra=RA, dec=Dec, unit='deg')
	c_gal = c.transform_to('galactic')

	l = c_gal.l.deg
	b = c_gal.b.deg

	Vlos = RV + 10 * np.cos(l * radcon) * np.cos(b * radcon) + 7.2 * np.sin(b * radcon) + (voS + 5.2) * (np.sin(l * radcon) * np.cos(b * radcon))

	distance = (1 / (par/1000))*units.pc

#Function:
#o= Orbit(vxvv=[RA,Dec,distance,pmRA,pmDec,Vlos],radec=True)
# RA and Dec are expressed in degrees
# the distance is expressed in kpc,
# proper motions are expressed in mas/yr (pmra = pmra' * cos[Dec] )
# Vlos is the heliocentric line-of-sight velocity given in km/s

#Initial conditions
	o=Orbit(vxvv=[RA,Dec,distance,pmra,pmdec,Vlos],radec=True,ro=roS,vo=voS)
	#o=o.flip()

#Integrate orbit over time
	o.integrate(ts,MWPotentialwBH)

#Plot the orbit
	o.plot(d1='R', d2='z')
	plt.plot([o.R()],[o.z()],'ro')
	plt.title(source_id)
	plt.axis('equal')
#plt.show()

#Save the plot
	plt.savefig('figures/zRsamescale/' + str(source_id) + '.pdf')
	plt.close()
	
	i += 1
	
	print(i)
