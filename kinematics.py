#calculates z, e, and jp (aka L or angm) for a sample of stars, with confidence intervals

import numpy as np
from astropy.io import ascii
from galpy.actionAngle import UnboundError
from astropy import units
from astropy.coordinates import SkyCoord
from astropy.coordinates import ICRS, Galactic, FK4, FK5
from astropy.coordinates import Angle, Latitude, Longitude
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

radcon = np.pi/180  # radian conversion factor

#loading bar
def tenstep(num):
	return num % 10 == 0

#Source file
#t = ascii.read("sdB.dat")
#t = ascii.read("koposov/data.csv")
#t = ascii.read("sdA_vt500.dat")
#t = ascii.read("raddi_novae/J1603.csv")
t = ascii.read("raddi_novae/J1825.csv")

#define potential
bp= PowerSphericalPotentialwCutoff(alpha=-1.8,rc=1.9*units.kpc, amp=0.5*10**10./bovy_conversion.mass_in_msol(220.,8.))
mp= MiyamotoNagaiPotential(a=3.*units.kpc,b=0.28*units.kpc, amp=6.8*10**10./bovy_conversion.mass_in_msol(220.,8.))
nfwp= NFWPotential(a=16*units.kpc, amp=1.5*10**12./bovy_conversion.mass_in_msol(220.,8.))
MWPotential= [bp,mp,nfwp]
MWPotentialwBH= [MWPotential,KeplerPotential(amp=4*10**6./bovy_conversion.mass_in_msol(220.,8.))]

#Define time and steps
ts= np.linspace(0.,5000.,20000)*units.Myr

#Monte Carlo Sample Size
n = 1000

#confidence interval boundaries
upper = int(n*0.84)
lower = int(n*0.16)

#arrays for saving
sid = [] #identification
ra = np.zeros(len(t))
dec = np.zeros(len(t))

e = np.zeros(len(t)) #eccentricity
e_err_u = np.zeros(len(t))
e_err_d = np.zeros(len(t))

angm = np.zeros(len(t)) #angular momentum
angm_err_u = np.zeros(len(t))
angm_err_d = np.zeros(len(t))

zh = np.zeros(len(t)) #z height
zh_err_u = np.zeros(len(t))
zh_err_d = np.zeros(len(t))

#arrays for calculation
e_calc = np.zeros(n)
L_calc = np.zeros(n)
orbit = np.zeros(n, dtype=Orbit)

#start at
i = 0
for line in t:

	#Read out table
	source_id = t['source_id'][i]

	RV = np.random.normal(t['RV'][i], t['RV_error'][i], n)

	mean = [t['ra'][i], t['dec'][i], t['parallax'][i], t['pmra'][i], t['pmdec'][i]]
	cov = [[t['ra_error'][i]**2, t['ra_dec_corr'][i]*t['ra_error'][i]*t['dec_error'][i], t['ra_parallax_corr'][i]*t['ra_error'][i]*t['parallax_error'][i] ,t['ra_pmra_corr'][i]*t['ra_error'][i]*t['pmra_error'][i], t['ra_pmdec_corr'][i]*t['ra_error'][i]*t['pmdec_error'][i]], \
		[t['ra_dec_corr'][i]*(t['ra_error'][i]*t['dec_error'][i]), t['dec_error'][i]**2, t['dec_parallax_corr'][i]*(t['dec_error'][i]*t['parallax_error'][i]), t['dec_pmra_corr'][i]*(t['dec_error'][i]*t['pmra_error'][i]), t['dec_pmdec_corr'][i]*(t['dec_error'][i]*t['pmdec_error'][i])], \
		[t['ra_parallax_corr'][i]*(t['ra_error'][i]*t['parallax_error'][i]), t['dec_parallax_corr'][i]*(t['dec_error'][i]*t['parallax_error'][i]), t['parallax_error'][i]**2, t['parallax_pmra_corr'][i]*t['parallax_error'][i]*t['pmra_error'][i], t['parallax_pmdec_corr'][i]*t['parallax_error'][i]*t['pmdec_error'][i]], \
		[t['ra_pmra_corr'][i]*(t['ra_error'][i]*t['pmra_error'][i]), t['dec_pmra_corr'][i]*(t['dec_error'][i]*t['pmra_error'][i]), t['parallax_pmra_corr'][i]*t['parallax_error'][i]*t['pmra_error'][i], t['pmra_error'][i]**2, t['pmra_pmdec_corr'][i]*t['pmra_error'][i]*t['pmdec_error'][i]], \
		[t['ra_pmdec_corr'][i]*(t['ra_error'][i]*t['pmdec_error'][i]), t['dec_pmdec_corr'][i]*(t['dec_error'][i]*t['pmdec_error'][i]), t['parallax_pmdec_corr'][i]*t['parallax_error'][i]*t['pmdec_error'][i], t['pmra_pmdec_corr'][i]*t['pmra_error'][i]*t['pmdec_error'][i], t['pmdec_error'][i]**2]]
	RA, Dec, pi, pmra, pmdec = np.random.multivariate_normal(mean, cov, size=n).T

	sid.append(source_id)
	ra[i] = t['ra'][i]
	dec[i] = t['dec'][i]

	c = SkyCoord(ra=RA, dec=Dec, unit='deg')
	c_gal = c.transform_to('galactic')

	l = c_gal.l.deg
	b = c_gal.b.deg

	#distance = np.random.normal(8884, 11, n)*units.pc #(1000./pi)*units.pc
	#distance = np.random.normal(1.77, 0.34, size=1000)*units.kpc #J1603
	distance = np.random.normal(0.93, 0.05, size=1000)*units.kpc #J1825

	#print(distance)

	Vlos = RV #+ 10 * np.cos(l * radcon) * np.cos(b * radcon) + 7.2 * np.sin(b * radcon) + (voS + 5.2) * (np.sin(l * radcon) * np.cos(b * radcon))


	#calculate z w/ conf int
	z = (distance/units.pc * np.sin(b * radcon))/1000
	zh[i] = np.median(z)
	z.sort()
	zh_err_u[i] = z[upper] - np.median(z)
	zh_err_d[i] = np.median(z) - z[lower]

	j = 0
	for j in range(n):
		#calculate analytic e estimate, catch any 'unbound' orbits
		try:
			orbit[j]=Orbit(vxvv=[RA[j],Dec[j],distance[j],pmra[j],pmdec[j],Vlos[j]],radec=True,ro=roS,vo=voS)
			e_calc[j] = orbit[j].e(analytic=True, type='adiabatic', pot=MWPotentialwBH, c=True)
		except UnboundError:
		#parameters cannot be estimated analytically
			e_calc[j] = np.nan

		#integrate the orbit and return the numerical e value
		#orbit[j].integrate(ts, MWPotentialwBH)
		#e_calc[j] = orbit[j].e(analytic=False)

		#print(orbit.L(ro=roS,vo=voS)[0,2])
		L_calc[j] = orbit[j].L(ro=roS,vo=voS)[0,2]		
		j = j+1
	
	#calculate median and conf int for e
	e[i] = np.median(e_calc)
	e_calc.sort()
	e_err_u[i] = e_calc[upper] - np.median(e_calc)
	e_err_d[i] = np.median(e_calc) - e_calc[lower]
	#calculate median and conf int for L
	angm[i] = np.median(L_calc)
	L_calc.sort()
	angm_err_u[i] = L_calc[upper] - np.median(L_calc)
	angm_err_d[i] = np.median(L_calc) - L_calc[lower]


	i=i+1
	if tenstep(i):
		print(i/len(t) * 100 , '%')

#save data
data = [sid, ra, dec, zh, zh_err_u, zh_err_d, e, e_err_u, e_err_d, angm, angm_err_u, angm_err_d]
ascii.write(data, 'raddi_novae/kine_J1825_2.dat', names=['source_id', 'ra', 'dec', 'z', 'z_u', 'z_d', 'e', 'e_u', 'e_d', 'L', 'L_u', 'L_d'], overwrite=True)

