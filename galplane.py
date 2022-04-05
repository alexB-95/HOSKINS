#find points where the star crosses the galactic plane as well as the crossing times
#can differentiate direction, but not between different crossing points and random scatter
#finding the time of crossing points does not work yet

# Importing packages
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import patches
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
from galpy.potential import MWPotential2014
from galpy.util import bovy_conversion

#Source file
#t = ascii.read("LP_40-365/lp40-365.csv")
#t = ascii.read("koposov/data.csv")
#t = ascii.read("PG1610/PG1610.csv")
#t = ascii.read("raddi_novae/J1603.csv")
#t = ascii.read("raddi_novae/J1825.csv")
t = ascii.read("J1211/J1211.csv")

#Choose row and number of random draws
#print("Choose row from 1 to " + str(len(t['ra'])))
i = 0 #int(input())-1
n = 1000

#should the axes be x-y or y-x?
yx = False

#solar values and other constants
#Schoenrich (2012): v_c = 238 +/- 9 km/s, R_0 = 8.27+/-0.29 kpc
voS= np.random.normal(238, 9, size=n)*(units.km/units.s) # 238.#
roS= np.random.normal(8.27, 0.29, size=n)*units.kpc # 8.27#

cf = 3.6 * 10 ** 6 #conversion factor mas - deg

radcon = np.pi/180 #radian conversion factor

#Define potential
bp= PowerSphericalPotentialwCutoff(alpha=-1.8,rc=1.9*units.kpc, amp=0.5*10**10./bovy_conversion.mass_in_msol(238.,8.27))
mp= MiyamotoNagaiPotential(a=3.*units.kpc,b=0.28*units.kpc, amp=6.8*10**10./bovy_conversion.mass_in_msol(238.,8.27))
nfwp= NFWPotential(a=16*units.kpc, amp=1.2*10**12./bovy_conversion.mass_in_msol(238.,8.27))
MWPotential= [bp,mp,nfwp]
MWPotentialwBH= [MWPotential,KeplerPotential(amp=4*10**6./bovy_conversion.mass_in_msol(238.,8.27))]

MWPotential2014wBH= [MWPotential2014,KeplerPotential(amp=4*10**6./bovy_conversion.mass_in_msol(238.,8.27))]

#Define time and steps
print("Choose integration end in Myr")
end = int(input())
steps = 5000
ts= np.linspace(0.,end,steps)*units.Myr


#Read out table and generate values
source_id = t['source_id'][i]

RV = np.random.normal(t['RV'][i], t['RV_error'][i], n)*(units.km / units.s)

#########
goodGaia = False

if goodGaia:
    mean = [t['ra'][i]*cf, t['dec'][i]*cf,
            t['parallax'][i],
            t['pmra'][i], t['pmdec'][i]]
    cov = [[t['ra_error'][i]**2,
            t['ra_dec_corr'][i]*t['ra_error'][i]*t['dec_error'][i],
            t['ra_parallax_corr'][i]*t['ra_error'][i]*t['parallax_error'][i],
            t['ra_pmra_corr'][i]*t['ra_error'][i]*t['pmra_error'][i],
            t['ra_pmdec_corr'][i]*t['ra_error'][i]*t['pmdec_error'][i]],
           [t['ra_dec_corr'][i]*(t['ra_error'][i]*t['dec_error'][i]),
            t['dec_error'][i]**2,
            t['dec_parallax_corr'][i]*(t['dec_error'][i]*t['parallax_error'][i]),
            t['dec_pmra_corr'][i]*(t['dec_error'][i]*t['pmra_error'][i]),
            t['dec_pmdec_corr'][i]*(t['dec_error'][i]*t['pmdec_error'][i])],
           [t['ra_parallax_corr'][i]*(t['ra_error'][i]*t['parallax_error'][i]),
            t['dec_parallax_corr'][i]*(t['dec_error'][i]*t['parallax_error'][i]),
            t['parallax_error'][i]**2,
            t['parallax_pmra_corr'][i]*t['parallax_error'][i]*t['pmra_error'][i],
            t['parallax_pmdec_corr'][i]*t['parallax_error'][i]*t['pmdec_error'][i]],
           [t['ra_pmra_corr'][i]*(t['ra_error'][i]*t['pmra_error'][i]),
            t['dec_pmra_corr'][i]*(t['dec_error'][i]*t['pmra_error'][i]),
            t['parallax_pmra_corr'][i]*t['parallax_error'][i]*t['pmra_error'][i],
            t['pmra_error'][i]**2,
            t['pmra_pmdec_corr'][i]*t['pmra_error'][i]*t['pmdec_error'][i]],
           [t['ra_pmdec_corr'][i]*(t['ra_error'][i]*t['pmdec_error'][i]),
            t['dec_pmdec_corr'][i]*(t['dec_error'][i]*t['pmdec_error'][i]),
            t['parallax_pmdec_corr'][i]*t['parallax_error'][i]*t['pmdec_error'][i],
            t['pmra_pmdec_corr'][i]*t['pmra_error'][i]*t['pmdec_error'][i],
            t['pmdec_error'][i]**2]]

    RA, Dec, pi, pmra, pmdec = np.random.multivariate_normal(mean, cov, size=n).T
    pi = pi*units.mas
   
else:
    mean = [t['ra'][i]*cf, t['dec'][i]*cf,
            t['pmra'][i], t['pmdec'][i]]
    cov = [[t['ra_error'][i]**2,
            t['ra_dec_corr'][i]*t['ra_error'][i]*t['dec_error'][i],
            t['ra_pmra_corr'][i]*t['ra_error'][i]*t['pmra_error'][i],
            t['ra_pmdec_corr'][i]*t['ra_error'][i]*t['pmdec_error'][i]],
           [t['ra_dec_corr'][i]*(t['ra_error'][i]*t['dec_error'][i]),
            t['dec_error'][i]**2,
            t['dec_pmra_corr'][i]*(t['dec_error'][i]*t['pmra_error'][i]),
            t['dec_pmdec_corr'][i]*(t['dec_error'][i]*t['pmdec_error'][i])],
           [t['ra_pmra_corr'][i]*(t['ra_error'][i]*t['pmra_error'][i]),
            t['dec_pmra_corr'][i]*(t['dec_error'][i]*t['pmra_error'][i]),
            t['pmra_error'][i]**2,
            t['pmra_pmdec_corr'][i]*t['pmra_error'][i]*t['pmdec_error'][i]],
           [t['ra_pmdec_corr'][i]*(t['ra_error'][i]*t['pmdec_error'][i]),
            t['dec_pmdec_corr'][i]*(t['dec_error'][i]*t['pmdec_error'][i]),
            t['pmra_pmdec_corr'][i]*t['pmra_error'][i]*t['pmdec_error'][i],
            t['pmdec_error'][i]**2]]

    RA, Dec, pmra, pmdec = np.random.multivariate_normal(mean, cov, size=n).T
#########

#ATTENTION: RA/Dec are in deg in the Gaia catalog but the errors are in mas, requiring some conversions back and forth:
RA = (RA/cf)*units.deg
Dec = (Dec/cf)*units.deg
pmra = pmra*(units.mas / units.yr)
pmdec = pmdec*(units.mas / units.yr)

#transform coordinates from radec to galactic
c = SkyCoord(ra=RA, dec=Dec, unit='deg')
c_gal = c.transform_to('galactic')
l = c_gal.l.deg*units.deg
b = c_gal.b.deg*units.deg

#print("l = ", np.median(l), " , b = ", np.median(b))

#calculate Vlos and distance
Vlos = RV #(RV/(units.km / units.s) + 10 * np.cos(l) * np.cos(b) + 7.2 * np.sin(b) + (voS/(units.km / units.s) + 5.2) * (np.sin(l) * np.cos(b)))*(units.km / units.s)

#distance = np.random.normal(0.632, 0.014, size=n)*units.kpc #LP 40-365
#distance = np.random.normal(8.884, 0.011, size=n)*units.kpc #Koposov
#distance = np.random.normal(17.30, 2.480, size=n)*units.kpc #PG 1610
#distance = np.random.normal(1.77, 0.34, size=n)*units.kpc #J1603
#distance = np.random.normal(0.93, 0.05, size=n)*units.kpc #J1825
distance = np.random.normal(5.4, 0.5, size=n)*units.kpc #J1211
#distance = (1000./pi)*units.pc #general

gcx = np.median(roS/(units.kpc) - distance/units.kpc*np.cos(b)*np.cos(l))
gcy = np.median(distance/units.kpc*np.cos(b)*np.sin(l))
gcz = np.median(distance/units.kpc*np.sin(b))

print(gcx, gcy, gcz)

#initiate arrays for the crossing points (and times)
cross_x_up = []
cross_x_down = []
cross_y_up = []
cross_y_down = []
cross_t_up = []
cross_t_down = []

i_x = []
i_y = []
i_z = []

#flip orbit for backwards integration?
#WARNING: if orbit is flipped, crosses from under the plane are actually from above and vice-versa, as we are integrating backwards in time!
flip = True

#"loading bar"
def tenstep(num):
	return num % 100 == 0

#Function:
#o= Orbit(vxvv=[RA,Dec,distance,pmRA,pmDec,Vlos],radec=True)
# RA and Dec are expressed in degrees
# the distance is expressed in kpc,
# proper motions are expressed in mas/yr (pmra = pmra' * cos[Dec] )
# Vlos is the heliocentric line-of-sight velocity given in km/s

#orbit integration
for j in range(n):
	o=Orbit(vxvv=[RA[j],Dec[j],distance[j],pmra[j],pmdec[j],Vlos[j]],radec=True,ro=roS[j],vo=voS[j])

	if flip==True:
		o = o.flip()

	o.integrate(ts,MWPotentialwBH)

	#o.plot([o.R()],[o.z()],'ro')
	#o.plot3d()	

	#look for the crossing points
	xs = o.x(ts)
	ys = o.y(ts)
	zs = o.z(ts)
	time = o.time(ts)

	shiftzs = np.roll(zs,-1)
	indx_up = (zs[:-1] < 0.)*(shiftzs[:-1] > 0.)
	indx_down = (zs[:-1] > 0.)*(shiftzs[:-1] < 0.)
	cross_x_up.extend(xs[:-1][indx_up])
	cross_x_down.extend(xs[:-1][indx_down])
	cross_y_up.extend(ys[:-1][indx_up])
	cross_y_down.extend(ys[:-1][indx_down])
	cross_t_up.extend(time[:-1][indx_up]/(1*units.Myr))
	cross_t_down.extend(time[:-1][indx_down]/(1*units.Myr))	

	i_x.append(xs[0])
	i_y.append(ys[0])
	i_z.append(zs[0])

	if tenstep(j):
		print(j/n * 100 , '%')

#print("Now:", np.median(i_x), np.median(i_y), np.median(i_z))
#print("Up:", np.median(cross_x_up), np.median(cross_y_up))
#print("Cross:", np.median(cross_x_down), np.median(cross_y_down))

#vup = 'up'
#vdown = 'down'

#if flip==False:
#	up = np.full(len(cross_t_up), 'up')
#	down = np.full(len(cross_t_down), 'down')
#else:
#	up = np.full(len(cross_t_up), 'down')
#	down = np.full(len(cross_t_down), 'up')

#ar_up = np.column_stack((cross_t_up, up))
#ar_down = np.column_stack((cross_t_down, down))

#data=np.concatenate((ar_up, ar_down), axis = 0)
#ascii.write(data, 'galplane' + str(source_id) + '.dat', names=['crossing_time', 'direction'], overwrite=True)
#ascii.write(time, 'time.dat', overwrite=True)

#PLOT:
#define solar circle, disk edge
sc = patches.Circle((0,0), 8)
sc.set_fill(False)
sc.set_linestyle(':')
dc = patches.Circle((0,0), 20)
dc.set_fill(False)
dc.set_linestyle(':')

fig, ax = plt.subplots(1)

if yx == True:
	ax.scatter(cross_y_up, cross_x_up, color="r", marker=".")
	ax.scatter(cross_y_down, cross_x_down, color="b", marker=".")
	ax.set(xlabel='y [kpc]', ylabel='x [kpc]')
else:
	ax.scatter(cross_x_up, cross_y_up, color="r", marker=".")
	ax.scatter(cross_x_down, cross_y_down, color="b", marker=".")
	ax.set(xlabel='x [kpc]', ylabel='y [kpc]')

ax.plot(0, 0, ".k")
ax.add_artist(sc)
ax.add_artist(dc)
ax.axis('equal')

plt.title(source_id)
plt.show()




