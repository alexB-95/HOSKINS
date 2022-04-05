import numpy as np
import matplotlib.pyplot as plt
#from matplotlib import patches
from matplotlib.patches import Ellipse, Polygon
from astropy import units
from astropy.io import ascii

#read out data
t = ascii.read("results/kine_err.dat")


jz = t['L']/1000
jz_u = t['L_u']/1000
jz_d = t['L_d']/1000
e = t['e']
e_u = t['e_u']
e_d = t['e_d']


#error bars
jz_e = [jz_d, jz_u]

e_e = [e_d, e_u]


#define region B
p1 = [0.27, 1100/1000]
p2 = [0.7, 400/1000]
p3 = [0.7, 1150/1000]
p4 = [0.27, 1800/1000]

fig, ax = plt.subplots(figsize=(7.97, 7.97))

ax.set_xlim(0, 1)
ax.set_ylim(-6, 7)

ax.errorbar(e, jz, xerr=e_e, yerr=jz_e, color='k', ecolor='0.5', fmt='D', elinewidth=1, capsize=1.5, markersize=2, zorder=-1) 
ax.add_patch(Polygon([p1, p2, p3, p4], closed=True, fill=False, lw=1, color='r'))

#labels
ax.set(ylabel='J$_z$ [1000 kpc km s$^{-1}$]')
ax.set(xlabel='e')
ax.set_xticks(np.arange(0, 1.1, step=0.1))
ax.set_yticks(np.arange(-6, 7.5, step=1))


#ax.legend()

#plt.show()

#save plot
plt.savefig('jze.pdf')
plt.close()
