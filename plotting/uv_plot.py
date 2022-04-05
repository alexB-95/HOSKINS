#Makes scatter plots of U,V,W values, including contours for the disk and halo, and error bars
#Two plots: V-U, V-W

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import patches
from astropy import units
from astropy.io import ascii


usol = 13.84 #-8.5
vsol = 238 #+ 13.38

#read out data
t = ascii.read("results/UVW.dat")

u = t['u_median'] #+ usol
v = t['v_median'] #+ vsol
uls = t['u_err_down']
uus = t['u_err_up']
vls = t['v_err_down']
vus = t['v_err_up']

#define ellipses
#thin disk, position
U_md = -18 + 13.84
V_md = 14 + 238
#thin disk, standard deviation
sU_md = 38
sV_md = 25
#thick disk, position
U_mh = -40 + 13.84
V_mh = 63 + 238
#thick disk, standard deviation
sU_mh = 56
sV_mh = 45

#thin disk in vu
eu1 = patches.Ellipse((V_md, U_md), 6*sV_md, 6*sU_md, color='r')
eu1.set_fill(False)
eu1.set_linestyle(':')

#thick disk in vu
eu2 = patches.Ellipse((V_mh, U_mh), 6*sV_mh, 6*sU_mh, color='g')
eu2.set_fill(False)
eu2.set_linestyle(':')


#plot
fig, axu = plt.subplots(1, sharex=True, squeeze=True, figsize=(7.97, 7.97))

axu.set(xlabel='V [km s$^{-1}$]', ylabel='U [km s$^{-1}$]')
axu.set_xlim(0, 500)
axu.set_ylim(-250, 250)
axu.set_xticks(np.arange(0, 550, step=100))
axu.set_yticks(np.arange(-200, 250, step=100))

#error bars
uerr = [uls, uus]
verr = [vls, vus]

axu.errorbar(v, u, xerr=[vls, vus], yerr=[uls, uus], color='k', ecolor='0.5', fmt='D', elinewidth=1, capsize=1.5, markersize=2, zorder=-1)
#axu.errorbar(vb, ub, xerr=[vlsb, vusb], yerr=[ulsb, uusb], color='g', ecolor='g', fmt='.', capthick=2, label='sdB')

#add ellipses
axu.add_artist(eu1)
axu.add_artist(eu2)

#axu.legend()


#save plot
plt.savefig('UV.pdf')
plt.close()
