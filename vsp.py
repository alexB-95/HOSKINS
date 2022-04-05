
import numpy as np
#from matplotlib import patches
from astropy import units
from astropy.io import ascii

tpm = ascii.read("data.dat")
tuv = ascii.read("UVW.dat")

v_uv = np.zeros(len(tuv))
v_uv_err = np.zeros(len(tuv))
v_pm = np.teros(len(tpm))
v_pm_err = np.teros(len(tpm))

sid = []

i = 0
for line in tpm:
    sid.append(tpm['source_id'][i])
    
    U = tuv['u_median'][i]
    U_err = (tuv['u_err_up'][i]+tuv['u_err_dn'][i])/2
    V = tuv['v_median'][i]
    V_err = (tuv['v_err_up'][i]+tuv['v_err_dn'][i])/2
    W = tuv['w_median'][i]
    W_err = (tuv['w_err_up'][i]+tuv['w_err_dn'][i])/2
    
    RV = tpm['RV'][i]
    RV_err = tpm['RV_err'][i]
    
    pi = np.random(tpm['parallax'][i], tpm['parallax_error'][i], 1000)
    dist = 1/pi
    dist.sort()
    d = np.median(dist)
    d_err_up = dist[840] - np.median(dist)
    d_err_down = np.median(dist) - dist[160]
    d_err = (d_err_up + d_err_down)/2
    
    pmra = tpm[''][i]
    pmra_err = tpm[''][i]
    pmdec = tpm[''][i]
    pmdec_err = tpm[''][i]
    
    v_uv[i] = np.sqrt(U**2 + V**2 + W**2)
    v_uv_err[i] = 1/v_uv * (abs(U)*U_err + abs(V)*V_err + abs(W)*W_err)

    v_pm[i] = np.sqrt(RV**2 + (pmra/d)**2 + (pmdec/d)**2)
    v_pm_err[i] = 1/v_pm * (abs(RV)*RV_err + 1/d*((1/d * pmra_err - abs(pmra)/d * d_err)+((1/d * pmdec_err - abs(pmdec)/d * d_err))))
    
    i += 1

data = [sid, v_uv, v_uv_err, v_pm, v_pm_err]
ascii.write(data, 'results/vsp.dat', names=['source_id', 'v_uvw', 'v_uvw_err', 'v_pm', 'v_pm_err'])

