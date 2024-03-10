import os
import sys
import csv
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
import astropy
from astropy.io import fits
import emcee
import corner
import time

start = time.time()

#Read the PANSTARRS and MESA_PANSTARRS data.
p1 = np.genfromtxt('tf2_oid_stack_sshah1502.csv', delimiter = ',', skip_header=1)
p2 = np.genfromtxt('interpolated_phoenix.txt')
p3 = np.genfromtxt('tf2-result.csv', delimiter=',', skip_header=1)

gaia_source_id = p3[:,0]; gaia_ra = p3[:,1]; gaia_ra_error = p3[:,2]; gaia_dec = p3[:,3]
gaia_dec_error = p3[:,4]; gaia_parallax = p3[:,5]; gaia_parallax_error = p3[:,6]; gaia_pm = p3[:,7]
gaia_pm_ra = p3[:,8]; gaia_pm_ra_error = p3[:,9]; gaia_pm_dec = p3[:,10]
gaia_pm_dec_error = p3[:,11]; gaia_ruwe = p3[:,12]



teffp2 = p2[:,0]
loggp2 = p2[:,2]
fehp2 = p2[:,1]
kurucz_sdss_gp2 = p2[:,3]
kurucz_sdss_rp2 = p2[:,4]
kurucz_sdss_ip2 = p2[:,5]
kurucz_sdss_zp2 = p2[:,6]
kurucz_sdss_yp2 = p2[:,7]
kurucz_2mass_jp2 = p2[:,8]
kurucz_2mass_hp2 = p2[:,9]
kurucz_2mass_kp2 = p2[:,10]

teff = []
logg = []
feh = []
kurucz_sdss_g = []
kurucz_sdss_r = []
kurucz_sdss_i = []
kurucz_sdss_z = []
kurucz_sdss_y = []
kurucz_2mass_j = []
kurucz_2mass_h = []
kurucz_2mass_k = []

indp2 = np.where((teffp2>2800)& (teffp2<5000) & (fehp2<-1.5) & (loggp2>3.0))[0]
print('a=', len(indp2), len(teffp2))
teff = np.append(teff, teffp2[indp2])
#teff = np.append(teff, teffp3[indp3])
logg = np.append(logg, loggp2[indp2])
#logg = np.append(logg, loggp3[indp3])
feh = np.append(feh, fehp2[indp2])
#feh = np.append(feh, fehp3[indp3])
kurucz_sdss_g = np.append(kurucz_sdss_g, kurucz_sdss_gp2[indp2])
#kurucz_sdss_g = np.append(kurucz_sdss_g, kurucz_sdss_gp3[indp3])
kurucz_sdss_r = np.append(kurucz_sdss_r, kurucz_sdss_rp2[indp2])
#kurucz_sdss_r = np.append(kurucz_sdss_r, kurucz_sdss_rp3[indp3])
kurucz_sdss_i = np.append(kurucz_sdss_i, kurucz_sdss_ip2[indp2])
#kurucz_sdss_i = np.append(kurucz_sdss_i, kurucz_sdss_ip3[indp3])
kurucz_sdss_z = np.append(kurucz_sdss_z, kurucz_sdss_zp2[indp2])
#kurucz_sdss_z = np.append(kurucz_sdss_z, kurucz_sdss_zp3[indp3])
kurucz_sdss_y = np.append(kurucz_sdss_y, kurucz_sdss_yp2[indp2])
#kurucz_sdss_y = np.append(kurucz_sdss_y, kurucz_sdss_yp3[indp3])
kurucz_2mass_j = np.append(kurucz_2mass_j, kurucz_2mass_jp2[indp2])
#kurucz_2mass_j = np.append(kurucz_2mass_j, kurucz_2mass_jp3[indp3])
kurucz_2mass_h = np.append(kurucz_2mass_h, kurucz_2mass_hp2[indp2])
#kurucz_2mass_h = np.append(kurucz_2mass_h, kurucz_2mass_hp3[indp3])
kurucz_2mass_k = np.append(kurucz_2mass_k, kurucz_2mass_kp2[indp2])
#kurucz_2mass_k = np.append(kurucz_2mass_h, kurucz_2mass_kp3[indp3])
ptsi=[]
for i in range(len(teff)):
    pti = np.where((teffp2==teff[i]) & (loggp2==logg[i])&(fehp2==feh[i]))[0]
    ptsi=np.append(ptsi,pti)

ptsi = np.int64(ptsi)

teffgp2 = np.delete(teffp2,ptsi)
logggp2 = np.delete(loggp2,ptsi)
fehgp2 = np.delete(fehp2, ptsi)
kurucz_sdss_gp2 = np.delete(kurucz_sdss_gp2,ptsi)
kurucz_sdss_rp2 = np.delete(kurucz_sdss_rp2,ptsi)
kurucz_sdss_ip2 = np.delete(kurucz_sdss_ip2,ptsi)
kurucz_sdss_zp2 = np.delete(kurucz_sdss_zp2,ptsi)
kurucz_sdss_yp2 = np.delete(kurucz_sdss_yp2,ptsi)
kurucz_2mass_jp2 = np.delete(kurucz_2mass_jp2,ptsi)
kurucz_2mass_hp2 = np.delete(kurucz_2mass_hp2,ptsi)
kurucz_2mass_kp2 = np.delete(kurucz_2mass_kp2,ptsi)

print('length of first set of Phoenix models=', len(teffgp2))

indp2 = np.where((teffgp2>2800)& (teffgp2<4000) & (fehgp2>-0.5) & (logggp2<3.0))[0]
print('b=', len(indp2), len(teff))
teff = np.append(teff, teffgp2[indp2])
#teff = np.append(teff, teffp3[indp3])
logg = np.append(logg, logggp2[indp2])
#logg = np.append(logg, loggp3[indp3])
feh = np.append(feh, fehgp2[indp2])
#feh = np.append(feh, fehp3[indp3])
kurucz_sdss_g = np.append(kurucz_sdss_g, kurucz_sdss_gp2[indp2])
#kurucz_sdss_g = np.append(kurucz_sdss_g, kurucz_sdss_gp3[indp3])
kurucz_sdss_r = np.append(kurucz_sdss_r, kurucz_sdss_rp2[indp2])
#kurucz_sdss_r = np.append(kurucz_sdss_r, kurucz_sdss_rp3[indp3])
kurucz_sdss_i = np.append(kurucz_sdss_i, kurucz_sdss_ip2[indp2])
#kurucz_sdss_i = np.append(kurucz_sdss_i, kurucz_sdss_ip3[indp3])
kurucz_sdss_z = np.append(kurucz_sdss_z, kurucz_sdss_zp2[indp2])
#kurucz_sdss_z = np.append(kurucz_sdss_z, kurucz_sdss_zp3[indp3])
kurucz_sdss_y = np.append(kurucz_sdss_y, kurucz_sdss_yp2[indp2])
#kurucz_sdss_y = np.append(kurucz_sdss_y, kurucz_sdss_yp3[indp3])
kurucz_2mass_j = np.append(kurucz_2mass_j, kurucz_2mass_jp2[indp2])
#kurucz_2mass_j = np.append(kurucz_2mass_j, kurucz_2mass_jp3[indp3])
kurucz_2mass_h = np.append(kurucz_2mass_h, kurucz_2mass_hp2[indp2])
#kurucz_2mass_h = np.append(kurucz_2mass_h, kurucz_2mass_hp3[indp3])
kurucz_2mass_k = np.append(kurucz_2mass_k, kurucz_2mass_kp2[indp2])
#kurucz_2mass_k = np.append(kurucz_2mass_h, kurucz_2mass_kp3[indp3])
print('c=',len(teff), len(teffgp2[indp2]))

p2 = np.genfromtxt('final_interpolated_ckurucz.txt')

teffp2 = p2[:,0]
loggp2 = p2[:,2]
fehp2 = p2[:,1]
kurucz_sdss_gp2 = p2[:,3]
kurucz_sdss_rp2 = p2[:,4]
kurucz_sdss_ip2 = p2[:,5]
kurucz_sdss_zp2 = p2[:,6]
kurucz_sdss_yp2 = p2[:,7]
kurucz_2mass_jp2 = p2[:,8]
kurucz_2mass_hp2 = p2[:,9]
kurucz_2mass_kp2 = p2[:,10]


indp2 = np.where((teffp2>4000))[0]
print('d=', len(indp2))
teff = np.append(teff, teffp2[indp2])
#teff = np.append(teff, teffp3[indp3])
logg = np.append(logg, loggp2[indp2])
#logg = np.append(logg, loggp3[indp3])
feh = np.append(feh, fehp2[indp2])
#feh = np.append(feh, fehp3[indp3])
kurucz_sdss_g = np.append(kurucz_sdss_g, kurucz_sdss_gp2[indp2])
#kurucz_sdss_g = np.append(kurucz_sdss_g, kurucz_sdss_gp3[indp3])
kurucz_sdss_r = np.append(kurucz_sdss_r, kurucz_sdss_rp2[indp2])
#kurucz_sdss_r = np.append(kurucz_sdss_r, kurucz_sdss_rp3[indp3])
kurucz_sdss_i = np.append(kurucz_sdss_i, kurucz_sdss_ip2[indp2])
#kurucz_sdss_i = np.append(kurucz_sdss_i, kurucz_sdss_ip3[indp3])
kurucz_sdss_z = np.append(kurucz_sdss_z, kurucz_sdss_zp2[indp2])
#kurucz_sdss_z = np.append(kurucz_sdss_z, kurucz_sdss_zp3[indp3])
kurucz_sdss_y = np.append(kurucz_sdss_y, kurucz_sdss_yp2[indp2])
#kurucz_sdss_y = np.append(kurucz_sdss_y, kurucz_sdss_yp3[indp3])
kurucz_2mass_j = np.append(kurucz_2mass_j, kurucz_2mass_jp2[indp2])
#kurucz_2mass_j = np.append(kurucz_2mass_j, kurucz_2mass_jp3[indp3])
kurucz_2mass_h = np.append(kurucz_2mass_h, kurucz_2mass_hp2[indp2])
#kurucz_2mass_h = np.append(kurucz_2mass_h, kurucz_2mass_hp3[indp3])
kurucz_2mass_k = np.append(kurucz_2mass_k, kurucz_2mass_kp2[indp2])
#kurucz_2mass_k = np.append(kurucz_2mass_h, kurucz_2mass_kp3[indp3])
print('e=', len(teff))

kurucz_sdss_gr = kurucz_sdss_g - kurucz_sdss_r
kurucz_sdss_ri = kurucz_sdss_r - kurucz_sdss_i
kurucz_sdss_gi = kurucz_sdss_g - kurucz_sdss_i
kurucz_sdss_gz = kurucz_sdss_g - kurucz_sdss_z
kurucz_sdss_gy = kurucz_sdss_g - kurucz_sdss_y
kurucz_sdss_rz = kurucz_sdss_r - kurucz_sdss_z
kurucz_sdss_ry = kurucz_sdss_r - kurucz_sdss_y
kurucz_sdss_iz = kurucz_sdss_i - kurucz_sdss_z
kurucz_sdss_iy = kurucz_sdss_i - kurucz_sdss_y
kurucz_sdss_zy = kurucz_sdss_z - kurucz_sdss_y


#Transformation equations to convert Kurucz SDSS magnitudes to PANSTARRS magnitudes. Adopted from Tonry et. al. 2012

kurucz_ps_g = kurucz_sdss_g# - 0.012*kurucz_sdss_gr - 0.139
kurucz_ps_r = kurucz_sdss_r# - 0.07*kurucz_sdss_gr
kurucz_ps_i = kurucz_sdss_i# - 0.014*kurucz_sdss_gr + 0.004
kurucz_ps_z = kurucz_sdss_z# - 0.07*kurucz_sdss_gr
kurucz_ps_y = kurucz_sdss_y# - 0.014*kurucz_sdss_gr + 0.004

#Transformation equations to convert Kurucz 2MASS magnitudes to UKIDSS magnitudes. Adopted from Hodgkin et. al. 2009
kurucz_ukidss_j = kurucz_2mass_j# - 0.065*(kurucz_2mass_j - kurucz_2mass_h)
kurucz_ukidss_h = kurucz_2mass_h
kurucz_ukidss_k = kurucz_2mass_k

kurucz_ps_gr = kurucz_sdss_gr
kurucz_ps_ri = kurucz_sdss_ri
kurucz_ps_gi = kurucz_sdss_gi
kurucz_ps_gz = kurucz_sdss_gz
kurucz_ps_gy = kurucz_sdss_gy
kurucz_ps_rz = kurucz_sdss_rz
kurucz_ps_ry = kurucz_sdss_ry
kurucz_ps_iz = kurucz_sdss_iz
kurucz_ps_iy = kurucz_sdss_iy
kurucz_ps_zy = kurucz_sdss_zy

#ebv is defined constant. But can be changed accoding to each field.
ebv = 0.0761
err_ebv = 0.022
ebv = 0.0405
aj = 0.062
ah = 0.039
ak = 0.026

#reading ra, dec and corresponding uncertainties from PANSTARRS file.
ps1_objid = p1[:,0]; ps_ra = p1[:,1]; ps_dec = p1[:,2];ps_ra_err = (p1[:,3]);ps_dec_err = (p1[:,4])

#We read the gmag, rmag and imag. Then we remove the nan values. PANSTARRS has -999 value while UKIDSS has -9.99e-08 value. kron magnitudes can also be used.

gmag = p1[:,5];gkron = p1[:,7];e_gmag = p1[:,6];e_gkron = p1[:,8];rmag = p1[:,9];rkron = p1[:,11]
e_rmag = p1[:,10];e_rkron = p1[:,12];imag = p1[:,13];ikron = p1[:,15];e_imag = p1[:,14];e_ikron = p1[:,16]

zmag = p1[:,17]; zkron = p1[:,19]; e_zmag = p1[:,18]; e_zkron = p1[:,20]; ymag = p1[:,21]
ykron = p1[:,23]; e_ymag = p1[:,22]; e_ykron = p1[:,24]

objinfoflag = p1[:,25]; qualityflag = p1[:,26]; ndetections = p1[:,27]; nstackdetections = p1[:,28]
ginfoflag = p1[:,29]; ginfoflag2 = p1[:,30]; ginfoflag3 = p1[:,31]
rinfoflag = p1[:,29]; rinfoflag2 = p1[:,30]; rinfoflag3 = p1[:,31]
iinfoflag = p1[:,29]; iinfoflag2 = p1[:,30]; iinfoflag3 = p1[:,31]
zinfoflag = p1[:,29]; zinfoflag2 = p1[:,30]; zinfoflag3 = p1[:,31]
yinfoflag = p1[:,29]; yinfoflag2 = p1[:,30]; yinfoflag3 = p1[:,31]

oid1 = [*set(ps1_objid)]
print(len(ps1_objid), len(oid1))
ptsf = []
for i in range(len(oid1)):
    ptsi = np.where(oid1[i]==ps1_objid)[0]
    if len(ptsi)>1.0:
        ptsi = ptsi[0]
    ptsf = np.append(ptsf, ptsi)
ptsf = np.int64(ptsf)
print(len(ptsf))
ps1_objid = ps1_objid[ptsf]; ps_ra = ps_ra[ptsf]; ps_dec = ps_dec[ptsf];ps_ra_err = ps_ra_err[ptsf];ps_dec_err = ps_dec_err[ptsf]

#We read the gmag, rmag and imag. Then we remove the nan values. PANSTARRS has -999 value while UKIDSS has -9.99e-08 value. kron magnitudes can also be used.

gmag = gmag[ptsf]; gkron = gkron[ptsf];e_gmag = e_gmag[ptsf];e_gkron = e_gkron[ptsf];rmag = rmag[ptsf];rkron = rkron[ptsf]
e_rmag = e_rmag[ptsf];e_rkron = e_rkron[ptsf];imag = imag[ptsf];ikron = ikron[ptsf];e_imag = e_imag[ptsf];e_ikron = e_ikron[ptsf]

zmag = zmag[ptsf]; zkron = zkron[ptsf]; e_zmag = e_zmag[ptsf]; e_zkron = e_zkron[ptsf]; ymag = ymag[ptsf]
ykron = ykron[ptsf]; e_ymag = e_ymag[ptsf]; e_ykron = e_ykron[ptsf]

objinfoflag = objinfoflag[ptsf]; qualityflag = qualityflag[ptsf]; ndetections = ndetections[ptsf]; nstackdetections = nstackdetections[ptsf]
ginfoflag = ginfoflag[ptsf]; ginfoflag2 = ginfoflag2[ptsf]; ginfoflag3 = ginfoflag3[ptsf]
rinfoflag = rinfoflag[ptsf]; rinfoflag2 = rinfoflag2[ptsf]; rinfoflag3 = rinfoflag3[ptsf]
iinfoflag = iinfoflag[ptsf]; iinfoflag2 = iinfoflag2[ptsf]; iinfoflag3 = iinfoflag3[ptsf]
zinfoflag = zinfoflag[ptsf]; zinfoflag2 = zinfoflag2[ptsf]; zinfoflag3 = zinfoflag3[ptsf]
yinfoflag = yinfoflag[ptsf]; yinfoflag2 = yinfoflag2[ptsf]; yinfoflag3 = yinfoflag3[ptsf]

print('actual no. of rows=', len(gmag))

#---------------------Here we remove the nan values and also add a condition to seperate stars and galaxies in the data. See Chambers et. al. 2016----------------------#
ind0 = np.where((gmag!= -999) & (imag!= -999) & (rmag != -999) & (zmag != -999) & (ymag != -999) & (gkron!= -999) & (ikron!= -999) & (zkron != -999) & (ykron != -999) & (rkron != -999) & (e_gmag != -999) & (e_rmag != -999) & (e_imag != -999) & (e_zmag != -999) & (e_ymag != -999) & (e_gkron != -999) & (e_rkron != -999) & (e_ikron != -999) & (e_zkron != -999) & (e_ykron != -999))[0]
plt.clf()
plt.figure(figsize=(10,10))
plt.errorbar(gmag[ind0], e_gmag[ind0], alpha=0.5, color = 'g', fmt='.',  label='$\sigma_{g}$')
plt.errorbar(rmag[ind0], e_rmag[ind0], alpha=0.5, color = 'r', fmt='.', label='$\sigma_{r}$')
plt.errorbar(imag[ind0], e_imag[ind0], alpha=0.5, color = 'c', fmt='.', label='$\sigma_{i}$')
plt.errorbar(zmag[ind0], e_zmag[ind0], alpha=0.5, color = 'm', fmt='.', label='$\sigma_{z}$')
plt.errorbar(ymag[ind0], e_ymag[ind0], alpha=0.5, color = 'b', fmt='.', label='$\sigma_{y}$')
plt.xlabel('PSF Magnitude', fontsize=18)
plt.ylabel('Errors', fontsize=18)
plt.legend(loc='best', fontsize=18)
plt.grid()
plt.title('Observed PS1 magnitudes vs uncertainties - TF1', fontsize=18)
#plt.tick_params(left = False)
#ax=plt.gca()
#ax.axes.yaxis.set_ticklabels([])
plt.savefig('stack_ps1_errors.png')
plt.clf()

print('ind0=', len(ind0))
plt.clf()
#plt.rcParams['axes.facecolor'] = 'silver'
plt.clf()
plt.figure(figsize=(10,10))
cm = plt.cm.get_cmap('cool')
cm = plt.cm.get_cmap('winter')
z_prim = imag[ind0]
sc_prim = plt.scatter(gmag[ind0] - imag[ind0], rmag[ind0] - imag[ind0], s=0.5, \
c=z_prim, cmap=cm, alpha = 0.3)
plt.xlim(-1,4)
cb = plt.colorbar(sc_prim)
cb.ax.set_ylabel('$i_{psf}$', fontsize=18)
#plt.ylim(12,23)
plt.title('CCD of the sources in the PS1 data', fontsize=15)
plt.savefig('cmd_gi_stacked_data.png')
plt.clf()

plt.clf()
plt.scatter(imag[ind0], imag[ind0] - ikron[ind0], s=0.5, color = 'g', alpha = 0.3)
#plt.xlim(-1,4)
#plt.ylim(12,23)
plt.grid()
plt.xlabel('$i_{psf}$')
plt.ylabel('$i_{psf}$ - $i_{kron}$')
plt.ylim(-1,2)
plt.savefig('ipsf_ikron_before_sgc.png')
plt.clf()

plt.clf()
plt.hist(imag[ind0], color = 'b', edgecolor = 'black', alpha = 0.5)
plt.xlabel('$I_{psf}$')
plt.ylabel('Bins')
plt.savefig('hist_imag_for_all_the_sources.png')
plt.clf()


ind = np.where((gmag!= -999) & (imag!= -999) & (rmag != -999) & (zmag != -999) & (ymag != -999) & (e_gmag != -999) & (e_rmag != -999) & (e_imag != -999) & (e_zmag != -999) & (e_ymag != -999) & (gkron!= -999) & (ikron!= -999) & (zkron != -999) & (ykron != -999) & (rkron != -999) & ((imag - ikron) < 0.05)  & ((gmag - gkron) < 0.05) & ((rmag - rkron) < 0.05) & ((zmag - zkron) < 0.05) & ((ymag - ykron) < 0.05)  & (e_gkron != -999) & (e_rkron != -999) & (e_ikron != -999) & (e_zkron != -999) & (e_ykron != -999)& (e_gmag < 0.2) & (e_rmag < 0.2) & (e_imag < 0.2) & (e_zmag < 0.2) & (e_ymag < 0.2))[0]
ind2 = np.where((gmag!= -999) & (imag!= -999) & (rmag != -999) & (zmag != -999) & (ymag != -999) & (e_gmag != -999) & (e_rmag != -999) & (e_imag != -999) & (e_zmag != -999) & (e_ymag != -999) & (gkron!= -999) & (ikron!= -999) & (zkron != -999) & (ykron != -999) & (rkron != -999) & ((imag - ikron) > 0.05)  & ((gmag - gkron) > 0.05) & ((rmag - rkron) > 0.05) & ((zmag - zkron) > 0.05) & ((ymag - ykron) > 0.05)  & (e_gkron != -999) & (e_rkron != -999) & (e_ikron != -999) & (e_zkron != -999) & (e_ykron != -999))[0]

plt.clf()
fig,ax = plt.subplots(1)
#plt.rcParams['axes.facecolor'] = 'silver'
plt.figure(figsize=(10,10))
cm = plt.cm.get_cmap('summer_r')
#cm = plt.cm.get_cmap('winter')
z_prim = imag[ind]; z_sec = imag[ind2]
sc_prim = plt.scatter(gmag[ind] - imag[ind], rmag[ind] - imag[ind], s=0.5, \
c='m', label='Stars')
sc_sec = plt.scatter(gmag[ind2] - imag[ind2], rmag[ind2] - imag[ind2], s=0.5,
c='k', label = 'Galaxies')
plt.xlabel('(g-i)', fontsize=18)
plt.ylabel('(r-i)', fontsize=18)
#cb = plt.colorbar(sc_prim)
#cb.ax.set_ylabel('$i_{psf}$', fontsize=18)
ax.set_yticklabels([])
ax.set_xticklabels([])
plt.legend(loc='best', fontsize=18)
plt.savefig('cmd_giri_after_sgc_stacked_data.png')
plt.clf()

plt.clf()
plt.scatter(imag[ind], imag[ind] - ikron[ind], s=0.5, color = 'm', alpha = 0.3)
plt.scatter(imag[ind2], imag[ind2] - ikron[ind2], s=0.5, color = 'k', alpha = 0.3)
plt.grid()
plt.xlabel('$i_{psf}$')
plt.ylabel('$i_{psf}$ - $i_{kron}$')
plt.ylim(-1,2)
plt.savefig('ipsf_ikron_after_sgc.png')
plt.clf()


print('Number of stars after first cut = ', len(ind), len(ind) -  len(ind0))

ps1_objid = ps1_objid[ind]
actual_ps_ra = ps_ra[ind]
actual_err_ps_ra = ps_ra_err[ind]
actual_ps_dec = ps_dec[ind]
actual_err_ps_dec = ps_dec_err[ind]

actual_gmag = gmag[ind]
actual_rmag = rmag[ind]
actual_imag = imag[ind]
actual_zmag = zmag[ind]
actual_ymag = ymag[ind]

plt.clf()
plt.hist(actual_imag, color = 'b', edgecolor = 'black', alpha = 0.5)
plt.xlabel('$I_{psf}$')
plt.ylabel('Bins')
plt.savefig('hist_imag_after_sgc_cut.png')
plt.clf()


actual_e_gmag = e_gmag[ind]
actual_e_rmag = e_rmag[ind]
actual_e_imag = e_imag[ind]
actual_e_zmag = e_zmag[ind]
actual_e_ymag = e_ymag[ind]

actual_gkron = gkron[ind]
actual_rkron = rkron[ind]
actual_ikron = ikron[ind]
actual_zkron = zkron[ind]
actual_ykron = ykron[ind]

actual_gmag = actual_gmag[~np.isnan(actual_gmag)]
actual_rmag = actual_rmag[~np.isnan(actual_rmag)]
actual_imag = actual_imag[~np.isnan(actual_imag)]
actual_zmag = actual_zmag[~np.isnan(actual_zmag)]
actual_ymag = actual_ymag[~np.isnan(actual_ymag)]

actual_gkron = actual_gkron[~np.isnan(actual_gkron)]
actual_rkron = actual_rkron[~np.isnan(actual_rkron)]
actual_ikron = actual_ikron[~np.isnan(actual_ikron)]
actual_zkron = actual_zkron[~np.isnan(actual_zkron)]
actual_ykron = actual_ykron[~np.isnan(actual_ykron)]

actual_e_gmag = actual_e_gmag[~np.isnan(actual_e_gmag)]
actual_e_rmag = actual_e_rmag[~np.isnan(actual_e_rmag)]
actual_e_imag = actual_e_imag[~np.isnan(actual_e_imag)]
actual_e_zmag = actual_e_zmag[~np.isnan(actual_e_zmag)]
actual_e_ymag = actual_e_ymag[~np.isnan(actual_e_ymag)]

objinfoflag = objinfoflag[ind]; qualityflag = qualityflag[ind];
ndetections = ndetections[ind]; nstackdetections = nstackdetections[ind]
ginfoflag = ginfoflag[ind]; ginfoflag2 = ginfoflag2[ind]; ginfoflag3 = ginfoflag3[ind]
rinfoflag = rinfoflag[ind]; rinfoflag2 = rinfoflag2[ind]; rinfoflag3 = rinfoflag3[ind]
iinfoflag = iinfoflag[ind]; iinfoflag2 = iinfoflag2[ind]; iinfoflag3 = iinfoflag3[ind]
zinfoflag = zinfoflag[ind]; zinfoflag2 = zinfoflag2[ind]; zinfoflag3 = zinfoflag3[ind]
yinfoflag = yinfoflag[ind]; yinfoflag2 = yinfoflag2[ind]; yinfoflag3 = yinfoflag3[ind]


quadrature_gr = np.sqrt(actual_e_gmag**2 + actual_e_rmag**2)
quadrature_ri = np.sqrt(actual_e_rmag**2 + actual_e_imag**2)
quadrature_gi = np.sqrt(actual_e_gmag**2 + actual_e_imag**2)
quadrature_gy = np.sqrt(actual_e_gmag**2 + actual_e_ymag**2)
quadrature_gz = np.sqrt(actual_e_gmag**2 + actual_e_zmag**2)
quadrature_ry = np.sqrt(actual_e_rmag**2 + actual_e_ymag**2)
quadrature_rz = np.sqrt(actual_e_rmag**2 + actual_e_zmag**2)
quadrature_iy = np.sqrt(actual_e_imag**2 + actual_e_ymag**2)
quadrature_iz = np.sqrt(actual_e_imag**2 + actual_e_zmag**2)
quadrature_zy = np.sqrt(actual_e_zmag**2 + actual_e_ymag**2)

ag = ebv*0.88*(3.613 - 0.0972*(actual_gmag - actual_imag) + 0.01*(actual_gmag - actual_imag)**2)
ar = ebv*0.88*(2.585 - 0.0315*(actual_gmag - actual_imag))
ai = ebv*0.88*(1.908 - 0.0152*(actual_gmag - actual_imag))
az = ebv*0.88*(1.499 - 0.0023*(actual_gmag - actual_imag))
ay = ebv*0.88*(1.251 - 0.0027*(actual_gmag - actual_imag))

de_reddened_actual_gmag = actual_gmag - ag#3.518*ebv
de_reddened_actual_rmag = actual_rmag - ar#2.617*ebv
de_reddened_actual_imag = actual_imag - ai#1.917*ebv
de_reddened_actual_zmag = actual_zmag - az#1.917*ebv
de_reddened_actual_ymag = actual_ymag - ay#1.917*ebv


de_reddened_gr = de_reddened_actual_gmag - de_reddened_actual_rmag
de_reddened_ri = de_reddened_actual_rmag - de_reddened_actual_imag
de_reddened_gi = de_reddened_actual_gmag - de_reddened_actual_imag
de_reddened_gy = de_reddened_actual_gmag - de_reddened_actual_ymag
de_reddened_gz = de_reddened_actual_gmag - de_reddened_actual_zmag
de_reddened_ry = de_reddened_actual_rmag - de_reddened_actual_ymag
de_reddened_rz = de_reddened_actual_rmag - de_reddened_actual_zmag
de_reddened_iy = de_reddened_actual_imag - de_reddened_actual_ymag
de_reddened_iz = de_reddened_actual_imag - de_reddened_actual_zmag
de_reddened_zy = de_reddened_actual_zmag - de_reddened_actual_ymag

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]

err_obs = np.sqrt((quadrature_gr)**2 + (quadrature_ri)**2 + (quadrature_gi)**2 + (quadrature_gy)**2 + (quadrature_gz)**2 + (quadrature_ry)**2 + (quadrature_rz)**2 + (quadrature_iy)**2 + (quadrature_iz)**2 + (quadrature_zy)**2)
#print('net observed error in colours=', err_obs)

def stdv(sfavg, v1, v2, v3, v4, v5):
    mu = sfavg
    n=5
    sig = (((mu - v1)**2 + (mu - v2)**2 + (mu-v3)**2 + (mu - v4)**2 + (mu - v5)**2)/5)**0.5
    return sig

header = ['#ps1_objid', 'ps1_ra', 'ps1_ra_error', 'ps1_dec', 'ps1_dec_error', 'ps1_gpsf', 'ps1_gpsf_error', \
'ps1_rpsf', 'ps1_rpsf_error', 'ps1_ipsf', 'ps1_ipsf_error', 'ps1_zpsf', 'ps1_zpsf_error', 'ps1_ypsf', 'ps1_ypsf_error',\
'teff', 'logg', 'feh', 'sam_g','sam_r','sam_i','sam_z','sam_y','sam_j','sam_h','sam_k', 'scale_factor', 'scale_factor_error',\
'chi2', 'computed_j', 'computed_j_error', 'computed_h', 'computed_h_error', 'computed_k', 'computed_k_error',\
'gaia_source_id', 'gaia_ra', 'gaia_ra_error', 'gaia_dec', 'gaia_dec_error', 'gaia_parallax', 'gaia_parallax_error', 'gaia_pm',\
'gaia_pm_ra', 'gaia_pm_ra_error', 'gaia_pm_dec', 'gaia_pm_dec_error', 'gaia_ruwe',\
'objinfoflag', 'qualityflag', 'ndetections', 'nstackdetections',\
'ginfoflag', 'ginfoflag2', 'ginfoflag3', 'rinfoflag', 'rinfoflag2',\
'rinfoflag3','iinfoflag', 'iinfoflag2', 'iinfoflag3','zinfoflag'\
'zinfoflag2', 'zinfoflag3', 'yinfoflag', 'yinfoflag2', 'yinfoflag3', 'SAM Flag']

sigma_sf = (1/5)*(actual_e_gmag**2+actual_e_rmag**2+actual_e_imag**2+actual_e_zmag**2\
                        +actual_e_ymag**2)

with open('match_parameters_kurucz.csv', 'w', encoding='UTF8') as file1:
    writer=csv.writer(file1)
    writer.writerow(header)
    for j in range(len(de_reddened_gr)):
                    d_min = []
                    d1 = (kurucz_ps_gr - de_reddened_gr[j])#/quadrature_gr[j]
                    d2 = (kurucz_ps_ri -  de_reddened_ri[j])#/quadrature_ri[j]
                    d3 = (kurucz_ps_gi - de_reddened_gi[j])#/quadrature_gi[j]
                    d4 = (kurucz_ps_gz - de_reddened_gz[j])#/quadrature_gz[j]
                    d5 = (kurucz_ps_gy - de_reddened_gy[j])#/quadrature_gy[j]
                    d6 = (kurucz_ps_rz - de_reddened_rz[j])#/quadrature_rz[j]
                    d7 = (kurucz_ps_ry - de_reddened_ry[j])#/quadrature_ry[j]
                    d8 = (kurucz_ps_iz - de_reddened_iz[j])#/quadrature_iz[j]
                    d9 = (kurucz_ps_iy - de_reddened_iy[j])#/quadrature_iy[j]
                    d10 = (kurucz_ps_zy - de_reddened_zy[j])#/quadrature_zy[j]
                    d_rms =  (d1**2 + d2**2 + d3**2 + d4**2 + d5**2 + d6**2 + d7**2 + d8**2 + d9**2 + d10**2)
                    minv = find_nearest(d_rms, np.min(d_rms))


                    if np.min(d_rms) > 0.0:
                        index_minv = np.where(minv == (d_rms))[0]
                        sf1 = ((de_reddened_actual_gmag[j] - kurucz_ps_g[index_minv]))
                        sf2 = ((de_reddened_actual_rmag[j] - kurucz_ps_r[index_minv]))
                        sf3 = ((de_reddened_actual_imag[j] - kurucz_ps_i[index_minv]))
                        sf4 = ((de_reddened_actual_zmag[j] - kurucz_ps_z[index_minv]))
                        sf5 = ((de_reddened_actual_ymag[j] - kurucz_ps_y[index_minv]))

                        sf_avg = (sf1 + sf2 + sf3 + sf4 + sf5)/5.0

                        computed_j = kurucz_ukidss_j[index_minv] + sf_avg + aj - 0.91
                        computed_h = kurucz_ukidss_h[index_minv] + sf_avg + ah - 1.39
                        computed_k = kurucz_ukidss_k[index_minv] + sf_avg + ak -1.85


                        computed_j_error = computed_h_error = computed_k_error = np.sqrt(sigma_sf**2 +err_ebv**2)

                        gamma_gaia = 3600*np.sqrt(((actual_ps_ra[j] - gaia_ra)*np.cos(np.radians(actual_ps_dec[j])))**2\
                         + (actual_ps_dec[j] - gaia_dec)**2)
                        indg = np.where(gamma_gaia<=1.0)[0]

                        if len(indg) > 1.0:
                            g_gaia1 = gamma_gaia[indg]
                            gf_gaia = gamma_gaia[np.where(np.min(g1)==gamma)[0]]

                            indgn = np.where(gf_gaia == gamma)[0]

                            data = ps1_objid[j], actual_ps_ra[j], actual_err_ps_ra[j], actual_ps_dec[j], actual_err_ps_dec[j], \
                            actual_gmag[j], actual_e_gmag[j], actual_rmag[j], actual_e_rmag[j], actual_imag[j], actual_e_imag[j],\
                            actual_zmag[j], actual_e_zmag[j], actual_ymag[j], actual_e_ymag[j], teff[index_minv][0], logg[index_minv][0], \
                            feh[index_minv][0], kurucz_ps_g[index_minv][0],  kurucz_ps_r[index_minv][0], kurucz_ps_i[index_minv][0], \
                            kurucz_ps_z[index_minv][0], kurucz_ps_y[index_minv][0], kurucz_ukidss_j[index_minv][0], kurucz_ukidss_h[index_minv][0],\
                             kurucz_ukidss_k[index_minv][0], sf_avg[0], sigma_sf[j], minv, computed_j[0], computed_j_error[0],\
                             computed_h[0], computed_h_error[0], computed_k[0], computed_k_error[0],\
                             gaia_source_id[indgn][0], gaia_ra[indgn][0], gaia_ra_error[indgn][0], gaia_dec[indgn][0], \
                             gaia_dec_error[indgn][0], gaia_parallax[indgn][0], \
                             gaia_parallax_error[indgn][0], gaia_pm[indgn][0], gaia_pm_ra[indgn][0], \
                             gaia_pm_ra_error[indgn][0], gaia_pm_dec[indgn][0], gaia_pm_dec_error[indgn][0], gaia_ruwe[indgn][0],\
                             objinfoflag[j], qualityflag[j], ndetections[j], nstackdetections[j],\
                             ginfoflag[j], ginfoflag2[j], ginfoflag3[j], rinfoflag[j], rinfoflag2[j],\
                             rinfoflag3[j], iinfoflag[j], iinfoflag2[j], iinfoflag3[j], zinfoflag[j],\
                             zinfoflag2[j], zinfoflag3[j], yinfoflag[j], yinfoflag2[j], yinfoflag3[j]

                            writer.writerow(data)
                        elif len(indg) == 0.0:
                            data = ps1_objid[j], actual_ps_ra[j], actual_err_ps_ra[j], actual_ps_dec[j], actual_err_ps_dec[j], \
                            actual_gmag[j], actual_e_gmag[j], actual_rmag[j], actual_e_rmag[j], actual_imag[j], actual_e_imag[j],\
                            actual_zmag[j], actual_e_zmag[j], actual_ymag[j], actual_e_ymag[j], teff[index_minv][0], logg[index_minv][0], \
                            feh[index_minv][0], kurucz_ps_g[index_minv][0],  kurucz_ps_r[index_minv][0], kurucz_ps_i[index_minv][0], \
                            kurucz_ps_z[index_minv][0], kurucz_ps_y[index_minv][0], kurucz_ukidss_j[index_minv][0], kurucz_ukidss_h[index_minv][0],\
                             kurucz_ukidss_k[index_minv][0], sf_avg[0], sigma_sf[j], minv, computed_j[0], computed_j_error[0],\
                             computed_h[0], computed_h_error[0], computed_k[0], computed_k_error[0],\
                             -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999,\
                             objinfoflag[j], qualityflag[j], ndetections[j], nstackdetections[j],\
                             ginfoflag[j], ginfoflag2[j], ginfoflag3[j], rinfoflag[j], rinfoflag2[j],\
                             rinfoflag3[j], iinfoflag[j], iinfoflag2[j], iinfoflag3[j], zinfoflag[j],\
                             zinfoflag2[j], zinfoflag3[j], yinfoflag[j], yinfoflag2[j], yinfoflag3[j]
                            #writer=csv.writer(file1)
                            #writer.writerow(header)
                            writer.writerow(data)
                        elif len(indg) == 1.0:
                            data = ps1_objid[j], actual_ps_ra[j], actual_err_ps_ra[j], actual_ps_dec[j], actual_err_ps_dec[j], \
                            actual_gmag[j], actual_e_gmag[j], actual_rmag[j], actual_e_rmag[j], actual_imag[j], actual_e_imag[j],\
                            actual_zmag[j], actual_e_zmag[j], actual_ymag[j], actual_e_ymag[j], teff[index_minv][0], logg[index_minv][0], \
                            feh[index_minv][0], kurucz_ps_g[index_minv][0],  kurucz_ps_r[index_minv][0], kurucz_ps_i[index_minv][0], \
                            kurucz_ps_z[index_minv][0], kurucz_ps_y[index_minv][0], kurucz_ukidss_j[index_minv][0], kurucz_ukidss_h[index_minv][0],\
                            kurucz_ukidss_k[index_minv][0], sf_avg[0], sigma_sf[j], minv, computed_j[0], computed_j_error[0],\
                            computed_h[0], computed_h_error[0], computed_k[0], computed_k_error[0],\
                            gaia_source_id[indg][0], gaia_ra[indg][0], \
                            gaia_ra_error[indg][0], gaia_dec[indg][0], gaia_dec_error[indg][0], gaia_parallax[indg][0], \
                            gaia_parallax_error[indg][0], gaia_pm[indg][0], gaia_pm_ra[indg][0], gaia_pm_ra_error[indg][0],\
                            gaia_pm_dec[indg][0], gaia_pm_dec_error[indg][0], gaia_ruwe[indg][0],\
                             objinfoflag[j], qualityflag[j], ndetections[j], nstackdetections[j],\
                             ginfoflag[j], ginfoflag2[j], ginfoflag3[j], rinfoflag[j], rinfoflag2[j],\
                             rinfoflag3[j], iinfoflag[j], iinfoflag2[j], iinfoflag3[j], zinfoflag[j],\
                             zinfoflag2[j], zinfoflag3[j], yinfoflag[j], yinfoflag2[j], yinfoflag3[j]
                            #writer=csv.writer(file1)
                            #writer.writerow(header)
                            writer.writerow(data)
                        #file1.write('%0.16f %0.8f %0.8f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f\n' %(sf_avg, actual_ps_ra[j], actual_ps_dec[j], actual_err_ps_ra[j], actual_err_ps_dec[j], kurucz_ukidss_j[index_minv], actual_gmag[j], actual_rmag[j], actual_imag[j], actual_zmag[j], actual_ymag[j], teff[index_minv], logg[index_minv], feh[index_minv], actual_e_gmag[j], actual_e_rmag[j], actual_e_imag[j], actual_e_zmag[j], actual_e_ymag[j], kurucz_ps_g[index_minv], kurucz_ps_r[index_minv], kurucz_ps_i[index_minv], kurucz_ps_z[index_minv], kurucz_ps_y[index_minv], minv, err_obs[j], stdv(sf_avg, sf1, sf2, sf3, sf4, sf5), computed_j, computed_h, computed_k, kurucz_ukidss_h[index_minv], kurucz_ukidss_k[index_minv]))
#-------------------Now we read the data from match_parameters_kurucz_t1 file which contains the kurucz_j magnitude of the stars from the PANSTARRS data which satisfy the cut in the previous two sections-------------------#

p3 = np.genfromtxt('match_parameters_kurucz.csv', delimiter=',')

ps1_objidf = p3[:,0]
actual_ps_raf = p3[:,1]
actual_ps_decf = p3[:,3]
actual_ps_ra_errf = p3[:,2]*0.0002777778
actual_ps_dec_errf = p3[:,4]*0.0002777778
actual_gmagf = p3[:,5]
actual_e_gmagf = p3[:,6]
actual_rmagf = p3[:,7]
actual_e_rmagf = p3[:,8]
actual_imagf = p3[:,9]
actual_e_imagf = p3[:,10]
actual_zmagf = p3[:,11]
actual_e_zmagf = p3[:,12]
actual_ymagf = p3[:,13]
actual_e_ymagf = p3[:,14]
tefff = p3[:,15]
loggf = p3[:,16]
fehf = p3[:,17]
kgf = p3[:,18]
krf = p3[:,19]
kif = p3[:,20]
kzf = p3[:,21]
kyf = p3[:,22]
kjf = p3[:,23]
khf = p3[:,24]
kkf = p3[:,25]
sf_avgf = p3[:,26]
sf_avg_errf = p3[:,27]
minvrf = p3[:,28]
computed_jf = p3[:,29]
computed_j_errf = p3[:,30]
computed_hf = p3[:,31]
computed_h_errf = p3[:,32]
computed_kf = p3[:,33]
computed_k_errf = p3[:,34]
gaia_source_idf = p3[:,35]
gaia_raf = p3[:,36]
gaia_ra_errf = p3[:,37]
gaia_decf = p3[:,38]
gaia_dec_errf = p3[:,39]
gaia_parallaxf = p3[:,40]
gaia_parallax_errorf = p3[:,41]
gaia_pmf = p3[:,42]
gaia_pm_raf = p3[:,43]
gaia_pm_ra_errorf = p3[:,44]
gaia_pm_decf = p3[:,45]
gaia_pm_dec_errorf = p3[:,46]
gaia_ruwef = p3[:,47]
objinfoflag = p3[:,48]
qualityflag = p3[:,49]
ndetections = p3[:,50]
nstackdetections = p3[:,51]
ginfoflag = p3[:,52]
ginfoflag2 = p3[:,53]
ginfoflag3 = p3[:,54]
rinfoflag = p3[:,55]
rinfoflag2 = p3[:,56]
rinfoflag3 = p3[:,57]
iinfoflag = p3[:,58]
iinfoflag2 = p3[:,59]
iinfoflag3 = p3[:,60]
zinfoflag = p3[:,61]
zinfoflag2 = p3[:,62]
zinfoflag3 = p3[:,63]
yinfoflag = p3[:,64]
yinfoflag2 = p3[:,65]
yinfoflag3 = p3[:,66]

print('faintest magnitude in J band =', np.max(computed_jf))

plt.clf()
plt.scatter(computed_jf, minvrf, s=3, color = 'b', edgecolor = 'black', alpha = 0.5)
plt.ylabel('$\chi^{2}$')
plt.xlabel('Computed J')
plt.ylim(0,10)
plt.savefig('hist_minvrf.png')

plt.clf()
plt.hist(actual_imagf, color = 'b', edgecolor = 'black', alpha = 0.5)
plt.xlabel('$I_{psf}$')
plt.ylabel('Bins')
plt.savefig('hist_imag_match_parameters.png')
plt.clf()

plt.clf()
plt.hist(computed_jf)
plt.xlabel('$J_{computed}$')
plt.ylabel('Bins')
plt.savefig('hist_jmag_computed.png')
plt.clf()


plt.clf()
plt.hist(computed_hf)
plt.xlabel('H_{computed}$')
plt.ylabel('Bins')
plt.savefig('hist_hmag_computed.png')
plt.clf()


plt.clf()
plt.hist(computed_kf)
plt.xlabel('$K_{computed}$')
plt.ylabel('Bins')
plt.savefig('hist_kmag_computed.png')
plt.clf()

print('Number of Stars after second cut =', len(kjf))

apparent_jf = computed_jf
apparent_hf = computed_hf
apparent_kf = computed_kf
err_computed_j = np.sqrt((sf_avg_errf)**2 + (err_ebv)**2)
err_computed_h = np.sqrt((sf_avg_errf)**2 + (err_ebv)**2)
err_computed_k = np.sqrt((sf_avg_errf)**2 + (err_ebv)**2)


plt.clf()
plt.scatter(apparent_jf, minvrf, s=5, alpha = 0.3)
plt.xlim(12,22)
plt.ylim(-0.05, 0.6)
plt.grid()
plt.xlabel('J magnitude')
plt.ylabel('$d_{quad}$')
plt.savefig('drms_vs_j.png')
plt.clf()


#indjn = np.where(apparent_jf<=18.0)[0]
#print('median observed error for stars brighter than 18 magnitude=', np.median(err_quad_obs[indjn]))


hdulist = fits.open('tf6.fits',  memmap=True)
p8 = hdulist[1].data
petro_j = p8['JPETROMAG']
e_petro_j = p8['jPetroMagErr']
petro_h = p8['HPETROMAG']
e_petro_h = p8['hPetroMagErr']
petro_k = p8['KPETROMAG']
e_petro_k = p8['kPetroMagErr']
t1_ra = (p8['RA'])*(180.0/np.pi)
t1_dec = p8['DEC']*(180.0/np.pi)

#Here we remove the nan values and use the filtered points.

ind3 = np.where((petro_j != -9.99999500e+08) & (e_petro_j != -9.99999500e+08) & (petro_h != -9.99999500e+08) \
& (e_petro_h != -9.99999500e+08) & (petro_k != -9.99999500e+08) & (e_petro_j < 0.2)& (e_petro_h < 0.2)& (e_petro_k < 0.2))[0]

print('Number of Stars in the UKIDSS data = ', len(ind3))

actual_petro_j = petro_j[ind3]
e_petro_h = e_petro_h[ind3]
actual_petro_h = petro_h[ind3]
e_petro_k = e_petro_k[ind3]
actual_petro_k = petro_k[ind3]
e_petro_j = e_petro_j[ind3]
actual_t1_ra = t1_ra[ind3]
actual_t1_dec = t1_dec[ind3]

plt.clf()
plt.hist(actual_petro_j)
plt.xlabel('$J_{petro}$')
plt.ylabel('Bins')
plt.savefig('hist_jpetro_ukidss.png')
plt.clf()

plt.clf()
plt.hist(actual_petro_h)
plt.xlabel('$H_{petro}$')
plt.ylabel('Bins')
plt.savefig('hist_hpetro_ukidss.png')
plt.clf()

plt.clf()
plt.hist(actual_petro_k)
plt.xlabel('$K_{petro}$')
plt.ylabel('Bins')
plt.savefig('hist_kpetro_ukidss.png')
plt.clf()
#Here we match the coordinates of UKIDSS and PANSTARRS. If they match within 1" (0.0002777778 degrees), we compute the difference between the observed and computed J magnitude for the matching stars. We convert the mesa J magnitude to apparent J magnitude by adding sf and reddening.

ob_j = []; ob_h = []; ob_k = []; diff_jf = []; diff_hf = []; diff_kf = []

header = ['ps1_objid', 'ps1_ra', 'ps1_ra_error', 'ps1_dec', 'ps1_dec_error', 'ps1_gpsf', 'ps1_gpsf_error', \
'ps1_rpsf', 'ps1_rpsf_error', 'ps1_ipsf', 'ps1_ipsf_error', 'ps1_zpsf', 'ps1_zpsf_error', 'ps1_ypsf', 'ps1_ypsf_error',\
'teff', 'logg', 'feh', 'sam_g','sam_r','sam_i','sam_z','sam_y','sam_j','sam_h','sam_k', 'scale_factor', 'scale_factor_error',\
'chi2', 'computed_j', 'computed_j_error', 'computed_h', 'computed_h_error', 'computed_k', 'computed_k_error',\
'gaia_source_id', 'gaia_ra', 'gaia_ra_error', 'gaia_dec', 'gaia_dec_error', 'gaia_parallax', 'gaia_parallax_error', 'gaia_pm',\
'gaia_pm_ra', 'gaia_pm_ra_error', 'gaia_pm_dec', 'gaia_pm_dec_error', 'gaia_ruwe', 'ob_j', 'ob_j_error', 'ob_h', \
'ob_h_error', 'ob_k', 'ob_k_error', 'diff_j', 'diff_h', 'diff_k', 'objinfoflag', 'qualityflag', 'ndetections',\
 'nstackdetections', 'ginfoflag', 'ginfoflag2', 'ginfoflag3', 'rinfoflag', 'rinfoflag2',\
'rinfoflag3','iinfoflag', 'iinfoflag2', 'iinfoflag3','zinfoflag',\
'zinfoflag2', 'zinfoflag3', 'yinfoflag', 'yinfoflag2', 'yinfoflag3']

with open('partial_catalogue_kurucz.csv', 'w', encoding = 'UTF8') as file4:
    writer=csv.writer(file4)
    writer.writerow(header)
    for i1 in range(len(actual_ps_raf)):
            gamma = 3600*np.sqrt(((actual_ps_raf[i1] - actual_t1_ra)*np.cos(np.radians(actual_ps_decf[i1])))**2 + (actual_ps_decf[i1] - actual_t1_dec)**2)
            ind4 = np.where(gamma<=1.0)[0]
            if len(ind4) > 1:
                g1 = gamma[ind4]
                gf = gamma[np.where(np.min(g1) == gamma)[0]]
                ind5 = np.where(gf == gamma)[0]
                ob_j = np.append(ob_j, actual_petro_j[ind5])
                ob_h = np.append(ob_h, actual_petro_h[ind5])
                ob_k = np.append(ob_k, actual_petro_k[ind5])
                diff_j = actual_petro_j[ind5] - apparent_jf[i1]
                diff_h = actual_petro_h[ind5] - apparent_hf[i1]
                diff_k = actual_petro_k[ind5] - apparent_kf[i1]
                diff_jf = np.append(diff_jf, diff_j)
                diff_hf = np.append(diff_hf, diff_h)
                diff_kf = np.append(diff_kf, diff_k)
                data = ps1_objidf[i1], actual_ps_raf[i1], actual_ps_ra_errf[i1], actual_ps_decf[i1], actual_ps_dec_errf[i1], \
                actual_gmagf[i1],actual_e_gmagf[i1], actual_rmagf[i1], actual_e_rmagf[i1], actual_imagf[i1], \
                actual_e_imagf[i1], actual_zmagf[i1], actual_e_zmagf[i1], actual_ymagf[i1], actual_e_ymagf[i1],\
                tefff[i1], loggf[i1], fehf[i1],kgf[i1], krf[i1], kif[i1], kzf[i1], kyf[i1], kjf[i1], khf[i1], kkf[i1], sf_avgf[i1], sf_avg_errf[i1],\
                minvrf[i1], apparent_jf[i1], err_computed_j[i1], apparent_hf[i1], err_computed_h[i1], apparent_kf[i1], \
                err_computed_k[i1], gaia_source_idf[i1], gaia_raf[i1], gaia_ra_errf[i1], gaia_decf[i1], gaia_dec_errf[i1], \
                gaia_parallaxf[i1], gaia_parallax_errorf[i1], gaia_pmf[i1], gaia_pm_raf[i1], gaia_pm_ra_errorf[i1],\
                gaia_pm_decf[i1], gaia_pm_dec_errorf[i1], gaia_ruwef[i1], actual_petro_j[ind5][0][0], e_petro_j[ind5][0][0], \
                actual_petro_h[ind5][0][0], e_petro_h[ind5][0][0], actual_petro_k[ind5][0][0], e_petro_k[ind5][0][0], diff_j[0][0], \
                diff_h[0][0], diff_k[0][0], objinfoflag[i1], qualityflag[i1], ndetections[i1], nstackdetections[i1],\
                ginfoflag[i1], ginfoflag2[i1], ginfoflag3[i1], rinfoflag[i1], rinfoflag2[i1],\
                rinfoflag3[i1], iinfoflag[i1], iinfoflag2[i1], iinfoflag3[i1], zinfoflag[i1],\
                zinfoflag2[i1], zinfoflag3[i1], yinfoflag[i1], yinfoflag2[i1], yinfoflag3[i1]
                writer.writerow(data)
            elif len(ind4) == 1:
                ob_j = np.append(ob_j, actual_petro_j[ind4])
                ob_h = np.append(ob_h, actual_petro_h[ind4])
                ob_k = np.append(ob_k, actual_petro_k[ind4])
                diff_j = actual_petro_j[ind4] - apparent_jf[i1]
                diff_h = actual_petro_h[ind4] - apparent_hf[i1]
                diff_k = actual_petro_k[ind4] - apparent_kf[i1]
                diff_jf = np.append(diff_jf, diff_j)
                diff_hf = np.append(diff_hf, diff_h)
                diff_kf = np.append(diff_kf, diff_k)
                data = ps1_objidf[i1], actual_ps_raf[i1], actual_ps_ra_errf[i1], actual_ps_decf[i1], actual_ps_dec_errf[i1], \
                actual_gmagf[i1], actual_e_gmagf[i1], actual_rmagf[i1], actual_e_rmagf[i1], actual_imagf[i1], \
                actual_e_imagf[i1], actual_zmagf[i1], actual_e_zmagf[i1], actual_ymagf[i1], actual_e_ymagf[i1],\
                tefff[i1], loggf[i1], fehf[i1], kgf[i1], krf[i1], kif[i1], kzf[i1], kyf[i1], kjf[i1], khf[i1], kkf[i1], \
                sf_avgf[i1], sf_avg_errf[i1],\
                minvrf[i1], apparent_jf[i1], err_computed_j[i1], apparent_hf[i1], err_computed_h[i1], apparent_kf[i1],\
                err_computed_k[i1], gaia_source_idf[i1], gaia_raf[i1], gaia_ra_errf[i1], gaia_decf[i1], gaia_dec_errf[i1], \
                gaia_parallaxf[i1], gaia_parallax_errorf[i1], gaia_pmf[i1], gaia_pm_raf[i1], gaia_pm_ra_errorf[i1],\
                gaia_pm_decf[i1], gaia_pm_dec_errorf[i1], gaia_ruwef[i1], actual_petro_j[ind4][0][0], e_petro_j[ind4][0][0],\
                actual_petro_h[ind4][0][0], e_petro_h[ind4][0][0], actual_petro_k[ind4][0][0], \
                e_petro_k[ind4][0][0], diff_j[0][0], diff_h[0][0], diff_k[0][0], objinfoflag[i1], qualityflag[i1], ndetections[i1],\
                nstackdetections[i1],\
                ginfoflag[i1], ginfoflag2[i1], ginfoflag3[i1], rinfoflag[i1], rinfoflag2[i1],\
                rinfoflag3[i1], iinfoflag[i1], iinfoflag2[i1], iinfoflag3[i1], zinfoflag[i1],\
                zinfoflag2[i1], zinfoflag3[i1], yinfoflag[i1], yinfoflag2[i1], yinfoflag3[i1]
                writer.writerow(data)

print('Number of Stars in the scatter plot =', len(ob_j))


indjp = np.where(np.abs(diff_jf)<0.2)[0]
indhp = np.where(np.abs(diff_hf)<0.2)[0]
indkp = np.where(np.abs(diff_kf)<0.2)[0]

bins2 = np.arange(diff_jf.min(), diff_jf.max()+.1, 0.05)
plt.clf()
fig = plt.figure(figsize=(10,10))
from matplotlib.gridspec import GridSpec
gs = GridSpec(4, 4)
ax_joint = fig.add_subplot(gs[1:4,0:3])
ax_marg_x = fig.add_subplot(gs[0,0:3])
ax_marg_y = fig.add_subplot(gs[1:4,3])
ax_joint.scatter(ob_j, diff_jf, alpha = 0.3, s = 5, color = 'g', label = 'No. of stars =' + str(len(ob_j)))
ax_joint.grid()
ax_joint.set_ylim(-2,2)
ax_marg_y.set_ylim(-2,2)
ax_joint.legend(fontsize=18, loc = 'best')
nx, bx, px = ax_marg_x.hist(ob_j, color = 'm', edgecolor = 'g', alpha = 0.5, label = 'Observed J')
ny, by, px = ax_marg_y.hist(diff_jf, bins = bins2, orientation="horizontal", edgecolor = 'g', alpha = 0.5, \
facecolor = 'orange', label = '$J_{o}$ - $J_{c}$')
biny_max = np.where(ny == ny.max())
print('maxbin', "{:.2f}".format(by[biny_max][0]), np.median(diff_jf), np.std(diff_jf))
#plt.text(500, -1.5, str('mode at '"{:.2f}".format(by[biny_max][0])), rotation = 270)
ax_joint.set_xlim(10,22.5)
ax_marg_x.grid()
ax_marg_x.legend(loc='best', fontsize=15)
ax_marg_y.grid()
ax_marg_y.legend(loc='best', fontsize=15)
ax_joint.set_title('Median and spread of the scatter =' + str("{:.3f}".format(np.median(diff_jf)))\
+'$\pm$'+str("{:.3f}".format(np.std(diff_jf))))
ax_marg_x.set_title('No. of sources lying in the range -0.2 < ($J_{o}$ - $J_{c}$) < 0.2 =' + str("{:.2f}".format(100*len(indjp)/len(diff_jf))+'%'))
# Turn off tick labels on marginals
plt.setp(ax_marg_x.get_xticklabels(), visible=False)
plt.setp(ax_marg_y.get_yticklabels(), visible=False)
# Set labels on joint
ax_joint.set_xlabel('Observed J magnitude')
ax_joint.set_ylabel('(Observed - Computed) J magnitude')
# Set labels on marginals
ax_marg_y.set_xlabel('N')
ax_marg_x.set_ylabel('N')
plt.savefig('scatter_whole_kurucz_j.png')
plt.clf()


bins2 = np.arange(diff_hf.min(), diff_hf.max()+.1, 0.05)
plt.clf()
fig = plt.figure(figsize=(10,10))
from matplotlib.gridspec import GridSpec
gs = GridSpec(4, 4)
ax_joint = fig.add_subplot(gs[1:4,0:3])
ax_marg_x = fig.add_subplot(gs[0,0:3])
ax_marg_y = fig.add_subplot(gs[1:4,3])
ax_joint.scatter(ob_h, diff_hf, alpha = 0.3, s = 5, color = 'g', label = 'No. of stars =' + str(len(ob_h)))
ax_joint.grid()
ax_joint.set_ylim(-2,2)
ax_marg_y.set_ylim(-2,2)
ax_joint.legend(fontsize=18, loc = 'best')
nx, bx, px = ax_marg_x.hist(ob_h, color = 'm', edgecolor = 'g', alpha = 0.5, label = 'Observed H')
ny, by, px = ax_marg_y.hist(diff_hf, bins = bins2, orientation="horizontal", edgecolor = 'g', alpha = 0.5, \
facecolor = 'orange', label = '$H_{o}$ - $H_{c}$')
biny_max = np.where(ny == ny.max())
print('maxbin', "{:.2f}".format(by[biny_max][0]), np.median(diff_hf), np.std(diff_hf))
ax_marg_x.set_title('No. of sources lying in the range -0.2 < ($H_{o}$ - $H_{c}$) < 0.2 =' + str("{:.2f}".format(100*len(indhp)/len(diff_hf))+'%'))
ax_joint.set_xlim(10,20)
ax_marg_x.grid()
ax_marg_x.legend(loc='best', fontsize=15)
ax_marg_y.grid()
ax_marg_y.legend(loc='best', fontsize=15)
ax_joint.set_title('Median and spread of the scatter =' + str("{:.3f}".format(np.median(diff_hf)))\
+'$\pm$'+str("{:.3f}".format(np.std(diff_hf))))
# Turn off tick labels on marginals
plt.setp(ax_marg_x.get_xticklabels(), visible=False)
plt.setp(ax_marg_y.get_yticklabels(), visible=False)
# Set labels on joint
ax_joint.set_xlabel('Observed H magnitude')
ax_joint.set_ylabel('(Observed - Computed) H magnitude')
# Set labels on marginals
ax_marg_y.set_xlabel('N')
ax_marg_x.set_ylabel('N')
plt.savefig('scatter_whole_kurucz_h.png')
plt.clf()


bins2 = np.arange(diff_kf.min(), diff_kf.max()+.1, 0.05)
plt.clf()
fig = plt.figure(figsize=(10,10))
from matplotlib.gridspec import GridSpec
gs = GridSpec(4, 4)
ax_joint = fig.add_subplot(gs[1:4,0:3])
ax_marg_x = fig.add_subplot(gs[0,0:3])
ax_marg_y = fig.add_subplot(gs[1:4,3])
ax_joint.scatter(ob_k, diff_kf, alpha = 0.3, s = 5, color = 'g', label = 'No. of stars =' + str(len(ob_k)))
ax_joint.grid()
ax_joint.set_ylim(-2,2)
ax_marg_y.set_ylim(-2,2)
ax_joint.legend(fontsize=18, loc = 'best')
nx, bx, px = ax_marg_x.hist(ob_k, color = 'm', edgecolor = 'g', alpha = 0.5, label = 'Observed K')
ny, by, px = ax_marg_y.hist(diff_kf, bins = bins2, orientation="horizontal", edgecolor = 'g', alpha = 0.5,\
facecolor = 'orange', label = '$K_{o}$ - $K_{c}$')
biny_max = np.where(ny == ny.max())
print('maxbin', "{:.2f}".format(by[biny_max][0]), np.median(diff_kf), np.std(diff_kf))
ax_joint.set_xlim(10,20)
ax_marg_x.grid()
ax_marg_x.legend(loc='best', fontsize=15)
ax_marg_y.grid()
ax_marg_y.legend(loc='best', fontsize=15)
ax_joint.set_title('Median and spread of the scatter =' + str("{:.3f}".format(np.median(diff_kf)))\
+'$\pm$'+str("{:.3f}".format(np.std(diff_kf))))
ax_marg_x.set_title('No. of sources lying in the range -0.2 < ($K_{o}$ - $K_{c}$) < 0.2 =' + str("{:.2f}".format(100*len(indkp)/len(diff_kf))+'%'))
# Turn off tick labels on marginals
plt.setp(ax_marg_x.get_xticklabels(), visible=False)
plt.setp(ax_marg_y.get_yticklabels(), visible=False)
# Set labels on joint
ax_joint.set_xlabel('Observed K magnitude')
ax_joint.set_ylabel('(Observed - Computed) K magnitude')
# Set labels on marginals
ax_marg_y.set_xlabel('N')
ax_marg_x.set_ylabel('N')
plt.savefig('scatter_whole_kurucz_k.png')
plt.clf()


end = time.time()
print(' time taken =' , end - start)


sys.exit(0)
