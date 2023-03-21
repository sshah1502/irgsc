import os
import sys
import csv
from datetime import date
import matplotlib.pylab as pylab
import numpy as np
from ._fitting import Generate_IRGSC
from ._read_data import Read_Data
from matplotlib import pyplot as plt
from ._get_data import Get_Data
from ._read_data import Read_Data

params = {'legend.fontsize': 'x-large',
          'figure.figsize': (10,10),
         'axes.labelsize': 'x-large',
         'axes.titlesize':'x-large',
         'xtick.labelsize':'x-large',
         'ytick.labelsize':'x-large'}
pylab.rcParams.update(params)
from matplotlib.offsetbox import AnchoredText
current_datetime = date.today()
home_dir = os.getcwd()

def find_nearest(array, value):
        array = np.asarray(array)
        idx = (np.abs(array - value)).argmin()
        return array[idx]

header = ['ps1_objid', 'ps1_ra', 'ps1_ra_error', 'ps1_dec', 'ps1_dec_error', 'ps1_gpsf', 'ps1_gpsf_error', \
'ps1_rpsf', 'ps1_rpsf_error', 'ps1_ipsf', 'ps1_ipsf_error', 'ps1_zpsf', 'ps1_zpsf_error', 'ps1_ypsf',\
    'ps1_ypsf_error', 'teff', 'logg', 'feh', 'sam_g','sam_r','sam_i','sam_z','sam_y','sam_j','sam_h',\
        'sam_k', 'scale_factor', 'scale_factor_error', 'chi2', 'computed_j', 'computed_j_error', 'computed_h',\
            'computed_h_error', 'computed_k', 'computed_k_error', 'gaia_source_id', 'gaia_ra', 'gaia_ra_error',\
                'gaia_dec', 'gaia_dec_error', 'gaia_parallax', 'gaia_parallax_error', 'gaia_pm', 'gaia_pm_ra',\
                    'gaia_pm_ra_error', 'gaia_pm_dec', 'gaia_pm_dec_error', 'gaia_ruwe', 'objinfoflag',\
                        'qualityflag', 'ndetections', 'nstackdetections', 'ginfoflag', 'ginfoflag2',\
                            'ginfoflag3', 'rinfoflag', 'rinfoflag2', 'rinfoflag3','iinfoflag', 'iinfoflag2',\
                                'iinfoflag3','zinfoflag', 'zinfoflag2', 'zinfoflag3', 'yinfoflag', 'yinfoflag2',\
                                    'yinfoflag3', 'diff_j', 'diff_h', 'diff_k', 'observed_j', 'obseved_j_error',\
                                        'observed_h', 'observed_h_error', 'observed_k', 'observed_k_error']
class Validate():
    def __init__(self, ra, dec):
        self.ra, self.dec = ra, dec
        self.rd = Read_Data(ra,dec)

    def read_irgsc(self):
        ra_name = str(self.ra).replace('.','_')
        dec_name = str(self.dec).replace('.', '_')
        try:
            irgsc_data = np.genfromtxt('IRGSC' + '_' + 'RA' + str(ra_name) + 'DEC' + str(dec_name) +\
                      str(current_datetime) + '.csv', delimiter=',', skip_header=1)
            ps1_objid = irgsc_data[:,0]; ps_ra = irgsc_data[:,1]; err_ps_ra = irgsc_data[:,2]; ps_dec = irgsc_data[:,3]; \
            err_ps_dec = irgsc_data[:,4]; ec_gmag = irgsc_data[:,5]; e_ec_gmag = irgsc_data[:,6]; ec_rmag = irgsc_data[:,7]; \
            e_ec_rmag = irgsc_data[:,8]; ec_imag = irgsc_data[:,9]; e_ec_imag = irgsc_data[:,10]; ec_zmag = irgsc_data[:,11]; \
            e_ec_zmag = irgsc_data[:,12]; ec_ymag = irgsc_data[:,13]; e_ec_ymag = irgsc_data[:,14]; teff = irgsc_data[:,15]; \
                logg = irgsc_data[:,16]; feh = irgsc_data[:,17]; sam_g = irgsc_data[:,18]; sam_r = irgsc_data[:,19]; sam_i = irgsc_data[:,20]; sam_z = irgsc_data[:,21];\
                sam_y = irgsc_data[:,22]; sam_j = irgsc_data[:,23]; sam_h = irgsc_data[:,24]; sam_k = irgsc_data[:,25]; sf_avg = irgsc_data[:,26]; sigma_sf = irgsc_data[:,27];\
                    min_dquad_element = irgsc_data[:,28]; computed_j = irgsc_data[:,29]; computed_j_error = irgsc_data[:,30]; computed_h = irgsc_data[:,31];\
                    computed_h_error = irgsc_data[:,32]; computed_k = irgsc_data[:,33]; computed_k_error = irgsc_data[:,34];\
                        gaia_source_id = irgsc_data[:,35]; gaia_ra = irgsc_data[:,36]; gaia_ra_error = irgsc_data[:,37];\
                        gaia_dec = irgsc_data[:,38]; gaia_dec_error = irgsc_data[:,39]; gaia_parallax = irgsc_data[:,40];\
                            gaia_parallax_error = irgsc_data[:,41]; gaia_pm = irgsc_data[:,42]; gaia_pm_ra = irgsc_data[:,43]; gaia_pm_ra_error = irgsc_data[:,44];\
                            gaia_pm_dec = irgsc_data[:,45]; gaia_pm_dec_error = irgsc_data[:,46]; gaia_ruwe = irgsc_data[:,47];\
                                objinfoflag = irgsc_data[:,48]; qualityflag = irgsc_data[:,49]; ndetections = irgsc_data[:,50];\
                                    nstackdetections = irgsc_data[:,51]; ginfoflag = irgsc_data[:,52]; ginfoflag2 = irgsc_data[:,53];\
                                        ginfoflag3 = irgsc_data[:,54]; rinfoflag = irgsc_data[:,55]; rinfoflag2 = irgsc_data[:,56];\
                                        rinfoflag3 = irgsc_data[:,57]; iinfoflag = irgsc_data[:,58]; iinfoflag2 = irgsc_data[:,59];\
                                            iinfoflag3 = irgsc_data[:,60]; zinfoflag = irgsc_data[:,61]; zinfoflag2 = irgsc_data[:,62];\
                                            zinfoflag3 = irgsc_data[:,63]; yinfoflag = irgsc_data[:,64]; yinfoflag2 = irgsc_data[:,65];\
                                            yinfoflag3 = irgsc_data[:,66]
            irgsc_data = ps1_objid, ps_ra, err_ps_ra, ps_dec, err_ps_dec, ec_gmag, e_ec_gmag, ec_rmag, e_ec_rmag,\
            ec_imag, e_ec_imag, ec_zmag, e_ec_zmag, ec_ymag, e_ec_ymag, teff, logg, feh, sam_g,\
                sam_r, sam_i, sam_z, sam_y, sam_j, sam_h, sam_k, sf_avg, sigma_sf, min_dquad_element,\
                computed_j, computed_j_error, computed_h, computed_h_error, computed_k, computed_k_error,\
                    gaia_source_id, gaia_ra, gaia_ra_error, gaia_dec, gaia_dec_error, gaia_parallax,\
                    gaia_parallax_error, gaia_pm, gaia_pm_ra, gaia_pm_ra_error, gaia_pm_dec,\
                        gaia_pm_dec_error, gaia_ruwe, objinfoflag, qualityflag, ndetections,\
                        nstackdetections, ginfoflag, ginfoflag2, ginfoflag3, rinfoflag, rinfoflag2,\
                            rinfoflag3, iinfoflag, iinfoflag2, iinfoflag3, zinfoflag, zinfoflag2,\
                            zinfoflag3, yinfoflag, yinfoflag2, yinfoflag3
        except FileNotFoundError:
            Generate_IRGSC.generate_irgsc(self.ra,self.dec)
            irgsc_data = np.genfromtxt('IRGSC' + '_' + 'RA' + str(ra_name) + 'DEC' + str(dec_name) +\
                      str(current_datetime) + '.csv', delimiter=',', skip_header=1)
            ps1_objid = irgsc_data[:,0]; ps_ra = irgsc_data[:,1]; err_ps_ra = irgsc_data[:,2]; ps_dec = irgsc_data[:,3]; \
            err_ps_dec = irgsc_data[:,4]; ec_gmag = irgsc_data[:,5]; e_ec_gmag = irgsc_data[:,6]; ec_rmag = irgsc_data[:,7]; \
            e_ec_rmag = irgsc_data[:,8]; ec_imag = irgsc_data[:,9]; e_ec_imag = irgsc_data[:,10]; ec_zmag = irgsc_data[:,11]; \
            e_ec_zmag = irgsc_data[:,12]; ec_ymag = irgsc_data[:,13]; e_ec_ymag = irgsc_data[:,14]; teff = irgsc_data[:,15]; \
                logg = irgsc_data[:,16]; feh = irgsc_data[:,17]; sam_g = irgsc_data[:,18]; sam_r = irgsc_data[:,19]; sam_i = irgsc_data[:,20]; sam_z = irgsc_data[:,21];\
                sam_y = irgsc_data[:,22]; sam_j = irgsc_data[:,23]; sam_h = irgsc_data[:,24]; sam_k = irgsc_data[:,25]; sf_avg = irgsc_data[:,26]; sigma_sf = irgsc_data[:,27];\
                    min_dquad_element = irgsc_data[:,28]; computed_j = irgsc_data[:,29]; computed_j_error = irgsc_data[:,30]; computed_h = irgsc_data[:,31];\
                    computed_h_error = irgsc_data[:,32]; computed_k = irgsc_data[:,33]; computed_k_error = irgsc_data[:,34];\
                        gaia_source_id = irgsc_data[:,35]; gaia_ra = irgsc_data[:,36]; gaia_ra_error = irgsc_data[:,37];\
                        gaia_dec = irgsc_data[:,38]; gaia_dec_error = irgsc_data[:,39]; gaia_parallax = irgsc_data[:,40];\
                            gaia_parallax_error = irgsc_data[:,41]; gaia_pm = irgsc_data[:,42]; gaia_pm_ra = irgsc_data[:,43]; gaia_pm_ra_error = irgsc_data[:,44];\
                            gaia_pm_dec = irgsc_data[:,45]; gaia_pm_dec_error = irgsc_data[:,46]; gaia_ruwe = irgsc_data[:,47];\
                                objinfoflag = irgsc_data[:,48]; qualityflag = irgsc_data[:,49]; ndetections = irgsc_data[:,50];\
                                    nstackdetections = irgsc_data[:,51]; ginfoflag = irgsc_data[:,52]; ginfoflag2 = irgsc_data[:,53];\
                                        ginfoflag3 = irgsc_data[:,54]; rinfoflag = irgsc_data[:,55]; rinfoflag2 = irgsc_data[:,56];\
                                        rinfoflag3 = irgsc_data[:,57]; iinfoflag = irgsc_data[:,58]; iinfoflag2 = irgsc_data[:,59];\
                                            iinfoflag3 = irgsc_data[:,60]; zinfoflag = irgsc_data[:,61]; zinfoflag2 = irgsc_data[:,62];\
                                            zinfoflag3 = irgsc_data[:,63]; yinfoflag = irgsc_data[:,64]; yinfoflag2 = irgsc_data[:,65];\
                                            yinfoflag3 = irgsc_data[:,66]
            irgsc_data = ps1_objid, ps_ra, err_ps_ra, ps_dec, err_ps_dec, ec_gmag, e_ec_gmag, ec_rmag, e_ec_rmag,\
            ec_imag, e_ec_imag, ec_zmag, e_ec_zmag, ec_ymag, e_ec_ymag, teff, logg, feh, sam_g,\
                sam_r, sam_i, sam_z, sam_y, sam_j, sam_h, sam_k, sf_avg, sigma_sf, min_dquad_element,\
                computed_j, computed_j_error, computed_h, computed_h_error, computed_k, computed_k_error,\
                    gaia_source_id, gaia_ra, gaia_ra_error, gaia_dec, gaia_dec_error, gaia_parallax,\
                    gaia_parallax_error, gaia_pm, gaia_pm_ra, gaia_pm_ra_error, gaia_pm_dec,\
                        gaia_pm_dec_error, gaia_ruwe, objinfoflag, qualityflag, ndetections,\
                        nstackdetections, ginfoflag, ginfoflag2, ginfoflag3, rinfoflag, rinfoflag2,\
                            rinfoflag3, iinfoflag, iinfoflag2, iinfoflag3, zinfoflag, zinfoflag2,\
                            zinfoflag3, yinfoflag, yinfoflag2, yinfoflag3
        return irgsc_data

    def validate(self, validate=True):
        """
            Function to compare the observed and computed NIR magnitudes for a given field.
            If this is set to True, the function first obtains the UKIDSS data for the given field.
            Raises an error if UKIDSS data is not available.
            Raises warning if J and H band data from UKIDSS is not available.
        """
        print("")
        print('###############################################################################')
        print('Now validating the generated IRGSC using the UKIDSS data file if it is present.')
        print("")
        print('###############################################################################')

        if validate is True:
            ukidss_data = self.rd.read_nir_data()
            irgsc_data = self.read_irgsc()
            ukidss_j, ukidss_h, ukidss_k, e_ukidss_j, e_ukidss_h, e_ukidss_k, ukidss_ra, ukidss_dec = ukidss_data
            ps1_objid, ps_ra, err_ps_ra, ps_dec, err_ps_dec, ec_gmag, e_ec_gmag, ec_rmag,\
                e_ec_rmag, ec_imag, e_ec_imag, ec_zmag, e_ec_zmag, ec_ymag, e_ec_ymag, teff,\
                logg, feh, sam_g,sam_r, sam_i, sam_z, sam_y, sam_j, sam_h, sam_k, sf_avg, sigma_sf,\
                    min_dquad_element, computed_j, computed_j_error, computed_h, computed_h_error,\
                        computed_k, computed_k_error, gaia_source_id, gaia_ra, gaia_ra_error,\
                            gaia_dec, gaia_dec_error, gaia_parallax, gaia_parallax_error,\
                            gaia_pm, gaia_pm_ra, gaia_pm_ra_error, gaia_pm_dec, gaia_pm_dec_error,\
                                gaia_ruwe, objinfoflag, qualityflag, ndetections, nstackdetections,\
                                ginfoflag, ginfoflag2, ginfoflag3, rinfoflag, rinfoflag2, rinfoflag3,\
                                    iinfoflag, iinfoflag2, iinfoflag3, zinfoflag, zinfoflag2,\
                                    zinfoflag3, yinfoflag, yinfoflag2, yinfoflag3 = irgsc_data
        
        validate_params=[]; ob_j = []; e_ob_j = []; ob_h = []; e_ob_h = []; ob_k = []; e_ob_k = []
        diff_jf = []; diff_hf = []; diff_kf = []
        ra_name = str(self.ra).replace('.','_'); dec_name = str(self.dec).replace('.', '_')
        with open('validated_IRGSC' + '_' + 'RA' + str(ra_name) + '_' + 'DEC' + str(dec_name) + '_' + str(current_datetime) + '.csv', 'w') as file2:
            writer=csv.writer(file2)
            writer.writerow(header)
            for i1 in range(len(ps_ra)):
                #positionally matching the sources in the UKIDSS within 1" to the PS1 sources in the catalogue
                gamma_ukidss = 3600*np.sqrt(((ps_ra[i1] - ukidss_ra)*np.cos(np.radians(ps_dec[i1])))**2\
                                      + (ps_dec[i1] - ukidss_dec)**2)
                index_ukidss_position = np.where(gamma_ukidss<=1.0)[0]
                if len(index_ukidss_position) > 1:
                    matched_positions = gamma_ukidss[index_ukidss_position]
                    minimum_seperation = gamma_ukidss[np.where(np.min(matched_positions) == gamma_ukidss)[0]]
                    index_minimum_seperation = np.where(minimum_seperation == gamma_ukidss)[0]
                    diff_j = ukidss_j[index_minimum_seperation] - computed_j[i1]
                    diff_h = ukidss_h[index_minimum_seperation] - computed_h[i1]
                    diff_k = ukidss_k[index_minimum_seperation] - computed_k[i1]
                    validate_params = ps1_objid[i1], ps_ra[i1], err_ps_ra[i1], ps_dec[i1], err_ps_dec[i1], ec_gmag[i1], e_ec_gmag[i1], ec_rmag[i1],\
                    e_ec_rmag[i1], ec_imag, e_ec_imag[i1], ec_zmag, e_ec_zmag[i1], ec_ymag, e_ec_ymag[i1], teff[i1], logg[i1], feh[i1], sam_g[i1],\
                    sam_r[i1], sam_i[i1], sam_z[i1], sam_y[i1], sam_j[i1], sam_h[i1], sam_k[i1], sf_avg[i1], sigma_sf[i1], min_dquad_element[i1],\
                    computed_j[i1], computed_j_error[i1], computed_h[i1], computed_h_error[i1], computed_k[i1], computed_k_error[i1], gaia_source_id[i1],\
                    gaia_ra[i1], gaia_ra_error[i1], gaia_dec[i1], gaia_dec_error[i1], gaia_parallax[i1], gaia_parallax_error[i1], gaia_pm[i1], gaia_pm_ra[i1],\
                    gaia_pm_ra_error[i1], gaia_pm_dec[i1], gaia_pm_dec_error[i1], gaia_ruwe[i1], objinfoflag[i1], qualityflag[i1], ndetections[i1], nstackdetections[i1],\
                    ginfoflag[i1], ginfoflag2[i1], ginfoflag3[i1], rinfoflag[i1], rinfoflag2[i1], rinfoflag3[i1], iinfoflag[i1], iinfoflag2[i1], iinfoflag3[i1],\
                    zinfoflag[i1], zinfoflag2[i1], zinfoflag3[i1], yinfoflag[i1], yinfoflag2[i1], yinfoflag3[i1], diff_j, diff_h, diff_k,\
                    ukidss_j[index_minimum_seperation], e_ukidss_j[index_minimum_seperation], ukidss_h[index_minimum_seperation], e_ukidss_h[index_minimum_seperation],\
                    ukidss_k[index_minimum_seperation], e_ukidss_k[index_minimum_seperation]
                    writer.writerow(validate_params)
                    ob_j = np.append(ob_j, ukidss_j[index_minimum_seperation])
                    e_ob_j = np.append(e_ob_j, e_ukidss_j[index_minimum_seperation])
                    ob_h = np.append(ob_h, ukidss_h[index_minimum_seperation])
                    e_ob_h = np.append(e_ob_h, e_ukidss_h[index_minimum_seperation])
                    ob_k = np.append(ob_k, ukidss_k[index_minimum_seperation])
                    e_ob_k = np.append(e_ob_k, e_ukidss_k[index_minimum_seperation])
                    diff_jf = np.append(diff_jf, diff_j)
                    diff_hf = np.append(diff_hf, diff_h)
                    diff_kf = np.append(diff_kf, diff_k)
                elif len(index_ukidss_position) == 1:
                    diff_j = ukidss_j[index_ukidss_position] - computed_j[i1]
                    diff_h = ukidss_h[index_ukidss_position] - computed_h[i1]
                    diff_k = ukidss_k[index_ukidss_position] - computed_k[i1]
                    validate_params = ps1_objid[i1], ps_ra[i1], err_ps_ra[i1], ps_dec[i1], err_ps_dec[i1], ec_gmag[i1], e_ec_gmag[i1], ec_rmag[i1],\
                    e_ec_rmag[i1], ec_imag, e_ec_imag[i1], ec_zmag, e_ec_zmag[i1], ec_ymag, e_ec_ymag[i1], teff[i1], logg[i1], feh[i1], sam_g[i1],\
                    sam_r[i1], sam_i[i1], sam_z[i1], sam_y[i1], sam_j[i1], sam_h[i1], sam_k[i1], sf_avg[i1], sigma_sf[i1], min_dquad_element[i1],\
                    computed_j[i1], computed_j_error[i1], computed_h[i1], computed_h_error[i1], computed_k[i1], computed_k_error[i1], gaia_source_id[i1],\
                    gaia_ra[i1], gaia_ra_error[i1], gaia_dec[i1], gaia_dec_error[i1], gaia_parallax[i1], gaia_parallax_error[i1], gaia_pm[i1], gaia_pm_ra[i1],\
                    gaia_pm_ra_error[i1], gaia_pm_dec[i1], gaia_pm_dec_error[i1], gaia_ruwe[i1], objinfoflag[i1], qualityflag[i1], ndetections[i1], nstackdetections[i1],\
                    ginfoflag[i1], ginfoflag2[i1], ginfoflag3[i1], rinfoflag[i1], rinfoflag2[i1], rinfoflag3[i1], iinfoflag[i1], iinfoflag2[i1], iinfoflag3[i1],\
                    zinfoflag[i1], zinfoflag2[i1], zinfoflag3[i1], yinfoflag[i1], yinfoflag2[i1], yinfoflag3[i1], diff_j, diff_h, diff_k,\
                    ukidss_j[index_ukidss_position], e_ukidss_j[index_ukidss_position], ukidss_h[index_ukidss_position], e_ukidss_h[index_ukidss_position],\
                    ukidss_k[index_ukidss_position], e_ukidss_k[index_ukidss_position]
                    writer.writerow(validate_params)
                    ob_j = np.append(ob_j, ukidss_j[index_ukidss_position])
                    e_ob_j = np.append(e_ob_j, e_ukidss_j[index_ukidss_position])
                    ob_h = np.append(ob_h, ukidss_h[index_ukidss_position])
                    e_ob_h = np.append(e_ob_h, e_ukidss_h[index_ukidss_position])
                    ob_k = np.append(ob_k, ukidss_k[index_ukidss_position])
                    e_ob_k = np.append(e_ob_k, e_ukidss_k[index_ukidss_position])
                    diff_jf = np.append(diff_jf, diff_j)
                    diff_hf = np.append(diff_hf, diff_h)
                    diff_kf = np.append(diff_kf, diff_k)
        
        plt.clf()
        plt.scatter(ob_j, e_ob_j, s=5, alpha = 0.5)
        plt.grid()
        plt.ylabel('Error in $J_{UKIDSS}$')
        plt.xlabel('$J_{UKIDSS}')
        plt.savefig('obj_vs_err_obj.png')
        plt.clf()

        plt.clf()
        plt.scatter(ob_h, e_ob_h, s=5, alpha = 0.5)
        plt.grid()
        plt.ylabel('Error in $H_{UKIDSS}$')
        plt.xlabel('$H_{UKIDSS}')
        plt.savefig('obh_vs_err_obh.png')
        plt.clf()

        plt.clf()
        plt.scatter(ob_k, e_ob_k, s=5, alpha = 0.5)
        plt.grid()
        plt.ylabel('Error in $K_{UKIDSS}$')
        plt.xlabel('$K_{UKIDSS}')
        plt.savefig('obk_vs_err_obk.png')
        plt.clf()


        bins2 = np.arange(np.min(diff_jf), np.max(diff_jf)+.1, 0.1)
        plt.clf()
        fig = plt.figure(figsize=(10,10))
        from matplotlib.gridspec import GridSpec
        gs = GridSpec(4,4)
        ax_joint = fig.add_subplot(gs[1:4,0:3])
        ax_marg_x = fig.add_subplot(gs[0,0:3])
        ax_marg_y = fig.add_subplot(gs[1:4,3])
        ax_joint.scatter(ob_j, diff_jf, alpha = 0.3, s = 5, color = 'g', label =\
                        'No. of stars =' + str(len(ob_j)))
        ax_joint.grid()
        ax_joint.set_ylim(-2,2)
        ax_joint.legend(loc = 'best')
        nx, bx, px = ax_marg_x.hist(ob_j, color = 'm', edgecolor = 'g', density =True,\
                                    alpha = 0.5, label = 'Observed J')
        ny, by, px = ax_marg_y.hist(diff_jf, bins = bins2, orientation="horizontal",\
                                    edgecolor = 'g', density=True, alpha = 0.5,\
                                    facecolor = 'orange', label = 'Difference')
        biny_max = find_nearest(ny, np.median(ny))
        ax_marg_y.set_title('Median at:%0.2f'%(by[np.where(ny==biny_max)[0][0]]))
        ax_marg_y.set_ylim(-2,2)
        ax_marg_x.grid()
        ax_marg_x.legend(loc='best')
        ax_marg_y.grid()
        ax_marg_y.legend(loc='best')
        # Turn off tick labels on marginals
        plt.setp(ax_marg_x.get_xticklabels(), visible=False)
        plt.setp(ax_marg_y.get_yticklabels(), visible=False)
        # Set labels on joint
        ax_joint.set_xlabel('Observed J magnitude')
        ax_joint.set_ylabel('$J_{UKIDSS}$ - $J_{Computed}$')
        # Set labels on marginals
        ax_marg_y.set_xlabel('N')
        ax_marg_x.set_ylabel('N')
        plt.savefig('validation_plot_j' + '_' + 'RA' + '_' + str(self.ra) + '_' + 'DEC' + str(self.dec)+'.png')
        plt.clf()

        bins2 = np.arange(diff_hf.min(), diff_hf.max()+.1, 0.1)
        plt.clf()
        fig = plt.figure(figsize=(10,10))
        from matplotlib.gridspec import GridSpec
        gs = GridSpec(4,4)
        ax_joint = fig.add_subplot(gs[1:4,0:3])
        ax_marg_x = fig.add_subplot(gs[0,0:3])
        ax_marg_y = fig.add_subplot(gs[1:4,3])
        ax_joint.scatter(ob_h, diff_hf, alpha = 0.3, s = 5, color = 'g',\
                            label = 'No. of stars =' + str(len(ob_j)))
        ax_joint.grid()
        ax_joint.set_ylim(-2,2)
        ax_joint.legend(fontsize=18, loc = 'best')
        nx, bx, px = ax_marg_x.hist(ob_h, color = 'm', edgecolor = 'g', alpha = 0.5,\
                                label = 'Observed J')
        ny, by, px = ax_marg_y.hist(diff_hf, bins = bins2, orientation="horizontal",\
                                 edgecolor = 'g', alpha = 0.5, facecolor = 'orange', label = 'Difference')
        biny_max = find_nearest(ny, np.median(ny))
        ax_marg_y.set_title('Median at: ' +str(by[np.where(ny==biny_max)[0]][0]))
        ax_marg_y.set_ylim(-2,2)
        ax_marg_x.grid()
        ax_marg_x.legend(loc='best')
        ax_marg_y.grid()
        ax_marg_y.legend(loc='best')
        # Turn off tick labels on marginals
        plt.setp(ax_marg_x.get_xticklabels(), visible=False)
        plt.setp(ax_marg_y.get_yticklabels(), visible=False)
        # Set labels on joint
        ax_joint.set_xlabel('Observed H magnitude')
        ax_joint.set_ylabel('$H_{UKIDSS}$ - $H_{Computed}$')
        # Set labels on marginals
        ax_marg_y.set_xlabel('N')
        ax_marg_x.set_ylabel('N')
        plt.savefig('validation_plot_h' + '_' + 'RA' + '_' + str(self.ra) + '_' + 'DEC' + str(self.dec)+'.png')
        plt.clf()

        bins2 = np.arange(np.min(diff_jf), np.max(diff_jf)+.1, 0.1)
        plt.clf()
        fig = plt.figure(figsize=(10,10))
        from matplotlib.gridspec import GridSpec
        gs = GridSpec(4,4)
        ax_joint = fig.add_subplot(gs[1:4,0:3])
        ax_marg_x = fig.add_subplot(gs[0,0:3])
        ax_marg_y = fig.add_subplot(gs[1:4,3])
        ax_joint.scatter(ob_k, diff_kf, alpha = 0.3, s = 5, color = 'g',\
                             label = 'No. of stars =' + str(len(ob_j)))
        ax_joint.grid()
        ax_joint.set_ylim(-2,2)
        ax_joint.legend(fontsize=18, loc = 'best')
        nx, bx, px = ax_marg_x.hist(ob_k, color = 'm', edgecolor = 'g', alpha = 0.5,\
                                        label = 'Observed J')
        ny, by, px = ax_marg_y.hist(diff_kf, bins = bins2, orientation="horizontal",\
                                    edgecolor = 'g', alpha = 0.5, facecolor = 'orange', label = 'Difference')
        biny_max = find_nearest(ny, np.median(ny))
        ax_marg_y.set_title('Median at T= ' + str(by[np.where(ny==biny_max)[0]]))
        ax_marg_y.set_ylim(-2,2)
        ax_marg_x.grid()
        ax_marg_x.legend(loc='best')
        ax_marg_y.grid()
        ax_marg_y.legend(loc='best')
        # Turn off tick labels on marginals
        plt.setp(ax_marg_x.get_xticklabels(), visible=False)
        plt.setp(ax_marg_y.get_yticklabels(), visible=False)
        # Set labels on joint
        ax_joint.set_xlabel('Observed K magnitude')
        ax_joint.set_ylabel('$K_{UKIDSS}$ - $K_{Computed}$')
        # Set labels on marginals
        ax_marg_y.set_xlabel('N')
        ax_marg_x.set_ylabel('N')
        plt.savefig('validation_plot_k' + '_' + 'RA' + '_' + str(self.ra) + '_' + 'DEC' + str(self.dec)+'.png')
        plt.clf()

