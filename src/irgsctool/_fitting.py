#pylint: disable=wrong-import-position
#pylint: disable=import-error
from datetime import date
import csv
import numpy as np
from ._sam import Models
from ._extinction_correction import ExtinctionCorrection as EC
from ._read_data import ReadData
current_datetime = date.today()


header = ['ps1_objid','ps1_ra','ps1_ra_error','ps1_dec','ps1_dec_error',\
'ps1_gpsf','ps1_gpsf_error','ps1_rpsf','ps1_rpsf_error','ps1_ipsf',\
'ps1_ipsf_error','ps1_zpsf','ps1_zpsf_error','ps1_ypsf','ps1_ypsf_error',\
'teff','logg','feh','sam_g','sam_r','sam_i','sam_z','sam_y','sam_j','sam_h',\
'sam_k','scale_factor','scale_factor_error','chi2','computed_j',\
'computed_j_error','computed_h','computed_h_error', 'computed_k',\
'computed_k_error','gaia_source_id','gaia_ra','gaia_ra_error','gaia_dec',\
'gaia_dec_error','gaia_parallax','gaia_parallax_error','gaia_pm','gaia_pm_ra',\
'gaia_pm_ra_error','gaia_pm_dec','gaia_pm_dec_error','gaia_ruwe','objinfoflag',\
'qualityflag','ndetections','nstackdetections','ginfoflag','ginfoflag2',\
'ginfoflag3','rinfoflag','rinfoflag2','rinfoflag3','iinfoflag','iinfoflag2',\
'iinfoflag3','zinfoflag','zinfoflag2','zinfoflag3','yinfoflag','yinfoflag2',\
'yinfoflag3', 'SAM Flag']

def stdv(sfavg, v1, v2, v3, v4, v5):
    mu = sfavg
    n=5
    sig = (((mu - v1)**2 + (mu - v2)**2 + (mu-v3)**2 + (mu - v4)**2 + (mu - v5)**2)/n)**0.5
    return sig

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]

def calc_sf(j, om, e_om, sm, index_min_ang_seperation, aj, ah, ak):
        ec_gmag, ec_rmag, ec_imag, ec_zmag, ec_ymag = om
        e_ec_gmag, e_ec_rmag, e_ec_imag, e_ec_zmag, e_ec_ymag = e_om
        sam_g, sam_r, sam_i, sam_z, sam_y, sam_j, sam_h, sam_k = sm

        sf_mean = (1/5.0)*((ec_gmag[j] - sam_g[index_min_ang_seperation])\
                           +(ec_rmag[j] - sam_r[index_min_ang_seperation])\
                                +(ec_imag[j] - sam_i[index_min_ang_seperation])\
                                    +(ec_zmag[j] - sam_z[index_min_ang_seperation])\
                                        + (ec_ymag[j] - sam_y[index_min_ang_seperation]))

        e_sf_mean = (1/5)*np.sqrt(e_ec_gmag[j]**2 + e_ec_rmag[j]**2 + e_ec_imag[j]**2\
                                + e_ec_zmag[j]**2  + e_ec_ymag[j]**2 )

        #0.91 is the conversion constant from J_AB to J_Vega
        #1.39 is the conversion constant from H_AB to H_Vega
        #1.85 is the conversion constant from K_AB to K_Vega

        cj = sf_mean + aj + sam_j[index_min_ang_seperation] - 0.91
        ch = sf_mean + ah + sam_h[index_min_ang_seperation] - 1.39
        ck = sf_mean + ak + sam_k[index_min_ang_seperation] - 1.85

        e_cj = np.sqrt(e_sf_mean**2)# + (self.e_aj)**2)
        e_ch = np.sqrt(e_sf_mean**2)# + (self.e_ah)**2)
        e_ck = np.sqrt(e_sf_mean**2)# + (self.e_ak)**2)

        return sf_mean, e_sf_mean, cj, e_cj, ch, e_ch, ck, e_ck


def compute_dquad(j, oc, mc):
        j = np.int64(j)
        obs_gr, obs_gi, obs_gz, obs_gy, obs_ri, obs_ry, obs_rz, obs_iz, obs_iy, obs_zy = oc
        sam_gr, sam_gi, sam_gz, sam_gy, sam_ri, sam_rz, sam_ry, sam_iz, sam_iy, sam_zy = mc
        dev_gr = obs_gr[j] - sam_gr
        dev_gi = obs_gi[j] - sam_gi
        dev_gz = obs_gz[j] - sam_gz
        dev_gy = obs_gy[j] - sam_gy
        dev_ri = obs_ri[j] - sam_ri
        dev_rz = obs_rz[j] - sam_rz
        dev_ry = obs_ry[j] - sam_ry
        dev_iz = obs_iz[j] - sam_iz
        dev_iy = obs_iy[j] - sam_iy
        dev_zy = obs_zy[j] - sam_zy
        dquad = dev_gr**2 + dev_gi**2 + dev_gz**2 + dev_gy**2 + dev_ri**2\
                 + dev_rz**2 + dev_ry**2 + dev_iz**2 + dev_iy**2 + dev_zy**2
        return dquad, np.min(dquad), dev_gr, dev_gi, dev_gz, dev_gy, dev_ri,\
            dev_rz, dev_ry, dev_iz, dev_iy, dev_zy


class GenerateIRGSC():
    """
        <justify> The *** GenerateIRGSC class *** hosts method to generate a catalog of
        probable stellar sources in the PANSTARRS data with their computed magnitudes, astrometric
        information from GAIA DR3 data, best fitted model parameters and flags.</justify>
    
    """

    def __init__(self, ra, dec):
         self.ra, self.dec = ra, dec
         self.rd = ReadData(ra,dec)
         self.ec = EC(ra, dec)

    def generate_irgsc(self, use_optimal_method=True):
        """
            `irgsctool.GenerateIRGSC.generate_irgsc(use_optimal_method=True)`
            
            <justify> This method finds 
            the best fitting model to the observed colors of the stellar source.
            The best fitting model is chosen from a combination of Kurucz/Castelli-Kurucz 
            and Phoenix synthetic spectra convolved with the PANSTARRS response function
            (or BANDPASS) which is integrated w.r.t. the wavelength and normalised to
            the product of the PANSTARRS response function and wavelength. This is also
            called as Effective Stimulus (ES). </justify>

            $$
            ES = \\frac{\\int{F_{\\lambda}P_{\\lambda}{\\lambda}
            d{\\lambda}}}{\\int{P_{\\lambda}{\\lambda}d{\\lambda}}}
            $$
            
            The spectra is obtained from pysynphot 
            ([More information here](https://pysynphot.readthedocs.io)).

            Returns:
                    irgsc_data: The generated IRGSC.

        """
        if use_optimal_method is True:
            print("")
            print('#########################################')
            print('Computing the NIR magnitudes for the sources using the optimal method')
            print('#########################################')
            print("")

            ps1_objid, ps_ra, err_ps_ra, ps_dec, err_ps_dec, ec_gmag, e_ec_gmag, gkron, e_gkron,\
            ec_rmag, e_ec_rmag, rkron, e_rkron, ec_imag, e_ec_imag, ikron, e_ikron, ec_zmag,\
            e_ec_zmag, zkron, e_zkron, ec_ymag, e_ec_ymag, ykron, e_ykron, objinfoflag, qualityflag,\
            ndetections, nstackdetections, ginfoflag, ginfoflag2, ginfoflag3, rinfoflag, rinfoflag2,\
            rinfoflag3, iinfoflag, iinfoflag2, iinfoflag3, zinfoflag, zinfoflag2, zinfoflag3,\
            yinfoflag, yinfoflag2, yinfoflag3 = self.ec.extinction_corrected_photometry()
            _,_,aj,ah,ak = self.ec.get_reddening()
            gaia_data = self.rd.read_gaia_data()
            gaia_source_id, gaia_ra, gaia_ra_error, gaia_dec, gaia_dec_error, gaia_parallax,\
            gaia_parallax_error, gaia_pm, gaia_pm_ra, gaia_pm_ra_error, gaia_pm_dec,\
            gaia_pm_dec_error, gaia_ruwe = gaia_data
            
            k0 = Models('Kurucz')
            k0.read_sam_file()
            
            c1 = Models('Phoenix')
            c1.read_sam_file()
            
            c2 = Models('Phoenix')
            c2.read_sam_file()

            model_params_k0 = k0.select_sam_range(teff_range=[4000,10000],
                                                logg_range=None, feh_range=None)
            model_params_c1 = c1.select_sam_range(teff_range=[2800,5000],
                                                logg_range=[3.0,5.5],
                                                feh_range=[-5.0,-1.5])
            model_params_c2 = c2.select_sam_range(teff_range=[2800,4000],
                                                logg_range=[0.0,3.0],
                                                feh_range=[-0.5,1.5])
            
            
            teff_c1, logg_c1, feh_c1, sam_g_c1, sam_r_c1, sam_i_c1, sam_z_c1, sam_y_c1, sam_j_c1,\
            sam_h_c1, sam_k_c1 = model_params_c1
            teff_c2, logg_c2, feh_c2, sam_g_c2, sam_r_c2, sam_i_c2, sam_z_c2, sam_y_c2, sam_j_c2,\
            sam_h_c2, sam_k_c2 = model_params_c2
            teff_k0, logg_k0, feh_k0, sam_g_k0, sam_r_k0, sam_i_k0, sam_z_k0, sam_y_k0, sam_j_k0,\
            sam_h_k0, sam_k_k0 = model_params_k0
            
            len_c1 = len(teff_c1)
            len_c2 = len(teff_c2)
            len_k0 = len(teff_k0)

            teff = np.concatenate((teff_c1, teff_c2, teff_k0), axis=0)
            logg = np.concatenate((logg_c1, logg_c2, logg_k0), axis=0)
            feh = np.concatenate((feh_c1, feh_c2, feh_k0), axis=0)
            sam_g = np.concatenate((sam_g_c1, sam_g_c2, sam_g_k0), axis=0)
            sam_r = np.concatenate((sam_r_c1, sam_r_c2, sam_r_k0), axis=0)
            sam_i = np.concatenate((sam_i_c1, sam_i_c2, sam_i_k0), axis=0)
            sam_z = np.concatenate((sam_z_c1, sam_z_c2, sam_z_k0), axis=0)
            sam_y = np.concatenate((sam_y_c1, sam_y_c2, sam_y_k0), axis=0)
            sam_j = np.concatenate((sam_j_c1, sam_j_c2, sam_j_k0), axis=0)
            sam_h = np.concatenate((sam_h_c1, sam_h_c2, sam_h_k0), axis=0)
            sam_k = np.concatenate((sam_k_c1, sam_k_c2, sam_k_k0), axis=0)

            sam_gr = sam_g - sam_r
            sam_ri = sam_r - sam_i
            sam_gi = sam_g - sam_i
            sam_gz = sam_g - sam_z
            sam_gy = sam_g - sam_y
            sam_ry = sam_r - sam_y
            sam_rz = sam_r - sam_z
            sam_iz = sam_i - sam_z
            sam_iy = sam_i - sam_y
            sam_zy = sam_z - sam_y

            observed_colours = (ec_gmag - ec_rmag), (ec_gmag - ec_imag),\
            (ec_gmag - ec_rmag), (ec_gmag - ec_ymag), (ec_rmag - ec_imag),\
            (ec_rmag - ec_ymag), (ec_rmag - ec_zmag), (ec_imag - ec_zmag),\
            (ec_imag - ec_ymag), (ec_zmag - ec_ymag)

            model_colours = sam_gr, sam_gi, sam_gz, sam_gy, sam_ri, sam_rz,\
            sam_ry, sam_iz, sam_iy, sam_zy

            observed_optical_magnitudes = ec_gmag, ec_rmag, ec_imag,\
                ec_zmag, ec_ymag
            e_observed_optical_magnitudes = e_ec_gmag, e_ec_rmag,\
                e_ec_imag, e_ec_zmag, e_ec_ymag
            sam_magnitudes = sam_g, sam_r, sam_i, sam_z, sam_y,\
                sam_j, sam_h, sam_k

            data = []
            ra_name=str(self.ra).replace('.','_');dec_name=str(self.dec)\
                .replace('.', '_')

            with open('IRGSC'+'_'+'RA'+str(ra_name)+'DEC'+str(dec_name)+\
                    str(current_datetime)+'.csv','w',encoding='UTF8') as file1:
                writer=csv.writer(file1)
                writer.writerow(header)
                for j in range(len(ec_gmag)):
                    dquad_arr, min_dquad, _, _, _, _, _, _, _, _, _, _ = \
                        compute_dquad(j, oc = observed_colours, mc = model_colours)
                    min_dquad_element = find_nearest(dquad_arr,min_dquad)
                    index_best_fit_sam = np.where(min_dquad_element==(dquad_arr))[0]
                    if index_best_fit_sam<=len_c1:
                        sam_model = 'c1'
                    elif index_best_fit_sam>len_c1 and index_best_fit_sam <=len_c2:
                        sam_model = 'c2'
                    elif index_best_fit_sam > len_c2:
                        sam_model = 'K0'
                        
                    sf_avg,sigma_sf,computed_j,computed_j_error,computed_h,computed_h_error,\
                        computed_k, computed_k_error=calc_sf(j, observed_optical_magnitudes,\
                                                            e_observed_optical_magnitudes,\
                                                            sam_magnitudes, index_best_fit_sam,\
                                                            aj, ah, ak)
                    gaia_angular_seperation = 3600*np.sqrt(((ps_ra[j]
                                                            -gaia_ra)*np.cos(np.radians(ps_dec[j])))**2
                                                            +(ps_dec[j] - gaia_dec)**2)
                    index_min_ang_seperation = np.where(gaia_angular_seperation<=1.0)[0]
                    if len(index_min_ang_seperation) > 1.0:
                            gaia_ang_seperation_selected = gaia_angular_seperation[index_min_ang_seperation]
                            min_gaia_ang_seperation = gaia_angular_seperation\
                                [np.where(np.min(gaia_ang_seperation_selected)\
                                          == gaia_angular_seperation)[0]]

                            index_min_ang_seperation = np.where(min_gaia_ang_seperation \
                                                                == gaia_angular_seperation)[0]

                            irgsc_data = ps1_objid[j], ps_ra[j], err_ps_ra[j], ps_dec[j], \
                                err_ps_dec[j], ec_gmag[j], e_ec_gmag[j], ec_rmag[j], \
                                e_ec_rmag[j], ec_imag[j], e_ec_imag[j], ec_zmag[j], \
                                e_ec_zmag[j], ec_ymag[j], e_ec_ymag[j], teff[index_best_fit_sam][0], \
                                logg[index_best_fit_sam][0], feh[index_best_fit_sam][0], \
                                sam_g[index_best_fit_sam][0], sam_r[index_best_fit_sam][0], sam_i[index_best_fit_sam][0], \
                                sam_z[index_best_fit_sam][0], sam_y[index_best_fit_sam][0], sam_j[index_best_fit_sam][0],\
                                sam_h[index_best_fit_sam][0], sam_k[index_best_fit_sam][0], sf_avg[0], sigma_sf,\
                                min_dquad_element, computed_j[0], computed_j_error, computed_h[0],\
                                computed_h_error, computed_k[0], computed_k_error,\
                                gaia_source_id[index_min_ang_seperation][0], gaia_ra[index_min_ang_seperation][0],\
                                gaia_ra_error[index_min_ang_seperation][0], gaia_dec[index_min_ang_seperation][0],\
                                gaia_dec_error[index_min_ang_seperation][0], gaia_parallax[index_min_ang_seperation][0],\
                                gaia_parallax_error[index_min_ang_seperation][0], gaia_pm[index_min_ang_seperation][0],\
                                gaia_pm_ra[index_min_ang_seperation][0], gaia_pm_ra_error[index_min_ang_seperation][0],\
                                gaia_pm_dec[index_min_ang_seperation][0], gaia_pm_dec_error[index_min_ang_seperation][0],\
                                gaia_ruwe[index_min_ang_seperation][0], objinfoflag[j], qualityflag[j], ndetections[j],\
                                nstackdetections[j], ginfoflag[j], ginfoflag2[j], ginfoflag3[j], rinfoflag[j], rinfoflag2[j],\
                                rinfoflag3[j], iinfoflag[j], iinfoflag2[j], iinfoflag3[j], zinfoflag[j], zinfoflag2[j],\
                                zinfoflag3[j], yinfoflag[j], yinfoflag2[j], yinfoflag3[j],sam_model
                            writer.writerow(irgsc_data)
                    elif len(index_min_ang_seperation) == 0.0:
                            irgsc_data = ps1_objid[j], ps_ra[j], err_ps_ra[j], ps_dec[j], err_ps_dec[j],\
                                ec_gmag[j], e_ec_gmag[j], ec_rmag[j], e_ec_rmag[j], ec_imag[j],\
                                e_ec_imag[j], ec_zmag[j], e_ec_zmag[j], ec_ymag[j], e_ec_ymag[j],\
                                teff[index_best_fit_sam][0], \
                                logg[index_best_fit_sam][0], feh[index_best_fit_sam][0], \
                                sam_g[index_best_fit_sam][0], sam_r[index_best_fit_sam][0], sam_i[index_best_fit_sam][0], \
                                sam_z[index_best_fit_sam][0], sam_y[index_best_fit_sam][0], sam_j[index_best_fit_sam][0],\
                                sam_h[index_best_fit_sam][0], sam_k[index_best_fit_sam][0], sf_avg[0], sigma_sf,\
                                min_dquad_element, computed_j[0], computed_j_error, computed_h[0],\
                                computed_h_error, computed_k[0], computed_k_error, -999, -999, -999,\
                                -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, objinfoflag[j],\
                                qualityflag[j], ndetections[j], nstackdetections[j], ginfoflag[j], ginfoflag2[j],\
                                ginfoflag3[j], rinfoflag[j], rinfoflag2[j], rinfoflag3[j], iinfoflag[j],\
                                iinfoflag2[j], iinfoflag3[j], zinfoflag[j], zinfoflag2[j], zinfoflag3[j],\
                                yinfoflag[j], yinfoflag2[j], yinfoflag3[j], sam_model
                            writer.writerow(irgsc_data)
                    elif len(index_min_ang_seperation) == 1.0:
                            irgsc_data = ps1_objid[j], ps_ra[j], err_ps_ra[j], ps_dec[j], err_ps_dec[j],\
                                ec_gmag[j], e_ec_gmag[j], ec_rmag[j], e_ec_rmag[j], ec_imag[j],\
                                e_ec_imag[j], ec_zmag[j], e_ec_zmag[j], ec_ymag[j], e_ec_ymag[j],\
                                teff[index_best_fit_sam][0], \
                                logg[index_best_fit_sam][0], feh[index_best_fit_sam][0], \
                                sam_g[index_best_fit_sam][0], sam_r[index_best_fit_sam][0], sam_i[index_best_fit_sam][0], \
                                sam_z[index_best_fit_sam][0], sam_y[index_best_fit_sam][0], sam_j[index_best_fit_sam][0],\
                                sam_h[index_best_fit_sam][0], sam_k[index_best_fit_sam][0], sf_avg[0], sigma_sf,\
                                min_dquad_element, computed_j[0], computed_j_error, computed_h[0],\
                                computed_h_error, computed_k[0], computed_k_error, gaia_source_id[index_min_ang_seperation][0],\
                                gaia_ra[index_min_ang_seperation][0], gaia_ra_error[index_min_ang_seperation][0],\
                                gaia_dec[index_min_ang_seperation][0], gaia_dec_error[index_min_ang_seperation][0],\
                                gaia_parallax[index_min_ang_seperation][0], gaia_parallax_error[index_min_ang_seperation][0],\
                                gaia_pm[index_min_ang_seperation][0], gaia_pm_ra[index_min_ang_seperation][0],\
                                gaia_pm_ra_error[index_min_ang_seperation][0], gaia_pm_dec[index_min_ang_seperation][0],\
                                gaia_pm_dec_error[index_min_ang_seperation][0], gaia_ruwe[index_min_ang_seperation][0],\
                                objinfoflag[j], qualityflag[j], ndetections[j], nstackdetections[j], ginfoflag[j],\
                                ginfoflag2[j], ginfoflag3[j], rinfoflag[j], rinfoflag2[j], rinfoflag3[j], iinfoflag[j],\
                                iinfoflag2[j], iinfoflag3[j], zinfoflag[j], zinfoflag2[j], zinfoflag3[j], yinfoflag[j],\
                                iinfoflag2[j], yinfoflag3[j], sam_model
                            writer.writerow(irgsc_data)
        return irgsc_data
