import os
import sys
from datetime import date
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import pylab
from ._get_data import GetData
current_datetime = date.today()
home_dir = os.getcwd()

params = {'legend.fontsize': 'x-large',
          'figure.figsize': (10,10),
         'axes.labelsize': 'x-large',
         'axes.titlesize':'x-large',
         'xtick.labelsize':'x-large',
         'ytick.labelsize':'x-large'}
pylab.rcParams.update(params)

class ReadData():
    """
    *** ReadData *** class contains methods to read the photometric data from PS1 DR2,
    GAIA DR3 and UKIDSS DR11.
    """
    def __init__(self, ra, dec):
         self.ra, self.dec = ra, dec
         self.gd = GetData(ra,dec)

    def read_optical_data(self):
        """
            `irgsctool.ReadData.read_optical_data()`
            ------------------------------------------

            <justify> 
            This function reads the input optical PANSTARRS data.
            The number of columns in the input file are 32.
            After reading the input data, this function filters it for
            nan values (if present) and restricts the data to the sources having
            detection in all the five bands and that have SNR atleast 5.
            This data is then fed to the Star-Galaxy classification routine to
            seperate stars and galaxies in the data. <\justify>

            Raises:
                FileNotFoundError: if the optical input data file is not available.
        """
        ra_name = str(self.ra).replace('.','_')
        dec_name = str(self.dec).replace('.', '_')
        file_name = 'PS1' + '_' + 'RA'+str(ra_name) + 'DEC' + str(dec_name)
        try:
            ps1_data = np.genfromtxt(str(file_name) +'.csv',\
                                     delimiter=',', skip_header=1)
 
        except FileNotFoundError:
            self.gd.get_panstarrs_data()
            ps1_data = np.genfromtxt(str(file_name) +'.csv',\
                                     delimiter=',', skip_header=1)
        ps1_objid = ps1_data[:,0]
        ps_ra = ps1_data[:,1]
        err_ps_ra = ps1_data[:,2]
        ps_dec = ps1_data[:,3]
        err_ps_dec = ps1_data[:,4]
        gmag = ps1_data[:,5]
        e_gmag = ps1_data[:,6]
        gkron = ps1_data[:,7]
        e_gkron = ps1_data[:,8]
        rmag = ps1_data[:,9]
        e_rmag = ps1_data[:,10]
        rkron = ps1_data[:,11]
        e_rkron = ps1_data[:,12]
        imag = ps1_data[:,13]
        e_imag = ps1_data[:,14]
        ikron = ps1_data[:,15]
        e_ikron = ps1_data[:,16]
        zmag = ps1_data[:,17]
        e_zmag = ps1_data[:,18]
        zkron = ps1_data[:,19]
        e_zkron = ps1_data[:,20]
        ymag = ps1_data[:,21]
        e_ymag = ps1_data[:,22]
        ykron = ps1_data[:,23]
        e_ykron = ps1_data[:,24]
        objinfoflag = ps1_data[:,25]
        qualityflag = ps1_data[:,26]
        ndetections = ps1_data[:,27]
        nstackdetections = ps1_data[:,28]
        ginfoflag = ps1_data[:,29]
        ginfoflag2 = ps1_data[:,30]
        ginfoflag3 = ps1_data[:,31]
        rinfoflag = ps1_data[:,32]
        rinfoflag2 = ps1_data[:,33]
        rinfoflag3 = ps1_data[:,34]
        iinfoflag = ps1_data[:,35]
        iinfoflag2 = ps1_data[:,36]
        iinfoflag3 = ps1_data[:,37]
        zinfoflag = ps1_data[:,38]
        zinfoflag2 = ps1_data[:,39]
        zinfoflag3 = ps1_data[:,40]
        yinfoflag = ps1_data[:,41]
        yinfoflag2 = ps1_data[:,42]
        yinfoflag3 = ps1_data[:,43]
            
        ps1_data = ps1_objid, ps_ra, err_ps_ra, ps_dec, err_ps_dec, gmag,\
            e_gmag, gkron, e_gkron, rmag, e_rmag, rkron, e_rkron, imag,\
            e_imag, ikron, e_ikron, zmag, e_zmag, zkron, e_zkron, ymag,\
            e_ymag, ykron, e_ykron, objinfoflag, qualityflag, ndetections,\
            nstackdetections, ginfoflag, ginfoflag2, ginfoflag3, rinfoflag,\
            rinfoflag2, rinfoflag3, iinfoflag, iinfoflag2, iinfoflag3,\
            zinfoflag, zinfoflag2, zinfoflag3, yinfoflag, yinfoflag2,\
            yinfoflag3

        print('Number of rows in the PANSTARRS file:', len(ps1_objid))

        oid1 = np.array([list(set(ps1_objid))])[0]
        oid1 = [*set(ps1_objid)]
        ptsf = []
        for i in range(len(oid1)):
            ptsi = np.where(oid1[i]==ps1_objid)[0]
            if len(ptsi)>1.0:
                ptsi = ptsi[0]
                ptsf = np.append(ptsf, ptsi)
                ptsf = np.int64(ptsf)
            else:
                 ptsf = np.append(ptsf, ptsi)

        ps1_objid = ps1_objid[ptsf]
        ps_ra = ps_ra[ptsf]
        ps_dec = ps_dec[ptsf]
        err_ps_ra = err_ps_ra[ptsf]
        err_ps_dec = err_ps_dec[ptsf]
        gmag = gmag[ptsf]
        gkron = gkron[ptsf]
        e_gmag = e_gmag[ptsf]
        e_gkron = e_gkron[ptsf]
        rmag = rmag[ptsf]
        rkron = rkron[ptsf]
        e_rmag = e_rmag[ptsf]
        e_rkron = e_rkron[ptsf]
        imag = imag[ptsf]
        ikron = ikron[ptsf]
        e_imag = e_imag[ptsf]
        e_ikron = e_ikron[ptsf]
        zmag = zmag[ptsf]
        zkron = zkron[ptsf]
        e_zmag = e_zmag[ptsf]
        e_zkron = e_zkron[ptsf]
        ymag = ymag[ptsf]
        ykron = ykron[ptsf]
        e_ymag = e_ymag[ptsf]
        e_ykron = e_ykron[ptsf]
        objinfoflag = objinfoflag[ptsf]
        qualityflag = qualityflag[ptsf]
        ndetections = ndetections[ptsf]
        nstackdetections = nstackdetections[ptsf]
        ginfoflag = ginfoflag[ptsf]
        ginfoflag2 = ginfoflag2[ptsf]
        ginfoflag3 = ginfoflag3[ptsf]
        rinfoflag = rinfoflag[ptsf]
        rinfoflag2 = rinfoflag2[ptsf]
        rinfoflag3 = rinfoflag3[ptsf]
        iinfoflag = iinfoflag[ptsf]
        iinfoflag2 = iinfoflag2[ptsf]
        iinfoflag3 = iinfoflag3[ptsf]
        zinfoflag = zinfoflag[ptsf]
        zinfoflag2 = zinfoflag2[ptsf]
        zinfoflag3 = zinfoflag3[ptsf]
        yinfoflag = yinfoflag[ptsf]
        yinfoflag2 = yinfoflag2[ptsf]
        yinfoflag3 = yinfoflag3[ptsf]
        
        print("")
        print('Now filtering the optical data for nan values')
        print("")
        print('Number of sources in the PANSTARRS data:',\
              len(ps1_objid))
        print("")

        indices_only_ifilered = np.where(imag!= -999)[0]
        binwidth=0.5
        bins = np.arange(np.min(imag[indices_only_ifilered]),\
                         np.max(imag[indices_only_ifilered]) +\
                            binwidth, binwidth)
        plt.clf()
        plt.hist(imag[indices_only_ifilered], bins=bins,\
                 facecolor='white', edgecolor='red',\
                    linestyle='--', label='$i_{mag}$')
        plt.xlabel('Only $i_{mag}$ observations')
        plt.legend(loc='best')
        plt.grid()
        plt.savefig('hist_only_imag_in_panstarrs_data.png')
        plt.clf()

        indices_all_filtered = np.where((gmag!= -999) & (imag!= -999) &\
                                        (rmag != -999) & (zmag != -999) &\
                                        (ymag != -999) & (e_gmag != -999) &\
                                        (e_rmag != -999) & (e_imag != -999) &\
                                        (e_zmag != -999) & (e_ymag != -999) &\
                                        (gkron!= -999) & (ikron!= -999) &\
                                        (zkron != -999) & (ykron != -999) &\
                                        (rkron != -999) & (e_gmag < 0.2) &\
                                        (e_rmag < 0.2) & (e_imag < 0.2) &\
                                        (e_zmag < 0.2) & (e_ymag < 0.2))[0]

        print('Number of sources having detections in five optical bands=',\
              len(indices_all_filtered))

        binwidth=0.5
        bins_g = np.arange(np.min(gmag[indices_all_filtered]),\
                           np.max(gmag[indices_all_filtered]) +\
                            binwidth, binwidth)
        bins_r = np.arange(np.min(rmag[indices_all_filtered]),\
                           np.max(rmag[indices_all_filtered]) +\
                            binwidth, binwidth)
        bins_i = np.arange(np.min(imag[indices_all_filtered]),\
                           np.max(imag[indices_all_filtered]) +\
                            binwidth, binwidth)
        bins_z = np.arange(np.min(zmag[indices_all_filtered]),\
                           np.max(zmag[indices_all_filtered]) +\
                            binwidth, binwidth)
        bins_y = np.arange(np.min(ymag[indices_all_filtered]),\
                           np.max(ymag[indices_all_filtered]) +\
                            binwidth, binwidth)
        
        plt.clf()
        plt.hist(gmag[indices_all_filtered], bins=bins_g,\
                 facecolor='white', edgecolor='red',\
                linestyle='--', label='$g_{mag}$')
        plt.hist(rmag[indices_all_filtered], bins=bins_r,\
                 facecolor='white', edgecolor='blue',\
                linestyle='--', label='$r_{mag}$')
        plt.hist(imag[indices_all_filtered], bins=bins_i,\
                 facecolor='white', edgecolor='green',\
                linestyle='--', label='$i_{mag}$')
        plt.hist(zmag[indices_all_filtered], bins=bins_z,\
                 facecolor='white', edgecolor='orange',\
                    linestyle='--', label='$z_{mag}$')
        plt.hist(ymag[indices_all_filtered], bins=bins_y,\
                 facecolor='white', edgecolor='purple',\
                    linestyle='--', label='$y_{mag}$')
        plt.xlabel('PANSTARRS observed data')
        plt.legend(loc='best')
        plt.grid()
        plt.savefig('hist_panstarrs_data.png')
        plt.clf()

        filtered_optical_data = ps1_objid[indices_all_filtered],\
            ps_ra[indices_all_filtered], err_ps_ra[indices_all_filtered],\
            ps_dec[indices_all_filtered], err_ps_dec[indices_all_filtered],\
            gmag[indices_all_filtered], e_gmag[indices_all_filtered],\
            gkron[indices_all_filtered], e_gkron[indices_all_filtered],\
            rmag[indices_all_filtered], e_rmag[indices_all_filtered],\
            rkron[indices_all_filtered], e_rkron[indices_all_filtered],\
            imag[indices_all_filtered], e_imag[indices_all_filtered],\
            ikron[indices_all_filtered], e_ikron[indices_all_filtered],\
            zmag[indices_all_filtered], e_zmag[indices_all_filtered],\
            zkron[indices_all_filtered], e_zkron[indices_all_filtered],\
            ymag[indices_all_filtered], e_ymag[indices_all_filtered],\
            ykron[indices_all_filtered], e_ykron[indices_all_filtered],\
            objinfoflag[indices_all_filtered], qualityflag[indices_all_filtered],\
            ndetections[indices_all_filtered], nstackdetections[indices_all_filtered],\
            ginfoflag[indices_all_filtered], ginfoflag2[indices_all_filtered],\
            ginfoflag3[indices_all_filtered], rinfoflag[indices_all_filtered],\
            rinfoflag2[indices_all_filtered], rinfoflag3[indices_all_filtered],\
            iinfoflag[indices_all_filtered], iinfoflag2[indices_all_filtered],\
            iinfoflag3[indices_all_filtered], zinfoflag[indices_all_filtered],\
            zinfoflag2[indices_all_filtered], zinfoflag3[indices_all_filtered],\
            yinfoflag[indices_all_filtered], yinfoflag2[indices_all_filtered],\
            yinfoflag3[indices_all_filtered]
        return filtered_optical_data

    def read_nir_data(self):
        """
        Reads the input UKIDSS NIR data. The number of columns are 8.

        Returns the input optical data with nan values removed (if present)
        and restricts the data to the sources having SNR >= 5.

        Some regions do not have J or H band data especially DXS or GCS
        surveys. For these regions, only K band data is imported.
        """

        ra_name = str(self.ra).replace('.','_')
        dec_name = str(self.dec).replace('.', '_')
        file_name = 'UKIDSS' + '_' + 'RA'+str(ra_name) + 'DEC' + str(dec_name) + '.csv'
        print('UKIDSS file name=', file_name)
        file_exists = os.path.exists(file_name)
        print("")
        print('Does UKIDSS observed NIR data file exist?', file_exists)
        print("")
        if file_exists is False:
            print('Validated catalogue does not exist')
            print("")
            print('############################################')
            print('Generating observed UKIDSS NIR data file...')
            print("")
            self.gd.get_ukidss_data()
            ukidss_data = np.genfromtxt(str(file_name), delimiter=',', skip_header=1)
            if len(ukidss_data)<9.0:
                raise ValueError('No observations in UKIDSS')
                sys.exit(0)
            petro_j = ukidss_data[:,2]; e_petro_j = ukidss_data[:,3]
            petro_h = ukidss_data[:,4]; e_petro_h = ukidss_data[:,5]
            petro_k = ukidss_data[:,6]; e_petro_k = ukidss_data[:,7]
            ukidss_ra = ukidss_data[:,0]; ukidss_dec = ukidss_data[:,1]

            nir_filter_index = np.where((petro_j != -999999488) & (e_petro_j != -999999488) &\
                    (petro_h != -999999488) & (e_petro_h != -999999488) &\
                        (petro_k != -999999488) & (e_petro_k != -999999488) &\
                            (e_petro_j < 0.2) & (e_petro_h < 0.2) & (e_petro_k < 0.2))[0]
            
            if len(nir_filter_index) == 0.0:

                    new_nir_filter_index = np.where((petro_j != -999999488)\
                                & (petro_k != -999999488) & (petro_h != -999999488) &\
                                (np.abs(e_petro_j) < 0.2) & (np.abs(e_petro_h) < 0.2)\
                                & (np.abs(e_petro_k) < 0.2))[0]
                    
                    e_petro_j = 0.005*petro_j[new_nir_filter_index]
                    e_petro_h = 0.05*petro_h[new_nir_filter_index]
                    e_petro_k = 0.05*petro_k[new_nir_filter_index]

                    print("")
                    print('Number of sources in the NIR data = ', len(new_nir_filter_index))
                    print("")
        
                    filtered_petro_j = petro_j[new_nir_filter_index]
                    filtered_petro_h = petro_h[new_nir_filter_index]
                    filtered_petro_k = petro_k[new_nir_filter_index]
                    filtered_ukidss_ra = ukidss_ra[new_nir_filter_index]
                    filtered_ukidss_dec = ukidss_dec[new_nir_filter_index]

                    binwidth=0.5
                    binsj = np.arange(np.min(filtered_petro_j),\
                                      np.max(filtered_petro_j)\
                                        + binwidth, binwidth)
                    binsh = np.arange(np.min(filtered_petro_h),\
                                      np.max(filtered_petro_h)\
                                      + binwidth, binwidth)
                    binsk = np.arange(np.min(filtered_petro_k),\
                                      np.max(filtered_petro_k)\
                                      + binwidth, binwidth)

                    plt.clf()
                    plt.hist(filtered_petro_j, bins = binsj,\
                             facecolor='white', edgecolor = 'red',\
                            linestyle='--', label='Observed J')
                    plt.hist(filtered_petro_h, bins = binsh,\
                             facecolor='white', edgecolor = 'blue',\
                            linestyle='--', label='Observed H')
                    plt.hist(filtered_petro_k, bins = binsk,\
                             facecolor='white',edgecolor = 'green',\
                            linestyle='--', label='Observed K')
                    plt.xlabel('petro magnitudes')
                    plt.ylabel('Bins')
                    plt.gird()
                    plt.legend(loc='best')
                    plt.savefig('hist_ukidss_nir.png')
                    plt.clf()

                    nir_data =  filtered_petro_j, filtered_petro_h, filtered_petro_k,\
                        e_petro_j, e_petro_h, e_petro_k, filtered_ukidss_ra, filtered_ukidss_dec
            else:
                    print("")
                    print('Number of sources in the NIR data = ', len(nir_filter_index))
                    print("")
        
                    filtered_petro_j = petro_j[nir_filter_index]
                    e_petro_h = e_petro_h[nir_filter_index]
                    filtered_petro_h = petro_h[nir_filter_index]
                    e_petro_k = e_petro_k[nir_filter_index]
                    filtered_petro_k = petro_k[nir_filter_index]
                    e_petro_j = e_petro_j[nir_filter_index]
                    filtered_ukidss_ra = ukidss_ra[nir_filter_index]
                    filtered_ukidss_dec = ukidss_dec[nir_filter_index]

                    binwidth=0.5
                    binsj = np.arange(np.min(filtered_petro_j),\
                                      np.max(filtered_petro_j)\
                                        + binwidth, binwidth)
                    binsh = np.arange(np.min(filtered_petro_h),\
                                      np.max(filtered_petro_h)\
                                      + binwidth, binwidth)
                    binsk = np.arange(np.min(filtered_petro_k),\
                                      np.max(filtered_petro_k)\
                                      + binwidth, binwidth)

                    plt.clf()
                    plt.hist(filtered_petro_j, bins = binsj,\
                             facecolor='white', edgecolor = 'red',\
                            linestyle='--', label='Observed J')
                    plt.hist(filtered_petro_h, bins = binsh,\
                             facecolor='white', edgecolor = 'blue',\
                            linestyle='--', label='Observed H')
                    plt.hist(filtered_petro_k, bins = binsk,\
                             facecolor='white',edgecolor = 'green',\
                            linestyle='--', label='Observed K')
                    plt.xlabel('petro magnitudes')
                    plt.ylabel('Bins')
                    plt.grid()
                    plt.legend(loc='best')
                    plt.savefig('hist_ukidss_nir.png')
                    plt.clf()
                    
                    nir_data =  filtered_petro_j, filtered_petro_h, filtered_petro_k,\
                        e_petro_j, e_petro_h, e_petro_k, filtered_ukidss_ra, filtered_ukidss_dec

        elif file_exists is True:
                print("")
                print('Reading the validated catalogue file:'+ str(file_name))
                print("")
                ukidss_data = np.genfromtxt(str(file_name), delimiter=',', skip_header=1)
                if len(ukidss_data)<9.0:
                    raise ValueError('No observations in UKIDSS')
                    sys.exit(0)
                petro_j = ukidss_data[:,2]; e_petro_j = ukidss_data[:,3]
                petro_h = ukidss_data[:,4]; e_petro_h = ukidss_data[:,5]
                petro_k = ukidss_data[:,6]; e_petro_k = ukidss_data[:,7]
                ukidss_ra = ukidss_data[:,0]; ukidss_dec = ukidss_data[:,1]

                nir_filter_index = np.where((petro_j != -999999488) & (e_petro_j != -999999488) &\
                    (petro_h != -999999488) & (e_petro_h != -999999488) &\
                        (petro_k != -999999488) & (e_petro_k != -999999488) &\
                            (e_petro_j < 0.2) & (e_petro_h < 0.2) & (e_petro_k < 0.2))[0]
                if len(nir_filter_index) == 0:

                    new_nir_filter_index = np.where((petro_j != -999999488)\
                                        & (petro_k != -999999488) & (petro_h != -999999488) &\
                                            (e_petro_j < 0.2) & (e_petro_h < 0.2) & (e_petro_k < 0.2))[0]
                    e_petro_j = 0.005*petro_j[new_nir_filter_index]
                    e_petro_h = 0.05*petro_h[new_nir_filter_index]
                    e_petro_k = 0.05*petro_k[new_nir_filter_index]

                    print("")
                    print('Number of sources in the NIR data = ',\
                          len(new_nir_filter_index))
                    print("")
        
                    filtered_petro_j = petro_j[new_nir_filter_index]
                    filtered_petro_h = petro_h[new_nir_filter_index]
                    filtered_petro_k = petro_k[new_nir_filter_index]
                    filtered_ukidss_ra = ukidss_ra[new_nir_filter_index]
                    filtered_ukidss_dec = ukidss_dec[new_nir_filter_index]

                    binwidth=0.5
                    binsj = np.arange(np.min(filtered_petro_j),\
                                      np.max(filtered_petro_j)\
                                        + binwidth, binwidth)
                    binsh = np.arange(np.min(filtered_petro_h),\
                                      np.max(filtered_petro_h)\
                                      + binwidth, binwidth)
                    binsk = np.arange(np.min(filtered_petro_k),\
                                      np.max(filtered_petro_k)\
                                      + binwidth, binwidth)

                    plt.clf()
                    plt.hist(filtered_petro_j, bins = binsj,\
                             facecolor='white', edgecolor = 'red',\
                            linestyle='--', label='Observed J')
                    plt.hist(filtered_petro_h, bins = binsh,\
                             facecolor='white', edgecolor = 'blue',\
                            linestyle='--', label='Observed H')
                    plt.hist(filtered_petro_k, bins = binsk,\
                             facecolor='white',edgecolor = 'green',\
                            linestyle='--', label='Observed K')
                    plt.xlabel('petro magnitudes')
                    plt.ylabel('Bins')
                    plt.grid()
                    plt.legend(loc='best')
                    plt.savefig('hist_ukidss_nir.png')
                    plt.clf()

                    nir_data =  filtered_petro_j, filtered_petro_h, filtered_petro_k,\
                        e_petro_j, e_petro_h, e_petro_k, filtered_ukidss_ra, filtered_ukidss_dec
                else:
                    print("")
                    print('Number of sources in the NIR data = ', len(nir_filter_index))
                    print("")
        
                    filtered_petro_j = petro_j[nir_filter_index]
                    e_petro_h = e_petro_h[nir_filter_index]
                    filtered_petro_h = petro_h[nir_filter_index]
                    e_petro_k = e_petro_k[nir_filter_index]
                    filtered_petro_k = petro_k[nir_filter_index]
                    e_petro_j = e_petro_j[nir_filter_index]
                    filtered_ukidss_ra = ukidss_ra[nir_filter_index]
                    filtered_ukidss_dec = ukidss_dec[nir_filter_index]

                    binwidth=0.5
                    binsj = np.arange(np.min(filtered_petro_j), np.max(filtered_petro_j) + binwidth, binwidth)
                    binsh = np.arange(np.min(filtered_petro_h), np.max(filtered_petro_h) + binwidth, binwidth)
                    binsk = np.arange(np.min(filtered_petro_k), np.max(filtered_petro_k) + binwidth, binwidth)

                    plt.clf()
                    plt.hist(filtered_petro_j, bins = binsj,\
                             facecolor='white', edgecolor = 'red',\
                                linestyle='--', label='Observed J')
                    plt.hist(filtered_petro_h, bins = binsh,\
                             facecolor='white', edgecolor = 'blue',\
                                linestyle='--', label='Observed H')
                    plt.hist(filtered_petro_k, bins = binsk,\
                             facecolor='white', edgecolor = 'green',\
                                linestyle='--', label='Observed K')
                    plt.xlabel('petro magnitudes')
                    plt.ylabel('Bins')
                    plt.grid()
                    plt.legend(loc='best')
                    plt.savefig('hist_ukidss_nir.png')
                    plt.clf()
                    
                    nir_data =  filtered_petro_j, filtered_petro_h, filtered_petro_k,\
                        e_petro_j, e_petro_h, e_petro_k, filtered_ukidss_ra, filtered_ukidss_dec
        return nir_data

    def read_gaia_data(self):
            """
                Reads the input GAIA DR3 data. The number of columns are 12.

            """
            header = ['source_id','ra','ra_error,dec','dec_error','parallax','parallax_error',\
                      'pm','pmra','pmra_error','pmdec','pmdec_error','ruwe']
            
            ra_name = str(self.ra).replace('.','_')
            dec_name = str(self.dec).replace('.', '_')    
            file_name = 'GAIA'+'_'+'RA'+str(ra_name)+'DEC'+str(dec_name)
            try:
                gaia_data = np.genfromtxt(str(file_name)+'.csv',delimiter=',',skip_header=1)
                gaia_source_id = gaia_data[:,0]
                gaia_ra = gaia_data[:,1]
                gaia_ra_error = gaia_data[:,2]
                gaia_dec = gaia_data[:,3]
                gaia_dec_error = gaia_data[:,4]
                gaia_parallax = gaia_data[:,5]
                gaia_parallax_error = gaia_data[:,6]
                gaia_pm = gaia_data[:,7]
                gaia_pm_ra = gaia_data[:,8]
                gaia_pm_ra_error = gaia_data[:,9]
                gaia_pm_dec = gaia_data[:,10]
                gaia_pm_dec_error = gaia_data[:,11]
                gaia_ruwe = gaia_data[:,12]
                gaia_data = gaia_source_id, gaia_ra, gaia_ra_error, gaia_dec,\
                    gaia_dec_error, gaia_parallax, gaia_parallax_error,\
                    gaia_pm, gaia_pm_ra, gaia_pm_ra_error, gaia_pm_dec,\
                    gaia_pm_dec_error, gaia_ruwe
            except FileNotFoundError:
                self.gd.get_gaia_data()
                gaia_data=np.genfromtxt(str(file_name)+'.csv',delimiter=',',\
                                        skip_header=1)
                gaia_source_id = gaia_data[:,0]
                gaia_ra = gaia_data[:,1]
                gaia_ra_error = gaia_data[:,2]
                gaia_dec = gaia_data[:,3]
                gaia_dec_error = gaia_data[:,4]
                gaia_parallax = gaia_data[:,5]
                gaia_parallax_error = gaia_data[:,6]
                gaia_pm = gaia_data[:,7]
                gaia_pm_ra = gaia_data[:,8]
                gaia_pm_ra_error = gaia_data[:,9]
                gaia_pm_dec = gaia_data[:,10]
                gaia_pm_dec_error = gaia_data[:,11]
                gaia_ruwe = gaia_data[:,12]
                gaia_data = gaia_source_id,gaia_ra,gaia_ra_error,gaia_dec,\
                    gaia_dec_error, gaia_parallax, gaia_parallax_error,\
                    gaia_pm, gaia_pm_ra, gaia_pm_ra_error, gaia_pm_dec,\
                    gaia_pm_dec_error, gaia_ruwe
            return gaia_data
