import os
import sys
import numpy as np
data_dir = os.getcwd()
class Models():
        def read_sam_file(use_sam=None):

                if use_sam == 'Kurucz':
                        print("")
                        print('Reading Interpolated Kurucz SAMs')
                        print("")
                        print('data_dir=', os.getcwd())

                        p2 = np.genfromtxt(str(os.getcwd())+'/'+'src/irgsctool/data/'+'interpolated_kurucz.txt')
                elif use_sam == 'Phoenix':
                        print("")
                        print('Reading Interpolated Phoenix SAMs')
                        print("")
                        p2 = np.genfromtxt(str(os.getcwd())+'/'+'src/irgsctool/data/'+'interpolated_phoenix.txt')
        
                teff = p2[:,0]; logg = p2[:,2]; feh = p2[:,1]; sam_g = p2[:,3]; sam_r = p2[:,4];\
                        sam_i = p2[:,5]; sam_z = p2[:,6]; sam_y = p2[:,7]; sam_j = p2[:,8];\
                        sam_h = p2[:,9]; sam_k = p2[:,10]
        
                sam_params =  teff, logg, feh, sam_g, sam_r, sam_i, sam_z, sam_y, sam_j, sam_h, sam_k
                return sam_params

        def select_sam(teff_range=None, logg_range=None, feh_range=None, use_sam =None, use_optimal_method=False):
                teff, logg, feh, sam_g, sam_r, sam_i, sam_z, sam_y, sam_j, sam_h, sam_k\
                                        = Models.read_sam_file(use_sam=use_sam)
                if use_optimal_method is True:
                        if teff_range is None and logg_range is None and feh_range is None:
                                raise  TypeError("Parameter range must be provided in order to use\
                                         the optimal method")
                        elif teff_range is not None and feh_range is not None and logg_range is not None:
                                teff_lower_limit = teff_range[0]; teff_upper_limit = teff_range[-1];\
                                logg_lower_limit = logg_range[0]; logg_upper_limit = logg_range[-1];\
                                feh_lower_limit = feh_range[0]; feh_upper_limit = feh_range[-1]

                                index_sam = np.where((teff>teff_lower_limit) & (teff<teff_upper_limit)\
                                   & (feh>feh_lower_limit) & (feh<feh_upper_limit)\
                                        & (logg>logg_lower_limit) & (logg<logg_upper_limit))[0]
                
                                sam_params =  teff[index_sam], logg[index_sam], feh[index_sam], sam_g[index_sam], sam_r[index_sam],\
                                                sam_i[index_sam], sam_z[index_sam], sam_y[index_sam], sam_j[index_sam],\
                                                        sam_h[index_sam], sam_k[index_sam]
                                
                        elif teff_range is not None and feh_range is None and logg_range is None:
                                teff_lower_limit = teff_range[0]; teff_upper_limit = teff_range[-1]
                                index_sam = np.where((teff>teff_lower_limit) & (teff<teff_upper_limit))[0]
                                sam_params =  teff[index_sam], logg[index_sam], feh[index_sam], sam_g[index_sam], sam_r[index_sam],\
                                                sam_i[index_sam], sam_z[index_sam], sam_y[index_sam], sam_j[index_sam],\
                                                        sam_h[index_sam], sam_k[index_sam]                                
                                
                                
                        return sam_params