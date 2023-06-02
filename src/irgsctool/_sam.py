#pylint: disable=wrong-import-position
#pylint: disable=import-error, C0103, R0914, W0311, C0114, C0301, R0903, W1401
from pathlib import Path
import numpy as np

data_dir = Path(__file__).parent.joinpath()

class Models():
        """

                ***Models*** class reads and selects the required
                Kurucz and Phoenix stellar atmospheric models in the
                generation of IRGSC.

        """
        def __init__(self, use_sam=None):
                self.use_sam = use_sam


        def read_sam_file(self):
                """
                        `irgsctool.Models.read_sam_file(use_sam=None)`
                
                        <justify> 
                        This method reads the model parameters and results from the
                        interpolated Kurucz or the Phoenix model files. use_sam is 
                        bool and decides which model file to be read.
                        To use the interpolated Kurucz models, set ***use_sam = Kurucz ***.
                        Similarly, *** use_sam = Phoenix *** to use interpolated Phoenix models.</justify>

                        Raises:
                                AttributeError: if use_sam is None.
                                FileNotFoundError: if the model files are not found.
                        Returns:
                            sam_params: A multi-dimwnsional array of the sam parameters
                                depending upon the model selected.

                """
                if self.use_sam is None:
                        raise AttributeError('Stellar Atmospheric Model to be used not specified.')
                if self.use_sam == 'Kurucz':
                        print("")
                        print('Reading Interpolated Kurucz SAMs')
                        print("")
                        try:
                                p2 = np.genfromtxt(str(data_dir) +'/data/interpolated_kurucz.txt')
                        except Exception:
                                FileNotFoundError('interpolated_kurucz.txt file not found')
                elif self.use_sam == 'Phoenix':
                        print("")
                        print('Reading Interpolated Phoenix SAMs')
                        print("")
                        try:
                                p2 = np.genfromtxt(str(data_dir)+'/data/interpolated_phoenix.txt')
                        except Exception:
                                FileNotFoundError('interpolated_phoenix.txt not found')
                teff = p2[:,0]
                logg = p2[:,2]
                feh = p2[:,1]
                sam_g = p2[:,3]
                sam_r = p2[:,4]
                sam_i = p2[:,5]
                sam_z = p2[:,6]
                sam_y = p2[:,7]
                sam_j = p2[:,8]
                sam_h = p2[:,9]
                sam_k = p2[:,10]

                sam_params =  teff, logg, feh, sam_g, sam_r, sam_i, sam_z, sam_y, sam_j, sam_h, sam_k
                return sam_params

        def select_sam_range(self, teff_range=None, logg_range=None, feh_range=None):
                """
                        
                        `irgsctool.Models.select_sam_range(teff_range=None, logg_range=None,
                        feh_range=None, use_optimal_method=False)`
                        
                        <justify> This method selects the range of the models to be
                        used in the generation of IRGSC.
                       
                        <justify> If *** use_optimal_method *** is set to True, the following 
                        range of model parameters is selected:</justify>

                        | Model Name    |               \(T_{eff}\) (K)   |    log(g) (dex)  |    [Fe/H] (dex)|
                        | :-------------| :-------------------------| :----------------| :--------------|
                        | Phoenix (C1)  |                2800 - 5000 |    3.0 - 5.5    |     -5.0 - -1.5| 
                        | Phoenix (C2)  |              2800 - 4000   |   0.0 - 3.0     |     -0.5 - 1.5 | 
                        | KuruczCastelli-Kurucz (K0) |  4000 - 10000  |   ---         |       ---      |

                        Raises:
                                TypreError: if range of parameters is not given.
                        Returns:
                                sam_params : A multi-dimwnsional array of the selected model parameters and model 
                                magnitudes to generate IRGSC.

                """

                teff, logg, feh, sam_g, sam_r, sam_i, sam_z, sam_y, sam_j, sam_h, sam_k\
                                        = self.read_sam_file()
                if teff_range is None and logg_range is None and feh_range is None:
                        raise  TypeError("At least one parameter range must be provided")
                if teff_range is not None:
                                teff_lower_limit = teff_range[0]
                                teff_upper_limit = teff_range[-1]
                                c1 = (teff>teff_lower_limit) & (teff<teff_upper_limit)
                else:
                                c1 = 1

                if logg_range is not None:
                                logg_lower_limit = logg_range[0]
                                logg_upper_limit = logg_range[-1]
                                c2 = (logg>logg_lower_limit) & (logg<logg_upper_limit)
                else:
                                c2 = 1

                if logg_range is not None:
                                feh_lower_limit = feh_range[0]
                                feh_upper_limit = feh_range[-1]

                                c3 = (feh>feh_lower_limit) & (feh<feh_upper_limit)
                else:
                                c3 = 1

                index_sam = np.where(c1 & c2 & c3)[0]

                sam_params =  teff[index_sam], logg[index_sam], feh[index_sam], sam_g[index_sam], sam_r[index_sam],\
                                        sam_i[index_sam], sam_z[index_sam], sam_y[index_sam], sam_j[index_sam],\
                                                sam_h[index_sam], sam_k[index_sam]
                return sam_params
