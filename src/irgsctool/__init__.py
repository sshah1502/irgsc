"""
    IRGSC, aa python pacakge to generate a NIR catalog of the stars in the TMT's observable sky.
    version - v1.0
    Author - Sarang Shah
    Please add the following acknowledgment if you use our package in your work.

            "This work has made use of Infrared Guide Star Catalog (IRGSC) developed 
            as a part of the Thirty Meter Telescope (TMT) project."

    If you have any questions or suggestions for improvements to this repo,
    please email: sarang.itcc@iiap.res.in.
    
"""
#pylint: disable=wrong-import-position
import sys
import os
from datetime import date
import numpy as np
from dustmaps.config import config
config['data_dir'] = os.getcwd()
import dustmaps.sfd
from ._fitting import Generate_IRGSC
from ._validate import Validate
from ._extinction_correction import Extinction_Correction
from ._read_data import Read_Data
from ._get_data import Get_Data
from ._sgc import Star_Galaxy_Classification
from ._sam import Models


home_dir = os.getcwd()
print('')
print('Your home directory is:', home_dir)
print('####################################################')
print('')
__author__ = "Sarang Shah"
__copyright__ = "Copyright 2023, TMT/DMS/IRGSC"
__credits__ = ["Sarang Shah"]
__license__ = "MIT"
__version__ = "0.1.1"
__maintainer__ = "Sarang Shah"
__email__ = "sarang.itcc@iiap.res.in"
__status__ = "Development"

"""MIT License

Copyright (c) [2023] [tmt-irgsc]

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARraNTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARraNTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTraCT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

"""

class Irgsc(Generate_IRGSC, Validate, Extinction_Correction,
            Read_Data, Get_Data, Star_Galaxy_Classification):
    """
    Initialisation of Irgsc Class
    """
    print ('##########################################################################')
    print("")
    print('Initializing')
    print("")
    print ('##########################################################################')
    print("")

    def __init__(self, ra, dec, validate=None):
        """
        __init__ function described using ra and decl.
        Checks whther the data from PANSTARRS, UKISS and
        Gaia can be obtained for the given set of ra and decl.

        """
        print('')
        print('##################################################')
        print('Checking the input coordinates')
        print('')
        print('##################################################')
        print('')
        self.ra = ra
        self.dec = dec
        self.validate=validate
        if self.ra < 0.0 or self.dec<-30.0:
            raise ValueError('Please check the input coordinates')
            sys.exit(0)

        ra_name = str(self.ra).replace('.','_')
        dec_name = str(self.dec).replace('.', '_')


        #Checking whether UKIDSS data is available for the given field.
        ###If yes then obtaining it.

        file_name = 'UKIDSS'+'_'+'ra'+str(ra_name)+'DEC'+str(dec_name)
        try:
            validating_data = np.genfromtxt(file_name+'.csv')
        except FileNotFoundError:
            Get_Data.get_ukidss_data(self.ra, self.dec)
            validating_data = Read_Data.read_nir_data(self.ra, self.dec)
            self.validate = True
            if len(validating_data) == 0.0:
                self.validate = False
                raise ValueError('UKIDSS Observed NIR data not available.\
                                 Validation of the generated IRGSC is not\
                                 possible for this field!!!')
            sys.exit(0)

        ###Obtaining the PANSTARRS data for the given field.

        file_name ='PS1'+'_'+'ra'+str(ra_name)+'DEC'+str(dec_name)
        try:
            optical_data = np.genfromtxt(file_name + '.csv')
        except FileNotFoundError:
            Get_Data.get_panstarrs_data(self.ra, self.dec)
            optical_data = np.genfromtxt(file_name + '.csv')
            if len(optical_data) == 0.0:
                raise ValueError('PANSTARRS optical data not available.\
                                 Please check the input coordinates!!!')
            sys.exit(0)

        #Obtaining the GAIA data for the given field

        file_name = 'GAIA' + '_' + 'ra'+str(ra_name) + 'DEC' + str(dec_name)
        try:
            gaia_data = np.genfromtxt(file_name + '.csv')
        except FileNotFoundError:
            Get_Data.get_gaia_data(self.ra, self.dec)
            gaia_data = np.genfromtxt(file_name + '.csv')
            if len(gaia_data) == 0.0:
                print('GAIA data not available for this field!!!')
