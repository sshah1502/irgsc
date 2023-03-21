#!/usr/bin/env python
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

from datetime import date
import os
from ._fitting import Generate_IRGSC
from ._validate import Validate
from ._extinction_correction import Extinction_Correction
from ._read_data import Read_Data
from ._get_data import Get_Data
from ._sgc import Star_Galaxy_Classification
from ._sam import Models
import numpy as np
from dustmaps.config import config
config['data_dir'] = Path(__file__).parent.joinpath()
import dustmaps.sfd
if "sfd" not in os.listdir(config['data_dir']):
    dustmaps.sfd.fetch()

home_dir = os.getcwd()
print('')
print('Your home directory is:', home_dir)
print('####################################################')
print('')
__author__ = "Sarang Shah"
__copyright__ = "Copyright 2023, TMT/DMS/IRGSC"
__credits__ = "Sarang Shah"
__license__ = "MIT"
__version__ = "0.0.1"
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

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

"""

class irgsc(Generate_IRGSC, Validate, Extinction_Correction,\
            Read_Data, Get_Data, Star_Galaxy_Classification):

    print ('##########################################################')
    print("")
    print('Initializing')
    print("")
    print ('##########################################################')
    print("")

    def __init__(self, ra, dec):
        print("")
        print('Checking whether UKIDSS data is available for this field')
        print('')
        print('##################################################')
        print('')
        ra_name = str(ra).replace('.','_'); dec_name = str(dec).replace('.', '_')
        file_name = 'UKIDSS' + '_' + 'RA'+str(ra_name) + 'DEC' + str(dec_name)
        try:
            validating_data = np.genfromtxt(str(home_dir)+ '/'+ str(file_name) + '.csv')
            if len(validating_data) == 0.0:
                print("""UKIDSS Observed NIR data not available. Validation of the 
                generated IRGSC not possible for this field!!!""")
            else:
                validating_data = validating_data
        except FileNotFoundError:
            Get_Data.get_ukidss_data(ra, dec)
            validating_data = np.genfromtxt(str(home_dir)+ '/'+ str(file_name) + '.csv')
            if len(validating_data) == 0.0:
                print("""UKIDSS Observed NIR data not available. Validation of the
                generated IRGSC not possible for this field!!!""")
            else:
                validating_data = validating_data

        file_name = 'PS1' + '_' + 'RA'+str(ra_name) + 'DEC' + str(dec_name)
       
        try:
            optical_data = np.genfromtxt(str(home_dir)+ '/'+ str(file_name) + '.csv')
            if len(optical_data) == 0.0:
                print('')
                print("""PANSTARRS optical data not available. 
                Please check the input coordinates!!!""")
                print('')
                Get_Data.get_panstarrs_data(ra, dec)
                optical_data = np.genfromtxt(str(home_dir)+ '/'+ str(file_name) + '.csv')
            else:
                optical_data = optical_data
        except FileNotFoundError:
            Get_Data.get_panstarrs_data(ra, dec)
            optical_data = np.genfromtxt(str(file_name) + '.csv')
            if len(optical_data) == 0.0:
                print("""PANSTARRS optical data not retrieved. Please check the 
                input coordinates or the ADSQL query!!!""")
            else:
                optical_data = optical_data

        file_name = 'GAIA' + '_' + 'RA'+str(ra_name) + 'DEC' + str(dec_name)
        try:
            gaia_data = np.genfromtxt(str(home_dir)+ '/'+ str(file_name) + '.csv')
            if len(gaia_data) == 0.0:
                print('GAIA data not available for this field!!!')
            else:
                gaia_data = gaia_data
        except FileNotFoundError:
            Get_Data.get_gaia_data(ra, dec)
            gaia_data = np.genfromtxt(str(home_dir)+ '/'+str(file_name) + '.csv')
            if len(gaia_data) == 0.0:
                print('GAIA data not available for this field!!!')
            else:
                gaia_data = gaia_data

