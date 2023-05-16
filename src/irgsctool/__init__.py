#pylint: disable=wrong-import-position
import sys
import os
from datetime import date
import warnings
import numpy as np
from dustmaps.config import config
config['data_dir'] = os.getcwd()

from ._read_data import ReadData
from ._get_data import GetData
from ._fitting import GenerateIRGSC
from ._validate import ValidateIRGSC
from ._extinction_correction import ExtinctionCorrection
from ._sgc import StarGalaxyClassification
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


class irgsc(GetData, ReadData, StarGalaxyClassification, ExtinctionCorrection, Models,
            GenerateIRGSC, ValidateIRGSC):
    """
    ------------------------------------------
    Initialisation of child *** irgsc class ***. This class has several child classes.
    """
    print ('##########################################################################')
    print("")
    print('Initializing')
    print("")
    print ('##########################################################################')
    print("")

    def __init__(self, ra, dec, validate=None):
        """
        This method described using input  ra and dec.
        It checks whther the data from PANSTARRS DR2, UKISS DR11 and
        Gaia DR3 can be obtained.
        Raises:
            ValueError: if the data is not available in UKIDSS or
                        PANSTARRS 3-pi survey. The code will not\
                             proceed further.
            Warning: if the data is not available in the Gaia DR3\
                        survey. In the IRGSC, the Gaia values will\
                            be replaces with -999.

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
        #gd = GetData(self.ra, self.dec)
        #rd = ReadData(self.ra, self.dec)
        if self.ra < 0.0 or self.dec<-30.0:
            raise ValueError('Please check the input coordinates')
            sys.exit(0)