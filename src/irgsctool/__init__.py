"""
irgsc class checks whether IRGSC can be generated for the given set of coordinates
"""
#python.linting.pylintArgs: [ "--disable=missing-module-docstring", \
# "--disable=missing-class-docstring", "--disable=missing-function-docstring"]
#pylint: disable=wrong-import-position
# pylint: disable=super-init-not-called

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
__version__ = "0.1"
__maintainer__ = "Sarang Shah"
__email__ = "sarang.itcc@iiap.res.in"
__status__ = "Development"


class irgsc(GetData, ReadData, StarGalaxyClassification, ExtinctionCorrection, Models,
            GenerateIRGSC, ValidateIRGSC):
    """
    ------------------------------------------
    Initialisation of *** irgsc class ***. This class is dependant on several other classes
    like GetData, ReadData, StarGalaxyClassification, ExtinctionCorrection, Models, GenerateIRGSC,
    ValidateIRGSC.
    """
    def __init__(self, ra, dec, validate=None):
        """
        This method describes using input ra and dec.
        It checks whther the IRGSC can be generated for a given set of input coordinates.
        Raises:
            ValueError: if the data is not available in UKIDSS or
                        PANSTARRS 3-pi survey. The code will not\
                             proceed further.

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
        