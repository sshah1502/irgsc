# irgsctool (mock for tmt-irgsc)

## Introduction
This is a python repository dedicated to the development of the Near-Infrared Guide Star Catalog (IRGSC) for the Adaptive Optics (AO) observations of the Thirty Meter Telescope (TMT) project. This package generates the catalog by computing the NIR magnitudes of the optical stellar sources in the PANSTARSS DR2 data.

# Packages
## tmt-irgsc:
A python package aimed to compute the NIR magnitudes of the optical sources in the PANSTARRS data by modelling them with the Kurucz, Castelli-Kurucz and Phoenix Stellar Atmospheric models. This package also validates the computed NIR magnitudes with the observed NIR data from UKIDSS (if it is available). The methodology implemented in this paython package is implemented on twenty test-fields across the TMT's observable sky. Most of the sources have the computed NIR magnitudes similar to the observed. The generated catalog contains astrometric information from GAIA DR3 as well.


### Installation
```
pip install irgsctool

```

# Usage
```
 class IRGSC
```
This class is defined by importing irgsctool module and passing the R.A. and Decl. arguments. In this package, the catalog is generated using the optimal method described in the work report (link). After initializing, this module alerts if there is no observed NIR UKIDSS data for the given field.

```
from irgsctool import Generate_IRGSC as GC
gc = GC(ra,dec)
gc.generate_irgsc()
```

The module Generate_IRGSC is the module that generates the catalog. Irrespective of whether UKIDSS data is available or not, this module (the command gc.generate_irgsc()) generates the catalog using the optical PANSTARRS data from 3pi steradian survey for given ra and decl.

```
import Validate
vd = Validate(ra,dec)
vd.validate()
```
The module Validate(ra,dec) is the module that validates the computed NIR magnitudes after importing the IRGSC library. If initialized without generating the catalog, this module independantly checks whether the UKIDSS observed NIR data can be obtained for the given field.

# Conclusion/Disclaimer

Please add the following acknowledgment if you use our package in your work.

"This work has made use of "tmt-irgsc" developed as part of the Thirty Meter Telescope (TMT) project."

If you have any questions or suggestions for improvements to this repo,
please contact the owners of the repository.
