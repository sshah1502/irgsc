# TMT-IRGSC

## Introduction
This is a python repository dedicated to the development of the Near-Infrared Guide Star Catalog (IRGSC) for the Adaptive Optics (AO) observations of the Thirty Meter Telescope (TMT) project. This package generates the catalog by computing the NIR magnitudes of the optical stellar sources in the PANSTARSS DR2 data.

# Packages
## irgsctool:
<p style="text-align: justify;">This is a python package aimed to compute the NIR magnitudes of the optical sources in the PANSTARRS stack-photometric data by modelling them with the Kurucz, Castelli-Kurucz and Phoenix stellar atmospheric models. This package also validates the computed NIR magnitudes with the observed NIR data from UKIDSS (if it is available). The methodology implemented in this python package is implemented on twenty test-fields across the TMT's observable sky and the generated as well as validated IRGSC is available on the GitHub homepage. Most of the sources have the computed NIR magnitudes similar to the observed. The generated catalog contains astrometric information from GAIA DR3 also. Read the section below to see the nature of IRGSC.

### Nature of the generated catalog
The IRGSC generated has various information about the sources shown in the following Table. This table describes the columns in the IRGSC generated for a particular test field. The details of the flags, e.g., infoflags, filterflags, and qualityflags can be found [here](https://outerspace.stsci.edu/display/PANSTARRS/PS1+StackObjectView+table+fields). These flags indicate various values assigned to
the source by the PANSTARRS team, which gives further information about the nature of the source
and the quality of its detection, which can help understand more about a particular object of interest.
It is to be noted that although this package relies on the PANSTARRS StackObjectView table, the Right
Ascension and Declination of the source is obtained from the mean photometric information as they are well calibrated using Gaia DR2.</p>

| Column Name | Description | Type  |
| :----------- |:------------:|:------|
| PS1_ObjID    | Object ID in the PANSTARRS data|
| PS1_ra       | R.A. of the source in PS1 weighted mean photometry


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
