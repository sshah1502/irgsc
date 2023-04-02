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

| Column Name | Description | Data Type  |
| :----------- |:------------|:------|
| PS1_ObjID    | Source identifier in the PANSTARRS data| float |
| PS1_ra       | Right Ascencion of the source in the PANSTARRS DR2 weighted mean photometry| float|
| PS1_ra_error | Uncertainty in PS1_ra| float|
| PS1_dec      | Declination of the source in the PANSTARRS DR2 weighted mean photometry| float|
| PS1_dec_error| Uncertainty in the PS1_dec| float|
| PS1_gpsf     | psf magnitude of the source in the g-band stacked photometry | float|
| PS1_gpsf     | Uncertainty in PS1_gpsf|
| PS1_rpsf     | psf magnitude of the source in the r-band stacked photometry | float|
| PS1_rpsf     | Uncertainty in PS1_rpsf | float|
| PS1_ipsf     | psf magnitude of the source in the i-band stacked photometry | float|
| PS1_ipsf     | Uncertainty in PS1_ipsf | float|
| PS1_zpsf     | psf magnitude of the source in the z-band stacked photometry | float|
| PS1_zpsf     | Uncertainty in PS1_zpsf | float|
| PS1_ypsf     | psf magnitude of the source in the y-band stacked photometry | float|
| PS1_ypsf     | Uncertainty in PS1_ypsf | float|
| SAM_Name     | Name of the best-fitted Stellar Atmospheric Model (SAM)| string|
| Teff         | Best-fitted model parameter: Teff| float|
| logg         | Best-fitted model parameter: log(g)| float|
| [Fe/H]       | Best-fitted model parameter: [Fe/H]| float|
|sam_g         | Best-fitted model magnitudes in PANSTARRS g-filter| float|
|sam_r         | Best-fitted model magnitudes in PANSTARRS r-filter| float|
|sam_i         | Best-fitted model magnitudes in PANSTARRS i-filter| float|
|sam_z         | Best-fitted model magnitudes in PANSTARRS z-filter| float|
|sam_y         | Best-fitted model magnitudes in PANSTARRS y-filter| float|
|sam_j         | Best-fitted model magnitudes in PANSTARRS j-filter| float|
|sam_h         | Best-fitted model magnitudes in PANSTARRS h-filter| float|
|sam_k         | Best-fitted model magnitudes in PANSTARRS k-filter| float|
| scale factor | The scale factor computed after fitting the SAM| float|
| scale factor error| Error in the computed scale factor| float|
| d_dev        | The parameter denoting the goodness-of-fit| float|
| Computed J   | The computed J magnitude in the Vega system| float|
| Computed J error| Error in computed J magnitude| float|
| Computed H   | The computed H magnitude in the Vega system| float|
| Computed H error| Error in computed H magnitude| float|
| Computed K   | The computed K magnitude in the Vega system| float|
| Computed K error| Error in computed K magnitude| float|
| gaia source id| Source identifier in Gaia DR3| float|
| gaia ra | Right Ascension of the source in Gaia DR3 catalog| float|
| gaia ra error| Uncertainty in gaia ra| float|
| gaia dec | Declination of the source in Gaia DR3 catalog| float|
| gaia dec error | Uncertainty in gaia dec | float|
| gaia parallax| Parallax (mas) of the source in the Gaia DR3 catalog| float|
| gaia parallax error| Uncertainty in gaia parallax| float|
| gaia pm| pm of the source (mas/yr) in Gaia DR3 catalog| float|
| gaia pm ra| pm of the source along R.A. axis in the Gaia DR3 catalog| float|
| gaia pm ra error| Uncertainty gaia pm ra | float|
| gaia pm dec| pm of the source along Dec. axis in the Gaia DR3 catalog| float|
| gaia pm dec error| Uncertainty gaia pm dec | float|
| gaia ruwe| Renormalised Unit Weight Error flag of the source in Gaia DR3| float|
| objinfoflag | These flag values of the source in PANSTARRS data specify whether the object is a QSO, transient, asteroid, extended, a known solar system object, etc. in nature| float|
| objqualityflag | These flag values denote if an object is real or a possible false positive | float|
| ndetections |The number of times something is detected from the individual exposures| float|
| nstackdetections | The number of stack detections after which the stack photometric measurements are done | float|
| ginfoflag | These flags indicate the details of the g filter stack photometry | float|
| ginfoflag2 | These flags indicate the details of the g filter stack photometry | float|
| ginfoflag3 | These flags indicate the details of the g filter stack photometry | float|
| rinfoflag | These flags indicate the details of the r filter stack photometry | float|
| rinfoflag2 | These flags indicate the details of the r filter stack photometry | float|
| rinfoflag3 | These flags indicate the details of the r filter stack photometry | float|
| iinfoflag | These flags indicate the details of the i filter stack photometry | float|
| iinfoflag2 | These flags indicate the details of the i filter stack photometry | float|
| iinfoflag3 | These flags indicate the details of the i filter stack photometry | float|
| zinfoflag | These flags indicate the details of the z filter stack photometry | float|
| zinfoflag2 | These flags indicate the details of the z filter stack photometry | float|
| zinfoflag3 | These flags indicate the details of the z filter stack photometry | float|
| yinfoflag | These flags indicate the details of the y filter stack photometry | float|
| yinfoflag2 | These flags indicate the details of the y filter stack photometry | float|
| yinfoflag3 | These flags indicate the details of the y filter stack photometry | float|

# Application of irgsctool on fields
The method developed for the generation of IRGSC has applied on twenty test fields across the sky. The generaed IRGSC is also valiated using the UKIDSS data available for those fields and the generated as well as validated catalog for these fields can be found in the 'generated_irgsc' directory.
In addition to the twenty test fields, additional ten catalogs are provided for the PANSTARRS Medium Deep Survey (MDS) Fields [more information available here](https://arxiv.org/abs/1612.05560). Since the MDS data is not publically released by the PANSTARRS, the optical data to generate the IRGSC for these fields is been taken from the 3-pi PANSTARRS survey.


# Requirements
This package is developed for Python versions above 3.6. It uses various other packages like:
astroquery,
astropy,
matplotlib,
astropy,
dustmaps,
numpy,
datetime,
requests, and 
pyvo.
### Note: It is recommended to install the package in a fresh environment.


# Installation

## 1. Using pip:
```
pip install irgsctool

```
## 2. Using .zip file from GitHub:
Download the .zip file from [here](https://github.com/sshah1502/irgsc) and unzip it. Then open the directory in terminal and type:
```
pip install .
```
## 3 Using the Development version from GitHub:
```
pip install git+https://github.com/sshah1502/irgsc@main
```

# Usage
```
 class GenerateIRGSC
```
This class is defined by importing irgsctool module and passing the R.A. and Decl. arguments. In this package, the catalog is generated using the optimal method described in the work report (link). After initializing, this module alerts the user if there is no observed NIR UKIDSS data for the given field.

```
from irgsctool import GenerateIRGSC as GC
gc = GC(ra(float),dec(float))
gc.generate_irgsc()
```

The module Generate_IRGSC is the module that generates the catalog in .csv format. Irrespective of whether UKIDSS data is available or not, this module (the command gc.generate_irgsc()) generates the catalog using the optical PANSTARRS data from 3pi steradian survey for given ra (float) and decl.(float). The name of the generated catalog has IRGSC prefix followed by RA, Dec., and the data of generation. eg. IRGSC_RA_0_0_DEC_0_0.csv for (ra,dec) = (0.0, 0.0)

```
class Validate
```

from irgssctool import Validate
vd = Validate(ra,dec)
vd.validate()
```
The module Validate(ra,dec) is the module that validates the computed NIR magnitudes after importing the IRGSC library. If initialized without generating the catalog, this module independantly checks whether the UKIDSS observed NIR data can be obtained for the given field.

# Conclusion/Disclaimer

Please add the following acknowledgment if you use our package in your work.

"This work has made use of "tmt-irgsc" developed as part of the Thirty Meter Telescope (TMT) project."

If you have any questions or suggestions for improvements to this repo,
please contact the owners of the repository.
