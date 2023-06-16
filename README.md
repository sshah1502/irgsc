# TMT-IRGSC

## Introduction
This is a python repository dedicated to the development of the Near-Infrared Guide Star Catalog (IRGSC) for the Adaptive Optics (AO) observations of the Thirty Meter Telescope (TMT) project. This package generates the catalog by computing the NIR magnitudes of the optical stellar sources in the PANSTARSS DR2 data.
The details of this work can be found in the __phase_III_report.pdf__ file.
# Packages
## irgsctool:
<p style="text-align: justify;">This is a python package aimed to compute the NIR magnitudes of the optical sources in the PANSTARRS stack-photometric data by modelling them with the Kurucz, Castelli-Kurucz and Phoenix stellar atmospheric models. This package also validates the computed NIR magnitudes with the observed NIR data from UKIDSS (if it is available). The methodology implemented in this python package is implemented on twenty test-fields across the TMT's observable sky and the generated as well as validated IRGSC is available in the generated_irgsc directory (seec details below). Most of the sources have the computed NIR magnitudes similar to the observed. The generated catalog contains astrometric information from GAIA DR3 also. Read the section below to see the nature of IRGSC.

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
| SAM | The name of the best-fitted Stellar Atmospheric Model (SAM) | string|

# Nature of the validated catalog using the UKIDSS Data
The computed NIR magnitudes of the sources in the IRGSC can also be validated using the readily available UKIDSS data (if any) for the given field. irgsctool first checks whether a validated IRGSC can be produced for a given field and alerts the user accordingly. The table below shows the additional columns in an validated IRGSC.
| Column Name | Description | Data Type  |
| :----------- |:------------|:------|
| diff_J    	| Difference in the observed and computed J| float |
| diff_H       | Difference in the observed and computed H| float|
| diff_K 		| Difference in the observed and computed K| float|
| J_UKIDSS      | Observed J| float|
| err_J_UKIDSS| Uncertainty in observed J| float|
| H_UKIDSS     | Observed H | float|
| err_H_UKIDSS     | Uncertainty in observed H | float |
| K_UKIDSS     | Observed K | float|
| err_K_UKIDSS     | Uncertainty in observed K | float|


# Application of irgsctool on fields
The method developed for the generation of IRGSC has applied on twenty test fields (see the following table) across the sky. The generaed IRGSC is also valiated using the UKIDSS data available for those fields and the generated as well as validated catalog for these fields can be found in the 'generated_irgsc' directory.
In addition to the twenty test fields, additional ten catalogs are provided for the PANSTARRS Medium Deep Survey (MDS) Fields ([more information available here](https://arxiv.org/abs/1612.05560)). Since the MDS data is not publically released by the PANSTARRS, the optical data to generate the IRGSC for these fields is been taken from the PANSTARRS $3\pi$-survey.
| R.A. | Decl. | l | b| E(B-V) |
| :----------- |:------------|:------|:------|:------|
|227.26	    |0.0		|359.27	    |47.24		|0.04|
|334.27		|0.38		|63.08	    |-43.84		|0.07|
|60.00		|1.25		|188.72	    |-36.53		|0.26|
|30.00		|0.50		|156.53	    |-57.82		|0.02|
|11.16	    |7.83		|120.00		|-55.00		|0.04|
|225.53	    |2.19		|0.00		|50.00		|0.04|
|269.93	    |-13.48	    |15.00		|5.00		|0.98|
|334.80	    |50.96	    |100.00		|-5.00		|0.28|
|324.09     |51.47	    |95.00		|-0.50		|2.48|
|298.02	    |34.02		|70.00		|3.00		|1.01|
|0.00   	|0.00		|96.33	    |-60.18		|0.02|
|34.50		|-5.16		|169.97	    |-59.87		|0.01|
|36.25		|-4.50		|171.65	    |-58.22		|0.02|
|164.25		|57.66		|148.39	    |53.43		|0.04|
|66.75		|15.86		|180.08	    |-22.32		|0.58|
|82.25		|-2.60		|205.62	    |-19.48		|0.62|
|189.83 	|0.00		|296.33	    |62.71		|0.01|
|150.25		|10.00		|227.71		|46.40		|0.03|
|15.0		|0.90		|127.47	    |-61.89		|0.02|
|35.0		|-3.50		|168.62 	|-58.28     |0.01|

The generated_irgsc directory contains following sub-directories:
### 1. irgsc_using_casjobs_30_arcmin_radius:
This directory contains the generated and validated irgsc using the PANSTARRS data for 30 arcmin radius downloaded from MAST CasJobs.

### 2. irgsc_using_casjobs_35_arcmin_radius:
This directory contains the generated and validated irgsc using the PANSTARRS data for 35 arcmin radius (1.06 sq. deg.) fields as the requirement was for 35 arcmin size.

### 3. irgsc_for_fields_close_to_galactic_plane:
This directory contains generated and validated irgsc for the fields located close to the galactic plane fields. Due to high source density in these fields, we have downloaded the PANSTARRS data for 5 arcmin size radius.

### 4. irgsc_using_pyvo_15_arcmin_radius:
This directory contains generated and validated irgsc for the PANSTARRS Medium Deep Survey fields with the optical data obtained using pyvo. Due to limitations of pyvo, the data can only be obtained for 15 arcmin radius area of the field.

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
### Note: It is recommended to install irgsctool in a fresh environment and requires a stable internet connection to fetch the data.


# Installation

## 1. Using .zip file from GitHub:
Download the .zip file from [here](https://github.com/sshah1502/irgsc) and unzip it. Then open the directory in terminal and type:
```
pip install .
```
## 2. Using the Development version from GitHub:
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

"This work has made use of "irgsctool" developed as part of the Thirty Meter Telescope (TMT) project."

If you have any questions or suggestions for improvements to this repo,
please email: sarang.itcc@iiap.res.in
