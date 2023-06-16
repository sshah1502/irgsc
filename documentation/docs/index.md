# irgsctool

## Getting Started
<p style="text-align: justify;">This package is a tool to generate the catalog of Natural Guide Stars (NGS) based on the requirements of the NFIRAOS, the facility Adaptive Optics system on the Thirty Meter Telescope (TMT). This tool obtains the optical data from PANSTARRS DR2 and generates a catalog of the computed Near-Infrared (NIR) magnitudes. There is also an option in this tool to validate the computed NIR magnitudes using the readily available NIR observed data from UKIDSS DR11 (in the regions where the UKIDSS data is available.)</p>

## Installation
This package can be installed in two ways:

### 1. Using pip
In a fresh environment, enter the following command
```bash
pip install irgsctool
```

### 2. Uing .zip file on GitHub [here](https://github.com/tmtsoftware/dms-irgsc)  

After downloading the .zip file, unzip it into a directory and type
```bash
pip install .
```
This will install irgsctool locally in your system.

### 3. Using the Development version from GitHub:
Open the terminal in your system and run:
```
pip install git+https://github.com/tmtsoftware/dms-irgsc@main

```
This will also install the package locally in your system.
