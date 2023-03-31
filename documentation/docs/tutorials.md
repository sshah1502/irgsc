## irgsctool - tutorials

## 1. Generating IRGSC for a given set of coordinates

``` python
from irgsctool import GenerateIRGSC
gc = GenerateIRGSC(ra=0.0, dec=0.0)
gc.generate_irgsc()
```
This results in a catalog file: 'IRGSC_R_A_0_0_DEC_0_0.csv' that contains probable stellar sources in 0.25 degrees radius around the input coordinates.