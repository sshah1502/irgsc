import os
import sys
from irgsctool import GenerateIRGSC as GC
from irgsctool import ValidateIRGSC as VG


gc = GC(227.26,0.0)
gc.generate_irgsc()
#vd = VG(53.1, -27.8)
#vd.validate(validate=True)