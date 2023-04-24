import os
import sys
from irgsctool import GenerateIRGSC as GC
from irgsctool import ValidateIRGSC as VG

#gc = VG(0.0, 0.0)
#gc.validate(validate=True)
ic = GC(334.27, 0.38)
ic.generate_irgsc(use_optimal_method=True)