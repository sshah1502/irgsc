"""Test for checking validate module
"""
from irgsctool import ValidateIRGSC
from irgsctool import GenerateIRGSC

g = GenerateIRGSC(0.0,0.0)
g.generate_irgsc()
v = ValidateIRGSC(0.0,0.0)

def test_read_irgsc():
  """
  Simple test to check table integrity as IRGSC 
  table is being read
  """
  
  irgsc_data = v.read_irgsc()
  
  assert len(irgsc_data)==68
  assert len(irgsc_data[0])>1
  
