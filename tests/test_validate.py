"""Test for checking validate module
"""
from irgsctool import ValidateIRGSC

v = ValidateIRGSC(0.0,0.0)

def test_read_irgsc():
  """
  Simple test to check table integrity as IRGSC 
  table is being read
  """
  
  tab = v.read_irgsc()
  
  assert len(tab)>1
  assert len(tab.keys())==68
