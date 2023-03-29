"""Test for checking SED fitting
"""
from irgsctool import GenerateIRGSC
from irgsctool._fitting import find_nearest, calc_sf, compute_dquad


def test_find_nearest():
  A = [1,2,3,4,5,6]
  B = 3.6
  
  out = find_nearest(A,B)  
  assert out == 4
   
