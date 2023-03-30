""" This modules tests the data acquisition class.
Checks if the table has non-zero rows and expected
number of columns"""

from irgsctool import GetData
import os

gd = GetData(0.0,0.0)


def test_PS():
  """PanSTARRS Test"""
  tab = gd.get_panstarrs_data()
  
  assert os.path.exists('PS1_RA0_0DEC0_0.csv')
  assert len(tab)>1 and len(tab.keys())== 44
  
def test_UKIDSS():
  """UKIDSS Test"""
  tab = gd.get_ukidss_data()
  
  assert os.path.exists('UKIDSS_RA0_0DEC0_0.csv') 
  assert len(tab)>1 and len(tab.keys()) == 9
  
def test_GAIA():
  """GAIA Test"""
  tab = gd.get_gaia_data()
  
  assert os.path.exists('GAIA_RA0_0DEC0_0.csv')
  assert len(tab)>1 and len(tab.keys()) == 14
