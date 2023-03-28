""" This modules tests the data acquisition class.
Checks if the table has non-zero rows and expected
number of columns"""

from irgsctool import Get_Data
from astropy.table import Table
import os


gd = Get_Data(227.26,0)


def test_PS():
  """PanSTARRS Test"""
  tab = gd.get_panstarrs_data()
  assert len(tab)>1 and len(tab.keys())== 44
  
def test_UKIDSS():
  """UKIDSS Test"""
  gd.get_ukidss_data()
  assert os.path.exists('UKIDSS_RA227_26DEC0.csv')
  
  tab = Table.read('UKIDSS_RA227_26DEC0.csv', format = 'csv')
  
  
  assert len(tab)>1 and len(tab.keys()) == 8
  
def test_GAIA():
  """GAIA Test"""
  gd.get_ukidss_data()
  assert os.path.exists('GAIA_RA0_0DEC0_0.csv')
  
  tab = Table.read('GAIA_RA227_26DEC0.csv', format = 'csv')
  
  assert len(tab)>1 and len(tab.keys()) == 14
