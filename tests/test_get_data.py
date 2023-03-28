""" This modules tests the data acquisition class.
Checks if the table has non-zero rows and expected
number of columns"""

from irgsctool import GetData
import pandas as pd
import os


gd = GetData(0.0,0.0)


def test_PS():
  """PanSTARRS Test"""
  tab = gd.get_panstarrs_data()
  assert len(tab)>1 and len(tab.keys())== 44
  
def test_UKIDSS():
  """UKIDSS Test"""
  gd.get_ukidss_data()
  assert os.path.exists('UKIDSS_RA0_0DEC0_0.csv')
  
  tab = pd.read_csv('UKIDSS_RA0_0DEC0_0.csv')
  
  assert len(tab)>1 and len(tab.keys()) == 8
  
def test_GAIA():
  """GAIA Test"""
  tab = gd.get_ukidss_data()
  
  assert os.path.exists('GAIA_RA0_0DEC0_0.csv')
  assert len(tab)>1 and len(tab.keys()) == 14
