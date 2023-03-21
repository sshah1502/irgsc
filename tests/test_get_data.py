from irgsctool import Get_Data
import os

def test_PS():
  tab = Get_Data.get_panstarrs_data(227.26,0)
  assert len(tab)>1 and len(tab.keys())== 44
  
def test_UKIDSS():
  Get_Data.get_ukidss_data(227.26,0)
  assert os.path.exists('UKIDSS_RA227_26DEC0.csv')
  
