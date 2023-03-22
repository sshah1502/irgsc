from irgsctool import Get_Data
import os

gd = Get_Data(227.26,0)
def test_PS():
  tab = gd.get_panstarrs_data()
  assert len(tab)>1 and len(tab.keys())== 44
  
def test_UKIDSS():
  gd.get_ukidss_data()
  assert os.path.exists('UKIDSS_RA227_26DEC0.csv')
 
