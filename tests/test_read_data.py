"""This modules tests if the downloaded data is
readable and processed properly"""

from irgsctool import ReadData

rd = ReadData(0.0,0.0)


def test_optical_data():
  """PanSTARRS"""
  opt_tab = rd.read_optical_data()
  
  assert len(opt_tab)==44  # No of Columns
  
def test_nir_data():
  """UKIDSS"""
  nir_tab = rd.read_nir_data()
  
  assert len(nir_tab)==8  # No of Columns
  
def test_gaia_data():
  """GAIA"""
  gaia_tab = rd.read_gaia_data()
  
  assert len(gaia_tab)==13   # No of Columns
  
