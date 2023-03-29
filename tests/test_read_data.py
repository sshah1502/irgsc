"""This modules tests if the downloaded data is
readable and processable"""

from irgsctool import ReadData

rd = ReadData(0.0,0.0)


def test_optical_data():
  opt_tab = rd.read_optical_data()
  
  assert len(opt_tab)>1
  
def test_nir_data():
  nir_tab = rd.read_nir_data()
  
  assert len(nir_tab)>1
  
def test_gaia_data():
  gaia_tab = rd.read_gaia_data()
  
  assert len(gaia_tab)>1
  
