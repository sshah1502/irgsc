"""Testing Stellar Atmospheric Data models which are part
of the package as well as Models Class"""

from irgsctool import Models
from pathlib import Path

data_dir = Path(__file__).parent.joinpath()

def test_kurucz():
  """ Testing Kurucz and Phoenix model data"""
  
  
  k = Models('Kurucz')
  p = Models('Phoenix')
  
  k_m = k.read_sam_file()
  p_m = p.read_sam_file()
  
  k_sm = k.select_sam_range(teff_range=[4000,10000])
  p_sm = p.select_sam_range(teff_range=[2800,5000], logg_range=[3.0,5.5],
                                       feh_range=[-5.0,-1.5])
  
  assert len(k_m)==len(p_m)==11
  assert k_m[0].min() == 3500.0 and k_sm[0].min() > 4000.0
  assert p_m[0].min() == 2000.0 and p_sm[0].min() > 2800.0
  assert p_sm[1].min() > 3.0 and p_sm[2].min() > -5.0
  
