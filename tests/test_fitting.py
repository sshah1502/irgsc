"""Test for checking SED fitting
"""
from irgsctool import GenerateIRGSC
from irgsctool._fitting import find_nearest, calc_sf, compute_dquad
from datetime import date
import pandas as pd
import os

def test_find_nearest():
  A = [1,2,3,4,5,6]
  B = 3.6
  
  out = find_nearest(A,B)  
  assert out == 4
 
def test_GenIRGSC():
  gc = GenerateIRGSC(0.0,0.0)
  gc.generate_irgsc()
  
  ra_name=str(gc.ra).replace('.','_')
  dec_name=str(gc.dec).replace('.', '_')
  
  current_datetime = date.today()
  file_name = 'IRGSC_'+'RA'+str(ra_name)+'DEC'+str(dec_name)+ str(current_datetime)+'.csv'
  
  df = pd.read_csv(file_name)
  
  assert os.path.exists(file_name)
  assert len(df)>1 and len(df.keys()) == 68
