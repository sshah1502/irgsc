from astroquery.ukidss import Ukidss
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.table import Table
import numpy as np
import csv
import pyvo as vo
from astroquery.mast import Mast
from astropy.table import Table
from astroquery.gaia import Gaia
import sys

"""
Functions to obtain PANSTARRS DR2 Stacked Photometry data and GAIA DR3 astrometry data.

The default search radius is 0.25 degrees due to the limitation of pyvo.

"""

class Get_Data():
    def __init__(self, ra, dec):
        self.ra, self.dec = ra, dec

    def get_panstarrs_data(self):
        ra_name = str(self.ra).replace('.','_'); dec_name = str(self.dec).replace('.', '_')
        file_name = 'PS1' + '_' + 'RA'+str(ra_name) + 'DEC'+str(dec_name) + '.csv'
        TAP_service = vo.dal.TAPService("https://vao.stsci.edu/PS1DR2/tapservice.aspx")
        #TAP_service.describe()
        TAP_tables = TAP_service.tables
        for tablename in TAP_tables.keys():
            if not "tap_schema" in tablename:  
                TAP_tables[tablename].describe()
                print("Columns={}".format(sorted([k.name for k in TAP_tables[tablename].columns ])))
                print("----")
        query = """
            SELECT objID, RAMean, RAMeanErr, DecMean, DecMeanErr, gPSFMag, gPSFMagErr, gKronMag, gKronMagErr, rPSFMag, rPSFMagErr, rKronMag, rKronMagErr, iPSFMag, iPSFMagErr, iKronMag, iKronMagErr, zPSFMag, zPSFMagErr, zKronMag, zKronMagErr, yPSFMag, yPSFMagErr, yKronMag, yKronMagErr, objInfoFlag, qualityFlag, nDetections, nStackDetections, ginfoFlag, ginfoFlag2, ginfoFlag3, rinfoFlag, rinfoFlag2, rinfoFlag3, iinfoFlag, iinfoFlag2, iinfoFlag3, zinfoFlag, zinfoFlag2, zinfoFlag3, yinfoFlag, yinfoFlag2, yinfoFlag3
            FROM dbo.StackObjectView
            WHERE 
            CONTAINS(POINT('ICRS', RAMean, DecMean),CIRCLE('ICRS',{},{},{}))=1
                """.format(self.ra,self.dec,0.25)
        try:
            job = TAP_service.search(query)
            TAP_results = job.to_table()
            np.savetxt(str(file_name),\
                TAP_results, delimiter=',', header = 'objid, RAMean, RAMeanErr, DecMean, DecMeanErr, gPSFMag, gPSFMagErr, gKronMag, gKronMagErr, rPSFMag, rPSFMagErr, rKronMag, rKronMagErr, iPSFMag, iPSFMagErr, iKronMag, iKronMagErr, zPSFMag, zPSFMagErr, zKronMag, zKronMagErr, yPSFMag, yPSFMagErr, yKronMag, yKronMagErr, objInfoFlag, qualityFlag, nDetections, nStackDetections, ginfoFlag, ginfoFlag2, ginfoFlag3, rinfoFlag, rinfoFlag2, rinfoFlag3, iinfoFlag, iinfoFlag2, iinfoFlag3, zinfoFlag, zinfoFlag2,zinfoFlag3, yinfoFlag, yinfoFlag2, yinfoFlag3')
        except:
            print('This field is outside the sky coverage of PANSTARRS')
            sys.exit(0)
        return TAP_results

    def get_gaia_data(self):
        ra_name = str(self.ra).replace('.','_'); dec_name = str(self.dec).replace('.', '_')
        file_name = 'GAIA' + '_' + 'RA'+str(ra_name)\
                    + 'DEC' + str(dec_name) + '.csv'
        tables = Gaia.load_tables(only_names=True)
        #for table in (tables):
            #print (table.get_qualified_name())
        coord = SkyCoord(ra=self.ra, dec=self.dec, unit=(u.degree, u.degree), frame='icrs')
        Gaia.MAIN_GAIA_TABLE = "gaiadr3.gaia_source"
        Gaia.ROW_LIMIT = -1
        try:
            job=Gaia.cone_search(coordinate=coord, radius=u.Quantity(1.0, u.deg),\
                        table_name="gaiadr3.gaia_source",\
                        output_file=file_name, output_format='csv', verbose=True,\
                            dump_to_file=True, columns=['source_id','ra', 'ra_error,dec',\
                                                        'dec_error', 'parallax', 'parallax_error',\
                                                        'pm', 'pmra', 'pmra_error', 'pmdec',\
                                                        'pmdec_error', 'ruwe'])
        except:
            print('No Gaia observations for this field')
            sys.exit(0)
        return job.get_results()


    def get_ukidss_data(self):
        catalogs = ['UDS', 'GCS', 'GPS', 'DXS', 'LAS']
        ra_name = str(self.ra).replace('.','_'); dec_name = str(self.dec).replace('.', '_')
        file_name = 'UKIDSS' + '_' + 'RA'+str(ra_name)\
                    + 'DEC' + str(dec_name) + '.csv'
        Ukidss.filters = {'H': 4,'J': 3, 'K': 5}
        for i in range(len(catalogs)):
            print('')
            print('Name of the catalog:', str(catalogs[i]))
            try:
                table = Ukidss.query_region(SkyCoord(self.ra, self.dec, unit=(u.deg, u.deg), frame='icrs'),\
                            radius = 0.25*u.deg, programme_id = str(catalogs[i]),\
                                database='UKIDSSDR11PLUS', attributes = ['ra', 'dec', 'jPetroMag', 'jPetroMagErr', 'hPetroMag', 'hPetroMagErr', 'kPetroMag', 'kPetroMagErr'], verbose=True)

                np.savetxt(str(file_name), table, delimiter=',', header='ra, dec, jPetroMag, jPetroMagErr, hPetroMag, hPetroMagErr, kPetroMag, kPetroMagErr')# sigRa, sigDec, epoch, muRa, muDec, sigMuRa, sigMuDec, chi2, nFrames, cx, cy, cz, htmID, l, b, lambda, eta, priOrSec, ymjPnt, ymjPntErr, jmhPnt, jmhPntErr, hmkPnt, hmkPntErr, j_1mhPnt, j_1mhPntErr, j_2mhPnt, j_2mhPntErr, ymj_1Pnt, ymj_1PntErr, ymj_2Pnt, ymj_2PntErr, ymjExt, ymjExtErr, jmhExt, jmhExtErr, hmkExt, hmkExtErr, j_1mhExt, j_1mhExtErr, j_2mhExt, j_2mhExtErr, ymj_1Ext, ymj_1ExtErr, ymj_2Ext, ymj_2ExtErr, mergedClassStat, mergedClass, pStar, pGalaxy, pNoise, pSaturated, eBV, aY, aJ, aH, aK, yHallMag, yHallMagErr, yPetroMag, yPetroMagErr, yPsfMag, yPsfMagErr, ySerMag2D, ySerMag2DErr, yAperMag3, yAperMag3Err, yAperMag4, yAperMag4Err, yAperMag6, yAperMag6Err, yGausig, yEll, yPA, yErrBits, yDeblend, yClass, yClassStat, yppErrBits, ySeqNum, yXi, yEta, jHallMag, jHallMagErr, jPetroMag, jPetroMagErr, jPsfMag, jPsfMagErr, jSerMag2D, jSerMag2DErr, jAperMag3, jAperMag3Err, jAperMag4, jAperMag4Err, jAperMag6, jAperMag6Err, jGausig, jEll, jPA, jErrBits, jDeblend, jClass, jClassStat, jppErrBits, jSeqNum, jXi, jEta, j_1HallMag, j_1HallMagErr, j_1PetroMag, j_1PetroMagErr, j_1PsfMag, j_1PsfMagErr, j_1SerMag2D, j_1SerMag2DErr, j_1AperMag3, j_1AperMag3Err, j_1AperMag4, j_1AperMag4Err,	j_1AperMag6, j_1AperMag6Err,  j_1Gausig, j_1Ell, j_1PA, j_1ErrBits,	j_1Deblend,	j_1Class, j_1ClassStat,	j_1ppErrBits, j_1SeqNum, j_1Xi, j_1Eta, j_2HallMag, 	j_2HallMagErr, j_2PetroMag,	j_2PetroMagErr, j_2PsfMag, j_2PsfMagErr, j_2SerMag2D, j_2SerMag2DErr, j_2AperMag3, j_2AperMag3Err, j_2AperMag4, j_2AperMag4Err, j_2AperMag6, j_2AperMag6Err, j_2Gausig, j_2Ell, j_2PA, j_2ErrBits, j_2Deblend, j_2Class, j_2ClassStat, j_2ppErrBits, j_2SeqNum, j_2Xi, j_2Eta, hHallMag, hHallMagErr, hPetroMag, hPetroMagErr, hPsfMag, hPsfMagErr, hSerMag2D, hSerMag2DErr, hAperMag3, hAperMag3Err, hAperMag4, hAperMag4Err, hAperMag6, hAperMag6Err, hGausig, hEll, hPA, hErrBits, hDeblend, hClass, hClassStat, hppErrBits, hSeqNum, hXi, hEta, kHallMag, kHallMagErr, kPetroMag, kPetroMagErr, kPsfMag, kPsfMagErr, kSerMag2D, kSerMag2DErr, kAperMag3, kAperMag3Err, kAperMag4, kAperMag4Err, kAperMag6, kAperMag6Err, kGausig, kEll, kPA, kErrBits, kDeblend, kClass, kClassStat, kppErrBits, kSeqNum, kXi, kEta, distance')
            except:
                print('No observations in'+' '+str(catalogs[i])+' ' + 'catalog of UKIDSS')

