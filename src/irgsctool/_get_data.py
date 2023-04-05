import sys
import numpy as np
import pyvo as vo
from astroquery.ukidss import Ukidss
from astroquery.gaia import Gaia
from astropy.coordinates import SkyCoord
import astropy.units as u

header = ['objid, RAMean, RAMeanErr, DecMean, DecMeanErr, gPSFMag, \
    gPSFMagErr, gKronMag, gKronMagErr, rPSFMag, rPSFMagErr, rKronMag, \
        rKronMagErr, iPSFMag, iPSFMagErr, iKronMag, iKronMagErr, zPSFMag, \
            zPSFMagErr, zKronMag, zKronMagErr, yPSFMag, yPSFMagErr, yKronMag, \
                yKronMagErr, objInfoFlag, qualityFlag, nDetections, nStackDetections, \
                    ginfoFlag, ginfoFlag2, ginfoFlag3, rinfoFlag, rinfoFlag2, rinfoFlag3, \
                        iinfoFlag, iinfoFlag2, iinfoFlag3, zinfoFlag, zinfoFlag2,zinfoFlag3, \
                            yinfoFlag, yinfoFlag2, yinfoFlag3']

class GetData():
    """
    <justify> *** GetData class *** contains methods to obtain 
    PANSTARRS DR2 stacked optical photometry
    data, UKIDSS NIR observed data and GAIA DR3 astrometry data.
    The default search radius is 0.25 degrees due to the
    limitation of pyvo.
    The data retrieved is stored in .csv format with the name
    of the survey + str(ra) + str(dec) </justify>

    Examples:
            >>> irgsctool.GetData.get_panstarrs_data(0.0,0.0)
            'PS1_RA_0_0_DEC_0_0.csv'
            >>> irgsctool.GetData.get_gaia_data(0.0,0.0)
            'GAIA_RA_0_0_DEC_0_0.csv'
            >>> irgsctool.GetData.get_ukidss_data(0.0,0.0)
            'UKIDSS_RA_0_0_DEC_0_0.csv'

    """

    def __init__(self, ra, dec):
        self.ra, self.dec = ra, dec

    def get_panstarrs_data(self):
        """
        `irgsctool.GetData.get_panstarrs_data()`

        <justify> This function sends a query to obtain the PANSTARRS 
        data from DR2 database.
        The query uses pyvo TAPService module for retrieving the data.
        The data is selected from StackObjectView db and the maximum
        search radius is 0.25 degrees.</justify>
        Raises:
                ValueError if there is no observed PANSTARRS DR2 data for the
                            given set of input coordinates.
        """
        ra_name = str(self.ra).replace('.','_')
        dec_name = str(self.dec).replace('.', '_')
        file_name = 'PS1'+'_'+'RA'+str(ra_name)+'DEC'+str(dec_name)+'.csv'
        Tap_Service = vo.dal.TAPService("https://vao.stsci.edu/PS1DR2/tapservice.aspx")
        Tap_Service.describe()
        Tap_Tables = Tap_Service.tables
        print('table keys=', Tap_Tables.keys())
        for tablename in Tap_Tables.keys():
            if not "TAP_schema" in tablename:  
                Tap_Tables[tablename].describe()
                print("Columns={}".format(sorted([k.name for k in\
                                                  Tap_Tables[tablename].columns ])))
                print("----")
        query = """
            SELECT objID, RAMean, RAMeanErr, DecMean, DecMeanErr, gPSFMag, gPSFMagErr, gKronMag, gKronMagErr, rPSFMag, rPSFMagErr, rKronMag, rKronMagErr, iPSFMag, iPSFMagErr, iKronMag, iKronMagErr, zPSFMag, zPSFMagErr, zKronMag, zKronMagErr, yPSFMag, yPSFMagErr, yKronMag, yKronMagErr, objInfoFlag, qualityFlag, nDetections, nStackDetections, ginfoFlag, ginfoFlag2, ginfoFlag3, rinfoFlag, rinfoFlag2, rinfoFlag3, iinfoFlag, iinfoFlag2, iinfoFlag3, zinfoFlag, zinfoFlag2, zinfoFlag3, yinfoFlag, yinfoFlag2, yinfoFlag3
            FROM dbo.StackObjectView
            WHERE 
            CONTAINS(POINT('ICRS', RAMean, DecMean),CIRCLE('ICRS',{},{},{}))=1
            AND
            primaryDetection>0
                """.format(self.ra,self.dec,0.25)

        try:
            #job = Tap_Service.search(query)
            job = Tap_Service.submit_job(query)
            job.run()
            job_url = job.url
            job = vo.dal.tap.AsyncTAPJob(job_url)
            #Tap_Results = job.to_table()
            Tap_Results = job.fetch_result()
            table = Tap_Results.table
            print(Tap_Results)
            #with open("result.csv") as f:
                #Tap_Results.to_table().write(output=f, header = 'objid, RAMean, RAMeanErr, DecMean, DecMeanErr, gPSFMag, gPSFMagErr, gKronMag, gKronMagErr, rPSFMag, rPSFMagErr, rKronMag, rKronMagErr, iPSFMag, iPSFMagErr, iKronMag, iKronMagErr, zPSFMag, zPSFMagErr, zKronMag, zKronMagErr, yPSFMag, yPSFMagErr, yKronMag, yKronMagErr, objInfoFlag, qualityFlag, nDetections, nStackDetections, ginfoFlag, ginfoFlag2, ginfoFlag3, rinfoFlag, rinfoFlag2, rinfoFlag3, iinfoFlag, iinfoFlag2, iinfoFlag3, zinfoFlag, zinfoFlag2,zinfoFlag3, yinfoFlag, yinfoFlag2, yinfoFlag3', format="csv")
            np.savetxt(str(file_name),\
                table, delimiter=',', header = 'objid, RAMean, RAMeanErr, DecMean, DecMeanErr, gPSFMag, gPSFMagErr, gKronMag, gKronMagErr, rPSFMag, rPSFMagErr, rKronMag, rKronMagErr, iPSFMag, iPSFMagErr, iKronMag, iKronMagErr, zPSFMag, zPSFMagErr, zKronMag, zKronMagErr, yPSFMag, yPSFMagErr, yKronMag, yKronMagErr, objInfoFlag, qualityFlag, nDetections, nStackDetections, ginfoFlag, ginfoFlag2, ginfoFlag3, rinfoFlag, rinfoFlag2, rinfoFlag3, iinfoFlag, iinfoFlag2, iinfoFlag3, zinfoFlag, zinfoFlag2,zinfoFlag3, yinfoFlag, yinfoFlag2, yinfoFlag3', fmt="%s")
        except Exception:
            raise ValueError('This field is outside the sky coverage of PANSTARRS')
        return Tap_Results

    def get_gaia_data(self):
        """
        `irgsctool.GetData.get_gaia_data()`

        <justify> This function sends a query to obtain GAIA DR3 data
        using the astroquery module.
        The ROW_LIMIT is set to -1 which implies that the query retriees all
        the rows in the given field.</justify>
        Raises:
                ValueError if there is no observed GAIA DR3 data for the
                            given set of input coordinates.
        """
        ra_name = str(self.ra).replace('.','_')
        dec_name = str(self.dec).replace('.', '_')
        file_name = 'GAIA' + '_' + 'RA'+str(ra_name)\
                    + 'DEC' + str(dec_name) + '.csv'
        #tables = Gaia.load_tables(only_names=True)
        #for table in (tables):
            #print (table.get_qualified_name())
        coord = SkyCoord(ra=self.ra, dec=self.dec, unit=(u.degree, u.degree), frame='icrs')
        Gaia.MAIN_GAIA_TABLE = "gaiadr3.gaia_source"
        Gaia.ROW_LIMIT = -1
        try:
            job=Gaia.cone_search(coordinate=coord, radius=u.Quantity(0.25, u.deg),\
                        table_name="gaiadr3.gaia_source",\
                        output_file=file_name, output_format='csv', verbose=True,\
                            dump_to_file=True, columns=['source_id','ra', 'ra_error','dec',\
                                                        'dec_error', 'parallax', 'parallax_error',\
                                                        'pm', 'pmra', 'pmra_error', 'pmdec',\
                                                        'pmdec_error', 'ruwe'])
        except Exception:
            raise ValueError('No Gaia observations for this field')
        return job.get_results()


    def get_ukidss_data(self):
        """
        `irgsctool.GetData.get_ukidss_data()`

        <justify> This function sends a query to obtain UKIDSS DR11
        NIR data using astroquery.
        UKIDSS consists of five sub-surveys viz. UDS, GPS, GCS,
        DXS and LAS. The query loops over this surveys and retrieves
        the data for the given coordinates.
        The surveys which do not contain J and H band data,
        the function sends an alert. <\justify>
        Raises:
                ValueError if there is no observed UKIDSS data for the
                            given set of input coordinates.


        """
        catalogs = ['UDS', 'GCS', 'GPS', 'DXS', 'LAS']
        ra_name = str(self.ra).replace('.','_')
        dec_name = str(self.dec).replace('.', '_')
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

                table.write(file_name, format = 'csv', overwrite=True)
                if len(table)>8.0:
                    break
            except Exception:
                raise ValueError('No observations in'+' '+str(catalogs[i])+' ' + 'catalog of UKIDSS')
        return table

