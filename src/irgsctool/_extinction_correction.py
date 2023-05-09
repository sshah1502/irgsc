#pylint: disable=wrong-import-position
#pylint: disable=import-error
import numpy as np
from astropy.coordinates import SkyCoord
from astropy import units as u
from dustmaps.sfd import SFDQuery
import dustmaps
from ._sgc import StarGalaxyClassification as sgc

class ExtinctionCorrection():
        """
                        *** ExtinctionCorrection class *** has two methods; one to obtain the reddening and 
                        NIR extinction coefficients, while the other to correct the PANSTARRS 
                        data for extinction in each optical filter.
        """
        def __init__(self, ra, dec):
                self.ra, self.dec = ra, dec
                self.sgc = sgc(ra,dec)

        def get_reddening(self):
                """
                        `irgsctool.ExtinctionCorrection.get_reddening()`
                
                        This method obtains the reddening value for a given set of 
                        input coordinates from Schelgel et.al. 1998 (sfd) reddening
                        map. irgsctool uses Schlafly & Finkbeiner 2011 (snf) reddening
                        map which is snf = 0.86*sfd.

                        Raises:
                                FileNotFoundError: if the sfd files are not present.
                        
                        Returns:
                                ebv: Reddening from Schlafly & Finkbeiner 2011
                                err_ebv: Uncertinty in reddening
                                aj: J-band extinction coefficient.
                                ah: H-band extinction coefficient.
                                ak: K-band extinction coefficient.

                """
                try:
                        coords = SkyCoord((self.ra)*u.degree, (self.dec)\
                                          *u.degree, frame='icrs')
                        sfd = SFDQuery()
                        #sfd reddening
                        sfd_ebv = sfd(coords)
                        #s&f reddening
                        snf_ebv = sfd_ebv*0.86
                        err_snf_ebv = 0.001*snf_ebv
                        aj = 0.709*snf_ebv
                        ah = 0.449*snf_ebv
                        ak = 0.302*snf_ebv
                except FileNotFoundError:
                        dustmaps.sfd.fetch()
                        coords = SkyCoord((self.ra)*u.degree, (self.dec)\
                                          *u.degree, frame='icrs')
                        sfd = SFDQuery()
                        #sfd reddening
                        sfd_ebv = sfd(coords)
                        #s&f reddening
                        snf_ebv = sfd_ebv*0.86
                        err_snf_ebv = 0.001*snf_ebv
                        aj = 0.709*snf_ebv
                        ah = 0.449*snf_ebv
                        ak = 0.302*snf_ebv
                return snf_ebv, err_snf_ebv, aj,ah,ak

        def extinction_corrected_photometry(self):
                """
                        `irgsctool.ExtinctionCorrection.extinction_correctdd_photometry()`

                        This method corrects the input optical PANSTARRS data for reddening
                        and extinction along the line of site.

                        Returns:
                                ndarray: Extinction corrected PANSTARRS optical photometry.

                """
                print("########################################")
                print("")
                print('Correcting the optical photometry of the probable stellar sources for extinction.')
                print("")
                print("########################################")
                print("")

                ebv, err_ebv,_,_,_ = self.get_reddening()
                print('ebv=', ebv)
                ps_phot = self.sgc.star_galaxy_classification()
                print('')
                print('Length of PS1 data before ec is:', len(ps_phot[0]))
                ps1_objid, ps_ra, err_ps_ra, ps_dec, err_ps_dec, gpsf, e_gpsf, gkron, e_gkron,\
                rpsf, e_rpsf, rkron, e_rkron, ipsf, e_ipsf, ikron, e_ikron, zpsf, e_zpsf,\
                zkron, e_zkron, ypsf, e_ypsf, ykron, e_ykron, objinfoflag, qualityflag,\
                ndetections, nstackdetections, ginfoflag, ginfoflag2, ginfoflag3, rinfoflag,\
                rinfoflag2, rinfoflag3, iinfoflag,iinfoflag2, iinfoflag3, zinfoflag, zinfoflag2,\
                zinfoflag3, yinfoflag, yinfoflag2, yinfoflag3 = ps_phot

                #extinction in ps1 filters taken from Tonry et.al. 2012
                ag = (ebv)*0.88*(3.613 - 0.0972*(gpsf - ipsf) + 0.01*(gpsf - ipsf)**2)
                ar = (ebv)*0.88*(2.585 - 0.0315*(gpsf - ipsf))
                ai = (ebv)*0.88*(1.908 - 0.0152*(gpsf - ipsf))
                az = (ebv)*0.88*(1.499 - 0.0023*(gpsf - ipsf))
                ay = (ebv)*0.88*(1.251 - 0.0027*(gpsf - ipsf))

                #error in extinction
                e_gi = np.sqrt(e_gpsf**2 + e_ipsf**2)
                e_ag = ((err_ebv)*(3.613 - 0.0972*(gpsf - ipsf) + 0.02*(gpsf - ipsf)**2))+\
                        ((ebv)*((-0.0972*e_gi)+(0.02*(gpsf - ipsf)*e_gi)))
                e_ar = (err_ebv)*(2.585 - 0.0315*(gpsf - ipsf)) + (ebv)*(-0.0315*e_gi)
                e_ai = (err_ebv)*(1.908 - 0.0152*(gpsf - ipsf)) + (ebv)*(-0.0152*e_gi)
                e_az = (err_ebv)*(1.499 - 0.0023*(gpsf - ipsf)) + (ebv)*(-0.0023*e_gi)
                e_ay = (err_ebv)*(1.251 - 0.0027*(gpsf - ipsf)) + (ebv)*(-0.0027*e_gi)

                #"ec_" stands for extinction corrected magnitudes

                ec_gmag = gpsf - ag
                ec_rmag = rpsf - ar
                ec_imag = ipsf - ai
                ec_zmag = zpsf - az
                ec_ymag = ypsf - ay

                #e_ec_ stands for error in extinction corrected magnitudes

                e_ec_gmag = np.sqrt((e_gpsf)**2 + (e_ag)**2)
                e_ec_rmag = np.sqrt((e_rpsf)**2 + (e_ar)**2)
                e_ec_imag = np.sqrt((e_ipsf)**2 + (e_ai)**2)
                e_ec_zmag = np.sqrt((e_zpsf)**2 + (e_az)**2)
                e_ec_ymag = np.sqrt((e_ypsf)**2 + (e_ay)**2)

                psf_phot = ps1_objid,ps_ra,err_ps_ra,ps_dec,err_ps_dec,\
                        ec_gmag,e_ec_gmag,gkron,e_gkron,ec_rmag,e_ec_rmag,\
                        rkron, e_rkron, ec_imag, e_ec_imag, ikron, e_ikron,\
                        ec_zmag,e_ec_zmag,zkron,e_zkron,ec_ymag,e_ec_ymag,\
                        ykron,e_ykron,objinfoflag,qualityflag,ndetections,\
                        nstackdetections,ginfoflag,ginfoflag2,ginfoflag3,\
                        rinfoflag,rinfoflag2,rinfoflag3,iinfoflag,iinfoflag2,\
                        iinfoflag3,zinfoflag,zinfoflag2,zinfoflag3,yinfoflag,\
                        yinfoflag2, yinfoflag3
                print('length of ps1 data after ec is:', len(psf_phot[0]))
                return psf_phot
