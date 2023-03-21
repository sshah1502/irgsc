from setuptools import setup, find_packages

setup(
    name          = 'irgsctool',
    version       = '0.1.1',    
    description   = 'Python tool to generate a catalog of NIR guide stars for the AO observations of the Thirty Meter Telescope.',
    url           = '',
    author        = 'Sarang Shah',
    author_email  = 'sshah1502@gmail.com',
    license       = 'BSD 2-clause', 
    package_dir   = {'':'src'},
    packages      = find_packages(where='src'),          
    install_requires =['astroquery', 'astropy', 'matplotlib', 'astropy', 'dustmaps',
                        'numpy', 'datetime', 'requests', 'pyvo'],
    include_package_data = True,  
    package_data       = {'': ['irgsctool/data/*']},
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: BSD License',  
        'Operating System :: OS Independent',        
        'Programming Language :: Python :: 3.6',
	    'Programming Language :: Python :: 3.9',
    ],
   
)
#import irgsctool
#from dustmaps.config import config
#config['data_dir'] = Path(__file__).parent.joinpath()

#import dustmaps.sfd
#dustmaps.sfd.fetch()
