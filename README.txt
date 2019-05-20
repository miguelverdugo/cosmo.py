README
======

cosmo.py: Simple and handy command line calculator of cosmological quantities. It relies in astropy.cosmology module for the actual calculations

usage: cosmo.py [-h] [-c COSMO] [-Om Omatter] [-H0 H0]
                redshifts [redshifts ...]

positional arguments:
  redshifts    list of redshifts

optional arguments:
  -h, --help   show this help message and exit
  -c COSMO     Specify the cosmological model
  -Om Omatter  Matter density at z=0
  -H0 H0       Hubble constant

Command line calculator of cosmological quantities.
For user provided redshift(s) and optionally a provided (flat) cosmology, the calculator returns

1) Distance modulus, DM 
2) Scale (kpc/arcsec), Scale
3) Luminosity Distance, DL
4) Angular Diameter Distance, DA
5) Age of the Universe, Age
6) Lookback time, tL 

Built-in cosmologies
--------------------

Name          Source                         H0      Om     Flat
WMAP5         Komatsu et al. 2009            70.2    0.277   Yes
WMAP7         Komatsu et al. 2011            70.4    0.272   Yes
WMAP9         Hinshaw et al. 2013            69.3    0.287   Yes
Planck13      Planck Collab 2013, Paper XVI  67.8    0.307   Yes
Planck15      Planck Collab 2015, Paper XIII 67.7    0.307   Yes
FlatLambdaCDM        (default)               70.0    0.3     Yes

Examples 
--------
cosmo.py 1.5 
cosmo.py 0.22 -c WMAP5 
cosmo.py 0.3 0.4 0.5 -Om 0.25 -H0 100 


