#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Command line calculator of cosmological quantities.
For user provided redshift(s) and optionally a provided (flat) cosmology, the calculator returns

1) Distance modulus, DM 
2) Scale (kpc/arcsec), Scale
3) Luminosity Distance, DL
4) Angular Diameter Distance, DA
5) Age of the Universe at the specified redshifts, Age
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

"""

import argparse
import importlib
import numpy as np
from astropy.table import Table, Column


cosmo_models = ["Planck15", "Planck13", "WMAP9", "WMAP7", "WMAP5", "FlatLambdaCDM"]


def get_arguments():
    parser = argparse.ArgumentParser(description="Command line calculator of cosmological quantities",
                                    formatter_class=argparse.RawDescriptionHelpFormatter)

    parser = argparse.ArgumentParser(prog="cosmo.py",
                                    formatter_class=argparse.RawDescriptionHelpFormatter,
                                    epilog=__doc__)

    parser.add_argument('redshifts', type=float, nargs='+', 
                        help='list of redshifts')

    parser.add_argument("-c", type=str, default="FlatLambdaCDM",
                        metavar="COSMO", dest="cosmo",
                        choices=cosmo_models, 
                        help="Specify the cosmological model")

    parser.add_argument("-Om", type=float, default=0.3, 
                        metavar="Omatter", dest="Om",
                        help="Matter density at z=0")

    parser.add_argument("-H0", type=float, default=70.0, 
                        metavar="H0", dest="H0",
                        help="Hubble constant")

    args = parser.parse_args()

    return args


def get_cosmo_model(cosmo_model="FlatLambdaCDM", omega_matter=0.3, hubble_const=70 ):

    if omega_matter<=0 or omega_matter>1:
        raise ValueError("Invalid value for Om(0)")
    if hubble_const<=0 or hubble_const>=200:
        raise ValueError("Invalid value for H0")
        
    print("Using", cosmo_model)

    if cosmo_model != "FlatLambdaCDM":
        cosmo = getattr(importlib.import_module('astropy.cosmology', package=cosmo_model), cosmo_model)

    else:
        cosmo = getattr(importlib.import_module('astropy.cosmology', package=cosmo_model), cosmo_model)
        cosmo = cosmo(H0=hubble_const, Om0=omega_matter)

    params = "with H0=" + str(cosmo.H(0)) + " and Om(0)=" + str(cosmo.Om0)
    print(params)  
    print("-"*len(params))
    
    return cosmo


def check_z(z):
    z = np.array(z)
    if (z<0).any():
        raise ValueError("redshifts must be possitive")
    else:
        return z


def create_table(z, cosmo):

    redshift = Column(data=z, name="z")
    DM    = Column(data=np.round(cosmo.distmod(z)), name="DM")
    scale = Column(data=np.round(1/cosmo.arcsec_per_kpc_proper(z),3), name="Scale")
    LD    = Column(data=np.round(cosmo.luminosity_distance(z),2), name="DL")
    Age   = Column(data=np.round(cosmo.age(z),2), name="Age")
    tL    = Column(data=np.round(cosmo.lookback_time(z),2), name="tL")

    table = Table([redshift, DM, scale, LD, Age, tL])

    return table


def main():

    args = get_arguments()

    z = args.redshifts
    cosmo_model = args.cosmo
    Om = args.Om
    H0 = args.H0

    z = check_z(z)
    cosmo = get_cosmo_model(cosmo_model, omega_matter=Om, hubble_const=H0)

    table_out = create_table(z, cosmo)

    print(table_out)


if __name__ == "__main__":
    main()
