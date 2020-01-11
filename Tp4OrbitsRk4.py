# -*- coding: utf-8 -*-
"""
Created on %(date)s

@author: %(username)s
"""
import re
import numpy as np
import gpstime as gpst
import math
import pip
import os

""" Import du module orbits """
#import sys
#sys.path.append("/home/beilin/progs/Prog_scientifique/Packages/yagnss_toolbox/python/src/pygnsstoolbox/gnsstoolbox")
#import orbits

import gnsstoolbox.orbits as orbits

class orbit_TP(orbits.orbit):
    """ Classe orbits à surcharger.
    Elle hérite de la classe orbit définie dans le module orbits.py """

    def print(self):
        """ Ecriture des attributs """
        for s in self.__dict__:
            print('%-35s : ' %(s), self.__dict__.get(s))

    def pos_sat_rk4(self,const,prn,mjd):
        """ Calcul de la postion du satellite "const/prn"
        à un instant donné mjd """
        X_ECEF, Y_ECEF, ZECEF, dte = 0,0,0,0

        """ A COMPLETER (TP4) """

        return X_ECEF, Y_ECEF, ZECEF, dte


if __name__ == "__main__":

    tic = gpst.gpstime()
#    for i in pip.get_installed_distributions(local_only=True):
#        print(i)

    Orb = orbit_TP()
    
    """ lecture du fichier BRDC """
    dir_orb = '../Guerledan'
    Orb.load_rinex_n(os.path.join(dir_orb, 'BRDC00IGS_R_20182850000_01D_MN.rnx'))

    """ Definition de l'instant pour lequel on cherche une position """
    t = gpst.gpstime(yyyy=2018,doy=285,dsec=5400)
    print(t)

    """ SATELLITE """
    constellation = 'R'
    prn = 1

    try:
        Eph = Orb.get_ephemeris(constellation,prn,t.mjd)
        print(Eph)
        print("TOC : ",Eph.tgps.st_iso_epoch())
    except:
        print("Unable to find satellite !!! ")


    """ Calcul avec la fonction integree a la toolbox """
    X,Y,Z, VX,VY,VZ,dte, deb = Orb.orb_from_RK4(constellation,prn,t.mjd)
    print("X = %.3f Y = %.3f Z = %.3f det = %.9f" % (X,Y,Z,dte))

    """ Calcul avec la fonction developpee lors du TP """
    X,Y,Z,dte = Orb.pos_sat_rk4(constellation,prn,t.mjd)
    print("X = %.3f Y = %.3f Z = %.3f dte = %.9f" % (X,Y,Z,dte))

    toc = gpst.gpstime()
    print ('%.3f sec elapsed ' % (toc-tic))

