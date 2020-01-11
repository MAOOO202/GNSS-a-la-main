# -*- coding: utf-8 -*-
"""
@author: BEILIN Jacques - IGN/ENSG
@date  : %(date)s
"""

import re
import numpy as np
import math
import os

import gpstime as gpst


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

#    def pos_sat_sp3(self,const, prn, mjd, ordre):
#        """ Calcul de la postion du satellite "const/prn"
#        à un instant donné mjd """
#        
    def pos_sat_sp3(self,const,prn,mjd,ordre):
        X_ECEF, Y_ECEF, Z_ECEF, dte = 0,0,0,0

        """ A COMPLETER (TP5) """
        (orb,nl) = self.get_sp3(const,prn)
        valeur = 0
        for i,line in enumerate(orb): 
            if line[0] < mjd:
                continue
            else: 
                valeur = i
                break
                
        
        extrait = orb[valeur-int((ordre+1)/2):valeur+int((ordre+1)/2)-1]
        print("\nN° de ligne:\n",valeur,"\nExtrait:\n",extrait)

        L = [1] * ordre
    
        for j in range(len(extrait)):
            for k in range(len(extrait)): 
                if j == k:
                    l = 1
                else:
                    l = (mjd - extrait[k][0]) /(extrait[j][0]-extrait[k][0])
                    L[j]*= l
               
        print("\nL\n",L)
        
        for i in range(len(L)):
            X_ECEF += L[i]*extrait[i,1]*1000
            Y_ECEF += L[i]*extrait[i,2]*1000
            Z_ECEF += L[i]*extrait[i,3]*1000
            dte += L[i]*extrait[i,4]*10**-6

        return X_ECEF, Y_ECEF, Z_ECEF, dte
    
    def vitSat():
        return 1

if __name__ == "__main__":

    tic = gpst.gpstime()

    mysp3 = orbit_TP()

    dir_orb = '../Guerledan'
    mysp3.load_sp3(os.path.join(dir_orb,'COM20225_5M.SP3'))

    """ Definition de l'instant pour lequel on cherche une position """
    t = gpst.gpstime(yyyy=2018,doy=285,dsec=5435.000)
#    print(t)

    """ SATELLITE """
    constellation = 'G'
    prn = 14

    """ Ordre du polynome """
    ordre = 9

    """ Calcul avec la fonction integree a la toolbox """
    (X_sp3,Y_sp3,Z_sp3,dte_sp3)	 = mysp3.orb_sp3_Lagrange(constellation,prn,t.mjd,ordre)
    print("X = %.3f m\nY = %.3f m\nZ = %.3f m\ndte = %.9f s" % (X_sp3,Y_sp3,Z_sp3,dte_sp3))

#    (X1,Y1,Z1,dte1)	 = mysp3.orb_sp3_Lagrange(constellation,prn,t.mjd-0.001/86400,ordre)
#    (X2,Y2,Z2,dte2)	 = mysp3.orb_sp3_Lagrange(constellation,prn,t.mjd+0.001/86400,ordre)
#
#    print("X = %.3f Y = %.3f Z = %.3f dte = %.9f" % (X1,Y1,Z1,dte1))
#    V = np.array([[X2-X1],[Y2-Y1],[Z2-Z1]])/0.002
#    print(np.linalg.norm(V))

#    print("X = %.3f Y = %.3f Z = %.3f dte = %.9f" % (X_sp3,Y_sp3,Z_sp3,dte_sp3))
#
#    """ Calcul avec la fonction developpee lors du TP """
    X,Y,Z,dte = mysp3.pos_sat_sp3(constellation,prn,t.mjd,ordre)
    print("X = %.3f Y = %.3f Z = %.3f dte = %.9f" % (X,Y,Z,dte))

    toc = gpst.gpstime()
    print ('%.3f sec elapsed ' % (toc-tic))

