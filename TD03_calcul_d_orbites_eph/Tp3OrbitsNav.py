# -*- coding: utf-8 -*-
"""
Created on %(date)s

@author: %(username)s
"""
import re
import os
import numpy as np
import gpstime as gpst
import math as m
import pip
import math

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


    def pos_sat_brdc(self,const,prn,mjd):
        """ Calcul de la postion du satellite "const/prn"
        à un instant donné mjd """
        X_ECEF, Y_ECEF, Z_ECEF, dte = 0,0,0,0

        """ A COMPLETER (TP3) """
        Eph = self.get_ephemeris(const,prn,mjd)
        
        # Etape 1
        t = gpst.gpstime(mjd=mjd)
        t = t.mjd
        TOC = Eph.TOC
        dt = t-TOC
        print("dt:\n",dt,"t:\n",t,"TOE",TOC)
        print("demi-grand axe a=%.3fm" % (Eph.sqrt_a**2))
        print(Eph.mjd, mjd)
        a = Eph.sqrt_a**2
        mu = 3.986005e14
        n_0 = np.sqrt(mu/a**3)
        n = n_0 + Eph.delta_n
        print("Moyen mouvement n:\n",n,"n0:\n",n_0,"delta n:\n",Eph.delta_n)
        M = Eph.M0 + n*dt
        print("Anomalie moyenne M:\n",M)
        
        e = Eph.e
        E0 = M
        print("E0",E0)
        Ek = M + e*np.sin(E0)
        print("Ek",Ek)
        compteur = 0
        while abs (Ek - E0) > 0.01/3e7:
            E0 = Ek
            Ek = M + e*np.sin(E0)
            compteur += 1
        E = Ek
        print("E final:\n",Ek)
        print("compteur:\n",compteur)
        v = 2* m.atan(np.sqrt((1+e)/(1-e))*np.tan(E/2))
        print("Anomalie vraie:\n",v)
        r = a * (1-e*np.cos(E))
        print("Rayon Terre-Satellite:\n",r)
        omega = Eph.omega
        cus = Eph.cus
        cuc = Eph.cuc
        crs = Eph.crs
        crc = Eph.crc
        phi = omega + v
        print("phi:\n",phi)
        dphi = cus * np.sin(2*phi) + cuc * np.cos(2*phi)
        dr =  crs* np.sin(2*phi) + crc * np.cos(2*phi)
        x_orb = (r+dr) * np.cos(phi + dphi)
        y_orb = (r+dr) * np.sin(phi + dphi)
        print("x_orb:\n",x_orb,"\ny_orb:\n",y_orb)
        # Etape 2
        OMEGA_0 = Eph.OMEGA
        i_DOT = Eph.IDOT
        OMEGA_DOT = Eph.OMEGA_DOT
        i0 = Eph.i0
        cis = Eph.cis
        cic = Eph.cic
        di = cis * np.sin(2*phi) + cic*np.cos(2*phi)
        i = i0 + i_DOT*dt + di
        OMEGA = OMEGA_0 + OMEGA_DOT*dt
        print("Inclinaison de l'orbite i:\n",i,"\nAscension de l'orbite w",OMEGA)
        A = np.array([[np.cos(OMEGA), np.sin(OMEGA),0.0],
                       [np.sin(OMEGA),np.cos(OMEGA),0.0],
                       [0.0,0.0,1.0]])
        B = np.array([[1.0,0.0,0.0],
                      [0.0,np.cos(i),np.sin(i)],
                      [0.0,np.sin(i), np.cos(i)]])
        C = np.array([[x_orb],[y_orb],[0.0]])
        N = np.dot(A,B)
        coor = np.dot(N,C)
        print("Résultats",coor)
        
        D = np.array([[np.cos(OMEGA_DOT*t),np.sin(OMEGA_DOT*t),0],
                       [np.sin(OMEGA_DOT*t),np.cos(OMEGA_DOT*t),0],
                       [0,0,1]])
        E = np.dot(D,coor)
        X_ECEF = E[0]
        Y_ECEF = E[1]
        Z_ECEF = E[2]
        return X_ECEF, Y_ECEF, Z_ECEF, dte
    
        
        # Tp 2

    
    


if __name__ == "__main__":

    tic = gpst.gpstime()
#    for i in pip.get_installed_distributions(local_only=True):
#        print(i)

    Orb = orbit_TP()

    """ lecture du fichier BRDC """
    dir_orb = '../Guerledan'
    mysp3 = orbit_TP()
    
    Orb.load_rinex_n(os.path.join(dir_orb, 'BRDC00IGS_R_20182850000_01D_MN.rnx'))

    """ Definition de l'instant pour lequel on cherche une position """
    t = gpst.gpstime(yyyy=2018,doy=285,dsec=5435)
    print(t.mjd)

    """ SATELLITE """
    constellation = 'G'
    prn = 14
    ordre = 9

#    try:
        #Eph = Orb.get_ephemeris(constellation,prn,t.mjd)
        #print("Eph\n",Eph)
        #print("TOC : ",Eph.tgps.st_iso_epoch())

             
#    except:
#        print("Unable to find satellite !!! ")


    #print(Eph.sqrt_a**2,Eph.e)

    """ Calcul avec la fonction integee a la toolbox """
    X,Y,Z,dte = Orb.orb_from_eph(constellation,prn,t.mjd) 
    print("X = %.3f m\nY = %.3f m\nZ = %.3f m\ndte = %.9f s" % (X,Y,Z,dte))
#    a= Orb.orb_from_eph(constellation,prn,t.mjd)
#    print(a, type(a))
#    print(Orb.debug.__dict__)
#
#    print(">>>>>>>>>>>>>",Orb.debug.X_ECI[0].squeeze())

    """ Calcul avec la fonction integee developpee lors du TP """
    X,Y,Z,dte = Orb.pos_sat_brdc(constellation,prn,t.mjd)
    print("X = %.3f Y = %.3f Z = %.3f dte = %.9f" % (X,Y,Z,dte))


    

#
###    print(Orb.debug.__dict__)

    toc = gpst.gpstime()
    print ('%.3f sec elapsed ' % (toc-tic))

