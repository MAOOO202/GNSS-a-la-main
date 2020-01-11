# -*- coding: utf-8 -*-
"""
@author: BEILIN Jacques - IGN/ENSG
@date  : %(date)s
"""

import numpy as np
import math
import gpstime as gpst
#import sys
#sys.path.append("/home/beilin/progs/Prog_scientifique/Packages/yagnss_toolbox/python/src/pygnsstoolbox/gnsstoolbox")
#import gnss_process as proc
import gnsstoolbox.gnss_process as proc

class Trilat():
    """ Classe trilateration (ou multilateration) """
    
    def __init__(self,PosSat,Dobs,SatIndex,ElevSat,X0=np.zeros((4,1))):
        self.PosSat = PosSat
        self.Dobs = Dobs
        self.X0 = X0
        self.SatIndex = SatIndex
        self.ElevSat = ElevSat
    
    def print(self):
        """Ecriture des attributs de la classe point()"""
        for s in self.__dict__:
            print('%-35s : ' %(s), self.__dict__.get(s))

    def estimation(self):
        """ Calcul d'une trilatération simple avec un offset
        correspondant à l'erreur d'horloge du récepteur
        Et 2 offsets correspondant aux décallages Glonass et Galileo 
        et ponderation selon d'elevation """
    
        self.Qxx = [] 
        self.V = []
        self.X = 0
        self.Y = 0
        self.Z = 0
        self.cdtr = 0
        self.cGGTO = 0
        self.cGPGL = 0
        self.sigma0_2 = 0
        return 0

if __name__ == "__main__":

    tic = gpst.gpstime()

    PosSat = np.genfromtxt("possat.txt",skip_header=1,delimiter=",")  
    Dobs = np.genfromtxt("dobs.txt",skip_header=1,delimiter=",")
    SatIndex = np.genfromtxt("sat_index.txt",skip_header=1,delimiter=",")
    ElevSat = np.genfromtxt("ElevSat.txt",skip_header=1,delimiter=",")
    X0=np.array([4202000,178000,4780000,0]) 
    
    X,Y,Z,cdtr,cGGTP, cGPGL, sigma0_2,V,Qxx = proc.TrilatGnssPonderationElev(PosSat,Dobs,X0,SatIndex,ElevSat)
    print(X,Y,Z,cdtr,cGGTP, cGPGL,sigma0_2,V,Qxx)

    T = Trilat(PosSat,Dobs,SatIndex,ElevSat,X0)
    T.estimation()
    print("X = %.3f Y = %.3f Z = %.3f cdtr = %.9f" % (T.X,T.Y,T.Z,T.cdtr))
    print("cGGTO = %.9f cGPGL = %.9f" % (T.cGGTO,T.cGPGL))
    print("sigma0_2 = %.3f" % (T.sigma0_2))
    print("V = ",T.V,"\nQxx = ",T.Qxx)
#    
#    toc = gpst.gpstime()
#    print ('%.3f sec elapsed ' % (toc-tic))
