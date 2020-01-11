# -*- coding: utf-8 -*-
"""
@author: BEILIN Jacques - IGN/ENSG
@date  : %(date)s
"""

import re
import numpy as np
import math
from numpy.linalg import inv

import gpstime as gpst

#import sys
#sys.path.append("/home/beilin/progs/Prog_scientifique/Packages/yagnss_toolbox/python/src/pygnsstoolbox/gnsstoolbox")
#import gnss_process as proc

import gnsstoolbox.gnss_process as proc

class Trilat():
    """ Classe trilateration (ou multilateration) """

    def __init__(self,PosSat,Dobs,ElevSat,X0=np.zeros((4,1))):
        self.PosSat = PosSat
        self.Dobs = Dobs
        self.X0 = X0
        self.ElevSat = ElevSat

    def print(self):
        """Ecriture des attributs de la classe point()"""
        for s in self.__dict__:
            print('%-35s : ' %(s), self.__dict__.get(s))


    def estimation(self):
        """ Calcul d'une trilatération simple avec un offset
        correspondant à l'erreur d'horloge du récepteur
        Et 2 offsets correspondant aux décallages Glonass et Galileo """
    
        #print("PosSat = ",self.PosSat)
        Nb = len(self.PosSat)
        c = 299792458.0
        P_obs = np.eye(Nb)
        for i in range(Nb):
            P_obs[i][i] = 2/np.sin(self.ElevSat[i])
        print("P_obs:",P_obs)
        p_approche = []
        for i in range(Nb):
            dist = np.sqrt((self.X0[0]-self.PosSat[i][0])**2+
                           (self.X0[1]-self.PosSat[i][1])**2+
                           (self.X0[2]-self.PosSat[i][2])**2)
            p_approche.append(dist)
            #print(p_approche)
        nb = len(X0)
        B = np.zeros(Nb)
        n = Nb-nb
            #print("\nn\n",n,"\nNb\n",Nb,"\nbn\n",nb)
        for i in range(Nb):
            B[i] = self.Dobs[i] - p_approche[i]
            #print("\nB:\n",B)
        A = np.zeros((Nb,4))
        A[:,0] = (self.X0[0]-self.PosSat[:,0])/p_approche
        A[:,1] = (self.X0[1]-self.PosSat[:,1])/p_approche
        A[:,2] = (self.X0[2]-self.PosSat[:,2])/p_approche
        A[:,3] = c
            #print("\nA:\n",A)
        a = np.dot(A.T,P_obs)
            
        N = np.dot(a,A)
        Ninv = inv(N)
        C = np.dot(a,B)
        X_sol = np.dot(Ninv,C)
            #print("\nSolution:\n",X_sol)
            #print("\nv:\n",V)
        self.X_com = X0+X_sol
            #print("\nX compensé:\n",X)
      
                
        self.Qxx = Ninv
        self.V = B - np.dot(A,X_sol)
        self.cdtr = c*X_sol[3]
        self.sigma0_2 = np.dot(V.T,V)/n
        self.X = self.X_com[0]
        self.Y = self.X_com[1]
        self.Z = self.X_com[2]
        self.cGGTO = 0
        self.cGPGL = 0
        self.sigma0_2 = 0
            
        return self.X,self.Y,self.Z,self.cdtr,self.V,self.sigma0_2,self.Qxx,self.X_com

if __name__ == "__main__":

    tic = gpst.gpstime()

    PosSat = np.genfromtxt("possat.txt",skip_header=1,delimiter=",")
    Dobs = np.genfromtxt("dobs.txt",skip_header=1,delimiter=",")
    ElevSat = np.genfromtxt("ElevSat.txt",skip_header=1,delimiter=",")
    X0=np.array([4202000,178000,4780000,0])
    satIndex = np.ones(Dobs.shape)

    X,Y,Z,cdtr,cGGTP, cGPGL, sigma0_2,V,SigmaX = proc.TrilatGnssPonderationElev(PosSat,Dobs,X0,satIndex,ElevSat)
    print(X,Y,Z,cdtr,sigma0_2,V)
    print("SigmaX", SigmaX)
    
    T = Trilat(PosSat,Dobs,ElevSat,X0)
    compteur = 1
    X,Y,Z,cdtr,V,sigma0_2,Qxx,X_com = T.estimation()
#    for i in range(10):
#        if all(abs(X_com-X0)) >0.001:
#            X0 = X_com
#            T = Trilat(PosSat,Dobs,X0)
#            X,Y,Z,cdtr,V,sigma0_2,Qxx = T.estimation()
#            compteur +=1
#            print(compteur)
#        else :
#            print("Les coordonnées compensées:\n",X_com[:3],
#                  "\nDistance du décalage de temps: \n",cdtr,
#                  "\nRésidus:\n",V,
#                  "\nEstimateur 0 carré:\n",sigma0_2)
#            break

    print("X = %.3f Y = %.3f Z = %.3f cdtr = %.9f" % (T.X,T.Y,T.Z,T.cdtr))
    print("cGGTO = %.9f cGPGL = %.9f" % (T.cGGTO,T.cGPGL))
    print("sigma0_2 = %.3f" % (T.sigma0_2))
    print("V = ",T.V,"\nQxx = ",T.Qxx)

    toc = gpst.gpstime()
    print ('%.3f sec elapsed ' % (toc-tic))
