
#!/usr/bin/python
# -*- coding: iso-8859-15 -*-

# GNSS observation class
# Beilin Jacques
# 2014-06-02


import re
import math
import numpy as np
import copy
import os

import gpstime as gpst
import gnsstoolbox.rinex_o as rx
import gnsstoolbox.orbits as orbits
import gnsstoolbox.gnss_const as const

class gnss_process():
    """GNSS proccess class"""


    def __init__(self):

        """Options par defaut"""
        self.process="spp" #"DGNSS","phase"
        self.freq="C1"
        self.cut_off=3*const.d2r
        self.X0 = np.zeros((3,1))
        self.iono = 'none'
        self.nav = 'sp3'
        self.const= 'GRE'
        self.type="o"
        self.constraint=0 # contraintes sur la position, utilisé pour solution DGNSS

        self.nb_sat = 0
        self.nb_GPS = 0
        self.nb_GLO = 0
        self.nb_GAL = 0

        self.nb_GPS_calc = 0
        self.nb_GLO_calc = 0
        self.nb_GAL_calc = 0

    def spp(self,epoch,brdc):
        """Process epoch using spp"""

        self.nb_sat=len(epoch.satellites)

        c = 299792458.0

        print(epoch.satellites[0].__dict__.keys())

        mjd = epoch.tgps.mjd
        for S in epoch.satellites:

            """ Test si le satellite S est bien dans une constellation qu'on souhaite utiliser """
            if not(re.search(S.const,self.const)):
                continue

            observable = 'C1'
            S.PR = S.obs.get(observable)

            """ Test si l'observation est bien cohérente """
            if S.PR<15e6:
                continue

            """ Dispose-t-on d'un message de navigation pour ce satellite/epoch ? """
            Eph = brdc.get_ephemeris(S.const,S.PRN,mjd)
            if not(hasattr(Eph,'mjd')):
                continue

            print(S.const,S.PRN)
            print("mjd = ",epoch.tgps.mjd,observable, "=", S.obs.get(observable))

            """ A ce stage on est certain de disposer d'une observation et d'un
            message de navigation utilisable pour le satellite S à l'époque qui nous interesse"""

             """ Calcul de la date de réception (mjd) """
            S.tr = epoch.tgps.mjd
            S.tr *= 86400

            """ Calcul du temps de vol et de la date d'émission """
            temp_vol = S.PR/c
            print("temp_vol",temp_vol,"Temps de recepteur:\n",S.tr)
            S.te = S.tr - temp_vol
            print("Temps d'émission:\n",S.te)


            """ Correction de la derive d'horloge du satellite """
            alpha0 = Eph.alpha0
            alpha1 = Eph.alpha1
            alpha2 = Eph.alpha2
            TOC = Eph.TOC
            TOC *= 86400
            dte = alpha0 + alpha1*(S.te-TOC*86400)+alpha2*(S.te-TOC*86400)**2
            dtes.append(dte)
            S.te -=dte
            S.PR += c*dte

            """ Calcul de l'effet relativiste """
            e = Eph.e
            a = Eph.sqrt_a**2
            dt = (t.mjd-Eph.TOC)*(24*60*60)
            F = -4.442807633e-10
            #print("dt:\n",dt,"t:\n",S.tr,"TOE",TOC)
            #print(Eph.mjd, mjd)

            mu = 3.986005e14
            n_0 = np.sqrt(mu/a**3)
            n = n_0 + Eph.delta_n
            #print("Moyen mouvement n:\n",n,"n0:\n",n_0,"delta n:\n",Eph.delta_n)
            M = Eph.M0 + n*dt
            #print("Anomalie moyenne M:\n",M)

            e = Eph.e
            E0 = M
            #print("E0",E0)
            Ek = M + e*np.sin(E0)
            print("Ek",Ek)
            compteur = 0
            while abs (Ek - E0) > 0.01/3e7:
                E0 = Ek
                Ek = M + e*np.sin(E0)
                compteur += 1
            E = Ek
            print("E final:\n",Ek)
            #print("compteur:\n",compteur)
            dtr = f*e*np.sqrt(a)*np.sin(E)
            dtrs.append(dtr)
            S.tr +=dtr
            S.PR -= c*dtr
            """ Calcul de la position des satellites a te """
            X,Y,Z,dte = brdc.pos_sat_brdc(S.const,S.PRN,S.te)
            print("X:\n",X,"Y:\n",Y,"Z:\n",Z)
            PosSat = [X,Y,Z]

            """ Et pourtant elle tourne """

            """ Sauvegarde des donnees """
            Prcorr = S.PR
            Xs = X
            Ys = Y
            Zs = Z
            Dobs.append(Prcorr)
            PosSat.append([Xs,Ys,Zs])
            PR_cor.append(Prcorr)

        Dobs = np.array(Dobs)
        PosSat = np.array(PosSat)
        trilat = Tp7.Trilat(PosSat, PR_cor, X0)
        X_com, cdtr, V, sigma0_2, Qxx = trilat.estimation()

        print("\nDécalage de temps recepteur:\n",dtes,"\nDécalage de temps relative:\n ",
              dtrs,"\nPR corrigée:\n",PR_cor,"\nPR:\n",PR)
        return X_com

def main():

    """ lecture du fichier BRDC """
    Orb = orbits.orbit()
    dir_orb = '../Guerledan'
    Orb.load_rinex_n(os.path.join(dir_orb, 'BRDC00IGS_R_20182850000_01D_MN.rnx'))


    """ lecture du fichier SP3 """
    sp3 = orbits.orbit()
    sp3.load_sp3(os.path.join(dir_orb,'COM20225.SP3'))

    """ Definition de l'instant pour lequel on cherche une position """
    t = gpst.gpstime(yyyy=2018,doy=285,dsec=5400)
    print(t)

    rnx = rx.rinex_o()
    filename = os.path.join(dir_orb,'edf1285b.18o')
    ret = rnx.load_rinex_o(filename)

    if ret<0:
        print(ret)
        return

    Ep = rnx.get_epoch_by_mjd(t.mjd)

    spp1 = gnss_process()
    spp1.const='GRE'
    spp1.constraint=0
    d2r = np.pi / 180
    spp1.cut_off = 10 * d2r
    spp1.X0[0]=rnx.headers[0].X
    spp1.X0[1]=rnx.headers[0].Y
    spp1.X0[2]=rnx.headers[0].Z
#    spp1.spp(Ep,Orb)

    spp2 = proc.gnss_process()
    spp2.const='GRE'
    Ep2 = spp2.spp(Ep,sp3)
    print("X = %.2fm Y = %.2fm Z = %.2fm " % (Ep2.X, Ep2.Y, Ep2.Z))

if __name__ == '__main__':


    tic = gpst.gpstime()
    main()
    toc = gpst.gpstime()
    print ('%.3f sec elapsed ' % (toc-tic))
