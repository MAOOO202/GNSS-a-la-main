
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

import gnsstoolbox.gnsstools as tools
import gpstime as gpst
import gnsstoolbox.rinex_o as rx
import gnsstoolbox.orbits as orbits
import gnsstoolbox.gnss_const as const
import gnsstoolbox.gnss_process as proc
sys.path.append('../TD10_trilat_GPS')
import Tp10_cqlcGPS as Tp10

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
            temp_vol = S.PR / c
            print("temp_vol", temp_vol, "Temps de recepteur:\n", S.tr)
            S.te = S.tr - temp_vol
            print("Temps d'émission:\n", S.te)

            """ Correction de la derive d'horloge du satellite """
            alpha0 = Eph.alpha0
            alpha1 = Eph.alpha1
            alpha2 = Eph.alpha2
            TOC = Eph.TOC
            TOC *= 86400
            dte = alpha0 + alpha1 * (S.te - TOC * 86400) + alpha2 * (S.te - TOC * 86400) ** 2
            dtes.append(dte)
            S.te -= dte
            S.PR += c * dte

            """ Calcul de l'effet relativiste """
            e = Eph.e
            a = Eph.sqrt_a ** 2
            dt = (t.mjd - Eph.TOC) * (24 * 60 * 60)
            f = -4.442807633e-10
            # print("dt:\n",dt,"t:\n",S.tr,"TOE",TOC)
            # print(Eph.mjd, mjd)

            mu = 3.986005e14
            n_0 = np.sqrt(mu / a ** 3)
            n = n_0 + Eph.delta_n
            # print("Moyen mouvement n:\n",n,"n0:\n",n_0,"delta n:\n",Eph.delta_n)
            M = Eph.M0 + n * dt
            # print("Anomalie moyenne M:\n",M)

            e = Eph.e
            E0 = M
            # print("E0",E0)
            Ek = M + e * np.sin(E0)
            print("Ek", Ek)
            compteur = 0
            while abs(Ek - E0) > 0.01 / 3e7:
                E0 = Ek
                Ek = M + e * np.sin(E0)
                compteur += 1
            E = Ek
            print("E final:\n", Ek)
            # print("compteur:\n",compteur)
            dtr = f * e * np.sqrt(a) * np.sin(E)
            dtrs.append(dtr)
            S.tr += dtr
            S.PR -= c * dtr

            """ Calcul de la position des satellites a te """
            X, Y, Z, dte = brdc.pos_sat_brdc(S.const, S.PRN, S.te)
            print("X:\n", X, "Y:\n", Y, "Z:\n", Z)


            """ Et pourtant elle tourne """

            """ Sauvegarde des donnees """
            Prcorr = S.PR
            Xs = X
            Ys = Y
            Zs = Z
            Dobs.append(Prcorr)
            PosSat.append([Xs, Ys, Zs])
            PR_cor.append(Prcorr)






        Dobs = np.array(Dobs)
        PosSat = np.array(PosSat)
        trilat = Tp7.Trilat(PosSat, Dobs, X0)
        X_com, cdtr, V, sigma0_2, Qxx = trilat.estimation()

        print("\nDécalage de temps recepteur:\n", dtes, "\nDécalage de temps relative:\n ",
              dtrs, "\nPR corrigée:\n", PR_cor, "\nPR:\n", S.PR)




        """ Estimation des elevations des satellites """
        elevations = []
        azimuts = []
        elevations_ref = []
        azimuts_ref = []
        hs_ref = []
        #[ATX_header, ATX_data] = load_antex('igs08.atx')
        #ATX_name = RNX_header.ANT_TYPE
        #[h] = get_antex(ATX_header, ATX_data, ATX_name, 'G01')
        for i in range(len(PosSat_local)):
            p = np.sqrt((PosSat_local[i][0]-X_com[0])**2 + (PosSat_local[i][1] - X_com[1])**2 + (PosSat_local[i][2] - X_com[2])**2)
            elevation = np.arcsin(PosSat_local[i][2]/p)
            azimut = np.arctan(PosSat_local[i][1]/PosSat_local[i][0])
            # hauteur =
            elevations.append(elevation)
            azimuts.append(azimut)

            [az, ele, h] = tools.tool_az_ele_h(X_com[i][0], X_com[i][1], X_com[i][2], PosSat_local[i][0],
                                               PosSat_local[i][1], PosSat_local[i][2])
            azimuts_ref.append(az)
            elevations_ref.append(ele)
            hs_ref.append(h)
        print(elevations)

        for j in elevations:           # Elimination des satellites dont l'elevation est < 7 degree
            if j < 7:
                print(j)
            else:
                continue

        print("elevations ref:\n",elevations_ref, "elevations:\n", elevations,
              "azimuts ref:\n", azimuts_ref, "azimuts", azimuts, "hauteur d'antenne", hs_ref)

        """ Correction du retard tropospheriaue """
        P = 1013.15
        T = 285
        Rh =0.5
        e = Rh/100*np.exp(-37.2465 + 0.213166*T - 0.000256908 * T**2)
        for i in range(len(elevations)):
            dtropo = 0.002277/np.cos(elevations[i])*(P + (1255/T + 0.05) * e)
            PR_cor[i] -= dtropo

        """ Correction des effets ionospheriques """


        return PosSat_local,dtropo

    def bifreq(self, C1, C2):
        f1 = 1574.42
        f2 = 1227.60
        gamma = (f1/f2)**2
        Pionofree = (C2 - gamma * C1)/(1-gamma)


def main():
    """ lecture du fichier BRDC """
    PosSat = np.genfromtxt("possat.txt", skip_header=1, delimiter=",")
    Dobs = np.genfromtxt("dobs.txt", skip_header=1, delimiter=",")

    #    Orb = orbits.orbit()
    Orb = Tp3.orbit_TP()

    dir_orb = '../Guerledan'
    Orb.load_rinex_n(os.path.join(dir_orb, 'BRDC00IGS_R_20182850000_01D_MN.rnx'))

    """ lecture du fichier SP3 """
    sp3 = orbits.orbit()
    sp3.load_sp3(os.path.join(dir_orb, 'COM20225_15M.SP3'))

    """ Definition de l'instant pour lequel on cherche une position """
    t = gpst.gpstime(yyyy=2018, doy=285, dsec=5400)
    # print(t)

    rnx = rx.rinex_o()
    filename = os.path.join(dir_orb, 'edf1285b.18o')
    ret = rnx.load_rinex_o(filename)

    if ret < 0:
        print(ret)
        return

    Ep = rnx.get_epoch_by_mjd(t.mjd)

    #    spp1 = gnss_process_TP()
    #    spp1.const='G'
    #    spp1.constraint=0
    #    spp1.cut_off = 10 * d2r
    #    spp1.X0[0]=rnx.headers[0].X
    #    spp1.X0[1]=rnx.headers[0].Y
    #    spp1.X0[2]=rnx.headers[0].Z
    #    spp1.spp(Ep,Orb)

    spp2 = proc.gnss_process()
    spp2.const = 'G'
    prn = 14
    ordre = 9

    X_com, cdtr, V, sigma0_2, Qxx = orb_eph.estimation()
    posrecep = spp2.spp(Ep, Orb)
    print(posrecep)


#    print("Calcul sur les orbites brdc")
#    #Ep2 = spp2.spp(Ep,Orb)
#    #print("X = %.2fm Y = %.2fm Z = %.2fm " % (Ep2.X, Ep2.Y, Ep2.Z))
##
##    print("PRN           %-16s %-16s %-16s %-16s" % ('Xs', 'Ys', 'Zs', 'PR'))
##    for s in Ep2.satellites:
##        print("%1s%02d %16.3f %16.3f %16.3f %16.3f" % (s.const,s.PRN, s.Xs, s.Ys, s.Zs, s.PR))
##    print("V",Ep2.V)
##
#    print("Calcul sur les orbites sp3")
#    Ep3 = spp2.spp(Ep,sp3)
#    print("X = %.2fm Y = %.2fm Z = %.2fm " % (Ep3.X, Ep3.Y, Ep3.Z))
#
#    print("PRN           %-16s %-16s %-16s %-16s" % ('Xs', 'Ys', 'Zs', 'PR'))
#    for s in Ep3.satellites:
#        print("%1s%02d %16.3f %16.3f %16.3f %16.3f" % (s.const,s.PRN, s.Xs, s.Ys, s.Zs, s.PR))
#    print("V",Ep3.V)
#    print(Ep2.__dict__)

#    [E,N,U]= tool.tool_cartloc_GRS80(float(spp2.X0[0]),float(spp2.X0[1]),float(spp2.X0[2]),E2.X,E2.Y,E2.Z)
#    print("dE = %.2fm dN = %.2fm dU = %.2fm " % (E,N,U))

if __name__ == '__main__':


    tic = gpst.gpstime()
    main()
    toc = gpst.gpstime()
    print ('%.3f sec elapsed ' % (toc-tic))
