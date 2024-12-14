#
#===IMPORT PACKAGES============================================================
#
import matplotlib.pyplot as plt
from CoolProp.CoolProp import PropsSI
from scipy import integrate as intgr
from scipy.optimize import fsolve
import numpy as np
from math import exp
import scipy.integrate as spi
from math import log
import scipy as Scipy


#=============Mes imports ============

from Tools.Pump import Pump1, Pump2



class ORC(object):

    def __init__(self,inputs,parameters,display):

        self.display             = display

        # HOT Fluid
        self.hot_fluid           = parameters['hot_fluid']
        self.T_8                 = parameters['T_8']
        self.T_7                 = parameters['T_7']
        self.T_9                 = parameters['T_9'] 
        self.T_pinch_ex_I        = parameters['T_pinch_ex_I']
        self.T_pinch_ex_II       = parameters['T_pinch_ex_II']   
        self.p_HF                = parameters['p_HF']
        self.dot_m_ex            = parameters['dot_m_ex']   # mass flow du fluid chaud  10 kg/s

        # Cold Fluid
        self.cold_fluid          = parameters['cold_fluid']
        self.m_dot_CF            = parameters['m_dot_CF']
        self.T_10                = parameters['T_10']
        self.T_11                = parameters['T_11']
        self.p_CF                = parameters['p_CF']
        self.T_pinch_cd         = parameters['T_pinch_cd']
        
        # Etat 3
        self.T_max               = parameters['T_max']

        # General 
        self.fluid               = parameters['fluid']
        self.k_cc_ex1            = parameters['k_cc_ex1']
        self.k_cc_ex2            = parameters['k_cc_ex2']

        # Pump
        self.eta_pump_1          = parameters['eta_pump_1']
        self.eta_pump_2          = parameters['eta_pump_2']
        self.T_cd_subcool        = parameters['T_cd_subcool']

        # Turbine
        self.eta_is_T            = parameters['eta_is_T']

        # ETAT REF
        self.p_ref, self.T_ref = parameters['p_ref'],   parameters['T_ref']
        self.h_ref = PropsSI('H','P',self.p_ref,'T',self.T_ref,self.fluid)
        self.s_ref = PropsSI('S','P',self.p_ref,'T',self.T_ref,self.fluid)

        # Contraintes
        self.p_7       = self.p_HF
        self.m_HF      = parameters['m_dot_HF']

        # Impose

        self.T_surchauffe = parameters['T_surchauffe']
        self.m_tot = parameters['m_tot']


        # ETAT REF
        self.p_ref, self.T_ref = parameters['p_ref'],   parameters['T_ref']
        self.h_ref = PropsSI('H','P',self.p_ref,'T',self.T_ref,self.fluid)
        self.s_ref = PropsSI('S','P',self.p_ref,'T',self.T_ref,self.fluid)

        # Graphes

        self.p1_guess_plot = []
        self.p2_guess_plot = []
        self.p5_guess_plot = []

    
    def afficher_matrice(self, matrice, lignes, colonnes):
        """
        Affiche une matrice joliment formatée avec des étiquettes pour les lignes et colonnes.

        Args:
            matrice (list of list): La matrice à afficher.
            lignes (list): Les étiquettes des lignes.
            colonnes (list): Les étiquettes des colonnes.
        """
        # Affiche les en-têtes des colonnes
        print("    | " + "  ".join(f"{colonne:>6}" for colonne in colonnes))
        print("-" * (8 + 7 * len(colonnes)))

        # Affiche chaque ligne avec son étiquette
        for etiquette, ligne in zip(lignes, matrice):
            print(f"{etiquette:>3} | " + "  ".join(f"{valeur:>6}" for valeur in ligne))
        
  

    def exergie(self,h_i,s_i):

        return ( h_i - self.h_ref ) - self.T_ref * ( s_i - self.s_ref ) 


    def CP_av(self, t1, t2, p1, p2, fluid):
        """
        Calcul de la capacité calorifique moyenne massique (cp) entre deux températures et deux pressions.

        :param t1: Température de départ (en K)
        :param t2: Température d'arrivée (en K)
        :param p1: Pression de départ (en Pa)
        :param p2: Pression d'arrivée (en Pa)
        :return: Capacité calorifique moyenne massique (en J/(kg.K))
        """
        if t1 == t2:
            raise ValueError("t1 et t2 ne peuvent pas être égaux, sinon division par zéro.")
        
        p = (p1 + p2) / 2  # Pression moyenne

        try:
            # Calcul de l'intégrale de cp entre t1 et t2
            cp_integrale = intgr.quad(
                lambda x: PropsSI('CPMASS', 'T', x, 'P', p, fluid),t1, t2)[0]
            
            # Capacité calorifique moyenne
            cp_av = cp_integrale / np.abs(t2 - t1)
            return cp_av
        except Exception as e:
            raise ValueError(f"Erreur lors du calcul de CP_av : {e}")

    def CP(self,T, P,fluid):
        cp = PropsSI('CPMASS','P', P ,'T',T,fluid) 
        return cp

    def evaluate(self):


        #region ETATS connnus

            #region ETAT 7
        # Etat 7 connu car données d'entrée du problème. 
        self.h_7 = PropsSI("H","T",self.T_7,"P",self.p_7,self.hot_fluid)
        self.s_7 = PropsSI("S","T",self.T_7,"P",self.p_7,self.hot_fluid)
        self.e_7 = self.exergie(self.h_7,self.s_7)
        self.x_7 = PropsSI("Q","T",self.T_7,"P",self.p_7,self.hot_fluid)
            #endregion etat 7


            #region ETAT 8
        self.p_8 = self.p_HF
        self.h_8 = PropsSI("H","T",self.T_8,"P",self.p_7,self.hot_fluid)
        self.s_8 = PropsSI("S","T",self.T_8,"P",self.p_7,self.hot_fluid)
        self.e_8 = self.exergie(self.h_8,self.s_8)
        self.x_8 = PropsSI("Q","T",self.T_8,"P",self.p_7,self.hot_fluid)

            #endregion etat 8

            #region ETAT 9
        self.p_9 = self.p_HF
        self.h_9 = PropsSI("H","T",self.T_9,"P",self.p_7,self.hot_fluid)
        self.s_9 = PropsSI("S","T",self.T_9,"P",self.p_7,self.hot_fluid)
        self.e_9 = self.exergie(self.h_9,self.s_9)
        self.x_9 = PropsSI("Q","T",self.T_9,"P",self.p_7,self.hot_fluid)
            #endregion etat 9

            #region ETAT 10
        self.p_10 = self.p_CF
        self.h_10 = PropsSI("H","T",self.T_10,"P",self.p_CF,self.cold_fluid)
        self.s_10 = PropsSI("S","T",self.T_10,"P",self.p_CF,self.cold_fluid)
        self.e_10 = self.exergie(self.h_10,self.s_10)
        self.x_10 = PropsSI("Q","T",self.T_10,"P",self.p_CF,self.cold_fluid)
            #endregion etat 10

            #region ETAT 11
        self.p_11 = self.p_CF
        self.h_11 = PropsSI("H","T",self.T_11,"P",self.p_CF,self.cold_fluid)
        self.s_11 = PropsSI("S","T",self.T_11,"P",self.p_CF,self.cold_fluid)
        self.e_11 = self.exergie(self.h_11,self.s_11)
        self.x_11 = PropsSI("Q","T",self.T_11,"P",self.p_CF,self.cold_fluid)
            #endregion etat 11
        
        #endregion etats connus


    
        # On pose les pressions en premier guess
        def cycle(T_I_guess, T_II_guess, T_C_guess) :
            self.etats = np.ones((8,6)) # p, T, s, h, x, e et etats 1,2,3,3',4,4',5,6
            if (T_I_guess < 0) or (T_II_guess < 0) or (T_C_guess < 0) :
                raise ValueError("Température négative")

            self.T_3 = self.T_7 - T_II_guess 
            self.T_4 = self.T_8 - T_I_guess 
            self.T_6 = self.T_11 + T_C_guess
            print("T_3 === ",self.T_3)
            print("T_4 === ",self.T_4)
            print("T_6 === ",self.T_6)
            
            # Pressions avec T_surchauffe et sous refroidissement
            self.p_4 = PropsSI("P","T",self.T_4 - self.T_surchauffe,"Q",1,self.fluid)
            self.p_3 = PropsSI("P","T",self.T_3 - self.T_surchauffe,"Q",1,self.fluid)
            self.p_6 = PropsSI("P","T",self.T_6 + self.T_cd_subcool,"Q",0,self.fluid)
            print("P_4 === ",self.p_4)
            print("P_3 === ",self.p_3)
            print("P_6 === ",self.p_6)
            self.p_5 = self.p_6
            self.p_1 = self.p_4
            self.p_2 = self.p_3

            self.h_3 = PropsSI("H","P",self.p_3,"T",self.T_3,self.fluid)
            self.h_4 = PropsSI("H","P",self.p_4,"T",self.T_4,self.fluid)
            self.h_6 = PropsSI("H","P",self.p_6,"T",self.T_6,self.fluid)

            # On trouve T2 via pinch.


            #region ETAT 5
            self.h_5 = self.m_dot_CF/self.m_tot * (self.h_10 - self.h_11) + self.h_6
            self.s_5 = PropsSI("S","P",p_5_guess,"H",self.h_5,self.fluid)
            self.e_5 = self.exergie(self.h_5,self.s_5)
            self.x_5 = PropsSI("Q","P",p_5_guess,"H",self.h_5,self.fluid)
            self.T_5 = PropsSI("T","P",p_5_guess,"H",self.h_5,self.fluid)
            self.etats[6] = [p_5_guess, self.T_5, self.s_5, self.h_5, self.x_5, self.e_5]
            #endregion etat 5

            # Via turbine : 
            #region ETAT 3_PRIME
            self.s_prim = PropsSI("S","P",p_5_guess,"Q", 0,self.fluid)
            self.s_prim_prim = PropsSI("S","P",p_5_guess,"Q", 1,self.fluid)
            self.x_3_prim_s = (self.s_3 - self.s_prim) / (self.s_prim_prim - self.s_prim)
            self.h_3_prime_s = (1 - self.x_3_prim_s) * PropsSI("H","P",p_5_guess,"Q",0,self.fluid) + self.x_3_prim_s * PropsSI("H","P",p_5_guess,"Q",1,self.fluid)
            # print("h_3_prime_s === ",self.h_3_prime_s)
            # self.h_3_prime_s = PropsSI("H","P",p_5_guess,"S",self.s_3,self.fluid)
            # print("h_3_prime_s_bis === ",self.h_3_prime_s)
            self.h_3_prime = self.h_3 - self.eta_is_T * (self.h_3 - self.h_3_prime_s)
            self.T_3_prime = PropsSI("T","P",p_5_guess,"H",self.h_3_prime,self.fluid)
            self.s_3_prime = PropsSI("S","P",p_5_guess,"H",self.h_3_prime,self.fluid)
            self.e_3_prime = self.exergie(self.h_3_prime,self.s_3_prime)
            self.x_3_prime = PropsSI("Q","P",p_5_guess,"H",self.h_3_prime,self.fluid)
            self.etats[3] = [p_5_guess, self.T_3_prime, self.s_3_prime, self.h_3_prime, self.x_3_prime, self.e_3_prime]
            #endregion etat 3_prime

            #region ETAT 4_PRIME
            #self.h_4_prime_s = PropsSI("H","P",p_5_guess,"S",self.s_4,self.fluid)
            self.x_4_prim_s = (self.s_4 - self.s_prim) / (self.s_prim_prim - self.s_prim)
            self.h_4_prime_s = (1 - self.x_4_prim_s) * PropsSI("H","P",p_5_guess,"Q",0,self.fluid) + self.x_4_prim_s * PropsSI("H","P",p_5_guess,"Q",1,self.fluid)

            self.h_4_prime = self.h_4 - self.eta_is_T * (self.h_4 - self.h_4_prime_s)
            self.T_4_prime = PropsSI("T","P",p_5_guess,"H",self.h_4_prime,self.fluid)
            self.s_4_prime = PropsSI("S","P",p_5_guess,"H",self.h_4_prime,self.fluid)
            self.e_4_prime = self.exergie(self.h_4_prime,self.s_4_prime)
            self.x_4_prime = PropsSI("Q","P",p_5_guess,"H",self.h_4_prime,self.fluid)
            self.etats[5] = [p_5_guess, self.T_4_prime, self.s_4_prime, self.h_4_prime, self.x_4_prime, self.e_4_prime]
            #endregion etat 4_prime

            #region débits massiques
            def equations(m) : 
                return [self.h_4_prime*m[0] + m[1]*self.h_3_prime - self.m_tot*self.h_5, m[0] + m[1] - self.m_tot]
            self.m_1, self.m_2 = fsolve(equations,[100,100])
            #endregion débits massiques

            #region ETAT 1
            def T_out_pumpI(T1_guess, arg) :
                T_6, p_5_guess, p_1_guess, fluid = arg
                p_av = (p_1_guess + p_5_guess) / 2
                rho_av = Scipy.integrate.quad(lambda x: PropsSI('D','P',p_av,'T',x,fluid),T_6,T1_guess)[0] / (T1_guess - T_6)

                cp_average = Scipy.integrate.quad(lambda x : self.CP(x,p_av,fluid),T_6,T1_guess)[0] / (T1_guess - T_6)
                return T1_guess - T_6 - (p_1_guess - p_5_guess) / (self.eta_pump_1 * cp_average*rho_av)

            self.args_1 = [self.T_6, p_5_guess, p_1_guess, self.fluid]
            self.T_1 = fsolve(T_out_pumpI, self.T_6*1.02, args=self.args_1)[0]
            self.h_1 = PropsSI("H","P",p_1_guess,"T",self.T_1,self.fluid)
            self.s_1 = PropsSI("S","P",p_1_guess,"T",self.T_1,self.fluid)
            self.e_1 = self.exergie(self.h_1,self.s_1)
            self.x_1 = PropsSI("Q","P",p_1_guess,"T",self.T_1,self.fluid)
            self.etats[0] = [p_1_guess, self.T_1, self.s_1, self.h_1, self.x_1, self.e_1]
            #endregion etat 1

            #region ETAT 2
            def T_out_pumpII(T2_guess, arg) :
                T_1, p_1_guess, p_2_guess, fluid = arg
                p_av = (p_1_guess + p_2_guess) / 2
                rho_av = Scipy.integrate.quad(lambda x: PropsSI('D','P',p_av,'T',x,fluid),T_1,T2_guess)[0] / (T2_guess - T_1)

                cp_average = Scipy.integrate.quad(lambda x : self.CP(x,p_av,fluid),T_1,T2_guess)[0] / (T2_guess - T_1)
                return T2_guess - T_1 - (p_2_guess - p_1_guess) / (self.eta_pump_2 * cp_average * rho_av)
            
            self.args_2 = [self.T_1, p_1_guess, p_2_guess, self.fluid]
            self.T_2 = fsolve(T_out_pumpII, self.T_1*1.05, args=self.args_2)[0]
            self.h_2 = PropsSI("H","P",p_2_guess,"T",self.T_2,self.fluid)
            self.s_2 = PropsSI("S","P",p_2_guess,"T",self.T_2,self.fluid)
            self.e_2 = self.exergie(self.h_2,self.s_2)
            self.x_2 = PropsSI("Q","P",p_2_guess,"T",self.T_2,self.fluid)
            self.etats[1] = [p_2_guess, self.T_2, self.s_2, self.h_2, self.x_2, self.e_2]
            #endregion etat 2

            return self.h_1, self.h_2, self.h_3, self.h_3_prime, self.h_4, self.h_4_prime, self.h_5, self.h_6, self.m_1, self.m_2, self.etats

        #region PITCH CONDENSEUR
        def Etat_i_condenseur(p5_guess, h_6):

            h6i = PropsSI('H','P',p5_guess,'Q',1,self.fluid)
            T6i = PropsSI('T','P',p5_guess,'Q',1,self.fluid)

            hcsi = (self.m_tot/self.m_dot_CF) * (h6i - h_6) + self.h_11
            Tcsi = PropsSI('T','H',hcsi,'P',self.p_10,self.cold_fluid)

            return T6i, Tcsi

        def Iter_function_condenseur(p5_guess, h_6) :

            T6i, Tcsi = Etat_i_condenseur(p5_guess, h_6)

            return T6i - Tcsi - self.T_pinch_cd
        #endregion PITCH CONDENSEUR
        
        #region PITCH ECHANGEUR I

        def Etat_i_evaporator_I(p3_guess, h_2, h_3, m_2):

            h_hs_ex = self.h_8
            h_hs_su = self.h_7
            h_cs_su = h_2
            h_cs_ex = h_3

            T_cs_i = PropsSI('T', 'P', p3_guess , 'Q', 0, self.fluid) 
            h_cs_i = PropsSI('H', 'P', p3_guess , 'Q', 0, self.fluid) 

            h_hs_i = h_hs_su - (self.m_HF/m_2) * (h_cs_ex - h_cs_su) 
            T_hs_i = PropsSI('T', 'H', h_hs_i, 'P', self.p_7, self.hot_fluid)

            return h_hs_i, h_cs_i, T_hs_i, T_cs_i
        
        def Iter_function_evaporator_I(p3_guess, h_2, h_3, m_2):
                
            h_hs_i, h_cs_i, T_hs_i, T_cs_i = Etat_i_evaporator_I(p3_guess, h_2, h_3, m_2)
        
            return T_hs_i - T_cs_i - self.T_pinch_ex_I

        #endregion PITCH ECHANGEUR I
        
        #region PITCH ECHANGEUR II
        def Etat_i_evaporator_II(p1_guess, h_1, h_4, m_1):

            h_hs_ex = self.h_9
            h_hs_su = self.h_8
            h_cs_su = h_1
            h_cs_ex = h_4

            T_cs_i = PropsSI('T', 'P', p1_guess , 'Q', 0, self.fluid) 
            h_cs_i = PropsSI('H', 'P', p1_guess , 'Q', 0, self.fluid) 

            h_hs_i = h_hs_su - (self.m_HF/m_1) * (h_cs_ex - h_cs_su) 
            T_hs_i = PropsSI('T', 'H', h_hs_i, 'P', self.p_7, self.hot_fluid)

            return h_hs_i, h_cs_i, T_hs_i, T_cs_i
        
        def Iter_function_evaporator_II(p1_guess, h_1, h_4, m_1):
                
            h_hs_i, h_cs_i, T_hs_i, T_cs_i = Etat_i_evaporator_II(p1_guess, h_1, h_4, m_1)
        
            return T_hs_i - T_cs_i - self.T_pinch_ex_II

        #endregion PITCH ECHANGEUR II


        #region ITERATION

        def check(p_guess) : 


            p1_guess = p_guess[0]
            p2_guess = p_guess[1]
            p5_guess = p_guess[2]

            h_1, h_2, h_3, h_3_prime, h_4, h_4_prime, h_5, h_6, m_1, m_2, matrice = cycle(p1_guess, p2_guess, p5_guess)

            delta_T_condenseur = Iter_function_condenseur(p5_guess, h_6)
            delta_T_evaporator_I = Iter_function_evaporator_I(p2_guess, h_2, h_3, m_2)
            delta_T_evaporator_II = Iter_function_evaporator_II(p1_guess, h_1, h_4, m_1)

            return delta_T_condenseur, delta_T_evaporator_I, delta_T_evaporator_II


        self.T_I_guess = 20
        self.T_II_guess = 10
        self.T_C_guess = 5
        self.T_I , self.T_II , self.T_C = fsolve(check, [self.T_I_guess, self.T_II_guess, self.T_C_guess])
        #endregion ITERATION


        

 


        #region ETATS

        self.h_1, self.h_2, self.h_3, self.h_3_prime, self.h_4, self.h_4_prime, self.h_5, self.h_6, self.m_1, self.m_2, self.matrice = cycle(self.T_I, self.T_II, self.T_C)

        #endregion etats

        print("m_1 === ",self.m_1)
        print("m_2 === ",self.m_2)
        print("m_tot === ",self.m_tot)

        #region Rendements
        self.e_4_prime = self.matrice[5][5]
        self.e_3_prime = self.matrice[3][5]
        self.e_4 = self.matrice[4][5]
        self.e_3 = self.matrice[2][5]

        Pm = self.m_1 * (self.h_4 - self.h_4_prime) + self.m_2 * (self.h_3 - self.h_3_prime) - ( (self.h_1 - self.h_6) * self.m_tot + (self.h_2 - self.h_1) * self.m_2 )
        Q = self.m_2 * (self.h_7 - self.h_8) + self.m_1 * (self.h_8 - self.h_9)  # h_7 et h_8 sont les enthalpies du fluid chaud. ??
        self.eta_cyclen = Pm / Q
        print("PM === ",Pm)
        print("Q === ",Q)
        
        
        self.eta_cyclen = (self.m_1*(self.h_4 - self.h_4_prime) + self.m_2*(self.h_3 - self.h_3_prime) ) / (self.m_2*(self.h_3 - self.h_2) + self.m_1*(self.h_4 - self.h_1))
        self.eta_cyclex = (self.m_1*(self.e_4 - self.e_4_prime) + self.m_2*(self.e_3 - self.e_3_prime) ) / (self.m_2*(self.h_3 - self.h_2) + self.m_1*(self.h_4 - self.h_1))
        #endregion rendements


        if self.display:
            print(self.fluid)
            print("P1 === ",self.p_1/1000,"kPa")
            print("P2 === ",self.p_2/1000,"kPa")
            print("P5 === ",self.p_5/1000,"kPa")
            print("Rendement cycle energie === ",self.eta_cyclen)
            print("Rendement cycle exergie === ",self.eta_cyclex)
            lignes = ["1", "2", "3", "3'", "4", "4'", "5", "6"]
            colonnes = ["p", "T", "s", "h", "x", "e"]
            print("Matrice :")
            self.afficher_matrice(self.matrice, lignes, colonnes)


            # Courbe de saturation du fluide
            s_liq = []
            T_liq = []
            s_vap = []
            T_vap = []
            for T in np.linspace(PropsSI('Ttriple', self.fluid), PropsSI('Tcrit', self.fluid), 500):
                s_liq.append(PropsSI('S', 'T', T, 'Q', 0, self.fluid) / 1000)  # Entropie du liquide saturé
                T_liq.append(T)
                s_vap.append(PropsSI('S', 'T', T, 'Q', 1, self.fluid) / 1000)  # Entropie de la vapeur saturée
                T_vap.append(T)
            # Tracé du diagramme T-s
            plt.figure(figsize=(10, 6))
            # Tracé de la courbe de saturation
            plt.plot(s_liq, T_liq, linestyle='--', color='blue', label="Saturation curve")
            plt.plot(s_vap, T_vap, linestyle='--', color='blue')

            # Ajout courbes cycle
            T_23 = np.linspace(self.T_2, self.T_3, 100)
            s_23 = []
            for T in T_23:
                if (T != PropsSI("T", "P", self.p_2, "Q", 1, self.fluid)):
                    s_23.append(PropsSI('S', 'T', T, 'P', self.p_2, self.fluid) / 1000)
                else:
                    s_23.append(PropsSI('S', 'T', T, 'Q', 0, self.fluid) / 1000)
            plt.plot(s_23, T_23, color='black')

            T_14 = np.linspace(self.T_1, self.T_4, 100)
            s_14 = []
            for T in T_14:
                if T != PropsSI("T", "P", self.p_1, "Q", 1, self.fluid):
                    s_14.append(PropsSI('S', 'T', T, 'P', self.p_1, self.fluid) / 1000)
                else :
                    s_14.append(PropsSI('S', 'T', T, 'Q', 0, self.fluid) / 1000)
            plt.plot(s_14, T_14, color='black')

            T_56 = np.linspace(self.T_5, self.T_6, 100)
            s_56 = []
            for T in T_56:
                if T != PropsSI("T", "P", self.p_5, "Q", 0, self.fluid):
                    s_56.append(PropsSI('S', 'T', T, 'P', self.p_5, self.fluid) / 1000)
                else:
                    s_56.append(PropsSI('S', 'T', T, 'Q', 0, self.fluid) / 1000)
            plt.plot(s_56, T_56, color='black')

            plt.plot([self.s_1/1000, self.s_2/1000], [self.T_1, self.T_2], color='black')
            plt.plot([self.s_6/1000, self.s_1/1000], [self.T_6, self.T_1], color='black')
            plt.plot([self.s_3/1000, self.s_4/1000], [self.T_3, self.T_4], color='black')
            plt.plot([self.s_4/1000, self.s_5/1000], [self.T_4, self.T_5], color='black')

            s = [self.s_1/1000, self.s_2/1000, self.s_3/1000, self.s_4/1000, self.s_5/1000, self.s_6/1000]
            T = [self.T_1, self.T_2, self.T_3, self.T_4, self.T_5, self.T_6]
            labels = ["1", "2", "3", "4", "5", "6"]
            for i in range(len(labels)):
                plt.scatter(s[i], T[i], color='gray')
                plt.text(s[i], T[i], labels[i], color='gray', fontsize=12, ha='right')

            # Ajout courbes exterieure
            T_78 = np.linspace(self.T_7, self.T_8, 100)
            s_78 = []
            for T in T_78:
                if T != PropsSI("T", "P", self.p_7, "Q", 1, self.hot_fluid):
                    s_78.append(PropsSI('S', 'T', T, 'P', self.p_7, self.hot_fluid) / 1000)
                else:
                    s_78.append(PropsSI('S', 'T', T, 'Q', 0, self.hot_fluid) / 1000)
            plt.plot(s_78, T_78, color='red', label="Hot Fluid Evaporator II", linestyle='--')

            T_89 = np.linspace(self.T_8, self.T_9, 100)
            s_89 = []
            for T in T_89:
                if T != PropsSI("T", "P", self.p_7, "Q", 1, self.hot_fluid):
                    s_89.append(PropsSI('S', 'T', T, 'P', self.p_7, self.hot_fluid) / 1000)
                else:
                    s_89.append(PropsSI('S', 'T', T, 'Q', 0, self.hot_fluid) / 1000)
            plt.plot(s_89, T_89, color='yellow', linestyle='--', label="Hot Fluid Evaporator I")

            T_1110 = np.linspace(self.T_11, self.T_10, 100)
            s_1110 = []
            for T in T_1110:
                if T != PropsSI("T", "P", self.p_CF, "Q", 1, self.cold_fluid):
                    s_1110.append(PropsSI('S', 'T', T, 'P', self.p_CF, self.cold_fluid) / 1000)
                else:
                    s_1110.append(PropsSI('S', 'T', T, 'Q', 0, self.cold_fluid) / 1000)
            plt.plot(s_1110, T_1110, color='g', linestyle='--', label="Cold Fluid Condenser")


            # Configuration du graphique
            plt.title("Diagramme T-s du Cycle ORC avec Courbe de Saturation")
            plt.xlabel("Entropie (s) [kJ/(kg.K)]")
            plt.ylabel("Température (T) [K]")
            plt.grid(True)
            plt.legend()
            plt.show()




   


