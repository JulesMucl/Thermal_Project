
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
from scipy.optimize import least_squares
from scipy.integrate import quad




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

        # Cold Fluid
        self.cold_fluid          = parameters['cold_fluid']
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
        self.eta_mec   = parameters['eta_mec']

        # Impose

        self.T_surchauffe_3 = parameters['T_surchauffe_3']
        self.T_surchauffe_4 = parameters['T_surchauffe_4']


        # ETAT REF
        self.p_ref, self.T_ref = parameters['p_ref'],   parameters['T_ref']
        self.h_ref = PropsSI('H','P',self.p_ref,'T',self.T_ref,self.fluid)
        self.s_ref = PropsSI('S','P',self.p_ref,'T',self.T_ref,self.fluid)


        # Guess
        self.p_1_guess = parameters['p_1_guess']
        self.p_2_guess = parameters['p_2_guess']
        self.p_5_guess = parameters['p_5_guess']

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


        #region ETATS

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


        # On pose les pressions en premier guess
        def cycle(p_guess) :
            p_1_guess, p_2_guess, p_5_guess = p_guess
            
            # Pressions avec T_surchauffe et sous refroidissement
            self.p1_guess_plot.append(p_1_guess)
            self.p2_guess_plot.append(p_2_guess)
            self.p5_guess_plot.append(p_5_guess)
            self.p_4 = p_1_guess
            self.p_3 = p_2_guess
            self.p_6 = p_5_guess

            self.T_3 = self.T_surchauffe_3 + PropsSI("T","P",p_2_guess,"Q",1,self.fluid) 
            self.T_4 = self.T_surchauffe_4 + PropsSI("T","P",p_1_guess,"Q",1,self.fluid) 
            self.T_6 = - self.T_cd_subcool + PropsSI("T","P",p_5_guess,"Q",0,self.fluid)
            print("P1 === ",p_1_guess)
            print("P2 === ",p_2_guess)
            print("P5 === ",p_5_guess)


            self.h_3 = PropsSI("H","P",self.p_3,"T",self.T_3,self.fluid)
            self.h_4 = PropsSI("H","P",self.p_4,"T",self.T_4,self.fluid)
            self.h_6 = PropsSI("H","P",self.p_6,"T",self.T_6,self.fluid)


            # Iteration sur T_1 via pitch
            def pitch_evap_I(T1) :
                h_1 = PropsSI("H","P",p_1_guess,"T",T1,self.fluid)
                T_1_prim = PropsSI("T","P",p_1_guess,"Q", 0,self.fluid)
                h_1_prim = PropsSI("H","P",p_1_guess,"Q", 0,self.fluid)

                m_1 = self.m_HF*(self.h_8 - self.h_9)/(self.h_4 - h_1)

                h_8_prim = self.h_8 - m_1/self.m_HF*(self.h_4 - h_1_prim)
                T_8_prim = PropsSI("T","P",self.p_7,"H",h_8_prim,self.hot_fluid)

                return T_8_prim - T_1_prim - self.T_pinch_ex_I
            
            self.T_1 = least_squares(pitch_evap_I, x0=10 + 273.15, bounds=(273.15, self.T_9)).x[0]
            self.h_1 = PropsSI("H","P",p_1_guess,"T",self.T_1,self.fluid)

            self.m_1 = self.m_HF*(self.h_8 - self.h_9)/(self.h_4 - self.h_1)
            print("m_1 === ",self.m_1)
            # Iteration sur T_2 via pitch

            def pitch_evap_II(T2) :
                h_2 = PropsSI("H","P",p_2_guess,"T",T2,self.fluid)
                T_2_prim = PropsSI("T","P",p_2_guess,"Q", 0,self.fluid)
                h_2_prim = PropsSI("H","P",p_2_guess,"Q", 0,self.fluid)

                m_2 = self.m_HF*(self.h_7 - self.h_8)/(self.h_3 - h_2)

                h_7_prim = self.h_7 - m_2/self.m_HF*(self.h_3 - h_2_prim)
                T_7_prim = PropsSI("T","P",self.p_7,"H",h_7_prim,self.hot_fluid)

                return T_7_prim - T_2_prim - self.T_pinch_ex_II

            self.T_2 = least_squares(pitch_evap_II, x0=10.1 + 273.15, bounds=(273.15, self.T_3 - self.T_surchauffe_3 - 1)).x[0] 
            self.h_2 = PropsSI("H","P",p_2_guess,"T",self.T_2,self.fluid)

            self.m_2 = self.m_HF*(self.h_7 - self.h_8)/(self.h_3 - self.h_2)
            self.m_tot = self.m_1 + self.m_2
            print("m_tot === ",self.m_tot)
            print("T6 === ",self.T_6)

            # Etat 5 via pinch
            def pitch_cond(T5) : 
                print("T5 === ",T5)
                T_6_sat = self.T_6 + self.T_cd_subcool
                h_5 = PropsSI("H","P",p_5_guess,"T",T5,self.fluid)
                T_5_prim = PropsSI("T","P",p_5_guess,"Q", 1,self.fluid)
                h_5_prim = PropsSI("H","P",p_5_guess,"Q", 1,self.fluid)

                m_CF = self.m_tot*(h_5 - self.h_6)/(self.h_10 - self.h_11)
                h_10_prim = self.h_11 + self.m_tot/m_CF*(h_5_prim - self.h_6)
                T_10_prim = PropsSI("T","P",self.p_CF,"H",h_10_prim,self.cold_fluid)
                return T_5_prim - T_10_prim - self.T_pinch_cd

            self.T_5 = least_squares(pitch_cond, x0=self.T_6 + 10, bounds=(self.T_6 + 0.1, 273.15 + 100)).x[0]
            self.h_5 = PropsSI("H","P",p_5_guess,"T",self.T_5,self.fluid)

            # Etat 3_prim
            self.s_3 = PropsSI("S","P",self.p_3,"T",self.T_3,self.fluid)
            self.s_prim = PropsSI("S","P",p_5_guess,"Q", 0,self.fluid)
            self.s_prim_prim = PropsSI("S","P",p_5_guess,"Q", 1,self.fluid)
            self.x_3_prim_s = (self.s_3 - self.s_prim) / (self.s_prim_prim - self.s_prim)
            self.h_3_prime_s = (1 - self.x_3_prim_s) * PropsSI("H","P",p_5_guess,"Q",0,self.fluid) + self.x_3_prim_s * PropsSI("H","P",p_5_guess,"Q",1,self.fluid)
            self.h_3_prime = self.h_3 - self.eta_is_T * (self.h_3 - self.h_3_prime_s)
            self.T_3_prime = PropsSI("T","P",p_5_guess,"H",self.h_3_prime,self.fluid)

            # Etat 4_prim
            self.s_4 = PropsSI("S","P",self.p_4,"T",self.T_4,self.fluid)
            self.x_4_prim_s = (self.s_4 - self.s_prim) / (self.s_prim_prim - self.s_prim)
            self.h_4_prime_s = (1 - self.x_4_prim_s) * PropsSI("H","P",p_5_guess,"Q",0,self.fluid) + self.x_4_prim_s * PropsSI("H","P",p_5_guess,"Q",1,self.fluid)

            self.h_4_prime = self.h_4 - self.eta_is_T * (self.h_4 - self.h_4_prime_s)
            self.T_4_prime = PropsSI("T","P",p_5_guess,"H",self.h_4_prime,self.fluid)

        
            ## FAIRE VIA ENTHALPIES
            def vI(P):
                return 1/PropsSI('D', 'H',(self.h_6 + self.h_1)/2,'P',P,self.fluid)
            
            vdp_variableI, _ = quad(vI, p_5_guess, p_1_guess)
            
            Eq1 = self.h_1 - self.h_6 - (vdp_variableI/self.eta_pump_1)


            def vII(P):
                return 1/PropsSI('D', 'H',(self.h_2 + self.h_1)/2,'P',P,self.fluid)
            
            vdp_variableII, _ = quad(vII, p_1_guess, p_2_guess)
            
            Eq2 = self.h_2 - self.h_1 - (vdp_variableII/self.eta_pump_2)
            

            # Etat 5
            First_eq = - self.h_5 + (self.h_4_prime*self.m_1 + self.m_2*self.h_3_prime)/self.m_tot

            return First_eq, Eq1, Eq2


        # Résolution
        self.p_1, self.p_2, self.p_5 = least_squares(cycle, x0=[self.p_1_guess, self.p_2_guess, self.p_5_guess], bounds = ([135323,5000,5000], [10**7,10**7,150865])).x
        print("p_1 === ",self.p_1)
        print("p_2 === ",self.p_2)
        print("p_5 === ",self.p_5)
        print(self.m_1)
        print(self.m_2)


