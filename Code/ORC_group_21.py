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


#=============Mes imports ============

from Tools.Pump import Pump1, Pump2



class ORC(object):

    def __init__(self,inputs,parameters,display):

        self.display             = display
        self.T_max               = parameters['T_max']
        self.p_3                 = parameters['p_3']
        self.r_pump_1            = parameters['r_pump_1'] # Pour moi on peut pas les imposer !!!!
        self.r_pump_2            = parameters['r_pump_2'] # Pour moi on peut pas les imposer !!!! Car redondant avec d'aures infos -> on va être bloqué ! 
        self.fluid               = parameters['fluid']
        self.k_cc_ex1            = parameters['k_cc_ex1']
        self.k_cc_ex2            = parameters['k_cc_ex2']
        self.T_cold_fluid_in     = parameters['T_cold_fluid_in']
        self.T_cold_fluid_out    = parameters['T_cold_fluid_out']
        self.cold_fluid          = parameters['cold_fluid']
        self.hot_fluid           = parameters['hot_fluid']
        self.dot_m_ex            = parameters['dot_m_ex']   # mass flow du fluid chaud  10 kg/s
        self.T_hot_fluid_in_II   = parameters['T_hot_fluid_in_II'] 
        self.p_hot_fluid         = parameters['p_hot_fluid']
        self.eta_pump_1          = parameters['eta_pump_1']
        self.eta_pump_2          = parameters['eta_pump_2']
        self.T_cd_subcool        = parameters['T_cd_subcool']
        self.T_surchauffe        = parameters['T_surchauffe']
        self.eta_is_T            = parameters['eta_is_T']
        self.m_tot               = parameters['m_tot']
        self.m_dot_CF            = parameters['m_dot_CF']
        self.T_10                = parameters['T_10']
        self.T_11                = parameters['T_11']
        self.p_CF                = parameters['p_CF']

        
        
       

        self.T_pinch_ex = 1


        # ETAT REF
        self.p_ref, self.T_ref = parameters['p_ref'],   parameters['T_ref']
        self.h_ref = PropsSI('H','P',self.p_ref,'T',self.T_ref,self.fluid)
        self.s_ref = PropsSI('S','P',self.p_ref,'T',self.T_ref,self.fluid)

       
    
        # Contraintes
        self.x5_lim = 0.88
        self.x6 = 0
        self.T_7 = self.T_hot_fluid_in_II
        self.p_7 = self.p_hot_fluid

       
    def T_pinch_cd(self,T5, T_cold_in, T_cold_out):
        """
        Retourne la température T_pinch du condensateur
        Ce fait via des itérations

        A faire 
        """
        return 1

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
                lambda x: PropsSI('CPMASS', 'T', x, 'P', p, fluid),
                t1, t2
            )[0]
            
            # Capacité calorifique moyenne
            cp_av = cp_integrale / np.abs(t2 - t1)
            return cp_av
        except Exception as e:
            raise ValueError(f"Erreur lors du calcul de CP_av : {e}")



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

            #endregion etat 8

        
        

        # On pose les pressions en premier guess
        def cycle(p_1_guess, p_2_guess, p_5_guess) : 

            self.T_3 = self.T_surchauffe + PropsSI("T","P",p_2_guess,"Q",1,self.fluid)
            self.T_4 = self.T_surchauffe + PropsSI("T","P",p_1_guess,"Q",1,self.fluid)
            self.T_6 = - self.T_cd_subcool + PropsSI("T","P",p_5_guess,"Q",0,self.fluid)

            #region ETAT 3
            self.p_3 = self.p_2_guess
            self.h_3 = PropsSI("H","P",self.p_2_guess,"T",self.T_3,self.fluid)
            self.s_3 = PropsSI("S","P",self.p_2_guess,"T",self.T_3,self.fluid)
            self.e_3 = self.exergie(self.h_3,self.s_3)
            self.x_3 = PropsSI("Q","P",self.p_2_guess,"T",self.T_3,self.fluid)
            #endregion etat 3

            #region ETAT 4
            self.p_4 = self.p_1_guess
            self.h_4 = PropsSI("H","P",self.p_1_guess,"T",self.T_4,self.fluid)
            self.s_4 = PropsSI("S","P",self.p_1_guess,"T",self.T_4,self.fluid)
            self.e_4 = self.exergie(self.h_4,self.s_4)
            self.x_4 = PropsSI("Q","P",self.p_1_guess,"T",self.T_4,self.fluid)
            #endregion etat 4

            #region ETAT 6
            self.p_6 = self.p_5_guess
            self.h_6 = PropsSI("H","P",self.p_5_guess,"T",self.T_6,self.fluid)
            self.s_6 = PropsSI("S","P",self.p_5_guess,"T",self.T_6,self.fluid)
            self.e_6 = self.exergie(self.h_6,self.s_6)
            self.x_6 = PropsSI("Q","P",self.p_5_guess,"T",self.T_6,self.fluid)
            #endregion etat 6

            #region ETAT 5
            self.h_5 = self.m_dot_CF/self.m_tot * (self.h_10 - self.h_11) + self.h_6
            self.s_5 = PropsSI("S","P",self.p_5_guess,"H",self.h_5,self.fluid)
            self.e_5 = self.exergie(self.h_5,self.s_5)
            self.x_5 = PropsSI("Q","P",self.p_5_guess,"H",self.h_5,self.fluid)
            self.T_5 = PropsSI("T","P",self.p_5_guess,"H",self.h_5,self.fluid)
            #endregion etat 5

            # Via turbine : 
            #region ETAT 3_PRIME
            self.h_3_prime_s = PropsSI("H","P",self.p_5_guess,"S",self.s_3,self.fluid)
            self.h_3_prime = self.h_3 - self.eta_is_T * (self.h_3 - self.h_3_prim_s)
            self.T_3_prime = PropsSI("T","P",self.p_5_guess,"H",self.h_3_prime,self.fluid)
            self.s_3_prime = PropsSI("S","P",self.p_5_guess,"H",self.h_3_prime,self.fluid)
            self.e_3_prime = self.exergie(self.h_3_prime,self.s_3_prime)
            self.x_3_prime = PropsSI("Q","P",self.p_5_guess,"H",self.h_3_prime,self.fluid)
            #endregion etat 3_prime

            #region ETAT 4_PRIME
            self.h_4_prime_s = PropsSI("H","P",self.p_5_guess,"S",self.s_4,self.fluid)
            self.h_4_prime = self.h_4 - self.eta_is_T * (self.h_4 - self.h_4_prim_s)
            self.T_4_prime = PropsSI("T","P",self.p_5_guess,"H",self.h_4_prime,self.fluid)
            self.s_4_prime = PropsSI("S","P",self.p_5_guess,"H",self.h_4_prime,self.fluid)
            self.e_4_prime = self.exergie(self.h_4_prime,self.s_4_prime)
            self.x_4_prime = PropsSI("Q","P",self.p_5_guess,"H",self.h_4_prime,self.fluid)
            #endregion etat 4_prime

            #region débits massiques
            def equations(m) : 
                return [self.h_4_prime*m[1] + m[0]*self.h_3_prime - self.m_tot*self.h_5, m[0] + m[1] - self.m_tot]
            self.m_1, self.m_2 = fsolve(equations,[1,1])
            #endregion débits massiques

            #region ETAT 1
            def T_out_pumpI(T1_guess) :
                cp_average = self.CP_av(T1_guess,self.T_2,self.p_5_guess,self.p_1_guess,self.fluid)
                
                return T1_guess - self.T_6 - (self.p_1_guess - self.p_5_guess) / (self.eta_pump_1 * cp_average)

            self.T_1 = fsolve(T_out_pumpI, self.T_6*1.02)[0]
            self.h_1 = PropsSI("H","P",self.p_1_guess,"T",self.T_1,self.fluid)
            self.s_1 = PropsSI("S","P",self.p_1_guess,"T",self.T_1,self.fluid)
            self.e_1 = self.exergie(self.h_1,self.s_1)
            self.x_1 = PropsSI("Q","P",self.p_1_guess,"T",self.T_1,self.fluid)
            #endregion etat 1

            #region ETAT 2
            def T_out_pumpII(T2_guess) :
                cp_average2 = self.CP_av(T2_guess,self.T_3,self.p_1_guess,self.p_2_guess,self.fluid)
                
                return T2_guess - self.T_1 - (self.p_2_guess - self.p_1_guess) / (self.eta_pump_2 * cp_average2)
            
            self.T_2 = fsolve(T_out_pumpII, self.T_1*1.02)[0]
            self.h_2 = PropsSI("H","P",self.p_2_guess,"T",self.T_2,self.fluid)
            self.s_2 = PropsSI("S","P",self.p_2_guess,"T",self.T_2,self.fluid)
            self.e_2 = self.exergie(self.h_2,self.s_2)
            self.x_2 = PropsSI("Q","P",self.p_2_guess,"T",self.T_2,self.fluid)
            #endregion etat 2

            #region CHECKS PINCH

            def Etat_i_condenseur(self,p_5_guess):

                h6 = self.h_6

                h6i = PropsSI('H','P',p_5_guess,'Q',1,self.fluid)
                T6i = PropsSI('T','P',p_5_guess,'Q',1,self.fluid)

                hcsi = (self.m_dot_tot/self.dot_m_cf) * (h6i - h6) + self.h_11
                Tcsi = PropsSI('T','H',hcsi,'P',self.p_10,self.cold_fluid)

                return T6i, Tcsi
    
            def Iter_function_condenseur(self,p_5_guess):

                T6i, Tcsi = self.Etat_i(p_5_guess)

                return T6i - Tcsi - self.T_pinch_cd

            # A CORRIGER !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            # def Etat_i_evaporator_I(self,p_3_guess):

            #     h_hs_ex = self.h_8
            #     h_cs_su = self.h_2

            #     T_hs_i = PropsSI('T', 'P', p_3_guess , 'Q', 0, self.fluid) # T_hs_sat_liq
            #     h_hs_i = PropsSI('H', 'P', p_3_guess , 'Q', 0, self.fluid) # h_hs_sat_liq

            #     h_cs_i = (1/self.ratio_I) * (h_hs_i- h_hs_ex) + h_cs_su
            #     T_cs_i = PropsSI('T', 'H', h_cs_i, 'P', self.p_7, self.hot_fluid)

            #     return h_hs_i, h_cs_i, T_hs_i, T_cs_i
            
            # def Iter_function_evaportor_I(self,p_3_guess):
                    
            #     h_hs_i, h_cs_i, T_hs_i, T_cs_i = self.Etat_i(p_3_guess)
            
            #     return T_cs_i - T_hs_i - self.T_pinch_ex_I
            
            
            #endregion checks pinch


        # Faire une matrice comme au HMW 3 : plus petite 

 




    

        #endregion etats

        #region massflowrates


        #endregion massflowrates


        #region Rendements
        
        #endregion rendements


        if self.display:
            print(self.fluid)
            print("P6 === ",self.p_6/1000,"kPa")
            print("T6 === ",self.T_6)
            print("T1 === ",self.T_1)

            print("P2 === ",self.p_2/1000,"kPa")
            print("T2 === ",self.T_2)
            print("e2 === ",self.e_2/1000,"kJ/kg")

            print("P3 === ",self.p_3/1000,"kPa")
            print("T3 === ",self.T_3)
            print("e3 === ",self.e_3/1000,"kJ/kg")
   


