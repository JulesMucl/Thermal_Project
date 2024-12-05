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

        # HOT Fluid
        self.hot_fluid           = parameters['hot_fluid']
        self.T_8                 = parameters['T_8']
        self.T_7                 = parameters['T_7']
        self.T_9                 = parameters['T_9'] 
        self.T_pinch_ex_I        = parameters['T_pinch_ex_I']
        self.T_pinch_ex_II       = parameters['T_pinch_ex_II']   
        self.p_HF                 = parameters['p_HF']
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
        self.x5_lim    = 0.88
        self.x6        = 0
        self.p_7       = self.p_HF

        # Impose

        self.T_surchauffe = parameters['T_surchauffe']
        self.m_tot = parameters['m_tot']

        # Guess
        self.p_1_guess = parameters['p_1_guess']
        self.p_2_guess = parameters['p_2_guess']
        self.p_5_guess = parameters['p_5_guess']

       


        # ETAT REF
        self.p_ref, self.T_ref = parameters['p_ref'],   parameters['T_ref']
        self.h_ref = PropsSI('H','P',self.p_ref,'T',self.T_ref,self.fluid)
        self.s_ref = PropsSI('S','P',self.p_ref,'T',self.T_ref,self.fluid)

       
    
        # Contraintes
        self.x5_lim = 0.88
        self.x6 = 0

       
  

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
        self.h_8 = PropsSI("H","T",self.T_8,"P",self.p_7,self.hot_fluid)
        self.s_8 = PropsSI("S","T",self.T_8,"P",self.p_7,self.hot_fluid)
        self.e_8 = self.exergie(self.h_8,self.s_8)
        self.x_8 = PropsSI("Q","T",self.T_8,"P",self.p_7,self.hot_fluid)

            #endregion etat 8

            #region ETAT 9

        self.h_9 = PropsSI("H","T",self.T_9,"P",self.p_7,self.hot_fluid)
        self.s_9 = PropsSI("S","T",self.T_9,"P",self.p_7,self.hot_fluid)
        self.e_9 = self.exergie(self.h_9,self.s_9)
        self.x_9 = PropsSI("Q","T",self.T_9,"P",self.p_7,self.hot_fluid)
            #endregion etat 9

            #region ETAT 10
        self.h_10 = PropsSI("H","T",self.T_10,"P",self.p_CF,self.cold_fluid)
        self.s_10 = PropsSI("S","T",self.T_10,"P",self.p_CF,self.cold_fluid)
        self.e_10 = self.exergie(self.h_10,self.s_10)
        self.x_10 = PropsSI("Q","T",self.T_10,"P",self.p_CF,self.cold_fluid)
            #endregion etat 10

            #region ETAT 11
        self.h_11 = PropsSI("H","T",self.T_11,"P",self.p_CF,self.cold_fluid)
        self.s_11 = PropsSI("S","T",self.T_11,"P",self.p_CF,self.cold_fluid)
        self.e_11 = self.exergie(self.h_11,self.s_11)
        self.x_11 = PropsSI("Q","T",self.T_11,"P",self.p_CF,self.cold_fluid)
            #endregion etat 11


        

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
            self.h_3_prime = self.h_3 - self.eta_is_T * (self.h_3 - self.h_3_prime_s)
            self.T_3_prime = PropsSI("T","P",self.p_5_guess,"H",self.h_3_prime,self.fluid)
            self.s_3_prime = PropsSI("S","P",self.p_5_guess,"H",self.h_3_prime,self.fluid)
            self.e_3_prime = self.exergie(self.h_3_prime,self.s_3_prime)
            self.x_3_prime = PropsSI("Q","P",self.p_5_guess,"H",self.h_3_prime,self.fluid)
            #endregion etat 3_prime

            #region ETAT 4_PRIME
            self.h_4_prime_s = PropsSI("H","P",self.p_5_guess,"S",self.s_4,self.fluid)
            self.h_4_prime = self.h_4 - self.eta_is_T * (self.h_4 - self.h_4_prime_s)
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
                cp_average = self.CP_av(self.T_6 , T1_guess , self.p_5_guess , self.p_1_guess , self.fluid)
                
                return T1_guess - self.T_6 - (self.p_1_guess - self.p_5_guess) / (self.eta_pump_1 * cp_average)

            self.T_1 = fsolve(T_out_pumpI, self.T_6*1.02)[0]
            print("T1 === ",self.T_1-273.15,"°C")
            self.h_1 = PropsSI("H","P",self.p_1_guess,"T",self.T_1,self.fluid)
            self.s_1 = PropsSI("S","P",self.p_1_guess,"T",self.T_1,self.fluid)
            self.e_1 = self.exergie(self.h_1,self.s_1)
            self.x_1 = PropsSI("Q","P",self.p_1_guess,"T",self.T_1,self.fluid)
            #endregion etat 1

            #region ETAT 2
            def T_out_pumpII(T2_guess) :
                cp_average2 = self.CP_av(self.T_1 , T2_guess , self.p_1_guess , self.p_2_guess , self.fluid)
                
                return T2_guess - self.T_1 - (self.p_2_guess - self.p_1_guess) / (self.eta_pump_2 * cp_average2)
            
            self.T_2 = fsolve(T_out_pumpII, self.T_1*1.05)[0]
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


            def Etat_i_evaporator_I(self,p_3_guess):


                h_hs_ex = self.h_8
                h_hs_su = self.h_7
                h_cs_su = self.h_2
                h_cs_ex = self.h_3


                T_cs_i = PropsSI('T', 'P', p_3_guess , 'Q', 0, self.fluid) 
                h_cs_i = PropsSI('H', 'P', p_3_guess , 'Q', 0, self.fluid) 

                h_hs_i = h_hs_su - (self.dot_m_cf/self.m_1) * (h_cs_ex - h_cs_su) 
                T_hs_i = PropsSI('T', 'H', h_hs_i, 'P', self.p_7, self.hot_fluid)


                return h_hs_i, h_cs_i, T_hs_i, T_cs_i
            
            def Iter_function_evaportor_I(self,p_3_guess):
                    
                h_hs_i, h_cs_i, T_hs_i, T_cs_i = self.Etat_i(p_3_guess)
            
                return T_hs_i - T_cs_i - self.T_pinch_ex_I
            
            def Etat_i_evaporator_II(self,p_1_guess):

                h_hs_ex = self.h_9
                h_hs_su = self.h_8
                h_cs_su = self.h_1
                h_cs_ex = self.h_4


                T_cs_i = PropsSI('T', 'P', p_1_guess , 'Q', 0, self.fluid) 
                h_cs_i = PropsSI('H', 'P', p_1_guess , 'Q', 0, self.fluid) 

                h_hs_i = h_hs_su - (self.dot_m_cf/self.m_1) * (h_cs_ex - h_cs_su) 
                T_hs_i = PropsSI('T', 'H', h_hs_i, 'P', self.p_7, self.hot_fluid)


                return h_hs_i, h_cs_i, T_hs_i, T_cs_i
            
            def Iter_function_evaportor_II(self,p_1_guess):
                    
                h_hs_i, h_cs_i, T_hs_i, T_cs_i = self.Etat_i_evaporator_II(p_1_guess)
            
                return T_hs_i - T_cs_i - self.T_pinch_ex_II
            
            #endregion checks pinch

            initial_guess = [self.p_1_guess, self.p_2_guess, self.p_5_guess]

            soluce = fsolve(cycle, initial_guess)
            p_1, p_2, p_5 = soluce[0], soluce[1], soluce[2]

            return p_1, p_2, p_5

        self.p_1, self.p_2, self.p_5 = cycle(self.p_1_guess, self.p_2_guess, self.p_5_guess)


            


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

            print("P1 === ",self.p_1/1000,"kPa")
            print("P2 === ",self.p_2/1000,"kPa")
            print("P5 === ",self.p_5/1000,"kPa")
   


