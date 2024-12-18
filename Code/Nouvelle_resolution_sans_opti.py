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
from scipy.optimize import minimize
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
        

        # General 
        self.fluid               = parameters['fluid']

        # Pump
        self.eta_pump_1          = parameters['eta_pump_1']
        self.eta_pump_2          = parameters['eta_pump_2']
        self.T_pinch_cd          = parameters['T_pinch_cd']
        self.T_surchauffe_3      = parameters['T_surchauffe_3']
        self.T_surchauffe_4      = parameters['T_surchauffe_4']
        self.T_cd_subcool        = parameters['T_cd_subcool']

        # Turbine
        self.eta_is_T            = parameters['eta_is_T']

        # ETAT REF
        self.p_ref, self.T_ref = parameters['p_ref'],   parameters['T_ref']
        self.h_ref = PropsSI('H','P',self.p_ref,'T',self.T_ref,self.fluid)
        self.s_ref = PropsSI('S','P',self.p_ref,'T',self.T_ref,self.fluid)

        # Contraintes
        self.m_HF      = parameters['m_dot_HF']
        self.eta_mec   = parameters['eta_mec']

        # Impose

        # Guess
        self.p_1_guess = parameters['p_1_guess']
        self.p_2_guess = parameters['p_2_guess']
        self.p_5_guess = parameters['p_5_guess']

        # ETAT REF
        self.p_ref, self.T_ref = parameters['p_ref'],   parameters['T_ref']
        self.h_ref = PropsSI('H','P',self.p_ref,'T',self.T_ref,self.fluid)
        self.s_ref = PropsSI('S','P',self.p_ref,'T',self.T_ref,self.fluid)

        # Graphes

        self.T1_guess_plot = []
        self.T2_guess_plot = []
        self.T5_guess_plot = []
        self.n_iterations = 0


    def exergie(self,h_i,s_i):

        return ( h_i - self.h_ref ) - self.T_ref * ( s_i - self.s_ref ) 



    def evaluate(self):


        #region ETATS

            #region ETAT 7
        # Etat 7 connu car données d'entrée du problème.
        self.p_7 = self.p_HF 
        self.h_7 = PropsSI("H","T",self.T_7,"P",self.p_7,self.hot_fluid)
        self.s_7 = PropsSI("S","T",self.T_7,"P",self.p_7,self.hot_fluid)
        self.e_7 = self.exergie(self.h_7,self.s_7)
        self.x_7 = PropsSI("Q","T",self.T_7,"P",self.p_7,self.hot_fluid)
            #endregion etat 7

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

            #region ETAT 8
        self.p_8 = self.p_HF
        self.h_8 = PropsSI("H","T",self.T_8,"P",self.p_7,self.hot_fluid)
        self.s_8 = PropsSI("S","T",self.T_8,"P",self.p_7,self.hot_fluid)
        self.e_8 = self.exergie(self.h_8,self.s_8)
        self.x_8 = PropsSI("Q","T",self.T_8,"P",self.p_7,self.hot_fluid)

        #endregion etat 8

        self.T_4 = self.T_8 - self.T_surchauffe_4
        self.T_3 = self.T_7 - self.T_surchauffe_3
        

        # On pose les pressions en premier guess
        def cycle(T) :
            T1, T2, T5 = T
            self.T1_guess_plot.append(T1)
            self.T2_guess_plot.append(T2)
            self.T5_guess_plot.append(T5)
            self.n_iterations += 1
            
            print("T1 === ",T1)
            print("T2 === ",T2)
            print("T5 === ",T5)

            # Iteration sur p_1 via pitch
            def pitch_evap_I(p1) :
                h_1 = PropsSI("H","P",p1,"T",T1,self.fluid)
                T_1_prim = PropsSI("T","P",p1,"Q", 0,self.fluid)
                h_1_prim = PropsSI("H","P",p1,"Q", 0,self.fluid)
                h_4 = PropsSI("H","T",self.T_4,"P",p1,self.fluid)

                m_1 = self.m_HF*(self.h_8 - self.h_9)/(h_4 - h_1)
                h_8_prim = self.h_8 - m_1/self.m_HF*(h_4 - h_1_prim)
                T_8_prim = PropsSI("T","P",self.p_7,"H",h_8_prim,self.hot_fluid)
                return T_8_prim - T_1_prim - self.T_pinch_ex_I

            solution = least_squares(pitch_evap_I, self.p_1_guess, bounds = (5000, 10**7))
            self.p_1 = solution.x[0]
            self.p_4 = self.p_1
            self.h_1 = PropsSI("H","P",self.p_1,"T",T1,self.fluid)
            self.h_4 = PropsSI("H","T",self.T_4,"P",self.p_4,self.fluid)
            self.m_1 = self.m_HF*(self.h_8 - self.h_9)/(self.h_4 - self.h_1)
            # Iteration sur T_2 via pitch

            def pitch_evap_II(p2) :
                h_2 = PropsSI("H","P",p2,"T",T2,self.fluid)
                T_2_prim = PropsSI("T","P",p2,"Q", 0,self.fluid)
                h_2_prim = PropsSI("H","P",p2,"Q", 0,self.fluid)
                h_3 = PropsSI("H","T",self.T_3,"P",p2,self.fluid)

                m_2 = self.m_HF*(self.h_7 - self.h_8)/(h_3 - h_2)

                h_7_prim = self.h_7 - m_2/self.m_HF*(h_3 - h_2_prim)
                T_7_prim = PropsSI("T","P",self.p_7,"H",h_7_prim,self.hot_fluid)

                return T_7_prim - T_2_prim - self.T_pinch_ex_II

            self.p_2 = fsolve(pitch_evap_II, self.p_2_guess)[0]
            self.p_3 = self.p_2
            self.h_2 = PropsSI("H","P",self.p_2,"T",T2,self.fluid)
            self.h_3 = PropsSI("H","T",self.T_3,"P",self.p_2,self.fluid)
            self.m_2 = self.m_HF*(self.h_7 - self.h_8)/(self.h_3 - self.h_2)
            self.m_tot = self.m_1 + self.m_2
            print("m_tot === ",self.m_tot)
            # Etat 5 via pinch
            def pitch_cond(T6) :
                T_6_sat = T6 + self.T_cd_subcool
                p6 = PropsSI("P","T",T_6_sat[0],"Q",0,self.fluid)
                T_5_prim = PropsSI("T","P",p6,"Q", 1,self.fluid)
                h_5_prim = PropsSI("H","P",p6,"Q", 1,self.fluid)
                h_6 = PropsSI("H","P",p6,"T", T6,self.fluid)
                h_5 = PropsSI("H","P",p6,"T", T5,self.fluid)
                m_CF = self.m_tot*(h_5 - h_6)/(self.h_10 - self.h_11)
                h_10_prim = self.h_11 + self.m_tot/m_CF*(h_5_prim - h_6)
                T_10_prim = PropsSI("T","P",self.p_CF,"H",h_10_prim,self.cold_fluid)
                return T_5_prim - T_10_prim - self.T_pinch_cd


            solution2 = least_squares(pitch_cond,  x0 = self.T_11+0.1, bounds=(self.T_11, T1))
            self.T_6 = solution2.x[0]
            self.p_6 = PropsSI("P","T",self.T_6 + self.T_cd_subcool,"Q",1,self.fluid)
            self.p_5 = self.p_6
            self.h_5 = PropsSI("H","P",self.p_5,"T",T5,self.fluid)
            self.h_6 = PropsSI("H","P",self.p_6,"T",self.T_6,self.fluid)
            self.m_dot_CF = self.m_tot*(self.h_5 - self.h_6)/(self.h_10 - self.h_11)
            # Etat 3_prim
            self.s_3 = PropsSI("S","P",self.p_2,"T",self.T_3,self.fluid)
            self.h_3_prime_s = PropsSI("H","P",self.p_6,"S",self.s_3,self.fluid)
            self.h_3_prime = self.h_3 - self.eta_is_T * (self.h_3 - self.h_3_prime_s)
            
            # Etat 4_prim
            self.s_4 = PropsSI("S","P",self.p_1,"T",self.T_4,self.fluid)
            self.h_4_prime_s = PropsSI("H","P",self.p_6,"S",self.s_4,self.fluid)
            self.h_4_prime = self.h_4 - self.eta_is_T * (self.h_4 - self.h_4_prime_s)
            

            ## FAIRE VIA ENTHALPIES
            def vI(P):
                return 1/PropsSI('D', 'H',(self.h_6 + self.h_1)/2,'P',P,self.fluid)
            
            vdp_variableI, _ = quad(vI, self.p_6, self.p_1)
            
            Eq1 = self.h_1 - self.h_6 - (vdp_variableI/self.eta_pump_1)


            def vII(P):
                return 1/PropsSI('D', 'H',(self.h_2 + self.h_1)/2,'P',P,self.fluid)
            
            vdp_variableII, _ = quad(vII, self.p_1, self.p_2)
            
            Eq2 = self.h_2 - self.h_1 - (vdp_variableII/self.eta_pump_2)
            
            
            First_eq = self.h_4_prime*self.m_1 + self.m_2*self.h_3_prime - self.m_tot*self.h_5
            return First_eq, Eq1, Eq2

        # Initial guess
        T_1_guess = 25.3 + 273.15
        T_2_guess = 25.6 + 273.15
        T_5_guess = 65 + 273.15

        # Résolution
        solution = least_squares(cycle, [T_1_guess, T_2_guess, T_5_guess], bounds = ([273.15 + 20,273.15 + 20,self.T_11],[self.T_8,self.T_7,273.15 + 120]), max_nfev=600)
        self.T1 = solution.x[0]
        self.T2 = solution.x[1]
        self.T5 = solution.x[2]


        # Calcul des états
        def pitch_evap_I_vrai(p1) :
            h_1 = PropsSI("H","P",p1,"T",self.T1,self.fluid)
            T_1_prim = PropsSI("T","P",p1,"Q", 0,self.fluid)
            h_1_prim = PropsSI("H","P",p1,"Q", 0,self.fluid)
            h_4 = PropsSI("H","T",self.T_4,"P",p1,self.fluid)

            m_1 = self.m_HF*(self.h_8 - self.h_9)/(h_4 - h_1)

            h_8_prim = self.h_8 - m_1/self.m_HF*(h_4 - h_1_prim)
            T_8_prim = PropsSI("T","P",self.p_7,"H",h_8_prim,self.hot_fluid)

            return T_8_prim - T_1_prim - self.T_pinch_ex_I
        
        self.p_1 = least_squares(pitch_evap_I_vrai, self.p_1_guess, bounds = (5000, 10**7)).x[0]
        self.h_1 = PropsSI("H","P",self.p_1,"T",self.T1,self.fluid) 
        self.s_1 = PropsSI("S","P",self.p_1,"T",self.T1,self.fluid)
        self.x_1 = PropsSI("Q","P",self.p_1,"T",self.T1,self.fluid)
        self.e_1 = self.exergie(self.h_1,self.s_1)

        self.p_4 = self.p_1
        self.h_4 = PropsSI("H","T",self.T_4,"P",self.p_4,self.fluid)
        self.s_4 = PropsSI("S","T",self.T_4,"P",self.p_4,self.fluid)
        self.x_4 = PropsSI("Q","T",self.T_4,"P",self.p_4,self.fluid)
        self.e_4 = self.exergie(self.h_4,self.s_4)

        self.m_1 = self.m_HF*(self.h_8 - self.h_9)/(self.h_4 - self.h_1)

        def pitch_evap_II_vrai(p2) :
            h_2 = PropsSI("H","P",p2,"T",self.T2,self.fluid)
            T_2_prim = PropsSI("T","P",p2,"Q", 0,self.fluid)
            h_2_prim = PropsSI("H","P",p2,"Q", 0,self.fluid)
            h_3 = PropsSI("H","T",self.T_3,"P",p2,self.fluid)

            m_2 = self.m_HF*(self.h_7 - self.h_8)/(h_3 - h_2)

            h_7_prim = self.h_7 - m_2/self.m_HF*(h_3 - h_2_prim)
            T_7_prim = PropsSI("T","P",self.p_7,"H",h_7_prim,self.hot_fluid)

            return T_7_prim - T_2_prim - self.T_pinch_ex_II

        self.p_2 = fsolve(pitch_evap_II_vrai, self.p_2_guess)[0]
        self.h_2 = PropsSI("H","P",self.p_2,"T",self.T2,self.fluid)
        self.s_2 = PropsSI("S","P",self.p_2,"T",self.T2,self.fluid)
        self.x_2 = PropsSI("Q","P",self.p_2,"T",self.T2,self.fluid)
        self.e_2 = self.exergie(self.h_2,self.s_2)

        self.p_3 = self.p_2
        self.h_3 = PropsSI("H","T",self.T_3,"P",self.p_3,self.fluid)
        self.s_3 = PropsSI("S","T",self.T_3,"P",self.p_3,self.fluid)
        self.x_3 = PropsSI("Q","T",self.T_3,"P",self.p_3,self.fluid)
        self.e_3 = self.exergie(self.h_3,self.s_3)

        self.m_2 = self.m_HF*(self.h_7 - self.h_8)/(self.h_3 - self.h_2)
        self.m_tot = self.m_1 + self.m_2

        def pitch_cond_vrai(T6) :
            T_6_sat = T6 + self.T_cd_subcool
            p6 = PropsSI("P","T",T_6_sat[0],"Q",0,self.fluid)
            T_5_prim = PropsSI("T","P",p6,"Q", 1,self.fluid)
            h_5_prim = PropsSI("H","P",p6,"Q", 1,self.fluid)
            h_6 = PropsSI("H","P",p6,"T", T6,self.fluid)
            h_5 = PropsSI("H","P",p6,"T", self.T5,self.fluid)
            m_CF = self.m_tot*(h_5 - h_6)/(self.h_10 - self.h_11)
            h_10_prim = self.h_11 + self.m_tot/m_CF*(h_5_prim - h_6)
            T_10_prim = PropsSI("T","P",self.p_CF,"H",h_10_prim,self.cold_fluid)
            return T_5_prim - T_10_prim - self.T_pinch_cd


        self.T_6 = least_squares(pitch_cond_vrai,  x0 = self.T_11+0.1, bounds=(self.T_11, self.T1)).x[0]
        self.p_6 = PropsSI("P","T",self.T_6 + self.T_cd_subcool,"Q",0,self.fluid)
        self.h_6 = PropsSI("H","P",self.p_6,"T",self.T_6,self.fluid)
        self.s_6 = PropsSI("S","P",self.p_6,"T",self.T_6,self.fluid)
        self.x_6 = PropsSI("Q","P",self.p_6,"T",self.T_6,self.fluid)
        self.e_6 = self.exergie(self.h_6,self.s_6)

        self.p_5 = self.p_6
        self.h_5 = PropsSI("H","P",self.p_5,"T",self.T5,self.fluid)
        self.s_5 = PropsSI("S","P",self.p_5,"T",self.T5,self.fluid)
        self.x_5 = PropsSI("Q","P",self.p_5,"T",self.T5,self.fluid)
        self.e_5 = self.exergie(self.h_5,self.s_5)
        self.m_dot_CF = self.m_tot*(self.h_5 - self.h_6)/(self.h_10 - self.h_11)

        # Etat 3_prim
        self.s_prim = PropsSI("S","P",self.p_5,"Q", 0,self.fluid)
        self.s_prim_prim = PropsSI("S","P",self.p_5,"Q", 1,self.fluid)
        self.x_3_prim_s = (self.s_3 - self.s_prim) / (self.s_prim_prim - self.s_prim)
        self.h_3_prime_s = (1 - self.x_3_prim_s) * PropsSI("H","P",self.p_5,"Q",0,self.fluid) + self.x_3_prim_s * PropsSI("H","P",self.p_5,"Q",1,self.fluid)
        
        self.h_3_prime = self.h_3 - self.eta_is_T * (self.h_3 - self.h_3_prime_s)
        self.T_3_prime = PropsSI("T","P",self.p_5,"H",self.h_3_prime,self.fluid)
        self.s_3_prime = PropsSI("S","P",self.p_5,"H",self.h_3_prime,self.fluid)
        self.x_3_prime = PropsSI("Q","P",self.p_5,"H",self.h_3_prime,self.fluid)
        self.e_3_prime = self.exergie(self.h_3_prime,self.s_3_prime)

        # Etat 4_prim
        self.x_4_prim_s = (self.s_4 - self.s_prim) / (self.s_prim_prim - self.s_prim)
        self.h_4_prime_s = (1 - self.x_4_prim_s) * PropsSI("H","P",self.p_5,"Q",0,self.fluid) + self.x_4_prim_s * PropsSI("H","P",self.p_5,"Q",1,self.fluid)

        self.h_4_prime = self.h_4 - self.eta_is_T * (self.h_4 - self.h_4_prime_s)
        self.T_4_prime = PropsSI("T","P",self.p_5,"H",self.h_4_prime,self.fluid)
        self.s_4_prime = PropsSI("S","P",self.p_5,"H",self.h_4_prime,self.fluid)
        self.x_4_prime = PropsSI("Q","P",self.p_5,"H",self.h_4_prime,self.fluid)
        self.e_4_prime = self.exergie(self.h_4_prime,self.s_4_prime)

        # Plot evolution des valeurs de T1, T2 et T5
        self.pas = np.linspace(1,self.n_iterations + 1,self.n_iterations)
        plt.figure(figsize=(10, 6))
        plt.plot(self.pas, self.T1_guess_plot, label="T1_guess_plot", linestyle='-', marker='o')
        plt.plot(self.pas, self.T2_guess_plot, label="T2_guess_plot", linestyle='--', marker='x')
        plt.plot(self.pas, self.T5_guess_plot, label="T5_guess_plot", linestyle='-.', marker='s')
        plt.title("Comparison of T1_guess_plot, T2_guess_plot, and T5_guess_plot")
        plt.xlabel("Index")
        plt.ylabel("Values")
        plt.legend()
        plt.show()

        # self.n_line = np.linspace(1,self.n+1, self.n)
        # plt.figure(figsize=(10, 6))
        # plt.plot(self.n_line, self.liste_Pe, label="Pe", linestyle='-', marker='o')
        # plt.title("Evolution of Pe")
        # plt.xlabel("Index")
        # plt.ylabel("Values")
        # plt.legend()
        # plt.show()


        print("m_1 === ",self.m_1)
        print("m_2 === ",self.m_2)
        print("m_tot === ",self.m_tot)

        #region Rendements
                
        self.Pm = self.m_1 * (self.h_4 - self.h_4_prime) + self.m_2 * (self.h_3 - self.h_3_prime) - self.m_tot*(self.h_1 - self.h_6) - self.m_2*(self.h_2 - self.h_1)
        self.Pe = self.eta_mec * self.Pm

        self.eta_cyclen = (self.m_1*(self.h_4 - self.h_4_prime) + self.m_2*(self.h_3 - self.h_3_prime) - self.m_tot*(self.h_1 - self.h_6) - self.m_2*(self.h_2 - self.h_1)) / (self.m_2*(self.h_3 - self.h_2) + self.m_1*(self.h_4 - self.h_1))
        self.eta_cyclex = (self.m_1*(self.h_4 - self.h_4_prime) + self.m_2*(self.h_3 - self.h_3_prime) - self.m_tot*(self.h_1 - self.h_6) - self.m_2*(self.h_2 - self.h_1)) / (self.m_2*(self.e_3 - self.e_2) + self.m_1*(self.e_4 - self.e_1))
        
        self.e_f = PropsSI("U", "T", self.T_ref, "P", self.p_ref, self.fluid)
        self.e_exh_I = self.e_7 - self.e_8
        self.e_exh_II = self.e_8 - self.e_9
        self.eta_transex = (self.m_2*(self.e_3 - self.e_2) + self.m_1*(self.e_4 - self.e_1))/(self.m_HF*(self.e_f - self.e_exh_II - self.e_exh_I))
        self.eta_conden = self.m_tot*(self.h_5 - self.h_6) / (self.m_2*(self.h_3 - self.h_2) + self.m_1*(self.h_4 - self.h_1))
        self.eta_condex = self.m_tot*(self.e_5 - self.e_6) / (self.m_2*(self.e_3 - self.e_2) + self.m_1*(self.e_4 - self.e_1))
        self.eta_rotex = (self.m_1*(self.h_4 - self.h_4_prime) + self.m_2*(self.h_3 - self.h_3_prime)- self.m_tot*(self.h_1 - self.h_6) - self.m_2*(self.h_2 - self.h_1))/(self.m_1*(self.e_4 - self.e_4_prime) + self.m_2*(self.e_3 - self.e_3_prime)- self.m_tot*(self.e_1 - self.e_6) - self.m_2*(self.e_2 - self.e_1))

        self.loss_mec = self.Pm - self.Pe
        self.loss_conden = self.m_tot*(self.h_6 - self.h_5)
        self.loss_exhen = self.m_HF*self.h_9
        self.loss_condex = self.m_tot*(self.e_6 - self.e_5)
        self.loss_rotex = (1 - self.eta_rotex) * self.Pm
        self.loss_transex_I = self.m_HF*(self.e_7 - self.e_8) - self.m_1*(self.e_4 - self.e_1)
        self.loss_transex_II = self.m_HF*(self.e_8 - self.e_9) - self.m_2*(self.e_3 - self.e_2)
        self.loss_transex_tot = self.loss_transex_I + self.loss_transex_II
        self.loss_exhex = self.m_HF*self.e_9
        #endregion rendements


        if self.display:
            print(self.fluid)
            print("P1 === ",self.p_1/1000,"kPa")
            print("P2 === ",self.p_2/1000,"kPa")
            print("P5 === ",self.p_5/1000,"kPa")
            print("Rendement cycle energie === ",self.eta_cyclen)
            print("Rendement cycle exergie === ",self.eta_cyclex)
            print("Pe === ",self.Pe/1000,"kW")

            #region Graphes T-S
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
            T_23 = np.linspace(self.T2, self.T_3, 100)
            s_23 = []
            for T in T_23:
                if (T != PropsSI("T", "P", self.p_2, "Q", 1, self.fluid)):
                    s_23.append(PropsSI('S', 'T', T, 'P', self.p_2, self.fluid) / 1000)
                else:
                    s_23.append(PropsSI('S', 'T', T, 'Q', 0, self.fluid) / 1000)
            plt.plot(s_23, T_23, color='black')

            T_14 = np.linspace(self.T1, self.T_4, 100)
            s_14 = []
            for T in T_14:
                if T != PropsSI("T", "P", self.p_1, "Q", 1, self.fluid):
                    s_14.append(PropsSI('S', 'T', T, 'P', self.p_1, self.fluid) / 1000)
                else :
                    s_14.append(PropsSI('S', 'T', T, 'Q', 0, self.fluid) / 1000)
            plt.plot(s_14, T_14, color='black')

            T_56 = np.linspace(self.T5, self.T_6, 100)
            s_56 = []
            for T in T_56:
                if T != PropsSI("T", "P", self.p_5, "Q", 0, self.fluid):
                    s_56.append(PropsSI('S', 'T', T, 'P', self.p_5, self.fluid) / 1000)
                else:
                    s_56.append(PropsSI('S', 'T', T, 'Q', 0, self.fluid) / 1000)
            plt.plot(s_56, T_56, color='black')

            plt.plot([self.s_1/1000, self.s_2/1000], [self.T1, self.T2], color='black')
            plt.plot([self.s_6/1000, self.s_1/1000], [self.T_6, self.T1], color='black')
            plt.plot([self.s_3/1000, self.s_4/1000], [self.T_3, self.T_4], color='black')
            plt.plot([self.s_4/1000, self.s_5/1000], [self.T_4, self.T5], color='black')

            s = [self.s_1/1000, self.s_2/1000, self.s_3/1000, self.s_4/1000, self.s_5/1000, self.s_6/1000]
            T = [self.T1, self.T2, self.T_3, self.T_4, self.T5, self.T_6]
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

            #endregion Graphes T-S


            #region Pies Charts
            plt.figure(1)
            self.losses = [np.abs(np.round(self.loss_mec,2)), np.abs(np.round(self.loss_conden,2)), np.abs(np.round(self.Pe,2)), np.abs(np.round(self.loss_exhen,2))]
            self.labels = ["Pertes mécaniques = "+ str(np.abs(np.round(self.loss_mec*10**(-6),2))) + "MW", "Pertes au condenseur = "+ str(np.abs(np.round(self.loss_conden*10**(-6),2))) + "MW", "Puissance électrique = " + str(np.abs(np.round(self.Pe*10**(-6),2))) + "MW", "Pertes à la sortie = " + str(np.abs(np.round(self.loss_exhen*10**(-6),2))) + "MW"]
            plt.pie(self.losses, labels=self.labels, autopct='%1.1f%%', startangle=90)
            plt.title("Répartition energie. Energie totale = " + str(self.m_HF*self.h_7) + " MW")
            self.pie_en = plt.figure(1)
            plt.show()

            plt.figure(2)
            self.losses = [np.abs(self.loss_exhex), np.abs(self.loss_mec), np.abs(self.loss_condex), np.abs(self.Pe), np.abs(self.loss_rotex), np.abs(self.loss_transex_tot)]
            self.labels = ["Pertes à la sortie = "+ str(np.abs(np.round(self.loss_exhex*10**(-6),2))) + "MW", "Pertes mécaniques = "+ str(np.abs(np.round(self.loss_mec*10**(-6),2))) + "MW", "Pertes exergétiques au condenseur = "+ str(np.abs(np.round(self.loss_condex*10**(-6),2))) + "MW", "Puissance électrique = " + str(np.abs(np.round(self.Pe*10**(-6),2))) + "MW", "Pertes aux turbines et pompes = " + str(np.abs(np.round(self.loss_condex*10**(-6),2))) + "MW", "Pertes aux échangeurs de chaleur = "+ str(np.abs(np.round(self.loss_transex_tot*10**(-6),2))) + "MW"]
            plt.pie(self.losses, labels=self.labels, autopct='%1.1f%%', startangle=90)
            plt.title("Répartition exergie. Exergie totale = " + str((self.m_HF*self.e_7)*10**(-6)) + " MW")
            self.pie_ex = plt.figure(2)
            plt.show()


            #endregion Pies Charts







            # Optimisation T_surchauffe et T_sous-refroidissement, T_8
            # Impact type de fluide (1 dans chaque catégorie), Impact température extérieure