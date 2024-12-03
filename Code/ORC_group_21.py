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
from Tools.Evaporator import Evaporator_I
from Tools.Condenser import Condenser




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
        self.p_9                 = parameters['p_9']
        self.dot_m_ex            = parameters['dot_m_ex']   # mass flow du fluid chaud  10 kg/s

        # Cold Fluid
        self.cold_fluid          = parameters['cold_fluid']
        self.m_dot_cf            = parameters['m_dot_cf']
        self.T_10                = parameters['T_10']
        self.T_11                = parameters['T_11']
        self.p_10                = parameters['p_10']
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

        # ETAT REF
        self.p_ref, self.T_ref = parameters['p_ref'],   parameters['T_ref']
        self.h_ref = PropsSI('H','P',self.p_ref,'T',self.T_ref,self.fluid)
        self.s_ref = PropsSI('S','P',self.p_ref,'T',self.T_ref,self.fluid)

        # Contraintes
        self.x5_lim    = 0.88
        self.x6        = 0
        self.p_7       = self.p_9

        # Impose
        self.m_dot_tot = parameters['m_dot_tot']
        self.X_I       = 0.5
        self.X_II      = 0.5

        self.dot_m_I = self.m_dot_tot * self.X_I
        self.dot_m_II = self.m_dot_tot * self.X_II

        # Guess
        self.p_3_guess = parameters['p_3_guess']
        self.p_6_guess = parameters['p_6_guess']



       
 
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

        self.h_10 = PropsSI("H","T",self.T_10,"P",self.p_10,self.cold_fluid)
        self.h_11 = PropsSI("H","T",self.T_11,"P",self.p_10,self.cold_fluid)


        #region ETATS

            #region ETAT 7 connu car données d'entrée du problème
        
        self.h_7 = PropsSI("H","T",self.T_7,"P",self.p_7,self.hot_fluid)
        self.s_7 = PropsSI("S","T",self.T_7,"P",self.p_7,self.hot_fluid)
        self.e_7 = self.exergie(self.h_7,self.s_7)
        self.x_7 = PropsSI("Q","T",self.T_7,"P",self.p_7,self.hot_fluid)
            #endregion etat 7


            # region ETAT 6
        self.condenser = Condenser(self)
        self.p_6 = self.condenser.P6()
        self.T_6 = PropsSI("T","P",self.p_6,"Q",0,self.fluid) - self.T_cd_subcool
        self.h_6 = PropsSI("H","P",self.p_6,"T", self.T_6 ,self.fluid)
        self.s_6 = PropsSI("S","P",self.p_6,"T", self.T_6 ,self.fluid)  
        self.e_6 = self.exergie(self.h_6,self.s_6)
        print('P6 === ',self.p_6/11e5,'bar')
        print('T6 === ',self.T_6-273.15,'°C')
        print('h6 === ',self.h_6/1000,'kJ/kg')
        print('s6 === ',self.s_6/1000,'kJ/kg.K')
        print('e6 === ',self.e_6/1000,'kJ/kg')

            #endregion etat6



        #     #region ETAT 2
        # self.p_2 = self.p_1 * self.r_pump_2
        # self.pump_2 = Pump2(self)
        # self.T_2 = self.pump_2.evaluate_T_out()
        # self.h_2 = PropsSI("H","T",self.T_2,"P",self.p_2,self.fluid)
        # self.s_2 = PropsSI("S","T",self.T_2,"P",self.p_2,self.fluid)
        # self.e_2 = self.exergie(self.h_2,self.s_2)

            #endregion etat2

             # region ETAT 3

        self.evaporator = Evaporator_I(self)
        self.p_3 = self.evaporator.Pressure()

        self.T_3 = PropsSI("T","P",self.p_3,"Q",1,self.fluid)
        print('T3 === ',self.T_3-273.15,'°C')
        print("p3 === ",self.p_3/1000,"kPa")
        self.h_3 = PropsSI("H","Q",0,"P",self.p_3,self.fluid)
        self.s_3 = PropsSI("S","T",self.T_3,"P",self.p_3,self.fluid)
        self.e_3 = self.exergie(self.h_3,self.s_3)

             # endregion etat 3

             # region ETAT 5
        self.T_5 = self.T_6 + self.T_cd_subcool

        # Via le rendement isentropique de la turbine : 
        self.s_5s = self.s_3
        self.x_5s = PropsSI("Q","T",self.T_5,"S",self.s_5s,self.fluid)
        self.h_5s = PropsSI("H","T",self.T_5,"Q",self.x_5s,self.fluid)
        self.h_5 = self.h_3 - (self.h_3 - self.h_5s) * self.eta_is_HP
        self.s_5 = PropsSI("S","T",self.T_5,"H",self.h_5,self.fluid)
        self.e_5 = self.exergie(self.h_5,self.s_5)
        self.x_5 = PropsSI("Q","T",self.T_5,"H",self.h_5,self.fluid)
        if(self.x_5 < self.x5_lim):
            print("Error : x5 < x5_lim. Too much liquid in the turbine")
        
            # endregion etat 5

            # region ETAT 4
        self.T_4 = self.T_max



        # Faire une matrice comme au HMW 3 : plus petite 

        
        # self.T_hot_fluid_out_I = self.T_1 + self.T_pinch_evap_I

        # self.dot_m_tot = self.dot_m_ex * self.CP_av(self.T_hot_fluid_in_II, self.T_hot_fluid_out_I,self.p_hot_fluid,self.p_hot_fluid,self.fluid) * (self.T_hot_fluid_in_II - self.T_hot_fluid_out_I) / (self.CP_av(self.T_2,self.T_3,self.p_2,self.p_3,self.fluid) * (self.T_3 - self.T_2) + self.CP_av(self.T_1,self.T_4,self.p_1,self.p_4,self.fluid) * (self.T_4 - self.T_1))
        # print('dot_m_tot === ',self.dot_m_tot,'[kg/s]')





    

        #endregion etats

        #region massflowrates

   
        #endregion massflowrates


        #region Rendements
        self.eta_toten = 0.4
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
   


