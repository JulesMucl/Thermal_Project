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


        # Doute de la méthodologie - Voir TODO : moi aussi je suis pas sure du pinch.
        self.T_pinch_evap_I   = parameters['T_pinch_evap_I']

        
        
       

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

            # region ETAT 6
        self.T_6 = self.T_cold_fluid_in + 1 #self.T_pinch_cd()  # Pas sure !!!!
        self.p_6 = PropsSI("P","T",self.T_6,"Q",self.x6,self.fluid)
        self.h_6 = PropsSI("H","T",self.T_6,"P",self.p_6,self.fluid)
        self.s_6 = PropsSI("S","T",self.T_6,"P",self.p_6,self.fluid)
        self.e_6 = self.exergie(self.h_6,self.s_6)
        #endregion etat6


        # Pour moi on peut pas définir les rapports d epression dans les pompes. 
        #     #region ETAT 1
        # self.p_1 = self.p_6 * self.r_pump_1
        # self.pump_1 = Pump1(self)  # self.p_in, self.p_out, self.T_in, self.eta_pump = params
        # self.T_1 = self.pump_1.evaluate_T_out()
        # self.h_1 = PropsSI("H","T",self.T_1,"P",self.p_1,self.fluid)
        # self.s_1 = PropsSI("S","T",self.T_1,"P",self.p_1,self.fluid)
        # self.e_1 = self.exergie(self.h_1,self.s_1)
        #     #endregion etat 1

        #     #region ETAT 2
        # self.p_2 = self.p_1 * self.r_pump_2
        # self.pump_2 = Pump2(self)
        # self.T_2 = self.pump_2.evaluate_T_out()
        # self.h_2 = PropsSI("H","T",self.T_2,"P",self.p_2,self.fluid)
        # self.s_2 = PropsSI("S","T",self.T_2,"P",self.p_2,self.fluid)
        # self.e_2 = self.exergie(self.h_2,self.s_2)

            #endregion etat2

            # region ETAT 3
        self.T_3 = self.T_max
        self.h_3 = PropsSI("H","T",self.T_3,"P",self.p_3,self.fluid)
        self.s_3 = PropsSI("S","T",self.T_3,"P",self.p_3,self.fluid)
        self.e_3 = self.exergie(self.h_3,self.s_3)
        self.x_3 = PropsSI("Q","T",self.T_3,"P",self.p_3,self.fluid)

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

        # Déterminé via un bilan énergétique (conservation énergie)
        self.m_dot_I = self.dot_m_ex * (self.h_8 - self.h_9)/(self.h_4 - self.h_1) 
        self.m_dot_II = self.dot_m_ex * (self.h_7 - self.h_8)/(self.h_3 - self.h_2)
        self.m_dot_tot = self.m_dot_I + self.m_dot_II
        self.m_dot_CF = self.m_dot_tot*(self.h_5 - self.h_6)/(self.h_10 - self.h_11)


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
   


