


#
#===IMPORT PACKAGES============================================================
#

import numpy as np
from scipy.optimize import fsolve
import CoolProp
import CoolProp.CoolProp as CP
from CoolProp.CoolProp import PropsSI
from scipy.integrate import quad
import math
import sympy as sp
import scipy as scp
from scipy.optimize import differential_evolution
from scipy.optimize import newton_krylov
from scipy.optimize import least_squares
import matplotlib.pyplot as plt
from tabulate import tabulate
import pandas as pd
from sympy import symbols, Eq, solve
from scipy.linalg import lstsq
from scipy.optimize import minimize_scalar


#### Rankine Organic Cycle ####

class ORC(object): 
    
    def __init__(self, inputs,parameters,display):
        """
        Create a steam turbine object.
        """                                       
        self.p_ref              = parameters['p_ref']        # [Pa] reference pressure for exergy computation
        self.T_ref              = parameters['T_ref']        # [K] reference temperature for exergy computation
        self.T_approach_II      = parameters['T_surchauffe_3']# [K] maximum steam temperature inlet 3
        self.T_approach_I       = parameters['T_surchauffe_4'] # [K] maximum steam temperature inlet 4
        self.T_cd_subcool       = parameters['T_cd_subcool'] # [K] condenser cold outlet temperature
        self.T_pinch_ex         = parameters['T_pinch_ex_I']   # [K] pinch temperature at heat exchangers
        self.T_pinch_cd         = parameters['T_pinch_cd']   # [K] pinch temperature at the condenser 
        self.eta_LP_is          = parameters['eta_is_T']    # [-] LP turbine isentropic efficiency (see text book pp. 55)
        self.eta_mec            = parameters['eta_mec']      # [-] shafts bearings mechanical efficiency
        self.eta_pump           = parameters['eta_pump_1']     # [-] internal efficiency of the pump (see text book pp. 53)
        self.fluid_recovery     = parameters['hot_fluid']    # [-] fluid recovery
        self.fluid              = parameters['fluid']        # [-] fluid
        self.fluid_condensor    = parameters['cold_fluid']   # [-] fluid condensor
        self.m_HTF              = parameters['m_dot_HF']        # [kg/s] flow rate of the fluid recovery
        self.T_inlet_rec        = parameters['T_7']     # [K] temperature at the inlet of the fluid recovery
        self.T_outlet_rec       = parameters['T_9']    # [K] temperature at the outlet of the fluid recovery
        self.p_rec              = parameters['p_HF']        # [Pa] pressure of the recovering fluid
        self.T_inlet_cond       = parameters['T_11']    # [K] temperature at the inlet of the condensor
        self.T_outlet_cond      = parameters['T_10']   # [K] temperature at the outlet of the condensor
        self.p_cond             = parameters['p_CF']       # [Pa] pressure of the condensor fluid 
        
        self.display            = display                    # if True, make the plots
   
    
    def evaluate(self):
        
        Dmolar = CP.PropsSI("Dmolar", "T", self.T_ref - 14.99, "P", self.p_ref, self.fluid)
        CP.set_reference_state(self.fluid, self.T_ref - 14.99, Dmolar,0, 0)

        self.h_ref = CP.PropsSI('H', 'P', self.p_ref,'T',self.T_ref, 'H2O')
        self.s_ref = CP.PropsSI('S', 'P', self.p_ref,'T',self.T_ref, 'H2O')
        
        def exergy (s,h) : 
            return (h - self.h_ref - self.T_ref*(s - self.s_ref))
        
        
        self.T_7 = self.T_inlet_rec
        self.p_7 = self.p_rec
        self.p_8 = self.p_7
        self.p_9 = self.p_7

        self.h_7 = CP.PropsSI('H', 'P', self.p_7, 'T', self.T_7, self.fluid_recovery)
        self.s_7 = CP.PropsSI('S', 'P', self.p_7, 'T', self.T_7, self.fluid_recovery)
        
        ##### HYPOTHÈSES #####
        
        self.T_9 = self.T_outlet_rec
        
        ####### Valeurs normées #######
        
        self.p_11 = self.p_cond
        self.p_10 = self.p_11
        self.T_11 = self.T_inlet_cond
        self.h_11 = CP.PropsSI('H', 'P', self.p_11, 'T', self.T_11, self.fluid_condensor)
        self.s_11 = CP.PropsSI('S', 'P', self.p_11, 'T', self.T_11, self.fluid_condensor)
        
        self.T_10 = self.T_outlet_cond
        self.h_10 = CP.PropsSI('H', 'P', self.p_11, 'T', self.T_10, self.fluid_condensor)
        self.s_10 = CP.PropsSI('S', 'P', self.p_11, 'T', self.T_10, self.fluid_condensor)
  
        ####### #######
        
        self.T_8 = (self.T_inlet_rec - self.T_outlet_rec)/2 + self.T_outlet_rec
        self.p_8 = self.p_7
        self.p_9 = self.p_7
        self.h_8 = CP.PropsSI('H', 'P', self.p_7, 'T', self.T_8, self.fluid_recovery)
        self.s_8 = CP.PropsSI('S', 'P', self.p_7, 'T', self.T_8, self.fluid_recovery)
        self.h_9 = CP.PropsSI('H', 'P', self.p_7, 'T', self.T_9, self.fluid_recovery)
        self.s_9 = CP.PropsSI('S', 'P', self.p_7, 'T', self.T_9, self.fluid_recovery)
        
        self.T_3 = self.T_7 - self.T_approach_II
        self.T_4 = self.T_8 - self.T_approach_I
        
        def equationsTOTALES(variables):
            print(variables)
            T_1, T_2, T_5 = variables

            ####### INTERLUDE OU INTERLUX p_2 ########
            
            initial_guess_equations_p_2 = [0.8e6]
            def equations_p_2(variable):
                p_2_i = variable
                
                T_2_p_i   = CP.PropsSI('T', 'Q', 0,'P', p_2_i, self.fluid)
                h_2_p_i   = CP.PropsSI('H', 'Q', 0,'P', p_2_i, self.fluid)
                
                h_2_i = CP.PropsSI('H', 'T', T_2,'P', p_2_i, self.fluid)
                h_3_i = CP.PropsSI('H', 'T', self.T_3,'P', p_2_i, self.fluid)
                m_WFII_i = self.m_HTF * (self.h_7 - self.h_8)/(h_3_i - h_2_i)
                
                h_PII_i = self.h_7 - (m_WFII_i/self.m_HTF) * (h_3_i - h_2_p_i)
                T_PII_i = CP.PropsSI('T', 'P', self.p_7, 'H', h_PII_i, self.fluid_recovery)
                
                Eq4 = T_PII_i - T_2_p_i - self.T_pinch_ex
                return Eq4
            
            p_2 = fsolve(equations_p_2, initial_guess_equations_p_2)[0]
            ##### Evaporator II #####
            h_2_p =   CP.PropsSI('H', 'Q', 0,'P', p_2, self.fluid) 
            h_2 = CP.PropsSI('H', 'T', T_2,'P', p_2, self.fluid)
            h_3 = CP.PropsSI('H', 'T', self.T_3,'P', p_2, self.fluid)
            m_WFII = self.m_HTF * (self.h_7 - self.h_8)/(h_3 - h_2)
            
            h_PII = self.h_7 - (m_WFII/self.m_HTF) * (h_3 - h_2_p)
            T_PII = CP.PropsSI('T', 'P', self.p_7, 'H', h_PII, self.fluid_recovery)
            
            ####### INTERLUDE OU INTERLUX p_1 ########
            
            initial_guess_equations_p_1 = 5e5
            def equations_p_1(variable):
                p_1_i = variable
                
                T_1_p_i   = CP.PropsSI('T', 'Q', 0,'P', p_1_i, self.fluid)
                h_1_p_i   = CP.PropsSI('H', 'Q', 0,'P', p_1_i, self.fluid)
                
                h_1_i = CP.PropsSI('H', 'T', T_1,'P', p_1_i, self.fluid)
                h_4_i = CP.PropsSI('H', 'T', self.T_4,'P', p_1_i, self.fluid)
                m_WFI_i = self.m_HTF * (self.h_8 - self.h_9)/(h_4_i - h_1_i)
                
                h_PI_i = self.h_8 - (m_WFI_i/self.m_HTF) * (h_4_i - h_1_p_i)
                T_PI_i = CP.PropsSI('T', 'P', self.p_7, 'H', h_PI_i, self.fluid_recovery)
                
                Eq5 = T_PI_i - T_1_p_i - self.T_pinch_ex
                return Eq5
            
            p_1 = fsolve(equations_p_1, initial_guess_equations_p_1)[0]
            print(p_1)
                        
            ##### Evaporator I #####
            h_1_p   = CP.PropsSI('H', 'Q', 0,'P', p_1, self.fluid) 
            h_1 = CP.PropsSI('H', 'T', T_1,'P', p_1, self.fluid)
            h_4 = CP.PropsSI('H', 'T', self.T_4,'P', p_1, self.fluid)
            m_WFI = self.m_HTF * (self.h_8 - self.h_9)/(h_4 - h_1)
            
            h_PI = self.h_8 - (m_WFI/self.m_HTF) * (h_4 - h_1_p)
            T_PI = CP.PropsSI('T', 'P', self.p_7, 'H', h_PI, self.fluid_recovery)
            
            m_WFTOT = m_WFI + m_WFII
            
            ####### INTERLUDE OU INTERLUX p_6 ########
            
            initial_guess_equations_T_6 = 28 + 273.15
            def equations_p_6(variable):                
                T_6_i = variable 
                
                T_sat_6_i = T_6_i[0]  + self.T_cd_subcool
                p_6_i =  CP.PropsSI('P', 'Q', 1,'T', T_sat_6_i, self.fluid)
                T_5_p_p_i = CP.PropsSI('T', 'Q', 1,'P', p_6_i, self.fluid)
                h_5_p_p_i = CP.PropsSI('H', 'Q', 1,'P', p_6_i, self.fluid)
                
                h_6_i = CP.PropsSI('H', 'P', p_6_i, 'T', T_6_i, self.fluid)
                h_5_i = CP.PropsSI('H', 'T', T_5,'P', p_6_i, self.fluid)
                m_C_i = m_WFTOT * (h_5_i - h_6_i)/ (self.h_10 - self.h_11)
                
                h_PC_i = (m_WFTOT/m_C_i) * (h_5_p_p_i - h_6_i) + self.h_11
                T_PC_i = CP.PropsSI('T', 'P', self.p_11, 'H', h_PC_i, self.fluid_condensor)
                
                Eq6 = T_5_p_p_i - T_PC_i - self.T_pinch_cd
                return Eq6
            
            T_6 = fsolve(equations_p_6, initial_guess_equations_T_6)[0]
            p_6 = CP.PropsSI('P', 'Q', 1,'T',  T_6 + self.T_cd_subcool , self.fluid)
            ##### Condensor #####
            
            h_5_p_p = CP.PropsSI('H', 'Q', 1,'P', p_6, self.fluid)  
            h_6 = CP.PropsSI('H', 'T', T_6,'P', p_6, self.fluid)
            h_5 = CP.PropsSI('H', 'T', T_5,'P', p_6, self.fluid)
            m_C = m_WFTOT * (h_5 - h_6)/ (self.h_10 - self.h_11)
            
            h_PC = (m_WFTOT/m_C) * (h_5_p_p - h_6) + self.h_11
            T_PC = CP.PropsSI('T', 'P', self.p_11, 'H', h_PC, self.fluid_condensor)
            
            ##### PUMP I #####
            def vI(P):
                return 1/CP.PropsSI('D', 'H',(h_6 + h_1)/2,'P',P,self.fluid)
            
            vdp_variableI, _ = quad(vI, p_6, p_1)
            
            Eq1 = h_1 - h_6 - (vdp_variableI/self.eta_pump)
            
            ##### PUMP II #####
            def vII(P):
                return 1/CP.PropsSI('D', 'H',(h_2 + h_1)/2,'P',P,self.fluid)
            
            vdp_variableII, _ = quad(vII, p_1, p_2)
            
            Eq2 = h_2 - h_1 - (vdp_variableII/self.eta_pump)
            
            ##### TURBINE 3-5 #####
            
            s_3 = CP.PropsSI('S', 'P', p_2, 'T', self.T_3, self.fluid)
            h_3_inter_is = CP.PropsSI('H', 'P', p_6, 'S', s_3, self.fluid)
            h_3_inter = h_3 - self.eta_LP_is * (h_3 - h_3_inter_is)
            
            ##### TURBINE 4-5 #####
            
            s_4 = CP.PropsSI('S', 'P', p_1, 'T', self.T_4, self.fluid)
            h_4_inter_is = CP.PropsSI('H', 'P', p_6, 'S', s_4, self.fluid)
            h_4_inter = h_4 - self.eta_LP_is * (h_4 - h_4_inter_is)
            
            ##### MELANGE #####
            
            Eq3 = m_WFI * h_3_inter + m_WFII * h_4_inter - m_WFTOT * h_5
            
    
            return Eq1, Eq2, Eq3

        initial_guess_equationsTOTALES = [25.3 + 273.15, 25.6 + 273.15, 100 + 273.15]


        self.solution = least_squares(
            equationsTOTALES,
            initial_guess_equationsTOTALES,
            bounds=([273.15 + 20, 273.15 + 20, self.T_11],
                    [273.15 + 50, 273.15 + 50, 273.15 + 120])
        )
        print(self.solution.fun)
        self.T_1 = self.solution.x[0]
        self.T_2 = self.solution.x[1]
        self.T_5 = self.solution.x[2]
        
        initial_guess_equations_p_2 = [8.59406195e5]
        def equations_p_2(variable):
            p_2 = variable
            
            T_2_p   = CP.PropsSI('T', 'Q', 0,'P', p_2, self.fluid)
            h_2_p   = CP.PropsSI('H', 'Q', 0,'P', p_2, self.fluid)
            
            h_2 = CP.PropsSI('H', 'T', self.T_2,'P', p_2, self.fluid)
            h_3 = CP.PropsSI('H', 'T', self.T_3,'P', p_2, self.fluid)
            m_WFII = self.m_HTF * (self.h_7 - self.h_8)/(h_3 - h_2)
            
            h_PII = self.h_7 - (m_WFII/self.m_HTF) * (h_3 - h_2_p)
            T_PII = CP.PropsSI('T', 'P', self.p_7, 'H', h_PII, self.fluid_recovery)
            
            Eq7 = T_PII - T_2_p - self.T_pinch_ex
            return Eq7

        self.p_2 = fsolve(equations_p_2, initial_guess_equations_p_2)[0]
        self.p_3 = self.p_2

        initial_guess_equations_p_1 = [5e5]
        def equations_p_1(variable):
            p_1 = variable
            
            T_1_p   = CP.PropsSI('T', 'Q', 0,'P', p_1, self.fluid)
            h_1_p   = CP.PropsSI('H', 'Q', 0,'P', p_1, self.fluid)
            
            h_1 = CP.PropsSI('H', 'T', self.T_1,'P', p_1, self.fluid)
            h_4 = CP.PropsSI('H', 'T', self.T_4,'P', p_1, self.fluid)
            m_WFI = self.m_HTF * (self.h_8 - self.h_9)/(h_4 - h_1)
            
            h_PI = self.h_8 - (m_WFI/self.m_HTF) * (h_4 - h_1_p)
            T_PI = CP.PropsSI('T', 'P', self.p_7, 'H', h_PI, self.fluid_recovery)
            
            Eq8 = T_PI - T_1_p - self.T_pinch_ex
            return Eq8
        
        self.p_1 = fsolve(equations_p_1, initial_guess_equations_p_1)[0]
        self.p_4 = self.p_1

        ##### Evaporator I #####
        self.T_2_p   = CP.PropsSI('T', 'Q', 0,'P', self.p_2, self.fluid)
        self.T_2_p_p = CP.PropsSI('T', 'Q', 1,'P', self.p_2, self.fluid)
        self.h_2_p_p = CP.PropsSI('H', 'Q', 1,'P', self.p_2, self.fluid)
        self.h_2_p   = CP.PropsSI('H', 'Q', 0,'P', self.p_2, self.fluid)
        
        self.h_2 = CP.PropsSI('H', 'T', self.T_2,'P', self.p_2, self.fluid)
        self.h_3 = CP.PropsSI('H', 'T', self.T_3,'P', self.p_2, self.fluid)
        self.m_WFII = self.m_HTF * (self.h_7 - self.h_8)/(self.h_3 - self.h_2)

        self.h_PII = self.h_7 - (self.m_WFII/self.m_HTF) * (self.h_3 - self.h_2_p)
        self.T_PII = CP.PropsSI('T', 'P', self.p_7, 'H', self.h_PII, self.fluid_recovery)

        ##### Evaporator I #####
        self.T_1_p   = CP.PropsSI('T', 'Q', 0,'P', self.p_1, self.fluid)
        self.T_1_p_p = CP.PropsSI('T', 'Q', 1,'P', self.p_1, self.fluid)
        self.h_1_p_p = CP.PropsSI('H', 'Q', 1,'P', self.p_1, self.fluid)
        self.h_1_p   = CP.PropsSI('H', 'Q', 0,'P', self.p_1, self.fluid)
        
        self.h_1 = CP.PropsSI('H', 'T', self.T_1,'P', self.p_1, self.fluid)
        self.h_4 = CP.PropsSI('H', 'T', self.T_4,'P', self.p_1, self.fluid)
        self.m_WFI = self.m_HTF * (self.h_8 - self.h_9)/(self.h_4 - self.h_1)

        self.h_PI = self.h_8 - (self.m_WFI/self.m_HTF) * (self.h_4 - self.h_1_p)
        self.T_PI = CP.PropsSI('T', 'P', self.p_7, 'H', self.h_PI, self.fluid_recovery)

        self.m_WFTOT = self.m_WFI + self.m_WFII

        initial_guess_equations_T_6 = [25 + 273.15]
        def equations_T_6(variable):
            T_6 = variable
            T_sat_6 = T_6[0] + self.T_cd_subcool
            p_6 = CP.PropsSI('P', 'Q', 1,'T',T_sat_6,self.fluid)
            T_5_p_p = CP.PropsSI('T', 'Q', 1,'P', p_6,self.fluid)
            h_5_p_p = CP.PropsSI('H', 'Q', 1,'P', p_6,self.fluid)
            
            h_6 = CP.PropsSI('H', 'T', T_6,'P', p_6, self.fluid)
            h_5 = CP.PropsSI('H', 'T', self.T_5,'P', p_6, self.fluid)
            m_C = self.m_WFTOT * (h_5 - h_6)/ (self.h_10 - self.h_11)
            
            h_PC = (self.m_WFTOT/m_C) * (h_5_p_p - h_6) + self.h_11
            T_PC = CP.PropsSI('T', 'P', self.p_11, 'H', h_PC, self.fluid_condensor)
            
            Eq9 = T_5_p_p - T_PC - self.T_pinch_cd
            return Eq9

        self.T_6 = fsolve(equations_T_6, initial_guess_equations_T_6)[0]
        self.T_sat_6 = self.T_6 + self.T_cd_subcool
        self.p_6 = CP.PropsSI('P', 'Q', 1,'T', self.T_sat_6,self.fluid)
        self.p_5 = self.p_6

        print("P1 :", self.p_1)
        print("P2 :", self.p_2)
        print("P6 :", self.p_6)

        ##### Condensor #####
        self.T_5_p   = CP.PropsSI('T', 'Q', 0,'P', self.p_6, self.fluid)
        self.T_5_p_p = CP.PropsSI('T', 'Q', 1,'P', self.p_6, self.fluid)  
        self.h_5_p_p = CP.PropsSI('H', 'Q', 1,'P', self.p_6, self.fluid)
        self.h_5_p   = CP.PropsSI('H', 'Q', 0,'P', self.p_6, self.fluid)
        
        self.h_6 = CP.PropsSI('H', 'T', self.T_6,'P', self.p_6, self.fluid)
        self.h_5 = CP.PropsSI('H', 'T', self.T_5,'P', self.p_6, self.fluid)
        self.m_C = self.m_WFTOT * (self.h_5 - self.h_6)/ (self.h_10 - self.h_11)

        self.h_PC = (self.m_WFTOT/self.m_C) * (self.h_5_p_p - self.h_6) + self.h_11
        self.T_PC = CP.PropsSI('T', 'P', self.p_11, 'H', self.h_PC, self.fluid_condensor)
        
        
        self.s_1 = CP.PropsSI('S', 'T',self.T_1,'P',self.p_1, self.fluid)
        self.e_1 = exergy(self.s_1, self.h_1)
        
        self.s_2 = CP.PropsSI('S', 'T',self.T_2,'P',self.p_2, self.fluid)
        self.e_2 = exergy(self.s_2, self.h_2)
        
        self.s_3 = CP.PropsSI('S', 'T',self.T_3,'P',self.p_3, self.fluid)
        self.h_3_inter_is = CP.PropsSI('H', 'P', self.p_6, 'S', self.s_3, self.fluid) 
        self.h_3_inter = self.h_3 - self.eta_LP_is * (self.h_3 - self.h_3_inter_is)
        self.s_3_inter = CP.PropsSI('S', 'P', self.p_6, 'H', self.h_3_inter, self.fluid) 
        
        self.e_3 = exergy(self.s_3, self.h_3)
        self.e_3_inter = exergy(self.s_3_inter,self.h_3_inter)
        
        self.s_4 = CP.PropsSI('S', 'T',self.T_4,'P',self.p_4, self.fluid)
        self.h_4_inter_is = CP.PropsSI('H', 'P', self.p_6, 'S', self.s_4, self.fluid) 
        self.h_4_inter = self.h_4 - self.eta_LP_is * (self.h_4 - self.h_4_inter_is)
        self.s_4_inter = CP.PropsSI('S', 'P', self.p_6, 'H', self.h_4_inter, self.fluid) 
        
        self.e_4_inter = exergy(self.s_4_inter,self.h_4_inter)
        self.e_4 = exergy(self.s_4, self.h_4)
        
        self.s_5 = CP.PropsSI('S', 'T',self.T_5,'P',self.p_5, self.fluid)
        self.e_5 = exergy(self.s_5, self.h_5)
        
        self.s_6 = CP.PropsSI('S', 'T',self.T_6,'P',self.p_6, self.fluid)
        self.e_6 = exergy(self.s_6, self.h_6)
        
        self.s_7 = CP.PropsSI('S', 'T',self.T_7,'P',self.p_7, self.fluid)
        self.e_7 = exergy(self.s_7, self.h_7)
        
        self.s_8 = CP.PropsSI('S', 'T',T_8,'P',self.p_8, self.fluid)
        self.e_8 = exergy(self.s_8, self.h_8)
        print("hihi :" ,T_8, self.e_8)
        
        self.s_9 = CP.PropsSI('S', 'T',self.T_9,'P',self.p_9, self.fluid)
        self.e_9 = exergy(self.s_9, self.h_9)
        
        self.s_10 = CP.PropsSI('S', 'T',self.T_10,'P',self.p_10, self.fluid)
        self.e_10 = exergy(self.s_10, self.h_10)
        
        self.s_11 = CP.PropsSI('S', 'T',self.T_11,'P',self.p_11, self.fluid)
        self.e_11 = exergy(self.s_11, self.h_11)
        
        self.W_3_5 = (self.h_3 - self.h_3_inter) * self.m_WFII
        self.W_4_5 = (self.h_4 - self.h_4_inter) * self.m_WFI
        
        self.W_6_1 = (self.h_1 - self.h_6) * self.m_WFTOT
        self.W_1_2 = (self.h_2 - self.h_1) * self.m_WFII
        
        self.W_tot = self.W_3_5 + self.W_4_5 - self.W_6_1 - self.W_1_2
        
        self.P_e = self.W_tot * self.eta_mec
        
        # self.Q_7_8 = (self.h_7 - self.h_8) * self.m_HTF
        # self.Q_8_9 = (self.h_8 - self.h_9) * self.m_HTF
        # self.Q_tot = self.Q_7_8 + self.Q_8_9
        
        # self.e_7_8 = (self.e_7 - self.e_8) * self.m_HTF
        # self.e_8_9 = (self.e_8 - self.e_9) * self.m_HTF
        # self.e_comb = self.e_7_8 + self.e_8_9
        
        self.Q_3_2 = (self.h_3 - self.h_2) * self.m_WFII
        self.Q_4_1 = (self.h_4 - self.h_1) * self.m_WFI
        self.Q_tot = self.Q_3_2 + self.Q_4_1
        
        self.e_3_2 = (self.e_3 - self.e_2) * self.m_WFII
        self.e_4_1 = (self.e_4 - self.e_1) * self.m_WFI
        self.e_comb = self.e_3_2 + self.e_4_1
        
        
        EVAP_I_HOT   = (self.e_8 - self.e_9) * self.m_HTF
        EVAP_I_COLD  = (self.e_4 - self.e_1) * self.m_WFI
        self.loss_transex_I = EVAP_I_HOT - EVAP_I_COLD
        self.eta_transex_I = EVAP_I_COLD / EVAP_I_HOT
        
        EVAP_II_HOT  = (self.e_7 - self.e_8) * self.m_HTF
        EVAP_II_COLD = (self.e_3 - self.e_2) * self.m_WFII
        self.loss_transex_II = EVAP_II_HOT - EVAP_II_COLD
        self.eta_transex_II = EVAP_II_COLD / EVAP_II_HOT
        
        self.loss_transex = self.loss_transex_I + self.loss_transex_II
        self.eta_transex = (self.eta_transex_I + self.eta_transex_II) / 2
        
        
        self.eta_condex = (self.m_WFTOT*(self.e_5 - self.e_6))/ (self.m_C *(self.e_10 - self.e_11))
        self.eta_rotex = (self.m_WFI * (self.h_4 - self.h_4_inter) + self.m_WFII * (self.h_3 - self.h_3_inter))/ (self.m_WFI * (self.e_4 - self.e_4_inter) + self.m_WFII * (self.e_3 - self.e_3_inter))
        
        self.loss_mec = (1 - self.eta_mec) * self.W_tot
        
        self.loss_conden = (self.h_5 - self.h_6) * self.m_WFTOT 
        self.loss_condex = (self.e_5 - self.e_6) * self.m_WFTOT
        
        self.loss_rotex = (1 - self.eta_rotex) * self.W_tot
        
        self.loss_exhaust = self.m_HTF*self.e_9
        print(self.m_HTF*self.e_7/1e6)
        
        
        
        
        self.p           = self.p_1, self.p_2, self.p_3,self.p_4,self.p_5,self.p_6,self.p_7,self.p_8,self.p_9,self.p_10,self.p_11 # [Pa]     tuple containing the pressure at each state
        self.T           = self.T_1, self.T_2, self.T_3,self.T_4,self.T_5,self.T_6,self.T_7,self.T_8,self.T_9,self.T_10,self.T_11 # [K]      temperature at each state
        self.s           = self.s_1, self.s_2, self.s_3,self.s_4,self.s_5,self.s_6,self.s_7,self.s_8,self.s_9,self.s_10,self.s_11 # [J/kg/K] entropy at each state
        self.h           = self.h_1, self.h_2, self.h_3,self.h_4,self.h_5,self.h_6,self.h_7,self.h_8,self.h_9,self.h_10,self.h_11 # [J/kg]   enthalpy at each state
        self.e           = self.e_1, self.e_2, self.e_3,self.e_4,self.e_5,self.e_6,self.e_7,self.e_8,self.e_9,self.e_10,self.e_11 # [J/kg]   exergy at each state (use ref conditions)
        self.DAT         = self.p,self.T,self.s,self.h,self.e
        
        pd.set_option("display.max_rows", None)
        pd.set_option("display.max_columns", None)
        pd.set_option("display.width", 700) 

        States = [
        "State 1", "State 2", "State 3", "State 4", "State 5", "State 6", "State 7", "State 8", 
        "State 9", "State 10", "State 11"]
        
        
        p = np.array(self.p)/1000
        T = np.array(self.T)-273.15
        s = np.array(self.s)/1000
        h = np.array(self.h)/1000
        e = np.array(self.e)/1000
        
        data = {
            "Etat": States,
            "Pression [kPa]": p,
            "Température [K]": T,
            "Entropie [kJ/kg/K]": s,
            "Enthalpie [kJ/kg]": h,
            "Exergie [kJ/kg]": e
        }   
        df = pd.DataFrame(data)
        #print(df)
        
        self.ETA         = self.eta_cyclen,self.eta_cyclex,self.eta_condex,self.eta_transex,self.eta_rotex
        # -> see text book pp. 53-94
        #      o eta_cyclen    [-]      cycle energy efficiency
        #      o eta_cyclex    [-]      cycle exergy efficiency
        #      o eta_condex    [-]      condenser exergy efficiency
        #      o eta_transex   [-]      bleedings heat exchangers overall exergy efficiency
        #      o eta_rotex     [-]      pumps and turbines exergy efficiency
        # Energy losses -------------------------------------------------------
        self.DATEN       = self.loss_mec,self.loss_conden
        #      o loss_mec      [W]      mechanical energy losses
        #      o loss_conden   [W]      condenser energy losses
        # Exergy losses -------------------------------------------------------
        self.DATEX       = self.loss_mec,self.loss_rotex,self.loss_transex,self.loss_condex, self.loss_exhaust
        #      o loss_mec      [W]      mechanical energy losses
        #      o loss_rotex    [W]      pumps and turbines exergy losses
        #      o loss_transex  [W]      bleedings heat exchangers overall exergy losses
        #      o loss_condex   [W]      condenser exergy losses
       
        
        
        if (self.display) :
            
            plt.figure(1)
            labels = ['Mechanical losses', 'Condensor losses', 'Electrical Power']
            sizes = [self.loss_mec, self.loss_conden, self.P_e]
            sizes_in_mw = [size / 1e6 for size in sizes]
            labels_with_values = [f'{label}\n{size:.3f} MW' for label, size in zip(labels, sizes_in_mw)]
            plt.pie(sizes, labels=labels_with_values, autopct='%1.1f%%', startangle=140)
            plt.axis('equal')
            self.fig_pie_en = plt.figure(1)
            
            plt.figure(2)
            labels = ['Mechanical losses','Condensor losses',  'Effective energy','Heat exchanger losses','Turbine and compressor irreversibilities',"Exhaust losses"]
            sizes = [ self.loss_mec,self.loss_condex,self.P_e, self.loss_transex, self.loss_rotex, self.loss_exhaust]
            sizes_in_mw = [size / 1e6 for size in sizes]
            
            labels_with_values = [f'{label}\n{size:.3f} MW' for label, size in zip(labels, sizes_in_mw)]
            plt.pie(sizes, labels=labels_with_values, autopct='%1.1f%%', startangle=180)
            plt.axis('equal')
            self.fig_pie_ex= plt.figure(2)
        
            self.Q_ECON_II = ( self.m_WFII * (self.h_2_p - self.h_2) ) / ( self.m_HTF * (self.h_7 - self.h_8) )
            #print("Chaleur échangée econ II:", Q_ECON_II)
            self.Q_SUPERH_II = ( self.m_WFII * (self.h_3 - self.h_2_p_p) ) / ( self.m_HTF * (self.h_7 - self.h_8) )
            #print("Chaleur échangée superII:",Q_SUPERH_II)
            
            plt.figure(3)
            plt.plot([0,100], [self.T_7 - 273.15, self.T_8 - 273.15], color="blue", marker="o")
            plt.plot([100 *(1 - self.Q_ECON_II),100], [self.T_2_p- 273.15, self.T_2- 273.15], color="red", marker="o")
            plt.plot([100 * (self.Q_SUPERH_II), 100 * (1 - self.Q_ECON_II)], [self.T_2_p_p- 273.15, self.T_2_p- 273.15], color="red", marker="o")
            plt.plot([0, self.Q_SUPERH_II * 100], [self.T_3- 273.15, self.T_2_p_p - 273.15], color="red", marker="o")
            plt.grid()
            plt.show()
    
            self.Q_ECON_I = ( self.m_WFI * (self.h_1_p - self.h_1) ) / ( self.m_HTF * (self.h_8 - self.h_9) )
            #print("Chaleur échangée econ I:",Q_ECON_I)
            self.Q_SUPERH_I = ( self.m_WFI * (self.h_4 - self.h_1_p_p) ) / ( self.m_HTF * (self.h_8 - self.h_9) )
            #print("Chaleur échangée super I:",Q_SUPERH_I)
            plt.figure(4)
            plt.plot([0,100], [self.T_8 - 273.15, self.T_9 - 273.15], color="blue", marker="o")
            plt.plot([100 *(1 - self.Q_ECON_I),100], [self.T_1_p- 273.15, self.T_1- 273.15], color="red", marker="o")
            plt.plot([100 * (self.Q_SUPERH_I), 100 * (1 - self.Q_ECON_I)], [self.T_1_p_p- 273.15, self.T_1_p- 273.15], color="red", marker="o")
            plt.plot([0, self.Q_SUPERH_I * 100], [self.T_4- 273.15, self.T_1_p_p - 273.15], color="red", marker="o")
            plt.grid()
            plt.show()
    
            self.Q_ECON_C = ( self.m_WFTOT * (self.h_5_p - self.h_6) ) / ( self.m_C * (self.h_10 - self.h_11) )
            #print("Chaleur échangée econ C:",Q_ECON_C)
            self.Q_SUPERH_C = ( self.m_WFTOT * (self.h_5 - self.h_5_p_p) ) / ( self.m_C * (self.h_10 - self.h_11) )
            #print("Chaleur échangée super C:",Q_SUPERH_C)
            #print(self.m_C, self.m_WFI, self.m_WFII, self.m_WFTOT)
            print(self.solution.fun)
            plt.figure(5)
            plt.plot([0,100], [self.T_11 - 273.15, self.T_10 - 273.15], color="blue", marker="o")
            plt.plot([100 *(1 - self.Q_SUPERH_C),100], [self.T_5_p_p- 273.15, self.T_5- 273.15], color="red", marker="o")
            plt.plot([100 * (self.Q_ECON_C), 100 * (1 - self.Q_SUPERH_C)], [self.T_5_p- 273.15, self.T_5_p_p- 273.15], color="red", marker="o")
            plt.plot([0, self.Q_ECON_C * 100], [self.T_6- 273.15, self.T_5_p - 273.15], color="red", marker="o")
            plt.grid()
            plt.show()
            
            num = 50
            
            
            T_range = np.linspace(135, 427.01, 200)  # Température entre 250 K (-23°C) et 430 K (157°C)
            
            s_liquid = []
            s_vapor = []
            T_liquid = []
            T_vapor = []
            
            # Calculer les propriétés pour chaque température
            for T in T_range:
                try:
                    # Entropie saturée liquide (Q = 0) et vapeur (Q = 1)
                    s_liq = CP.PropsSI('S', 'T', T, 'Q', 0, self.fluid)   # kJ/kg·K
                    s_vap = CP.PropsSI('S', 'T', T, 'Q', 1, self.fluid)   # kJ/kg·K
            
                    # Enregistrer les valeurs
                    s_liquid.append(s_liq)
                    s_vapor.append(s_vap)
                    T_liquid.append(T)  # Convertir T en °C
                    T_vapor.append(T )
                except:
                    # Éviter les erreurs près du point critique
                    continue
                
            T_1_4 = np.linspace(self.T_1,self.T_4,num)
            s_1_4 = CP.PropsSI('S', 'T',T_1_4, 'P',self.p_1, self.fluid)
            
            T_2_3 = np.linspace(self.T_2,self.T_3,num)
            s_2_3 = CP.PropsSI('S', 'T',T_2_3, 'P',self.p_2, self.fluid)
            
            T_6_5 = np.linspace(self.T_5,self.T_6,num)
            s_6_5 = CP.PropsSI('S', 'T',T_6_5, 'P',self.p_6, self.fluid)
            
            T_8_7 = np.linspace(self.T_8,self.T_7,num)
            s_8_7 = CP.PropsSI('S', 'T',T_8_7, 'P',self.p_7, self.fluid_recovery)
            
            T_8_9 = np.linspace(self.T_7,self.T_9,num)
            s_8_9 = CP.PropsSI('S', 'T',T_8_9, 'P',self.p_7, self.fluid_recovery)
            
     
            p_2_6 = np.linspace(self.p_3, self.p_5,num)
            
            h_3_inter_is = CP.PropsSI('H', 'P', p_2_6, 'S', self.s_3, self.fluid)
            h_3_inter = self.h_3 - self.eta_LP_is * (self.h_3 - h_3_inter_is)
            
            T_3_3_p = CP.PropsSI('T', 'P', p_2_6, 'H', h_3_inter, self.fluid)
            s_3_3_p = CP.PropsSI('S', 'P', p_2_6, 'H', h_3_inter, self.fluid)
            
            ##### TURBINE 4-5 #####
            
            p_1_6 = np.linspace(self.p_1, self.p_6,num)
            h_4_inter_is = CP.PropsSI('H', 'P', p_1_6, 'S', self.s_4, self.fluid)
            h_4_inter = self.h_4 - self.eta_LP_is * (self.h_4 - h_4_inter_is)
            
            T_4_4_p = CP.PropsSI('T', 'P', p_1_6, 'H', h_4_inter, self.fluid)
            s_4_4_p = CP.PropsSI('S', 'P', p_1_6, 'H', h_4_inter, self.fluid)
            
            
            ##### MELANGE #####
            
            h_5_v =(self.m_WFI * h_3_inter + self.m_WFII * h_4_inter)/self.m_WFTOT 
            
            plt.figure(6)
            plt.plot(s_liquid, T_liquid, label='Saturated Liquid Line', color='lightblue')
            plt.plot(s_vapor, T_vapor, label='Saturated Vapor Line', color='lightblue')
            plt.plot(s_1_4,T_1_4)
            plt.plot(s_2_3,T_2_3)
            plt.plot(s_6_5,T_6_5)
            plt.plot(s_3_3_p,T_3_3_p)
            plt.plot(s_4_4_p,T_4_4_p)
            plt.grid()
            plt.show()
           
            
          
            
          


        