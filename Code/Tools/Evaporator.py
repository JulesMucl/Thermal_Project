from CoolProp.CoolProp import PropsSI
from scipy.optimize import fsolve

class Evaporator_I(object):

    def __init__(self, orc_instance):

        self.orc            = orc_instance

        self.hot_fluid      = self.orc.hot_fluid
        self.dot_m_ex       = self.orc.dot_m_ex
        self.p_7           = self.orc.p_7

        self.T_9            = self.orc.T_9
        self.T_8            = self.orc.T_8
        self.T_7            = self.orc.T_7
        self.T_10           = self.orc.T_10
        self.T_11           = self.orc.T_11

        self.T_pinch_ex_I   = self.orc.T_pinch_ex_I
        self.T_pinch_ex_II  = self.orc.T_pinch_ex_II

        self.fluid          = self.orc.fluid
        self.dot_m_I        = self.orc.dot_m_I
        self.dot_m_II       = self.orc.dot_m_II

        self.p_3_guess      = self.orc.p_3_guess

        self.ratio_I = self.dot_m_ex / self.dot_m_I

    def H_3(self, p_3_guess):
        """

        Trouver h_3 = h_hs_ex grace à la balance energétique

        """
        
        h_hs_su = PropsSI('H', 'P', p_3_guess, 'Q', 0 , self.fluid)  # Après la pompe c'est encore full liquide, non ?
       
        h_cs_su = PropsSI('H', 'P', self.p_7, 'T', self.T_7, self.hot_fluid)  
        h_cs_ex = PropsSI('H', 'P', self.p_7, 'T', self.T_8, self.hot_fluid)  

        h_hs_ex = h_hs_su - self.ratio_I * ( h_cs_ex - h_cs_su )

        return h_hs_ex
    


    def Etat_i(self,p_3_guess):

        h_hs_ex = self.H_3(p_3_guess)
        h_cs_su = PropsSI('H', 'P', self.p_7, 'T', self.T_7, self.hot_fluid)  

        T_hs_i = PropsSI('T', 'P', p_3_guess , 'Q', 0, self.fluid) # T_hs_sat_liq
        h_hs_i = PropsSI('H', 'P', p_3_guess , 'Q', 0, self.fluid) # h_hs_sat_liq

        h_cs_i = self.ratio_I * (h_hs_i- h_hs_ex) + h_cs_su
        T_cs_i = PropsSI('T', 'H', h_cs_i, 'P', self.p_7, self.hot_fluid)

        return h_hs_i, h_cs_i, T_hs_i, T_cs_i
    
    def T_equ(self,p_3_guess):
            
        h_hs_i, h_cs_i, T_hs_i, T_cs_i = self.Etat_i(p_3_guess)
    
        return T_cs_i - T_hs_i - self.T_pinch_ex_I
           
    
    def Pressure(self):

        pressure = fsolve(self.T_equ, self.p_3_guess)[0]
       
        return pressure