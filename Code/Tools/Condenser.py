from CoolProp.CoolProp import PropsSI
from scipy.optimize import fsolve

class Condenser(object):

    def __init__(self, orc_instance):

        self.orc            = orc_instance

        self.cold_fluid     = self.orc.cold_fluid
        self.dot_m_cf       = self.orc.m_dot_cf

        self.fluid          = self.orc.fluid
        self.m_dot_tot   = self.orc.m_dot_tot

        self.T_10           = self.orc.T_10
        self.T_11           = self.orc.T_11
        self.p_10           = self.orc.p_10
        self.h_10           = self.orc.h_10
        self.h_11           = self.orc.h_11

        self.p_6_guess      = self.orc.p_6_guess
        self.T_pinch_cd  = self.orc.T_pinch_cd
        self.T_cd_subcool = self.orc.T_cd_subcool


    def Etat_6_guess(self,p_6_guess):

        T6 = PropsSI('T','P',p_6_guess,'Q',0,self.fluid) - self.T_cd_subcool
        h6 = PropsSI('H','P',p_6_guess,'T',T6,self.fluid)
        return h6
    
    def Etat_i(self,p_6_guess):

        h6 = self.Etat_6_guess(p_6_guess)

        h6i = PropsSI('H','P',p_6_guess,'Q',1,self.fluid)
        T6i = PropsSI('T','P',p_6_guess,'Q',1,self.fluid)

        hcsi = (self.m_dot_tot/self.dot_m_cf) * (h6i - h6) + self.h_11
        Tcsi = PropsSI('T','H',hcsi,'P',self.p_10,self.cold_fluid)

        return T6i, Tcsi
    
    def Iter_function(self,p_6_guess):

        T6i, Tcsi = self.Etat_i(p_6_guess)

        return T6i - Tcsi - self.T_pinch_cd
    
    def P6(self):

        p_6 = fsolve(self.Iter_function,self.p_6_guess)

        return p_6[0]
    

    





    


        



