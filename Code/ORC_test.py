
"""
LELME2150 - Thermal cycles
Homework 3 - Steam turbine

Test code for your function

@author: Antoine Laterre
@date: October 30, 2022
"""

#
#===IMPORT PACKAGES============================================================
#
import matplotlib.pyplot as plt
import pandas as pd
import os
import sys

from Nouvelle_resolution_sans_opti import ORC

# Charger le fichier Excel
excel_path = os.path.join(os.path.dirname(__file__), '..', 'Param', 'Projet_contraintes.xlsx')
data = pd.read_excel(excel_path)

# Créer le dictionnaire à partir des données Excel
params_excel = {row['Name']: row['Value'] for _, row in data.iterrows()}

params = {}
params.update(params_excel)
inputs = 0



my_ORC = ORC(inputs,params,False)
try:
    my_ORC.evaluate()
except ValueError as e:
    print(f"Une erreur est survenue dans my_ORC.evaluate(): {e}")

