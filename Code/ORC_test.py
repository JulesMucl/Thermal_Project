
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

from ORC_group_21 import ORC

# Charger le fichier Excel
excel_path = os.path.join(os.path.dirname(__file__), '..', 'Param', 'Projet_contraintes.xlsx')
data = pd.read_excel(excel_path)

# Créer le dictionnaire à partir des données Excel
params_excel = {row['Name']: row['Value'] for _, row in data.iterrows()}

params = {}
params.update(params_excel)
inputs = 0



my_ORC = ORC(inputs,params,True)
try:
    my_ORC.evaluate()
except ValueError as e:
    print(f"Une erreur est survenue dans my_ORC.evaluate(): {e}")

#eta_en = my_ORC.eta_toten
p1 = my_ORC.p1_guess_plot
p2 = my_ORC.p2_guess_plot
p5 = my_ORC.p5_guess_plot



plt.figure(figsize=(10, 6))
plt.plot(p1, label="p1_guess_plot", linestyle='-', marker='o')
plt.plot(p2, label="p2_guess_plot", linestyle='--', marker='x')
plt.plot(p5, label="p5_guess_plot", linestyle='-.', marker='s')

# Ajout des légendes et titres
plt.title("Comparison of p1_guess_plot, p2_guess_plot, and p5_guess_plot")
plt.xlabel("Index")
plt.ylabel("Values")
plt.legend()
plt.grid(True)

# Affichage du graphe
plt.show()
