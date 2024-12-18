

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
import numpy as np

from Nouvelle_resolution_sans_opti import ORC

# region BON 


# Définir le répertoire de sortie (vous pouvez le personnaliser)
output_directory = "Graphes_R245ca_Ts=10"  # Dossier pour sauvegarder les fichiers de sortie

# Charger le fichier Excel
excel_path = os.path.join(os.path.dirname(__file__), '..', 'Param', 'Projet_contraintes.xlsx')
data = pd.read_excel(excel_path)

# Créer le dictionnaire à partir des données Excel
params_excel = {row['Name']: row['Value'] for _, row in data.iterrows()}

params = {}
params.update(params_excel)
inputs = 0

my_ORC = ORC(inputs,params,display=True, output_dir=output_directory)

my_ORC.evaluate()

print('Pe =', my_ORC.Pe)
print('m_tot =', my_ORC.m_tot)
print('eta_cyclex =', my_ORC.eta_cyclex)

#endregion BON


# # region T3
# # Charger le fichier Excel
# excel_path = os.path.join(os.path.dirname(__file__), '..', 'Param', 'Projet_contraintes.xlsx')
# data = pd.read_excel(excel_path)

# # Créer le dictionnaire à partir des données Excel
# params_excel = {row['Name']: row['Value'] for _, row in data.iterrows()}

# params = {}
# params.update(params_excel)
# inputs = 0

# # Initialisation des variables pour T_11
# i_values = np.linspace(5, 30, 10)  # Exemple : linspace de 5°C à 30°C avec 10 points
# T4_surchauffe_values = np.linspace(5, 30, 10)  # Exemple : linspace de 5°C à 30°C avec 10 points
# results = []

# # Boucle sur T_11
# for i in i_values:
#     print(f"Running ORC model for T_11 = {i} °C")
    
#     # Mettre à jour le paramètre dans le dictionnaire
#     params['T_11'] = i  # Mise à jour dynamique
#     params['T_surchauffe_4'] = i  # Mise à jour dynamique
    
#     # Initialiser et évaluer le modèle
#     my_ORC = ORC(inputs, params, display=False)
#     my_ORC.evaluate()
    
#     # Stocker les résultats (vous pouvez choisir ce que vous voulez stocker)
#     results.append({
#         'T_11': i,
#         'eta_cycle': my_ORC.eta_cyclen,  # Rendement énergétique
#         'Pe': my_ORC.Pe,  # Puissance électrique
#         'm_tot': my_ORC.m_tot  # Débit massique total
#     })

# # Convertir les résultats en DataFrame pour analyse
# results_df = pd.DataFrame(results)

# # Afficher les résultats sous forme de tableau
# print(results_df)

# # Graphique des résultats
# plt.figure(figsize=(10, 6))
# plt.plot(results_df['T_11'], results_df['eta_cycle'], label="Rendement énergétique η", marker='o')
# #plt.plot(results_df['T_11'], results_df['Pe'], label="Puissance électrique (kW)", marker='x')
# plt.xlabel("T_11 (°C)")
# plt.ylabel("Valeurs")
# plt.title("Impact de T_11 sur le modèle ORC")
# plt.legend()
# plt.grid()


# # Sauvegarder le graphe final
# output_dir = "graphs"  # Dossier pour sauvegarder le fichier
# os.makedirs(output_dir, exist_ok=True)  # Créer le dossier s'il n'existe pas
# final_graph_filename = os.path.join(output_dir, "Rendement_energétique_T_11.png")
# plt.savefig(final_graph_filename, dpi=300)
# print(f"Graphique final sauvegardé sous : {final_graph_filename}")



# plt.show()

# # # endregion T3


# # region COPS

# # Charger le fichier Excel
# excel_path = os.path.join(os.path.dirname(__file__), '..', 'Param', 'Projet_contraintes.xlsx')
# data = pd.read_excel(excel_path)

# # Créer le dictionnaire à partir des données Excel
# params_excel = {row['Name']: row['Value'] for _, row in data.iterrows()}

# params = {}
# params.update(params_excel)
# inputs = 0

# # Initialisation des variables pour T_11
# variable_1 = np.linspace(274.15, 295.13, 10)  # Attention 0 degré pas bon
# results = []

# # Boucle sur T_11
# for i in variable_1:
#     print(f"Running ORC model for T_11 = {i} °C")
    
#     # Mettre à jour le paramètre dans le dictionnaire
#     params['T_11'] = i  # Mise à jour dynamique
    
#     # Initialiser et évaluer le modèle
#     my_ORC = ORC(inputs, params, display=False)
#     my_ORC.evaluate()
    
#     # Stocker les résultats (vous pouvez choisir ce que vous voulez stocker)
#     results.append({
#         'T_11': i,
#         'eta_cycle': my_ORC.eta_cyclen,  # Rendement énergétique
#         'Pe': my_ORC.Pe,  # Puissance électrique
#         'm_tot': my_ORC.m_tot  # Débit massique total
#     })

# # Convertir les résultats en DataFrame pour analyse
# results_df = pd.DataFrame(results)

# # Afficher les résultats sous forme de tableau
# print(results_df)

# # Graphique des résultats
# plt.figure(figsize=(10, 6))
# plt.plot(results_df['T_11'], results_df['eta_cycle'], label="Rendement énergétique η", marker='o')
# #plt.plot(results_df['T_11'], results_df['Pe'], label="Puissance électrique (kW)", marker='x')
# plt.xlabel("T_11 (°C)")
# plt.ylabel("Valeurs")
# plt.title("Impact de T_11 sur le modèle ORC")
# plt.legend()
# plt.grid()


# # Sauvegarder le graphe final
# output_dir = "Graphs_R245ca_COPS"  # Dossier pour sauvegarder le fichier
# os.makedirs(output_dir, exist_ok=True)  # Créer le dossier s'il n'existe pas
# final_graph_filename = os.path.join(output_dir, "Rendement_energétique_T_11.png")
# plt.savefig(final_graph_filename, dpi=300)
# print(f"Graphique final sauvegardé sous : {final_graph_filename}")



# plt.show()

# # endregion COPS

























#eta_en = my_ORC.eta_toten
# p1 = my_ORC.p1_guess_plot
# p2 = my_ORC.p2_guess_plot
# p5 = my_ORC.p5_guess_plot



# plt.figure(figsize=(10, 6))
# plt.plot(p1, label="p1_guess_plot", linestyle='-', marker='o')
# plt.plot(p2, label="p2_guess_plot", linestyle='--', marker='x')
# plt.plot(p5, label="p5_guess_plot", linestyle='-.', marker='s')

# # Ajout des légendes et titres
# plt.title("Comparison of p1_guess_plot, p2_guess_plot, and p5_guess_plot")
# plt.xlabel("Index")
# plt.ylabel("Values")
# plt.legend()
# plt.grid(True)

# # Affichage du graphe
# plt.show()



















# """
# LELME2150 - Thermal cycles
# Homework 3 - Steam turbine

# Test code for your function

# @author: Antoine Laterre
# @date: October 30, 2022
# """

# #
# #===IMPORT PACKAGES============================================================
# #
# import matplotlib.pyplot as plt
# import pandas as pd
# import os
# import sys

# from Nouvelle_resolution_sans_opti import ORC

# # Charger le fichier Excel
# excel_path = os.path.join(os.path.dirname(__file__), '..', 'Param', 'Projet_contraintes.xlsx')
# data = pd.read_excel(excel_path)

# # Créer le dictionnaire à partir des données Excel
# params_excel = {row['Name']: row['Value'] for _, row in data.iterrows()}

# params = {}
# params.update(params_excel)
# inputs = 0



# my_ORC = ORC(inputs,params,False)
# try:
#     my_ORC.evaluate()
# except ValueError as e:
#     print(f"Une erreur est survenue dans my_ORC.evaluate(): {e}")

