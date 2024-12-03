import numpy as np
import matplotlib.pyplot as plt
from CoolProp.CoolProp import PropsSI

# Définir les limites de pression pour tracer la courbe de saturation
pressures = np.linspace(PropsSI('P_MIN', 'Pentane'), PropsSI('P_CRITICAL', 'Pentane'), 500)
print("P_min ===",PropsSI('P_MIN', 'Pentane')/1e5,'bar')
print("P_critical ===",PropsSI('P_CRITICAL', 'Pentane')/1e5,'bar')


# Listes pour stocker la température et l'entropie pour les phases liquide et vapeur
T_liq = []
T_vap = []
S_liq = []
S_vap = []

# Boucle sur chaque pression pour obtenir la température de saturation, l'entropie du liquide et de la vapeur
for P in pressures:
    T = PropsSI('T', 'P', P, 'Q', 0, 'Pentane')  # Température de saturation (Q=0 pour liquide)
    S_f = PropsSI('S', 'P', P, 'Q', 0, 'Pentane')  # Entropie du liquide saturé
    S_g = PropsSI('S', 'P', P, 'Q', 1, 'Pentane')  # Entropie de la vapeur saturée

    T_liq.append(T - 273.15)  # Convertir en Celsius
    T_vap.append(T - 273.15)  # Convertir en Celsius
    S_liq.append(S_f / 1000)  # Convertir en kJ/kg.K
    S_vap.append(S_g / 1000)  # Convertir en kJ/kg.K

# Tracer le diagramme T-s
plt.figure(figsize=(10, 6))
plt.plot(S_liq, T_liq, label='Ligne de saturation - Liquide', color='blue')
plt.plot(S_vap, T_vap, label='Ligne de saturation - Vapeur', color='red')

# Tracer les courbes isobares pour différentes pressions
isobare_pressures = np.linspace(pressures[0], pressures[-1], 8)  # en Pascal
for P in isobare_pressures:
    S_vals = []
    T_vals = []
    for Q in np.linspace(0, 1, 100):  # Parcours de la qualité de 0 à 1
        T = PropsSI('T', 'P', P, 'Q', Q, 'Pentane')
        S = PropsSI('S', 'P', P, 'Q', Q, 'Pentane')
        T_vals.append(T - 273.15)  # Convertir en Celsius
        S_vals.append(S / 1000)  # Convertir en kJ/kg.K
    plt.plot(S_vals, T_vals, '--', label=f'Isobare P={P/1e5:.1f} bar')

# Configurations du graphique
plt.xlabel('Entropie (kJ/kg.K)')
plt.ylabel('Température (°C)')
plt.title('Diagramme T-s du Pentane avec Isobares')
plt.legend()
plt.grid(True)
plt.show()


