"""
TD4
Gabriel Maisonneuve
Samuel Joly   
git version
"""

# packages calcul numérique 
import scipy as sp
import numpy as np

#package graphique
from matplotlib import pyplot as plt
import math


# Caracteristique de l'etoile
Ts = 5800;          # Temperature surface K
Rs = 6.96E5*1000;   # rayon m
A_etoile = 4*np.pi*Rs**2; # surface etoile 
epsilon = 1;        # Emossivite

# Caracteristique de la planete Hades
r_ha = 6400*1000;    # m
d_ha = 1.5E11;       # m
e_ha = 0.6;
a_ha = 0.39;
A_etoile_hades = 4*np.pi*d_ha**2;
A_disque_hades = np.pi*r_ha**2;

# Caracteristique de la planete Hephaistos
r_he = 8000000;      # m
d_he = 2.2E11;       # m
e_he = 0.82;
a_he = 0.75;
A_etoile_hep = 4*np.pi*d_he**2;
A_disque_hep = np.pi*r_he**2;

# donnees utiles
sigma = 5.67E-8;    # Constante de Stefan‐Boltzmann
h = 6.63E-34;       # Planck
k = 1.38E-23;       # Constante de Boltzmann
c = 3E8;            # Vitesse lumiere



# ********** Question 1 **********
# ***** 1a) *****
flux_surf_etoile = epsilon*sigma*Ts**4;
flux_therm_etoile = flux_surf_etoile*A_etoile;

# ***** 1b) *****
# Planete Hades
flux_surf_hades = flux_therm_etoile/A_etoile_hades;
flux_therm_hades = flux_surf_hades*A_disque_hades;

# Planete Hephaistos
flux_surf_hep = flux_therm_etoile/A_etoile_hep;
flux_therm_hep = flux_surf_hep*A_disque_hep;

# ***** 1c) *****
lambda_max = h*c/(5*k*Ts);


print("Tianlang")
print("Flux surfacique emis = ", flux_surf_etoile, "W/m^2")
print("Flux thermique emis = ", flux_therm_etoile, "J/s")
print("Hades")
print("Flux surfacique emis = ", flux_surf_hades, "W/m^2")
print("Flux thermique emis = ", flux_therm_hades, "J/s")
print("Hephaistos")
print("Flux surfacique emis = ", flux_surf_hep, "W/m^2")
print("Flux thermique emis = ", flux_therm_hep, "J/s")
print("Longueur d'onde max = ", lambda_max, "m")

# ********** Question 2 **********
# ***** 2a) *****
flux_surf_hades_emis = flux_surf_hades/4;
flux_surf_hep_emis = flux_surf_hep/4;

# ***** 2b) *****
T_hades = (flux_surf_hades_emis/(epsilon*sigma))**(0.25);
T_hep = (flux_surf_hep_emis/(epsilon*sigma))**(0.25);

# ***** 2c) *****
# avec atmosphere
# hades
Js_ha = flux_surf_hades;
Js_capte = (1-a_ha)*Js_ha; 
Jrs = Js_capte/4;
Jrs = Js_capte/(1-1.5*e_ha);
Ts_surface = (Jrs/(e_ha*sigma*epsilon))**(0.25);


# avec atmosphere
# hep
Js_he = flux_surf_hep;
Js_capte_he = (1-a_he)*Js_he;
Jrs_he = Js_capte_he/4;




