"""
TD4
Gabriel Maisonneuve 2133674
Samuel Joly 2116135
"""

# packages calcul numérique 
import scipy as sp
import numpy as np

#package graphique
from matplotlib import pyplot as plt
import math
from sympy import symbols, Eq, solve

# Caracteristique de l'etoile
Ts = 5800;          # Temperature surface K
Rs = 6.96E5*1000;   # rayon m
A_etoile = 4*np.pi*Rs**2; # surface etoile 
epsilon = 1;        # Emissivite

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
# flux sans atmosphere
flux_surf_hades_emis = flux_surf_hades/4;
flux_surf_hep_emis = flux_surf_hep/4;

# ***** 2b) *****
# temperature sans atmosphere
T_hades = (flux_surf_hades_emis/(epsilon*sigma))**(0.25); # valider
T_hep = (flux_surf_hep_emis/(epsilon*sigma))**(0.25); # valider


# ***** 2c) *****
# avec atmosphere
# hades
Js_ha = flux_surf_hades_emis;
Js_capte = (1-a_ha)*Js_ha; 
Jrs = Js_capte/(1-e_ha/2);
Ts_surface_hades = (Jrs/(sigma*epsilon))**(0.25);

# temperature de l'atmosphere de hades
Jra_ha = Jrs*e_ha/2;
Ta_hades = (Jra_ha/(e_ha*sigma))**(0.25);

# avec atmosphere
# hep
Js_he = flux_surf_hep_emis;
Js_capte_he = (1-a_he)*Js_he;
Jrs_he = Js_capte_he/(1-e_he/2);
Ts_surface_hep = (Jrs_he/(sigma*epsilon))**(0.25);

# temperature de l'atmosphere de hep
Jra_he = Jrs_he*e_he/2;
Ta_hep = (Jra_he/(e_he*sigma))**(0.25);

# ***** 2d) *****
Ta_hep_c = (((((1-a_he)*(((flux_surf_etoile*A_etoile)/(4*np.pi*d_he**2))/4))/(1-e_he/2))*e_he/2)/(e_he*sigma))**(0.25);

# Define the variable
new_d_hep = symbols('new_d_hep')

# Define the equation
equation = Eq(Ts_surface_hades, (((((1-a_he)*(((flux_surf_etoile*A_etoile)/(4*np.pi*new_d_hep**2))/4))/(1-e_he/2))*e_he/2)/(e_he*sigma))**(0.25))

# Solve the equation for x
solution = solve(equation, new_d_hep)

# ***** reponses *****
print("Flux surfacique émis par la surface d'Hadès = ", flux_surf_hades_emis, "W/m^2")
print("Flux surfacique émis par la surface d'Héphaïstos = ", flux_surf_hep_emis, "J/s")

print("Température de surface d'Hadès sans atmosphère = ", T_hades, "K")
print("Température de surface d'Héphaïstos sans atmosphère = ", T_hep, "K")

print("Hades")
print("Température de surface avec atmosphère = ", Ts_surface_hades, "K")
print("Température de l'atmosphère = ", Ta_hades, "K")

print("Héphaïstos")
print("Température de surface avec atmosphère = ", Ts_surface_hep,"K")
print("Température de l'atmosphère = ", Ta_hep, "K")

print("Le nouveau paramètre est la distance entre l'étoile et la planète")
print("La nouvelle distance est", solution[1], "m")
