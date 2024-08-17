# *** Constantes physiques ***

# Période de révolution de la terre (SI)
T_terre = 24*3600

# Rayon terrestre (SI)
R_terre = 6371e3

# Paramètre gravitationnel standard terrestre (SI)
μ = 3.986_004_418e14 

# Altitudes minimales et maximales (SI) pour un satellite sur une orbite L.E.O.
h_min = 400e3
h_max = 1400e3 

# Rayons d'orbites minimales et maximales (SI) pour un satellite sur une orbite
# L.E.O. 
a_min = R_terre + h_min
a_max = R_terre + h_max
