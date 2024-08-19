# *** Constantes physiques ***

# Période de révolution de la terre (SI)
T_terre = 24*3600

# Rayon terrestre (SI)
R_terre = 6378136.3

# Paramètre gravitationnel standard terrestre (SI)
μ = 3.986_004_418e14 

# Altitudes minimales et maximales (SI) pour un satellite sur une orbite L.E.O.
h_min = 400e3
h_max = 1400e3 

# Rayons d'orbites minimales et maximales (SI) pour un satellite sur une orbite
# L.E.O. 
a_min = R_terre + h_min
a_max = R_terre + h_max

# Circumeference of a sphere (Earth)
CircumferenceEarth = R_terre*2*pi # Km

# Pas de temps de la simulation
dt = 180.0# sec

# Number of interval periodic in a Earth rotation
n_0 = 12

# Nombre de pas de temps dans la simulation
NUM_TIME = convert(Int,round(T_terre/dt))

# Number of period 
NUM_PERIOD = convert(Int,T_terre/n_0/dt)

# Earth's rotation rate in rad/s
we= 7.2921e-5

# Parameters for Satellite Toolbox
Modele = Val(:twobody)

# Obtention ensemble de valeurs admis pour hauteur satellite en fonction du nombre de N revolutions du satellite
N_test=1:24
altitudeSet = Float64[]
PERIOD_SAT = Float64[]
for i in N_test
    altitudeSetVal = ( ( μ*((T_terre)/N_test[i])^2 ) / 4π^2 )^(1/3) - R_terre
    if altitudeSetVal >= 400000 && altitudeSetVal <= 1400000
        append!(altitudeSet, [altitudeSetVal])
        append!(PERIOD_SAT, [T_terre/N_test[i]])
    end
end
#println(altitudeSet)

PERIOD_SAT_MIN = PERIOD_SAT[end]
altitudeLength = length(altitudeSet)
time_zero_simulation = date_to_jd(2000,1,1,12,0,0)#1970
	
# Compute Theta and alpha of a satellite
theta_min = ((2*pi*dt)/(2*PERIOD_SAT_MIN))*1.2
numbOfDiscretize = Int(ceil(pi/theta_min))
alphaHalf = atan((sin(theta_min))/(((R_terre+altitudeSet[end])/R_terre)-cos(theta_min)))

Theta = fill(0.0,length(altitudeSet))
for j in eachindex(altitudeSet)
    if ((R_terre+altitudeSet[j])/R_terre)*sin(alphaHalf) > 1
        alpha_lim = asin(R_terre / (R_terre + altitudeSet[j]))
        Theta[j] = (-alpha_lim + asin(round(((R_terre+altitudeSet[j])/R_terre)*sin(alpha_lim),digits=6)))
    else	
        Theta[j] = (-alphaHalf + asin(((R_terre+altitudeSet[j])/R_terre)*sin(alphaHalf)))
    end
end

#Number total sat, discretization od Keplerian's parameters and the related cardinality
global InclinaisonSet = LinRange(0.0,pi,numbOfDiscretize)
global RAANSet = LinRange(0.0,2*pi,numbOfDiscretize)
global MeanAnomalySet = LinRange(0.0,2*pi,numbOfDiscretize)
global N_inclinaison = size(InclinaisonSet)[1]
global N_RAAN = size(RAANSet)[1]
global N_meanAnomaly = size(MeanAnomalySet)[1]
global NbrTOTALsat = N_inclinaison*N_RAAN*N_meanAnomaly*altitudeLength
