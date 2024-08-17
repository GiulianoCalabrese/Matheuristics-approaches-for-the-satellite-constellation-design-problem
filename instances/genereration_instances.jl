# Générations d'instances 
using Random 
using FileIO
using Dates
using StaticArrays
using SatelliteToolbox
push!(LOAD_PATH, "$(pwd())/../SCDP")
include("../SCDP/SCDP/src/constantes_physiques.jl")

using SCDP

Random.seed!(0)


dt = 180
N_test=1:24
global altitudeSet = Float64[]
global PERIOD_SAT = Float64[]
for i in N_test
	altitudeSetVal = ( ( μ*((T_terre)/N_test[i])^2 ) / 4π^2 )^(1/3) - R_terre
	if altitudeSetVal >= 400000 && altitudeSetVal <= 1400000
		append!(altitudeSet, [altitudeSetVal])
		append!(PERIOD_SAT, [T_terre/N_test[i]])
	end
end
# println(altitudeSet)

global PERIOD_SAT_MIN = PERIOD_SAT[end]
global theta_min = ((2*pi*dt)/(2*PERIOD_SAT_MIN))*1.2
global numbOfDiscretize = Int(ceil(pi/theta_min))
global alphaHalf = atan((sin(theta_min))/(((R_terre+altitudeSet[end])/R_terre)-cos(theta_min)))

Theta = fill(0.0,length(altitudeSet))
for j in 1:length(altitudeSet)
	if ((R_terre+altitudeSet[j])/R_terre)*sin(alphaHalf) > 1
		alpha_lim = asin(R_terre / (R_terre + altitudeSet[j]))
		global Theta[j] = (-alpha_lim + asin(round(((R_terre+altitudeSet[j])/R_terre)*sin(alpha_lim),digits=6)))
	else	
		global Theta[j] = (-alphaHalf + asin(((R_terre+altitudeSet[j])/R_terre)*sin(alphaHalf)))
	end
end

println(PERIOD_SAT)
println(alphaHalf*180/pi)
println(Theta*180/pi)

# angleOuverture = 20π/180
instantInitial = DateTime(1970, 1, 1, 0) 
m_simu = 1 # constante dans la relation de périodicité
n_0 = 3 #nombre d'intervals

nbCibles = 10
# inst = instanceAleatoire(nbCibles, angleOuverture, instantInitial, m_simu, n_0)
# save("instance_10_cibles.jld2", "inst", inst)
# ecrireInst(inst)

nbCibles = 5
# inst = instanceAleatoire(nbCibles, angleOuverture, instantInitial, m_simu, n_0)
# save("instance_5_cibles.jld2", "inst", inst)
# ecrireInst(inst)

nbCibles = 40
# inst = instanceAleatoire(nbCibles, angleOuverture, instantInitial, m_simu, n_0)
# save("instance_40_cibles.jld2", "inst", inst)
# ecrireInst(inst)

nbCibles = 3
# inst = instanceAleatoire(nbCibles, angleOuverture, instantInitial, m_simu, n_0)
# save("instance_3_cibles.jld2", "inst", inst)
# ecrireInst(inst)

Latitude = [14,-37,-57]*pi/180
Longitude = [206, 18, 56]*pi/180 
coordsCibles = [SA[Latitude[i], Longitude[i]] for i in 1:nbCibles]
inst = DonneesSCDP(coordsCibles, alphaHalf, instantInitial, m_simu, n_0)
println(coordsCibles)

println(inst.dt)
save("instance_3_OPTE_cibles.jld2", "inst", inst)
# ecrireInst(inst)

# println(inst)
