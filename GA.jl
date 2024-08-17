push!(LOAD_PATH, "$(pwd())/SCDP")
using FileIO 
using Tables
using CSV
using StaticArrays
using JLD2
using SatelliteToolbox
using SCDP
include("SCDP/SCDP/src/constantes_physiques.jl")
include("Parameters.jl")

#Hyper parameter
instantInitial = SatelliteToolbox.DateTime(1970, 1, 1, 0) 
m_simu = 1 # constante dans la relation de périodicité
n_0 = 4 #nombre d'intervals

nbCibles = 1
# L'instance du problème que l'on considère . 
# inst = load("instances/instance_"*string(nbCibles)*"_cibles.jld2", "inst")
# inst = instanceAleatoire(nbCibles, alphaHalf, nbCibles*n_0, instantInitial, m_simu, n_0)
# save("instances/instance_"*string(nbCibles)*"_cibles.jld2", "inst", inst)

Latitude = [22]*pi/180
Longitude = [33]*pi/180 
coordsCibles = [SA[Latitude[i], Longitude[i]] for i in 1:nbCibles]
inst = DonneesSCDP(coordsCibles, alphaHalf, nbCibles*n_0, instantInitial, m_simu, n_0)
save("instances/instance_"*string(nbCibles)*"_cibles.jld2", "inst", inst)

ecrireInst(inst)

na = [CartesianIndex()]
nbMaxGen = 500 # MODIFFF
ϵ = 0.005

# Nombre de répétitions pour chaque valeur des paramètres. 
nbRep = 1  # MODIFFF
dureesCalcul = Vector{Float64}(undef, nbRep)

probasMut = [0.85]
probasModifNbSats = [0.01]

stats = Vector{Float64}[]
sizehint!(stats, nbMaxGen)
function fDiag(nbIndiv; ite, pop, donneesCrit)
   fCVMoy = sum(getfield.(pop[1:nbIndiv], :CV))/nbIndiv
   pourcentIndiv = 100sum(getfield.(pop[1:nbIndiv], :CV) .≈ 0)/nbIndiv
   nbSatMoy = sum( getindex.(getfield.(pop[1:nbIndiv], :y), 1) )/nbIndiv
   nbSatMin = pop[1].y[1]
   push!(stats, SA[fCVMoy, pourcentIndiv, nbSatMoy, nbSatMin])
end
nomCols = ["CV moyenne", "Pourcent indiv respectant contrainte", "Nb sat moyen", 
           "Nb sat min"]



# GA
nbIndiv = 1000 # MODIFFFF
for pMut in probasMut
   for pModif in probasModifNbSats 
	  TitleDir = "$nbIndiv-$pMut-$pModif"
	  if !isdir(TitleDir) 
		 mkdir(TitleDir)*"/"
      end
	  for ind_rep in 1:nbRep
         empty!(stats)
         dureesCalcul[ind_rep] = @elapsed local pop = 
            resoudre_SCDP(inst, nbIndiv, nbMaxGen, ϵ; probaMut = pMut, 
                            probaModifNbSat = pModif, 
                            fDiag = (;kwargs...) -> fDiag(nbIndiv; kwargs...)) 
         
         matStat = vcat([l[na, :] for l in stats]...)

         # Enregistrement de la population final et des statistiques 
         save(TitleDir *"/"*string(nbCibles)*"targets/pop$(ind_rep).jld2", "pop", pop)
         CSV.write(TitleDir *"/"*string(nbCibles)*"targets/stats$(ind_rep).csv", Tables.table(matStat, header = nomCols))
      end
      # CSV.write(TitleDir *"/"*string(nbCibles)*"targets/tps_calcul.csv", Tables.table(dureesCalcul), 
                # header = ["Durées de calcul (s)"])
   end
end
if !isdir(string(nbIndiv)*"-"*string(probasMut[1])*"-"*string(probasModifNbSats[1])*"\\"*string(nbCibles)*"targets")
	mkdir(string(nbIndiv)*"-"*string(probasMut[1])*"-"*string(probasModifNbSats[1])*"\\"*string(nbCibles)*"targets")
end



# Affichage résultats 
@load string(nbIndiv)*"-"*string(probasMut[1])*"-"*string(probasModifNbSats[1])*"\\"*string(nbCibles)*"targets\\pop1.jld2" pop
# on selectionne meilleur élement, conste est ordonné de façon décroissante 
indiv = pop[1]
conste = consteVar2Conste(indiv.x.p[1:4indiv.x.nbBlocs[]], inst)
latTarget = fill(0.0,nbCibles)
LongTarget = fill(0.0,nbCibles)
for i in 1:nbCibles
	latTarget[i] = inst.coordsCibles[i][1]*(180/pi)
	LongTarget[i] = inst.coordsCibles[i][2]*(180/pi)
end
println(latTarget)
println(LongTarget)
println("NUM_SATELLITE = "*string(conste.nbSats))
println("Altitude = "*string(conste.as))
println("Inc = "*string(conste.is))
println("RAAN = "*string(conste.Ωs))
println("MeanAnomaly = "*string(conste.Ms))

# angleOuverture = inst.angleOuverture
# println(angleOuverture*(180/pi))






