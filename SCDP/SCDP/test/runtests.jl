"""
Tester le module SCDP.
"""

using SDCP
using Test
using Random 
using Dates
using FileIO
using SatelliteToolbox
using Plots
plotlyjs()
Random.seed!(0)

# Instance de référence pour les tests
inst = load("../../../instances/instance_3_cibles.jld2", "inst")
ecrireInst(inst)

# Test de la discrétisation. Ce test est visuel : les zones de couvertures aux instants
# successifs doivent se superposer d'environ 20%.
println(@test begin 
   prop = init_orbit_propagator(Val(:twobody), 
                             KeplerianElements(inst.instantInitialJulian, 
                                          DCP.calculRayonOrbite(inst.n_max, inst.m_simu), 0, 
                                          0, 0, 0, 0))
   nbInstants = 10

   title = "Trajectoire d'un satellite<br>et zones couvertes aux $nbInstants \
   premiers instants considérés"
   plt = plot(title = title, size = (1000, 1000),
              xlabel = "x", ylabel = "y", zlabel = "z")
   DCP.dessinerTerre!(plt)
   DCP.dessinerTrajECEF!(plt, prop, inst.angleOuverture, inst.ts[1:nbIndiv], 
                         inst.mats_ECI_ECEF[1:nbInstants], Colors.parse(Colorant, :blue), 
                         "Trajectoire", 200)
   display(plt)
   println("Test ok ? [o/n]")
   resTest = readline()
   resTest == "o" ? true : false 
end
)

# Test du calcul du nombre d'intervalles où une cible n'est pas survolé 
# (fonction calculNbIntervsSansSurvol) . Les vrais résultats
# sont connus grace à la fonction dessinerIntervallesSurvol qui elle-même repose sur la
# fonction calculIntervallesSurvol. Les fonctions calculNbIntervsSansSurvol et
# calculIntervallesSurvol sont indépendantes mais ont toute les deux été écrites par la même 
# personne. 

# println(@test begin 
#   ns1 = [14, 

# On veut vérifier que l'objectif décroit au cours de l'évolution de NSGAII. Cela est vrai à
# condition d'assigner une valeur d'objectif infinie aux individus qui ne respectent pas la
# contrainte. 
calculStat(pop, nbIndiv, stats) = push!(stats, pop[1].CV ≈ 0 ? pop[1].y[1] : Inf)

nbIndiv = 10
nbGen = 10

println(@test begin 
   stats = Float64[]
   resoudre_DCP(inst, nbIndiv, nbGen, GrpSat; periodeStat = 1, 
                calculStat = (pop, nbIndiv) -> calculStat(pop, nbIndiv, stats))
   all(stats[2:end] .≤ stats[1:end-1])
end
)

println(@test begin 
   stats = Float64[]
   resoudre_DCP(inst, nbIndiv, nbGen, GrpType; periodeStat = 1, 
                calculStat = (pop, nbIndiv) -> calculStat(pop, nbIndiv, stats))
   all(stats[2:end] .≤ stats[1:end-1])
end
)

println(@test begin 
   stats = Float64[]
   resoudre_DCP(inst, nbIndiv, nbGen; periodeStat = 1, 
                calculStat = (pop, nbIndiv) -> calculStat(pop, nbIndiv, stats))
   all(stats[2:end] .≤ stats[1:end-1])
end
)
