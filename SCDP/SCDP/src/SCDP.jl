"""
Résoudre le SCDP (Satellites Constellation Design Problem) avec un algorithme génétique. 

Ce module implémente plusieurs variations de l'algorithme NSGAII afin de résoudre le SCDP.
Il contient aussi des fonctions de visualisation. 
"""
module SCDP 

export resoudre_SCDP, Conste, GrpSat, GrpType, instanceAleatoire, ecrireInst,
       ecrireIndiv, ecrirePop, DonneesSCDP,
       ecrireIntervallesSurvol, dessinerIntervallesSurvol, dessinerConsteECI!,
       dessinerSolECEF!, extraireDeConsteFixe, extraireDeConsteVar, 
       consteFixe2Conste, consteVar2Conste

using NSGAII 
using LinearAlgebra
using StaticArrays
using Dates
using SatelliteToolbox
using Random 
using Plots 
using PrettyTables
using StatsBase

import ColorSchemes.darkrainbow

# Backend pour les plots 
pyplot()
# plotlyjs()

# newaxis
na = [CartesianIndex()]

include("constantes_physiques.jl")
include("structs.jl")
include("simu.jl")
include("visu.jl")

# import Metaheuristics.PerformanceIndicators.hypervolume 


""" 
Utiliser l'algorithme génétique NSGAII pour résoudre l'instance `inst` du problème. Renvoie
la population obtenue à la fin de l'exécution de l'algorithme. 

# Arguments 
- `inst` : une instance du problème. 
- `nbIndiv` : nombre d'individus dans la population. 
- `nbGen` : nombre de générations. 
- `periodeStat`
"""

function resoudre_SCDP(inst::DonneesSCDP, nbIndiv::Integer, nbMaxGen::Integer, 
      Grp::Union{Type{GrpSat}, Type{GrpType}}, ϵ::Union{Nothing, Real}; 
      probaMut::Real = 0.05, nbGenGlissant::Int = 60, 
      periodeDiag::Int = 1, fDiag::Function = (; vargs...) -> nothing)

   # Codage
   # Types des éléments orbitaux pour un satellite. 
   typesElemsOrbitaux = [:Bin, :Int, :Cont, :Cont, :Cont]
   # Bornes des éléments orbitaux (Activation,nombre revolution, inclinaison,longtiude,anomalie)
   bornesInf = [0, inst.n_min, 0, -π, 0]
   bornesSup = [1, inst.n_max, π, π, 2π]

   if Grp == GrpSat 
      bc = BinaryCoding(inst.nbCS, repeat(typesElemsOrbitaux, outer = inst.nbMaxSat),
                           repeat(bornesInf, outer = inst.nbMaxSat),
                           repeat(bornesSup, outer = inst.nbMaxSat))
   elseif Grp == GrpType 
      bc = BinaryCoding(inst.nbCS, repeat(typesElemsOrbitaux, inner = inst.nbMaxSat),
                           repeat(bornesInf, inner = inst.nbMaxSat),
                           repeat(bornesSup, inner = inst.nbMaxSat))
   end
   
   # Objectif 
   # z doit renvoyer un tuple
   # (on minimse lorsque l'on appelle nsga)
   z(conste) = (calculNbSatsActifsConsteFixe(Grp, conste, inst), )

   # Fonction de violation de contrainte. Les variables entières sont codés en binaire d'une
   # telle manière qu'elles peuvent exécéder leur bornes supérieures après croisement et
   # mutation. On exclut ces cas de la population grâce à la fonction de violation de
   # contrainte. 
   CV(conste) = begin 
       _, ns, is, Ωs, Ms = extraireDeConsteFixe(Grp, conste, inst)

       nb_ns_depasse = sum(ns .> inst.n_max) 
       nb_ns_depasse > 0 && return inst.nbCibles*inst.n_0 + nb_ns_depasse

       as = Float64[calculRayonOrbite(n, inst.m_simu) for n in ns]
       return calculNbIntervsSansSurvol(inst, as, is, Ωs, Ms)
    end
   # CV(conste) = begin 
   #    _, ns, is, Ωs, Ms = extraireDeConsteFixe(Grp, conste, inst)
   #    any(ns .> inst.n_max) && return Inf 
   #    as = Float64[calculRayonOrbite(n, inst.m_simu) for n in ns]
   #    return calculNbIntervsSansSurvol(inst, as, is, Ωs, Ms)
   # end
   
   # Critère d'arrêt
   ca = ϵ != nothing ? (donneesCrit::Dict, pop) -> 
      critArretECGlissant(donneesCrit::Dict, pop, ϵ, nbGenGlissant, inst.nbMaxSat) : 
         (donneesCrit::Dict, pop) -> false

   popFinal = nsga(nbIndiv, nbMaxGen, z, bc; fCV = CV, pmut = probaMut, critArret = ca,
                   fDiag = fDiag, periodeDiag = periodeDiag)

   return popFinal
end

function resoudre_SCDP(inst::DonneesSCDP, nbIndiv::Integer, nbMaxGen::Integer,
      ϵ::Union{Nothing, Real}; probaModifNbSat::Real = 0.05, probaMut::Real = 0.05, 
      nbGenGlissant::Int = 60, fDiag::Function = (; vargs...) -> nothing, periodeDiag::Int = 1)

   # Codage pour des phénotypes et génotypes de longueur variable
   typesElemsOrbitaux = [:Int, :Cont, :Cont, :Cont]
   # Bornes des éléments orbitaux
   # bornesInf = [inst.n_min, 0, -π, 0]
   # bornesSup = [inst.n_max, π, π, 2π]
   bornesInf = [inst.n_min, 0, 0, 0]
   bornesSup = [inst.n_max, π, 2π, 2π]
   bc = BinaryCoding(inst.nbCS, typesElemsOrbitaux, bornesInf, bornesSup)
   bvc = BinaryVariableCoding(bc, inst.nbMaxSat)

   # Objectif 
   # z doit renvoyer un tuple
   # (on minimse lorsque l'on appelle nsga)
   z(x) = (x.nbBlocs[],)

   # Fonction de violation de contrainte. Les variables entières sont codés en binaire d'une
   # telle manière qu'elles peuvent exécéder leur bornes supérieures après croisement et
   # mutation. On exclut ces cas de la population grâce à la fonction de violation de
   # contrainte. 
   CV(x) = begin 
      _, ns, is, Ωs, Ms = extraireDeConsteVar(view( x.p, 1:(bvc.bc.nbvar*x.nbBlocs[]) ))

      nb_ns_depasse = sum(ns .> inst.n_max) 
      nb_ns_depasse > 0 && return inst.nbCibles*inst.n_0 + nb_ns_depasse

      as = Float64[calculRayonOrbite(n, inst.m_simu) for n in ns]
      return calculNbIntervsSansSurvol(inst, as, is, Ωs, Ms)
   end

   # Critère d'arrêt
   ca = ϵ != nothing ? (donneesCrit::Dict, pop) -> 
      critArretECGlissant(donneesCrit::Dict, pop, ϵ, nbGenGlissant, inst.nbMaxSat) : 
         (donneesCrit::Dict, pop) -> false

   popFinal = nsga(nbIndiv, nbMaxGen, z, bvc; fCV = CV, pmut = probaMut, 
                   probaModifNbBlocs = probaModifNbSat, critArret = ca,
                   fDiag = fDiag, periodeDiag = periodeDiag)

   return popFinal
end

end

# *** Test de resoudre_NSGAII ***
# Vérification avec SatelliteToolbox
if false

   # conste = [true, 13, 1.57080, 1.88495, 0.94248, 
   #          false, 0, 0, 0, 0, 
   #          true, 16, 0.3, 0.8, 0.1,  
   #          true, 15, 0.4, 0.3, 0.13, 
   #          false, 0, 0, 0, 0, 
   #          false, 0, 0, 0, 0, 
   #          false, 0, 0, 0, 0, 
   #          true, 0.466, 0.6, 0.144, 0.78]

   coordsCibles = [SA[0, 0], SA[0, π/2], SA[0, π], SA[0, -π/2]]
   angleOuverture = π/6 
   instantInitial = DateTime(2022, 3, 20, 12, 0)
   m_simu = 1 
   n_0 = 6 
   Δ = 30 

   inst = DonneesProb(coordsCibles, angleOuverture, instantInitial, m_simu, n_0, Δ)
   ecrireInst(inst)

   nbGen = 1
   nbIndiv = 10

   pop, _ = resoudre_NSGAII(inst, nbIndiv, nbGen)

   # On trie la population finale selon le rang
   # sort!(pop, by = ind -> ind.rank)

   ecrirePop(inst, pop)
   conste = pop[1].x.p
   ecrireIntervallesNonSurvol(inst, conste)

   # On récupère les paramètres orbitaux de chaque satellite actif de la constellation. 
   nbMaxSat = inst.nbMaxSat
   actifs = Bool.(conste[1:nbMaxSat])
   nbSat = sum(actifs)
   ns = conste[(nbMaxSat+1):(2nbMaxSat)][actifs]
   as = @. (μ*(m_simu*T_terre/(2π*ns))^2)^(1/3)
   is = conste[(2nbMaxSat+1):(3nbMaxSat)][actifs]
   Ωs = conste[(3nbMaxSat+1):(4nbMaxSat)][actifs]
   Ms = conste[(4nbMaxSat+1):(5nbMaxSat)][actifs]

   # Propagateurs des satellites actifs 
   props = init_orbit_propagator.(Val(:twobody), 
      KeplerianElements.(inst.instantInitialJulian, as, 0, is, Ωs, 0, Ms))

   angleElevation = @. π/2 - asin(as[1]*sin(inst.angleOuverture)/R_terre)

   intervallesNonSurvol = ground_station_gaps(props[1], inst.coordsCiblesECEF[1], T_terre, TOD(), 
                                              PEF(), θ = angleElevation)

   intervallesNonSurvolSec = @. 1e-3*getfield(intervallesNonSurvol - instantInitial, :value)
end


