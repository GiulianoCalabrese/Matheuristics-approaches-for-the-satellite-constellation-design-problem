"""
Utiliser l'algorithme NSGAII pour résoudre une problème d'optimisation multiobjectif. 

# Arguments
- _ : le type concret d'un individu sur lequel travaille `_nsga`. C'est un sous-type de 
   `Indiv{G, Y} where {G, Y}`.
- `sense` : un des types constants `Max` ou `Min`. Détermine si on maximise ou minimise 
   l'objectif. 
- `popSize` : le nombre d'individus qui composent la population. 
- `nbMaxGen` : l'algorithme fait au maximum `nbMaxGen` itérations. 
- `init` : fonction qui génère des solutions, utilisée pour initialiser la 
   population. De prototype `() -> ::G`.
- `z` : la fonction objectif. De prototype (::G) -> (::Y).
- `fCV` : la fonction de dépassement de contrainte. De prototype (::G) -> Float64.
- `pmut` : proba qu'un individu mute. Lorsque la population enfant est générée, l'opérateur 
   de mutation est appliqué à chaque enfant avec une proba `pmut`.
- `fmut` : opérateur de mutation. De prototype `(x::G) -> ()`. 
- `fcross` : opérateur de croisement. De prototype `(pa::G, pb::G, ca::G, cb::G) -> ()', il 
   croise les solutions `pa` et `pb`, les solutions enfants obtenues sont stockées dans `ca`
   et `cb`.
- `seed` : Array de solutions (de type `G` ou qui peuvent être convertit en `G` via convert) 
   pour initialiser une partie de la population. 
- `critarret` : un critère d'arrêt pour l'algorithme. De prototype
   `critarret(donneesCrit::Dict, pop::AbstractVector{<:Indiv{G, Y}} where {G, Y}) 
       -> Bool` 
   où `donneesCrit` est un dictionnaire qui sera créé par `_nsga` dans lequel peuvent être 
   stockée des données nécessaires pour le fonctionnement du critère et `pop` est la 
   population courante. 
   `criArret` est appelé à chaque itération. `_nsga` s'arrête quand il renvoie `true`.
- `fdiag` : Fonction qui est appelée toute les `periodediag` générations. Les objects 
   `donneesCrit`, le vecteur de la population courante `pop` l'itération courante `ite` sont
   passé comme en arguments positionnels.

# Return 
Un vecteur d'individus (de type `Indiv{G, Y}`, où `Y` est le type d'une
variable renvoyée par la fonction `z`) représentant la population obtenue à la fin de
l'exécution de l'algorithme est renvoyé. 
"""
function _nsga(::Indiv{G, Y}, sense, popSize, nbmaxgen, init, z, 
    fCV, pmut, fmut, fcross, seed, critarret, fdiag, periodediag, 
    refreshtime)::Vector{Indiv{G, Y}} where {G,Y}

   # Données pour le critère d'arrêt. 
   donneescrit = Dict()

   popSize = max(popSize, length(seed))
   isodd(popSize) && (popSize += 1)
   P = Vector{Indiv{G, Y}}(undef, 2 * popSize)
   P[1:(popSize - length(seed))] .= [createIndiv(init(), z, fCV) for _ = 1:(popSize - length(seed))]
   for i = 1:length(seed)
      P[popSize - length(seed) + i] = createIndiv(convert(G, seed[i]), z, fCV)
      if fCV(P[popSize - length(seed) + i].x) > 0
         @warn "element $i of the seed is unfeasible"
      end
   end
   # Ici, on obtient Q_0 en faisant une copie de P_0. Dans [Deb K. et al., 2000], Q_0 est 
   # obtenue en faisant évoluer P_0.
   for i = 1:popSize
      P[popSize+i] = deepcopy(P[i])
   end
   fast_non_dominated_sort!(view(P, 1:popSize), sense)

   boolarret = critarret(donneescrit, P)

   fdiag(;ite = 0, pop = P, donneesCrit = donneescrit)

   @showprogress refreshtime for gen = 1:nbmaxgen
      for i = 1:2:popSize

         pa = tournament_selection(P)
         pb = tournament_selection(P)

         crossover!(pa, pb, fcross, P[popSize + i], P[popSize + i + 1])

         rand() < pmut && mutate!(P[popSize + i], fmut)
         rand() < pmut && mutate!(P[popSize + i + 1], fmut)

         eval!(P[popSize + i], z, fCV)
         eval!(P[popSize + i + 1], z, fCV)
      end

      fast_non_dominated_sort!(P, sense)
      sort!(P, by = x -> x.rank, alg = Base.Sort.QuickSort)
        
      let f::Int = 1
         ind = 0
         indnext = findlast(x -> x.rank == f, P)
         while 0 < indnext <= popSize
            ind = indnext
            f += 1
            indnext = findlast(x -> x.rank == f, P)
         end
         indnext == 0 && (indnext = length(P))
         crowding_distance_assignment!(view(P, ind+1:indnext))
         sort!(view(P, (ind + 1):indnext), by = x -> x.crowding, rev = true, 
               alg = PartialQuickSort(popSize - ind))
      end

      boolarret = critarret(donneescrit, P) 
      boolarret && break

      gen % periodediag == 0 && fdiag(;ite = gen, pop = P, donneesCrit = donneescrit)

   end

   !boolarret && @warn "Le nombre maximum d'itérations ($nbmaxgen) a été atteint
   sans que le critère de convergence ne soit satisfait"

   # On renvoie uniquement les individus non dominés
   filter(x -> x.rank == 1, view(P, 1:popSize))
end

""" 
Etant donnée une population d'individu de type `Indiv`, calculer le niveau de non-domination
de chacun des ces individus avec l'algorithme décrit dans [Deb K. et al., 2000]. Le niveau
est stocké dans le champs `rank`. Les champs `dom_count` et `dom_list` sont aussi utilisé
pour stoker des données utiles lors de l'exécution (voir [Deb K. et al., 2000]).
"""
function fast_non_dominated_sort!(pop::AbstractVector{T}, sense) where {T}
    n = length(pop)

    for p in pop
        empty!(p.dom_list)
        p.dom_count = 0
        p.rank = 0
    end

    @inbounds for i in 1:n
        for j in i+1:n
            if dominates(sense, pop[i], pop[j])
                push!(pop[i].dom_list, j)
                pop[j].dom_count += 1
            elseif dominates(sense, pop[j], pop[i])
                push!(pop[j].dom_list, i)
                pop[i].dom_count += 1
            end
        end
        if pop[i].dom_count == 0
            pop[i].rank = 1
        end
    end

    k = UInt16(2)
    @inbounds while any(==(k-one(UInt16)), (p.rank for p in pop)) #ugly workaround for #15276
        for p in pop 
            if p.rank == k-one(UInt16)
                for q in p.dom_list
                    pop[q].dom_count -= one(UInt16)
                    if pop[q].dom_count == zero(UInt16)
                        pop[q].rank = k
                    end
                end
            end
        end
        k += one(UInt16)
    end
    nothing
end

function crowding_distance_assignment!(pop::AbstractVector{Indiv{X, NTuple{N, T}}}) where {X, N, T}
   # Question : le code semble sous-entendre que lorsque l'on trie par ordre croissant sur 
   # le premier objectif, on trie par ordre décroissant sur le deuxième. C'est vraie si les
   # objectifs sont contradictoires et que l'on a des solutions proches du front de Pareto
   # mais pas toujours vrai en général. Ainsi, dans le cas N = 2, il me semble que la
   # fonction ci-dessous peut renvoyer un valeur approximative. Est-ce un problème ? 
   # Est-ce que cela change quelque chose ?
    if N == 2
        sort!(pop, by = x -> x.y[1])
        pop[1].y[1] == pop[end].y[1] && return # Don't waste time if all indivs are the same
        pop[1].crowding = pop[end].crowding = Inf

        width_y1 = (pop[end].y[1] - pop[1].y[1])
        width_y2 = (pop[1].y[2] - pop[end].y[2])
        @inbounds for i = 2:length(pop)-1
            pop[i].crowding = (pop[i+1].y[1] - pop[i-1].y[1]) / width_y1 + (pop[i-1].y[2] - pop[i+1].y[2]) / width_y2
        end
    else
        for ind in pop
            ind.crowding = 0.
        end
        @inbounds for j = 1:length(first(pop).y) # Foreach objective
            let j = j #https://github.com/JuliaLang/julia/issues/15276
                sort!(pop, by = x -> x.y[j]) #sort by the objective value
            end
            pop[1].crowding = pop[end].crowding = Inf #Assign infinite value to extremas
            if pop[1].y[j] != pop[end].y[j]
                for i = 2:length(pop)-1
                    pop[i].crowding += (pop[i+1].y[j] - pop[i-1].y[j]) / (pop[end].y[j] - pop[1].y[j])
                end
            end
        end
    end
end

"""
Sélection par tournoi : deux individus sont tirés au hasard dans la population P. Le 
meilleur selon l'ordre "Crowded comparison operator " dans [Deb K. et al., 2000] est 
renvoyé. (dans ce package, cet ordre est défini dans le fichier `indiv.jl`.)
"""
function tournament_selection(P)
    a, b = rand(1:length(P)÷2), rand(1:length(P)÷2)
    P[a] < P[b] ? P[a] : P[b]
end


"""
Arrêter _nsga2 lorsque l'écart-type de l'objectif divisé par `nbMaxSat` est inférieur à un
certain seuil.  Ce critère n'est pertinent que si l'objectif est unidimensionnel. 

Ce critère calcule l'écart type de la valeur de l'objectif du meilleur individu la population
prise entre l'itération où la réalisabilité est atteinte (cad il existe au moins un individu
réalisable dans la population) et l'itération courante. Il renvoie true si cette valeur 
divisée par `nbMaxSat` est inférieur à `ϵ` et si `nbMinGen` générations se sont écoulées
depuis l'atteinte de la réalisabilité. 
"""
function critArretEC(donneesCrit::Dict, pop, ϵ::Real, nbMinGen::Int, nbMaxSat::Int)
   if donneesCrit == Dict()
      if pop[1].CV ≈ 0 
         donneesCrit[:iteDepuisRea] = 0
         donneesCrit[:moyObj] = pop[1].y[1]
         donneesCrit[:moyCarreObj] = (pop[1].y[1])^2 
      end
      return false
   end
   donneesCrit[:iteDepuisRea] += 1
   n = donneesCrit[:iteDepuisRea] 
   donneesCrit[:moyObj] = (n*donneesCrit[:moyObj] + pop[1].y[1])/(n+1) 
   donneesCrit[:moyCarreObj] = 
      (n*donneesCrit[:moyCarreObj] + pop[1].y[1]^2)/(n+1) 

   donneesCrit[:iteDepuisRea] < nbMinGen && return false 

   return (√(donneesCrit[:moyCarreObj] - donneesCrit[:moyObj]^2))/nbMaxSat < ϵ
end

"""
Arrêter _nsga2 lorsque l'écart-type de l'objectif divisé par `nbMaxSat`, pris sur les
`nbGenGlissant` dernières générations, est inférieur à ϵ. Ce critère n'est pertinent que si
l'objectif est unidimensionnel. On ne considère les valeurs de l'objectif qu'une fois que la 
réalisabilité a été atteinte. 
"""
function critArretECGlissant(donneesCrit::Dict, pop, ϵ::Real, nbGenGlissant::Int, 
      nbMaxSat::Int)
   
   # Il y a trois états possibles qui seront atteints successivement
   # 1- la réalisabilité n'a pas encore été atteinte.
   # 2- la réalisabilité est atteinte depuis moins de nbGenGlissant générations
   # 3- la réalisabilité est atteinte depuis plus de nbGenGlissant générations

   # Si on est dans l'état 1
   if donneesCrit == Dict()
      if pop[1].CV ≈ 0
         donneesCrit[:valsObj] = [pop[1].y[1]]
         sizehint!(donneesCrit[:valsObj], nbGenGlissant) 
      end
      return false 
   # Sinon, si on est dans l'état 2
   elseif length(donneesCrit[:valsObj]) < nbGenGlissant 
      push!(donneesCrit[:valsObj], pop[1].y[1])
      return false 
   # Sinon, on est dans l'état 3
   else 
      popfirst!(donneesCrit[:valsObj])
      push!(donneesCrit[:valsObj], pop[1].y[1])
      return std(donneesCrit[:valsObj])/nbMaxSat < ϵ 
   end
end





