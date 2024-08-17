__precompile__()
module NSGAII

export nsga, critArretECGlissant, BinaryCoding, BinaryCodedIndiv, BinaryVariableCoding, 
BinaryVariableCodedIndiv

using ProgressMeter
using Random
using StatsBase
using LinearAlgebra
using Statistics

include("binarycoding.jl")
include("indivs.jl")
include("functions.jl")
include("crossover.jl")
include("mutation.jl")

# include("vOptWrapper.jl")

function nsga(popSize::Integer, nbMaxGen::Integer, z::Function, init::Function ; 
    fCV = x -> 0., pmut = 0.05, fmut = default_mutation!, fcross = default_crossover!, 
    seed = typeof(init())[], critArret = (donneesCrit, P) -> false, fDiag = x -> nothing, 
    periodeDiag = 1, showprogress = true)

    X = createIndiv(init(), z, fCV)
    
    return _nsga(X, Min(), popSize, nbMaxGen, init, z, fCV, pmut, fmut, fcross, seed, 
                 critArret, fDiag, periodeDiag, showprogress ? 0.5 : Inf)
end

"""
Utiliser l'algorithme NSGAII pour résoudre un problème d'optimisation multiobjectif avec des 
solutions codées en binaires et de longueur fixée.
Une solution est représenté par une structure de type `BinaryCodedIndiv` et le codage en
binaire est défini par le paramètre `bc`.

# Arguments 
- `popSize`
- `nbGen`
- `z` : fonction objectif. De prototype `p::Vector{Float64} -> Y` où `p` est le phénotype 
   d'une solution. Voir la structure `BinaryCodedIndiv` qui est utilisée dans ce cas pour 
   représenter une solution. 
- `bc` : Le codage en binaire des solutions. 

- `fCV` : Fonction de dépassement de contrainte. De prototype 
   `p::Vector{Float64} -> Float64` où `p` est le phénotype d'une solution.
- `pmut`
- `fmut` : opérateur de mutation. De prototype `(g::BitVector) -> ()` où `g` est le 
   génotype d'une solution. Voir la structure `BinaryCodedIndiv`.
- `fcross` : opérateur de croisement. De prototype 
   `(pa::BitVector, pb::BitVector, ca::BitVector, cb::BitVector) -> ()', il 
   croise les génotypes `pa` et `pb`, les solutions enfants obtenues sont stockées dans les 
   génotypes `ca` et `cb`.
- `seed` : Vecteur de phénotypes utilisées pour initialiser la population. 

# Return 
Un vecteur d'individus (de type `Indiv{BinaryCodedIndiv, Y}`, où `Y` est le type d'une
variable renvoyée par la fonction `z`) représentant la population obtenue à la fin de
l'exécution de l'algorithme est renvoyé. 
"""
function nsga(popSize::Integer, nbMaxGen::Integer, z::Function, bc::BinaryCoding ; 
    fCV = x -> 0., pmut = 0.05, fmut = indiv -> default_mutation!(indiv.x), 
    fcross = (indivs...) -> default_crossover!(getproperty.(indivs, :x)...), 
    seed = Vector{Float64}[], critArret = (donneesCrit, P) -> false, 
    fDiag = (;vargs...) -> nothing, periodeDiag = 1, showprogress = true)

    init = () -> BinaryCodedIndiv(bitrand(bc.nbbitstotal), zeros(bc.nbvar))
    
    # MODIF 
    # C'est la fonction _fCV qui met à jour le phénotype indiv.p, et non la fonction _z. 
    # C'est pour être cohérent avec l'appel de eval! 
    # _z = indiv -> (decode!(indiv, bc) ; z(indiv.p))
    _z = indiv -> z(indiv.p)

    # MODIF  
    # La fonction fCV que l'on passe à nsga prend pour argument un phénotype. Lorsque l'on
    # utilise BinaryCoding, un phénotype est codée sous la forme d'un array contenant les
    # variables. La fonction fCV (et z) que l'on passe à _nsga prend en argument un
    # structure de type Indiv.G
    _fCV = indiv -> (decode!(indiv, bc) ; fCV(indiv.p))

    X = createIndiv(init(), _z, _fCV)
    return _nsga(X, Min(), popSize, nbMaxGen, init, _z, _fCV , pmut, fmut, fcross,
                 encode.(seed, Ref(bc)), critArret, fDiag, periodeDiag, 
                 showprogress ? 0.5 : Inf)
end

"""
Utiliser l'algorithme NSGAII pour résoudre un problème d'optimisation multiobjectif avec des 
solutions codées en binaires et de longueur variable.
Une solution est représenté par une structure de type `BinaryVariableCodedIndiv` et le
codage en binaire est défini par le paramètre `bvc`.

# Arguments 
- `popSize`
- `nbGen`
- `z` : fonction objectif. De prototype `::BinaryVariableCodedIndiv -> Y` 
- `bvc` : Le codage en binaire des solutions. 

- `fCV` : Fonction de dépassement de contrainte. De prototype 
   `BinaryVariableCodedIndiv -> Float64`.
- `pmut`
- `fmut` : opérateur de mutation. De prototype `(::BinaryVariableCodedIndiv) -> ()`.
- `fcross` : opérateur de croisement. De prototype 
   `(pa::BinaryVariableCodedIndiv, pb::BinaryVariableCodedIndiv, 
   ca::BinaryVariableCodedIndiv, cb::BinaryVariableCodedIndiv) -> ()', il 
   croise les solutions `pa` et `pb`, les solutions enfants obtenues sont stockées dans `ca`
   et `cb`.
- `seed` : tuple dont la première composante est un vecteur de phénotype, la deuxième le
   nombre de bloc que possèdent chacun de ces phénotypes. Les phénotypes seront utilisés
   pour initialiser la population. 

# Return 
Un vecteur d'individus (de type `Indiv{BinaryVariableCodedIndiv, Y}`, où `Y` est le type
d'une variable renvoyée par la fonction `z`) représentant la population obtenue à la fin de
l'exécution de l'algorithme est renvoyé. 
"""
function nsga(popSize::Integer, nbMaxGen::Integer, z::Function, bvc::BinaryVariableCoding; 
      probaModifNbBlocs = 0.05,
      fmut = ind -> rand_flip!(ind, bvc, probaModifNbBlocs),
      fcross = (pa, pb, ca, cb) -> one_point_crossover!(pa, pb, ca, cb, bvc),
      fCV = x -> 0., pmut = 0.05, seed = (Vector{Float64}[], Int[]), 
      critArret = (donneesCrit, P) -> false,
      fDiag = (;vargs...) -> nothing, 
      periodeDiag = 1, showprogress = true)

   function init() 
      nbBlocs = rand(0:bvc.nbMaxBlocs)
      binaryIndiv = BinaryVariableCodedIndiv(
                                 BitVector(undef, bvc.nbMaxBlocs*bvc.bc.nbbitstotal), 
                                 Vector{Float64}(undef, bvc.nbMaxBlocs*bvc.bc.nbvar),
                                 Ref(nbBlocs))
      binaryIndiv.x[1:nbBlocs*bvc.bc.nbbitstotal] .= bitrand(nbBlocs*bvc.bc.nbbitstotal)
      return binaryIndiv 
   end
    
   # MODIF  
   # La fonction fCV que l'on passe à nsga prend pour argument un phénotype. Lorsque l'on
   # utilise BinaryCoding, un phénotype est codée sous la forme d'un array contenant les
   # variables. La fonction fCV (et z) que l'on passe à _nsga prend en argument un
   # structure de type Indiv.G
   _fCV = x -> (decode!(x, bvc) ; fCV(x))

   X = createIndiv(init(), z, _fCV)
   return _nsga(X, Min(), popSize, nbMaxGen, init, z, _fCV , pmut, fmut, fcross, 
                encode.(seed..., Ref(bvc)), critArret, fDiag, periodeDiag, 
                showprogress ? 0.5 : Inf)
end

end # module
