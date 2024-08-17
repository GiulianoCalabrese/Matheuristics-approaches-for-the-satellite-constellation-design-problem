"""
Représentation d'un individu sur laquelle travaille la fonction `_nsga`.

# Champs 
- `x` : la solution que représente l'individu. Le type `G` peut être quelconque. Cependant, 
   il faudra alors fournir à _nsga des opérateurs adaptée. 
- `y` : la valeur de l'objectif associé à la solution `y`. Le type `Y` doit être itérable.
   Chaque élément de `y` représente une composante de l'objectif. Il doit donc exister une 
   méthode `<` pour le type des éléments de `y`. En pratique, `Y` est le type d'un tuple 
   dont les composantes sont réelles. 
- `CV` : La valeur de la fonction de dépassement de contrainte associée à la solution `x`. 
- `rank` : le niveau de non-domination de la solution `x` relativement à la population sur 
   laquelle travaille `_nsga`.
- `crowding` : la distance de peuplement de la solution `x` relativement à la population sur
   laquelle travaille `_nsga`.
- `dom_count` et `dom_list` : des variables utilisées par la fonction 
   `fast_non_dominated_sort!`. Voir [Deb K. et al., 2000] pour leur utilité précise. 
"""
mutable struct Indiv{G, Y} # Genotype, Type(Y_N)
    x::G
    y::Y
    CV::Float64
    rank::UInt16
    crowding::Float64
    dom_count::UInt16
    dom_list::Vector{UInt16}
    Indiv(x::G, y::Y, cv) where {G, Y} = new{G, Y}(x, y, cv, zero(UInt16), 0., zero(UInt16), UInt16[])
end

"""
Initialiser un individu à partir d'une solution `x`. 

# Arguments
- `x` : une solution. 
- `z` : la fonction objective. Si `G` est le type de `x`, le prototype de `z` doit être
   (x::G) -> (y::Y) où `Y` est un type quelconque. 
- `fCV` : la fonction de dépassement de contrainte. La prototype de `fCV` doit être 
   (x::G) -> Float64.

# Return 
   Un individu de type `Indiv{G, Y}`
"""
function createIndiv(x, z, fCV)
    # MODIF : j'appelle d'abord la fonction de dépassement de contrainte car cela appelle
    # decode!, ce qui met à jour le phénotype x.p. (Dans le cas où l'on travaille avec des 
    # solutions codées avec BinaryCoding ou BinaryVariableCoding.)
    # y = z(x)
    # cv = fCV(x)
    cv = fCV(x)
    y = z(x)
    Indiv(x, y, cv)
end

struct Max end
struct Min end

"""
Renvoie true ssi l'individu `a` domine l'individu `b` au sens de la minimisation. La notion
de domination utilisée est définie dans le rapport. 
"""
function dominates(::Min, a::Indiv, b::Indiv)
    a.CV != b.CV && return a.CV < b.CV
    res = false
    for i in eachindex(a.y)
        @inbounds a.y[i] > b.y[i] && return false
        @inbounds a.y[i] < b.y[i] && (res = true)
    end
    res
end

"""
Renvoie true ssi l'individu `a` domine l'individu `b` au sens de la maximisation. La notion
de domination utilisée est définie dans le rapport. 
"""
function dominates(::Max, a::Indiv, b::Indiv)
    a.CV != b.CV && return a.CV < b.CV
    res = false
    for i in eachindex(a.y)
        @inbounds a.y[i] < b.y[i] && return false
        @inbounds a.y[i] > b.y[i] && (res = true)
    end
    res
end

Base.:(==)(a::Indiv, b::Indiv) = a.x == b.x
Base.hash(a::Indiv) = hash(a.x)
Base.isless(a::Indiv, b::Indiv) = a.rank < b.rank || a.rank == b.rank && a.crowding >= b.crowding #Comparison operator for tournament selection
Base.show(io::IO, ind::Indiv) = print(io, "Indiv($(repr_pheno(ind.x)) : $(ind.y) | rank : $(ind.rank))")
repr_pheno(x) = repr(x)
function repr_pheno(x::Union{BitVector, Vector{Bool}}) 
    res = map(x -> x ? '1' : '0', x)
    if length(res) <= 40
        return "["*String(res)*"]"
    else
        return "["*String(res[1:15])*"..."*String(res[end-14:end])*"]"
    end
end

function eval!(indiv::Indiv, z::Function, fCV::Function)
    indiv.CV = fCV(indiv.x)
    indiv.CV ≈ 0 && (indiv.y = z(indiv.x))
    indiv
end
