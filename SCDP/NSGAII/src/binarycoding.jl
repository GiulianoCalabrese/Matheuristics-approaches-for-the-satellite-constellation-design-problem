struct BinaryCoding
    nbvar::Int
    types::Vector{Symbol}
    lb::Vector{Float64}
    ub::Vector{Float64}
    nbbits::Vector{Int}
    nbbitstotal::Int
end

struct BinaryCodedIndiv
    x::BitVector
    p::Vector{Float64}
end

struct BinaryVariableCoding 
   # Codage d'un bloc
   bc::BinaryCoding
   # Nombre maximal de blocs qu'une solution peut contenir.
   nbMaxBlocs::Int
end

"""
Une structure pour représenter des individus codés en binaire avec un nombre non constant 
de variable contenu dans l'individu
Bloc de base est représenté avec BinaryCoding.
Il existe un nombre maximal de blocs donc `x` et `p` ont une taille maximale fixée. Seul
l'information correspondant aux `nbBlocs` premiers blocs dans `x` et `p` a un sens. 
"""
struct BinaryVariableCodedIndiv
   x::BitVector
   p::Vector{Float64}
   nbBlocs::Base.RefValue{Int64}
end

function BinaryCoding(ϵ::Int, types::Vector{Symbol}, lb, ub)
    @assert length(types) == length(lb) == length(ub)
    @assert all(lb .< ub)
    @assert UInt128(10)^ϵ <= UInt128(2)^127

    nbvar = length(lb)
    nbbits = ones(Int, nbvar)
    for i = 1:nbvar
        if types[i] == :Cont
            while ((UInt128(10)^ϵ)*(ub[i]-lb[i]) >= UInt128(2)^nbbits[i])
                nbbits[i] += 1
            end
        elseif types[i] == :Int
            while 2^nbbits[i] <= ub[i]-lb[i]
                nbbits[i] += 1
            end
        end
    end
    BinaryCoding(nbvar, types, lb, ub, nbbits, sum(nbbits))
end
BinaryCoding(ϵ::Int, lb, ub) = BinaryCoding(ϵ, fill(:Cont, length(lb)), lb, ub)

"""
Cette fonction décode dans `p` le tronçon de `x` codant le `ind_bloc`ième bloc de `p` selon 
le codage binaire `bc`.
"""
function decode!(p::AbstractVector, x::BitVector, ind_bloc::Int, bc::BinaryCoding)
   ind_p = (ind_bloc-1)*bc.nbvar
   ind_x = (ind_bloc-1)*bc.nbbitstotal

   for i = 1:bc.nbvar
      ind_x += bc.nbbits[i]
      if bc.types[i] == :Bin
         p[ind_p + i] = x[ind_x] == 1 ? 1. : 0.
      else
         val = zero(UInt128)
         puis = one(UInt128)
         for ind = ind_x:-1:ind_x-bc.nbbits[i]+1
             x[ind] && (val += puis)
             puis *= 2
         end

         if bc.types[i] == :Cont
           p[ind_p + i] = bc.lb[i] + val * (bc.ub[i] - bc.lb[i]) / (UInt128(2)^bc.nbbits[i] - 1)
         else
            p[ind_p + i] = val + bc.lb[i]
         end

      end
   end
end

"""
Mettre à jour le phénotype de l'individu `indiv` en décodant son génotype.
"""
decode!(indiv::BinaryCodedIndiv, bc::BinaryCoding) = decode!(indiv.p, indiv.x, 1, bc)

"""
Mettre à jour le phénotype de l'individu `indiv` en décodant son génotype.
"""
decode!(indiv::BinaryVariableCodedIndiv, bvc::BinaryVariableCoding) = 
   decode!.(Ref(indiv.p), Ref(indiv.x), 1:indiv.nbBlocs[], Ref(bvc.bc))

"""
Coder dans `x` selon `bc` les variables de `p` correspondant au `ind_bloc`ième bloc.
"""
function encode!(p::AbstractVector, x::BitVector, bc::BinaryCoding, ind_bloc::Int)
   ind_p = (ind_bloc-1)*bc.nbvar
   ind_x = (ind_bloc-1)*bc.nbbitstotal+1

   for i = 1:bc.nbvar
      if bc.types[i] == :Int
         tab = reverse(digits(Bool, round(Int, p[ind_p + i] - bc.lb[i]), base = 2, 
                             pad = bc.nbbits[i]))
      elseif bc.types[i] == :Bin
         tab = p[ind_p + i] != 0
      else
         t = (p[ind_p + i] - bc.lb[i]) / (bc.ub[i] - bc.lb[i]) * (UInt128(2)^bc.nbbits[i] - 1)
         target = round(UInt128, t)
         if target == UInt128(2)^bc.nbbits[i] - 1
            target -= 1
         end
         tab = reverse(digits(Bool, target, base = 2, pad = bc.nbbits[i]))
      end
      x[ind_x : ind_x+bc.nbbits[i]-1] .= tab
      ind_x += bc.nbbits[i]
   end
end

"""
Coder, telque définie par le code binaire,  le phénotype `p` en binaire dans le cas où 
l'on travaille avec des individus définis par un nombre fixe de variables. 
"""
function encode(p::AbstractVector, bc::BinaryCoding)::BinaryCodedIndiv
   binaryIndiv = BinaryCodedIndiv(BitVector(undef, bc.nbbitstotal), p)
   encode!(p, binaryIndiv.x, bc, 1)
   return binaryIndiv
end
 
# encode(p, d) = encode([p], d)

"""
Coder en binaire le phénotype `p` dans le cas où l'on travaille avec des individus définis
par un nombre non fixe de variables. Le codage est défini par `bvc`.
"""
function encode(p::AbstractVector, nbBlocs::Int, bvc::BinaryVariableCoding)::BinaryVariableCodedIndiv

   @assert nbBlocs ≤ bvc.nbMaxBlocs
   @assert length(p) == nbBlocs*bvc.bc.nbvar

   binaryIndiv = BinaryVariableCodedIndiv(
                                       BitVector(undef, bvc.nbMaxBlocs*bvc.bc.nbbitstotal), 
                                       similar(p, bvc.nbMaxBlocs*bvc.bc.nbvar),
                                       Ref(nbBlocs))

   binaryIndiv.p[1:bvc.bc.nbvar*nbBlocs] .= p
   encode!.(Ref(p), Ref(binaryIndiv.x), Ref(bvc.bc), 1:nbBlocs)

   return binaryIndiv
end

#  ----------------------------------------------------------------------------------

"""
Fonction `encode` d'origine de nsga2. Je l'ai modifiée simplement pour éviter les 
redondances de code. 
"""
function encode_old(p::AbstractVector, d::BinaryCoding)::BinaryCodedIndiv
    res = BinaryCodedIndiv(BitVector(), p)
    sizehint!(res.x, d.nbbitstotal)
    for i = 1:d.nbvar
        if d.types[i] == :Int
            tab = reverse(digits(Bool, round(Int, res.p[i] - d.lb[i]), base = 2, pad = d.nbbits[i]))
            append!(res.x, tab)
        elseif d.types[i] == :Bin
            push!(res.x, p[i]!=0)
        else
            t = (p[i] - d.lb[i]) / (d.ub[i] - d.lb[i]) * (UInt128(2)^d.nbbits[i] - 1)
            target = round(UInt128, t)
            if target == UInt128(2)^d.nbbits[i] - 1
                target -= 1
            end
            tab = reverse(digits(Bool, target, base = 2, pad = d.nbbits[i]))
            append!(res.x, tab)
        end
    end
    res
end

"""
Fonction `decode!` d'origine de nsga2. Je l'ai modifiée simplement pour éviter les 
redondances de code. 
"""
function decode_old!(indiv::BinaryCodedIndiv, d::BinaryCoding)
    j = 0
    for i = 1:d.nbvar
        j += d.nbbits[i]
        if d.types[i] == :bin
            indiv.p[i] = indiv.x[j] == 1 ? 1. : 0.
        else
            val = zero(uint128)
            puis = one(uint128)
            for ind = j:-1:j-d.nbbits[i]+1
                indiv.x[ind] && (val += puis)
                puis *= 2
            end

            if d.types[i] == :cont
                indiv.p[i] = d.lb[i] + val * (d.ub[i] - d.lb[i]) / (uint128(2)^d.nbbits[i] - 1)
            else
                indiv.p[i] = val + d.lb[i]
            end
        end
    end
    indiv
end

"""
Plus besoin ... 
Ajouter au vecteur de bits `bitVector` la représentation binaire du phénotype `p` telque 
définie par le codage binaire `d`. 
"""
function phenToBitVector!(p::AbstractVector, x::BitVector, d::BinaryCoding)
    for i = 1:d.nbvar
        if d.types[i] == :Int
            tab = reverse(digits(Bool, round(Int, p[i] - d.lb[i]), base = 2, pad = d.nbbits[i]))
            append!(bitVector, tab)
        elseif d.types[i] == :Bin
            push!(bitVector, p[i]!=0)
        else
            t = (p[i] - d.lb[i]) / (d.ub[i] - d.lb[i]) * (UInt128(2)^d.nbbits[i] - 1)
            target = round(UInt128, t)
            if target == UInt128(2)^d.nbbits[i] - 1
                target -= 1
            end
            tab = reverse(digits(Bool, target, base = 2, pad = d.nbbits[i]))
            append!(bitVector, tab)
        end
    end
 end

# Tests
# types = [:Bool, :Cont, :Cont, :Int]
# lb = [0, 0, 3.3, 4]
# ub = [1, 2.1, 5.7, 10]
# bcVar = BinaryVariableCoding(BinaryCoding(3, types, lb, ub), 4)
# bc = BinaryCoding(3, repeat(types, outer = 3), repeat(lb, outer = 3), 
#                      repeat(ub, outer = 3))
# p = [0, 1.2, 3.4, 5,
#      1, 0.2, 4., 5,
#      1, 2., 5.1, 10]

# ind1 = encode(p, bc)
# ind2 = encode(p, 3, bcVar)
# test1 = ind2.x[1:bcVar.bc.nbbitstotal*3] == ind1.x

# for i in 1:bcTot.nbbitstotal
#    rand() < 0.6 && (ind1.x[i] = !ind1.x[i]; ind2.x[i] = !ind2.x[i])
# end

# decode!(ind1, bc)
# decode!(ind2, bcVar)

