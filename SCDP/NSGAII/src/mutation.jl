mutate!(ind::Indiv, fmut!) = fmut!(ind.x)
default_mutation!(p::Vector{Int}) = rand_swap!(p)
default_mutation!(b::T) where T<:AbstractVector{Bool} = rand_flip!(b)

function rand_flip!(bits)
    nb = length(bits)
    for i = 1:nb
        if rand() < 1/nb
            @inbounds bits[i] = 1 - bits[i]
        end
    end
end

"""
Opérateur de mutation pour les génotypes codée en longueur variable. 
Cet opérateur est décrit dans mon rapport. 

Attention, j'ai mal choisi mes noms de variables : ind signifie à la fois individu et 
indice...
"""
function rand_flip!(ind::BinaryVariableCodedIndiv, bvc::BinaryVariableCoding, 
      probaModifNbBlocs::Real)
   
   #
   # println("probaModifNbBlocs = $probaModifNbBlocs")

   # Nombre de bits par blocs
   nbBB = bvc.bc.nbbitstotal 

   # Le nombre de bloc d'un individu mute avec une probabilité probaModifNbBlocs
   if rand() < probaModifNbBlocs 
      nouvNbBlocs = rand(0:bvc.nbMaxBlocs)
      
      # Si le nouveau génotype a plus de bloc, on génère aléatoirement des blocs que l'on
      # insère à des endroits choisis aléatoirement sur le génotype jusqu'à atteindre
      # nouvNbBlocs 
      if nouvNbBlocs > ind.nbBlocs[]
         while ind.nbBlocs[] < nouvNbBlocs 
            # Choisir un endroit aléatoire sur le génotype 
            indBloc = rand(0:ind.nbBlocs[])
            # Faire de la place sur le génotype
            copyto!(ind.x, (indBloc+1)*nbBB+1, ind.x, indBloc*nbBB+1, 
                   (ind.nbBlocs[]-indBloc)*nbBB)
            # Ajouter le nouveau bloc
            copyto!(ind.x, indBloc*nbBB+1, bitrand(nbBB), 1, nbBB)

            ind.nbBlocs[] += 1 
         end
      # Si le nouveau génotype a moins de blocs, on retire des blocs choisit aléatoirement.
      elseif nouvNbBlocs < ind.nbBlocs[]
         # Les numéros des blocs que l'on conserve 
         numsBlocs = sort(sample(1:ind.nbBlocs[], nouvNbBlocs, replace = false))
         for indBloc in 1:nouvNbBlocs 
            copyto!(ind.x, (indBloc-1)*nbBB + 1, ind.x, (numsBlocs[indBloc]-1)*nbBB + 1, nbBB)
         end
         ind.nbBlocs[] = nouvNbBlocs
      end
   end

   rand_flip!(view(ind.x, 1:(nbBB*ind.nbBlocs[])))
end

function rand_flip_old!(ind::BinaryVariableCodedIndiv, bvc::BinaryVariableCoding, 
      probaModifNbBlocs::Real)
   
   # Le nombre de bloc d'un individu mute avec une probabilité probaModifNbBlocs
   if rand() < probaModifNbBlocs 
      nouvNbBlocs = rand(0:bvc.nbMaxBlocs)
      if nouvNbBlocs > ind.nbBlocs[]
         ind.x[ind.nbBlocs[]*bvc.bc.nbbitstotal+1 : nouvNbBlocs*bvc.bc.nbbitstotal] .=
         bitrand((nouvNbBlocs-ind.nbBlocs[])*bvc.bc.nbbitstotal)
      end
      ind.nbBlocs[] = nouvNbBlocs
   end

   rand_flip!(view(ind.x, 1:(bvc.bc.nbbitstotal*ind.nbBlocs[])))
end

function rand_swap!(perm::Vector{Int})
    i = j = rand(1:length(perm))
    while j == i
        j = rand(1:length(perm))
    end
    @inbounds perm[i], perm[j] = perm[j], perm[i]
end

# Test 
# types = [:Bool, :Cont, :Int]
# lb = [0, 3.3, 4]
# ub = [1, 5.7, 10]
# bvc = BinaryVariableCoding(BinaryCoding(2, types, lb, ub), 3)

# p1 = [0, 3.4, 5,
#     1, 4., 5,
#     1, 5.1, 10]

#p2 = [1, 3.9, 4,
#      0, 4.8, 7]

#ind1 = encode(p1, 3, bvc)
#ind2 = encode(p2, 2, bvc)
#ch1 = encode([], 0, bvc)
# ch2 = encode([], 0, bvc)
# one_point_crossover!(ind1, ind2, ch1, ch2, bvc)

# mat = hcat(ind1.x, ind2.x, ch1.x, ch2.x)
# pretty_table(mat, show_row_number = true)
