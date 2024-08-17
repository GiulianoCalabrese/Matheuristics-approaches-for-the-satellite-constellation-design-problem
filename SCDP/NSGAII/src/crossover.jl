function crossover!(ind_a, ind_b, fcross, child_a, child_b)
    fcross(ind_a.x, ind_b.x, child_a.x, child_b.x)
end
(default_crossover!(pa::T, pb::T, ca, cb)) where T<:AbstractVector{Bool} = two_point_crossover!(pa, pb, ca, cb)
(default_crossover!(pa::T, pb::T, ca, cb)) where T<:AbstractVector{Int} = PMX_crossover!(pa, pb, ca, cb)

function two_point_crossover!(bits_a, bits_b, child1, child2)
    cut_a = cut_b = rand(2:length(bits_a)-1)
    while(cut_b == cut_a)
        cut_b = rand(2:length(bits_a))
    end
    cut_a,cut_b = minmax(cut_b,cut_a)

    copyto!(child1, 1, bits_a, 1, cut_a-1)
    copyto!(child1, cut_a, bits_b, cut_a, cut_b-cut_a+1)
    copyto!(child1, cut_b+1, bits_a, cut_b+1, length(bits_a)-cut_b)

    copyto!(child2, 1, bits_b, 1, cut_a-1)
    copyto!(child2, cut_a, bits_a, cut_a, cut_b-cut_a+1)
    copyto!(child2, cut_b+1, bits_b, cut_b+1, length(bits_a)-cut_b)
end

"""
Croisement à un point pour les génotypes de longueur variable. 
Cet opérateur est décrit dans mon rapport. 
"""
function one_point_crossover!(pa::BinaryVariableCodedIndiv, pb::BinaryVariableCodedIndiv, 
      ca::BinaryVariableCodedIndiv, cb::BinaryVariableCodedIndiv, bvc::BinaryVariableCoding)

   # Nombre de bits par bloc
   nbBB = bvc.bc.nbbitstotal

   cut_a = rand(0:pa.nbBlocs[])
   cut_b = rand(0:pb.nbBlocs[])


   # Construction du premier enfant 
   nbBlocs_ca = cut_a + pb.nbBlocs[] - cut_b
   if nbBlocs_ca ≤ bvc.nbMaxBlocs
      ca.nbBlocs[] = nbBlocs_ca
      copyto!(ca.x, 1, pa.x, 1, cut_a*nbBB)
      copyto!(ca.x, cut_a*nbBB+1, pb.x, cut_b*nbBB+1, (ca.nbBlocs[] - cut_a)*nbBB)
   # Si le premier enfant a un nombre de blocs strictement plus grand que le nombre maximal
   # de bloc, on lui supprime des blocs jusqu'à ce qu'il comporte le nombre maximal de bloc. 
   else 
      ca.nbBlocs[] = bvc.nbMaxBlocs
      # Sélection des blocs que l'on conserve 
      numsBloc = sort(sample(1:nbBlocs_ca, bvc.nbMaxBlocs, replace = false))
      for indBloc in 1:bvc.nbMaxBlocs
         if numsBloc[indBloc] ≤ cut_a 
            copyto!(ca.x, (indBloc-1)*nbBB + 1, pa.x, (numsBloc[indBloc]-1)*nbBB + 1, nbBB)
         else
            copyto!(ca.x, (indBloc-1)*nbBB + 1, pb.x, (numsBloc[indBloc]-cut_a-1)*nbBB + 1,
                    nbBB)
         end
      end
   end

   # Construction du deuxième enfant 
   nbBlocs_cb = cut_b + pa.nbBlocs[] - cut_a
   if nbBlocs_cb ≤ bvc.nbMaxBlocs
      cb.nbBlocs[] = nbBlocs_cb
      copyto!(cb.x, 1, pb.x, 1, cut_b*nbBB)
      copyto!(cb.x, cut_b*nbBB+1, pa.x, cut_a*nbBB+1, (cb.nbBlocs[] - cut_b)*nbBB)
   # Si le premier enfant a un nombre de blocs strictement plus grand que le nombre maximal
   # de bloc, on lui supprime des blocs jusqu'à ce qu'il comporte le nombre maximal de bloc. 
   else 
      cb.nbBlocs[] = bvc.nbMaxBlocs
      # Sélection des blocs que l'on conserve 
      numsBloc = sort(sample(1:nbBlocs_cb, bvc.nbMaxBlocs, replace = false))
      for indBloc in 1:bvc.nbMaxBlocs
         if numsBloc[indBloc] ≤ cut_b 
            copyto!(cb.x, (indBloc-1)*nbBB + 1, pb.x, (numsBloc[indBloc]-1)*nbBB + 1, nbBB)
         else
            copyto!(cb.x, (indBloc-1)*nbBB + 1, pa.x, (numsBloc[indBloc]-cut_b-1)*nbBB + 1,
                    nbBB)
         end
      end
   end
end

function one_point_crossover_old!(pa::BinaryVariableCodedIndiv, pb::BinaryVariableCodedIndiv, 
      ca::BinaryVariableCodedIndiv, cb::BinaryVariableCodedIndiv, bvc::BinaryVariableCoding)
   cut_a = rand(0:pa.nbBlocs[])
   cut_b = rand(0:pb.nbBlocs[])

   ca.nbBlocs[] = min(bvc.nbMaxBlocs, cut_a + pb.nbBlocs[] - cut_b)
   copyto!(ca.x, 1, pa.x, 1, cut_a*bvc.bc.nbbitstotal)
   copyto!(ca.x, cut_a*bvc.bc.nbbitstotal+1, pb.x, cut_b*bvc.bc.nbbitstotal+1, 
           (ca.nbBlocs[] - cut_a)*bvc.bc.nbbitstotal)

   cb.nbBlocs[] = min(bvc.nbMaxBlocs, cut_b + pa.nbBlocs[] - cut_a)
   copyto!(cb.x, 1, pb.x, 1, cut_b*bvc.bc.nbbitstotal)
   copyto!(cb.x, cut_b*bvc.bc.nbbitstotal+1, pa.x, cut_a*bvc.bc.nbbitstotal+1, 
           (cb.nbBlocs[] - cut_b)*bvc.bc.nbbitstotal)

end

function PMX_crossover!(pa, pb, ca, cb)
    cut_a = cut_b = rand(1:length(pa))
    while(cut_b == cut_a)
        cut_b = rand(1:length(pa))
    end
    cut_a, cut_b = minmax(cut_a, cut_b)

    copyto!(ca, 1, pb, 1, cut_a-1)
    copyto!(ca, cut_a, pa, cut_a, cut_b-cut_a+1)
    copyto!(ca, cut_b+1, pb, cut_b+1, length(pa)-cut_b)

    copyto!(cb, 1, pa, 1, cut_a-1)
    copyto!(cb, cut_a, pb, cut_a, cut_b-cut_a+1)
    copyto!(cb, cut_b+1, pa, cut_b+1, length(pa)-cut_b)

    @inbounds for i = cut_a:cut_b
        if pa[i] ∉ view(pb, cut_a:cut_b)
            j = findfirst(isequal(pb[i]), pa)
            while j ∈ cut_a:cut_b
                j = findfirst(isequal(pb[j]), pa)
            end
            cb[j] = pa[i]
        end

        if pb[i] ∉ view(pa, cut_a:cut_b)
            j = findfirst(isequal(pa[i]), pb)
            while j ∈ cut_a:cut_b
                j = findfirst(isequal(pa[j]), pb)
            end
            ca[j] = pb[i]
        end
    end
end

# Test 
# types = [:Bool, :Cont, :Int]
# lb = [0, 3.3, 4]
# ub = [1, 5.7, 10]
# bvc = BinaryVariableCoding(BinaryCoding(2, types, lb, ub), 3)

# p1 = [0, 3.4, 5,
#      1, 4., 5,
#      1, 5.1, 10]

# p2 = [1, 3.9, 4,
#       0, 4.8, 7]

# ind1 = encode(p1, 3, bvc)
# ind2 = encode(p2, 2, bvc)
# ch1 = encode([], 0, bvc)
# ch2 = encode([], 0, bvc)
# print(one_point_crossover!(ind1, ind2, ch1, ch2, bvc))

# mat = hcat(ind1.x, ind2.x, ch1.x, ch2.x)
# pretty_table(mat, show_row_number = true)
