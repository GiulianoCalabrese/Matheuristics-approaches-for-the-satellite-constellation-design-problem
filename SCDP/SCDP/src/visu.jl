"""
Visualisation des données intervenant dans le SCDP
"""

alphabet = string.(collect('a':'z'))

"""
Afficher dans la console les données de l'instance `inst`.
"""
function ecrireInst(inst::DonneesSCDP)
   matCoordsCibles = vcat((180coord[na, :]/π for coord in inst.coordsCibles)...)
   printstyled("\n*** Données de l'instance ***\n\n", color = :red)
   printstyled("Coordonnées des $(inst.nbCibles) cibles\n", color =:blue)
   pretty_table(matCoordsCibles, header = (["Lattitude", "Longitude"], ["deg", "deg"]), 
                row_names = inst.nbCibles ≤ 26 ? alphabet[1:inst.nbCibles] : nothing)
   print("Angle d'ouverture : $(180inst.angleOuverture/π) degrés\n")
   print("Instant initial : $(inst.instantInitial)\n")
   print("Nombre de jours de simulation : $(inst.m_simu)\n")
   print("Nombre d'intervalles de surveillance sur la durée de simulation : $(inst.n_0)\n")
   print("Pas de simulation : $(inst.dt) secondes\n")
   print("Nombre de chiffres significatifs (après la virgule) : $(inst.nbCS) \n")
   print("Indice de rayon d'orbite minimal : $(inst.n_min)\n")
   print("Indice de rayon d'orbite maximal : $(inst.n_max)\n")
end

"""
Mettre en forme les données de l'individu `indiv`, solution de l'instance `inst` dans une
matrice. 

# Arguments 
   - extraireParamsOrbitaux : fonction de prototype
   ::NSGAII.Indiv, ::DonneesSCDP -> ::Tuple(Int,
                                            AbstractVector(Float64), 
                                            AbstractVector(Float64), 
                                            AbstractVector(Float64),
                                            AbstractVector(Float64)) 
   qui renvoie les paramètres orbitaux des satellites d'une constellation représentée par un 
   individu. 
"""
function indivToMat(inst::DonneesSCDP, indiv::NSGAII.Indiv,
      extraireParamsOrbitaux::Function)

   nbSats, ns, is, Ωs, Ms = extraireParamsOrbitaux(indiv, inst)
   as = calculRayonOrbite.(ns, inst.m_simu)

   as_km = @. as/1000
   Ts_min =  @. (inst.m_simu*T_terre/ns)/60
   is_deg = @. 180is/π
   Ωs_deg = @. 180Ωs/π
   Ms_deg = @. 180Ms/π

   nums = collect(1:nbSats)
   gauche = vcat([nbSats, indiv.CV, indiv.rank][na, :],
                 fill(".", nbSats-1, 3))
   droite = hcat(nums, ns, as_km, Ts_min, is_deg, Ωs_deg, Ms_deg)

   return hcat(gauche, droite)
end

"""
Afficher dans la console les données de l'individu `indiv`, solution du problème `inst`.

# Arguments 
   - extraireParamsOrbitaux : fonction de prototype
   ::NSGAII.Indiv, ::DonneesSCDP -> ::Tuple(Int,
                                            AbstractVector(Float64), 
                                            AbstractVector(Float64), 
                                            AbstractVector(Float64),
                                            AbstractVector(Float64)) 
   qui renvoie les paramètres orbitaux des satellites d'une constellation représentée par un 
   individu. 
"""
function ecrireIndiv(inst::DonneesSCDP, indiv::NSGAII.Indiv, 
      extraireParamsOrbitaux::Function; kwargs...)

   header = (["Nombre de satellites", "Dépassement contrainte", "Rang", "Numéro", 
              "Indice du rayon d'orbite", "Rayon d'orbite", "Période de rotation", 
              "Inclinaison", "Longitude du noeud ascendant Ω", "Anomalie moyenne initiale"],
            ["", "min", "", "", "", "km", "min", "deg", "deg", "deg"])

   highlighter = Highlighter((data, i, j) -> i==1 && j ≤ 4, crayon"yellow bold")

   pretty_table(indivToMat(inst, indiv, extraireParamsOrbitaux); header = header, 
                crop = :none, columns_width = 20, highlighters = highlighter,
                kwargs...)
end
   
"""
Afficher dans la console la population d'invidus `pop`, chaque individu étant solution du
problème `inst`.

# Arguments 
   - extraireParamsOrbitaux : fonction de prototype
   ::NSGAII.Indiv, ::DonneesProb -> ::Tuple(Int,
                                            AbstractVector(Float64), 
                                            AbstractVector(Float64), 
                                            AbstractVector(Float64),
                                            AbstractVector(Float64)) 
   qui renvoie les paramètres orbitaux des satellites d'une constellation représentée par un 
   individu. 
   - `kwargs` : arguments nommés additionnels qui seront passés à pretty_table.
"""
function ecrirePop(inst::DonneesSCDP, pop::AbstractVector{ <:NSGAII.Indiv}, 
      extraireParamsOrbitaux::Function; kwargs...)

   tab = vcat((indivToMat(inst, indiv, extraireParamsOrbitaux) for indiv in pop)...)

   header = (["Nombre de satellites", "Dépassement contrainte", "Rang", "Numéro", 
              "Indice du rayon d'orbite", "Rayon d'orbite", "Période de rotation", 
              "Inclinaison", "Longitude du noeud ascendant Ω", "Anomalie moyenne initiale"],
            ["", "min", "", "", "", "km", "min", "deg", "deg", "deg"])

   highlighter = Highlighter((data, i, j) -> j ≤ 4 && data[i, j] ≠ ".", crayon"yellow bold")

   printstyled("\n*** Population ***\n", color = :red)
   
   pretty_table(tab; header = header, crop = :none, columns_width = 20, 
                highlighters = highlighter, 
                kwargs...)
end

"""
Ecrire dans la console les intervalles où les cibles sont survolées par les 
satellites de la constellation `conste`. 

Le tableau affiché correspond à la concaténation verticale sur les cibles de la matrice 
renvoyé par `matIntervallesSurvol`.

# Arguments
- `format` : le format des bornes des intervalles. Si `format = :sec`, les bornes des
   intervalles sont en secondes écoulées depuis `inst.instantInitial`. Si `format = :DT`, 
   bornes sont au format `DateTime`.
"""
function ecrireIntervallesSurvol(inst::DonneesSCDP, conste::Conste;
   format::Symbol = :DT)

   @assert inst.instantInitialJulian ≈ conste.date 

   table = vcat((matIntervallesSurvol(inst, conste, indCible, format = format) 
                 for indCible in 1:inst.nbCibles)...)

   header = (vcat(["Cibles"], (["Sat $indSat", ""] for indSat in 1:conste.nbSats)..., 
                  ["Union", ""]),)

   format == :sec && (header = (header..., ["", fill("sec", 2(conste.nbSats+1))...]))

   printstyled("\n*** Intervalles de survol ***\n", color = :red)
   pretty_table(table, header = header, crop = :none)
end

"""
Dessiner les intervalles sur lesquelles les cibles sont survolées par au moins un satellite
de la constellation `conste`.

# Arguments 
- `kwargs` : arguments nommés qui seront passés à plot. 
"""
function dessinerIntervallesSurvol(inst::DonneesSCDP, conste::Conste, δ::Float64 = 1e-3, 
      nbMaxIter::Int = 100; kwargs...)

   @assert inst.instantInitialJulian ≈ conste.date 

   # Propagateurs des satellites de la constellation 
   props = init_orbit_propagator.(Val(:twobody), 
      KeplerianElements.(inst.instantInitialJulian, conste.as, 0, conste.is, conste.Ωs, 0, 
                         conste.Ms))

   θs = calcul_θ.(conste.as, inst.angleOuverture)

   # Vecteur contenant intervalles de survols de chaque ville.
   intervsSurvol = calculIntervallesSurvol.(inst.coordsCiblesECEF, Ref(props), Ref(inst.ts), 
                                           Ref(θs), δ, nbMaxIter) 
   
   yTicks = inst.nbCibles ≤ 26 ? (1:1:inst.nbCibles, alphabet[1:inst.nbCibles]) : nothing
   plt = plot(;size = (1300, 300), grid = false, xlabel = "Temps (s)", ylabel = "Cibles", 
              title = "Intervalles de survol", yticks = yTicks, kwargs...)

   # Dessiner, pour chaque cible, les intervalles où elle est survolée par au moins un
   # satellite. 
   # label = inst.nbCibles ≤ 26 ? alphabet[1:inst.nbCibles][na, :] : nothing
   # Vecteurs des couleurs pour les cibles
   cols = get(darkrainbow, inst.nbCibles == 1 ? 0.5 : LinRange(0, 1, inst.nbCibles))

   for ind_cible in 1:inst.nbCibles
      for inter in intervsSurvol[ind_cible]
         plot!(plt, inter, [ind_cible, ind_cible], color = cols[ind_cible], label = "", 
              linewidth = 15)
      end
   end

   # Dessiner des traits verticaux pour les bornes des intervalles de surveillance. 
   for k in 0:inst.n_0
      plot!(plt, fill(inst.ts[1 + k*inst.nbInstants], 2), [0, inst.nbCibles+1], 
            color = :silver, label = "")
   end

   plot!(plt, ylim = (0.5, inst.nbCibles+0.5))

   return plt 
end

"""
Dessiner sur le plot `plt` la surface de la terre. 

# Arguments
   - kwargs : arguments nommés additionnels qui seront passé à plot. 
"""
function dessinerTerre!(plt; kwargs...)
   surf = calculSphere(R_terre/1000)
   plot!(plt, getindex.(surf, 1), getindex.(surf, 2), getindex.(surf, 3),
              seriestype = :surface, label = "", colorbar = false, 
              color = :gainsboro, opacity = 0.5, kwargs...)
end

"""
Dessiner sur le plot `plt` la trajectoire du satellite de propagateur `prop` dans le
référentiel ECI. Si `longueurFleche` est différent de `nothing`, dessiner une flèche à la
position initiale du satellite indiquant le sens de parcours de la trajectoire.

# Arguments
- `ts` : instants (en secondes écoulées depuis l'instant initial de `prop`) pour lequels on
   calcule la position du satellite. 
- `couleur` : couleur de la trajectoire
- `label` : étiquette de la trajectoire. 
"""
function dessinerTrajECI!(plt::Plots.Plot, prop::OrbitPropagator, ts::AbstractVector{<:Real},
      couleur::Colorant, label::Union{String, Int}, 
      longueurFleche::Union{Nothing, Real})

   trajECI, _ = propagate!(prop, ts)
   trajECI_km = trajECI/1000

   # Dessiner la trajectoire. 
   plot!(plt, getindex.(trajECI_km, 1), getindex.(trajECI_km, 2), 
         getindex.(trajECI_km, 3); seriestype = :path3d, seriescolor = couleur,
         label = label)

   # Dessiner la flèche si requis
   if longueurFleche != nothing 
      pos1 = trajECI_km[1]
      pos2 = trajECI_km[2]
      Δpos = pos2 - pos1
      surfCote, surfBase = calculCone(pos1, pos1 + longueurFleche*Δpos/norm(Δpos), 
                                      longueurFleche/3)
      plot!(plt, getindex.(surfCote, 1), getindex.(surfCote, 2), getindex.(surfCote, 3), 
            seriestype = :surface, color = couleur, opacity = 1, colorbar = false)
      plot!(plt, getindex.(surfBase, 1), getindex.(surfBase, 2), getindex.(surfBase, 3), 
            seriestype = :surface, color = couleur, opacity = 1, colorbar = false)
   end
end

"""
Dessiner sur le plot `plt` la trajectoire du satellite de propagateur `prop` dans le
référentiel ECEF. Si `angleOuverture` est différent de `nothing`, dessiner aussi les zones
terrestres couvertes par le satellite aux instants de `ts` (en secondes écoulées puis
`get_epoch(prop)`. Si `longueurFleche` est différent de `nothing`, dessiner une flèche à la
position initiale du satellite indiquant le sens de parcours de la trajectoire.

# Arguments
- `ts` : instants (en secondes écoulées depuis l'instant initial de `prop`) pour lequels on
   calcule la position du satellite. 
- `mats_ECI_ECEF` : Vecteur de matrice de rotation. La multiplication par mats_ECI_ECEF[i] 
   permet de passer des coordonnées ECI aux coordonnées ECEF à l'instant ts[i]. 
- `couleur` : couleur de la trajectoire
- `label` : étiquette de la trajectoire. 
"""
function dessinerTrajECEF!(plt::Plots.Plot, prop::OrbitPropagator, 
      angleOuverture::Union{Real, Nothing}, ts::AbstractVector{<:Real},
      mats_ECI_ECEF::AbstractVector{<:AbstractMatrix{<:Real}}, couleur::Colorant,
      label::Union{String, Int}, longueurFleche::Union{Real, Nothing})

   trajECI, _ = propagate!(prop, ts)
   trajECEF = @. (mats_ECI_ECEF * trajECI)
   trajECEF_km = trajECEF/1000

   # Dessiner la trajectoire. 
   plot!(plt, getindex.(trajECEF_km, 1), getindex.(trajECEF_km, 2), 
         getindex.(trajECEF_km, 3); seriestype = :path3d, seriescolor = couleur,
         label = label)

   # Dessiner la surface terrestre couverte 
   if angleOuverture != nothing
      a = get_a(get_mean_elements(prop))

      θ_max = calcul_θ(a, angleOuverture)

      trajGeo = ecef_to_geodetic.(trajECEF)
      ϕs = getindex.(trajGeo, 1)
      λs = getindex.(trajGeo, 2)

      surfs = calculDisqueGeomSpherique.(R_terre/1000, ϕs, λs, θ_max) 
      for surf in surfs
         plot!(plt, getindex.(surf, 1), getindex.(surf, 2), getindex.(surf, 3),
            seriestype = :surface, opacity = 0.5, seriescolor = couleur, colorbar = false, 
            label = label)
      end
   end

   # Dessiner la flèche si requis
   if longueurFleche != nothing 
      pos1 = trajECEF_km[1]
      pos2 = trajECEF_km[2]
      Δpos = pos2 - pos1
      surfCote, surfBase = calculCone(pos1, pos1 + longueurFleche*Δpos/norm(Δpos), 
                                      longueurFleche/3)
      plot!(plt, getindex.(surfCote, 1), getindex.(surfCote, 2), getindex.(surfCote, 3), 
            seriestype = :surface, color = couleur, opacity = 1, colorbar = false)
      plot!(plt, getindex.(surfBase, 1), getindex.(surfBase, 2), getindex.(surfBase, 3), 
            seriestype = :surface, color = couleur, opacity = 1, colorbar = false)
   end
end

"""
Dessiner sur `plt` les trajectoires des satellites de la constellation `conste` dans le
référentiel ECI. Si `longueurFleche` est différent de `nothing`, dessiner une flèche à la
position initiale des satellites indiquant le sens de parcours des trajectoires. Dessiner 
aussi la surface terreste. 

# Arguments 
- `ts` : instants (en secondes écoulées depuis l'instant initial de `conste`) pour lequels
   on calcule la position du satellite. 
"""
function dessinerConsteECI!(plt::Plots.Plot, conste::Conste, ts::AbstractVector{<:Real}; 
      longueurFleche::Union{Real, Nothing} = 700)

   # Propagateurs des satellites actifs 
   props = init_orbit_propagator.(Val(:twobody), 
      KeplerianElements.(conste.date, conste.as, 0, conste.is, conste.Ωs, 0, 
                         conste.Ms))
   
   plot!(plt, foreground_color_grid = :black, gridalpha = 0.3, size = (1000, 1000), xlabel =
         "x", ylabel = "y", zlabel = "z", title = " Trajectoires (ECI)")
   
   dessinerTerre!(plt)

   # Dessiner les trajectoire des satellites 
   # Vecteurs des couleurs pour les orbites
   cols = get(darkrainbow, conste.nbSats == 1 ? 0.5 : LinRange(0, 1, conste.nbSats))

   dessinerTrajECI!.(Ref(plt), props, Ref(ts), cols, collect(1:conste.nbSats), 
                     longueurFleche)
end

"""
Dessiner sur le plot `plt` les trajectoires des satellites de la constellation `conste`,
solution de l'instance `inst` dans le référentiel ECEF. Si `longueurFleche` est différent de
nothing, dessinez des flèches indiquant les sens de parcours des trajectoires. 
"""
function dessinerSolECEF!(plt::Plots.Plot, inst::DonneesSCDP, conste::Conste;
      longueurFleche::Union{Real, Nothing} = 700)

   @assert inst.instantInitialJulian ≈ conste.date

   # Propagateurs des satellites de la constellation
   props = init_orbit_propagator.(Val(:twobody), 
      KeplerianElements.(inst.instantInitialJulian, conste.as, 0, conste.is, conste.Ωs, 0, 
                         conste.Ms))
   
   plot!(plt, foreground_color_grid = :black, gridalpha = 0.3, size = (1000, 1000), xlabel =
         "x", ylabel = "y", zlabel = "z", title = " Trajectoires (ECEF)")

   dessinerTerre!(plt)

   # Dessiner les trajectoire des satellites 
   # Vecteurs des couleurs pour les orbites
   cols = get(darkrainbow, conste.nbSats == 1 ? 0.5 : LinRange(0, 1, conste.nbSats))

   dessinerTrajECEF!.(Ref(plt), props, nothing, Ref(inst.ts), Ref(inst.mats_ECI_ECEF),
                  cols, collect(1:conste.nbSats), longueurFleche)

   # Dessiner les cibles
   coordsCiblesECEF_km = inst.coordsCiblesECEF/1000
   label = inst.nbCibles ≤ 26 ? alphabet[1:inst.nbCibles][na, :] : nothing
   # Vecteurs des couleurs pour les cibles
   cols = get(darkrainbow, inst.nbCibles == 1 ? 0.5 : LinRange(0, 1, inst.nbCibles))
   scatter!(plt, getindex.(coordsCiblesECEF_km, 1)[na, :], 
            getindex.(coordsCiblesECEF_km, 2)[na, :], 
            getindex.(coordsCiblesECEF_km, 3)[na, :],
            markershape = :circle, markersize = 1, color = cols[na, :], 
            markerstrokecolor = cols[na, :], label = label)
end

"""
Fonction acotan à valeur dans ]0, π[. (ce qui n'est pas le cas de acot de Julia)
"""
acotan(y) = π/2 - atan(y)

""" 
Etant donné un point à la surface de la sphère de rayon `R`, de coords géographiques
lattitude, longitude `ϕ`, `λ` (en radian), renvoie le vecteur des coords cartésiennes de ce
point (dans la même unité que `R`). 
"""
function geog_to_carte(R::Number, ϕ::Number, λ::Number)
   return SA[R*cos(ϕ)*cos(λ), R*cos(ϕ)*sin(λ), R*sin(ϕ)]
end

"""
Renvoie la surface de la sphère de rayon `R` paramétrée en lattitude, longitude. Les points
sont exprimés en coordonnéees cartésiennes dans la même unité que `R`. 
"""
function calculSphere(R::Number, dens::Int = 3)
   lats = LinRange(-π/2, π/2, Int(floor(π*dens)))
   longs = LinRange(-π, π, Int(floor(2π*dens)))
   return [geog_to_carte(R, lat, long) for lat in lats, long in longs] 
end

"""
On travaille en géométrie sphérique sur la sphère de rayon `R` centrée en 0.  Renvoie la
surface correspondant à un disque de centre  M0 := (`ϕ_0`, `λ_0`) (lattitude, longitude en
radian) et rayon angulaire `θ_max` (en radian) paramétrée en coordonnées polaires et dont
les points sont exprimés en coordonnées cartésiennes dans le même unité que `R`. 
Les coordonnées polaires de centre M0 de la géométrie sphérique d'un point M sont noté (η,
θ). η est l'angle entre le méridien passant par M_0 et le grand cercle passant par M0 et M.
θ est la distance angulaire entre le centre du disque et M. 
"""
function calculDisqueGeomSpherique(R::Number, ϕ_0::Number, λ_0::Number, θ_max::Number, 
      nbPoints_η::Int = 8, nbPoints_θ::Int = 3)

   ηs = LinRange(-π, π, nbPoints_η)
   θs = LinRange(0, θ_max, nbPoints_θ)

   # Paramétrisation de la surface
   function paramSurf(η, θ)

      # Trigonométrie sphérique
      ϕ = asin(sin(ϕ_0)*cos(θ) + cos(ϕ_0)*sin(θ)*cos(η))
      Δλ = acotan((cot(θ)*cos(ϕ_0) - sin(ϕ_0)*cos(η))/sin(η)) 
      η < 0 && (Δλ = Δλ - π)

      λ = λ_0 + Δλ
      # Conserver une longitude dans ]-π, π]
      λ > π && (λ = λ - 2π)
      λ ≤ -π && (λ = λ + 2π) 

      return geog_to_carte(R, ϕ, λ)
   end

   return [paramSurf(η, θ) for η in ηs, θ in θs]
end

"""
Calculer la surface d'un cône. 

Si A et B sont deux points de R³, calculer la surface du cône d'axe de révolution (AB), de
sommet B, de plan de base passant par A (et parallèle à (AB)) et de rayon à la base r_0.
Renvoie un couple de matrice, la première correspond au paramétrage du coté du cône, les
éléments de la matrice étant les points de la surface exprimés en coordonnées cartésiennes
et la deuxième correspond à la base du cône. 

# Arguments
   `A` : Array de longueur 3 représentant un point de R³.
   `B` : Array de longueur 3 représentant un point de R³.
   `r_0` : Rayon à l'origine du cône.  
"""
function calculCone(A::AbstractVector{<:Real}, B::AbstractVector{<:Real}, r_0::Real, 
   nbPoints_h::Int = 5, nbPoints_η::Int = 8, nbPoints_s::Int = 3, nbPoints_α = 8)

   H = norm(B .- A)
   θ = acos((B[3] - A[3])/H)
   ϕ = atan(B[2] - A[2], B[1] - A[1])
   
   cos_θ = cos(θ)
   sin_θ = sin(θ)
   R_θ = [cos_θ  0 sin_θ
          0      1 0
          -sin_θ 0 cos_θ]

   cos_ϕ = cos(ϕ)
   sin_ϕ = sin(ϕ)
   S_ϕ = [cos_ϕ -sin_ϕ 0
          sin_ϕ  cos_ϕ 0
          0      0     1]
   
   # Calcul de la surface correspondant au côté du cône
   hs = LinRange(0, H, nbPoints_h)
   ηs = LinRange(0, 2π, nbPoints_η)
   paramSurfCote(h, η) = begin 
      r = r_0*(1 - h/H)
      return A .+ S_ϕ*R_θ*[r*cos(η), r*sin(η), h]
   end

   # Calcul de la surface correspondant à la base du cône
   ss = LinRange(0, r_0, nbPoints_s)
   αs = LinRange(0, 2π, nbPoints_α)
   paramSurfBase(s, α) =  A .+ S_ϕ*R_θ*[s*cos(α), s*sin(α), 0]

   return [paramSurfCote(h, η) for h in hs, η in ηs], 
         [paramSurfBase(s, α) for s in ss, α in αs]
end

