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
   pretty_table(matCoordsCibles, header = (["Latitude", "Longitude"], ["deg", "deg"]))#, 
   #             row_names = inst.nbCibles ≤ 26 ? alphabet[1:inst.nbCibles] : nothing)
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

"""
Position pixel en fonction du temps (inspiré de page Wiki https://fr.wikipedia.org/wiki/Coordonn%C3%A9es_sph%C3%A9riques)
"""
function TargetInTime(NUM_PIXEL,x,y,z)
	global x_t = fill(0.0, NUM_PIXEL, NUM_TIME+1)
	global y_t = fill(0.0, NUM_PIXEL, NUM_TIME+1)
	global z_t = fill(0.0, NUM_PIXEL, NUM_TIME+1)
	for p in 1:NUM_TIME
	# vit:(Rad/s * tps:s) #https://fr.wikipedia.org/wiki/Temps_sid%C3%A9ral
		local dist_ang_pix = (2*pi/T_terre)*(p-1)*dt
		for j in 1:NUM_PIXEL
			atan_longitude_p = atan(y[j],x[j]) + dist_ang_pix
		# Radian : latitude (inclinaison par rapport au plan equatoriale)
			local LatPix = pi/2 - acos(z[j]/R_terre) 
			x_t[j,p] = R_terre * cos(atan_longitude_p) * cos(LatPix)
			y_t[j,p] = R_terre * sin(atan_longitude_p) * cos(LatPix)
			z_t[j,p] = z[j]
		end
	end
	
	return x_t,y_t,z_t
end

"""
Position pixel en fonction du temps (inspiré de page Wiki https://fr.wikipedia.org/wiki/Coordonn%C3%A9es_sph%C3%A9riques)
"""
function TargetInTimeECI(NUM_PIXEL,x,y,z)
	global x_t = fill(0.0, NUM_PIXEL, NUM_TIME )
	global y_t = fill(0.0, NUM_PIXEL, NUM_TIME)
	global z_t = fill(0.0, NUM_PIXEL, NUM_TIME)
	for p in 1:NUM_TIME
	# vit:(Rad/s * tps:s) #https://fr.wikipedia.org/wiki/Temps_sid%C3%A9ral
		local dist_ang_pix = (2*pi/T_terre)*(p-1)*dt
		global PrimoGiorno = DatetoJD(1961,4,12,0,0,0) #Date initiale de propagation
		local matrix_rotation = r_ecef_to_eci(PEF(), J2000(), PrimoGiorno+(dt*p))
		
		for j in 1:NUM_PIXEL
			atan_longitude_p = atan(y[j],x[j]) + dist_ang_pix
		# Radian : latitude (inclinaison par rapport au plan equatoriale)
			local LatPix = pi/2 - acos(z[j]/R_terre) 
			x_temp = R_terre * cos(atan_longitude_p) * cos(LatPix)
			y_temp = R_terre * sin(atan_longitude_p) * cos(LatPix)
			z_temp = z[j]
			local ECEFtoECI = matrix_rotation*([x_temp,y_temp,z_temp])
			x_t[j,p] = ECEFtoECI[1]
			y_t[j,p] = ECEFtoECI[2]
			z_t[j,p] = ECEFtoECI[3]
		end
	end
	
	return x_t,y_t,z_t
end

"""
"""
function TargetSatPlot(numbOfDiscretize, j, k, l, s)
	# Discretization Keplerian parameters
	inclinaisonSet = Discretization(numbOfDiscretize,pi)
	noeudAscendantSet = Discretization(numbOfDiscretize,2*pi)
	meanAnomalySet = Discretization(numbOfDiscretize,2*pi)

	x_projSat1 = fill(0.0,NUM_TIME)
	y_projSat1 = fill(0.0,NUM_TIME)
	z_projSat1 = fill(0.0,NUM_TIME)

	for p=1:NUM_TIME

		x_projSat1[p] = (((cos(meanAnomalySet[k])*cos(sqrt(μ/((R_terre+altitudeSet[j])^3))*(p-1)*dt) - 
		sin(meanAnomalySet[k])*sin(sqrt(μ/((R_terre+altitudeSet[j])^3))*(p-1)*dt)) * cos(noeudAscendantSet[l])*(R_terre+altitudeSet[j])) - 
		((sin(meanAnomalySet[k])*cos(sqrt(μ/((R_terre+altitudeSet[j])^3))*(p-1)*dt) + cos(meanAnomalySet[k])*sin(sqrt(μ/((R_terre+altitudeSet[j])^3))*(p-1)*dt))
		* cos(inclinaisonSet[s])*sin(noeudAscendantSet[l])*(R_terre+altitudeSet[j]))) 

		y_projSat1[p] = (((cos(meanAnomalySet[k])*cos(sqrt(μ/((R_terre+altitudeSet[j])^3))*(p-1)*dt) - 
		sin(meanAnomalySet[k])*sin(sqrt(μ/((R_terre+altitudeSet[j])^3))*(p-1)*dt)) * 
		sin(noeudAscendantSet[l])*(R_terre+altitudeSet[j])) + ((sin(meanAnomalySet[k])*cos(sqrt(μ/((R_terre+altitudeSet[j])^3))*(p-1)*dt) + 
		cos(meanAnomalySet[k])*sin(sqrt(μ/((R_terre+altitudeSet[j])^3))*(p-1)*dt)) * cos(noeudAscendantSet[l])*cos(inclinaisonSet[s])*(R_terre+altitudeSet[j]))) 
			
		z_projSat1[p] = ((sin(meanAnomalySet[k])*cos(sqrt(μ/((R_terre+altitudeSet[j])^3))*(p-1)*dt) + 
		cos(meanAnomalySet[k])*sin(sqrt(μ/((R_terre+altitudeSet[j])^3))*(p-1)*dt)) * sin(inclinaisonSet[s])*(R_terre+altitudeSet[j])) 

	end
	return x_projSat1,y_projSat1,z_projSat1
end

"""
Temps de passage
"""
function CountRevisitTime(TimeSeen_track,inst)
	for index_ville in 1:size(TimeSeen_track)[1]
		figure()
		PyPlot.title("Target :"*string(index_ville))
		PyPlot.xlabel("Time SIMULATION in seconds")
		PyPlot.ylabel("Visibility")
		PyPlot.yticks([])
		PyPlot.plot(collect(0:dt:T_terre),TimeSeen_track[index_ville,:])
		DELTA_T = T_terre/inst.n_0
		println("number revisit time ",string(sum(TimeSeen_track[index_ville,:])),"; number minimal revisit time ",string(T_terre/DELTA_T))
		Time_intervall_revisit = collect(0:convert(Int64,(DELTA_T)/dt):convert(Int64,T_terre/dt))
		for time_inter in Time_intervall_revisit
			interval_time = 0*collect(0:dt:T_terre)
			interval_time[time_inter+1] = 1.25
			PyPlot.plot(collect(0:dt:T_terre),interval_time, color = "red")
		end
	end
end

"""
"""
function CreateSfera(R_terre_sfera,x_init=0.0,y_init=0.0,z_init=0.0,N=32,colorSfera="gray",Axes="True")
	# Create a sfera
	phi = range(0, stop=2*pi, length=N)
	theta = range(0, stop=pi, length=N)
	x = R_terre_sfera*cos.(phi) .* sin.(theta)'
	y = R_terre_sfera*sin.(phi) .* sin.(theta)'
	z = R_terre_sfera*repeat(cos.(theta)',outer=[N, 1])
	plot_wireframe(x_init.+x,y_init.+y,z_init.+z, color=colorSfera,alpha=0.1)
	if Axes == "True"
		zerozero = zeros(1000)
		lunghezza_asse = range(-12*1e6,stop = 12*1e6,length = 1000)
		PyPlot.plot(zerozero, zerozero, lunghezza_asse, color = "black",alpha=0.2)
		PyPlot.plot(lunghezza_asse,zerozero, zerozero, color = "black",alpha=0.1)
		PyPlot.plot(zerozero,lunghezza_asse, zerozero, color = "black",alpha=0.1)
	end
	PyPlot.xlabel("X (m)")
	PyPlot.ylabel("Y (m)")
	PyPlot.zlabel("Z (m)")
end

"""
Create Piano equatoriale
"""
function CreateEquatorialPlane(N=8000000)
	dist = range(-N, stop=N, step=100000)
	len_discret = length(dist)
	x_eq = zeros(len_discret,len_discret)
	z_eq = zeros(len_discret,len_discret)
	global count = 1
	@views for row in dist
		global count
		b = x_eq[count, :]
		b[:] .= row
		count= count+1
	end
	y_eq = x_eq'
	PyPlot.plot_surface(x_eq,y_eq,z_eq, color="yellow",alpha=0.2 )
end

"""
"""
function plot_target(x_target,y_target,z_target)
	CreateEquatorialPlane()
	CreateSfera()
	
	PyPlot.title("NUMBER of TARGET : "*string(length(x_target)))
	PyPlot.scatter3D(x_target,y_target,z_target, color = "blue")#, marker = "o", markersize = 5)	
end

"""
Plot orbitals
"""
function OrbitalPlanes(R_terre_SAT,x_init,y_init,z_init,x_sat,y_sat,z_sat,N = 10)

	for i in 1:length(x_init)
		# SPhère couverture sat
		#CreateSfera(R_terre_SAT,x_init[i],y_init[i],z_init[i],N,"orange","False")		
	end
	# Position sat en ECEF
	PyPlot.scatter3D(x_sat,y_sat,z_sat,color = "violet", label = "Optimal Satellite_ECEF")
	# projection sat sur surface en ECEF
	PyPlot.scatter3D(x_init,y_init,z_init,color = "green", label = "PROJECTION Satellite_ECEF")
end

"""
Output a vector binary with 1 if the target is seen at time t
"""
function Validation(altitudeSet,InclinaisonSet,RAANSet,MeanAnomalySet,IndexMeanAnom,IndexRaan,Indexinclination,IndexAltitude,inst)

   NUM_PIXEL=inst.nbCibles
   global TimeSeen_track = zeros(NUM_PIXEL,length(collect(0:dt:T_terre)))   

   for i in 1:length(IndexAltitude)	
      
      local LatitudeSat,LongitudeSat = ProjSatPlotPOLARDominique(R_terre+altitudeSet[IndexAltitude[i]], 
                                                   InclinaisonSet[Indexinclination[i]], 
                                                   RAANSet[IndexRaan[i]], 
                                                   MeanAnomalySet[IndexMeanAnom[i]])
      for j in 1:length(LatitudeSat)
         for indexTarget in 1:(inst.nbCibles)    
            Latitude = inst.coordsCibles[indexTarget][1]
            Longitude = inst.coordsCibles[indexTarget][2]
            if 	((Theta[IndexAltitude[i]] >= abs(LatitudeSat[j]-Latitude)) && (Theta[IndexAltitude[i]] >= abs(LongitudeSat[j]-Longitude)*cos(Latitude)))			
               TimeSeen_track[indexTarget,j] = 1.0
            end
         end
      end
   end	

   return TimeSeen_track
end

"""
Possibility to show orbital plan in 3D with a sphere (Earth) using Keplerian's parameters as inputs
"""
function GraphicalValidation(altitudeSet,InclinaisonSet,RAANSet,MeanAnomalySet,IndexMeanAnom,IndexRaan,Indexinclination,IndexAltitude,inst)

   CreateSfera(R_terre)
   CreateEquatorialPlane()

   NUM_PIXEL=inst.nbCibles
   global TimeSeen_track = zeros(NUM_PIXEL,length(collect(0:dt:T_terre)))   

   for i in 1:NUM_PIXEL
      x_temp,y_temp,z_temp = PolarToCartesian(R_terre,pi/2-inst.coordsCibles[i][1],inst.coordsCibles[i][2])
      PyPlot.scatter3D(x_temp,y_temp,z_temp, label = "Target "*string(i))
   end

   for i in 1:length(IndexAltitude)	
      
      local LatitudeSat,LongitudeSat = ProjSatPlotPOLARDominique(R_terre+altitudeSet[IndexAltitude[i]], 
                                                   InclinaisonSet[Indexinclination[i]], 
                                                   RAANSet[IndexRaan[i]], 
                                                   MeanAnomalySet[IndexMeanAnom[i]])
      println("Theta :"*string(Theta[IndexAltitude[i]]))
      println("a,i,nA,mA :"*string([R_terre+altitudeSet[IndexAltitude[i]],InclinaisonSet[Indexinclination[i]],RAANSet[IndexRaan[i]],MeanAnomalySet[IndexMeanAnom[i]]]))

      local orbp = init_orbit_propagator(Modele, time_zero_simulation, R_terre+altitudeSet[IndexAltitude[i]], 0.0, InclinaisonSet[Indexinclination[i]], RAANSet[IndexRaan[i]],0.0, 
                                 MeanAnomalySet[IndexMeanAnom[i]])
      local r,v = propagate!(orbp, collect(0:dt:T_terre)) 

      local x_sat=fill(0.0, length(r)) 
      local y_sat=fill(0.0, length(r)) 
      local z_sat=fill(0.0, length(r))
      
      for j in 1:(length(r))
         local matrix_rotation = r_eci_to_ecef(TOD(), PEF(), time_zero_simulation+(dt*(j-1))/T_terre)
         local ECItoECEF = matrix_rotation*(r[j][:])
         
         local x_sat[j]=ECItoECEF[1]
         local y_sat[j]=ECItoECEF[2]
         local z_sat[j]=ECItoECEF[3]

         for indexTarget in 1:(inst.nbCibles)    
            Latitude = inst.coordsCibles[indexTarget][1]
            Longitude = inst.coordsCibles[indexTarget][2]
            if 	((Theta[IndexAltitude[i]] >= abs(LatitudeSat[j]-Latitude)) && (Theta[IndexAltitude[i]] >= abs(LongitudeSat[j]-Longitude)*cos(Latitude)))			
               TimeSeen_track[indexTarget,j] = 1.0
            end
         end
      end

      # Position sat en ECEF
      PyPlot.plot(x_sat,y_sat,z_sat,color = "magenta")#, label = legend_plot[end-1])

   end	

   legend()	
   CountRevisitTime(TimeSeen_track,inst)
end