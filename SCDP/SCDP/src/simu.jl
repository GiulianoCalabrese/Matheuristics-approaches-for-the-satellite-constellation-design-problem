"""
Déterminer la couverture de cibles terrestres que réalise une constellation de
satellites. 
"""

""" 
Calculer θ (rad) le demi-angle de la zone couverte sur Terre par le satellite rayon d'orbite a.
(voir définition dans le rapport).

# Arguments
- `a` : rayon d'orbite du satellite (m).  
- `α` : demi-angle d'ouvertue du satellite (rad). 
"""
function calcul_θ(a::Real, α::Real)
   sin_α = sin(α)
   return sin_α < (R_terre/a) ? asin(a*sin_α/R_terre) - α : π/2 - asin(R_terre/a)
end

"""
Calculer α (rad) le demi-angle d'ouverture d'un satellite de rayon d'orbite `a` (m) à partir
de `θ` (rad) le demi-angle de la zone terrestre couverte par ce satellite. A `a` fixé, cette
fonction est la réciproque de la précédente `calcul_θ`.
"""
calcul_α(a::Real, θ::Real) = atan(R_terre*sin(θ)/(a-R_terre*cos(θ))) 

"""
Obtenir coordonnées polaires des cartésiennes

# Arguments
- `x,y,z` : coordonnées cartésiennes 
"""
function CartesianToPolar(x,y,z)
	# position of projection of satellite on Earth (ground track)
	ρ = sqrt(x^2+y^2+z^2)
	θ = acos(z/ρ)
	ϕ = atan(y,x)
	return ρ,θ,ϕ
end

"""
Renvoie  true ssi la cible de coordonnées `coordsCible` est survolée par le satellite de 
coordonnées `coordsSat` couvrant la surface terrestre de demi-angle `θ`.

# Arguments
- `coordsSat` : coordonnées cartésiennes du satellite en m dans le référentiel ECEF. 
- `coordscible` : coordonnées cartésiennes de la cible en m dans le référentiel ECEF. 
- `θ` : demi-angle de la surface couverte par le satellite (rad).
"""
function detSiSurvol(coordsSat::AbstractVector{<:Real}, coordsCible::AbstractVector{<:Real}, 
      θ::Real)
	  
    # cos_η = dot(coordsSat/norm(coordsSat), coordsCible/norm(coordsCible))
    # return cos_η ≥ cos(θ)
   
	ρ_temp,colat_temp,long_temp = CartesianToPolar(coordsSat[1],coordsSat[2],coordsSat[3])
	LatSat = pi/2-colat_temp
	LongSat = long_temp+pi
   
    lat_η = (abs.(LatSat - coordsCible[1]) <= θ)
    # long_η = (abs.(LongSat - coordsCible[2]) <= θ)
 
    long_η = ((abs.(LongSat - coordsCible[2]).*cos.(coordsCible[1])) <= θ)
    
	return (lat_η && long_η) ≥ 1


end

"""
Etant donné une instance du SCDP `inst` et une constellation de satellites, renvoie le
nombre d'intervalles de surveillance sur lesquelles les cibles ne sont survolées par aucun
satellite de la constellation. 

Plus précisément, le résultat renvoyé s'obtient en calculant, pour chaque cible, le nombre 
d'intervalles de surveillance sur lesquelle cette cible n'est survolée par aucun satellite 
puis en sommant ces nombres.

# Arguments
- `inst` : Une instance du SCDP. 
- `as` : rayons d'orbite des satellites de la constellation (m). 
- `is` : inclinaisons des satellites de la constellation (rad). 
- `Ωs` : longitudes des noeuds ascendants des satellites de la constellation (rad). 
- `Ms` : anomalies moyennes initiales des satellites de la constellation (rad). 
"""
function calculNbIntervsSansSurvol(inst::DonneesSCDP, as::AbstractVector{<:Real}, 
      is::AbstractVector{<:Real}, Ωs::AbstractVector{<:Real}, Ms::AbstractVector{<:Real})

   nbIntervsSansSurvol = inst.n_0*inst.nbCibles
   
   # Si la constellation ne possède aucun satellite, la fonction peut s'arrêter là.
   as == [] && return nbIntervsSansSurvol

   # Propagateurs des satellites actifs 
   props = init_orbit_propagator.(Val(:twobody), 
      KeplerianElements.(inst.instantInitialJulian, as, 0, is, Ωs, 0, Ms))

   θs = calcul_θ.(as, inst.angleOuverture)
   
   # Si ind_cible est un indice de cible, boolsCiblesSurvolees[ind_cible] vaut true ssi la
   # cible d'indice ind_cible a été survolée par un des satellites actifs entre l'instant de
   # début de l'intervalle de surveillance courant et l'instant courant. 
   boolsCiblesSurvolees =  Vector{Bool}(undef, inst.nbCibles)

   # Parcours des intervalles de surveillance
   for ind_interv in 1:inst.n_0

      boolsCiblesSurvolees .= false
      # Indice de temps depuis le début de l'intervalle courant.
      ind_t_loc = 1 
      
      # Parcours des instants de l'intervalle de surveillance n°ind_interv.
      while ind_t_loc ≤ inst.nbInstants && !all(boolsCiblesSurvolees) 

         ind_t = inst.nbInstants*(ind_interv-1) + ind_t_loc
         # Coordonnées des satellites (actifs) à l'instant actuel.
         coordsSatsECEF = [inst.mats_ECI_ECEF[ind_t]*propagate!(prop, inst.ts[ind_t])[1] 
                           for prop in props]
         
         # Parmis les cibles non survolée depuis l'instant de début de l'intervalle de
         # surveillance courant, chercher s'il en existe qui sont survolées à l'instant
         # actuel.
         nouvellesCiblesSurvolees = dropdims(any(
               # detSiSurvol.(coordsSatsECEF[na, :], 
                            # inst.coordsCiblesECEF[.!boolsCiblesSurvolees, na], 
                            # θs[na, :]), 
							
			   detSiSurvol.(coordsSatsECEF[na, :], 
                            inst.coordsCibles[.!boolsCiblesSurvolees, na], 
                            θs[na, :]), 	
							
               dims = 2), dims = 2)
         
         # Mettre à jour les cibles survolée depuis le début de l'intervalle courant. 
         boolsCiblesSurvolees[.!boolsCiblesSurvolees] = nouvellesCiblesSurvolees 
         
         nbIntervsSansSurvol -= sum(nouvellesCiblesSurvolees)
         ind_t_loc += 1 
      end
   end

   return nbIntervsSansSurvol
end

"""
Sur la période de durée `ts[end]` secondes depuis l'instant initial de `props[1]`, calculer
les intervalles de temps sur lesquelles la cible de coordonnées `coordCibleECEF` est
survolée par au moins un satellite dont le propagateur est dans `props`. Renvoie un array
d'arrays statiques, chacun des ces arrays statiques étant un intervalle sur lequel la cible
est survolée (les bornes des intervalles sont en secondes écoulées depuis `get_epoch(prop)`).
Il est supposé que les propagateurs dans `props` ont même instant initial. 

# Arguments 
- `coordCibleECEF` : vecteur contenant les coordonnées de la cible (en m) dans le
   référentiel ECEF PEF. 
- `props` : vecteur des propagateurs des satellites.
- `ts` : vecteur de temps sur lequel l'analyse est conduite, en seconde écoulées depuis
   l'instant initial de prop.
- `θs` : Demi-angles des surfaces terrestres couvertes par les satellites de `props`. 
- `δ` : définit la précision selon laquelle les bornes des intervalles de temps (cad les 
   instants où la cible passe de l'état survolé à l'état non survolé et vice-versa) sont 
   calculés. 
- `nbMaxIter` : nombre max d'itérations réalisées lorsque l'on détermine les instants où un
   une cible passe de l'état survolée à l'état non-survolée et vice-versa. 

Il s'agit d'une adaptation de la fonction `ground_station_accesses` de SatelliteToolbox.
"""
function calculIntervallesSurvol(
      coordCibleECEF::AbstractVector{<:Real}, 
      props::AbstractVector{<:OrbitPropagator},
      ts::AbstractVector{<:Real},
      θs::AbstractVector{<:Real}, 
      δ::Real,
      nbMaxIter::Int)

   intervallesSurvol = SVector{2, Float64}[]

   instantInitial = get_epoch(props[1])

   # Renvoie true si la cible est survolée par un satellite de props à l'instant t (en
   # seconde depuis get_epoch(props[1])), false sinon. 
   calculCibleSurvolee(t) = begin 
      mat_ECI_ECEF = r_eci_to_ecef(TOD(), PEF(), instantInitial + t/T_terre)
      coordsSatsECEF_t  = [mat_ECI_ECEF*propagate!(prop, t)[1] for prop in props]
      return any(detSiSurvol.(coordsSatsECEF_t, Ref(coordCibleECEF), θs))
   end

   # Initialisation avec l'instant initial
   instantDebutSurvol = 0

   # Initialisation
   etat = calculCibleSurvolee(0) ? :survolee : :nonSurvolee

   for ind_t in 2:length(ts)
      boolCibleSurvolee = calculCibleSurvolee(ts[ind_t])
      if etat == :survolee && !boolCibleSurvolee
         instantFinSurvol = SatelliteToolbox.find_crossing(calculCibleSurvolee, 
                              ts[ind_t-1], ts[ind_t], true, false; Δ = δ, max = nbMaxIter)
         etat = :nonSurvolee
         push!(intervallesSurvol, SA[instantDebutSurvol, instantFinSurvol])
      elseif etat == :nonSurvolee && boolCibleSurvolee
         instantDebutSurvol = SatelliteToolbox.find_crossing(calculCibleSurvolee, 
                              ts[ind_t-1], ts[ind_t], false, true; Δ = δ, max = nbMaxIter)
         etat = :survolee
      end
   end
   
   # Si l'analyse s'achève alors que la cible est survolée, l'instant de fin de
   # l'analyse est considéré comme l'instant de fin de survol.
   if etat == :survolee
      push!(intervallesSurvol, SA[instantDebutSurvol, ts[end]])
   end

   return intervallesSurvol
end

"""
Calculer la matrice des intervalles où les satellites de la constellation `conste` survolent
la cible de coordonnées `inst.coordsCibles[ind_cible]`.

Si `conste` contient `nbSat` satellites, la matrice renvoyée possède 2(nbSat+1) colonnes.
Pour s = 1, ..., nbSat, la (2s-1)ième colonne de la matrice renvoyée contient les bornes
inférieures des intervalles où le sième satellite voit la cible et la (2s)ième colonne les
bornes supérieures de ces intervalles. 
La (2nbSat+1)ième colonne contient les bornes inférieures des intervalles où la cible est 
survolée par au moins un satellite de la constellation et la (2(nbSat+1))ième colonnes les
bornes supérieures de ces intervalles. 
Afin que les colonnes aient même longueur, on ajoute aux colonnes des chaines vides jusqu'à 
ce qu'elles aient la longueur de la plus longue colonne. 

# Arguments
- `inst` : une instance du SCDP.
- `conste` : constellation solution de l'instance `inst`.
- `ind_cible` : indice d'une cible de l'instance `inst`.
- `format` : le format des bornes des intervalles. Si `format = :sec`, les bornes des
   intervalles sont en secondes écoulées depuis `inst.instantInitial`. Si `format = :DT`,
   les bornes sont au format `DateTime`.
- `δ` : définit la précision selon laquelle les bornes des intervalles de temps (cad les 
   instants où la cible passe de l'état survolé à l'état non survolé et vice-versa) sont 
   calculés. Ces bornes sont connues avec une précision inférieure à `δ` secondes (sauf si 
   le nombre maximal d'itérations est atteint). 
- `nbMaxIter` : nombre max d'itérations réalisées lorsque l'on détermine les instants où un
   une cible passe de l'état survolée à l'état non-survolée et vice-versa. 
"""
function matIntervallesSurvol(inst::DonneesSCDP, conste::Conste, 
      ind_cible::Int; format::Symbol = :DT, δ::Float64 = 1e-3, nbMaxIter::Int=100)
   
   @assert inst.instantInitialJulian ≈ conste.date 

   # Propagateurs des satellites de la constellation. 
   props = init_orbit_propagator.(Val(:twobody), 
      KeplerianElements.(inst.instantInitialJulian, conste.as, 0, conste.is, conste.Ωs, 0, 
                         conste.Ms))

   θs = calcul_θ.(conste.as, inst.angleOuverture)

   # Vecteur contenant les matrices des intervalles de survol pour chaque satellite.
   intervsSurvolSats = map(1:conste.nbSats) do ind_sat 
      intervs = calculIntervallesSurvol(inst.coordsCiblesECEF[ind_cible],
                                       props[ind_sat:ind_sat], 
                                       inst.ts, 
                                       θs[ind_sat:ind_sat], 
                                       δ, 
                                       nbMaxIter)
      hcat(getindex.(intervs, 1), getindex.(intervs, 2))
   end

   # Ajout des intervalles où la cible est survolée par au moins un satellite que l'on
   # convertit en min.
   unionIntervsSurvol = calculIntervallesSurvol(inst.coordsCiblesECEF[ind_cible], 
                                                props, 
                                                inst.ts,
                                                θs, 
                                                δ, 
                                                nbMaxIter)
   intervsSurvol = [intervsSurvolSats..., 
                       hcat(getindex.(unionIntervsSurvol, 1), 
                            getindex.(unionIntervsSurvol, 2))]

   # On convertit les bornes en DateTime si requis
   format == :DT && (
      intervsSurvol = 
      [inst.instantInitial .+ Dates.Second.(floor.(intervs)) for intervs in intervsSurvol]
     )

   # Complétion de chaque matrice afin d'obtenir des matrices qui ont même nombre de lignes.
   nbMaxInterv = maximum(size.(intervsSurvol, 1))
   intervsSurvol = map(intervsSurvol) do intervs 
      nbInterv = size(intervs, 1)
      vcat(intervs, fill("", nbMaxInterv - nbInterv, 2))
   end

   colCibles = [inst.nbCibles ≤ 26 ? alphabet[ind_cible] : "*", fill(".", nbMaxInterv - 1)...]
   return hcat(colCibles, intervsSurvol...)
end

"""
Renvoie l'angle en radians dans ]-π, π] de la rotation entre le référentiel ECI et le
référentiel ECEF au jour julien date. 
"""
function calcul_angle_ECI_ECEF(date::Real)
   mat_ECI_ECEF = r_eci_to_ecef(TOD(), PEF(), date)
   return atan(-mat_ECI_ECEF[2, 1], mat_ECI_ECEF[1, 1])
end

"""
Renvoie la street of coverage du satellite en fonction de son plan orbitale.
INPUT :  4 paramètres képleriens de 1 satellite
OUTPUT : street of coverage en coordonnée polaires dans le temps
"""
function ProjSatPlotPOLARDominique(Altezza, inclinaison, noeudAscendant, meanAnomaly)

   LatitudeSat = fill(0.0,NUM_TIME+1)
   LongitudeSat = fill(0.0,NUM_TIME+1)
   t_p = (sqrt(μ/(Altezza^3)))
   t_u = sqrt(μ/(Altezza))
   t_GM = sqrt(Altezza/μ)

   for p in 1:(NUM_TIME+1)
      
      LatitudeSat[p] = asin(round(((sin(inclinaison)*Altezza*sin(meanAnomaly))*cos(t_p*((p-1)*dt))/Altezza) 
      + ((sin(inclinaison)*t_u*cos(meanAnomaly))*sin(t_p*((p-1)*dt))*t_GM), digits=8))

      LongitudeSat[p] = (-(calcul_angle_ECI_ECEF(time_zero_simulation) + (we*((p-1)*dt))) + 
      atan((((sin(noeudAscendant)*Altezza*cos(meanAnomaly)) + (cos(noeudAscendant)*cos(inclinaison)*
      Altezza*sin(meanAnomaly)))*cos(t_p*((p-1)*dt))/Altezza)	+ ((-(sin(noeudAscendant)*t_u*sin(meanAnomaly))
      + (cos(noeudAscendant)*cos(inclinaison)*t_u*cos(meanAnomaly)))*sin(t_p*((p-1)*dt))*t_GM),
      ((((cos(noeudAscendant)*Altezza*cos(meanAnomaly)) - (sin(noeudAscendant)*cos(inclinaison)*
      Altezza*sin(meanAnomaly)))* cos(t_p*((p-1)*dt))/Altezza) + ((-(cos(noeudAscendant)*t_u*sin(meanAnomaly))
      -(sin(noeudAscendant)*cos(inclinaison)*t_u*cos(meanAnomaly)))*sin((t_p*((p-1)*dt)))*t_GM))))%(2*pi)
      
      if LongitudeSat[p] <= 0 
         LongitudeSat[p] += 2*pi
      end
   end

   return LatitudeSat,LongitudeSat
end

"""
Fonction qui permet la discretization d'un interval continu à discrétiser avec borne inf et sup en fonction d'un nombre de discrétisation à donner.
INPUT : nombre d'éléements dans l'ensemble discret, borne inf et sup de l'interval 
OUTPUT : ensemble discret de l'interval continu
"""
function Discretization(numbOfDiscretize,lowerBound,upperBound)
   angle_discret = LinRange(lowerBound,upperBound,numbOfDiscretize)
   return angle_discret
end

"""
Obtenir coordonnées cartésiennes des polaires
"""
function PolarToCartesian(ρ,θ,ϕ)
	# position of projection of satellite on Earth (ground track)
	x = ρ*cos(ϕ)*sin(θ)
	y = ρ*sin(ϕ)*sin(θ)
	z = ρ*cos(θ)
	return x,y,z
end

"""
Coordonnées cartésiennes en 3D des points à surveiller périodiquement
"""
function generateTargetsCartesians(NUM_PIXEL)
	
	local x = fill(0.0, NUM_PIXEL)
	local y = fill(0.0, NUM_PIXEL)
	local z = fill(0.0, NUM_PIXEL)
    
	# https://mathworld.wolfram.com/SpherePointPicking.html
	for i in 1:NUM_PIXEL
		phi  = 2. * pi * rand()
		csth = 1.0 - 2.0 * rand()
		snth = sqrt(1.0 - csth*csth)
		x[i]=R_terre * snth * cos(phi)
		y[i]=R_terre * snth * sin(phi)
		z[i]=R_terre * csth
	end	
	return x,y,z
end

"""
Get hyper matrix with all the distance between sat anche target at each time step for eache possibl orbital
"""
function getHyperMatrix()

	NbrTOTALsat,N_inclinaison,N_RAAN,N_meanAnomaly,InclinaisonSet,RAANSet,MeanAnomalySet = getCardinalKeplerianParam()

	if !isfile("HyperMatrixLat_"*string(dt)*"sec_"*string(NUM_PIXEL)*"Targets.jld2")
		CoverageSatLat = fill(0.0,NbrTOTALsat,NUM_TIME,NUM_PIXEL)
		CoverageSatLong = fill(0.0,NbrTOTALsat,NUM_TIME,NUM_PIXEL)
		Thetamax = zeros(NbrTOTALsat,NUM_TIME,NUM_PIXEL)
		for j in 0:(altitudeLength-1)
			for s in 0:(N_inclinaison-1)
				for l in 0:(N_RAAN-1)
					for k in 1:N_meanAnomaly
						global lat_Sat,long_Sat = ProjSatPlotPOLARDominique(R_terre+altitudeSet[j+1], InclinaisonSet[s+1],
						RAANSet[l+1], MeanAnomalySet[k])
						for p in 1:NUM_TIME
							for t in 1:NUM_PIXEL
								CoverageSatLat[(k+(l*N_RAAN)+(s*N_RAAN*N_inclinaison)+(j*N_RAAN*N_inclinaison*N_meanAnomaly)),p,t] =
									abs(Latitude[t] - lat_Sat[p])
								CoverageSatLong[(k+(l*N_RAAN)+(s*N_RAAN*N_inclinaison)+(j*N_RAAN*N_inclinaison*N_meanAnomaly)),p,t] =
									(abs(Longitude[t] - long_Sat[p])*cos(Latitude[t]))
								Thetamax[(k+(l*N_RAAN)+(s*N_RAAN*N_inclinaison)+(j*N_RAAN*N_inclinaison*N_meanAnomaly)),p,t] = Theta[j+1]
							end
						end
					end
				end
			end
		end

		save("HyperMatrixLat_"*string(dt)*"sec_"*string(NUM_PIXEL)*"Targets.jld2","CoverageSatLat",CoverageSatLat)
		save("HyperMatrixLong_"*string(dt)*"sec_"*string(NUM_PIXEL)*"Targets.jld2","CoverageSatLong",CoverageSatLong)
		save("HyperMatrixThetamax_"*string(dt)*"sec_"*string(NUM_PIXEL)*"Targets.jld2","Thetamax",Thetamax)
	else
		CoverageSatLat = load("HyperMatrixLat_"*string(dt)*"sec_"*string(NUM_PIXEL)*"Targets.jld2","CoverageSatLat")
		CoverageSatLong = load("HyperMatrixLong_"*string(dt)*"sec_"*string(NUM_PIXEL)*"Targets.jld2","CoverageSatLong")
		Thetamax = load("HyperMatrixThetamax_"*string(dt)*"sec_"*string(NUM_PIXEL)*"Targets.jld2","Thetamax")
	end

	return CoverageSatLat,CoverageSatLong,Thetamax
end

"""
Get distance in lat and long from satellite to target for all time for a specific orbital indicated by INDEX. It give also ThetaMax the limit coverage of satellite.
"""
function getCoverageSat(INDEX,NUM_PIXEL,Latitude,Longitude)
	IndexAltitude,Indexinclination,IndexRaan,IndexMeanAnom = ExtractKepParamFromINDEX(INDEX)
	LenIndex = length(INDEX)
	CoverageSatLat = fill(0.0,LenIndex,NUM_TIME,NUM_PIXEL)
	CoverageSatLong = fill(0.0,LenIndex,NUM_TIME,NUM_PIXEL)
	Thetamax = zeros(LenIndex,NUM_TIME,NUM_PIXEL)

	for i in 1:LenIndex
		lat_Sat,long_Sat = ProjSatPlotPOLARDominique(R_terre+altitudeSet[IndexAltitude[i]],InclinaisonSet[Indexinclination[i]],RAANSet[IndexRaan[i]],MeanAnomalySet[IndexMeanAnom[i]])
		for p in 1:NUM_TIME
			for t in 1:NUM_PIXEL
				CoverageSatLat[i,p,t] =	abs(Latitude[t] - lat_Sat[p])
				CoverageSatLong[i,p,t] = (abs(Longitude[t] - long_Sat[p])*cos(Latitude[t]))
				Thetamax[i,p,t] = Theta[IndexAltitude[i]]
			end
		end
	end

	return CoverageSatLat,CoverageSatLong,Thetamax
end

"""
Output all the Keplerian's set possible, the cardinality of the sets and the total number of possible orbitals
"""
function getCardinalKeplerianParam()

	InclinaisonSet = Discretization(numbOfDiscretize,0.0,pi)
	RAANSet = Discretization(numbOfDiscretize,0.0,2*pi)
	MeanAnomalySet = Discretization(numbOfDiscretize,0.0,2*pi)
	N_inclinaison = size(InclinaisonSet)[1]
	N_RAAN = size(RAANSet)[1]
	N_meanAnomaly = size(MeanAnomalySet)[1]
	NbrTOTALsat = N_inclinaison*N_RAAN*N_meanAnomaly*altitudeLength

	return NbrTOTALsat,N_inclinaison,N_RAAN,N_meanAnomaly,InclinaisonSet,RAANSet,MeanAnomalySet
end

"""
From a i-est orbital indicated by INDEX, the function send the index of Kepelrrian's sets corresponding to the i-est orbital
"""
function ExtractKepParamFromINDEX(INDEX)

	IndexMeanAnom = fill(1,length(INDEX))
	IndexRaan = fill(1,length(INDEX))
	Indexinclination = fill(1,length(INDEX))
	IndexAltitude = fill(1,length(INDEX))
	#AltProduct=N_meanAnomaly*N_RAAN*N_inclinaison
	#IncProduct=N_meanAnomaly*N_RAAN

	global indtemp = 0
	for i in INDEX
		global indtemp +=1
		IndexMeanAnom[indtemp] = (i%N_meanAnomaly)
		IndexAltitude[indtemp] = div(i,(N_meanAnomaly*N_RAAN*N_inclinaison)) + 1
		Indexinclination[indtemp] = div((i-((IndexAltitude[indtemp]-1)*(N_meanAnomaly*N_RAAN*N_inclinaison))),
											(N_meanAnomaly*N_RAAN)) + 1
		IndexRaan[indtemp] = div((i-((IndexAltitude[indtemp]-1)*(N_meanAnomaly*N_RAAN*N_inclinaison))) -
											(Indexinclination[indtemp]-1)*(N_meanAnomaly*N_RAAN),N_meanAnomaly) + 1
		if IndexMeanAnom[indtemp]==0
			IndexMeanAnom[indtemp]=N_meanAnomaly
			IndexRaan[indtemp] -= 1
		end

		#if i != (((IndexAltitude[indtemp]-1)*AltProduct)+
		#	((Indexinclination[indtemp]-1)*IncProduct)+
		#	((IndexRaan[indtemp]-1)*N_meanAnomaly)+
		#	IndexMeanAnom[indtemp])
		#	println("ERROR to take Kep. params from INDEX (RMP.jl)")
		#end
	end

	return IndexAltitude,Indexinclination,IndexRaan,IndexMeanAnom
end