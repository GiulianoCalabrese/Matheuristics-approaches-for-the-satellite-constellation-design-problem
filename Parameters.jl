using SatelliteToolbox
include("SCDP/SCDP/src/constantes_physiques.jl")
#################### LIST of PARAMETERS ######################

# RAYON de la Terre
	global RAYON = 6378136.3 # m
	
# Circumeference of a sphere (Earth)
	global CircumferenceEarth = RAYON*2*pi # Km

# hai [1261.9898537417057, 893.6996722009753, 566.8052002788218] come differenti altitudini per 24h di simulazione
# ottieni questo con apertura 40 gradi come lato in KM del quadrato del satellite : 930.1114958405926 km   655.7710767836845 km    414.35089323796853 km

# Temps total de simulation
	global TEMPS_SIMU = 24.0*3600 # sec
	# sideral time in seconds 23*3600 + 56*60 + 4

# Pas de temps de la simulation
	global dt = 180.0# sec
	
# Nombre de pas de temps dans la simulation
	global NUM_TIME = convert(Int,round(TEMPS_SIMU/dt))

# Earth's rotation rate in rad/s
	global we= 7.2921e-5

# Nombre de points à observer sur la Terre (nous rappelons que la Terre est une sphère discrétisée)
	global RANGE_NUM_PIXEL = 1:1#1*fill(1, 13)#1:2
	
# Obtention ensemble de valeurs admis pour hauteur satellite en fonction du nombre de N revolutions du satellite
	N_test=1:24
	global altitudeSet = Float64[]
	global PERIOD_SAT = Float64[]
	for i in N_test
		altitudeSetVal = ( ( μ*((TEMPS_SIMU)/N_test[i])^2 ) / 4π^2 )^(1/3) - RAYON
		if altitudeSetVal >= 400000 && altitudeSetVal <= 1400000
			append!(altitudeSet, [altitudeSetVal])
			append!(PERIOD_SAT, [TEMPS_SIMU/N_test[i]])
		end
	end
	# println(altitudeSet)
	
	global PERIOD_SAT_MIN = PERIOD_SAT[end]
	
# Initialize variable in model
	Altitude =  [7.640227192785597e6, 7.640227192785597e6]
	Inc =  [1.2161003820347587, 2.026833970057931]
	RAAN =  [1.1469306513105595, 1.6455961518803681]
	MeanAnomaly =   [5.6847867064958155, 1.296530301481502]
	
	global start_variable = 0.0
	global start_altitude = [[0,1,0],[0,1,0]]
	global start_activation = [2,2,2,1,1,1,1,1,2]
	global start_sin_inclinaison = sin.(Inc)
	global start_cos_inclinaison = cos.(Inc)
	global start_sin_noeudAscendant = sin.(RAAN)
	global start_cos_noeudAscendant = cos.(RAAN)
	global start_sin_meanAnomaly = sin.(MeanAnomaly)
	global start_cos_meanAnomaly = cos.(MeanAnomaly)
	
	global AltitudeVariable = altitudeSet[1]
	global PERIODVariable = PERIOD_SAT[1]
	global altitudeLength = length(altitudeSet)

	global Modele = Val(:twobody)
	global color = ["brown","green","red","blue","violet","black"]
	global M_limit = 2*pi
	global M_limit_1 = 2*pi
	global M_limit_lat = 2*pi
	global M_limit_long = 2*pi
	global time_zero_simulation = date_to_jd(2000,1,1,12,0,0)#1970
	
# Compute Theta and alpha of a satellite
	global Theta_min = ((2*pi*dt)/(2*PERIOD_SAT_MIN))*1.2
	global numbOfDiscretize = Int(ceil(pi/Theta_min))
	global alphaHalf = atan((RAYON/(RAYON+altitudeSet[end]))*tan(Theta_min))

	global Theta = fill(0.0,length(altitudeSet))
	for j in 1:altitudeLength
		if ((RAYON+altitudeSet[j])/RAYON)*sin(alphaHalf) > 1
			alpha_lim = asin(RAYON / (RAYON + altitudeSet[j]))
			global Theta[j] = -alpha_lim + asin(((RAYON+altitudeSet[j])/RAYON)*sin(alpha_lim))
		else
			global Theta[j] = (-alphaHalf + asin(((RAYON+altitudeSet[j])/RAYON)*sin(alphaHalf)))
		end
	end
	
# Discretization 
	function Discretization(numbOfDiscretize,lowerBound,upperBound)
		angle_discret = LinRange(lowerBound,upperBound,numbOfDiscretize)
		return angle_discret
	end
	
	global inclinaisonSet = Discretization(numbOfDiscretize,Inc[1],1*pi+Inc[1])
	global noeudAscendantSet = Discretization(numbOfDiscretize,0.0,2*pi)
	global meanAnomalySet = copy(noeudAscendantSet)

	probasMut = [0.95]
	probasModifNbSats = [0.05]
	nbMaxGen = 500 # MODIFFF
	ϵ = 0.005
	nbIndiv = 1500 # MODIFFFF
	stats = Vector{Float64}[]
	
# Target Position
	# global lat = [22,31,45,51,60]*(pi/180)
	# global long = [206, 18, 56, 120, 342]*(pi/180)
	nbCibles = 1
	# @load "instances/instance_"*string(nbCibles)*"_cibles.jld2" inst
	# global lat = fill(0.0,nbCibles)
	# global long = fill(0.0,nbCibles)
	# for i in 1:nbCibles
		# lat[i] = inst.coordsCibles[i][1]*(180/pi)
		# long[i] = inst.coordsCibles[i][2]*(180/pi)
	# end
	# inst = load("instances/instance_"*string(nbCibles)*"_cibles.jld2", "inst")
	Latitude = [30.1289733044805, 11.202345525890758] *pi/180
	Longitude = [21.921097190662763, 134.10212177710707]*pi/180 
	coordsCibles = [SA[Latitude[i], Longitude[i]] for i in 1:nbCibles]
	n_0 = 2
	nbMaxSat = 4
	inst = DonneesSCDP(coordsCibles,alphaHalf,nbMaxSat, SatelliteToolbox.DateTime(1970, 1, 1, 0),1,n_0,dt)
	# inst = instanceAleatoire(nbCibles,alphaHalf,nbMaxSat,SatelliteToolbox.DateTime(1970, 1, 1, 0),1,n_0)
	# save("instances/instance_Revisit_"*string(Int(TEMPS_SIMU/(3600*n_0)))*"_cibles_"*string(nbCibles)*".jld2", "inst", inst)
	ecrireInst(inst)
	
# Function used in NSGAII	
	function fDiag(nbIndiv; ite, pop, donneesCrit)
	   fCVMoy = sum(getfield.(pop[1:nbIndiv], :CV))/nbIndiv
	   pourcentIndiv = 100sum(getfield.(pop[1:nbIndiv], :CV) .≈ 0)/nbIndiv
	   nbSatMoy = sum( getindex.(getfield.(pop[1:nbIndiv], :y), 1) )/nbIndiv
	   nbSatMin = pop[1].y[1]
	   push!(stats, SA[fCVMoy, pourcentIndiv, nbSatMoy, nbSatMin])
	end

"""
Renvoie l'angle dans ]-π, π] de la rotation entre le référentiel ECI et le référentiel ECEF
au jour julien date. 
"""
	function calcul_angle_ECI_ECEF(date::Real)
	   mat_ECI_ECEF = r_eci_to_ecef(TOD(), PEF(), date)
	   return atan(-mat_ECI_ECEF[2, 1], mat_ECI_ECEF[1, 1])
	end

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

# parameters for GA.jl
	N_test=1:24
	global altitudeSet = Float64[]
	global PERIOD_SAT = Float64[]
	for i in N_test
		altitudeSetVal = ( ( μ*((T_terre)/N_test[i])^2 ) / 4π^2 )^(1/3) - R_terre
		if altitudeSetVal >= 400000 && altitudeSetVal <= 1400000
			append!(altitudeSet, [altitudeSetVal])
			append!(PERIOD_SAT, [T_terre/N_test[i]])
		end
	end

	global PERIOD_SAT_MIN = PERIOD_SAT[end]
	global theta_min = ((2*pi*dt)/(2*PERIOD_SAT_MIN))*1.2
	global numbOfDiscretize = Int(ceil(pi/theta_min))
	global alphaHalf = atan((sin(theta_min))/(((R_terre+altitudeSet[end])/R_terre)-cos(theta_min)))

	Theta = fill(0.0,length(altitudeSet))
	for j in 1:length(altitudeSet)
		if ((R_terre+altitudeSet[j])/R_terre)*sin(alphaHalf) > 1
			alpha_lim = asin(R_terre / (R_terre + altitudeSet[j]))
			global Theta[j] = (-alpha_lim + asin(round(((R_terre+altitudeSet[j])/R_terre)*sin(alpha_lim),digits=6)))
		else	
			global Theta[j] = (-alphaHalf + asin(((R_terre+altitudeSet[j])/R_terre)*sin(alphaHalf)))
		end
	end
