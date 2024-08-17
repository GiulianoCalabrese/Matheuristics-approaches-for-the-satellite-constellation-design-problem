using LinearAlgebra, PyPlot, SatelliteToolbox
include("Parameters.jl")

# Obtenir coordonnées polaires des cartésiennes
function CartesianToPolar(x,y,z)
	# position of projection of satellite on Earth (ground track)
	ρ = sqrt(x^2+y^2+z^2)
	θ = acos(z/ρ)
	ϕ = atan(y,x)
	return ρ,θ,ϕ
end

# Obtenir coordonnées cartésiennes des polaires
function PolarToCartesian(ρ,θ,ϕ)
	# position of projection of satellite on Earth (ground track)
	x = ρ*cos(ϕ)*sin(θ)
	y = ρ*sin(ϕ)*sin(θ)
	z = ρ*cos(θ)
	return x,y,z
end

function greenwitchsideraltime(PrimoGiorno)
	julianDate = PrimoGiorno
	dut1 = 0.0
	julianDateUt1= (julianDate + dut1 / 86400.0);

	julianCenturiesUt1 = (julianDateUt1 - 2451545.0) / 36525.0;
	greenwichSideralTime = -6.2e-6 * julianCenturiesUt1 * julianCenturiesUt1 * julianCenturiesUt1
							+ 0.093104 * julianCenturiesUt1 * julianCenturiesUt1
							+ (876600.0 * 3600 + 8640184.812866) * julianCenturiesUt1
							+ 67310.54841;  # sec
	greenwichSideralTime = greenwichSideralTime * pi/180 / 240.0% 2pi # //360/86400 = 1/240 : sec to deg. And then deg to rad

	return greenwichSideralTime
end

# Coordonnées cartésiennes en 3D des points à surveiller périodiquement
function generateTargetsCartesians(NUM_PIXEL)
	
	local x = fill(0.0, NUM_PIXEL)
	local y = fill(0.0, NUM_PIXEL)
	local z = fill(0.0, NUM_PIXEL)
    
	# https://mathworld.wolfram.com/SpherePointPicking.html
	for i in 1:NUM_PIXEL
		phi  = 2. * pi * rand()
		csth = 1.0 - 2.0 * rand()
		snth = sqrt(1.0 - csth*csth)
		x[i]=RAYON * snth * cos(phi)
		y[i]=RAYON * snth * sin(phi)
		z[i]=RAYON * csth
	end	
	"""
	if NUM_PIXEL == 1
		x = [5.405520928000631e6]
		y = [2.99698027511728e6]
		z = [1.574507983111815e6]
	elseif NUM_PIXEL == 2
		x = [5.405520928000631e6, -2.4786455107910903e6]
		y = [2.99698027511728e6, 5.464405512810984e6]
		z = [1.574507983111815e6, -2.1626861734365863e6]
	end
	"""
	return x,y,z
end

# Position pixel en fonction du temps (inspiré de page Wiki https://fr.wikipedia.org/wiki/Coordonn%C3%A9es_sph%C3%A9riques)
function TargetInTime(NUM_PIXEL,x,y,z)
	global x_t = fill(0.0, NUM_PIXEL, NUM_TIME+1)
	global y_t = fill(0.0, NUM_PIXEL, NUM_TIME+1)
	global z_t = fill(0.0, NUM_PIXEL, NUM_TIME+1)
	for p in 1:NUM_TIME
	# vit:(Rad/s * tps:s) #https://fr.wikipedia.org/wiki/Temps_sid%C3%A9ral
		local dist_ang_pix = (2*pi/TEMPS_SIMU)*(p-1)*dt
		for j in 1:NUM_PIXEL
			atan_longitude_p = atan(y[j],x[j]) + dist_ang_pix
		# Radian : latitude (inclinaison par rapport au plan equatoriale)
			local LatPix = pi/2 - acos(z[j]/RAYON) 
			x_t[j,p] = RAYON * cos(atan_longitude_p) * cos(LatPix)
			y_t[j,p] = RAYON * sin(atan_longitude_p) * cos(LatPix)
			z_t[j,p] = z[j]
		end
	end
	
	return x_t,y_t,z_t
end

# Position pixel en fonction du temps (inspiré de page Wiki https://fr.wikipedia.org/wiki/Coordonn%C3%A9es_sph%C3%A9riques)
function TargetInTimeECI(NUM_PIXEL,x,y,z)
	global x_t = fill(0.0, NUM_PIXEL, NUM_TIME )
	global y_t = fill(0.0, NUM_PIXEL, NUM_TIME)
	global z_t = fill(0.0, NUM_PIXEL, NUM_TIME)
	for p in 1:NUM_TIME
	# vit:(Rad/s * tps:s) #https://fr.wikipedia.org/wiki/Temps_sid%C3%A9ral
		local dist_ang_pix = (2*pi/TEMPS_SIMU)*(p-1)*dt
		global PrimoGiorno = DatetoJD(1961,4,12,0,0,0) #Date initiale de propagation
		local matrix_rotation = r_ecef_to_eci(PEF(), J2000(), PrimoGiorno+(dt*p))
		
		for j in 1:NUM_PIXEL
			atan_longitude_p = atan(y[j],x[j]) + dist_ang_pix
		# Radian : latitude (inclinaison par rapport au plan equatoriale)
			local LatPix = pi/2 - acos(z[j]/RAYON) 
			x_temp = RAYON * cos(atan_longitude_p) * cos(LatPix)
			y_temp = RAYON * sin(atan_longitude_p) * cos(LatPix)
			z_temp = z[j]
			local ECEFtoECI = matrix_rotation*([x_temp,y_temp,z_temp])
			x_t[j,p] = ECEFtoECI[1]
			y_t[j,p] = ECEFtoECI[2]
			z_t[j,p] = ECEFtoECI[3]
		end
	end
	
	return x_t,y_t,z_t
end

# Discretization 
function Discretization(numbOfDiscretize,lowerBound,upperBound)
	angle_discret = LinRange(lowerBound,upperBound,numbOfDiscretize)
	return angle_discret
end

function TargetSatPlot(numbOfDiscretize, j, k, l, s)
	# Discretization Keplerian parameters
	inclinaisonSet = Discretization(numbOfDiscretize,pi)
	noeudAscendantSet = Discretization(numbOfDiscretize,2*pi)
	meanAnomalySet = Discretization(numbOfDiscretize,2*pi)

	x_projSat1 = fill(0.0,NUM_TIME)
	y_projSat1 = fill(0.0,NUM_TIME)
	z_projSat1 = fill(0.0,NUM_TIME)

	for p=1:NUM_TIME

		x_projSat1[p] = (((cos(meanAnomalySet[k])*cos(sqrt(μ/((RAYON+altitudeSet[j])^3))*(p-1)*dt) - 
		sin(meanAnomalySet[k])*sin(sqrt(μ/((RAYON+altitudeSet[j])^3))*(p-1)*dt)) * cos(noeudAscendantSet[l])*(RAYON+altitudeSet[j])) - 
		((sin(meanAnomalySet[k])*cos(sqrt(μ/((RAYON+altitudeSet[j])^3))*(p-1)*dt) + cos(meanAnomalySet[k])*sin(sqrt(μ/((RAYON+altitudeSet[j])^3))*(p-1)*dt))
		* cos(inclinaisonSet[s])*sin(noeudAscendantSet[l])*(RAYON+altitudeSet[j]))) 

		y_projSat1[p] = (((cos(meanAnomalySet[k])*cos(sqrt(μ/((RAYON+altitudeSet[j])^3))*(p-1)*dt) - 
		sin(meanAnomalySet[k])*sin(sqrt(μ/((RAYON+altitudeSet[j])^3))*(p-1)*dt)) * 
		sin(noeudAscendantSet[l])*(RAYON+altitudeSet[j])) + ((sin(meanAnomalySet[k])*cos(sqrt(μ/((RAYON+altitudeSet[j])^3))*(p-1)*dt) + 
		cos(meanAnomalySet[k])*sin(sqrt(μ/((RAYON+altitudeSet[j])^3))*(p-1)*dt)) * cos(noeudAscendantSet[l])*cos(inclinaisonSet[s])*(RAYON+altitudeSet[j]))) 
			
		z_projSat1[p] = ((sin(meanAnomalySet[k])*cos(sqrt(μ/((RAYON+altitudeSet[j])^3))*(p-1)*dt) + 
		cos(meanAnomalySet[k])*sin(sqrt(μ/((RAYON+altitudeSet[j])^3))*(p-1)*dt)) * sin(inclinaisonSet[s])*(RAYON+altitudeSet[j])) 

	end
	return x_projSat1,y_projSat1,z_projSat1
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

# Temps de passage
function CountRevisitTime(TimeSeen_track)
	for index_ville in 1:size(TimeSeen_track)[1]
		figure()
		PyPlot.title("Target :"*string(index_ville))
		PyPlot.xlabel("Time SIMULATION in seconds")
		PyPlot.ylabel("Visibility")
		PyPlot.yticks([])
		PyPlot.plot(collect(0:dt:TEMPS_SIMU),TimeSeen_track[index_ville,:])
		println("number revisit time ",string(sum(TimeSeen_track[index_ville,:])),"; number minimal revisit time ",string(TEMPS_SIMU/DELTA_T))
		Time_intervall_revisit = collect(0:convert(Int64,(DELTA_T)/dt):convert(Int64,TEMPS_SIMU/dt))
		for time_inter in Time_intervall_revisit
			interval_time = 0*collect(0:dt:TEMPS_SIMU)
			interval_time[time_inter+1] = 1.25
			PyPlot.plot(collect(0:dt:TEMPS_SIMU),interval_time, color = "red")
		end
	end
end
############### PLOT ###############

function CreateSfera(RAYON_sfera,x_init=0.0,y_init=0.0,z_init=0.0,N=32,colorSfera="gray",Axes="True")
	# Create a sfera
	phi = range(0, stop=2*pi, length=N)
	theta = range(0, stop=pi, length=N)
	x = RAYON_sfera*cos.(phi) .* sin.(theta)'
	y = RAYON_sfera*sin.(phi) .* sin.(theta)'
	z = RAYON_sfera*repeat(cos.(theta)',outer=[N, 1])
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

# Create Piano equatoriale
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

function plot_target(x_target,y_target,z_target)
	CreateEquatorialPlane()
	CreateSfera()
	
	PyPlot.title("NUMBER of TARGET : "*string(length(x_target)))
	PyPlot.scatter3D(x_target,y_target,z_target, color = "blue")#, marker = "o", markersize = 5)	
end

# Plot orbitals
function OrbitalPlanes(RAYON_SAT,x_init,y_init,z_init,x_sat,y_sat,z_sat,N = 10)

	for i in 1:length(x_init)
		# SPhère couverture sat
		#CreateSfera(RAYON_SAT,x_init[i],y_init[i],z_init[i],N,"orange","False")		
	end
	# Position sat en ECEF
	PyPlot.scatter3D(x_sat,y_sat,z_sat,color = "violet", label = "Optimal Satellite_ECEF")
	# projection sat sur surface en ECEF
	PyPlot.scatter3D(x_init,y_init,z_init,color = "green", label = "PROJECTION Satellite_ECEF")
end