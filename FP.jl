using JuMP, CPLEX, Ipopt, Distributions
# include("../Functions.jl")
# include("Parameters.jl")

function FpModel(nbMaxSat,NUM_PIXEL,lat,long,PERIOD,NUM_PERIOD)

	global angle0 = calcul_angle_ECI_ECEF(time_zero_simulation)
	global altitude_final
	global sin_inclinaison_final
	global cos_inclinaison_final
	global sin_noeudAscendant_final
	global cos_noeudAscendant_final
	global sin_meanAnomaly_final
	global cos_meanAnomaly_final
	global start_variable
	start_variable = 0.0
	global num_local = 1
	global objective
	global lat_cos = cos.(lat)
	global instance = 0

	global instance += 1

	status = "nothing"
	println("nombre target "*string(NUM_PIXEL))

	global NUM_SATELLITE = 0
	flag = 1
	global cpu_time = 0

	@time begin
		while NUM_SATELLITE < nbMaxSat && flag==1
			NUM_SATELLITE = NUM_SATELLITE+1 
			println("nombre satellites "*string(NUM_SATELLITE))
			global iter_max = 20 * NUM_SATELLITE
			# Constante de linéarisation (MINLP)
			M_LIMIT = NUM_PERIOD * NUM_SATELLITE
			M_LIMIT_1 = 3*pi

			global model = Model(Ipopt.Optimizer)
			model1 = Model(CPLEX.Optimizer)
			set_optimizer_attribute(model1, "CPX_PARAM_TILIM", 3600)
			set_optimizer_attribute(model1, "CPX_PARAM_SCRIND", 0)
			#set_optimizer_attribute(model1, "CPX_PARAM_EPRHS", 10e-2)
			set_optimizer_attribute(model, "max_cpu_time", 3600.0)
			set_optimizer_attribute(model, "max_iter", 1000000)
			set_optimizer_attribute(model, "tol", 10e-2)
			set_optimizer_attribute(model, "print_level", 0)
			JuMP.register(model, :atan, 2, atan, autodiff=true)
			JuMP.register(model, :sin, 1, sin, autodiff=true)
			JuMP.register(model, :%, 2, %; autodiff = true)

			# Altitude de l'orbite [Km]
			@variable(model, altitude[1:NUM_SATELLITE] >= RAYON)
			@variable(model, 0 <= activation_altitude[1:NUM_SATELLITE,1:length(altitudeSet)] <= 1)
			@variable(model, activation_altitude_old[1:NUM_SATELLITE,1:length(altitudeSet)])

			# inclinaison de l'orbite
			@variable(model, 0 <= sin_inclinaison[1:NUM_SATELLITE] <= 1)
			@variable(model, -1 <= cos_inclinaison[1:NUM_SATELLITE] <= 1)

			# Noeud ascendant
			@variable(model, -1 <= sin_noeudAscendant[1:NUM_SATELLITE] <= 1)
			@variable(model, -1 <= cos_noeudAscendant[1:NUM_SATELLITE] <= 1)

			# Anomalie moyenne du plan orbit
			@variable(model, -1 <= sin_meanAnomaly[1:NUM_SATELLITE] <= 1)
			@variable(model, -1 <= cos_meanAnomaly[1:NUM_SATELLITE] <= 1)

			# Identité trigonométrique
			@constraint(model, [i=1:NUM_SATELLITE], sin_inclinaison[i]^2 + cos_inclinaison[i]^2 == 1)
			@constraint(model, [i=1:NUM_SATELLITE], sin_noeudAscendant[i]^2 + cos_noeudAscendant[i]^2 == 1)
			@constraint(model, [i=1:NUM_SATELLITE], sin_meanAnomaly[i]^2 + cos_meanAnomaly[i]^2 == 1)

			# Compute the geocentric distance.
			@constraint(model, [i=1:NUM_SATELLITE], sum(activation_altitude[i,j] for j = 1:length(altitudeSet))==1)
			@constraint(model, [i=1:NUM_SATELLITE], altitude[i] == RAYON + sum(altitudeSet[j]*activation_altitude[i,j] for j=1:length(altitudeSet)))

			# Thetamax permet de déterminer la surface couverte par le satellite à partir du demi angle d'ouverture des capteurs, de l'altitude et de l'excentricité
			@variable(model, Thetamax[1:NUM_SATELLITE] >= theta_min)
			@constraint(model, [i=1:NUM_SATELLITE], Thetamax[i] == (-alphaHalf + sum(activation_altitude[i,j]*asin(((RAYON+altitudeSet[j])/RAYON)*sin(alphaHalf))
			for j=1:length(altitudeSet))))

			# Variables de linéarisation (MINLP)
			@variable(model, 0 <= Xi[1:NUM_SATELLITE, 1:NUM_TIME, 1:NUM_PIXEL] <= 1)
		
			@NLexpression(model, LatitudeSat[i=1:NUM_SATELLITE, p=1:NUM_TIME], ((asin(((sin_inclinaison[i]*altitude[i]*sin_meanAnomaly[i])*
			cos(sqrt(μ/(altitude[i]^3))*((p-1)*dt))/altitude[i]) + ((sin_inclinaison[i]*sqrt(μ/(altitude[i]))*cos_meanAnomaly[i])*sin(sqrt(μ/(altitude[i]^3))*((p-1)*dt))*sqrt(altitude[i]/μ))))))

				@NLconstraint(model, [p=1:NUM_TIME,j=1:NUM_PIXEL,i=1:NUM_SATELLITE], abs(- (asin(((sin_inclinaison[i]*altitude[i]*sin_meanAnomaly[i])*
				cos(sqrt(μ/(altitude[i]^3))*((p-1)*dt))/altitude[i]) + ((sin_inclinaison[i]*sqrt(μ/(altitude[i]))*cos_meanAnomaly[i])*sin(sqrt(μ/(altitude[i]^3))*((p-1)*dt))*sqrt(altitude[i]/μ)))) + lat[j]) <= Thetamax[i] + (1-Xi[i,p,j])*M_limit)


			# @NLconstraint(model, [p=1:NUM_TIME,j=1:NUM_PIXEL,i=1:NUM_SATELLITE], abs(- (asin(((sin_inclinaison[i]*altitude[i]*sin_meanAnomaly[i])*
			# cos(sqrt(μ/(altitude[i]^3))*((p-1)*dt))/altitude[i]) + ((sin_inclinaison[i]*sqrt(μ/(altitude[i]))*cos_meanAnomaly[i])*sin(sqrt(μ/(altitude[i]^3))*((p-1)*dt))*sqrt(altitude[i]/μ)))) + lat[j]) <= Thetamax[i] + (1-Xi[i,p,j])^2*M_limit_lat*10)

			# @NLconstraint(model, [p=1:NUM_TIME,j=1:NUM_PIXEL,i=1:NUM_SATELLITE], abs(- (asin(((sin_inclinaison[i]*altitude[i]*sin_meanAnomaly[i])*
			# cos(sqrt(μ/(altitude[i]^3))*((p-1)*dt))/altitude[i]) + ((sin_inclinaison[i]*sqrt(μ/(altitude[i]))*cos_meanAnomaly[i])*sin(sqrt(μ/(altitude[i]^3))*((p-1)*dt))*sqrt(altitude[i]/μ)))) + lat[j]) >= Thetamax[i] - Xi[i,p,j]^2*M_limit_lat*10)

			# @NLexpression(model, LongitudeSat[i=1:NUM_SATELLITE, p=1:NUM_TIME], ((((((-(angle0 + (we*((p-1)*dt))) +
			# atan(((sin_noeudAscendant[i]*altitude[i]*cos_meanAnomaly[i])
			# + (cos_noeudAscendant[i]*cos_inclinaison[i]*altitude[i]*sin_meanAnomaly[i]))*cos(((sqrt(μ/(altitude[i]^3)))*((p-1)*dt))/altitude[i])
			# + ((-(sin_noeudAscendant[i]*sqrt(μ/(altitude[i]))*sin_meanAnomaly[i])
			# + (cos_noeudAscendant[i]*cos_inclinaison[i]*sqrt(μ/(altitude[i]))*cos_meanAnomaly[i]))*sin(((sqrt(μ/(altitude[i]^3)))*((p-1)*dt))*sqrt(altitude[i]/μ))),
			# ((cos_noeudAscendant[i]*altitude[i]*cos_meanAnomaly[i])
			# - (sin_noeudAscendant[i]*cos_inclinaison[i]*altitude[i]*sin_meanAnomaly[i]))*
			# cos(((sqrt(μ/(altitude[i]^3)))*((p-1)*dt))/altitude[i])
			# + ((-(cos_noeudAscendant[i]*sqrt(μ/(altitude[i]))*sin_meanAnomaly[i])-(sin_noeudAscendant[i]*cos_inclinaison[i]*sqrt(μ/(altitude[i]))*
			# cos_meanAnomaly[i]))*sin(((sqrt(μ/(altitude[i]^3)))*((p-1)*dt))*sqrt(altitude[i]/μ))))) %(2*pi)) + (2*pi))%(2*pi)))))

			# @NLconstraint(model, [p=1:NUM_TIME,j=1:NUM_PIXEL,i=1:NUM_SATELLITE], abs(- (((((-(angle0 + (we*((p-1)*dt))) +
			# atan(((sin_noeudAscendant[i]*altitude[i]*cos_meanAnomaly[i])
			# + (cos_noeudAscendant[i]*cos_inclinaison[i]*altitude[i]*sin_meanAnomaly[i]))*cos(((sqrt(μ/(altitude[i]^3)))*((p-1)*dt))/altitude[i])
			# + ((-(sin_noeudAscendant[i]*sqrt(μ/(altitude[i]))*sin_meanAnomaly[i])
			# + (cos_noeudAscendant[i]*cos_inclinaison[i]*sqrt(μ/(altitude[i]))*cos_meanAnomaly[i]))*sin(((sqrt(μ/(altitude[i]^3)))*((p-1)*dt))*sqrt(altitude[i]/μ))),
			# ((cos_noeudAscendant[i]*altitude[i]*cos_meanAnomaly[i])
			# - (sin_noeudAscendant[i]*cos_inclinaison[i]*altitude[i]*sin_meanAnomaly[i]))*
			# cos(((sqrt(μ/(altitude[i]^3)))*((p-1)*dt))/altitude[i])
			# + ((-(cos_noeudAscendant[i]*sqrt(μ/(altitude[i]))*sin_meanAnomaly[i])-(sin_noeudAscendant[i]*cos_inclinaison[i]*sqrt(μ/(altitude[i]))*
			# cos_meanAnomaly[i]))*sin(((sqrt(μ/(altitude[i]^3)))*((p-1)*dt))*sqrt(altitude[i]/μ))))) ) %(2*pi)) + (2*pi))%(2*pi)) + long[j]) <= Thetamax[i] + (1-Xi[i,p,j])^2*M_limit_long*10)

			# @NLconstraint(model, [p=1:NUM_TIME,j=1:NUM_PIXEL,i=1:NUM_SATELLITE], abs(- (((((-(angle0 + (we*((p-1)*dt))) +
			# atan(((sin_noeudAscendant[i]*altitude[i]*cos_meanAnomaly[i])
			# + (cos_noeudAscendant[i]*cos_inclinaison[i]*altitude[i]*sin_meanAnomaly[i]))*cos(((sqrt(μ/(altitude[i]^3)))*((p-1)*dt))/altitude[i])
			# + ((-(sin_noeudAscendant[i]*sqrt(μ/(altitude[i]))*sin_meanAnomaly[i])
			# + (cos_noeudAscendant[i]*cos_inclinaison[i]*sqrt(μ/(altitude[i]))*cos_meanAnomaly[i]))*sin(((sqrt(μ/(altitude[i]^3)))*((p-1)*dt))*sqrt(altitude[i]/μ))),
			# ((cos_noeudAscendant[i]*altitude[i]*cos_meanAnomaly[i])
			# - (sin_noeudAscendant[i]*cos_inclinaison[i]*altitude[i]*sin_meanAnomaly[i]))*
			# cos(((sqrt(μ/(altitude[i]^3)))*((p-1)*dt))/altitude[i])
			# + ((-(cos_noeudAscendant[i]*sqrt(μ/(altitude[i]))*sin_meanAnomaly[i])-(sin_noeudAscendant[i]*cos_inclinaison[i]*sqrt(μ/(altitude[i]))*
			# cos_meanAnomaly[i]))*sin(((sqrt(μ/(altitude[i]^3)))*((p-1)*dt))*sqrt(altitude[i]/μ))))) ) %(2*pi)) + (2*pi))%(2*pi)) + long[j]) >= Thetamax[i] - Xi[i,p,j]^2*M_limit_long*10)
			
				
				@NLexpression(model, LongitudeSat[i=1:NUM_SATELLITE, p=1:NUM_TIME], ((((((-(angle0 + (we*((p-1)*dt))) + 
				atan(((sin_noeudAscendant[i]*altitude[i]*cos_meanAnomaly[i])
				+ (cos_noeudAscendant[i]*cos_inclinaison[i]*altitude[i]*sin_meanAnomaly[i]))*cos(((sqrt(μ/(altitude[i]^3)))*((p-1)*dt))/altitude[i])
				+ ((-(sin_noeudAscendant[i]*sqrt(μ/(altitude[i]))*sin_meanAnomaly[i])
				+ (cos_noeudAscendant[i]*cos_inclinaison[i]*sqrt(μ/(altitude[i]))*cos_meanAnomaly[i]))*sin(((sqrt(μ/(altitude[i]^3)))*((p-1)*dt))*sqrt(altitude[i]/μ))),
				((cos_noeudAscendant[i]*altitude[i]*cos_meanAnomaly[i])
				- (sin_noeudAscendant[i]*cos_inclinaison[i]*altitude[i]*sin_meanAnomaly[i]))*
				cos(((sqrt(μ/(altitude[i]^3)))*((p-1)*dt))/altitude[i])
				+ ((-(cos_noeudAscendant[i]*sqrt(μ/(altitude[i]))*sin_meanAnomaly[i])-(sin_noeudAscendant[i]*cos_inclinaison[i]*sqrt(μ/(altitude[i]))*
				cos_meanAnomaly[i]))*sin(((sqrt(μ/(altitude[i]^3)))*((p-1)*dt))*sqrt(altitude[i]/μ))))) %(2*pi)) + (2*pi))%(2*pi)))))

				@NLconstraint(model, [p=1:NUM_TIME,j=1:NUM_PIXEL,i=1:NUM_SATELLITE], abs(- (((((-(angle0 + (we*((p-1)*dt))) + 
				atan(((sin_noeudAscendant[i]*altitude[i]*cos_meanAnomaly[i])
				+ (cos_noeudAscendant[i]*cos_inclinaison[i]*altitude[i]*sin_meanAnomaly[i]))*cos(((sqrt(μ/(altitude[i]^3)))*((p-1)*dt))/altitude[i])
				+ ((-(sin_noeudAscendant[i]*sqrt(μ/(altitude[i]))*sin_meanAnomaly[i])
				+ (cos_noeudAscendant[i]*cos_inclinaison[i]*sqrt(μ/(altitude[i]))*cos_meanAnomaly[i]))*sin(((sqrt(μ/(altitude[i]^3)))*((p-1)*dt))*sqrt(altitude[i]/μ))),
				((cos_noeudAscendant[i]*altitude[i]*cos_meanAnomaly[i])
				- (sin_noeudAscendant[i]*cos_inclinaison[i]*altitude[i]*sin_meanAnomaly[i]))*
				cos(((sqrt(μ/(altitude[i]^3)))*((p-1)*dt))/altitude[i])
				+ ((-(cos_noeudAscendant[i]*sqrt(μ/(altitude[i]))*sin_meanAnomaly[i])-(sin_noeudAscendant[i]*cos_inclinaison[i]*sqrt(μ/(altitude[i]))*
				cos_meanAnomaly[i]))*sin(((sqrt(μ/(altitude[i]^3)))*((p-1)*dt))*sqrt(altitude[i]/μ))))) %(2*pi)) + (2*pi))%(2*pi))) + long[j]) <= Thetamax[i] + (1-Xi[i,p,j])*M_limit)

				# @NLconstraint(model, [p=1:NUM_TIME,j=1:NUM_PIXEL,i=1:NUM_SATELLITE], abs(- (((((-(angle0 + (we*((p-1)*dt))) +
				# atan(((sin_noeudAscendant[i]*altitude[i]*cos_meanAnomaly[i])
				# + (cos_noeudAscendant[i]*cos_inclinaison[i]*altitude[i]*sin_meanAnomaly[i]))*cos(((sqrt(μ/(altitude[i]^3)))*((p-1)*dt))/altitude[i])
				# + ((-(sin_noeudAscendant[i]*sqrt(μ/(altitude[i]))*sin_meanAnomaly[i])
				# + (cos_noeudAscendant[i]*cos_inclinaison[i]*sqrt(μ/(altitude[i]))*cos_meanAnomaly[i]))*sin(((sqrt(μ/(altitude[i]^3)))*((p-1)*dt))*sqrt(altitude[i]/μ))),
				# ((cos_noeudAscendant[i]*altitude[i]*cos_meanAnomaly[i])
				# - (sin_noeudAscendant[i]*cos_inclinaison[i]*altitude[i]*sin_meanAnomaly[i]))*
				# cos(((sqrt(μ/(altitude[i]^3)))*((p-1)*dt))/altitude[i])
				# + ((-(cos_noeudAscendant[i]*sqrt(μ/(altitude[i]))*sin_meanAnomaly[i])-(sin_noeudAscendant[i]*cos_inclinaison[i]*sqrt(μ/(altitude[i]))*
				# cos_meanAnomaly[i]))*sin(((sqrt(μ/(altitude[i]^3)))*((p-1)*dt))*sqrt(altitude[i]/μ))))) ) %(2*pi)) + (2*pi))%(2*pi)) + long[j]) >= Thetamax[i] - Xi[i,p,j]^2*M_limit_long*10)
				
			
			# @variable(model,L5>=0)
			# @constraint(model, [j=1:NUM_PIXEL,k=1:NUM_PERIOD], sum(Xi[i,p,j]^2 for p in (k-1)*PERIOD+1:(PERIOD*k) for i in 1:NUM_SATELLITE) >= L5)

					@constraint(model, [j=1:NUM_PIXEL,k=1:PERIOD], sum(Xi[i,p,j] for p in (k-1)*NUM_PERIOD+1:(NUM_PERIOD*k) for i in 1:NUM_SATELLITE) >= 1)

			@variable(model, Xi_old[1:NUM_SATELLITE, 1:NUM_TIME, 1:NUM_PIXEL])
			@variable(model, LongitudeSat_old[1:NUM_SATELLITE, 1:NUM_TIME])
			@variable(model, LatitudeSat_old[1:NUM_SATELLITE, 1:NUM_TIME])
			@variable(model, nlp_obj)

			# @NLobjective(model, Min, (sum(abs(Xi[i,p,j]-Xi_old[i,p,j]) for i=1:NUM_SATELLITE for j=1:NUM_PIXEL for p=1:NUM_TIME))
			# + (sum(abs(activation_altitude[i,j]-activation_altitude_old[i,j]) for i=1:NUM_SATELLITE for j=1:length(altitudeSet)))- 100*L5)

				@NLobjective(model, Min, sum((Xi[i,p,j]-Xi_old[i,p,j])^2 for i=1:NUM_SATELLITE for j=1:NUM_PIXEL for p=1:NUM_TIME)
				+ sum((activation_altitude[i,j]-activation_altitude_old[i,j])^2 for i=1:NUM_SATELLITE for j=1:length(altitudeSet)))


			@variable(model1, activation_altitudeA1[1:NUM_SATELLITE,1:length(altitudeSet)])
			@variable(model1, altitude1[1:NUM_SATELLITE])
			@variable(model1, activation_altitude1[1:NUM_SATELLITE,1:length(altitudeSet)], Bin)
			@constraint(model1, [i=1:NUM_SATELLITE], sum(activation_altitude1[i,j] for j = 1:length(altitudeSet))==1)
			@constraint(model1, [i=1:NUM_SATELLITE], altitude1[i] == RAYON + sum(altitudeSet[j]*activation_altitude1[i,j] for j=1:length(altitudeSet)))

			# Variables de linéarisation (MINLP)
			@variable(model1, Xi1[1:NUM_SATELLITE, 1:NUM_TIME, 1:NUM_PIXEL], Bin)
			@variable(model1, Xi1_old[1:NUM_SATELLITE, 1:NUM_TIME, 1:NUM_PIXEL])

			@variable(model1, Thetamax1[1:NUM_SATELLITE])
			@constraint(model1, [i=1:NUM_SATELLITE], Thetamax1[i] == (-alphaHalf + sum(activation_altitude1[i,j]*asin(((RAYON+altitudeSet[j])/RAYON)*sin(alphaHalf))
			for j=1:length(altitudeSet))))

			@variable(model1, -pi/2 <= LatitudeSat1[1:NUM_SATELLITE,1:NUM_TIME] <= pi/2) ## in [-pi/2,pi/2]
			@variable(model1, 0 <= LongitudeSat1[1:NUM_SATELLITE,1:NUM_TIME] <= 2*pi) ## in [0,2pi]

			@constraint(model1, [p=1:NUM_TIME,j=1:NUM_PIXEL,i=1:NUM_SATELLITE], Xi1[i,p,j] => {LatitudeSat1[i,p] - lat[j] <= Thetamax1[i]})
			@constraint(model1, [p=1:NUM_TIME,j=1:NUM_PIXEL,i=1:NUM_SATELLITE], Xi1[i,p,j] => {- LatitudeSat1[i,p] + lat[j] <= Thetamax1[i]})
			@constraint(model1, [p=1:NUM_TIME,j=1:NUM_PIXEL,i=1:NUM_SATELLITE], Xi1[i,p,j] => {(LongitudeSat1[i,p] - long[j]) <= Thetamax1[i]})
			@constraint(model1, [p=1:NUM_TIME,j=1:NUM_PIXEL,i=1:NUM_SATELLITE], Xi1[i,p,j] => {(- LongitudeSat1[i,p] + long[j]) <= Thetamax1[i]})

			@constraint(model1, [j=1:NUM_PIXEL,k=1:PERIOD], sum(Xi1[i,p,j] for p in (k-1)*NUM_PERIOD+1:(NUM_PERIOD*k) for i in 1:NUM_SATELLITE) >= 1)
			@variable(model1, abs[1:NUM_SATELLITE, 1:NUM_TIME, 1:NUM_PIXEL] >= 0)

			@constraint(model1, [i=1:NUM_SATELLITE, p=1:NUM_TIME, j=1:NUM_PIXEL], abs[i,p,j] >= Xi1[i,p,j] - Xi1_old[i,p,j])
			@constraint(model1, [i=1:NUM_SATELLITE, p=1:NUM_TIME, j=1:NUM_PIXEL], abs[i,p,j] >= Xi1_old[i,p,j] - Xi1[i,p,j])

			@variable(model1, LongitudeSat1_old[1:NUM_SATELLITE, 1:NUM_TIME])
			@variable(model1, LatitudeSat1_old[1:NUM_SATELLITE, 1:NUM_TIME])
			@variable(model1, absLat[1:NUM_SATELLITE, 1:NUM_TIME] >= 0)
			@variable(model1, absLong[1:NUM_SATELLITE, 1:NUM_TIME] >= 0)

			@constraint(model1, [i=1:NUM_SATELLITE, p=1:NUM_TIME], absLat[i,p] >= LatitudeSat1[i,p] - LatitudeSat1_old[i,p])
			@constraint(model1, [i=1:NUM_SATELLITE, p=1:NUM_TIME], absLat[i,p] >= LatitudeSat1_old[i,p] - LatitudeSat1[i,p])

			@constraint(model1, [i=1:NUM_SATELLITE, p=1:NUM_TIME], absLong[i,p] >= LongitudeSat1[i,p] - LongitudeSat1_old[i,p])
			@constraint(model1, [i=1:NUM_SATELLITE, p=1:NUM_TIME], absLong[i,p] >= LongitudeSat1_old[i,p] - LongitudeSat1[i,p])

			@variable(model1, obj_milp)
			@variable(model1, absAlt[1:NUM_SATELLITE, 1:length(altitudeSet)] >= 0)
			@objective(model1, Min, 100*sum(absLat[i,p] for i=1:NUM_SATELLITE for p=1:NUM_TIME) + 100*sum(absLong[i,p] for i=1:NUM_SATELLITE for p=1:NUM_TIME)
			+ sum(absAlt[i,k] for i=1:NUM_SATELLITE for k=1:length(altitudeSet)) + sum(abs[i,p,j] for p=1:NUM_TIME for j=1:NUM_PIXEL for i=1:NUM_SATELLITE))

			@variable(model1, abs3[1:NUM_SATELLITE, 1:NUM_TIME, 1:NUM_PIXEL, 1:iter_max], Bin)
			@variable(model1, zz[1:NUM_SATELLITE, 1:NUM_TIME, 1:NUM_PIXEL, 1:iter_max], Bin)

			@constraint(model1, [i=1:NUM_SATELLITE, k=1:length(altitudeSet)], absAlt[i,k] >= activation_altitude1[i,k] - activation_altitudeA1[i,k])
			@constraint(model1, [i=1:NUM_SATELLITE, k=1:length(altitudeSet)], absAlt[i,k] >= activation_altitudeA1[i,k] - activation_altitude1[i,k])
			
			value_nlp = 1
			value_milp = 1

			for i in 1:NUM_SATELLITE
				for p in 1:NUM_TIME
					for j in 1:NUM_PIXEL
						A = rand(0:1)
						fix(Xi_old[i,p,j], 1)
						fix(Xi1_old[i,p,j], 1)
					end
				end
			end

			for i in 1:NUM_SATELLITE
				for k in 1:length(altitudeSet)
					A = rand(0:1)
					# if k == 1
						# fix(activation_altitude_old[i,k], 1)
						# fix(activation_altitudeA1[i,k], 1)
					# else
						# fix(activation_altitude_old[i,k], 0)
						# fix(activation_altitudeA1[i,k], 0)
					# end
						fix(activation_altitude_old[i,k], A)
						fix(activation_altitudeA1[i,k], A)
				end
			end

			global A = fill(0.0, NUM_SATELLITE, NUM_PIXEL, NUM_TIME)
			global D = fill(0.0, NUM_SATELLITE, NUM_TIME, NUM_PIXEL)
			global theta_max = fill(0.0, NUM_SATELLITE)

			global A1 = fill(0.0, NUM_SATELLITE, NUM_PIXEL, NUM_TIME)
			global D1 = fill(0.0, NUM_SATELLITE, NUM_TIME, NUM_PIXEL, iter_max)

			global Alt = fill(1.0, NUM_SATELLITE, length(altitudeSet))
			global Alt1 = fill(1.0, NUM_SATELLITE, length(altitudeSet))

			global LongitudeSat_final = fill(long[1], NUM_SATELLITE, NUM_TIME)
			global LatitudeSat_final = fill(lat[1], NUM_SATELLITE, NUM_TIME)
			global Thetamax_final = fill(0.0, NUM_SATELLITE)

			for i in 1:NUM_SATELLITE
				for p in 1:NUM_TIME
					fix(LongitudeSat1_old[i,p], LongitudeSat_final[i,p])
					fix(LatitudeSat1_old[i,p], LatitudeSat_final[i,p])
				end
			end

			global iteration = 0
			global value_nlp_min = 10e7
			fix(obj_milp, 10e7)

			while (value_nlp>=10e-2 || value_milp>=10e-2) && flag==1

				global iteration += 1
				println("iteration: ", iteration)

				if iteration >= iter_max
					break
				end

				if iteration < iter_max 

					for i in 1:NUM_SATELLITE
						for p in 1:NUM_TIME
							for j in 1:NUM_PIXEL
								fix(Xi1_old[i,p,j], D[i,p,j])
							end
						end
					end

					for i in 1:NUM_SATELLITE
						for k in 1:length(altitudeSet)
							fix(activation_altitudeA1[i,k], Alt[i,k])
						end
					end

					@constraint(model1, sum(abs3[i,p,j,iteration] for p=1:NUM_TIME for j=1:NUM_PIXEL for i=1:NUM_SATELLITE) >= 0.5)
					@constraint(model1, [i=1:NUM_SATELLITE, p=1:NUM_TIME, j=1:NUM_PIXEL], abs3[i,p,j,iteration] >= Xi1[i,p,j] - D1[i,p,j,iteration])
					@constraint(model1, [i=1:NUM_SATELLITE, p=1:NUM_TIME, j=1:NUM_PIXEL], abs3[i,p,j,iteration] >= D1[i,p,j,iteration] - Xi1[i,p,j])
					@constraint(model1, [i=1:NUM_SATELLITE, p=1:NUM_TIME, j=1:NUM_PIXEL], zz[i,p,j,iteration] => {abs3[i,p,j,iteration] <= Xi1[i,p,j] - D1[i,p,j,iteration]})
					@constraint(model1, [i=1:NUM_SATELLITE, p=1:NUM_TIME, j=1:NUM_PIXEL], !zz[i,p,j,iteration] => {abs3[i,p,j,iteration] <= D1[i,p,j,iteration] - Xi1[i,p,j]})
					@constraint(model1, [i=1:NUM_SATELLITE, p=1:NUM_TIME, j=1:NUM_PIXEL], !zz[i,p,j,iteration] => {Xi1[i,p,j] - D1[i,p,j,iteration] <= D1[i,p,j,iteration] - Xi1[i,p,j]})
					@constraint(model1, [i=1:NUM_SATELLITE, p=1:NUM_TIME, j=1:NUM_PIXEL], zz[i,p,j,iteration] => {D1[i,p,j,iteration] - Xi1[i,p,j] <= Xi1[i,p,j] - D1[i,p,j,iteration]})

					print("solving unfixed milp...")
					optimize!(model1)
					print(termination_status(model1))

					cpu_time = cpu_time + solve_time(model1)

					if (string(termination_status(model1)) == "INFEASIBLE" && iteration==1)
						break
					end

					if (string(termination_status(model1)) == "ALMOST_LOCALLY_SOLVED" || string(termination_status(model1)) == "LOCALLY_SOLVED" || string(termination_status(model1)) == "OPTIMAL")
						value_milp = objective_value(model1)
						print(": ")
						println(value_milp)

						for it in 1:iteration
							if D1[:,:,:,it] == value.(model1[:Xi1])
								println("ERROR: ", it)
							end
						end

						global D1[:,:,:,iteration+1] = value.(model1[:Xi1])
						global Alt1 = value.(model1[:activation_altitude1])
						global abs3a = value.(model1[:abs3])

						for i in 1:NUM_SATELLITE
							for p in 1:NUM_TIME
								for j in 1:NUM_PIXEL
									fix(Xi_old[i,p,j], D1[i,p,j,iteration+1])
								end
							end
						end

						for i in 1:NUM_SATELLITE
							global theta_max[i] = value.(model1[:Thetamax1])[i]
							for k in 1:length(altitudeSet)
								fix(activation_altitude_old[i,k], Alt1[i,k])
							end
						end

						for i in 1:NUM_SATELLITE
							for p in 1:NUM_TIME
								fix(LongitudeSat_old[i,p], value.(model1[:LongitudeSat1])[i,p])
								fix(LatitudeSat_old[i,p], value.(model1[:LatitudeSat1])[i,p])
							end
						end

						global obj_milp_current = 0
						for i in 1:NUM_SATELLITE
							for p=1:NUM_TIME
								for j=1:NUM_PIXEL
									global obj_milp_current += value.(model1[:abs])[i,p,j]
								end
							end
						end
						if obj_milp_current <= value.(model1[:obj_milp])
							fix(obj_milp, obj_milp_current)
						end
						println("obj_milp: ", obj_milp_current)

					else
						break
					end
				end

				if (iteration == 1)
					objective = 10e7
					fix(nlp_obj,objective)
				end

				termination = "empty"
				objective_min = 10e7
				for iter in 1:5

					for i in 1:NUM_SATELLITE
						for p in 1:NUM_TIME
							for j in 1:NUM_PIXEL
								# set_start_value(Xi[i,p,j], 1)
								set_start_value(Xi[i,p,j], rand(0:1))
							end
						end
					end

					for i in 1:NUM_SATELLITE
						for k in 1:length(altitudeSet)
							# if k == 1
								# set_start_value(activation_altitude[i,k], 1)
							# else
								# set_start_value(activation_altitude[i,k], 0)
							# end
							set_start_value(activation_altitude[i,k], rand(0:1))
						end
					end


					# for i=1:NUM_SATELLITE
						# start_variable_inc = rand(Uniform(0,1))
						# coin_flip =  rand(0:1)*2 - 1
						# set_start_value(sin_inclinaison[i], start_variable_inc)
						# set_start_value(cos_inclinaison[i], coin_flip*sqrt(1-start_variable_inc^2))

						# coin_flip =  rand(0:1)*2 - 1
						# start_variable_noeud = rand(Uniform(-1,1))
						# set_start_value(sin_noeudAscendant[i], start_variable_noeud)
						# set_start_value(cos_noeudAscendant[i], coin_flip*sqrt(1-start_variable_noeud^2))

						# coin_flip =  rand(0:1)*2 - 1
						# start_variable_mean = rand(Uniform(-1,1))
						# set_start_value(sin_meanAnomaly[i], start_variable_mean)
						# set_start_value(cos_meanAnomaly[i], coin_flip*sqrt(1-start_variable_mean^2))
						# set_start_value(altitude[i], RAYON+altitudeSet[1])
					# end
					
					if (iter > 1 || iteration == 1)
						for i=1:NUM_SATELLITE
							start_variable_inc = rand(Uniform(0,1))
							set_start_value(sin_inclinaison[i], start_variable_inc)
							set_start_value(cos_inclinaison[i], sqrt(1-start_variable_inc^2))
							start_variable_noeud = rand(Uniform(-1,1))
							set_start_value(sin_noeudAscendant[i], start_variable_noeud)
							set_start_value(cos_noeudAscendant[i], sqrt(1-start_variable_noeud^2))
							start_variable_mean = rand(Uniform(-1,1))
							set_start_value(sin_meanAnomaly[i], start_variable_mean)
							set_start_value(cos_meanAnomaly[i], sqrt(1-start_variable_mean^2))
							set_start_value(altitude[i], RAYON+altitudeSet[1])
						end
					end

					print("solving unfixed nlp...")
					optimize!(model)
					value_nlp = objective_value(model)
					termination = string(termination_status(model))
					print(termination)
					print(": ")
					println(value_nlp)
					cpu_time = cpu_time + solve_time(model)
					objective_min = 10e7
					global Inclinaison = asin.(round.(value.(model[:sin_inclinaison]);digits=4))
					global RAAN = asin.(round.(value.(model[:sin_noeudAscendant]);digits=4))
					global MeanAnomaly = asin.(round.(value.(model[:sin_meanAnomaly]);digits=4))
					global SemiAxis = (round.(value.(model[:altitude]);digits=4))
					if ((string(termination_status(model)) == "ALMOST_LOCALLY_SOLVED" || string(termination_status(model)) == "LOCALLY_SOLVED" || string(termination_status(model)) == "OPTIMAL") && value_nlp <= value_nlp_min)
						objective = 0
						for i=1:NUM_SATELLITE
							for j=1:NUM_PIXEL
								for p=1:NUM_TIME
									objective += (value.(model[:Xi])[i,p,j]-value.(model[:Xi_old])[i,p,j])^2
								end
							end
						end
						for i=1:NUM_SATELLITE
							for j=1:length(altitudeSet)
								objective += (value.(model[:activation_altitude])[i,j]-value.(model[:activation_altitude_old])[i,j])^2
							end
						end
						if objective <= objective_min
							D = value.(model[:Xi])
							Alt = value.(model[:activation_altitude])
							LongitudeSat_final = value.(model[:LongitudeSat])
							LatitudeSat_final = value.(model[:LatitudeSat])
							Thetamax_final = value.(model[:Thetamax])
							objective_min = objective
						end
						
						global observability = fill(0.0, 1:NUM_PIXEL, 1:PERIOD)
						global observability_pix = fill(1.0, 1:NUM_PIXEL)
						global observability_tot = 1.0
						println("minimum: ", minimum(D))
						# println(value.(model[:L5]))

						# println("observed period(s):")
						for j in 1:NUM_PIXEL
							for k in 1:PERIOD
								for p in (k-1)*NUM_PERIOD+1:(NUM_PERIOD*k)
									for i in 1:NUM_SATELLITE
										if((Thetamax_final[i] >= LatitudeSat_final[i,p]-lat[j]) && (Thetamax_final[i] >= -LatitudeSat_final[i,p]+lat[j]) 
										&& (Thetamax_final[i] >= (LongitudeSat_final[i,p]-long[j])*cos(lat[j])) && (Thetamax_final[i] >= (-LongitudeSat_final[i,p]+long[j])*cos(lat[j])))
											global observability[j,k] += 1
											# println(i, " ", j ," ", p, ": ", D[i,p,j])
										end
									end
								end
							end
						end

						for j in 1:NUM_PIXEL
							for k in 1:PERIOD
								if(observability[j,k] < 1.0)
									observability_pix[j] = 0.0
									break
								end
							end
						end

						for j in 1:NUM_PIXEL
							if observability_pix[j] < 1.0
								observability_tot = 0.0
								break
							end
						end

						if(observability_tot == 1.0)
							flag = 0
							break
						end
					end
				end
				
				println("objective_nlp: ", objective_min)
				fix(nlp_obj,objective)
				global observability = fill(0.0, 1:NUM_PIXEL, 1:PERIOD)
				global observability_pix = fill(1.0, 1:NUM_PIXEL)
				global observability_tot = 1.0

				# println("observed period(s):")
				for j in 1:NUM_PIXEL
					for k in 1:PERIOD
						for p in (k-1)*NUM_PERIOD+1:(NUM_PERIOD*k)
							for i in 1:NUM_SATELLITE
								if((Thetamax_final[i] >= LatitudeSat_final[i,p]-lat[j]) && (Thetamax_final[i] >= -LatitudeSat_final[i,p]+lat[j]) 
								&& (Thetamax_final[i] >= (LongitudeSat_final[i,p]-long[j])*cos(lat[j])) && (Thetamax_final[i] >= (-LongitudeSat_final[i,p]+long[j])*cos(lat[j])))
									global observability[j,k] += 1
									# println(i, " ", j ," ", p, ": ", D[i,p,j])
								end
							end
						end
					end
				end

				for j in 1:NUM_PIXEL
					for k in 1:PERIOD
						if(observability[j,k] < 1.0)
							observability_pix[j] = 0.0
							break
						end
					end
				end

				for j in 1:NUM_PIXEL
					if observability_pix[j] < 1.0
						observability_tot = 0.0
						break
					end
				end

				if(observability_tot == 1.0)
					flag = 0
					println(" status : ", string(termination_status(model)))
					println(" time CPU :", string(cpu_time))
					println(" alpha : ",string(alphaHalf*180/pi))
					println(" NUM SATELLITES : ", string(NUM_SATELLITE))
					println(" demi-axe : ", string((value.(model[:altitude]))))
					println(" sin_inclinaison : ",string(value.(model[:sin_inclinaison])))
					println(" cos_inclinaison : ",string(value.(model[:cos_inclinaison])))
					println(" sin_noeudAscendant : ",string(value.(model[:sin_noeudAscendant])))
					println(" cos_noeudAscendant : ",string(value.(model[:cos_noeudAscendant])))
					println(" sin_meanAnomaly : ",string(value.(model[:sin_meanAnomaly])))
					println(" cos_meanAnomaly : ",string(value.(model[:cos_meanAnomaly])))
					break
				end

				for i in 1:NUM_SATELLITE
					for p in 1:NUM_TIME
						fix(LongitudeSat1_old[i,p], LongitudeSat_final[i,p])
						fix(LatitudeSat1_old[i,p], LatitudeSat_final[i,p])
					end
				end
			end
		end
	end
	return Inclinaison,RAAN,MeanAnomaly,SemiAxis,NUM_SATELLITE
end