using JuMP, CPLEX
include("Parameters.jl")

function MilpModel(model,NUM_SATELLITE,NUM_PIXEL,inclinaisonSet,RAANSet,meanAnomalySet,lat,long,PERIOD,NUM_PERIOD)
	
	n_inclinaison = size(inclinaisonSet)[1]
	n_RAAN = size(RAANSet)[1]
	n_meanAnomaly = size(meanAnomalySet)[1]

	CoverageSatLat = fill(0.0,NUM_PIXEL,NUM_TIME,n_inclinaison,n_RAAN,n_meanAnomaly,length(altitudeSet))
	CoverageSatLong = fill(0.0,NUM_PIXEL,NUM_TIME,n_inclinaison,n_RAAN,n_meanAnomaly,length(altitudeSet))
	
	for j in 1:altitudeLength
		for s in 1:n_inclinaison
			for l in 1:n_RAAN
				for k in 1:n_meanAnomaly
					global lat_Sat,long_Sat = ProjSatPlotPOLARDominique(RAYON+altitudeSet[j], inclinaisonSet[s],
					RAANSet[l], meanAnomalySet[k])

					for p in 1:NUM_TIME
						for t in 1:NUM_PIXEL
							CoverageSatLat[t,p,s,l,k,j] = abs(lat[t] - lat_Sat[p])
							CoverageSatLong[t,p,s,l,k,j] = abs(long[t] - long_Sat[p])*cos(lat[t])
						end
					end
				end
			end
		end
	end
	
	# Altitude de l'orbite [Km]
	@variable(model, altitude[i=1:NUM_SATELLITE])
	@variable(model, activation_altitude[1:NUM_SATELLITE,1:altitudeLength], Bin)

	# Inclinaison de l'orbite
	@variable(model, inclinaison[i=1:NUM_SATELLITE], start = inclinaisonSet[i])
	@variable(model, activation_inclinaison[1:NUM_SATELLITE,1:n_inclinaison], Bin)

	# Noeud ascendant
	@variable(model, RAAN[i=1:NUM_SATELLITE], start = RAANSet[i])
	@variable(model, activation_RAAN[1:NUM_SATELLITE,1:n_RAAN], Bin)

	# Anomalie moyenne du plan orbit
	@variable(model, meanAnomaly[i=1:NUM_SATELLITE], start = meanAnomalySet[i])
	@variable(model, activation_meanAnomaly[1:NUM_SATELLITE,1:n_meanAnomaly], Bin)

	# Variables de linéarisation
	@variable(model, Xi[1:NUM_SATELLITE, 1:NUM_TIME, 1:nbCibles], Bin)
	@variable(model, zeta[1:NUM_PIXEL,1:PERIOD], Bin)

	# Fonction objectif : Minimize the number of active satellite (choice between the next 2:)
	# objective(model, Min, sum(zeta[i] for i=1:NUM_SATELLITE))
	# @objective(model, Min, 0)
	# @objective(model, Max, sum(zeta[j,k] for j in 1:NUM_PIXEL for k in 1:PERIOD))
	@objective(model, Max, 100*NUM_SATELLITE*NUM_PERIOD*sum(zeta[j,k] for j in 1:NUM_PIXEL for k in 1:PERIOD) + sum(Xi[i,p,j] for i in 1:NUM_SATELLITE for p in 1:NUM_TIME for j in 1:NUM_PIXEL))

	# Compute the geocentric distance.
	@constraint(model, [i=1:NUM_SATELLITE], sum(activation_altitude[i,j] for j = 1:altitudeLength)==1)
	@constraint(model, [i=1:NUM_SATELLITE], altitude[i] == RAYON + sum(altitudeSet[j]*activation_altitude[i,j] for j=1:length(altitudeSet)))
	# Verification sinus/ cosinus relationship for angle inclinaison
	@constraint(model, [i=1:NUM_SATELLITE], sum(activation_inclinaison[i,j] for j = 1:n_inclinaison)==1)
	@constraint(model, inc_act[i=1:NUM_SATELLITE], inclinaison[i] == sum(inclinaisonSet[j]*activation_inclinaison[i,j] for j=1:n_inclinaison))
	# Verification sinus/ cosinus relationship for angle RAAN
	@constraint(model, [i=1:NUM_SATELLITE], sum(activation_RAAN[i,j] for j = 1:n_RAAN)==1)
	@constraint(model, noed_act[i=1:NUM_SATELLITE], RAAN[i] == sum(RAANSet[j]*activation_RAAN[i,j] for j=1:n_RAAN))
	# Verification sinus/ cosinus relationship for angle Mean Anomaly
	@constraint(model, [i=1:NUM_SATELLITE], sum(activation_meanAnomaly[i,j] for j = 1:n_meanAnomaly)==1)
	@constraint(model, mean_act[i=1:NUM_SATELLITE], meanAnomaly[i] == sum(meanAnomalySet[j]*activation_meanAnomaly[i,j] for j=1:n_meanAnomaly))

	# Pour la linéarisation des variables en utilisant des variables binaires
	@variable(model, act4[i=1:NUM_SATELLITE, j=1:altitudeLength, k=1:n_meanAnomaly, l=1:n_RAAN, s=1:n_inclinaison], Bin)
	@constraint(model, [i=1:NUM_SATELLITE, j=1:altitudeLength, k=1:n_meanAnomaly, l=1:n_RAAN, s=1:n_inclinaison], act4[i,j,k,l,s] <= activation_altitude[i,j])
	@constraint(model, [i=1:NUM_SATELLITE, j=1:altitudeLength, k=1:n_meanAnomaly, l=1:n_RAAN, s=1:n_inclinaison], act4[i,j,k,l,s] <= activation_meanAnomaly[i,k])
	@constraint(model, [i=1:NUM_SATELLITE, j=1:altitudeLength, k=1:n_meanAnomaly, l=1:n_RAAN, s=1:n_inclinaison], act4[i,j,k,l,s] <= activation_inclinaison[i,s])
	@constraint(model, [i=1:NUM_SATELLITE, j=1:altitudeLength, k=1:n_meanAnomaly, l=1:n_RAAN, s=1:n_inclinaison], act4[i,j,k,l,s] <= activation_RAAN[i,l])
	@constraint(model, [i=1:NUM_SATELLITE, j=1:altitudeLength, k=1:n_meanAnomaly, l=1:n_RAAN, s=1:n_inclinaison], act4[i,j,k,l,s] >=
	activation_altitude[i,j] + activation_meanAnomaly[i,k] + activation_inclinaison[i,s] + activation_RAAN[i,l] - 3)

	# Thetamax permet de déterminer la surface couverte par le satellite à partir du demi angle d'ouverture des capteurs, de l'altitude et de l'excentricité
	@variable(model, Thetamax[1:NUM_SATELLITE])
	@constraint(model, [i=1:NUM_SATELLITE], Thetamax[i] == sum(Theta[j]*activation_altitude[i,j] for j=1:length(altitudeSet)))
			
	@constraint(model, [p=1:NUM_TIME,j=1:NUM_PIXEL,i=1:NUM_SATELLITE], Xi[i,p,j] => {Thetamax[i] >= (sum(act4[i,a,k,l,s]*CoverageSatLat[j,p,s,l,k,a]
	for s=1:n_inclinaison for l=1:n_RAAN for k in 1:n_meanAnomaly for a=1:altitudeLength))})
	# @constraint(model, [p=1:NUM_TIME,j=1:NUM_PIXEL,i=1:NUM_SATELLITE], Xi[i,p,j] => {Thetamax[i] >= (sum(act4[i,a,k,l,s]*(-1)*CoverageSatLat[j,p,s,l,k,a]
	# for s=1:n_inclinaison for l=1:n_RAAN for k in 1:n_meanAnomaly for a=1:altitudeLength))})
	@constraint(model, [p=1:NUM_TIME,j=1:NUM_PIXEL,i=1:NUM_SATELLITE], Xi[i,p,j] => {Thetamax[i] >= (sum(act4[i,a,k,l,s]*CoverageSatLong[j,p,s,l,k,a]
	for s=1:n_inclinaison for l=1:n_RAAN for k in 1:n_meanAnomaly for a=1:altitudeLength))})
	# @constraint(model, [p=1:NUM_TIME,j=1:NUM_PIXEL,i=1:NUM_SATELLITE], Xi[i,p,j] => {Thetamax[i] >= (sum(act4[i,a,k,l,s]*(-1)*CoverageSatLong[j,p,s,l,k,a]
	# for s=1:n_inclinaison for l=1:n_RAAN for k in 1:n_meanAnomaly for a=1:altitudeLength))})

	# Cut the simmetry of problem. To help branch-and-bound
	#@constraint(model, [i=1:NUM_SATELLITE-1], zeta[i] >= zeta[i+1])

	# Linearisation du modèle
	#@constraint(model, [j=1:NUM_PIXEL,k=1:PERIOD], sum(Xi[i,p,j] for p in (((k-1)*NUM_PERIOD)+1):(NUM_PERIOD*k) for i in 1:NUM_SATELLITE) >= 1)
	@constraint(model, [j=1:NUM_PIXEL,k=1:PERIOD], zeta[j,k] => {sum(Xi[i,p,j] for p in (((k-1)*NUM_PERIOD)+1):(NUM_PERIOD*k) for i in 1:NUM_SATELLITE) >= 1})

	#delete(model, noed_act)
	optimize!(model)
	status = termination_status(model)	
	return model
end	

function printing(Output,flag_infeasible,model,status,solvedTime,dt,alphaHalf,Theta_min,value, NUM_SATELLITE)
	write(Output," status : ", string(status))
	write(Output,";  time CPU :", string(solvedTime))
	write(Output,"\n  dt : ",string(dt))
	write(Output,";  alpha : ",string(alphaHalf*180/pi))
	write(Output,";  theta_min : ",string(Theta_min*180/pi))
	write(Output,";  theta : ",string(value.(model[:Thetamax])*(180/pi)))		
	write(Output,"\n NUM_SATELLITE = ", string(NUM_SATELLITE))
	write(Output,"\n  demi-axe = ", string((value.(model[:altitude]))))
	write(Output,"\n  inclinaison = ",string(value.(model[:inclinaison])))
	write(Output,"\n  RAAN = ",string(value.(model[:RAAN])))
	write(Output,"\n  meanAnomaly = ",string(value.(model[:meanAnomaly])))
	# write(Output,"\n  index_semiAxis =  ",string(findall(a->a==1,value.(model[:activation_altitude])[:,:][1])))
	# write(Output,"\n  k =  ",string(findall(a->a==1,value.(model[:activation_inclinaison])[:,:][1])))
	# write(Output,"\n  l =  ",string(findall(a->a==1,value.(model[:activation_RAAN])[:,:][1])))
	# write(Output,"\n  s =  ",string(findall(a->a==1,value.(model[:activation_meanAnomaly])[:,:][1])))
	write(Output,"\n"," ")
	#println(flag_infeasible)
end

function PrintingNotFeasible(Xi_new,NUM_SATELLITE,NUM_PIXEL,iter,NUM_PERIOD,PERIOD)
	flag_infeasible = fill(0,NUM_PIXEL)
	for j=1:NUM_PIXEL
		for k=1:PERIOD
			if(sum(Xi_new[i,p,j] for p in (((k-1)*NUM_PERIOD)+1):(NUM_PERIOD*k) for i in 1:NUM_SATELLITE) < 1)
				println("SOLUTION NOT FEASIBLE FOR TARGET ", j)
				flag_infeasible[j] = 1
				break
			end
		end
	end
	return flag_infeasible
end

function DetectingCycle(flag_infeasible,iter,Xi_new,Xi_old,NUM_PIXEL,NUM_SATELLITE,KepParam,UppBound)
	global cycle = false
	for j in 1:NUM_PIXEL
		if flag_infeasible[j] == 1 && cycle == false
			for k in 1:iter
				if isequal(Xi_new[:,:,j],Xi_old[:,:,j,k])
					println("CYCLE DETECTED")
					KepParam = UppBound*pi*rand(NUM_SATELLITE)
					global cycle = true
					break
				end
			end
		end
	end
	return KepParam
end

function PrintingObserving(model,Xi_new,dt,NUM_PIXEL,NUM_SATELLITE)		
	for i=1:NUM_SATELLITE
		for j=1:NUM_PIXEL
			for p=1:NUM_TIME
				if(Xi_new[i,p,j][1] > 0)
					println("observing ",j," at time: ",p*dt)
				end
			end
		end
	end
end