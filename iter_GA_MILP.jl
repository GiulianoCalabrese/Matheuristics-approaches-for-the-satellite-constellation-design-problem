push!(LOAD_PATH, "$(pwd())/SCDP")
using FileIO 
using Tables
using StaticArrays
using JLD2
using SCDP
include("MILPmodel.jl")

global name_file = "GA-MILP_"*string(TEMPS_SIMU/(3600*inst.n_0))*"h_Nb_Cibles"*string(nbCibles)*".txt"
println(name_file)
open(name_file, "w") do Output
	global max_iter = 1*100
	
	discretizaionSet2pi = Discretization(numbOfDiscretize,0.0,2*pi)
	discretizaionSet1pi = Discretization(numbOfDiscretize,0.0,pi)
	
	global NUM_PIXEL = inst.nbCibles
	global NUM_SATELLITE = 1	
	global Xi_old = fill(1.0,1:NUM_SATELLITE,1:NUM_TIME,1:NUM_PIXEL,1:max_iter)
	global Xi_new = fill(0.0,1:NUM_SATELLITE,1:NUM_TIME,1:NUM_PIXEL)
	lat = []
	long = []
	write(Output," coordonnées cibles : ")
	for index in 1:NUM_PIXEL
		append!(lat,inst.coordsCibles[index][1])
		append!(long,inst.coordsCibles[index][2])
	end
	write(Output," lat : ",  string(lat*(180/pi)))
	write(Output," long : ",  string(long*(180/pi)))
	write(Output,"\n"," ")	
	status = "nothing"
	println("NUM_SATELLITE = "*string(NUM_SATELLITE))
	
	j = 0
	iter = 1
	flag_infeasible = PrintingNotFeasible(Xi_new,NUM_SATELLITE,NUM_PIXEL,iter,inst.nbInstants,inst.n_0)
	
	for pMut in probasMut
	   for pModif in probasModifNbSats 
			 global popul = resoudre_SCDP(inst, nbIndiv, nbMaxGen, ϵ; probaMut = pMut, 
								probaModifNbSat = pModif, 
								fDiag = (;kwargs...) -> fDiag(nbIndiv; kwargs...)) 
	   end
	end

	indiv = popul[1]
	conste = consteVar2Conste(indiv.x.p[1:4indiv.x.nbBlocs[]], inst)
	global Inclinaison = conste.is
	global RAAN = conste.Ωs
	global MeanAnomaly = conste.Ms
	println("GA__NUM_SATELLITE = "*string(conste.nbSats))
	solvedTime = 0

	if conste.nbSats == NUM_SATELLITE
		println("GA_SOLUTION")
		printing(Output,flag_infeasible,model,status,solvedTime,dt,alphaHalf,Theta_min,value,NUM_SATELLITE)
	end
	
	while (NUM_SATELLITE <= nbMaxSat) && (sum(flag_infeasible[p] for p in 1:NUM_PIXEL) != 0)			
		while (iter <= max_iter) && (sum(flag_infeasible[p] for p in 1:NUM_PIXEL) != 0)		
			
			println(iter)
			global model = Model(CPLEX.Optimizer)
			set_optimizer_attribute(model, "CPX_PARAM_TILIM", 100)
			set_optimizer_attribute(model, "CPX_PARAM_SCRIND", 0)
			
			if j == 0
				inclinaisonSet = sort(Inclinaison)
				RAANSet = sort(RAAN)
				meanAnomalySet = union(MeanAnomaly,discretizaionSet2pi)			
				global model = MilpModel(model,NUM_SATELLITE,NUM_PIXEL,inclinaisonSet,RAANSet,meanAnomalySet,lat,long,inst.n_0,inst.nbInstants)
				MeanAnomaly = value.(model[:meanAnomaly])		
			elseif j == 1
				inclinaisonSet = union(Inclinaison,discretizaionSet1pi)
				RAANSet = sort(RAAN)
				meanAnomalySet = sort(MeanAnomaly)	
				global model = MilpModel(model,NUM_SATELLITE,NUM_PIXEL,inclinaisonSet,RAANSet,meanAnomalySet,lat,long,inst.n_0,inst.nbInstants)
				Inclinaison = value.(model[:inclinaison])
			else
				inclinaisonSet = sort(Inclinaison)
				RAANSet = union(RAAN,discretizaionSet2pi)
				meanAnomalySet = sort(MeanAnomaly)	
				global model = MilpModel(model,NUM_SATELLITE,NUM_PIXEL,inclinaisonSet,RAANSet,meanAnomalySet,lat,long,inst.n_0,inst.nbInstants)
				RAAN = value.(model[:RAAN])			
			end	
			solvedTime += solve_time(model)
			Xi_old[:,:,:,iter] = Xi_new
			Xi_new = value.(model[:Xi])		
			flag_infeasible = PrintingNotFeasible(Xi_new,NUM_SATELLITE,NUM_PIXEL,iter,inst.nbInstants,inst.n_0)
			# PrintingObserving(Xi_new,dt,NUM_PIXEL,NUM_SATELLITE)
			if j == 0
				println("MILP_MEAN_ANOMALY")
				MeanAnomaly = sort(value.(model[:meanAnomaly]))
				MeanAnomaly = DetectingCycle(flag_infeasible,iter,Xi_new,Xi_old,NUM_PIXEL,NUM_SATELLITE,MeanAnomaly,2)
			elseif j == 1
				println("MILP_INCLINAISON")
				Inclinaison = sort(value.(model[:inclinaison]))
				Inclinaison = DetectingCycle(flag_infeasible,iter,Xi_new,Xi_old,NUM_PIXEL,NUM_SATELLITE,Inclinaison,1)
			else
				println("MILP_NOEUD_ASCENDANT")
				RAAN = sort(value.(model[:RAAN]))
				RAAN = DetectingCycle(flag_infeasible,iter,Xi_new,Xi_old,NUM_PIXEL,NUM_SATELLITE,RAAN,2)
			end				
			j = (j+1)%3
			iter += 1						
		end	
		
		if (sum(flag_infeasible[p] for p in 1:NUM_PIXEL) != 0)	
			println("NO SOLUTION")
			j = 0
			iter = 1
			NUM_SATELLITE = NUM_SATELLITE + 1
			global Inclinaison = conste.is
			global RAAN = conste.Ωs
			global MeanAnomaly = conste.Ms
			Xi_new = fill(0.0,1:NUM_SATELLITE,1:NUM_TIME,1:NUM_PIXEL)
			Xi_old = fill(1.0,1:NUM_SATELLITE,1:NUM_TIME,1:NUM_PIXEL,1:max_iter)
			flag_infeasible = PrintingNotFeasible(Xi_new,NUM_SATELLITE,NUM_PIXEL,iter,inst.nbInstants,inst.n_0)
			empty!(model)
		else 
			status = termination_status(model)
			println("RESOLVED")
			printing(Output,flag_infeasible,model,status,solvedTime,dt,alphaHalf,Theta_min,value,NUM_SATELLITE)
		end
	end
end