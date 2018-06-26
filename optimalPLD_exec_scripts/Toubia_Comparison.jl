# Script to run the four formulations -- Bertsimas and Mišić, Belloni et al., McBride and Zufryden, and the Utility formulation -- in their integer and relaxed forms, 
# for the Toubia data set.
# This script produces the results needed for Tables EC.1 and EC.2 in Section EC.3.1 of the electronic companion of 
# Bertsimas and Mišić.
# 
# For each value of n and numPermutations, the script randomly selects a set of n products out of the universe
# of 3584, and numPermutations customers out of the set of 330, and solves the corresponding PLD problem with 
# those parameters.
# The script takes two command line arguments that specify the "rep range". To run 20 repetitions for 
# each n and numPermutations value (as per the paper), type 
#
# > julia Toubia_Comparison.jl 1 20
#
# in the terminal.

using JuMP, Gurobi, Stats, Distributions #, HDF5, JLD

include("../optimalPLD_code/optimalPLD_createPath.jl")
include("../optimalPLD_code/optimalPLD_solveSAO_constraints.jl")
include("../optimalPLD_code/optimalPLD_solveSAO_BFSS_constraints.jl")
include("../optimalPLD_code/optimalPLD_solveSAO_McBZ_constraints.jl")
include("../optimalPLD_code/optimalPLD_solveSAO_Utility1_constraints.jl")


srand(1204)
machine = 2;

n_vec = [10, 20, 50, 100]
numPerm_vec = [10, 100] #, 330]

numReps = 100;

counter = 0;

LR = parse(ARGS[1]);
UR = parse(ARGS[2]);

expdirpath = optimalPLD_createPath(machine, "exp_dir", "Toubia_Comparison_reps$(LR)to$(UR)");
mkdir(expdirpath)

b_vec = [5; 10; 20; 100];


MIP_timeLimit = 1*3600;
isRelaxation = false;


data_prefix = "Toubia2003"

outcsvfilepath = string(expdirpath, "out.csv");
outcsvhandle = open(outcsvfilepath, "w")
row_string = string("n", ",",
					"numPermutations", ",",
					"rep", ",",
					"b_RHS", ",",
					"sumxi", ",",
					"BM_expRevenue", ",",
					"BM_upperBound", ",",
					"BM_solutiontime", ",",
					"BM_MIPgap", ",",
					"BM_LP_expRevenue", ",",
					"BM_LP_upperBound", ",",
					"BM_LP_solutiontime", ",",
					"BM_LP_MIPgap", ",",
					"BM_LP_isIntegral", ",",
					"BFSS_expRevenue", ",",
					"BFSS_upperBound", ",",
					"BFSS_solutiontime", ",",
					"BFSS_MIPgap", ",",
					"BFSS_LP_expRevenue", ",",
					"BFSS_LP_upperBound", ",",
					"BFSS_LP_solutiontime", ",",
					"BFSS_LP_MIPgap", ",",
					"BFSS_LP_isIntegral", ",",
					"Utility1_expRevenue", ",",
					"Utility1_upperBound", ",",
					"Utility1_solutiontime", ",",
					"Utility1_MIPgap", ",",
					"Utility1_LP_expRevenue", ",",
					"Utility1_LP_upperBound", ",",
					"Utility1_LP_solutiontime", ",",
					"Utility1_LP_MIPgap", ",",
					"Utility1_LP_isIntegral", ",",
					"McBZ_expRevenue", ",",
					"McBZ_upperBound", ",",
					"McBZ_solutiontime", ",",
					"McBZ_MIPgap", ",",
					"McBZ_LP_expRevenue", ",",
					"McBZ_LP_upperBound", ",",
					"McBZ_LP_solutiontime", ",",
					"McBZ_LP_MIPgap", ",",
					"McBZ_LP_isIntegral" )
row_string = string(row_string, "\n");

print(outcsvhandle, row_string);
flush(outcsvhandle);


rep_range = LR:UR

counter = 0;
for nind = 1:length(n_vec)
	n = n_vec[nind]
	for nPind = 1:length(numPerm_vec)
		numPermutations = numPerm_vec[nPind];

		println("**** Reading from data directory...")
		data_name = string(data_prefix,"_neq3584_Keq330_v2")
		data_dir = optimalPLD_createPath(machine, "data_dir", data_name)
		revenues_mat = readcsv(string(data_dir, "revenues_mat.csv"))
		lambda_mat = readcsv(string(data_dir, "lambda_mat.csv"))
		orderings_mat = readcsv(string(data_dir, "orderings_mat.csv"))
		utilities_mat = readcsv(string(data_dir, "utilities_mat.csv"))

		orderings_mat = convert(Array{Int64}, orderings_mat);

		for rep in rep_range

			prod_universe = randperm(3584)[1:n];
			prod_universe = sort(prod_universe);
			push!(prod_universe, n+1);
			cust_universe = randperm(330)[1:numPermutations];
			cust_universe = sort(cust_universe);

			revenues = revenues_mat[1,:][prod_universe]
			lambda = 1/numPermutations *ones(numPermutations); #[1:numPermutations]
			orderings_raw = orderings_mat[cust_universe,:]

			orderings = zeros(Int64, numPermutations, n+1);
			for k=1:numPermutations
				p2r = sortperm(orderings_raw[k,:]);
				p2r = p2r[prod_universe];
				orderings[k,:] = sortperm(p2r);
			end

			utilities = utilities_mat[cust_universe, prod_universe];
			for k=1:numPermutations
				utilities[k,:] -= minimum(utilities[k,:]);
			end


			writecsv(string(expdirpath, "prod_universe_neq", n, "_Keq", numPermutations, "rep_", rep, ".csv"), prod_universe)
			writecsv(string(expdirpath, "cust_universe_neq", n, "_Keq", numPermutations, "rep_", rep, ".csv"), cust_universe)

			b_vec_2 = filter(i -> i < n, b_vec)
			push!(b_vec_2, n); 


			for b_feas_2 in b_vec_2

				
				b_feas = b_feas_2;
				if (b_feas_2 >= n)
					b_feas = Int64[];
				end
			
				A_feas = ones(1,n); 
				x_initial = [];

				println("**** Solving n = $n, nP = $numPermutations, rep = $rep ****")
				println("**** Solving Bertsimas and Misic MIO ****")
				counter = counter + 1;
				isRelaxation = false;
				x_val, y_val, BM_expRevenue, BM_upperBound, BM_solutiontime, BM_MIPgap = optimalPLD_solveSAO_constraints( n, revenues, numPermutations, lambda, orderings, MIP_timeLimit, isRelaxation, x_initial, A_feas, b_feas)

				sumxi = sum(x_val); 

				println("**** Solving Bertsimas and Misic LP ****")
				isRelaxation = true;
				BM_LP_x_val, BM_LP_y_val, BM_LP_expRevenue, BM_LP_upperBound, BM_LP_solutiontime, BM_LP_MIPgap = optimalPLD_solveSAO_constraints( n, revenues, numPermutations, lambda, orderings, MIP_timeLimit, isRelaxation, x_initial, A_feas, b_feas )
				BM_LP_isIntegral = isinteger(BM_LP_x_val) ? 1 : 0; 



				println("**** Solving Beloni MIO ****")
				isRelaxation = false;
				x_val, y_val, BFSS_expRevenue, BFSS_upperBound, BFSS_solutiontime, BFSS_MIPgap = optimalPLD_solveSAO_BFSS_constraints( n, revenues, numPermutations, lambda, orderings, MIP_timeLimit, isRelaxation, x_initial, A_feas, b_feas )

				println("**** Solving Beloni LP ****")
				isRelaxation = true;
				BFSS_LP_x_val, BFSS_y_val, BFSS_LP_expRevenue, BFSS_LP_upperBound, BFSS_LP_solutiontime, BFSS_LP_MIPgap = optimalPLD_solveSAO_BFSS_constraints( n, revenues, numPermutations, lambda, orderings, MIP_timeLimit, isRelaxation, x_initial, A_feas, b_feas )
				BFSS_LP_isIntegral = isinteger(BFSS_LP_x_val) ? 1 : 0; 


				println("**** Solving McBride Zufryden MIO ****")
				isRelaxation = false;
				x_val, y_val, McBZ_expRevenue, McBZ_upperBound, McBZ_solutiontime, McBZ_MIPgap = optimalPLD_solveSAO_McBZ_constraints( n, revenues, numPermutations, lambda, utilities, MIP_timeLimit, isRelaxation, x_initial, A_feas, b_feas )
				

				println("**** Solving McBride Zufryden LP ****")
				isRelaxation = true;
				McBZ_LP_x_val, McBZ_LP_y_val, McBZ_LP_expRevenue, McBZ_LP_upperBound, McBZ_LP_solutiontime, McBZ_LP_MIPgap = optimalPLD_solveSAO_McBZ_constraints( n, revenues, numPermutations, lambda, utilities, MIP_timeLimit, isRelaxation, x_initial, A_feas, b_feas )
				McBZ_LP_isIntegral = isinteger(McBZ_LP_x_val) ? 1 : 0; 


				println("**** Solving Utility1 MIO ****")
				isRelaxation = false;
				x_val, y_val, Utility1_expRevenue, Utility1_upperBound, Utility1_solutiontime, Utility1_MIPgap = optimalPLD_solveSAO_Utility1_constraints( n, revenues, numPermutations, lambda, utilities, MIP_timeLimit, isRelaxation, x_initial, A_feas, b_feas )			

				println("**** Solving Utility1 LP ****")
				isRelaxation = true;
				Utility1_LP_x_val, Utility1_LP_y_val, Utility1_LP_expRevenue, Utility1_LP_upperBound, Utility1_LP_solutiontime, Utility1_LP_MIPgap = optimalPLD_solveSAO_Utility1_constraints( n, revenues, numPermutations, lambda, utilities, MIP_timeLimit, isRelaxation, x_initial, A_feas, b_feas )			
				Utility1_LP_isIntegral = isinteger(Utility1_LP_x_val) ? 1 : 0; 




				row_string = string(n, ",",
						numPermutations, ",",
						rep, ",",
						b_feas_2, ",",
						sumxi, ",",
						BM_expRevenue, ",",
						BM_upperBound, ",",
						BM_solutiontime, ",",
						BM_MIPgap, ",",
						BM_LP_expRevenue, ",",
						BM_LP_upperBound, ",",
						BM_LP_solutiontime, ",",
						BM_LP_MIPgap, ",",
						BM_LP_isIntegral, ",",
						BFSS_expRevenue, ",",
						BFSS_upperBound, ",",
						BFSS_solutiontime, ",",
						BFSS_MIPgap, ",",
						BFSS_LP_expRevenue, ",",
						BFSS_LP_upperBound, ",",
						BFSS_LP_solutiontime, ",",
						BFSS_LP_MIPgap, ",",
						BFSS_LP_isIntegral, ",",
						Utility1_expRevenue, ",",
						Utility1_upperBound, ",",
						Utility1_solutiontime, ",",
						Utility1_MIPgap, ",",
						Utility1_LP_expRevenue, ",",
						Utility1_LP_upperBound, ",",
						Utility1_LP_solutiontime, ",",
						Utility1_LP_MIPgap, ",",
						Utility1_LP_isIntegral, ",",
						McBZ_expRevenue, ",",
						McBZ_upperBound, ",",
						McBZ_solutiontime, ",",
						McBZ_MIPgap, ",",
						McBZ_LP_expRevenue, ",",
						McBZ_LP_upperBound, ",",
						McBZ_LP_solutiontime, ",",
						McBZ_LP_MIPgap, ",",
						McBZ_LP_isIntegral )
				row_string = string(row_string, "\n");

				print(outcsvhandle, row_string);
				flush(outcsvhandle);
			end
		end
	end
end




