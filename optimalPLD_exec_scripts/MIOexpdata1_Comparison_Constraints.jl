# Script to run the four formulations -- Bertsimas and Mišić, Belloni et al., McBride and Zufryden, and the Utility formulation -- 
# in their integer and relaxed forms, for synthetic data, with cardinality constraints (sum x_i <= b).
# This script produces part of the results needed for Tables 5 and 6 of the paper.
#
# The script takes two command line arguments that specify the "rep range". The MIOexpdata1 synthetic instances
# were generated in such a way that for each n and numPermutations combination, there are 100 instances; only the 
# first 20 are used in Bertsimas and Mišić. To run the first twenty, one would type:
#
# > julia MIOexpdata1_Comparison_Constraints.jl 1 20
#
# in the terminal. 

using JuMP, Gurobi, Stats, Distributions 

include("../optimalPLD_code/optimalPLD_createPath.jl")
include("../optimalPLD_code/optimalPLD_solveSAO_constraints.jl")
include("../optimalPLD_code/optimalPLD_solveSAO_BFSS_constraints.jl")
include("../optimalPLD_code/optimalPLD_solveSAO_McBZ_constraints.jl")
include("../optimalPLD_code/optimalPLD_solveSAO_Utility1_constraints.jl")


srand(1204)
machine = 2;

n_vec = [20, 50, 100]
numPerm_vec = [100,200,500,1000]

numReps = 100;

counter = 0;

LR = parse(ARGS[1]);
UR = parse(ARGS[2]);

expdirpath = optimalPLD_createPath(machine, "exp_dir", "MIOexpdata1_constrained_comparison_reps$(LR)to$(UR)");
mkdir(expdirpath)


b_vec = [5; 10]; # test sum x_i <= b for b = 5, 10.


MIP_timeLimit = 1*3600; # Impose a time limit of 1 hour on every MIO formulation
isRelaxation = false;



data_prefix = "MIOexpdata1"

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
		data_name = string(data_prefix,"_neq",n,"_Keq",numPermutations)
		data_dir = optimalPLD_createPath(machine, "data_dir", data_name)
		revenues_mat = readcsv(string(data_dir, "revenues_mat.csv"))
		lambda_mat = readcsv(string(data_dir, "lambda_mat.csv"))
		orderings_mat = readcsv(string(data_dir, "orderings_mat.csv"))

		orderings_mat = convert(Array{Int64}, orderings_mat);

		for rep in rep_range
			revenues = revenues_mat[rep,:][:]
			lambda = lambda_mat[rep,:] #[1:numPermutations]
			orderings = orderings_mat[numPermutations*(rep-1)+1:rep*numPermutations,:]

			for b_feas in b_vec


				if (b_feas >= n)
					continue;
				end
			
				A_feas = ones(1,n); 

				x_initial = [];

				println("**** Solving n = $n, nP = $numPermutations, rep = $rep ****")
				println("**** Solving Bertsimas and Misic MIO ****")
				counter = counter + 1;
				isRelaxation = false;
				x_val, y_val, BM_expRevenue, BM_upperBound, BM_solutiontime, BM_MIPgap = optimalPLD_solveSAO_constraints( n, revenues, numPermutations, lambda, orderings, MIP_timeLimit, isRelaxation, x_initial, A_feas, b_feas)

				sumxi = sum(x_val); 

				println("**** Solving Bertsimas and Misic LO ****")
				isRelaxation = true;
				BM_LP_x_val, BM_LP_y_val, BM_LP_expRevenue, BM_LP_upperBound, BM_LP_solutiontime, BM_LP_MIPgap = optimalPLD_solveSAO_constraints( n, revenues, numPermutations, lambda, orderings, MIP_timeLimit, isRelaxation, x_initial,  A_feas, b_feas )
				BM_LP_isIntegral = isinteger(BM_LP_x_val) ? 1 : 0; 


				println("**** Solving Belloni et al. MIO ****")
				isRelaxation = false;
				x_val, y_val, BFSS_expRevenue, BFSS_upperBound, BFSS_solutiontime, BFSS_MIPgap = optimalPLD_solveSAO_BFSS_constraints( n, revenues, numPermutations, lambda, orderings, MIP_timeLimit, isRelaxation, x_initial, A_feas, b_feas )

				println("**** Solving Belloni et al. LO ****")
				isRelaxation = true;
				BFSS_LP_x_val, BFSS_y_val, BFSS_LP_expRevenue, BFSS_LP_upperBound, BFSS_LP_solutiontime, BFSS_LP_MIPgap = optimalPLD_solveSAO_BFSS_constraints( n, revenues, numPermutations, lambda, orderings, MIP_timeLimit, isRelaxation, x_initial, A_feas, b_feas )
				BFSS_LP_isIntegral = isinteger(BFSS_LP_x_val) ? 1 : 0; 



				utilities = zeros(numPermutations, n+1)
				for k=1:numPermutations
					utilities[k,:] = n - (sortperm(orderings[k,:][:]) - 1)
				end

				println("**** Solving McBride & Zufryden MIO ****")
				isRelaxation = false;
				x_val, y_val, McBZ_expRevenue, McBZ_upperBound, McBZ_solutiontime, McBZ_MIPgap = optimalPLD_solveSAO_McBZ_constraints( n, revenues, numPermutations, lambda, utilities, MIP_timeLimit, isRelaxation, x_initial, A_feas, b_feas )
				

				println("**** Solving McBride Zufryden LO ****")
				isRelaxation = true;
				McBZ_LP_x_val, McBZ_LP_y_val, McBZ_LP_expRevenue, McBZ_LP_upperBound, McBZ_LP_solutiontime, McBZ_LP_MIPgap = optimalPLD_solveSAO_McBZ_constraints( n, revenues, numPermutations, lambda, utilities, MIP_timeLimit, isRelaxation, x_initial, A_feas, b_feas )
				McBZ_LP_isIntegral = isinteger(McBZ_LP_x_val) ? 1 : 0; 


				println("**** Solving Utility1 MIO ****")
				isRelaxation = false;
				x_val, y_val, Utility1_expRevenue, Utility1_upperBound, Utility1_solutiontime, Utility1_MIPgap = optimalPLD_solveSAO_Utility1_constraints( n, revenues, numPermutations, lambda, utilities, MIP_timeLimit, isRelaxation, x_initial, A_feas, b_feas )			

				println("**** Solving Utility1 LO ****")
				isRelaxation = true;
				Utility1_LP_x_val, Utility1_LP_y_val, Utility1_LP_expRevenue, Utility1_LP_upperBound, Utility1_LP_solutiontime, Utility1_LP_MIPgap = optimalPLD_solveSAO_Utility1_constraints( n, revenues, numPermutations, lambda, utilities, MIP_timeLimit, isRelaxation, x_initial, A_feas, b_feas )			
				Utility1_LP_isIntegral = isinteger(Utility1_LP_x_val) ? 1 : 0; 




			

				row_string = string(n, ",",
						numPermutations, ",",
						rep, ",",
						b_feas, ",",
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


close(outcsvhandle)

