# Script to run the MIO formulation from Bertsimas and Mišić for n = 500, K = 1000, with and without cardinality constraints.
# This script will only solve the MIO directly. A time limit of 1 hour is imposed on each solve. 
# This script produces part of the results needed for Table 4 of the paper.
#
# The script takes command line arguments to determine which instances from MIOexpdata1 to do. To run the first 20 instances
# (as per the paper), type 
#
# > julia MIOexpdata1_BMonly_n500_K1000.jl 1 20
#
# into the command line.



using JuMP, Gurobi, Stats, Distributions 

include("../optimalPLD_code/optimalPLD_createPath.jl")
include("../optimalPLD_code/optimalPLD_solveSAO_constraints.jl")

srand(1204)
machine = 2;

n_vec = [500]
numPerm_vec = [1000]

numReps = 100;

counter = 0;

LR = parse(ARGS[1]);
UR = parse(ARGS[2]);

expdirpath = optimalPLD_createPath(machine, "exp_dir", "MIOexpdata1_BMonly_reps$(LR)to$(UR)_n500_K1000");
mkdir(expdirpath)

b_vec = [5; 10];

MIP_timeLimit = 60 #3600
isRelaxation = false;





data_prefix = "MIOexpdata1"

outcsvfilepath = string(expdirpath, "out.csv");
outcsvhandle = open(outcsvfilepath, "w")
row_string = string("n", ",",
					"numPermutations", ",",
					"rep", ",",
					"b_feas", ",",
					"sumxi", ",",
					"MIP_timeLimit", ",",
					"BM_expRevenue", ",",
					"BM_upperBound", ",",
					"BM_solutiontime", ",",
					"BM_MIPgap") #, ",",

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

		# if (n < 100)
		# 	continue;
		# end

		

		for rep in rep_range
			# if ( (n == 500) & (numPermutations == 500) & (rep < 11))
			# 	continue;
			# end

			revenues = revenues_mat[rep,:][:]
			lambda = lambda_mat[rep,:] #[1:numPermutations]
			orderings = orderings_mat[numPermutations*(rep-1)+1:rep*numPermutations,:]

			#subdirpath = optimalPLD_createPath(machine, "exp_subdir", "pe1", "neq$(n)_nPeq$(numPermutations)")
			#mkdir(subdirpath)

			b_vec_actual = [b_vec; n]


			for b_feas in b_vec_actual

				# counter += 1;

				# if ( (counter >= 2) & (nind == 1) )
				# 	continue;
				# end

				A_feas = ones(1,n); 

				println("**** Solving n = $n, nP = $numPermutations, rep = $rep ****")
				
				counter = counter + 1;
				isRelaxation = false;
				x_initial = [];

				if (b_feas >= n)
					println("**** Solving Bertsimas and Misic MIO (unconstrained) ****")
         			x_val, y_val, BM_expRevenue, BM_upperBound, BM_solutiontime, BM_MIPgap = optimalPLD_solveSAO_constraints( n, revenues, numPermutations, lambda, orderings, MIP_timeLimit, isRelaxation, x_initial, [], [])
				else
					println("**** Solving Bertsimas and Misic MIO (constrained) ****")
					x_val, y_val, BM_expRevenue, BM_upperBound, BM_solutiontime, BM_MIPgap = optimalPLD_solveSAO_constraints( n, revenues, numPermutations, lambda, orderings, MIP_timeLimit, isRelaxation, x_initial, A_feas, b_feas)
				end

				sumxi = sum(x_val); 


				row_string = string(n, ",",
						numPermutations, ",",
						rep, ",",
						b_feas, ",",
						sumxi, ",",
						MIP_timeLimit, ",",
						BM_expRevenue, ",",
						BM_upperBound, ",",
						BM_solutiontime, ",",
						BM_MIPgap) #, ",",
						
				row_string = string(row_string, "\n");

				print(outcsvhandle, row_string);
				flush(outcsvhandle);
			end
		end
	end
end




