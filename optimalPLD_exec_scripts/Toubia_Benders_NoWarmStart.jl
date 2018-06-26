# Script to run the cold-started Benders algorithm (LO relaxation phase and integer phase; divide and conquer is not executed).
# This script produces the numbers needed for Table EC.3 in Section EC.3.2 of the electronic companion of the paper.
#
# This script does not take any command line arguments. To run it,
# type 
#
# > julia Toubia_Benders_NoWarmStart.jl 
#
# in the terminal. 


using JuMP, Gurobi

include("../optimalPLD_code/optimalPLD_createPath.jl")
include("../optimalPLD_code/optimalPLD_solveSAO_Benders.jl")


include("../optimalPLD_code/optimalPLD_solveSAO_LP_Benders.jl")



srand(1683)
machine = 2;


n_vec = [3584]
numPerm_vec = [330]

numReps = 1;

counter = 0;

data_prefix = "toubia2003"

expdirpath = optimalPLD_createPath(machine, "exp_dir", "Toubia_Benders_NoWarmStart");
mkdir(expdirpath)


numInst = 8;
inst_vec = [1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8]
b_feas_vec = [1, -1, 2, -2, 3, -3, 4, -4, 5, -5, 10, -10, 20, -20, 50, -50]
n = 3584;
A_feas_mat = [ones(1,n); 
			 -ones(1,n);
			 ones(1,n); 
			 -ones(1,n);
			 ones(1,n); 
			 -ones(1,n);
			 ones(1,n); 
			 -ones(1,n);
			 ones(1,n); 
			 -ones(1,n);
			 ones(1,n); 
			 -ones(1,n);
			 ones(1,n); 
			 -ones(1,n);
			 ones(1,n); 
			 -ones(1,n)]


text_description_by_inst = ["sum x_i == 1", "sum x_i == 2", "sum x_i == 3", 
							"sum x_i == 4", "sum x_i == 5", "sum x_i == 10",
							"sum x_i = 20", "sum x_i = 50"];



max_removals = Inf;
MIP_timeLimit = 6*3600;
isRelaxation = false;



outcsvfilepath = string(expdirpath, "out.csv");
outcsvhandle = open(outcsvfilepath, "w")
row_string = string("n", ",",
					"numPermutations", ",",
					"rep", ",",
					"instance", ",",
					"text_description", ",",
					"LPB_expRevenue", ",",
					"LPB_solutiontime", ",",
					"LPB_numIterations", ",", 
					"LPB_numCutsGenerated", ",",
					"MIPB_expRevenue", ",",
					"MIPB_upperBound", ",",
					"MIPB_solutiontime", ",",
					"MIPB_gap", ",",
					"MIPB_total_num_lazy_checked", ",",
					"MIPB_total_num_lazy_generated")
row_string = string(row_string, "\n");

# print_unescaped(outcsvhandle, row_string);
print(outcsvhandle, row_string);
flush(outcsvhandle);



counter = 0;
for nind = 1:length(n_vec)
	n = n_vec[nind]
	for nPind = 1:length(numPerm_vec)
		numPermutations = numPerm_vec[nPind];

		println("**** n = $n, numPermutations = $numPermutations ****")



		println("**** Reading from data directory...")
		data_name = string(data_prefix,"_neq",n,"_Keq",numPermutations, "_v2")
		data_dir = optimalPLD_createPath(machine, "data_dir", data_name)
		revenues_mat = readcsv(string(data_dir, "revenues_mat.csv"))
		lambda_mat = readcsv(string(data_dir, "lambda_mat.csv"))
		orderings_mat = readcsv(string(data_dir, "orderings_mat.csv"))

		for instance=1:numInst
			rep = 1;
			# numPermutations = 10;
			revenues = revenues_mat[rep,:][:]
			lambda = lambda_mat[rep,:][1:numPermutations]
			orderings = orderings_mat[numPermutations*(rep-1)+1:rep*numPermutations,:]
			
			relinds = find(inst_vec .== instance);
			@show relinds
			A_feas = A_feas_mat[relinds,:];
			b_feas = b_feas_vec[relinds];
			text_description = text_description_by_inst[instance];

			# --- No warm start here ---
			x_initial = [];

			

			println("**** Solving n = $n, nP = $numPermutations, rep = $rep ****")
			counter = counter + 1;

			best_x_val, LPB_expRevenue, LPB_solutiontime,
				rankingid_list, alpha_list, beta_list, gamma_list, LPB_numIterations = optimalPLD_solveSAO_LP_Benders( n, revenues, numPermutations, lambda, orderings, MIP_timeLimit, x_initial, A_feas, b_feas)


			writecsv( string(expdirpath, "rankingid_list_instance", instance, ".csv"), rankingid_list)
			writecsv( string(expdirpath, "alpha_list_instance", instance, ".csv"), alpha_list)
			writecsv( string(expdirpath, "beta_list_instance", instance, ".csv"), beta_list)
			writecsv( string(expdirpath, "gamma_list_instance", instance, ".csv"), gamma_list)

			LPB_numCutsGenerated = length(rankingid_list);
			
			x_val, MIPB_expRevenue, MIPB_upperBound, MIPB_solutiontime, MIPB_gap, MIPB_total_num_lazy_checked, MIPB_total_num_lazy_generated,  = optimalPLD_solveSAO_Benders( n, revenues, numPermutations, lambda, orderings, MIP_timeLimit, isRelaxation, x_initial, 
																			A_feas, b_feas,
																			rankingid_list, alpha_list, beta_list, gamma_list)


			@show MIPB_total_num_lazy_checked
			@show MIPB_total_num_lazy_generated


			row_string = string(n, ",",
								numPermutations, ",",
								rep, ",",
								instance, ",",
								text_description, ",",
								LPB_expRevenue, ",",
								LPB_solutiontime, ",",
								LPB_numIterations, ",", 
								LPB_numCutsGenerated, ",",
								MIPB_expRevenue, ",",
								MIPB_upperBound, ",",
								MIPB_solutiontime, ",",
								MIPB_gap, ",",
								MIPB_total_num_lazy_checked, ",",
								MIPB_total_num_lazy_generated)
			row_string = string(row_string, "\n");
			# print_unescaped(outcsvhandle, row_string);
			print(outcsvhandle, row_string);
			flush(outcsvhandle);


		end
	end
end



close(outcsvhandle)

