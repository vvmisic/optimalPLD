"""
`optimalPLD_solveSAO_LP_Benders( n, revenues, numPermutations, lambda, orderings, MIP_timeLimit, x_initial, A_feas, b_feas)`

Solve the integer version of the Bertsimas and Mišić first-choice PLD formulation using Benders decomposition. 

This function will formulate the Benders reformulation in Gurobi by specifying a lazy constraint callback to generate the Benders cuts. Cuts are generated using the closed-form expressions for the Benders cuts when the x_i's are integer (see Section 4.1 of Bertsimas and Mišić).

`n` = number of candidate products.

`revenues` = array of n numbers representing marginal revenue/profit of each of the n products. (These are the pi_i's in the Bertsimas and Mišić paper.)

`numPermutations` = number of rankings / customer types.

`lambda` = array specifying the probability of each ranking / customer type; should be nonnegative and sum to 1.

`orderings` = 2D array of integers encoding each ranking. orderings[k,j] contains the index of ranking k's jth top product; (n+1) encodes the no-purchase option. For example, if n = 5 and

	orderings[k,:] = [3,5,2,6,1,4]

then the customer prefers 3 < 5 < 2 < 0 (no-purchase option).

`MIP_timeLimit` = maximum time limit for Gurobi to solve the integer problem. If set to 0, no time limit is imposed. 

`x_initial` = array with with n values (between 0 and 1) indicating the initial values of the x_i variables to use to warm-start the Benders decomposition method. If empty, no warm-start solution will

`A_feas`, `b_feas` = 2D array and 1D array specifying the system of linear inequalities Ax <= b for the PLD problem. This can be used to impose side constraints (this is the Cx <= d constraint system in the formulation in Bertsimas and Mišić 2018). If b_feas is empty, then no constraints will be added (i.e., it will solve the unconstrained problem).

`rankingid_list`, `alpha_list`, `beta_list`, `gamma_list` = 1D, 2D, 2D and 1D arrays containing information on cuts generated in Benders LP solution process. These cuts are added to the Benders integer formulation so as to warm-start it. 


Returns:

`x_val` = array of optimal x_i variables for the integer problem.

`expRevenue` = expected per-customer revenue/profit of the integer solution.

`upperBound` = best upper bound on the optimal per-customer revenue/profit.

`solutiontime` = time required to solve the problem.

`MIPgap` = optimality gap of the best integer solution found.

`total_num_lazy_checked` = number of lazy constraints that are checked (incremented by 1 for each customer type and each invocation of the lazy constraint callback)

`total_num_lazy_generated` = number of lazy constraints added to the formulation (incremented by 1 every time a lazy constraint for any customer type is found to be violated and is added to the formulation)


"""
function optimalPLD_solveSAO_Benders( n, revenues, numPermutations, lambda, orderings, MIP_timeLimit, isRelaxation, x_initial,
								A_feas, b_feas,
								rankingid_list, alpha_list, beta_list, gamma_list )

	maxRev = maximum(revenues);

	if (MIP_timeLimit > 0)
		m = Model( solver = GurobiSolver(TimeLimit = MIP_timeLimit))
	else
		m = Model( solver = GurobiSolver())
	end



	if (!isRelaxation)
		@variable(m, x[1:n], Bin)
	else
		@variable(m, 0 <= x[1:n] <= 1)
	end

	if (!isRelaxation)
		if (!isempty(x_initial))
			for i=1:n
				setvalue(x[i], x_initial[i]);
			end
		end
	end

	@variable(m, 0 <= subobj[1:numPermutations] <= maxRev)

	for t=1:length(b_feas)
		@constraint(m, sum( A_feas[t,i]*x[i] for i=1:n) <= b_feas[t])
	end

	for t=1:length(rankingid_list)
		k = rankingid_list[t];
		@constraint(m, subobj[k] <= gamma_list[t] + sum( alpha_list[t,i]*x[i] for i=1:n) + sum(beta_list[t,i]*(1-x[i]) for i=1:n) )
	end

	@objective(m, Max, sum( lambda[k] * subobj[k] for k=1:numPermutations))

	const total_num_lazy_generated = 0;
	const total_num_lazy_checked = 0;

	orderings = convert( Array{Int64}, orderings);
	p2r_mat = zeros(Int64, numPermutations, n+1);
	temp_rev = [revenues;0];
	revenues_ordered_mat = zeros(numPermutations, n+1);
	sorted_rev_inds_mat = zeros(Int64, numPermutations, n+1);

	for k in 1:numPermutations
		single_ordering = convert(Array{Int64}, orderings[k,:][:]); 
		revenues_ordered =  temp_rev[single_ordering[:]]
		revenues_ordered_mat[k,:] = revenues_ordered
		sorted_rev_inds_mat[k,:] = sortperm(revenues_ordered[:], rev= true);
		p2r_mat[k,:] = sortperm(single_ordering[:]); # single_ordering is r2p
	end

	tol_integrality = 1e-4;
	# Compute current lower bound.
	function benders_cb(cb)
		x_val = getvalue(x)[:]

		sol_is_integral = all( i->( (i>1-tol_integrality) &  (i<tol_integrality)), x_val )
		
		# println("candSet = $candSet")
		subobj_val = getvalue(subobj)[:]

		candSet_ponly = find(x_val .> 0.5)
		candSet = [candSet_ponly; n+1];
		for k=1:numPermutations
			total_num_lazy_checked = total_num_lazy_checked + 1;
			#p2r = sortperm(orderings[k,:][:])
			gamma_val = 0;
			alpha_val = zeros(n);
			beta_val = zeros(n+1);
			#temp = indmin(p2r_mat[k,candSet[:]])
			temp = indmin(p2r_mat[k,candSet])
			if (candSet[temp] <= n)
				gamma_val = revenues[candSet[temp]];
			else
				gamma_val = 0.0
			end
			# println("\tgamma_val is ")
			# println(gamma_val)

			if (!isempty(candSet_ponly))
				maxrev = maximum(revenues[candSet_ponly[:]]);
				beta_val[candSet[temp]] = maxrev - gamma_val
			end
			# println("\tbeta_val[temp] is (temp = $temp) ")
			# println(beta_val[candSet[temp]])

			#alpha_val = zeros(size(alpha_val))
			for i2=1:(n+1)
				if (orderings[k,i2] <= n)
					alpha_val[orderings[k,i2]] = max(revenues[orderings[k,i2]] - gamma_val - sum( beta_val[orderings[k,1:(i2-1)]]), 0 )
				end
			end


			if (gamma_val < subobj_val[k])
				# println("\tperm $k: gamma_val = $gamma_val, subobj_val[k] = $(subobj_val[k])")
				@lazyconstraint(cb, gamma_val + sum( alpha_val[i] * x[i] for i=1:n) + sum( beta_val[i] * (1 - x[i]) for i=1:n) >= subobj[k])
				total_num_lazy_generated = total_num_lazy_generated+1;
			end
		end
	end

	addlazycallback(m, benders_cb)

	solve(m)

	println(getvalue(x))

	x_val = getvalue(x)[:]

	solutiontime = Gurobi.get_runtime(m.internalModel.inner);
	if (!isRelaxation)
		MIPgap = Gurobi.get_dblattr(m.internalModel.inner, "MIPGap");
	else
		MIPgap = -1;
	end

	expRevenue = getobjectivevalue(m)

	upperBound = Gurobi.get_dblattr(m.internalModel.inner, "ObjBound");


	return x_val, expRevenue, upperBound, solutiontime, MIPgap, total_num_lazy_checked, total_num_lazy_generated
end