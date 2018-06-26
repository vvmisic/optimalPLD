"""
`optimalPLD_solveSAO_LP_Benders( n, revenues, numPermutations, lambda, orderings, MIP_timeLimit, x_initial, A_feas, b_feas)`

Solve the LP relaxation of the Bertsimas and Mišić first-choice PLD formulation using Benders decomposition.

`n` = number of candidate products.

`revenues` = array of n numbers representing marginal revenue/profit of each of the n products. (These are the pi_i's in the Bertsimas and Mišić paper.)

`numPermutations` = number of rankings / customer types.

`lambda` = array specifying the probability of each ranking / customer type; should be nonnegative and sum to 1.

`orderings` = 2D array of integers encoding each ranking. orderings[k,j] contains the index of ranking k's jth top product; (n+1) encodes the no-purchase option. For example, if n = 5 and

	orderings[k,:] = [3,5,2,6,1,4]

then the customer prefers 3 < 5 < 2 < 0 (no-purchase option).

`MIP_timeLimit` = maximum time limit on each solve of the Benders LP master problem. If set to 0, no time limit is imposed. (This argument is not used in any experiment scripts, since master LP generally solves quickly.)

`x_initial` = array with with n values (between 0 and 1) indicating the initial values of the x_i variables to use to warm-start the Benders decomposition method. If empty, no warm-start solution will be used. 

`A_feas`, `b_feas` = 2D array and 1D array specifying the system of linear inequalities Ax <= b for the PLD problem. This can be used to impose side constraints (this is the Cx <= d constraint system in the formulation in Bertsimas and Mišić). If b_feas is empty, then no constraints will be added (i.e., it will solve the unconstrained problem).

Returns:

`x_val` = array of optimal x_i variables for the LP relaxation.

`expRevenue` = expected per-customer revenue/profit of the LP solution.

`solutiontime` = time required to solve the problem.

`rankingid_list`, `alpha_list`, `beta_list`, `gamma_list` = 1D, 2D, 2D and 1D arrays containing information on cuts generated in Benders LP solution process. Each row of these arrays corresponds to a cut for a particular customer type. rankingid_list[i] is the customer type of the ith cut; alpha_list[i,:], beta_list[i,:] and gamma_list[i] are the dual variable values defining the ith cut (see Bertsimas and Mišić for more details.)


"""
function optimalPLD_solveSAO_LP_Benders( n, revenues, numPermutations, lambda, orderings, MIP_timeLimit, x_initial, A_feas, b_feas)

	tic()
	maxRev = maximum(revenues);

	if (MIP_timeLimit > 0)
		m = Model( solver = GurobiSolver(TimeLimit = MIP_timeLimit, Method = 2, Crossover = 0))
	else
		m = Model( solver = GurobiSolver(Method = 2, Crossover = 0))
	end

	@variable(m, 0 <= x[1:n] <= 1)
	@variable(m, 0 <= subobj[1:numPermutations] <= maxRev)

	for t=1:length(b_feas)
		@constraint(m, sum( A_feas[t,i]*x[i] for i=1:n) <= b_feas[t])
	end

	@objective(m, Max, sum( lambda[k] * subobj[k] for k=1:numPermutations))

	solve(m)

	

	# BEGIN preprocessing for inline 
	temp_rev = [revenues;0];
	revenues_ordered_mat = zeros(numPermutations, n+1);
	sorted_rev_inds_mat = zeros(Int64, numPermutations, n+1);
	p2r_mat = zeros(Int64, numPermutations, n+1);
	for k in 1:numPermutations
		single_ordering = convert(Array{Int64}, orderings[k,:][:]); 
		revenues_ordered =  temp_rev[single_ordering[:]]
		revenues_ordered_mat[k,:] = revenues_ordered
		sorted_rev_inds_mat[k,:] = sortperm(revenues_ordered[:], rev= true);
		p2r_mat[k,:] = sortperm(single_ordering[:]); # single_ordering is r2p
	end

	# END preprocessing for inline test

	bestUB = maxRev;
	bestLB = 0;
	x_val = getvalue(x)[:];
	if (!isempty(x_initial))
		oracle_subobj_vec = zeros(numPermutations)
		for k=1:numPermutations
			temp_ordering = convert(Array{Int64}, orderings[k,:][:]); 
			temp_subobj_val, alpha_val, beta_val, gamma_val = oracle_inline_fn( n, x_initial, temp_ordering, 
																				revenues_ordered_mat[k,:][:], temp_rev, 
																				sorted_rev_inds_mat[k,:][:], p2r_mat[k,:][:] ) 


			oracle_subobj_vec[k] = temp_subobj_val;
		end
		x_val = x_initial
		bestLB = sum( oracle_subobj_vec .* lambda);
	end
	best_x_val = x_val;


	rankingid_list = zeros(Int64,0);
	alpha_list = zeros(Float64,0,n+1);
	beta_list = zeros(Float64,0,n+1);
	gamma_list = zeros(Float64,0);


	tol_UB_LB = 0.01;
	tol_subobj = 0.01;
	foundAcceptableSolution = false;
	itercounter = 0;

	while ( !foundAcceptableSolution )
		itercounter = itercounter + 1;
		bestUB = getobjectivevalue(m);
		solution_subobj_vec = getvalue(subobj)[:];
		oracle_subobj_vec = zeros(numPermutations)
		if (itercounter != 1)
			x_val = getvalue(x)[:];
		end
		
		loop_start = time();
		alpha_mat = zeros(Float64, numPermutations,n+1);
		beta_mat = zeros(Float64, numPermutations,n+1);
		gamma_mat = zeros(Float64, numPermutations);
		oracle_subobj_vec = zeros(Float64, numPermutations)

		for k=1:numPermutations
			temp_ordering = convert(Array{Int64}, orderings[k,:][:]); 
			temp_subobj_val, alpha_val, beta_val, gamma_val = oracle_inline_fn( n, x_val, temp_ordering, 
																				revenues_ordered_mat[k,:][:], temp_rev, 
																				sorted_rev_inds_mat[k,:][:], p2r_mat[k,:][:] ) 
			alpha_mat[k,:] = alpha_val
			beta_mat[k,:] = beta_val
			gamma_mat[k] = gamma_val
			oracle_subobj_vec[k] = temp_subobj_val;
		end

		

		for k=1:numPermutations
			if (oracle_subobj_vec[k] < solution_subobj_vec[k] - tol_subobj)
				alpha_val = alpha_mat[k,:][:]
				beta_val = beta_mat[k,:][:]
				gamma_val = gamma_mat[k];
				@constraint(m,  subobj[k] <= gamma_val + sum(alpha_val[i] * x[i] for i=1:n) + sum(beta_val[i] * (1-x[i]) for i=1:n) )

				push!(rankingid_list, k)
				alpha_list = vcat(alpha_list, alpha_val' )
				beta_list = vcat(beta_list, beta_val')
				push!(gamma_list, gamma_val)
			end
		end
		loop_etime = time() - loop_start
		
		@show loop_etime
		@show oracle_subobj_vec
		@show solution_subobj_vec
		

		currentLB = sum( oracle_subobj_vec .* lambda);
		if (currentLB > bestLB)
			bestLB = currentLB
			best_x_val = x_val;
		end
		gap = round( (bestUB - bestLB) / bestUB * 100, 2);
		println("Iter $itercounter : bestLB = $bestLB, bestUB = $bestUB, gap = $gap")

		if (gap < tol_UB_LB)
			foundAcceptableSolution = true;
		else
			solve(m);
		end
	end

	solutiontime = toc()

	println(best_x_val)

	expRevenue = getobjectivevalue(m)


	return best_x_val, expRevenue, solutiontime, rankingid_list, alpha_list, beta_list, gamma_list, itercounter
end



function oracle_inline_fn( n::Int64, x_val::Array{Float64}, single_ordering::Array{Int64}, revenues_ordered::Array{Float64}, temp_rev::Array{Float64}, sorted_rev_inds::Array{Int64}, p2r::Array{Int64} ) 
#function oracle_inline_fn( n, x_val, single_ordering, revenues_ordered, temp_rev, sorted_rev_inds, p2r ) 
	temp_x = [x_val; 1.0];
	x_ordered = temp_x[single_ordering];

	remainingbudget = 1.0;
	ub = 1 - x_ordered;


	#minimizing_inds = zeros(Int32, size(y))
	minimizing_ind = 0;

	is_alpha_nonneg = zeros(Int64, n+1)
	is_beta_nonneg = zeros(Int64, n+1)
	is_beta_link = zeros(Int64, n+1);
	#is_gamma = zeros(Int64, n+1);
	gamma_ind = -1;
	# prev_remainingbudget = remainingbudget;

	# tic()
	# Primal phase
	guess_revenue = 0.0;
	y_rev_ind = 0.0;
	@inbounds for cur_prod=1:(n+1)
		# println("Loop body")
		rev_ind = sorted_rev_inds[cur_prod];

		if (rev_ind == 1)
			min_ub = Inf;
		else
			min_ub = minimum(ub[1:rev_ind-1]);
		end

		y_rev_ind = 0.0;
		if ( remainingbudget <= min_ub )
			if (remainingbudget <= x_ordered[rev_ind])
				minimizing_ind = 1
				y_rev_ind = remainingbudget
			else
				minimizing_ind = 2
				y_rev_ind = x_ordered[rev_ind]
			end
		else
			if (x_ordered[rev_ind] <= min_ub) ###
				minimizing_ind = 2
				y_rev_ind = x_ordered[rev_ind]
			else
				minimizing_ind = 2 + indmin(ub[1:rev_ind-1]);
				y_rev_ind = min_ub;
			end
		end

		guess_revenue += y_rev_ind * revenues_ordered[rev_ind];

		# update 
		# prev_remainingbudget = remainingbudget
		remainingbudget = remainingbudget - y_rev_ind; #y[rev_ind];
		#ub[1:rev_ind-1] = ub[1:rev_ind-1] - y[rev_ind];
		ub[1:rev_ind-1] -=  y_rev_ind; #y[rev_ind];

		# println("Before dual conditionals")
		if (minimizing_ind == 2)
			is_alpha_nonneg[rev_ind] = 1
		elseif (minimizing_ind > 2)
			temp = minimizing_ind - 2; #minimizing_inds[rev_ind] - 2;
			
			if (is_beta_nonneg[temp] == 0)
				is_beta_nonneg[temp] = 1
				is_beta_link[rev_ind] = 1;
				# @show temp
				# @show rev_ind
			end
		else
			#is_gamma[rev_ind] = 1;
			gamma_ind = rev_ind;
			break;
		end
		#cur_prod = cur_prod+1;
	end

	#guess_revenue = sum( y .* revenues_ordered);
	#guess_revenue = dot(y, revenues_ordered);
	# println("Guess revenue = $guess_revenue")


	# Dual phase
	# println("inline: Start of dual phase")
	#gamma_ind =  find(is_gamma)
	if (gamma_ind == -1) #(isempty(gamma_ind))
		gamma_val = 0.0;
	else
		gamma_val = revenues_ordered[gamma_ind];
	end
	
	beta_inds = find(is_beta_nonneg .> 0)
	beta_link_inds = find(is_beta_link .> 0)
	beta_val = zeros(Float64, n+1);

	if (!isempty(beta_inds))
		beta_val[beta_inds[1]] = max(revenues_ordered[beta_link_inds[1]] - gamma_val, 0)
		for j=2:length(beta_inds)
			beta_val[beta_inds[j]] = max(revenues_ordered[beta_link_inds[j]] - gamma_val - sum(beta_val[1:(beta_link_inds[j]-1)]) , 0);
		end
		# println("here")
	end

	alpha_inds = find(is_alpha_nonneg)
	alpha_val = zeros(Float64,n+1);
	for i=1:length(alpha_inds)
		# alpha_val[alpha_inds[i]] = maximum( [ revenues_ordered[alpha_inds[i]] - gamma_val - sum( beta_val[1:(alpha_inds[i]-1)] ), 0])
		alpha_val[alpha_inds[i]] = max( revenues_ordered[alpha_inds[i]] - gamma_val - sum( beta_val[1:(alpha_inds[i]-1)] ), 0)
	end

	# guess_dual_rev = gamma_val + sum(alpha_val .* x_ordered[:]) + sum(beta_val .* (1 - x_ordered)[:])

	alpha_val = alpha_val[p2r[:]];
	beta_val = beta_val[p2r[:]];

	#guess_dual_rev_unsorted = gamma_val + sum(alpha_val .* temp_x[:]) + sum(beta_val .* (1 - temp_x)[:])
	subobj_val = guess_revenue;

	return subobj_val, alpha_val, beta_val, gamma_val  #, y
end



