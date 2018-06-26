"""
`optimalPLD_solveSAO_McBZ_constraints( n, revenues, numPermutations, lambda, utilities, MIP_timeLimit, isRelaxation, x_initial, A_feas, b_feas )`

Formulate and solve the McBride and Zufryden first-choice PLD formulation (see McBride and Zufryden 1988). 

`n` = number of candidate products.

`revenues` = array of n numbers representing marginal revenue/profit of each of the n products. (These are the pi_i's in the Bertsimas and Mišić paper.)

`numPermutations` = number of rankings / customer types.

`lambda` = array specifying the probability of each ranking / customer type; should be nonnegative and sum to 1.

`utilities` = 2D array of real numbers specifying the utility of each product for each customer type (i.e., u[k,j] is the utility of option j for customer type k); (n+1) encodes the no-purchase option.

`MIP_timeLimit` = maximum time to run Gurobi for. If set to 0, then no time limit is set.

`isRelaxation` = boolean to indicate whether we are solving LP relaxation (true) or the integer problem (false)

`x_initial` = array with with n values (either 0 or 1) indicating the initial values of the x_i variables to use to warm-start Gurobi. If empty, no warm-start solution will be given to Gurobi. Only allowed for integer problem (isRelaxation = false).

`A_feas`, `b_feas` = 2D array and 1D array specifying the system of linear inequalities Ax <= b for the PLD problem. This can be used to impose side constraints (this is the Cx <= d constraint system in the formulation in Bertsimas and Mišić 2018). If b_feas is empty, then no constraints will be added (i.e., it will solve the unconstrained problem).

Returns:

`x_val`, `y_val` = arrays of optimal x_i and y^k_i variables for either the integer problem or the relaxation.

`expRevenue` = expected per-customer revenue/profit of the best integer solution; this is the optimal objective if the problem is solved to full optimality.

`upperbound` = best upper bound on the optimal expected per customer revenue/profit found at termination; equal to the optimal objective (modulo Gurobi's default tolerance) if the problem is solved to full optimality.

`solutiontime` = time required to solve the problem.

`MIPgap` = optimality gap upon termination.


# References

R. D. McBride and F. S. Zufryden. An integer programming approach to the optimal product line selection problem. Marketing Science, 7(2):126–140, 1988.

"""
function optimalPLD_solveSAO_McBZ_constraints( n, revenues, numPermutations, lambda, utilities, MIP_timeLimit, isRelaxation, x_initial, A_feas, b_feas )

	if (MIP_timeLimit > 0)
		m = Model( solver = GurobiSolver(TimeLimit = MIP_timeLimit))
	else
		m = Model( solver = GurobiSolver())
	end

	#weights = rand(Gamma())

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

	Umax = maximum(utilities,2)
	
	# @variable(m, y[1:(n+1), 1:numPermutations] >= 0 )

	if (!isRelaxation)
		@variable(m, y[1:(n+1), 1:numPermutations], Bin )
	else
		@variable(m, 0 <= y[1:(n+1), 1:numPermutations] <= 1 )
	end

	for t=1:length(b_feas)
		@constraint(m, sum( A_feas[t,i]*x[i] for i=1:n) <= b_feas[t])
	end


	@objective(m, Max, sum( lambda[k] * revenues[i] * y[i,k] for i=1:n, k=1:numPermutations) )
	# for t=1:20
	# 	@constraint(m, sum{ x[i], i=(5*(t-1)+1):(5*t)} <= 1)
	# end
	# println("Here")
	for k=1:numPermutations
		for i=1:n
			@constraint(m, y[i,k] <= x[i])
		end
		for i=1:n
			for i2=1:n
				if (i != i2)
					@constraint(m, utilities[k,i] * x[i] >= utilities[k,i2] * x[i2] - Umax[k] * (1 - y[i,k]) )
				end
			end
			@constraint(m, utilities[k,i]*x[i] >= utilities[k,n+1] - Umax[k] * (1 - y[i,k]) )
		end

		for i2=1:n
			@constraint(m, utilities[k,n+1] >= utilities[k,i2] * x[i2] - Umax[k] * (1 - y[n+1,k]) )
		end

		@constraint(m, sum( y[i,k] for i=1:(n+1)) == 1)
	end


	solve(m)

	println(getvalue(x))

	x_val = getvalue(x)[:]
	y_val = getvalue(y)[:,:]

	solutiontime = Gurobi.get_runtime(m.internalModel.inner);
	if (!isRelaxation)
		upperBound = Gurobi.get_dblattr(m.internalModel.inner, "ObjBound");
		MIPgap = Gurobi.get_dblattr(m.internalModel.inner, "MIPGap");
	else
		upperBound = -1;
		MIPgap = -1;
	end

	expRevenue = getobjectivevalue(m)


	return x_val, y_val, expRevenue, upperBound, solutiontime, MIPgap
end