"""
`optimalPLD_createRanking_rfn(n, revenues, numPermutations, orderings, lambda)`

Create a function that will compute the expected per-customer profit/revenue of a product line according to the first-choice model specified by the input arguments.

`n` = number of candidate products.

`revenues` = array of n numbers representing marginal revenue/profit of each of the n products. (These are the pi_i's in the Bertsimas and Mišić paper.)

`numPermutations` = number of rankings / customer types.

`orderings` = 2D array of integers encoding each ranking. orderings[k,j] contains the index of ranking k's jth top product; (n+1) encodes the no-purchase option. For example, if n = 5 and

	orderings[k,:] = [3,5,2,6,1,4]

then the customer prefers 3 < 5 < 2 < 0 (no-purchase option).

`lambda` = array specifying the probability of each ranking / customer type; should be nonnegative and sum to 1.

"""
function optimalPLD_createRanking_rfn(n, revenues, numPermutations, orderings, lambda)
	unpermed_orderings = zeros(Int64,size(orderings))
	for k=1:numPermutations
		unpermed_orderings[k,:] = sortperm(orderings[k,:][:])
	end
	revenue_fn = function ( s::Array{Int64} )
		s0 = [s;n+1];
		temp = 0;
		# @show revenues
		for k=1:numPermutations
			# @show s 
			s_ranks = unpermed_orderings[k,s0[:]]
			# @show s_ranks
			ind = indmin(s_ranks)
			temp = temp + revenues[s0[ind]] * lambda[k]
			# @show k
			# @show 
		end

		return temp
	end

	return revenue_fn

end