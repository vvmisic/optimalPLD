"""
`optimalPLD_runDNCHeuristicObjFn(P, obj_fn, total_num_products, maxIter, heuristic_seed)`

Run the divide and conquer heuristic (see Green and Krieger 1993, Belloni et al. 2008). 

`P` = number of products in the product line.

`obj_fn` = arbitrary function that will compute expected revenue of a product line (= an array of product indices). Should only include actual products (not no-purchase option).

`total_num_products` = total number of candidate products from which the product line may be chosen.

`maxIter` = integer specifying number of random starting points to use.

`heuristic_seed` = integer to set the seed of the random number generator (to be used for generating the random starting points).



# References

P. E. Green and A. M. Krieger. Conjoint analysis with product-positioning applications. In J. Eliashberg and G. L. Lilien, editors, Handbooks in Operations Research and Management Science, volume 5, pages 467–515. Elsevier, 1993.

A. Belloni, R. Freund, M. Selove, and D. Simester. Optimizing product line designs: Efficient methods and comparisons. Management Science, 54(9):1544–1552, 2008.

"""
function optimalPLD_runDNCHeuristicObjFn(P, obj_fn, total_num_products, maxIter, heuristic_seed) 

tic()


heuristicRevenues = zeros(maxIter,1);
heuristicSolns = zeros(maxIter,P);

srand(heuristic_seed)

for iter=1:maxIter
	currentProduct = 1;
	chosenProducts = randperm(total_num_products)[1:P]

	currentRevenue = obj_fn(chosenProducts)
	initialRevenue = currentRevenue;

	numExhaustedProducts = 0;

	while (numExhaustedProducts < P)
		# Keep going

		# Evaluate
		bestExpectedCandRevenue = 0;
		bestCandProduct = -1;
		for p2=1:total_num_products
			if (!in(p2, chosenProducts))
				cand_chosenProducts = copy(chosenProducts);
				cand_chosenProducts[currentProduct] = p2;
                expectedCandRevenue = obj_fn(cand_chosenProducts); #sum(lambdaFull .* revenues_sorted); #sum(sum(lambda .* revenues_sorted));
                if (bestExpectedCandRevenue < expectedCandRevenue)
                	bestCandProduct = p2;
                	bestExpectedCandRevenue = expectedCandRevenue;
                end
			end
		end

		if (bestExpectedCandRevenue > currentRevenue)
			numExhaustedProducts = 1;
			chosenProducts[currentProduct] = bestCandProduct;
			currentRevenue = bestExpectedCandRevenue;
		else
			numExhaustedProducts = numExhaustedProducts + 1;
		end

		currentProduct = currentProduct + 1;
		if (currentProduct > P)
			currentProduct = 1;
		end
	end

	print(currentRevenue)
	print("\t")
	println(chosenProducts)
	heuristicRevenues[iter] = currentRevenue
	heuristicSolns[iter,:] = chosenProducts

end

bestSolnInd = indmax(heuristicRevenues)
bestHeuristicSoln = heuristicSolns[bestSolnInd,:][:]
bestHeuristicRevenue = maximum(heuristicRevenues)

heuristicTime = toc()

return bestHeuristicRevenue, bestHeuristicSoln, heuristicTime

end