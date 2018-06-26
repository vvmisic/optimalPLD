# Exact First-Choice Product Line Optimization

This repository contains all the code used in the numerical experiments in the paper

> D. Bertsimas and V. V. Mišić (2018). Exact first-choice product line optimization. Working paper; available at SSRN: https://ssrn.com/abstract=3020502. 

## Citation

If you use the code and/or data in this repository in your own research, please cite the above paper as follows:

```
@article{bertsimas2018exact,
	title={Exact first-choice product line optimization},
	author={Bertsimas, Dimitris and Mi\v{s}i\'{c}, Velibor V.}},
	journal={Working paper},
	year={2018},
	note={Available at SSRN: \url{https://ssrn.com/abstract=3020502}}
```

*In addition, if you use the processed data files in the `optimalPLD_data/toubia2003_neq3584_Keq330_v2/` directory*: 
1. **Please cite** the paper of Toubia et al. (2003). Full reference: 
  > O. Toubia, D. I. Simester, J. R. Hauser and E. Dahan (2003). Fast polyhedral adaptive conjoint estimation. *Marketing Science*, 22(3):273-303. 

  For convenience, the following .bib entry is provided:
  ```
  @article{toubia2003fast,
	title={Fast polyhedral adaptive conjoint estimation},
	author={Toubia, O. and Simester, D. I. and Hauser, J. R. and Dahan, E.},
	journal={Marketing Science},
	volume={22},
	number={3},
	pages={273--303},
	year={2003},
	publisher={INFORMS}
  }
  ```
  In addition, the use of the processed data in `optimalPLD_data/toubia2003_neq3584_Keq330_v2/` is subject to the same limitations included in the README file of the data set provided by the authors of Toubia et al. (2003). Please refer to the README file of that data for further details.

2. In addition, **please also cite** the paper of Belloni et al. (2008). Full reference: 
  > A. Belloni, R. Freund, M. Selove, and D. Simester (2008). Optimizing product line designs: Efficient methods and comparisons. *Management Science*, 54(9):1544–1552. 

  For convenience, the .bib entry is provided below:
  ```
@article{belloni2008optimizing,
	title={Optimizing product line designs: Efficient methods and comparisons},
	author={Belloni, A. and Freund, R. and Selove, M. and Simester, D.},
	journal={Management Science},
	volume={54},
	number={9},
	pages={1544--1552},
	year={2008},
	publisher={INFORMS}
}
  ```

## License 

This code is available under the MIT License.

Copyright (C) 2018 Velibor Misic

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

## Software Requirements

To run the code, you need to have:
+ Julia version 0.6.0 or later
+ JuMP version 0.18.0 or later
+ Gurobi version 7.5 or later

The code should be compatible with later versions. Future updates may be made to make the code compatible with newer releases of Julia/JuMP/Gurobi.


## Repository Structure

This repository contains several directories:

+ `optimalPLD_data/`: Contains the data files used in the numerical experiments. There are two types of subdirectories:
  + `MIOexpdata1_neq*_Keq**/`: Contains data for synthetic data instances (see Sections 5.1 and 5.2 of Bertsimas and Mišić), where `*` and `**` are the values of *n* and *K* respectively. Each directory contains data for 20 randomly generated instances.
  + `toubia2003_neq3584_Keq330_v2/`: Contains data for the real conjoint data set from Toubia et al. used in Section 5.3 of the main paper and Section EC.3 of the electronic companion.

  Each directory contains several files; the most important are:
  + `orderings_mat.csv`: Customer rankings (`\sigma^k` in the paper).
  + `lambda_mat.csv`: Customer type probabilities (`\lambda^k` in the paper).
  + `revenues_mat.csv`: Marginal product profits/revenues (`\pi_i` in the paper).

+ `optimalPLD_code/`: Contains functions needed to formulate and solve the various optimization problems in the paper, as well as helper functions (such as `optimalPLD_createPath.jl`, which creates directory strings systematically).

  Each function is formatted to include documentation that can be accessed using the `?` functionality in Julia. For example, running the following (assuming current directory is `optimalPLD/optimalPLD_code/`)
  ```
  > include("optimalPLD_solveSAO_constraints.jl")
  > ? optimalPLD_solveSAO_constraints
  ```
  will output some information on how to use `optimalPLD_solveSAO_constraints` (which formulates and solves the main formulation of Bertsimas and Mišić). 

+ `optimalPLD_exec_scripts/`: Contains scripts that can be executed to run the functions in `optimalPLD_code/` on a large swathe of instances. Each script contains instructions on how to run it. For example, 

  ``` 
  > julia MIOexpdata1_Comparison.jl 1 20 
  ```

  will solve all four formulations and their relaxations directly, without any side constraints, for *n* = 100, 200, 500, 1000, *K* = 100, 200, 500, 1000.
+ `optimalPLD_testing/`: Upon execution of scripts in `optimalPLD_exec_scripts/`, this  directory will contain directories with results (on, e.g., objective values and run times). 



