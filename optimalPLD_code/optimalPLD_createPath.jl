"""
`optimalPLD_createPath(machine, typeofpath, varargin...)`

Function to create a path to a particular folder; returns a string with the specified path.

`machine` = integer specifying computer where code is running. A value of 2 means to use the immediate subdirectory of the current working directory.

`typeofpath` = string specifying whether to create path to results, data, etc.

`varargin` = additional arguments that specify the directory (see code definition).

# Example

If the current directory is `/Users/example/Documents/optimalPLD/optimalPLD_exec_scripts/`, then

`optimalPLD_createPath(2, "data_dir", "testDataSet")`

will return

`/Users/example/Documents/optimalPLD/optimalPLD_data/testDataSet/`



"""
function optimalPLD_createPath(machine, typeofpath, varargin...)
	if (machine == 1)
		# VVM's macbook pro
		root = "/Users/vvmisic/Documents/optimalPLD/"
	elseif (machine == 2)
		# Catch-all -- go one directory up from the directory where .jl script
		# was executed.
		root = "../"
	else
		error("Unknown machine argument; use either 1 or 2")
	end

	if (typeofpath == "exp_dir")
		expname = varargin[1]
		path = "$root/optimalPLD_testing/$expname/";
	elseif (typeofpath == "data_dir")
		dataname = varargin[1]
		path = "$root/optimalPLD_data/$dataname/";
	else
		error("Unknown typeofpath! See optimalPLD_createPath code.")
	end

	return path

end