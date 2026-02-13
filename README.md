# STEQ

## Create .out file from .l file
* flex **_filename.l_** && g++ lex.yy.c -lfl -o **_filename.out_** && rm lex.yy.c  


## To run this
* cd STEQ
* ./Summarizer/LexToOut
* **./Summarizer/STEQ.out** (**input_file_path**) (**Central_Tendency_Type**) (**OutputMatrix_Path**) (**OutputTree_Path**) (**FASTME**) [--nqd_method baseline|optimized|optimized_all_lca] [--validate_nqd] [--benchmark_nqd] [--validate_pairs N]
	* Central_Tendency: 0 --> Mean, 1 --> Median, 2 --> Mode, 3 --> Normal Distribution
	* FastME: 0 --> FastME_BAL, 1 --> FastME_NJ
	* NQD method:
		* `baseline` (default): original per-path traversal
		* `optimized`: LCA + prefix-sum query after per-tree preprocessing
		* `optimized_all_lca`: same optimized formula + precomputed all leaf-pair LCA lookup table
	* `--validate_nqd`: compares baseline vs optimized on sampled/all pairs and on the final distance matrix
	* `--benchmark_nqd`: prints preprocessing/fill/total timing for both methods and speedup
