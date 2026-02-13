# STEQ

## Create .out file from .l file
* flex **_filename.l_** && g++ lex.yy.c -lfl -o **_filename.out_** && rm lex.yy.c  


## To run this
* cd STEQ
* ./Summarizer/LexToOut
* **./Summarizer/STEQ.out** (**input_file_path**) (**Central_Tendency_Type**) (**OutputMatrix_Path**) (**OutputTree_Path**) (**FASTME**) [--dist-calc-method all-lca|rmq|path-traversal]
	* Central_Tendency: 0 --> Mean, 1 --> Median, 2 --> Mode, 3 --> Normal Distribution
	* FastME: 0 --> FastME_BAL, 1 --> FastME_NJ
	* Distance calculation method:
		* `all-lca` (default): precompute all leaf-pair LCAs + O(1) query
		* `rmq`: Euler-tour RMQ-based LCA + O(1) query
		* `path-traversal`: original naive path traversal
