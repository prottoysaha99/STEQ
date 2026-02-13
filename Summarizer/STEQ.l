%option noyywrap

%{
#include "GeneTree.cpp" 
#include "Central_Tendency.cpp"

ofstream token;

map<string, int> taxaNum;
map<string, int> istaxaExists;

vector<Node> geneTree;
vector<GeneTree> geneTrees;
vector<string> taxalist;

stack<Node> st;

int nodeId = -1;
int leafCnt = 0;

double geneTreeCnt = 0.0;
bool isBinary;

enum NQDMethod {
	NQD_BASELINE = 0,
	NQD_OPTIMIZED = 1,
	NQD_OPTIMIZED_ALL_LCA = 2
};

NQDMethod nqdMethod = NQD_BASELINE;
bool validateNQD = false;
bool benchmarkNQD = false;
int validatePairs = 20000;
double validateTolerance = 1e-12;

void deallocate(double **matrix);

void Initialize_GeneTree() {
	st.pop();
	nodeId = -1;
	leafCnt = 0;
	geneTree.clear();
	taxaNum.clear();
	isBinary = true;
}

int subtreeId(int id) {
	// returns the subtreeNo. 
	// 0 --> leftSubtree, 1 --> rightSubtree, 2 --> middleSubtree, -1 -->  Indicates non-binary Tree

	if (geneTree[id].maxId[0] == id) return 0;
	else if (geneTree[id].maxId[1] == id) return 1;
	else if (geneTree[id].maxId[2] == id) return 2;
	else return -1;
}

double getDistance(vector <double> distance, int type) {
	if (type == 0) return getMean(distance);
	else if (type == 1) return getMedian(distance);
	else if (type == 2) return getMode(distance);
	else return getFromND(distance, 2.0);
}

bool PopulateDistanceMatrix(double **matrix, int type, NQDMethod method) {
	int taxaSize = taxalist.size();
	bool isCompleteDistanceMatrix = true;
	
	for(int i = 0; i < taxaSize - 1; i++){
		for(int j = i + 1; j < taxaSize; j++){
			vector<double> distances;
			for(int treeNo=0;treeNo<geneTreeCnt;treeNo++) {
				if (geneTrees[treeNo].isPresent(taxalist[i],taxalist[j])) {
					if (method == NQD_OPTIMIZED) {
						distances.push_back(geneTrees[treeNo].findDistOptimized(taxalist[i],taxalist[j]));
					}
					else if (method == NQD_OPTIMIZED_ALL_LCA) {
						distances.push_back(geneTrees[treeNo].findDistOptimizedAllLCA(taxalist[i],taxalist[j]));
					}
					else {
						distances.push_back(geneTrees[treeNo].findDistBaseline(taxalist[i],taxalist[j]));
					}
				}
			}

			if (distances.size()) {
				matrix[i][j] =  matrix[j][i] = getDistance(distances, type);
			}

			else {
				matrix[i][j] = matrix[j][i] = -1.0;
				isCompleteDistanceMatrix = false;			
			}
		}
	}

	return isCompleteDistanceMatrix;
}

double** AllocateMatrix() {
	int taxaSize = taxalist.size();
	double** matrix = new double*[taxaSize];
	for (int i = 0; i < taxaSize; i++) {
		matrix[i] = new double[taxaSize];
		for (int j = 0; j < taxaSize; j++) matrix[i][j] = 0.0;
	}
	return matrix;
}

bool CompareMatrices(double **a, double **b, double tolerance, double &worstDiff, string &worstCell) {
	int taxaSize = taxalist.size();
	worstDiff = 0.0;
	worstCell.clear();

	for (int i = 0; i < taxaSize; i++) {
		for (int j = 0; j < taxaSize; j++) {
			double diff = fabs(a[i][j] - b[i][j]);
			if (diff > worstDiff) {
				worstDiff = diff;
				worstCell = taxalist[i] + "," + taxalist[j];
			}
			if (diff > tolerance) return false;
		}
	}
	return true;
}

bool ValidateMethods(int centralTendency) {
	bool ok = true;
	double globalWorstOptimized = 0.0;
	double globalWorstAllLCA = 0.0;
	string worstTreePairOptimized = "";
	string worstTreePairAllLCA = "";

	for (int treeNo = 0; treeNo < geneTreeCnt; treeNo++) {
		unsigned int seed = 1337 + treeNo;

		double worstOpt = 0.0;
		string pairOpt;
		bool treeOkOpt = geneTrees[treeNo].validateOptimized(validatePairs, seed, validateTolerance, worstOpt, pairOpt);
		if (worstOpt > globalWorstOptimized) {
			globalWorstOptimized = worstOpt;
			worstTreePairOptimized = "tree=" + to_string(treeNo) + ", pair=(" + pairOpt + ")";
		}
		if (!treeOkOpt) {
			ok = false;
			cout << "[NQD-VALIDATION] optimized mismatch in " << worstTreePairOptimized << ", diff=" << setprecision(15) << globalWorstOptimized << endl;
			break;
		}

		double worstAll = 0.0;
		string pairAll;
		bool treeOkAll = geneTrees[treeNo].validateOptimizedAllLCA(validatePairs, seed, validateTolerance, worstAll, pairAll);
		if (worstAll > globalWorstAllLCA) {
			globalWorstAllLCA = worstAll;
			worstTreePairAllLCA = "tree=" + to_string(treeNo) + ", pair=(" + pairAll + ")";
		}
		if (!treeOkAll) {
			ok = false;
			cout << "[NQD-VALIDATION] optimized_all_lca mismatch in " << worstTreePairAllLCA << ", diff=" << setprecision(15) << globalWorstAllLCA << endl;
			break;
		}
	}

	double **baseline = AllocateMatrix();
	double **optimized = AllocateMatrix();
	double **optimizedAllLCA = AllocateMatrix();
	PopulateDistanceMatrix(baseline, centralTendency, NQD_BASELINE);
	PopulateDistanceMatrix(optimized, centralTendency, NQD_OPTIMIZED);
	PopulateDistanceMatrix(optimizedAllLCA, centralTendency, NQD_OPTIMIZED_ALL_LCA);

	double matrixWorstOpt = 0.0;
	string matrixCellOpt;
	bool matrixOkOpt = CompareMatrices(baseline, optimized, validateTolerance, matrixWorstOpt, matrixCellOpt);

	double matrixWorstAll = 0.0;
	string matrixCellAll;
	bool matrixOkAll = CompareMatrices(baseline, optimizedAllLCA, validateTolerance, matrixWorstAll, matrixCellAll);

	deallocate(baseline);
	deallocate(optimized);
	deallocate(optimizedAllLCA);

	if (!matrixOkOpt) {
		ok = false;
		cout << "[NQD-VALIDATION] optimized matrix mismatch at (" << matrixCellOpt << "), diff=" << setprecision(15) << matrixWorstOpt << endl;
	}
	if (!matrixOkAll) {
		ok = false;
		cout << "[NQD-VALIDATION] optimized_all_lca matrix mismatch at (" << matrixCellAll << "), diff=" << setprecision(15) << matrixWorstAll << endl;
	}

	if (ok) {
		cout << "[NQD-VALIDATION] passed, worst optimized diff=" << setprecision(15) << globalWorstOptimized
			 << ", worst optimized_all_lca diff=" << setprecision(15) << globalWorstAllLCA << endl;
	}

	return ok;
}

void BenchmarkMethods(int centralTendency) {
	using namespace chrono;

	auto t0 = high_resolution_clock::now();
	double **baseline = AllocateMatrix();
	auto t1 = high_resolution_clock::now();
	bool completeBaseline = PopulateDistanceMatrix(baseline, centralTendency, NQD_BASELINE);
	auto t2 = high_resolution_clock::now();

	auto prepStart = high_resolution_clock::now();
	size_t memoryBytesOptimized = 0;
	for (int treeNo = 0; treeNo < geneTreeCnt; treeNo++) {
		geneTrees[treeNo].preprocessOptimized();
		memoryBytesOptimized += geneTrees[treeNo].optimizedMemoryBytes();
	}
	auto prepEnd = high_resolution_clock::now();

	double **optimized = AllocateMatrix();
	auto t3 = high_resolution_clock::now();
	bool completeOptimized = PopulateDistanceMatrix(optimized, centralTendency, NQD_OPTIMIZED);
	auto t4 = high_resolution_clock::now();

	auto prepAllStart = high_resolution_clock::now();
	size_t memoryBytesAllLCA = 0;
	for (int treeNo = 0; treeNo < geneTreeCnt; treeNo++) {
		geneTrees[treeNo].preprocessAllLeafPairLCA();
		memoryBytesAllLCA += geneTrees[treeNo].allLeafPairLCAMemoryBytes();
	}
	auto prepAllEnd = high_resolution_clock::now();

	double **optimizedAllLCA = AllocateMatrix();
	auto t5 = high_resolution_clock::now();
	bool completeOptimizedAllLCA = PopulateDistanceMatrix(optimizedAllLCA, centralTendency, NQD_OPTIMIZED_ALL_LCA);
	auto t6 = high_resolution_clock::now();

	double matrixWorstOpt = 0.0;
	string matrixCellOpt;
	bool matrixOkOpt = CompareMatrices(baseline, optimized, validateTolerance, matrixWorstOpt, matrixCellOpt);

	double matrixWorstAll = 0.0;
	string matrixCellAll;
	bool matrixOkAll = CompareMatrices(baseline, optimizedAllLCA, validateTolerance, matrixWorstAll, matrixCellAll);

	deallocate(baseline);
	deallocate(optimized);
	deallocate(optimizedAllLCA);

	double allocBaselineMs = duration<double, milli>(t1 - t0).count();
	double fillBaselineMs = duration<double, milli>(t2 - t1).count();
	double prepMs = duration<double, milli>(prepEnd - prepStart).count();
	double allocOptimizedMs = duration<double, milli>(t3 - prepEnd).count();
	double fillOptimizedMs = duration<double, milli>(t4 - t3).count();
	double prepAllLCAMs = duration<double, milli>(prepAllEnd - prepAllStart).count();
	double allocOptimizedAllLCAMs = duration<double, milli>(t5 - prepAllEnd).count();
	double fillOptimizedAllLCAMs = duration<double, milli>(t6 - t5).count();

	double totalBaselineMs = allocBaselineMs + fillBaselineMs;
	double totalOptimizedMs = prepMs + allocOptimizedMs + fillOptimizedMs;
	double totalOptimizedAllLCAMs = prepAllLCAMs + allocOptimizedAllLCAMs + fillOptimizedAllLCAMs;
	double ratioOptimized = (totalOptimizedMs > 0.0) ? (totalBaselineMs / totalOptimizedMs) : 0.0;
	double ratioOptimizedAllLCA = (totalOptimizedAllLCAMs > 0.0) ? (totalBaselineMs / totalOptimizedAllLCAMs) : 0.0;

	cout << "[NQD-BENCHMARK]" << endl;
	cout << "  baseline_fill_ms         : " << fillBaselineMs << endl;
	cout << "  baseline_total_ms        : " << totalBaselineMs << endl;
	cout << "  optimized (RMQ LCA):" << endl;
	cout << "  optimized_preprocess_ms  : " << prepMs << endl;
	cout << "  optimized_fill_ms        : " << fillOptimizedMs << endl;
	cout << "  optimized_total_ms       : " << totalOptimizedMs << endl;
	cout << "  optimized_speedup        : " << ratioOptimized << endl;
	cout << "  optimized_memory_bytes~  : " << memoryBytesOptimized << endl;
	cout << "  optimized_all_lca:" << endl;
	cout << "  all_lca_precompute_ms    : " << prepAllLCAMs << endl;
	cout << "  all_lca_fill_ms          : " << fillOptimizedAllLCAMs << endl;
	cout << "  all_lca_total_ms         : " << totalOptimizedAllLCAMs << endl;
	cout << "  all_lca_speedup          : " << ratioOptimizedAllLCA << endl;
	cout << "  all_lca_memory_bytes~    : " << memoryBytesAllLCA << endl;

	if (!matrixOkOpt || completeBaseline != completeOptimized) {
		cout << "[NQD-BENCHMARK] warning: optimized mismatch at (" << matrixCellOpt << "), diff=" << setprecision(15) << matrixWorstOpt << endl;
	}
	if (!matrixOkAll || completeBaseline != completeOptimizedAllLCA) {
		cout << "[NQD-BENCHMARK] warning: optimized_all_lca mismatch at (" << matrixCellAll << "), diff=" << setprecision(15) << matrixWorstAll << endl;
	}
}

bool ParseOptionalArgs(int argc, char *argv[]) {
	for (int i = 6; i < argc; i++) {
		string arg = argv[i];
		if (arg == "--validate_nqd" || arg == "--debug_compare_methods") {
			validateNQD = true;
		}
		else if (arg == "--benchmark_nqd") {
			benchmarkNQD = true;
		}
		else if (arg == "--nqd_method") {
			if (i + 1 >= argc) {
				cout << "Missing value for --nqd_method (baseline|optimized|optimized_all_lca)" << endl;
				return false;
			}
			string method = argv[++i];
			if (method == "baseline") nqdMethod = NQD_BASELINE;
			else if (method == "optimized") nqdMethod = NQD_OPTIMIZED;
			else if (method == "optimized_all_lca") nqdMethod = NQD_OPTIMIZED_ALL_LCA;
			else {
				cout << "Invalid --nqd_method value: " << method << endl;
				return false;
			}
		}
		else if (arg == "--validate_pairs") {
			if (i + 1 >= argc) {
				cout << "Missing value for --validate_pairs" << endl;
				return false;
			}
			validatePairs = atoi(argv[++i]);
		}
		else {
			cout << "Unknown argument: " << arg << endl;
			return false;
		}
	}

	return true;
}

void write_to_File(double **matrix) {
	int taxaSize = taxalist.size();

	token << taxaSize << endl;
	
	for(int i = 0; i < taxaSize; i++){
		token<< taxalist[i] << " ";

		for(int j = 0; j < taxaSize; j++) token<< matrix[i][j] << " ";

		token << endl;
	}
}

void InferTreeFromMatrix(bool isComplete, string matrixPath, string treePath, int fastmeType) {
	string str;

	if (isComplete){
		string fastme_exec = string("../Binaries/fastme-2.1.5-linux64") + string(" -i ") + matrixPath + string(" -o ") + treePath;

		if (fastmeType == 0){ // fastme_BAL
			str = fastme_exec + string(" -n -m 'BAL' > wastage.txt && 'rm' wastage.txt ") + matrixPath;
		}
		else if (fastmeType == 1){ //fastme_NJ
			str = fastme_exec + string(" -n -m 'NJ' > wastage.txt && 'rm' wastage.txt ") + matrixPath;
		}
		else if (fastmeType == 2){ // fastme_UNJ
			str = fastme_exec + string(" -m 'UNJ' > wastage.txt && 'rm' wastage.txt ") + matrixPath;
		}
	}

	else{
		str = string("../Binaries/BIONJ ") +  matrixPath + string(" ") + treePath + string(" > wastage.txt && 'rm' wastage.txt ") + matrixPath;
	}

	system(str.c_str());
}

void deallocate(double **matrix) {
	for (int i = 0; i < taxalist.size(); i++) delete[] matrix[i];
	delete[] matrix;
}

%}

WHITESPACE		 [ \t\f\r\v\n]+ 
LETTER 			 [a-zA-Z_]
DIGIT 			 [0-9] 
ALPHANUMERIC	 [a-zA-Z0-9_.-]
STRING           {LETTER}+
INTEGER 		 {DIGIT}+
FLOAT			 [-]?{DIGIT}*(".")?{DIGIT}+(e[+-]?({DIGIT}+))? 

%%

{WHITESPACE} {
		}
		
"("		{
			nodeId++;

			if (st.empty()){
				Node root(nodeId, -1);
				geneTree.push_back(root);
				st.push(root);
			}	
			else{
				Node parent = st.top();
				Node internalNode(nodeId, parent.id);
				geneTree.push_back(internalNode);
				st.push(internalNode);
			}
		}
		
")"({ALPHANUMERIC}*)	{
			int childId = st.top().id;			
			st.pop();
			int id = st.top().id;

			int subtreeNo = subtreeId(id);

			if (subtreeNo == -1) {
				isBinary = false;
				cout << "Non-Binary Tree" << endl;
			}
			else {
				geneTree[id].maxId[subtreeNo] = max(geneTree[childId].maxId[0], geneTree[childId].maxId[1]);
				for (int i = 0; i < 3; i++) geneTree[id].descendants[subtreeNo] += geneTree[childId].descendants[i];
			}
		}	

","		{
			int childId = st.top().id;			
			st.pop();
			int id = st.top().id; // id of current node id

			int subtreeNo = subtreeId(id);

			if (subtreeNo == -1) {
				isBinary = false;
				cout << "Non-Binary Tree" << endl;
			}
			else {
				geneTree[id].maxId[subtreeNo] = max(geneTree[childId].maxId[0], geneTree[childId].maxId[1]);
				for (int i = 0; i < 3; i++) geneTree[id].descendants[subtreeNo] += geneTree[childId].descendants[i];
			}
		}
		
":"{WHITESPACE}*({FLOAT}|{INTEGER})		{ 
		}
		
";"		{
			if (isBinary){
				GeneTree gt(geneTree, taxaNum, leafCnt);
				geneTrees.push_back(gt);
				geneTreeCnt += 1.0;
			}
			Initialize_GeneTree();
		}

{ALPHANUMERIC}+	{		
			string taxa(yytext); 
			nodeId++;
			
			Node parent = st.top();
			Node leafNode(nodeId, parent.id);
			leafNode.descendants[0] = 1;

			geneTree.push_back(leafNode);
			taxaNum[taxa] = nodeId;
			st.push(leafNode);

			leafCnt += 1;

			if (!istaxaExists[taxa]){
				taxalist.push_back(taxa);
				istaxaExists[taxa] = 1;
			}	   
		}
	
. {
}	
	
%%

int main(int argc,char *argv[]){
	/*
		argv[1] = input file path
		argv[2] = type of central tendency; 0: mean, 1: median, 2: mode, 3: mean from considering values from mean - sigma * sd to mean + sigma * sd
		argv[3] = output matrix file path
		argv[4] = output tree file path
		argv[5] = fastme method; 0: fastme_BAL, 1: fastme_NJ, 2: fastme_UNJ
	*/
	if(argc < 6){
		printf("Usage: STEQ.out <input> <central_tendency> <output_matrix> <output_tree> <fastme_method> [--nqd_method baseline|optimized|optimized_all_lca] [--validate_nqd] [--benchmark_nqd] [--validate_pairs N]\n");
		return 0;
	} 

	if (!ParseOptionalArgs(argc, argv)) return 0;
	
	FILE *fin=fopen(argv[1],"r");

	if(fin == NULL){
		printf("Cannot open specified file\n");
		return 0;
	}
	yyin= fin;
	
	token.open(argv[3]); 
	token << std::setprecision(4) << std::fixed;

	geneTreeCnt = 0.0;
	isBinary = true;
	int type = *argv[2] - '0';
	int fastmeType = *argv[5] - '0'; 
	
	yylex(); 
	
	int taxaSize = taxalist.size();

	// matrix to store distance between any pair of taxa
	double** matrix = AllocateMatrix();

	if (validateNQD) {
		bool ok = ValidateMethods(type);
		if (!ok) {
			deallocate(matrix);
			fclose(yyin);
			token.close();
			return 1;
		}
	}

	if (benchmarkNQD) {
		BenchmarkMethods(type);
	}
	
	// populate the distance matrix where each entry represents the distance of two species and returns if the matrix is complete
	bool isCompleteDistanceMatrix = PopulateDistanceMatrix(matrix, type, nqdMethod);
	
	write_to_File(matrix);

	InferTreeFromMatrix(isCompleteDistanceMatrix, argv[3], argv[4], fastmeType);

	deallocate(matrix);

	fclose(yyin);
	token.close();

	return 0;
}
