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

bool PopulateDistanceMatrix(double **matrix, int type) {
	int taxaSize = taxalist.size();
	bool isCompleteDistanceMatrix = true;
	
	for(int i = 0; i < taxaSize - 1; i++){
		for(int j = i + 1; j < taxaSize; j++){
			vector<double> distances;
			for(int treeNo=0;treeNo<geneTreeCnt;treeNo++) {
				if (geneTrees[treeNo].isPresent(taxalist[i],taxalist[j])) {
					distances.push_back(geneTrees[treeNo].findDist(taxalist[i],taxalist[j]));
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
	if(argc != 6){
		printf("Please provide input file name and try again\n");
		return 0;
	} 
	
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
	double** matrix = new double*[taxaSize];
	for(int i = 0; i < taxaSize; i++){
		matrix[i] = new double[taxaSize];

		for(int j = 0; j < taxaSize; j++) matrix[i][j] = 0.0;
	}
	
	// populate the distance matrix where each entry represents the distance of two species and returns if the matrix is complete
	bool isCompleteDistanceMatrix = PopulateDistanceMatrix(matrix, type);
	
	write_to_File(matrix);

	InferTreeFromMatrix(isCompleteDistanceMatrix, argv[3], argv[4], fastmeType);

	deallocate(matrix);

	fclose(yyin);
	token.close();

	return 0;
}
