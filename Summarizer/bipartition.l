%option noyywrap

%{
#include<bits/stdc++.h>
using namespace std;


class Node{
public:
    //Node *leftSubtree, *rightSubtree, *middleSubtree;
    Node *subTree[3];
    bool leaf;
    string name;

    Node(){
        //leftSubtree = rightSubtree = middleSubtree = NULL;
        for (int i = 0; i < 3; i++) subTree[i] = NULL;
        leaf = false;
    }

    Node (string taxaName){
        leaf = true;
        name = taxaName;
    }
};


ofstream token;
ofstream outfile;

int leafCnt = 0;
Node *root;
stack<Node*> st;
int subsetsCnt = 1;


void printName(Node *current, int subsetId){
    string filename = "../subsets/subset_" + to_string(subsetId) + ".tre";
    outfile.open(filename);
    stack<Node*> st;
    st.push(current); 

    while(!st.empty()){
        Node *top = st.top();
        st.pop();

        if (top->leaf) {
            outfile << top->name << endl;
        }
        else {
            for (int i = 0; i < 3; i++){
                if (top->subTree[i] != NULL) st.push(top->subTree[i]);
            }
        }
    }

    outfile.close();
}


int getBipartition(Node *current, Node *parent, int subtreeId, int maxTaxa){
    if (current == NULL) return 0; 
    if (current->leaf) return 1;

    int l = 0, m = 0, r = 0;
    
    l = getBipartition(current->subTree[0], current, 0, maxTaxa);
    r = getBipartition(current->subTree[1], current, 1, maxTaxa);
    m = getBipartition(current->subTree[2], current, 2, maxTaxa);

    if (leafCnt - l - r - m < maxTaxa) {
        if (parent == NULL) {
            printName(current, subsetsCnt++);
            outfile << endl;

            return 0;
        }
        return l + r + m;
    }

    if (l + r + m >= maxTaxa) { 
        printName(current, subsetsCnt++);
        outfile << endl;

        leafCnt -= l + r + m; 
        parent->subTree[subtreeId] = NULL;
        current = NULL;

        return 0;
    }
    return l + r + m;
}


%}

WHITESPACE		 [ \t\f\r\v\n]+ 
LETTER 			 [a-zA-Z_]
DIGIT 			 [0-9] 
STRING           {LETTER}+
INTEGER 		 {DIGIT}+
FLOAT			 [-]?{DIGIT}*(".")?{DIGIT}+(e[+-]?({DIGIT}+))? 

%%

{WHITESPACE} {
		}
 
		
"("		{
			if (st.empty()) st.push(root);
			
			else {
				Node *parent = st.top();
				Node *internalNode = new Node();

                if (parent->subTree[0] == NULL) parent->subTree[0] = internalNode;
                else if (parent->subTree[1] == NULL) parent->subTree[1] = internalNode;
                else parent->subTree[2] = internalNode;

				st.push(internalNode);
			}
		}
")"{FLOAT}* {
		    st.pop();
		}	
","         {
            st.pop();
        }
		
":"{WHITESPACE}*({FLOAT}|{INTEGER})		{ 
		}

";" {
        st.pop();
}

({INTEGER}|{FLOAT}|{STRING})	{		
			string taxa(yytext); 
			Node *parent = st.top();
			Node *leafNode = new Node(taxa);

            if (parent->subTree[0] == NULL) parent->subTree[0] = leafNode;
            else if (parent->subTree[1] == NULL) parent->subTree[1] = leafNode;
            else parent->subTree[2] = leafNode;

			st.push(leafNode);
			leafCnt += 1;   
		}
	
. {
}	
%%

int main(int argc,char *argv[]){
	FILE *fin=fopen(argv[1],"r");

    if(fin == NULL){
        printf("Cannot open specified file\n");
        return 0;
	}

    root = new Node();
    int threshold = 10;
	yyin= fin;

    yylex();
	
    int _ = getBipartition(root, NULL, -1, threshold);

    cout << subsetsCnt << endl;

	fclose(yyin);
	return subsetsCnt;
}


