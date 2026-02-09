#include<bits/stdc++.h>
using namespace std;


class Node{
    public:
        int id;
        int parentId;

        /*  
            maxId[x]: 
            x = maxId of [0(left subtree), 1(right subtree), 2(middle subtree)]
        */
        int maxId[3];
        int descendants[3];

        Node(int Id, int p){ 
            maxId[0] = maxId[1] = maxId[2] = id = Id;
            parentId = p;
            descendants[0] = descendants[1] = descendants[2] = 0;
        } 

        void PrintNode(){
            printf("id = %d, parentId = %d, child = (%d, %d, %d), descendants = (%d, %d, %d)\n", id, parentId, maxId[0], maxId[1], maxId[2], descendants[0], descendants[1], descendants[2]);
        }
};

class GeneTree{
    public:
    vector<Node> tree;  
    map<string,int>taxaMap;
    int leafCnt;
    
    GeneTree(vector<Node> genetree, map<string,int>tNum, int leafcnt){
        tree =  genetree;
        taxaMap = tNum;
        leafCnt = leafcnt;

        for (int i = 1; i < tree.size(); i++) {
            tree[i].descendants[2] = leafCnt - tree[i].descendants[0] - tree[i].descendants[1];
        }
    }

    bool isPresent(string taxa1, string taxa2) {
        return taxaMap[taxa1] && taxaMap[taxa2];
    }

    int subtreeNo(int current_nodeId, int leafId) {
        if (leafId > current_nodeId && leafId <= tree[current_nodeId].maxId[0]) return 0;
        else if (leafId > tree[current_nodeId].maxId[0] && leafId <= tree[current_nodeId].maxId[1]) return 1;
        else return 2;
    }

    double findDist(string taxa1, string taxa2){
        int m = taxaMap[taxa1];
        int n = taxaMap[taxa2];

        int x = min(m,n);
        int y = max(m,n);

        double dist = 0.0;

        int currentId = tree[x].parentId;
        int subtreeX = subtreeNo(currentId, x);
        int subtreeY = subtreeNo(currentId, y);

        int a, b, R, lca;

        while (subtreeX != subtreeY) {
            a = tree[currentId].descendants[subtreeX];
            b = tree[currentId].descendants[subtreeY];
            R = leafCnt - a - b;

            dist += calc(a, b, R);

            lca = currentId;

            if (currentId == 0) break;
            currentId = tree[currentId].parentId;
            subtreeX = subtreeNo(currentId, x);
            subtreeY = subtreeNo(currentId, y);
        }

        currentId = tree[y].parentId;

        while(currentId != lca) {
            subtreeX = subtreeNo(currentId, x);
            subtreeY = subtreeNo(currentId, y);

            a = tree[currentId].descendants[subtreeX];
            b = tree[currentId].descendants[subtreeY];
            R = leafCnt - a - b;

            dist += calc(a, b, R);
            currentId = tree[currentId].parentId;
        }

        return dist / 2;
    }

    int calc(int aStar, int bStar, int r){
        return (aStar + bStar - 2);
    }
};
