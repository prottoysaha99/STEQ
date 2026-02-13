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

    private:
    bool optimizedPreprocessed;
    int rootId;

    vector<vector<int> > children;
    vector<int> depth;
    vector<long long> subtreeLeafCount;
    vector<long long> parentSideLeafCount;
    vector<long long> downPrefix;
    vector<long long> upPrefix;

    vector<int> euler;
    vector<int> eulerDepth;
    vector<int> firstOccurrence;
    vector<int> lg2;
    vector<vector<int> > rmq;

    public:
    
    GeneTree(vector<Node> genetree, map<string,int>tNum, int leafcnt){
        tree =  genetree;
        taxaMap = tNum;
        leafCnt = leafcnt;
        optimizedPreprocessed = false;
        rootId = 0;

        for (int i = 1; i < tree.size(); i++) {
            tree[i].descendants[2] = leafCnt - tree[i].descendants[0] - tree[i].descendants[1];
        }
    }

    bool isPresent(string taxa1, string taxa2) {
        map<string, int>::const_iterator it1 = taxaMap.find(taxa1);
        map<string, int>::const_iterator it2 = taxaMap.find(taxa2);
        return (it1 != taxaMap.end() && it1->second > 0 && it2 != taxaMap.end() && it2->second > 0);
    }

    int subtreeNo(int current_nodeId, int leafId) {
        if (leafId > current_nodeId && leafId <= tree[current_nodeId].maxId[0]) return 0;
        else if (leafId > tree[current_nodeId].maxId[0] && leafId <= tree[current_nodeId].maxId[1]) return 1;
        else return 2;
    }

    double findDistBaseline(string taxa1, string taxa2){
        if (taxa1 == taxa2) return 0.0;
        int m = taxaMap[taxa1];
        int n = taxaMap[taxa2];

        int x = min(m,n);
        int y = max(m,n);

        double dist = 0.0;

        int currentId = tree[x].parentId;
        int subtreeX = subtreeNo(currentId, x);
        int subtreeY = subtreeNo(currentId, y);

        int a, b, R;
        int lca = rootId;

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

    double findDist(string taxa1, string taxa2){
        return findDistBaseline(taxa1, taxa2);
    }

    private:
    void buildEulerTour(int nodeId) {
        if (firstOccurrence[nodeId] == -1) firstOccurrence[nodeId] = (int)euler.size();
        euler.push_back(nodeId);
        eulerDepth.push_back(depth[nodeId]);

        for (int i = 0; i < (int)children[nodeId].size(); i++) {
            int child = children[nodeId][i];
            buildEulerTour(child);
            euler.push_back(nodeId);
            eulerDepth.push_back(depth[nodeId]);
        }
    }

    int lca(int u, int v) const {
        int left = firstOccurrence[u];
        int right = firstOccurrence[v];
        if (left > right) swap(left, right);

        int len = right - left + 1;
        int k = lg2[len];
        int idx1 = rmq[k][left];
        int idx2 = rmq[k][right - (1 << k) + 1];
        return (eulerDepth[idx1] <= eulerDepth[idx2]) ? euler[idx1] : euler[idx2];
    }

    public:
    void preprocessOptimized() {
        if (optimizedPreprocessed) return;

        int nodeCount = (int)tree.size();
        if (nodeCount == 0) {
            optimizedPreprocessed = true;
            return;
        }

        rootId = 0;
        for (int i = 0; i < nodeCount; i++) {
            if (tree[i].parentId == -1) {
                rootId = i;
                break;
            }
        }

        children.assign(nodeCount, vector<int>());
        depth.assign(nodeCount, 0);
        subtreeLeafCount.assign(nodeCount, 0);
        parentSideLeafCount.assign(nodeCount, 0);
        downPrefix.assign(nodeCount, 0);
        upPrefix.assign(nodeCount, 0);

        for (int i = 0; i < nodeCount; i++) {
            int parent = tree[i].parentId;
            if (parent != -1) children[parent].push_back(i);
        }

        vector<int> preorder;
        preorder.reserve(nodeCount);
        stack<int> dfsStack;
        dfsStack.push(rootId);
        depth[rootId] = 0;

        while (!dfsStack.empty()) {
            int u = dfsStack.top();
            dfsStack.pop();
            preorder.push_back(u);

            for (int i = (int)children[u].size() - 1; i >= 0; i--) {
                int v = children[u][i];
                depth[v] = depth[u] + 1;
                dfsStack.push(v);
            }
        }

        for (int i = (int)preorder.size() - 1; i >= 0; i--) {
            int u = preorder[i];
            if (children[u].empty()) subtreeLeafCount[u] = 1;
            else {
                long long subtotal = 0;
                for (int j = 0; j < (int)children[u].size(); j++) {
                    subtotal += subtreeLeafCount[children[u][j]];
                }
                subtreeLeafCount[u] = subtotal;
            }
        }

        parentSideLeafCount[rootId] = 0;
        downPrefix[rootId] = 0;
        upPrefix[rootId] = 0;

        for (int i = 1; i < (int)preorder.size(); i++) {
            int u = preorder[i];
            int parent = tree[u].parentId;
            parentSideLeafCount[u] = (long long)leafCnt - subtreeLeafCount[u];
            downPrefix[u] = downPrefix[parent] + subtreeLeafCount[u];
            upPrefix[u] = upPrefix[parent] + parentSideLeafCount[u];
        }

        euler.clear();
        eulerDepth.clear();
        firstOccurrence.assign(nodeCount, -1);
        buildEulerTour(rootId);

        int m = (int)euler.size();
        lg2.assign(m + 1, 0);
        for (int i = 2; i <= m; i++) lg2[i] = lg2[i / 2] + 1;

        int maxLog = lg2[m] + 1;
        rmq.assign(maxLog, vector<int>(m, 0));
        for (int i = 0; i < m; i++) rmq[0][i] = i;

        for (int k = 1; k < maxLog; k++) {
            int span = 1 << k;
            int half = span >> 1;
            for (int i = 0; i + span <= m; i++) {
                int idx1 = rmq[k - 1][i];
                int idx2 = rmq[k - 1][i + half];
                rmq[k][i] = (eulerDepth[idx1] <= eulerDepth[idx2]) ? idx1 : idx2;
            }
        }

        optimizedPreprocessed = true;
    }

    double findDistOptimized(string taxa1, string taxa2){
        if (taxa1 == taxa2) return 0.0;
        preprocessOptimized();

        int x = taxaMap[taxa1];
        int y = taxaMap[taxa2];

        int w = lca(x, y);
        int parentX = tree[x].parentId;
        int parentY = tree[y].parentId;

        long long upToParentY = (parentY == -1) ? 0LL : upPrefix[parentY];
        long long upToParentX = (parentX == -1) ? 0LL : upPrefix[parentX];

        long long sumCx = (downPrefix[x] - downPrefix[w]) + (upToParentY - upPrefix[w]);
        long long sumCy = (downPrefix[y] - downPrefix[w]) + (upToParentX - upPrefix[w]);

        long long internalCount = (long long)depth[x] + depth[y] - 2LL * depth[w] - 1LL;
        if (internalCount < 0) internalCount = 0;

        return 0.5 * (double)(sumCx + sumCy - 2LL * internalCount);
    }

    bool validateOptimized(int pairSampleLimit, unsigned int seed, double tolerance, double &worstDiff, string &worstPair) {
        vector<string> taxa;
        taxa.reserve(taxaMap.size());
        for (map<string, int>::const_iterator it = taxaMap.begin(); it != taxaMap.end(); ++it) {
            if (it->second > 0) taxa.push_back(it->first);
        }

        int n = (int)taxa.size();
        if (n < 2) return true;

        long long totalPairs = 1LL * n * (n - 1) / 2;
        bool checkAll = (pairSampleLimit <= 0 || totalPairs <= pairSampleLimit);

        worstDiff = 0.0;
        worstPair.clear();
        mt19937 rng(seed);

        if (checkAll) {
            for (int i = 0; i < n - 1; i++) {
                for (int j = i + 1; j < n; j++) {
                    double b = findDistBaseline(taxa[i], taxa[j]);
                    double o = findDistOptimized(taxa[i], taxa[j]);
                    double diff = fabs(b - o);
                    if (diff > worstDiff) {
                        worstDiff = diff;
                        worstPair = taxa[i] + "," + taxa[j];
                    }
                    if (diff > tolerance) return false;
                }
            }
            return true;
        }

        unordered_set<long long> sampled;
        sampled.reserve(pairSampleLimit * 2);

        while ((int)sampled.size() < pairSampleLimit) {
            int i = rng() % n;
            int j = rng() % n;
            if (i == j) continue;
            if (i > j) swap(i, j);
            long long key = 1LL * i * n + j;
            if (sampled.find(key) != sampled.end()) continue;
            sampled.insert(key);

            double b = findDistBaseline(taxa[i], taxa[j]);
            double o = findDistOptimized(taxa[i], taxa[j]);
            double diff = fabs(b - o);
            if (diff > worstDiff) {
                worstDiff = diff;
                worstPair = taxa[i] + "," + taxa[j];
            }
            if (diff > tolerance) return false;
        }

        return true;
    }

    size_t optimizedMemoryBytes() {
        preprocessOptimized();

        size_t bytes = 0;
        bytes += children.size() * sizeof(vector<int>);
        for (int i = 0; i < (int)children.size(); i++) bytes += children[i].capacity() * sizeof(int);
        bytes += depth.capacity() * sizeof(int);
        bytes += subtreeLeafCount.capacity() * sizeof(long long);
        bytes += parentSideLeafCount.capacity() * sizeof(long long);
        bytes += downPrefix.capacity() * sizeof(long long);
        bytes += upPrefix.capacity() * sizeof(long long);
        bytes += euler.capacity() * sizeof(int);
        bytes += eulerDepth.capacity() * sizeof(int);
        bytes += firstOccurrence.capacity() * sizeof(int);
        bytes += lg2.capacity() * sizeof(int);
        bytes += rmq.size() * sizeof(vector<int>);
        for (int i = 0; i < (int)rmq.size(); i++) bytes += rmq[i].capacity() * sizeof(int);

        return bytes;
    }

    int calc(int aStar, int bStar, int r){
        return (aStar + bStar - 2);
    }
};
