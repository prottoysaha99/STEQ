#include <bits/stdc++.h>
#include "GeneTree.cpp"

using namespace std;

struct RawNode {
    string label;
    RawNode *left;
    RawNode *right;
    RawNode(const string &lab) : label(lab), left(NULL), right(NULL) {}
    RawNode(RawNode *l, RawNode *r) : label(""), left(l), right(r) {}
};

struct BuildInfo {
    int maxId;
    int leafCount;
};

RawNode* buildRandomTree(vector<string> labels, mt19937 &rng) {
    if (labels.size() == 1) return new RawNode(labels[0]);

    shuffle(labels.begin(), labels.end(), rng);
    int split = 1 + (rng() % (labels.size() - 1));

    vector<string> leftLabels(labels.begin(), labels.begin() + split);
    vector<string> rightLabels(labels.begin() + split, labels.end());

    RawNode *left = buildRandomTree(leftLabels, rng);
    RawNode *right = buildRandomTree(rightLabels, rng);
    return new RawNode(left, right);
}

void deleteTree(RawNode *node) {
    if (!node) return;
    deleteTree(node->left);
    deleteTree(node->right);
    delete node;
}

BuildInfo assignIdsAndPopulate(
    RawNode *node,
    int parentId,
    int &nextId,
    vector<Node> &nodes,
    map<string, int> &taxaMap
) {
    int id = nextId++;
    nodes.push_back(Node(id, parentId));

    if (!node->left && !node->right) {
        nodes[id].descendants[0] = 1;
        taxaMap[node->label] = id;
        return {id, 1};
    }

    BuildInfo leftInfo = assignIdsAndPopulate(node->left, id, nextId, nodes, taxaMap);
    BuildInfo rightInfo = assignIdsAndPopulate(node->right, id, nextId, nodes, taxaMap);

    nodes[id].maxId[0] = leftInfo.maxId;
    nodes[id].maxId[1] = rightInfo.maxId;
    nodes[id].descendants[0] = leftInfo.leafCount;
    nodes[id].descendants[1] = rightInfo.leafCount;

    return {max(leftInfo.maxId, rightInfo.maxId), leftInfo.leafCount + rightInfo.leafCount};
}

GeneTree buildGeneTreeFromLabels(const vector<string> &labels, mt19937 &rng) {
    RawNode *root = buildRandomTree(labels, rng);
    vector<Node> nodes;
    map<string, int> taxaMap;
    int nextId = 0;
    assignIdsAndPopulate(root, -1, nextId, nodes, taxaMap);
    deleteTree(root);
    return GeneTree(nodes, taxaMap, (int)labels.size());
}

vector<string> mapKeys(const map<string, int> &m) {
    vector<string> keys;
    keys.reserve(m.size());
    for (map<string, int>::const_iterator it = m.begin(); it != m.end(); ++it) keys.push_back(it->first);
    return keys;
}

void assertClose(double a, double b, const string &msg) {
    double diff = fabs(a - b);
    if (diff > 1e-12) {
        cerr << "ASSERTION FAILED: " << msg << " baseline=" << setprecision(15) << a
             << " optimized=" << b << " diff=" << diff << endl;
        exit(1);
    }
}

void compareTreeAllPairs(GeneTree &gt) {
    vector<string> taxa = mapKeys(gt.taxaMap);
    for (int i = 0; i < (int)taxa.size(); i++) {
        for (int j = i; j < (int)taxa.size(); j++) {
            double b = gt.findDistBaseline(taxa[i], taxa[j]);
            double o = gt.findDistOptimized(taxa[i], taxa[j]);
            double a = gt.findDistOptimizedAllLCA(taxa[i], taxa[j]);
            assertClose(b, o, "tree pair " + taxa[i] + "," + taxa[j]);
            assertClose(b, a, "tree pair(all_lca) " + taxa[i] + "," + taxa[j]);
        }
    }
}

vector<vector<double> > aggregateMatrix(vector<GeneTree> &trees, const vector<string> &allTaxa, int mode) {
    int n = (int)allTaxa.size();
    vector<vector<double> > matrix(n, vector<double>(n, 0.0));

    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            double sum = 0.0;
            int cnt = 0;
            for (int t = 0; t < (int)trees.size(); t++) {
                if (!trees[t].isPresent(allTaxa[i], allTaxa[j])) continue;
                double d = trees[t].findDistBaseline(allTaxa[i], allTaxa[j]);
                if (mode == 1) d = trees[t].findDistOptimized(allTaxa[i], allTaxa[j]);
                else if (mode == 2) d = trees[t].findDistOptimizedAllLCA(allTaxa[i], allTaxa[j]);
                sum += d;
                cnt++;
            }
            matrix[i][j] = matrix[j][i] = (cnt == 0 ? -1.0 : (sum / cnt));
        }
    }

    return matrix;
}

void compareMatrices(const vector<vector<double> > &a, const vector<vector<double> > &b) {
    for (int i = 0; i < (int)a.size(); i++) {
        for (int j = 0; j < (int)a[i].size(); j++) {
            assertClose(a[i][j], b[i][j], "matrix cell " + to_string(i) + "," + to_string(j));
        }
    }
}

void testFixedSmallTree() {
    mt19937 rng(7);
    vector<string> labels;
    labels.push_back("A");
    labels.push_back("B");
    labels.push_back("C");
    labels.push_back("D");
    labels.push_back("E");
    labels.push_back("F");

    GeneTree gt = buildGeneTreeFromLabels(labels, rng);
    compareTreeAllPairs(gt);
}

void testRandomTrees() {
    mt19937 rng(99);
    for (int t = 0; t < 40; t++) {
        int n = 8 + (rng() % 10);
        vector<string> labels;
        for (int i = 0; i < n; i++) labels.push_back("T" + to_string(t) + "_" + to_string(i));
        GeneTree gt = buildGeneTreeFromLabels(labels, rng);
        compareTreeAllPairs(gt);
    }
}

void testMissingTaxaAggregation() {
    mt19937 rng(123);

    vector<string> allTaxa;
    allTaxa.push_back("A");
    allTaxa.push_back("B");
    allTaxa.push_back("C");
    allTaxa.push_back("D");
    allTaxa.push_back("E");
    allTaxa.push_back("F");

    vector<GeneTree> trees;

    vector<string> l1;
    l1.push_back("A"); l1.push_back("B"); l1.push_back("C"); l1.push_back("D");
    vector<string> l2;
    l2.push_back("B"); l2.push_back("C"); l2.push_back("E"); l2.push_back("F");
    vector<string> l3;
    l3.push_back("A"); l3.push_back("C"); l3.push_back("D"); l3.push_back("E");

    trees.push_back(buildGeneTreeFromLabels(l1, rng));
    trees.push_back(buildGeneTreeFromLabels(l2, rng));
    trees.push_back(buildGeneTreeFromLabels(l3, rng));

    vector<vector<double> > baseline = aggregateMatrix(trees, allTaxa, 0);
    vector<vector<double> > optimized = aggregateMatrix(trees, allTaxa, 1);
    vector<vector<double> > optimizedAllLCA = aggregateMatrix(trees, allTaxa, 2);
    compareMatrices(baseline, optimized);
    compareMatrices(baseline, optimizedAllLCA);
}

int main() {
    testFixedSmallTree();
    testRandomTrees();
    testMissingTaxaAggregation();
    cout << "All NQD tests passed." << endl;
    return 0;
}
