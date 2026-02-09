#!/lusr/bin/python
'''
Created on Jun 3, 2011

@author: smirarab, modified by Bayzid
'''
import dendropy
import sys
import os
import os.path

if __name__ == '__main__':
    treeName = sys.argv[1]
    inputs = sys.argv[2]
    files = inputs.split(',')

    
    for file in files:
        sample = open(file)
        included = [s[:-1] for s in sample.readlines()]

        resultsFile = "%s%s%s" % (os.path.dirname(sample.name), "/in.", os.path.basename(sample.name))
        
        trees = dendropy.TreeList.get_from_path(treeName, 'newick')
        filt = lambda node: True if (node is not None and node.label not in included) else False
        for tree in trees:
            tree.retain_taxa_with_labels(included)
           
        trees.write(file=open(resultsFile,'w'),schema='newick') 
        os.remove(file)
