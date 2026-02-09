#! /usr/bin/python

import sys
from PRD import PRDWrapper

if ( __name__ == '__main__' ) :
    filepath = sys.argv[1]
    subsetFile = open(sys.argv[2], "w+")#file path where the subsets will be stored
    subset = int(sys.argv[3])
    overlap = int(sys.argv[4])
    dec = PRDWrapper(subset, overlap)
    output = dec.decompose_dataset(filepath)
    print(len(output))
    for i, subset in enumerate(output) :
        subsetFile.write('subset {0:2d}: {1}\n'.format(i, subset))
