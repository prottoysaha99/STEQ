#!/usr/bin/env python

###########################################################################
##    Copyright 2010 Rahul Suri and Tandy Warnow.
##    This file is part of ReUP.
##
##    ReUP is free software: you can redistribute it and/or modify
##    it under the terms of the GNU General Public License as published by
##    the Free Software Foundation, either version 3 of the License, or
##    (at your option) any later version.
##
##    ReUP is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##    GNU General Public License for more details.
##
##    You should have received a copy of the GNU General Public License
##    along with ReUP.  If not, see <http://www.gnu.org/licenses/>.
###########################################################################


# make sure you have dendropy 3.8.1 for this script.

from optparse import OptionParser, OptionGroup
from reup.reup import supplyReup
import os

def parseOptions (commandLine = None):
    '''Parse command line for options'''

    desc = ' '.join(['This script runs the SuperFine algorithm on a set of input trees',  
                     'given in a file, in Newick format.'])

    parser = OptionParser(usage = "usage: %prog [options] input_trees_file > output", 
                          version = "%prog 1.1", 
                          description = desc)

    parser.set_defaults(reconciler = "tnt",
                        numIters = 100, 
                        outputSCM = None)

    group4InfoString = ' '.join(["These options enable selection of the supertree", 
                                 "algorithm to be used as a subroutine within", 
                                 "SuperFine for resolving polytomies.  The 'qmc'",
                                 "option requires that Quartets MaxCut be", 
                                 "installed; the 'gmrp' (greedy consensus of MP", 
                                 "trees) and 'rmrp' (random MP tree) options", 
                                 "require that PAUP* be installed, the 'tnt' (MRP(TNT)) option requires TNT be installed;  and the 'rml' (MRL(RAxML)) option",
                 "requires that RAxML be installed.  Note that",
                                 "the selected subroutine's binary must be in", 
                                 "the system's executable search path."])
    #group4InfoString = ' '.join(["These options enable selection of the supertree", 
                                 #"algorithm to be used as a subroutine within", 
                                 #"SuperFine for resolving polytomies.  The 'qmc'",
                                 #"option requires that Quartets MaxCut be", 
                                 #"installed; the 'gmrp' (greedy consensus of MP", 
                                 #"trees) and 'rmrp' (random MP tree) options", 
                                 #"require that PAUP* be installed, the 'tnt' (MRP(TNT)) option requires TNT be installed; the 'fml' option",
                 #"requires that FastTree be installed; and the 'rml' option",
                 #"requires that RAxML be installed.  Note that",
                                 #"the selected subroutine's binary must be in", 
                                 #"the system's executable search path."])
                                 

    group5InfoString = ' '.join(["This option causes the SCM tree", 
                                 "to be written to disk."])

    group6DirectoryString = ' '.join(["This option allows the user to set the 					temporary", 
                                "working directory."])

    group7LogString = ' '.join(["This option allows the user to set the log output file."])
    
    group8SCMString = ' '.join(["This option allows the user give the SCM tree."])


    group4 = OptionGroup(parser, "Tree Reconciliation Options".upper(), group4InfoString)
    group5 = OptionGroup(parser, "Data Output Options".upper(), group5InfoString)
    group6 = OptionGroup(parser, "Working Directory Options".upper(), group6DirectoryString)
    group7 = OptionGroup(parser, "Log Output File".upper(), group7LogString)
    group8 = OptionGroup(parser, "SCM Tree File".upper(), group8SCMString)

    group4.add_option("-r", "--reconcile", choices = ("qmc", "gmrp", "rmrp", "rml", "tnt", "none"), 
                      dest = "reconciler", metavar = "ALG", 
                      help = "use ALG to reconcile relabeled trees, "
                             "where ALG is one of {gmrp, rml, tnt} "
                             "[default: %default]")

    group4.add_option("-n", "--numIters", type = "int", 
                      dest = "numIters", metavar = "N", 
                      help = "use N ratchet iterations when resolving with MRP "
                             "[default: %default]")

    group5.add_option("-w", "--write", dest = "writeData", metavar = "SUFFIX",  
                      help = "write merger tree and final tree to disk in same "
                             "directory as source trees file, append .SUFFIX "
                             "to written file names "
                             "[default: %default]")

    group5.add_option("-s", "--scmFile", dest = "outputSCM",  
                      help = "write merger tree to outputSCM file"
                             "[default: %default]")

    group5.add_option("-o", "--output", dest = "outputFile", metavar = "OUTPUT",
                      help = "write final tree to disk to OUTPUT file"
                             "[default: %default]")

    group6.add_option("-d", "--directory", dest = "tempDirectory",
                      help = "write files to temporary directory"
                             "[default: %default]")

    group7.add_option("-l", "--log", dest = "logFile",  
                      help = "write to log file"
                             "[default: %default]")

    group8.add_option("-t", "--tree", dest = "treeFile",
                      help = "SCM Tree"
                             "[default: %default]")

    parser.add_option_group(group4)
    parser.add_option_group(group5)
    parser.add_option_group(group6)
    parser.add_option_group(group7)
    parser.add_option_group(group8)

    if commandLine:
         (options, args) = parser.parse_args(commandLine)
    else:
        (options, args) = parser.parse_args()

    if len(args) != 1:
        parser.error("Incorrect number of arguments. Try the -h flag for help.")

    input = args[0]

    return (input, options)


# MAIN
if __name__ == '__main__':
    #Check to see if environment variable properly set
    path = os.getenv("FASTMRP")	    
    if (path is None or not os.path.isfile(path + "/mrp.jar")):
        print ("$FASTMRP variable must be defined and point to directory containing mrp.jar")

    #Increase recursion limit, can hit on very large datasets
    os.sys.setrecursionlimit(1000000)
    
    (input, options) = parseOptions()
    supplyReup(input, options)
