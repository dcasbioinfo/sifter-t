#! /usr/bin/env python
#  -*- coding: iso-8859-1 -*-  

########## ########## ########### ########## ########## ########## ##########
#  Sifter-T - Sifter framework for large scale Functional Annotation.       #
#                                                                           #
#  Copyright 2013 Almeida-e-Silva, D.C.; Vêncio, R.Z.N.                     #
#  All rights reserved.                                                     #
#                                                                           #
#  If you use this work or any portion thereof in published work,           #
#  please cite it as:                                                       #
#                                                                           #
#     Almeida-e-Silva D.C. and Vêncio R.Z.N. 2013. Sifter-T: A functional   #
#     framework for large-scale probabilistic protein domain annotation.    #
#                                                                           #
########## ########## ########### ########## ########## ########## ##########

"""
 * This script is used to replace builded hundred families files for the 
   original hundred families.
"""

import os
import sys
import shutil
from optparse import OptionParser

if __name__ == '__main__':
    usage = "\n     %prog -i DIR -o DIR"
    description = "Replace builded hundred families files for the original "   \
                  "hundred families."
    #Defines input variables
    parser = OptionParser(usage=usage, version="%prog 0.1", description=description)
    parser.add_option("-i",
        dest="input_dir", 
        help="Full path for input directory [with original hundred families].",
        metavar="/DIRECTORY/DIR/")
    parser.add_option("-o",
        dest="output_dir", 
        help="Full path for output directory [with builded hundred families].",
        metavar="/DIRECTORY/DIR/")

    (options, args) = parser.parse_args()

    if len(args) > 0:
        print "\n# Extra arguments. Wrong usage. Exiting... \n"
        sys.exit(1)

    if not options.input_dir or not options.output_dir:
        print "Not all required parameters were specified. "                   \
              "Type \"-h\" for help. \nExiting..."
        sys.exit(1)

    options.input_dir = os.path.abspath(options.input_dir)+"/"
    options.output_dir = os.path.abspath(options.output_dir)+"/"
    options.stdir = os.path.abspath(os.getcwd())+"/"

    families = set()

    if os.path.exists(options.output_dir+"PF00007/"):
        with open(options.stdir+'hundred_families_paper.txt', 'r') as handle:
            for line in handle:
                families.add(line.strip())
    else:
        with open(options.stdir+'hundred_families_site.txt', 'r') as handle:
            for line in handle:
                families.add(line.strip())


    for pfam in families:
        print "Removing "+pfam+" files..."
        os.remove(options.output_dir+pfam+"/"+pfam+".pli")
        os.remove(options.output_dir+pfam+"/"+pfam+".nhx")
        os.remove(options.output_dir+pfam+"/"+pfam+".rdata")
        os.remove(options.output_dir+pfam+"/"+pfam+"-sifterout.txt")
        os.remove(options.output_dir+pfam+"/alpha-"+pfam.lower()+".fx")
        os.remove(options.output_dir+pfam+"/scale-"+pfam.lower()+".fx")
        os.remove(options.output_dir+pfam+"/infer-"+pfam.lower()+".fx")
        print "Replacing \".pli\" and \".nhx\" files..."
        shutil.copyfile(options.input_dir+pfam+".pli", options.output_dir+pfam+"/"+pfam+".pli")
        shutil.copyfile(options.input_dir+pfam+".nhx", options.output_dir+pfam+"/"+pfam+".nhx")


    shutil.move(options.output_dir+"REPORT.txt", options.output_dir+"old_REPORT.txt")
    shutil.move(options.output_dir+"REPORT.tab", options.output_dir+"old_REPORT.tab")
    sys.exit()

