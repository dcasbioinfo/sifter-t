#! /usr/bin/env python

########## ########## ########### ########## ########## ########## ##########
#  Sifter-T - Sifter framework for large scale Functional Annotation.       #
#                                                                           #
#  Copyright 2013 Almeida-e-Silva, D.C.; Vencio, R.Z.N.                     #
#  All rights reserved.                                                     #
#                                                                           #
#  If you use this work or any portion thereof in published work,           #
#  please cite it as:                                                       #
#                                                                           #
#     Almeida-e-Silva D.C. and Vencio R.Z.N. Sifter-T: A functional         #
#     framework for large-scale probabilistic protein domain annotation.    #
#     (In preparation...)                                                   #
#                                                                           #
########## ########## ########### ########## ########## ########## ##########

"""
 * This script automate the inference engine execution.

    run_sifter1 (options, java_options, scodes)
    run_sifter0 (options, java_options, scodes)
    sifter_multi(options, familylist, scodes)
    _main()

"""

from optparse import OptionParser
from multiprocessing import JoinableQueue
from multiprocessing import Process
from Queue import Empty
from time import sleep
import pickle
import os
import sys
q = JoinableQueue()    
n = JoinableQueue()   

def run_sifter0 (options, java_options, scodes):
    '''
    Run sifter2.0 using "-generate" parameter.
    '''
    global q
    while 1:
        try:
            fam = q.get(block=True, timeout=0.1)
        except Empty: 
            break
        else:
            os.system("java "+java_options+" -jar " \
    ""+options.sdir.replace(" ","\ ")+"sifter.jar " \
    ""+fam+" --output "+options.outdir.replace(" ","\ ")+"" \
    ""+fam.upper()+"/"+fam+".rdata --protein " \
    ""+options.outdir.replace(" ","\ ")+fam.upper()+"/" \
    ""+fam+".pli --reconciled " \
    ""+options.outdir.replace(" ","\ ")+"" \
    ""+fam.upper()+"/"+fam+".nhx -v -g --em --ontology " \
    ""+options.dbdir.replace(" ","\ ")+"function.ontology"\
    " --familyfile "+options.outdir.replace(" ","\ ")+"" \
    ""+fam.upper()+"/infer-"+fam.lower()+".fx --scale " \
    ""+options.outdir.replace(" ","\ ")+fam.upper()+"" \
    "/scale-"+fam.lower()+".fx --alpha " \
    ""+options.outdir.replace(" ","\ ")+fam.upper()+"" \
    "/alpha-"+fam.lower()+".fx "+scodes+" > "  \
    ""+options.outdir.replace(" ","\ ")+fam.upper()+"/"+fam+"-sifterout.txt")
            q.task_done() 


def run_sifter1 (options, java_options, scodes):
    '''
    Run sifter2.0 for probability function annotation.
    '''
    global q
    global n
    while 1:
        try:
            fam = q.get(block=True, timeout=0.1)
        except Empty: 
            if n.qsize() > 0:
                for _ in range(n.qsize()):
                    print n.get(),
                print ""
            break
        else:
            n.put(fam)
            if n.qsize() >= 10:
                for _ in range(10):
                    print n.get(),
                print ""
            os.system("java "+java_options+" -jar " \
    ""+options.sdir.replace(" ","\ ")+"sifter.jar "+fam+"" \
    " --output "+options.outdir.replace(" ","\ ")+"" \
    ""+fam.upper()+"/"+fam+".rdata --protein " \
    ""+options.outdir.replace(" ","\ ")+fam.upper()+"/" \
    ""+fam+".pli --reconciled " \
    ""+options.outdir.replace(" ","\ ")+fam.upper()+"/" \
    ""+fam+".nhx --scale " \
    ""+options.outdir.replace(" ","\ ")+fam.upper()+"" \
    "/scale-"+fam.lower()+".fx --familyfile " \
    ""+options.outdir.replace(" ","\ ")+fam.upper()+"" \
    "/infer-"+fam.lower()+".fx --alpha " \
    ""+options.outdir.replace(" ","\ ")+fam.upper()+"" \
    "/alpha-"+fam.lower()+".fx --truncation 1 --xvalidation --folds 20 -v " \
    "--ontology "+options.dbdir.replace(" ","\ ")+"" \
    "function.ontology "+scodes+" > " \
    ""+options.outdir.replace(" ","\ ")+fam.upper()+"/" \
    ""+fam+"-sifterout.txt ")
            q.task_done()


def sifter_multi(options, familylist, scodes):
    '''
    Run "run_sifter0" and "run_sifter1" on multiple threads.
    '''
    # Probability engine
    global q
    java_options = "-XX:+UseFastAccessorMethods -Xmx2048m"
    for fam in familylist:
        q.put(fam)
    print "# Generating required files...\n"
    for i in range(options.threads):
        p0 = Process(target=run_sifter0, name='%i' % (i+1), 
                      args = (options, java_options, scodes))
        p0.start()
    sleep(options.threads*0.05)
    q.join()
    sleep(options.threads*0.05)
    if p0.is_alive() and q.empty():
        sleep(options.threads*0.2)
        if p0.is_alive() and q.empty():
            p0.terminate()
    sleep(options.threads*0.05)
    numfunc_familylist = list()
    for fam in familylist:
        handle = open(options.outdir+fam.upper()+"/alpha-"+fam.lower()+".fx","r")
        numfunc_familylist.append((len(handle.readlines()), fam))
        handle.close()
    numfunc_familylist.sort()
    numfunc_familylist.reverse()
    for item in numfunc_familylist:
        q.put(item[1])
    global n

    print "# Running Sifter for the following families: \n"
    for i in range(options.threads):
        p1 = Process(target=run_sifter1, name='%i' % (i+1), 
                     args = (options, java_options, scodes))
        p1.start()
    sleep(options.threads*0.05)
    q.join()
    sleep(options.threads*0.05)
    if p1.is_alive() and q.empty():
        sleep(options.threads*0.2)
        if p1.is_alive() and q.empty():
            p1.terminate()


    
def _main():
    '''
    Main function for standalone usage.
    '''
    #Defines usage and help
    usage = "\n     %prog -d DIR"
    description = "This script automate the inference engine execution."
    #Defines input variables
    parser = OptionParser(usage=usage, 
        version="%prog 0.6.1", 
        description=description)
    parser.add_option("-d", "--directory", 
        dest="dir", 
        help="full path for the working directory.",     
        metavar="/DIRECTORY/DIR/")
    (options, args) = parser.parse_args()

    if len(args) > 0:
        print "\n# Extra arguments. Wrong usage. Exiting... \n"
        sys.exit(1)

    ### Exception care ###
    if not options.dir:
        print "Not all needed parameters were specified. " \
              "Type \"-h\" for help. \nExiting..."
        sys.exit(1)

    if not os.path.exists(options.dir+"options.pk"):
        print "\"options.pk\" not found.\nExiting..."
        sys.exit(1)

    options = pickle.load(file(options.dir+"options.pk", "r"))

    print "\n###  --------------------------------------------------  ###" 
    print "###           Posterior probability calculation          ###" 
    print "###  --------------------------------------------------  ###\n" 
    
    #Defines sifter arguments related to evidence codes
    scodes = "--with-iea --with-ic --with-iep --with-igc --with-igi " \
             "--with-ipi --with-iss --with-rca --with-tas --with-nas"

    familylist = list()
    if options.coverage:
        handle = open(options.outdir+"useful_pfam3.txt","r")
        for line in handle:
            d = line.strip().split()
            familylist.append(d[0])
        handle.close()
    else:
        handle = open(options.outdir+"useful_pfam.txt","r")
        for line in handle:
            d = line.strip().split()
            familylist.append(d[0])
        handle.close()

    sifter_multi(options, familylist, scodes)
    sys.exit()


if __name__ == '__main__':
    _main()

