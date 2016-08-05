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
#     Almeida-e-Silva D.C. and Vencio R.Z.N. (2015) SIFTER-T: A scalable    #
#     and optimized framework for the SIFTER phylogenomic method of         #
#     probabilistic protein domain annotation. BioTechniques, Vol. 58,      #
#     No. 3, March 2015, pp. 140-142                                        #
#                                                                           #
########## ########## ########### ########## ########## ########## ##########

"""
 * Gene tree preparation for SIFTER pipeline usage.

    clean_tree (options, fam, gene_sp)
    nseq(filename)
    pli2tree (options, gene_sp)
    gentree_multi(options, familylist, gene_sp)
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

def clean_tree (options, fam, gene_sp):
    ''' Clean_tree.py code '''
    if os.path.exists(options.outdir+fam.upper()+"/"+fam.upper()+".nex"):
        inf = open(options.outdir+fam.upper()+"/"+fam.upper()+".nex", 'r')
        outf = open(options.outdir+fam.upper()+"/"+fam.upper()+".nhx", 'w')
        species_set = set()
        for line in inf:
            d = line.strip().split(',')
            newline = ''
            for di in d:
                if len(newline) > 0:
                    outf.write(newline+',')
                while '/' in di:
                    di = di[0:di.find('/')]+di[di.find(':'):]
                if "query" in di:
                    spe = str(options.input_species)
                else:
                    spe_di = di.replace("(", "")
                    if options.type == "pf":
                        spe = "131567"
                    else:
                        spe = gene_sp[spe_di[:spe_di.find(":")]]
                    species_set.add(spe)
                if ')' in di:
                    direst = di[(di.find(")")+1):]
                    di = di[0:di.find(")")]+"[&&NHX:S="+spe+"])"
                    while ')' in direst:
                        di = di+direst[direst.find(":"):direst.find(")")]
                        if not options.reconciliation:
                            di = di+"[&&NHX:D=N]"
                        di = di+")"
                        direst = direst[(direst.find(")")+1):]
                    di = di+direst[direst.find(":"):]
                    if not options.reconciliation:
                        di = di+"[&&NHX:D=N]"
                else:
                    di = di+"[&&NHX:S="+spe+"]"
                newline = di
            outf.write(newline[0:newline.find(";")]+";\n")
        inf.close()
        outf.close()
        with open(options.outdir+fam.upper()+"/"+fam.upper()+".sp_list", "w") as handle:
            for item in species_set:
                handle.write(item+"\n")
    return None

def nseq(filename):
    '''
    Count number of sequences in a multifasta file.
    '''
    with open(filename, "r") as handle:
        n_seq = 0
        for line in handle:
            if line[0] == '>':
                n_seq += 1
    return n_seq


def pli2tree (options, gene_sp):
    '''
    Build trees and alignments.
    '''
    global q
    global n
    to_screen = set()
    while 1:
        try:
            fam = q.get(block=True, timeout=0.05)
        except: 
            for _ in range(len(to_screen)):
                n.put(to_screen.pop())
            break
        to_screen.add(fam)
        if len(to_screen) > 10:
            for _ in range(10):
                print to_screen.pop(),
            print ""
        if nseq(options.outdir+fam.upper()+"/aligned.fasta") > 1500:
            fastest = " -fastest "
        else:
            fastest = " "
        os.system(options.stdir.replace(" ","\ ")+"FastTree -nopr -quiet"  \
            ""+fastest+options.outdir.replace(" ","\ ")+fam.upper()+""   \
            "/aligned.fasta > "+options.outdir.replace(" ","\ ")+""        \
            ""+fam.upper()+"/"+fam.upper()+".nex")
        clean_tree(options, fam, gene_sp)
        q.task_done()
    return None


def gentree_multi(options, familylist, gene_sp):
    '''
    Run "pli2tree" on multiple threads.
    '''
    global q
    for fam in familylist:
        q.put(fam)
    global n
    print "# Building tree from multiple alignment for the following families: \n"
    for i in range(options.threads):
        p = Process(target=pli2tree, name='%i' % (i+1), args = (options, gene_sp))
        p.start()
    sleep(options.threads*0.05)
    q.join()
    while n.qsize() > 0:
        for _ in range(10):
            if n.qsize() > 0:
                print n.get(),
        print ""
    sleep(options.threads*0.05)
    if p.is_alive() and q.empty():
        sleep(options.threads*0.2)
        if p.is_alive() and q.empty():
            p.terminate()
    return None


def _main():
    '''
    Main function for standalone usage.
    '''
    #Defines usage and help
    usage = "\n    v %prog -d DIR"
    description = "Gene tree preparation for SIFTER pipeline usage."
    #Defines input variables
    parser = OptionParser(usage=usage, 
        version="%prog 0.2.2", 
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
        print "Not all needed parameters were specified. Type \"-h\" for help. \nExiting..."
        sys.exit(1)

    if not os.path.exists(options.dir+"options.pk"):
        print "\"options.pk\" not found.\nExiting..."
        sys.exit(1)

    options = pickle.load(file(options.dir+"options.pk", "r"))
    print "\n###  --------------------------------------------------  ###" 
    print "###                Building Gene Trees                   ###" 
    print "###  --------------------------------------------------  ###\n" 

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

    gene_sp = dict()
    with open(options.dbdir+"gene_sp.list", 'r') as handle:
        for line in handle:
            d = line.strip().split()
            gene_sp[d[0]] = d[1]

    gentree_multi(options, familylist, gene_sp)

    sys.exit()

    
if __name__ == '__main__':

    _main()

