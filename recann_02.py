#! /usr/bin/env python
# -*- coding: iso-8859-1 -*-  

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
 * Annotation recovery and \".pli\" files preparation for Sifter-T usage.

    pfam2pli (options, families, uniprotskeys_set, uniprots)
    get_forbidden_genes(options)
    get_uniprots(options)
    get_familylist(options)
    get_families(options, familylist)
    pfam2pli_multi(options)
    _main()

"""

from optparse import OptionParser
from multiprocessing import Process
from multiprocessing import JoinableQueue
from Queue import Empty
from time import sleep
import pickle
import os
import sys
from Bio import SeqIO
q = JoinableQueue() 

def pfam2pli (options, families, uniprotskeys_set, uniprots):
    '''
    Generates ".pli" files for each protein family.
    '''
    global q
    while 1:
        try:
            fam = q.get(block=True, timeout=0.1)
        except Empty: 
            break
        else:
            if not os.path.exists(options.outdir+fam.upper()+"/"+fam.upper()+".pli") or options.force:
                outf = open(options.outdir+fam.upper()+"/"+fam.upper()+".pli", 'w')
                outf.write("<?xml version=\"1.0\"?>\n")
                outf.write("<Family>\n")
                outf.write("  <FamilyID>"+fam+"</FamilyID>\n")
                for uni in families[fam]:
                    unip = uni
                    outf.write("  <Protein>\n")
                    if uni in uniprotskeys_set or uni[0:(uni.find('_'))] in uniprotskeys_set:
                        if uni[0:(uni.find('_'))] in uniprotskeys_set:
                            unip = uni[0:(uni.find('_'))]
                        outf.write("    <ProteinName>"+uni+"</ProteinName>\n")
                        outf.write("    <ProteinNumber>"+uni[0:(uni.find('_'))]+"</ProteinNumber>\n")
                        outf.write("    <GONumber>[")
                        for upfi in range(len(uniprots[unip])):
                            upf = uniprots[unip][upfi]
                            if upfi == 0:
                                outf.write(upf[0][3:])
                            else:
                                outf.write(", "+upf[0][3:])
                        outf.write("]</GONumber>\n")
                        outf.write("    <MOC>[")
                        for upfi in range(len(uniprots[unip])):
                            upf = uniprots[unip][upfi]
                            if upfi == 0:
                                outf.write(upf[1])
                            else:
                                outf.write(", "+upf[1])
                        outf.write("]</MOC>\n")
                    else:
                        outf.write("    <ProteinName>"+uni+"</ProteinName>\n")
                        outf.write("    <ProteinNumber>"+uni[0:(uni.find('_'))]+"</ProteinNumber>\n")
                    outf.write("  </Protein>\n")
                outf.write("</Family>")
                outf.close()
            q.task_done()


def get_forbidden_genes(options):
    '''
    Load set with the forbidden genes due to specie's annotation removal.
    '''
    forbidden_genes = set()
    if os.path.exists(options.outdir+"forbidden_sp_genes.txt"):
        handle = open(options.outdir+"forbidden_sp_genes.txt","r")
        for line in handle:
            forbidden_genes.add(line.strip().split()[0])
        handle.close()
    return forbidden_genes


def get_uniprots(options):
    #Sifter2.0 bug fix - evidence codes probability
    '''
    # Sifter2.0 accepts the following evidence codes:
    [EVCODE]	[probability a priori]
    "GOR"	0.9
    "IDA"	0.9
    "TAS"	0.9
    "IGI"	0.8*
    "IMP"	0.8
    "IPI"	0.8
    "E" 	0.6
    "IC"	0.4
    "IEP"	0.4
    "IGC"	0.4*
    "ISS"	0.4
    "RCA"	0.4
    "NAS"	0.3
    "ND"	0.3
    "NR"	0.3
    "P" 	0.2
    "IEA"	0.2

    * Bug - Misattribution of probabilities
    
    # Besides, Sifter2.0 accepts the following evidence codes:
    sift = ["E","GOR","IC","IDA","IEA","IEP","IGC","IGI","IMP","IPI","ISS",
            "NAS","ND","NR","P","RCA","TAS"]

    # However, according to GO, these are the full list of evidence codes:
    codes = ["IPI", "IDA", "ISS", "TAS", "IMP", "NAS", "ND", "IGI", "IC", "EXP",
             "IEP", "RCA", "IGC", "ISO", "ISA", "ISM", "IEA", "IBA", "IBD"]

    # We decided to cover the newer evidence codes changing their names for 
    # other evidence codes with appropriate probability
    codes - sift = ['ISM', 'IBA', 'IBD', 'ISO', 'EXP', 'ISA']


    New IEA prior probability attributions (Sifter2.0 source):
    IMP = 0.9
    IPI = 0.8
    ISS = 0.7
    IDA = 0.6
    IEP = 0.5
    TAS = 0.4
    NAS = 0.3
    RCA = 0.2
    IEA = 0.1
    
    Mappings


    '''
    ec_values = {}
    ec_values["0.9"] = "IMP"
    ec_values["0.8"] = "IPI"
    ec_values["0.7"] = "ISS"
    ec_values["0.6"] = "IDA"
    ec_values["0.5"] = "IEP"
    ec_values["0.4"] = "TAS"
    ec_values["0.3"] = "NAS"
    ec_values["0.2"] = "RCA"
    ec_values["0.1"] = "IEA"

    codes = ["IPI", "IDA", "ISS", "TAS", "IMP", "NAS", "ND", "IGI", "IC", "EXP",
             "IEP", "RCA", "IGC", "ISO", "ISA", "ISM", "IEA", "IBA", "IBD"]

    ec_mapping = {}
    handle = open(options.stdir+"prob_conversion.txt","r")
    for line in handle:
        d = line.strip().split()
        if d[1] not in ["0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9"]:
            print "Invalid value "+str(d[1])+"in \"prob_conversion.txt\" file. \nExiting... "
            sys.exit(1)
        else:
            ec_mapping[d[0]] = d[1]

    handle.close()

    ec_convert = {}
    for code in codes:
        ec_convert[code] = ec_values[ec_mapping[code]]


    # Loading GOA annotations
    uniprots = {}
    handle = open(options.dbdir+"summary_gene_association.goa_uniprot","r")
    forbidden_genes = get_forbidden_genes(options)
    for line in handle:
        d = line.strip().split()
        if d[2] in options.experimental and d[3] not in forbidden_genes:
            if d[2] in ec_convert:
                d[2] = ec_convert[d[2]]
            try:
                uniprots[d[3]].append((d[1], d[2]))
            except:
                uniprots[d[3]] = list()
                uniprots[d[3]].append((d[1], d[2]))
    handle.close()
    return uniprots


def get_uniprots_cov(options, familylist):
    #Sifter2.0 bug fix - evidence codes probability
    '''
    # Sifter2.0 accepts the following evidence codes:
    [EVCODE]	[probability a priori]
    "GOR"	0.9
    "IDA"	0.9
    "TAS"	0.9
    "IGI"	0.8*
    "IMP"	0.8
    "IPI"	0.8
    "E" 	0.6
    "IC"	0.4
    "IEP"	0.4
    "IGC"	0.4*
    "ISS"	0.4
    "RCA"	0.4
    "NAS"	0.3
    "ND"	0.3
    "NR"	0.3
    "P" 	0.2
    "IEA"	0.2

    * Bug - Misattribution of probabilities
    
    # Besides, Sifter2.0 accepts the following evidence codes:
    sift = ["E","GOR","IC","IDA","IEA","IEP","IGC","IGI","IMP","IPI","ISS",
            "NAS","ND","NR","P","RCA","TAS"]

    # However, according to GO, these are the full list of evidence codes:
    codes = ["IPI", "IDA", "ISS", "TAS", "IMP", "NAS", "ND", "IGI", "IC", "EXP",
             "IEP", "RCA", "IGC", "ISO", "ISA", "ISM", "IEA", "IBA", "IBD"]

    # We decided to cover the newer evidence codes changing their names for 
    # other evidence codes with appropriate probability
    codes - sift = ['ISM', 'IBA', 'IBD', 'ISO', 'EXP', 'ISA']


    New IEA prior probability attributions (Sifter2.0 source):
    IMP = 0.9
    IPI = 0.8
    ISS = 0.7
    IDA = 0.6
    IEP = 0.5
    TAS = 0.4
    NAS = 0.3
    RCA = 0.2
    IEA = 0.1
    
    Mappings


    '''
    ec_values = {}
    ec_values["0.9"] = "IMP"
    ec_values["0.8"] = "IPI"
    ec_values["0.7"] = "ISS"
    ec_values["0.6"] = "IDA"
    ec_values["0.5"] = "IEP"
    ec_values["0.4"] = "TAS"
    ec_values["0.3"] = "NAS"
    ec_values["0.2"] = "RCA"
    ec_values["0.1"] = "IEA"

    codes = ["IPI", "IDA", "ISS", "TAS", "IMP", "NAS", "ND", "IGI", "IC", "EXP",
             "IEP", "RCA", "IGC", "ISO", "ISA", "ISM", "IEA", "IBA", "IBD"]

    ec_mapping = {}
    handle = open(options.stdir+"prob_conversion.txt","r")
    for line in handle:
        d = line.strip().split()
        if d[1] not in ["0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9"]:
            print "Invalid value "+str(d[1])+"in \"prob_conversion.txt\" file. \nExiting... "
            sys.exit(1)
        else:
            ec_mapping[d[0]] = d[1]

    handle.close()

    ec_convert = {}
    for code in codes:
        ec_convert[code] = ec_values[ec_mapping[code]]

    # Load families genes
    full_genes_set = set()
    for fam in familylist:
        with open(options.dbdir+"align/gene_list/"+fam+".gene") as handle:
            for line in handle:
                d0 = line.strip().split()[0]
                d1 = d0[0:d0.find("/")]
                full_genes_set.add(d1)

    # Loading GOA annotations
    uniprots = {}
    handle = open(options.dbdir+"summary_gene_association.goa_uniprot","r")
    forbidden_genes = get_forbidden_genes(options)
    for line in handle:
        d = line.strip().split()
        if d[2] == "IEA" and d[3] not in forbidden_genes and d[3] in full_genes_set:
            try:
                uniprots[d[3]].append((d[1], d[2]))
            except:
                uniprots[d[3]] = list()
                uniprots[d[3]].append((d[1], d[2]))
    handle.close()
    return uniprots


def get_familylist(options):
    '''
    Load list of families to be treated.
    '''
    familylist = list()
    handle = open(options.outdir+"useful_pfam.txt","r")
    for line in handle:
        d = line.strip().split()
        familylist.append(d[0])
    handle.close()
    return familylist


def get_familylist2(options):
    '''
    Load list of families to be treated in extended coverage.
    '''
    familylist = list()
    handle = open(options.outdir+"useful_pfam2.txt","r")
    for line in handle:
        d = line.strip().split()
        familylist.append(d[0])
    handle.close()
    return familylist


def get_families(options, familylist):
    '''
    Load gene list for each selected family.
    '''
    families = {}
    for pf in familylist:
        families[pf] = set()
        handle = open(options.outdir+pf+"/aligned.fasta", "r")
        for nuc_rec in SeqIO.parse(handle, "fasta"):
            families[pf].add(nuc_rec.id[0:nuc_rec.id.find("/")])
        handle.close()
    return families


def pfam2pli_multi(options):
    '''
    Run "pfam2pli" on multiple threads.
    '''
    global q
    familylist = get_familylist(options)
    for pf in familylist:
        q.put(pf)
    if options.force:
        for pf in familylist:
            if os.path.exists(options.outdir+pf.upper()+"/"+pf.upper()+".pli"):
                os.remove(options.outdir+pf.upper()+"/"+pf.upper()+".pli")
    uniprots = get_uniprots(options)
    uniprotskeys_set = set(uniprots)

    print "# Building \".pli\" files for the following families: \n"
    n = 0
    for pf in familylist:
        print pf,
        n = n+1
        if n == 10:
            n = 0
            print ""

    families = get_families(options, familylist)
    #Start a pool of workers
    for i in range(options.threads):
        p = Process(target=pfam2pli, name='%i' % (i+1), 
            args = (options, families, uniprotskeys_set, uniprots))
        p.start()
    sleep(options.threads*0.05)
    q.join()
    sleep(options.threads*0.05)
    if p.is_alive() and q.empty():
        sleep(options.threads*0.2)
        if p.is_alive() and q.empty():
            p.terminate()
    if options.coverage:
        print ""
        familylist = get_familylist2(options)
        for pf in familylist:
            q.put(pf)
        if options.force:
            for pf in familylist:
                if os.path.exists(options.outdir+pf.upper()+"/"+pf.upper()+".pli"):
                    os.remove(options.outdir+pf.upper()+"/"+pf.upper()+".pli")
        del(uniprots)
        del(uniprotskeys_set)
        uniprots = get_uniprots_cov(options, familylist)
        uniprotskeys_set = set(uniprots)
    
        n = 0
        for pf in familylist:
            print pf,
            n = n+1
            if n == 10:
                n = 0
                print ""
    
        families = get_families(options, familylist)
        #Start a pool of workers
        for i in range(options.threads):
            p = Process(target=pfam2pli, name='%i' % (i+1), 
                args = (options, families, uniprotskeys_set, uniprots))
            p.start()
        sleep(options.threads*0.05)
        q.join()
        sleep(options.threads*0.05)
        if p.is_alive() and q.empty():
            sleep(options.threads*0.2)
            if p.is_alive() and q.empty():
                p.terminate()
        print "\n"



def _main():
    '''
    Main function for standalone usage.
    '''
    #Defines usage and help
    usage = "\n 	%prog -d DIR"
    description = "Annotation recovery and \".pli\" files preparation for " \
                  "SIFTER pipeline usage."
    
    #Defines input variables
    parser = OptionParser(usage=usage, version="%prog 0.3.0", 
        description=description)
    parser.add_option("-d", "--directory", 
        dest="dir", 
        help="full path for the working directory.", 
        metavar="/DIRECTORY/DIR/")
    parser.add_option("-f", "--force", 
        help="(Optional) Force file substitution. (default False)", 
        action="store_true", 
        dest="force", 
        default = False)
    (options, args) = parser.parse_args()

    if len(args) > 0:
        print "\n# Extra arguments. Wrong usage. Exiting... \n"
        sys.exit(1)

    ### Exception care ###
    if not options.dir:
        print "Not all needed parameters were specified. Type \"-h\" for help."\
              " \nExiting..."
        sys.exit(1)

    if not os.path.exists(options.dir+"options.pk"):
        print "\"options.pk\" not found.\nExiting..."
        sys.exit(1)

    force = options.force
    options = pickle.load(file(options.dir+"options.pk", "r"))
    options.force = force
    del(force)
    print "\n###  --------------------------------------------------  ###" 
    print "###                  Annotation Recovery                 ###" 
    print "###  --------------------------------------------------  ###\n" 
    pfam2pli_multi(options)
    sys.exit()

if __name__ == '__main__':
    _main()

