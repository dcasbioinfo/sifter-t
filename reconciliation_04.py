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
 * Species tree creation and Protein Family reconciliation for Protein Families 
 * with suitable size.

    get_species_set(options, fam, desc_anc)
    buid_species_tree(options, fam, desc_anc, anc_desc)
    build_reconciled_tree(options)
    build_rec_multi(options, familylist)
    _main()

"""


from sys import stdout
try: 
    import dendropy
except:
    raise Exception("Can't load Dendropy. Exiting...")
from optparse import OptionParser
from multiprocessing import Process
from multiprocessing import JoinableQueue
from Queue import Empty
from time import sleep
import pickle
import os
import shutil
import sys
from Bio import SeqIO
q = JoinableQueue() 
n = JoinableQueue() 


def get_species_set(options, fam, desc_anc):
    '''
    Given a set of selected species of interest, return the set of the selected
    species of interest plus their ancestors.
    '''
    fam_sp = set()

    with open(options.outdir+fam.upper()+"/"+fam.upper()+".sp_list") as handle:
        for line in handle:
            fam_sp.add(line.strip())

    fam_sp = fam_sp | set([str(options.input_species)])

    temp_sp_set = fam_sp.copy()
    last_len_sp_set = 0
    
    sp_not = set()
    
    while last_len_sp_set != len(temp_sp_set):
        last_len_sp_set = len(temp_sp_set)
        for sp in set(temp_sp_set):
            try:
                temp_sp_set.add(desc_anc[sp])
            except:
                if not sp in desc_anc:
                    sp_not.add(sp)
    return temp_sp_set


def buid_species_tree(options, fam, desc_anc, anc_desc):
    '''
    Builds the species tree for a given Protein Family.
    '''
    print fam,
    tree1 = dendropy.Tree()
    
    temp_dict = dict()
    
    all_nodes = get_species_set(options, fam, desc_anc)

    for item in all_nodes:
        temp_dict[item] = {}
    
    for item in all_nodes:
        temp_dict[item]["oid"] = tree1.taxon_set.require_taxon(oid=item)
        temp_dict[item]["taxon"] = tree1.taxon_set.require_taxon(label=item)

    temp_dict["1"]["taxon"] = tree1.taxon_set.require_taxon(label="1")
    temp_dict["1"]["node"] = tree1.seed_node.new_child(taxon=temp_dict["1"]["taxon"])
    temp_dict["1"]["node"] = tree1.seed_node.new_child(taxon=temp_dict["1"]["oid"])

    i = 0
    all_nodes = all_nodes - set(["-"])
    sp_level = set(["1"])
    past_sp = set()
    x = (all_nodes - past_sp)

    while len(all_nodes - past_sp) > 0:
        temp_sp_level = set()
        for sp in set(sp_level):
            past_sp.add(sp)
            if sp in anc_desc:
                for desc in ((anc_desc[sp] & all_nodes) - past_sp):
                    if "taxon" in temp_dict[desc] and "oid" in temp_dict[desc]:
                        temp_dict[desc]["node"] = temp_dict[sp]["node"].new_child(taxon=temp_dict[desc]["taxon"])
                        temp_dict[desc]["node"] = temp_dict[sp]["node"].new_child(taxon=temp_dict[desc]["oid"])
                    else:
                        temp_dict[desc]["node"] = temp_dict[sp]["node"].new_child(taxon=temp_dict[desc]["oid"])
                    temp_sp_level.add(desc)
        sp_level = temp_sp_level.copy()
        if len(all_nodes -past_sp) > 0 and x == set(all_nodes -past_sp):
            i += 1
            if i > 5: 
                print "Problem in Protein Family", fam, x
                break
        x = (all_nodes - past_sp)
    
    tree1.reroot_at_node(temp_dict["1"]["node"])
    tree1.resolve_polytomies()
    tree1.update_splits(delete_outdegree_one=False)
    
    with open(options.outdir+fam.upper()+"/"+fam.upper()+".sptree", "w") as handle:
        handle.write(tree1.as_string("newick", suppress_edge_lengths=True))
    
    
def build_reconciled_tree(options):
    '''
    Runs Notung to reconciliate all Protein Families.
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
        os.system("java -jar "+options.stdir.replace(" ","\ ")+"notung/"       \
                  "Notung.jar -s "+options.outdir.replace(" ","\ ")+""     \
                  ""+fam.upper()+"/"+fam.upper()+".sptree -g "                 \
                  ""+options.outdir.replace(" ","\ ")+fam.upper()+"/"          \
                  ""+fam.upper()+".nhx --reconcile --speciestag nhx "          \
                  "--treeoutput nhx --nolosses --usegenedir --silent")
        shutil.move(options.outdir+fam.upper()+"/"+fam.upper()+".nhx",
            options.outdir+fam.upper()+"/"+fam.upper()+".nhx.notreconciled")
        if os.path.exists(options.outdir+fam.upper()+"/"+fam.upper()+".nhx.reconciled.0"):
            shutil.move(options.outdir+fam.upper()+"/"+fam.upper()+".nhx.reconciled.0",
                options.outdir+fam.upper()+"/"+fam.upper()+".nhx")
        elif os.path.exists(options.outdir+fam.upper()+"/"+fam.upper()+".nhx.reconciled"):
            shutil.move(options.outdir+fam.upper()+"/"+fam.upper()+".nhx.reconciled",
                options.outdir+fam.upper()+"/"+fam.upper()+".nhx")
        q.task_done()
    return None


def build_rec_multi(options, familylist):
    '''
    Run "build_reconciled_tree" on multiple threads.
    '''
    global q
    for fam in familylist:
        q.put(fam)
    global n

    for i in range(options.threads):
        p = Process(target=build_reconciled_tree, name='%i' % (i+1), args = (options,))
        p.start()
    sleep(options.threads*0.05)
    q.join()
    while n.qsize() > 0:
        for _ in range(10):
            if n.qsize() > 0:
                print n.get(),
        print ""
    sleep(options.threads*0.1)
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
    usage = "\n 	%prog -d DIR"
    description = "Species tree creation and Protein Family reconciliation for"\
                  " Protein Families with suitable size."
    
    #Defines input variables
    parser = OptionParser(usage=usage, version="%prog 0.2.0", 
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

    if not options.reconciliation:
        sys.exit()

    print "\n###  --------------------------------------------------  ###" 
    print "###                  Building Species Tree               ###" 
    print "###  --------------------------------------------------  ###\n" 

    desc_anc = dict()
    with open(options.dbdir+"summary_ncbi_taxonomy.obo") as handle:
        for line in handle:
            d = line.strip().split()
            desc_anc[d[0]] = d[1]

    desc_anc_alt = dict()
    with open(options.dbdir+"summary_taxonomy.txt") as handle:
        for line in handle:
            d = line.strip().split()
            if d[0] in desc_anc:
                if desc_anc[d[0]] != d[1]:
                    desc_anc_alt[d[0]] = d[1]
            else:
                desc_anc[d[0]] = d[1]

    anc_desc = dict()
    for item in desc_anc.values():
        anc_desc[item] = set()
    
    for item in desc_anc:
        anc_desc[desc_anc[item]].add(item)

    familylist = list()
    with open(options.outdir+"useful_pfam.txt","r") as handle:
        for line in handle:
            d = line.strip().split()
            familylist.append(d[0])

    i = 0
    for fam in familylist:
        buid_species_tree(options, fam, desc_anc, anc_desc)
        i +=1
        if i >= 10:
            print ""
            i = 0
    if i > 0:
            print ""
            i = 0

    print "\n###  --------------------------------------------------  ###" 
    print "###                  Reconciling Trees                   ###" 
    print "###  --------------------------------------------------  ###\n" 

    build_rec_multi(options, familylist)

    sys.exit()


if __name__ == '__main__':
    _main()

