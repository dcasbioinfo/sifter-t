#! /usr/bin/env python
# -*- coding: iso-8859-1 -*-  

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
 * Protein Family multiple alignment generation for Sifter-T usage.

    nt_prepare_input(options)
    aa_prepare_input(options)
    pfam_scan_mp(options)
    split_query_fasta(options, n_seq)
    pfam_scan_multi(options)
    get_input_pfam_genes(options, num_files)
    get_pfam_genes(options, input_pfam_genes)
    get_annot_genes_all(options)
    get_no_annot(options, pfam_genes)
    get_no_annot2(options, pfam_genes, no_annot, annot_genes_all)
    get_sp_gene(options)
    get_sp_desc_anc(options)
    get_sp_branch_set(options, sp_desc_anc)
    get_forbidden_sp_gene(options, sp_branch_set, sp_gene, annot_genes_all)
    get_input_genes_pfam__handle_pf(options, useful_pfam, input_pfam_genes)
    write_input_ntaa(options, input_genes_pfam, handle_pf)
    write_fasta_pf(options, useful_pfam)
    write_selected_pfam_genes(options, annot_genes_all)
    multi_write_selected_pfam_genes(options, useful_pfam, annot_genes_all)
    clean_useful_pfam(options, useful_pfam)
    align_sequences(options, useful_pfam)
    write_pfam_list(options, useful_pfam)
    clean_annot_genes_all(annot_genes_all, pfam_genes, useful_pfam)
    _main()

"""

from optparse import OptionParser
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from time import sleep
import pickle
import os
import sys
import shutil
from multiprocessing import Process
from multiprocessing import JoinableQueue
from Queue import Empty
q = JoinableQueue()    
n = JoinableQueue()    

def nt_prepare_input(options):
    '''
    Converts nucleotide to 6 frames, translate each, convert names to
    Sifter-T intermediate names

    Intermediate fasta file: options.outdir+"query.fasta"
    Names conversion: options.outdir+"input_names.txt"
    '''
    if options.type == "nt":
        if not os.path.exists(options.outdir+"query.fasta") or options.force:
            handle_table = open(options.outdir+"input_names.txt","w")
            j = 0
            handle = open(options.file, "rU")
            out = open(options.outdir+"query.fasta", "w")
            for nuc_rec in SeqIO.parse(handle, "fasta"):
                if nuc_rec.id == nuc_rec.description:
                    nuc_rec.description = ""
                handle_table.write("query%s\t%s %s\n" % (str(j), 
                                              nuc_rec.id, nuc_rec.description))
                rev = nuc_rec.seq.reverse_complement()
                strand = "c"
                for i in range(0, 3):
                    SeqIO.write(SeqRecord(seq = nuc_rec.seq[(abs(i)):].translate(cds=False,
                        table=options.translation), 
                        id = "query%s_%s%s" % (str(j), 
                        strand, str(i+1)), 
                        description = ""), out, "fasta")
                strand = "w"
                for i in range(0, 3):
                    SeqIO.write(SeqRecord(seq = rev[(abs(i)):].translate(cds=False,
                        table=options.translation), 
                        id = "query%s_%s%s" % (str(j), strand, str(i+1)), 
                        description = ""), out, "fasta")
                j += 1
            handle.close()
            out.close()
            handle_table.close()


def aa_prepare_input(options):
    '''
    Convert names to Sifter-T intermediate names
 
    Intermediate fasta file: options.outdir+"query.fasta"
    Names conversion: options.outdir+"input_names.txt"
    '''
    if options.type == "aa":
        if not os.path.exists(options.outdir+"query.fasta") or options.force:
            handle_table = open(options.outdir+"input_names.txt","w")
            i = 0
            handle = open(options.file, "rU")
            out = open(options.outdir+"query.fasta", "w")
            for nuc_rec in SeqIO.parse(handle, "fasta"):
                handle_table.write("query%s\t%s %s\n" % (str(i), nuc_rec.id, 
                    nuc_rec.description))
                SeqIO.write(SeqRecord(seq = nuc_rec.seq, 
                    id = "query%s" % str(i), description = ""), out, "fasta")
                i += 1
            handle.close()
            out.close()
            handle_table.close()


def pfam_scan_mp(options):
    '''
    Run Pfam_Scan over input aminoacids or translated nucleotides. 
    '''
    global q
    global n
    while True:
        try:
            i = q.get(block=True, timeout=0.05)
        except Empty: 
            break
        else:
            os.system("perl " \
                ""+os.path.abspath(options.stdir).replace(" ","\ ")+"" \
                "/pfam_scan.pl -e_dom "+str(options.pcut)+" -e_seq " \
                ""+str(options.pcut)+" -cpu 4 -fasta " \
                ""+os.path.abspath(options.outdir).replace(" ","\ ")+"" \
                "/query_temp"+str(i)+".fasta -outfile " \
                ""+os.path.abspath(options.outdir).replace(" ","\ ")+"" \
                "/query_temp"+str(i)+".pfam -d " \
                ""+os.path.abspath(options.dbdir).replace(" ","\ "))
            q.task_done()


def split_query_fasta(options, n_seq):
    '''
    Split "query.fasta" for multithreading use of PfamScan
    '''
    handle_in = open(options.outdir+"query.fasta","r")
    n_threads = int(options.threads/2)
    n_seq_split = (n_seq/n_threads)+1
    i = 0
    n_seq_temp = 0
    handle_out = open(options.outdir+"query_temp"+str(i)+".fasta","w")
    for record in SeqIO.parse(handle_in, "fasta"):
        if n_seq_temp <= n_seq_split:
            SeqIO.write(SeqRecord(seq = record.seq, id = record.id, 
                description = record.description), handle_out, "fasta")
            n_seq_temp = n_seq_temp + 1
        else:
            handle_out.close()
            i = i + 1
            n_seq_temp = 0
            handle_out = open(options.outdir+"query_temp"+str(i)+".fasta","w")
            SeqIO.write(SeqRecord(seq = record.seq, id = record.id, 
                description = record.description), handle_out, "fasta")
            n_seq_temp = n_seq_temp + 1
    try:
        handle_out.close()
    except:
        pass
    handle_in.close()
    return i, n_threads


def pfam_scan_multi(options):
    '''
    Run pfam_scan to recover the domains/protein families associated to
    each sequence.
    '''
    print "# Running PfamScan over the input sequences...\n"
    handle = open(options.outdir+"query.fasta", "r")
    n_seq = 0
    for line in handle:
        if line[0] == ">":
            n_seq = n_seq+1
    handle.close()
    if options.force and os.path.exists(options.outdir+"query_temp.pfam"):
        os.remove(options.outdir+"query_temp.pfam")
    if n_seq < 100 or options.threads == 1:
        os.system("perl "+os.path.abspath(options.stdir).replace(" ","\ ")+"" \
            "/pfam_scan.pl -e_dom "+str(options.pcut)+" -e_seq " \
            ""+str(options.pcut)+" -cpu "+str(options.threads)+" -fasta " \
            ""+os.path.abspath(options.outdir).replace(" ","\ ")+"/query.fasta"\
            " -outfile "+os.path.abspath(options.outdir).replace(" ","\ ")+"" \
            "/query_temp.pfam -d " \
            ""+os.path.abspath(options.dbdir).replace(" ","\ "))
        return 1
    elif n_seq >= 100 or options.threads > 1:
        i, n_threads = split_query_fasta(options, n_seq)
        global q
        num_files = i+1
        if options.force:
            for j in range(num_files):
                if os.path.exists(options.outdir+"query_temp"+str(j)+".pfam"):
                    os.remove(options.outdir+"query_temp"+str(j)+".pfam")
        for j in range(num_files):
            q.put(j)

        sleep(options.threads*0.05)

        for j in range(n_threads):
            p = Process(target=pfam_scan_mp, name='%i' % (j+1), 
                args = (options,))
            p.start()

        sleep(options.threads*0.05)

        q.join()            

        sleep(options.threads*0.05)

        if p.is_alive() and q.empty():
            p.terminate()
        return num_files
    

def get_input_pfam_genes(options, num_files):
    '''
    Build input_pfam_genes, containing the families found and 
    (when aa or nt) the input genes belonging to each family
    '''
    if (options.type == "aa" or options.type == "nt") and n < 100:
        input_pfam_genes = {}
        out = open(options.outdir+"query.pfam","w")
        handle = open(options.outdir+"query_temp.pfam","r")
        for line in handle:
            d = line.strip().split()
            if len(d) > 0 and d[0][0:5] == "query":
                pf = d[5][0:d[5].find(".")]
                gene = d[0]
                start = d[1]
                end = d[2]
                out.write("%s\t%s\t%s\t%s\n" % (gene, start, end, pf))
                try:
                    input_pfam_genes[pf].add("%s|%s-%s" % (gene, start, end))
                except:
                    input_pfam_genes[pf] = set()
                    input_pfam_genes[pf].add("%s|%s-%s" % (gene, start, end))
        out.close()
        handle.close()
    elif (options.type == "aa" or options.type == "nt") and n >= 100:
        input_pfam_genes = {}
        out = open(options.outdir+"query.pfam","w")
        for j in range(num_files):
            handle = open(options.outdir+"query_temp"+str(j)+".pfam","r")
            for line in handle:
                d = line.strip().split()
                if len(d) > 0 and d[0][0:5] == "query":
                    pf = d[5][0:d[5].find(".")]
                    gene = d[0]
                    start = d[1]
                    end = d[2]
                    out.write("%s\t%s\t%s\t%s\n" % (gene, start, end, pf))
                    try:
                        input_pfam_genes[pf].add("%s|%s-%s" % (gene,
                                                               start, end))
                    except:
                        input_pfam_genes[pf] = set()
                        input_pfam_genes[pf].add("%s|%s-%s" % (gene, 
                                                               start, end))
            handle.close()
        out.close()
    elif options.type == "pf":
        input_pfam_genes = set()
        handle = open(options.file,"r")
        for line in handle:
            input_pfam_genes.add(line.strip().split()[0])
        handle.close()
    return input_pfam_genes


def get_pfam_genes(options, input_pfam_genes):
    '''
    Load selected protein familie's genes.
    '''
    print "# Loading Protein Familie's genes...\n"
    pfam_genes = {}
    for pf in input_pfam_genes:
        handle = open(options.dbdir+"align/gene_list/"+pf.upper()+".gene","r")
        pfam_genes[pf] = set()
        for line in handle:
            pfam_genes[pf].add(line.strip()[0:line.find("/")])
        handle.close()    
    return pfam_genes


def get_annot_genes_all(options):
    '''
    Load GOA annotations.
    '''
    print "# Loading GOA annotations...\n"
    annot_genes_all = set()
    handle = open(options.dbdir+"summary_gene_association.goa_uniprot","r")
    for line in handle:
        d = line.strip().split()
        if d[2] in options.experimental:
            annot_genes_all.add(d[3])
    handle.close()
    return annot_genes_all


def get_no_annot(options, pfam_genes):
    '''
    Find wich family have no annotations at all.
    '''
    print "# Checking for untractable families...\n"
    handle = open(options.dbdir+"pf_noannot.list","r")
    pf_noannot = set()
    for line in handle:
        pf_noannot.add(line.strip())
    handle.close() 
    no_annot = set(pfam_genes) & pf_noannot
    # write no_annot
    handle = open(options.outdir+"without_annotations_goa.txt","w")
    if len(no_annot) > 0:
        print "# Sifter-T will not treat the following families due to " \
              "lack of GOA annotations: \n"
        i = 0
        for pf in no_annot:
            print pf,
            handle.write(pf+"\n")
            i = i+1
            if i == 10:
                i = 0
                print ""
        print "\n"
    handle.close()
    return no_annot


def get_no_annot2(options, pfam_genes, no_annot, annot_genes_all):
    ''' 
    Find wich family have no annotations due to incomplete evidence
    codes selection.
    '''
    no_annot2 = set()
    for pf in (set(pfam_genes) - no_annot):
        if len(pfam_genes[pf] & annot_genes_all) == 0:
            no_annot2.add(pf)
    # write no_annot2
    handle = open(options.outdir+"without_annotations_ec.txt","w")
    if len(no_annot2) > 0:
        i = 0
        print "# Sifter-T will not treat the following families due to " \
              "incomplete selection of\n# evidence codes annotations: \n"
        for pf in no_annot2:
            print pf,
            handle.write(pf+"\n")
            i = i+1
            if i == 10:
                i = 0
                print ""
        print "\n"
    handle.close()
    return no_annot2


def get_sp_gene(options):
    '''
    Load dictionary with the corresponding species (NCBITaxID) for each gene for
    all selected families
    '''
    print "# Loading PFam genes and their species.\n"
    sp_gene = {}
    if len(options.branch) > 0 or len(options.species) > 0 or options.reconciliation:
        handle = open(options.dbdir+"gene_sp.list","r")
        for line in handle:
            d = line.strip().split()
            if len(d) > 1:
                try:
                    sp_gene[d[1]].add(d[0])
                except:
                    sp_gene[d[1]] = set()
                    sp_gene[d[1]].add(d[0])
        handle.close()
    return sp_gene


def get_sp_desc_anc(options):
    '''
    Load full specie's tree.
    '''
    print "# Loading species tree...\n"
    sp_desc_anc = dict()
    handle = open(options.dbdir+"summary_ncbi_taxonomy.obo","r")
    for line in handle:
        d = line.strip().split()
        sp_desc_anc[d[0]] = d[1]
    handle.close()

    handle = open(options.dbdir+"summary_taxonomy.txt","r")
    for line in handle:
        d = line.strip().split()
        if d[0] not in sp_desc_anc:
            sp_desc_anc[d[0]] = d[1]
    handle.close()
    return sp_desc_anc


def get_sp_branch_set(options, sp_desc_anc):
    ''' 
    Load full species tree branch.
    '''
    sp_branch_set = set()
    for sp in options.branch:
        sp_branch_set.add(str(sp))
    sp_branch_set_temp = set()
    for sp2 in sp_desc_anc:
        if sp_desc_anc[sp2] in set(options.branch):
            sp_branch_set_temp.add(sp2)
    while len(sp_branch_set_temp - sp_branch_set) > 0:
        sp_branch_set = sp_branch_set | sp_branch_set_temp
        sp_branch_set_temp = set()
        for sp2 in sp_desc_anc:
            if sp_desc_anc[sp2] in sp_branch_set:
                sp_branch_set_temp.add(sp2)
    return sp_branch_set


def get_forbidden_sp_gene(options, sp_branch_set, sp_gene, annot_genes_all):
    '''
    According to forbidden species annotation input, return a set with genes to 
    be removed from the annotated set.
    '''
    forbidden_sp_gene = "set() "
    i = 0
    temp_set = set()
    print "# Removing annotations from forbidden species.\n"
    for sp in ((set(options.species) | sp_branch_set | set(["-"])) & set(sp_gene)):
        if i >= 500:
            temp_set = (temp_set | eval(forbidden_sp_gene)) & annot_genes_all
            forbidden_sp_gene = "set() "
            i = 0
        forbidden_sp_gene = "%s| sp_gene[\"%s\"] " % (forbidden_sp_gene, sp)
        i = i + 1
    forbidden_sp_gene = (temp_set | eval(forbidden_sp_gene)) & annot_genes_all
    handle_sp = open(options.outdir+"forbidden_sp_gene.txt","w")
    for item in forbidden_sp_gene:
        handle_sp.write(item+"\n")
    handle_sp.close()
    return forbidden_sp_gene


def get_input_genes_pfam__handle_pf(options, useful_pfam, input_pfam_genes):
    '''
    Create a dir for each useful pfam and a variable with information about
    genes, recognized pfam and aminoacid region of a given pfam
    '''
    handle_pf = {}
    print "# Preparing fasta files... \n"
    for pf in useful_pfam:
        if not os.path.exists(options.outdir+pf):
            os.mkdir(options.outdir+pf)
        handle_pf[pf] = open(options.outdir+pf+"/input.fasta", "w")
        handle_pf[pf].close()
        handle_pf[pf] = options.outdir+pf+"/input.fasta"
    input_genes_pfam = {}
    for pf in useful_pfam:
        for seq in input_pfam_genes[pf]:
            gene = seq[0:seq.find("|")]
            input_genes_pfam[gene] = set()
    for pf in useful_pfam:
        for seq in input_pfam_genes[pf]:
            gene = seq[0:seq.find("|")]
            start = seq[seq.find("|")+1:seq.find("-")]
            end = seq[seq.find("-")+1:]
            input_genes_pfam[gene].add((pf, start, end))
    return input_genes_pfam, handle_pf


def write_input_ntaa(options, input_genes_pfam, handle_pf):
    '''
    For each Protein Family creates a file with input aminoacid sequences of 
    gene with regions identified as belonging to this Protein Family
    '''
    handle = open(options.outdir+"query.fasta","r")
    for nuc_rec in SeqIO.parse(handle, "fasta"):
        if nuc_rec.id.strip().split()[0] in input_genes_pfam:
            gene = nuc_rec.id.strip().split()[0]
            for item in input_genes_pfam[gene]:
                pf = item[0]
                start = int(item[1])
                end = int(item[2])
                handle2 = open(handle_pf[pf], "a")
                SeqIO.write(SeqRecord(seq = nuc_rec.seq[start:end], 
                    id = nuc_rec.id+"|"+str(start)+"-"+str(end), 
                    description = ""), handle2, "fasta")
                handle2.close()
    handle.close()
    handle_pf = {}


def write_fasta_pf(options, useful_pfam):
    '''
    Copy the full Protein Family multiple alignment to a Protein Family 
    subfolder in outdir.
    '''
    for pf in useful_pfam:
        if not os.path.exists(options.outdir+pf):
            os.mkdir(options.outdir+pf)
            shutil.copy(options.dbdir+"align/"+pf.upper()+".fasta", 
                        options.outdir+pf+"/"+pf+".fasta")


def write_selected_pfam_genes(options, annot_genes_all):
    '''
    For each Protein Family write a second multiple alignment file with just the
    annotated genes.
    '''
    global q 
    while 1:
        try:
            pf = q.get(block=True, timeout=0.1)
        except Empty: 
            break
        else:
            handle = open(options.dbdir+"align/"+pf.upper()+".fasta","r")
            handle_out = open(options.outdir+pf+"/"+pf.upper()+".fasta", "w")
            for nuc_rec in SeqIO.parse(handle, "fasta"):
                if nuc_rec.id[0:nuc_rec.id.find("/")] in annot_genes_all:
                    SeqIO.write(SeqRecord(seq = nuc_rec.seq, id = nuc_rec.id, 
                        description = ""), handle_out, "fasta")
            handle_out.close()
            handle.close()
            q.task_done()


def multi_write_selected_pfam_genes(options, useful_pfam, annot_genes_all):
    '''
    Run "write_selected_pfam_genes" on multiple threads. 
    '''
    global q
    q = JoinableQueue() 
    for fam in useful_pfam:
        q.put(fam)
    for i in range(options.threads):
        p = Process(target = write_selected_pfam_genes, name = '%i' % (i+1), 
                    args = (options, annot_genes_all))
        p.start()
    sleep(options.threads*0.05)
    q.join()
    sleep(options.threads*0.05)
    if p.is_alive() and q.empty():
        p.terminate()


def clean_useful_pfam(options, useful_pfam):
    ''' 
    Remove protein families without annotations after all previous removals.
    '''
    for pf in list(useful_pfam):
        i = 0
        handle = open(options.outdir+pf+"/"+pf.upper()+".fasta", "r")
        for line in handle:
            if line[0] == ">":
                i = i+1
        handle.close()
        if i == 0:
            useful_pfam.remove(pf)
    return useful_pfam


def align_sequences(options, useful_pfam):
    '''
    Align input sequences with the processed Protein Family.
    '''
    if options.type == "aa" or options.type == "nt":
        print "# Aligning sequences from the following protein families: \n"
    i = 0
    for pf in useful_pfam:
        if options.type == "aa" or options.type == "nt":
            os.system("mafft --auto --amino --anysymbol --quiet --thread " \
                "%s --add %s%s/input.fasta %s%s/%s.fasta > %s%s/aligned.fasta" \
                "" % (str(options.threads), options.outdir.replace(" ","\ "),
                      pf, options.outdir.replace(" ","\ "), pf, pf.upper(),
                      options.outdir.replace(" ","\ "), pf))
        elif options.type == "pf":
            shutil.move(options.outdir+pf+"/"+pf+".fasta", 
                options.outdir+pf+"/aligned.fasta")
        print pf,
        i = i + 1
        if i == 10:
            i = 0
            print ""
    print "\n"


def write_pfam_list(options, useful_pfam):
    '''
    Sort useful_pfam from largest to smallest family and store in 
    "useful_pfam.txt".
    '''
    useful_pfam_list = list()
    for pf in useful_pfam:
        i = 0
        handle = open(options.outdir+pf+"/aligned.fasta", "r")
        for line in handle:
            if line[0] == ">":
                i = i + 1
        handle.close()
        if i > 0:
            useful_pfam_list.append((i, pf))

    useful_pfam_list.sort()
    useful_pfam_list.reverse()

    handle = open(options.outdir+"useful_pfam.txt","w")
    for item in useful_pfam_list:
        handle.write(item[1]+"\n")
    handle.close()


def clean_annot_genes_all(annot_genes_all, pfam_genes, useful_pfam):
    '''
    Reduce "annot_genes_all" dimension.
    '''
    annot_genes_all_clean = "set() "
    i = 0
    temp_set = set()
    for pf in useful_pfam:
        if i >= 250:
            temp_set = temp_set | eval(annot_genes_all_clean)
            annot_genes_all_clean = "set() "
            i = 0
        annot_genes_all_clean = "%s| (annot_genes_all & pfam_genes[\"%s\"]) " % (annot_genes_all_clean, pf)
        i = i + 1
    annot_genes_all_clean = temp_set | eval(annot_genes_all_clean)

    return annot_genes_all_clean



def _main():
    ''' 
    Main function for standalone usage.
    '''
    #Defines usage and help
    usage = "\n 	%prog -d DIR"
    description = "Protein Family multiple alignment generation for " \
                  "Sifter-T usage."
    #Defines input variables
    parser = OptionParser(usage=usage, version="%prog 0.7.1", 
        description=description)
    parser.add_option("-d", "--directory", 
        dest="dir", help="full path for the working directory.", 
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
              "\nExiting..."
        sys.exit(1)
    
    if not os.path.exists(options.dir+"options.pk"):
        print "\"options.pk\" not found.\nExiting..."
        sys.exit(1)

    force = options.force
    options = pickle.load(file(options.dir+"options.pk", "r"))
    options.force = force
    del(force)

    print "\n###  --------------------------------------------------  ###" 
    print "###      Sequence and Multiple Alignment Preparation     ###" 
    print "###  --------------------------------------------------  ###\n" 

    if options.type == "aa" or options.type == "nt":
        print "# Preparing the input sequences...\n"
        if options.type == "nt":
            nt_prepare_input(options)
        elif options.type == "aa":
            aa_prepare_input(options)
        num_files = pfam_scan_multi(options)
    elif options.type == "pf":
        num_files = 1

    input_pfam_genes = get_input_pfam_genes(options, num_files)

    no_annot = get_no_annot(options, set(input_pfam_genes))

    pfam_genes = get_pfam_genes(options, (set(input_pfam_genes) - no_annot))

    annot_genes_all = get_annot_genes_all(options)

    annot_genes_all = annot_genes_all - get_forbidden_sp_gene(options, 
                           get_sp_branch_set(options, get_sp_desc_anc(options)),
                           get_sp_gene(options), annot_genes_all)

    no_annot2 = get_no_annot2(options, pfam_genes, no_annot, annot_genes_all)

    useful_pfam = (set(input_pfam_genes) - (no_annot | no_annot2))

    annot_genes_all = clean_annot_genes_all(annot_genes_all, pfam_genes, useful_pfam)

    del(no_annot)
    del(no_annot2)
    del(pfam_genes)

    if options.type == "nt" or options.type == "aa":
        input_genes_pfam, handle_pf = get_input_genes_pfam__handle_pf(options,
                                                  useful_pfam, input_pfam_genes)
        write_input_ntaa(options, input_genes_pfam, handle_pf)
        del(input_genes_pfam)
        del(handle_pf)
    elif options.type == "pf":
        write_fasta_pf(options, useful_pfam)
    multi_write_selected_pfam_genes(options, useful_pfam, annot_genes_all)
    useful_pfam = clean_useful_pfam(options, useful_pfam)
    align_sequences(options, useful_pfam)
    write_pfam_list(options, useful_pfam)
    sys.exit()

if __name__ == '__main__':
    _main()

