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
 * This script prepare input databases for Sifter-T usage.

    remove_former_db(options)
    recover_stockholm(options)
    write_pfamlist(options, fullfamilyset)
    get_memsize()
    stockfastaconvert (options)
    stockfasta_multi(options, fullfamilyset)
    get_acession_numbers(options)
    get_fullgeneset(gene_ac)
    write_gene_ac(options, gene_ac)
    get_ac_gene(gene_ac)
    summary_goa(options,fullgeneset,ac_gene)
    get_annot_gene_set(options)
    write_annot_gene_list(options,annot_gene_set)
    write_gene_sp(options, annot_gene_set)
    write_pf_noannot(options, annot_gene_set, fullfamilyset)
    summ_tax1(options)
    is_number(s)
    summ_tax2(options)
    obo2_functionontology(options)
    write_go_names(options)
    prepare_hmm(options)
    _main()

"""

try:
    from Bio import SeqIO
except:
    raise Exception("Can't load BioPython. Exiting...")
from optparse import OptionParser
import os
import sys
import gzip
import zlib
import bz2
from multiprocessing import JoinableQueue
from multiprocessing import Process
import multiprocessing
from Queue import Empty
from time import sleep
q = JoinableQueue()
n = JoinableQueue()


def remove_former_db(options):
    '''
    Clean old files in case of "force".
    '''
    if options.force:
        if os.path.exists(options.dbdir+"summary_gene_association.goa_uniprot"):
            os.remove(options.dbdir+"summary_gene_association.goa_uniprot")
        if os.path.exists(options.dbdir+"align"):
            import shutil
            shutil.rmtree(os.path.abspath(options.dbdir+"align"))
        if os.path.exists(options.dbdir+"Pfam-A.hmm.h3f"):
            os.remove(options.dbdir+"Pfam-A.hmm.h3f")
        if os.path.exists(options.dbdir+"Pfam-A.hmm.h3i"):
            os.remove(options.dbdir+"Pfam-A.hmm.h3i")
        if os.path.exists(options.dbdir+"Pfam-A.hmm.h3m"):
            os.remove(options.dbdir+"Pfam-A.hmm.h3m")
        if os.path.exists(options.dbdir+"Pfam-A.hmm.h3p"):
            os.remove(options.dbdir+"Pfam-A.hmm.h3p")
        if os.path.exists(options.dbdir+"function.ontology"):
            os.remove(options.dbdir+"function.ontology")
        if os.path.exists(options.dbdir+"annot_gene.list"):
            os.remove(options.dbdir+"annot_gene.list")
        if os.path.exists(options.dbdir+"gene_sp.list"):
            os.remove(options.dbdir+"gene_sp.list")
        if os.path.exists(options.dbdir+"pfam.list"):
            os.remove(options.dbdir+"pfam.list")
        if os.path.exists(options.dbdir+"pfam_gene_ac.list"):
            os.remove(options.dbdir+"pfam_gene_ac.list")
        if os.path.exists(options.dbdir+"summary_taxonomy.txt"):
            os.remove(options.dbdir+"summary_taxonomy.txt")
        if os.path.exists(options.dbdir+"summary_ncbi_taxonomy.obo"):
            os.remove(options.dbdir+"summary_ncbi_taxonomy.obo")


def recover_stockholm(options):
    '''
    Multiple alignment recover for all protein families - Stockholm format
    '''
    print "# Recovering stockholm multiple alignments from PFam Database...\n"
    if not os.path.exists(options.dbdir+"align"):
        os.mkdir(options.dbdir+"align")
    if not os.path.exists(options.dbdir+"align/stockholm"):
        os.mkdir(options.dbdir+"align/stockholm")
    if not os.path.exists(options.dbdir+"align/gene_list"):
        os.mkdir(options.dbdir+"align/gene_list")

    infpf = open(options.dbdir+"Pfam-A.full", 'r')
    outf = open(options.dbdir+"temp.fasta", 'w')
    printout = False
    pfamily = ''
    i = 0
    fullfamilyset = set()
    for line in infpf:
        if '#=GF AC' in line:
            pfamily = line[line.find('PF'):line.find('.')]
            fullfamilyset.add(pfamily)
            if not os.path.exists(options.dbdir+"align/stockholm/"+pfamily+"." \
                    "stockholm") and not \
                    os.path.exists(options.dbdir+"align/"+pfamily+".fasta") \
                    or options.force:
                print pfamily,
                i += 1
                if i == 10:
                    print ""
                    i = 0
                outf.close()
                outf = open(options.dbdir+"align/stockholm/"+pfamily+".stockholm", 'w')
                outf.write("# STOCKHOLM 1.0\n#=GF ID   "+pfamily+"\n")
                printout = True
        if line[0:2] == '//':
            printout = False
        elif printout:
            outf.write(line)
    outf.close()
    infpf.close()
    if i > 0:
        print ""
    print ""
    return fullfamilyset


def write_pfamlist(options, fullfamilyset):
    '''
    Write family list.
    '''
    print "# Loading family list...\n"
    handle = open(options.dbdir+"pfam.list", "w")
    fullfamilylist = list(fullfamilyset)
    fullfamilylist.sort()
    for item in fullfamilylist:
        handle.write(item+"\n")
    handle.close()
    del fullfamilylist


def get_memsize():
    '''
    Recover system memory size.
    '''
    handle = open("/proc/meminfo", "r")
    memsize = int(handle.readlines()[0].strip().split()[1])
    handle.close()
    return memsize


def stockfastaconvert(options):
    '''
    Conversion of multiple alignment - Stockholm to Multifasta
    '''
    while True:
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
            if os.path.exists(options.dbdir+"align/stockholm/"+fam+".stockholm"):
                if os.path.exists(options.dbdir+"align/"+fam+".fasta"):
                    os.remove(options.dbdir+"align/"+fam+".fasta", "fasta")
                SeqIO.convert(options.dbdir+"align/stockholm/"+fam+".stockholm",
                    "stockholm", options.dbdir+"align/"+fam+".fasta", "fasta")
                handle = open(options.dbdir+"align/"+fam.upper()+".fasta", "r")
                temp = set()
                for nuc_rec in SeqIO.parse(handle, "fasta"):
                    temp.add(nuc_rec.id)
                handle.close()
                handle = open(options.dbdir+"align/gene_list/"+fam.upper()+".gene", "w")
                for gene in temp:
                    handle.write(gene+"\n")
                handle.close()
                os.remove(options.dbdir+"align/stockholm/"+fam+".stockholm")
            q.task_done()


def stockfasta_multi(options, fullfamilyset):
    '''
    Conversion of multiple alignment - Stockholm to Multifasta.
    '''
    print "# Conversion of multiple alignment - Stockholm to Multifasta\n"
    fam_size = list()
    for pfam in fullfamilyset:
        if os.path.exists(options.dbdir+"align/stockholm/"+pfam+".stockholm"):
            fam_size.append((os.path.getsize(options.dbdir+"align/stockholm/"+pfam+".stockholm"), pfam))
    fam_size.sort()
    fam_size.reverse()

    for item in fam_size:
        if os.path.exists(options.dbdir+"align/stockholm/"+item[1]+".stockholm") \
        and not os.path.exists(options.dbdir+"align/"+item[1]+".fasta"):
            q.put(item[1])
    lastline = True
    if q.qsize() == 0:
        lastline = False

    memsize = get_memsize()
    if (memsize/1572864.) <= 2.:
        num_proc = 1
    elif (memsize/1572864.) < options.threads:
        num_proc = int(memsize/1572864)
    else:
        num_proc = options.threads

    for i in range(num_proc):
        proc = Process(target=stockfastaconvert, name='%i' % (i+1), args=(options,))
        proc.start()

    sleep(options.threads*0.2)
    q.join()
    sleep(options.threads*0.05)
    if proc.is_alive() and q.empty():
        sleep(options.threads*0.2)
        if proc.is_alive() and q.empty():
            proc.terminate()
    if lastline:
        print ""
    del fam_size


def get_acession_numbers(options):
    '''
    Get genes and accession numbers from Pfam-A.full
    '''
    print "# Loading genes and accession numbers (1st phase)...\n"
    #gene_ac = {}
    pf = ""
    infpf = open(options.dbdir+"Pfam-A.full", 'r')
    handle = open(options.dbdir+"gene.list", "w")
    for line in infpf:
        d = line.split()
        if len(d) > 1:
            if d[0] == '#=GS':
                gene = d[1][0:d[1].find("/")]
                if d[2] == "AC":
                    ac = d[3][0:d[3].find(".")]
                handle.write("%s\n" % (gene))
                gene = ""
                ac = ""
            elif d[1] == "AC" and d[0] == '#=GF':
                pf = d[2][0:d[2].find(".")]
                #gene_ac[pf] = {}
    infpf.close()
    handle.close()
    return None


def get_fullgeneset(options):
    '''
    Get the complete set of genes present on Pfam database.
    '''
    print "# Loading genes and accession numbers (2nd phase)...\n"
    handle = open(options.dbdir+"gene.list", "r")
    fullgeneset = set()
    for line in handle:
        d = line.strip().split("\t")
        if len(d) >= 1:
            fullgeneset.add(d[0])
    return fullgeneset


def write_gene_ac(options, gene_ac):
    '''
    Write gene and their acession numbers to "pfam_gene_ac.list".
    '''
    handle = open(options.dbdir+"pfam_gene_ac.list", "w")
    for pf in gene_ac:
        for gene in gene_ac[pf]:
            handle.write("%s\t%s\t%s\n" % (pf, gene, gene_ac[pf][gene]))
    handle.close()


def get_ac_gene(options):
    '''
    Get acession numbers correspondence to their genes.
    '''
    ac_gene = {}
    handle = open(options.dbdir+"pfam_gene_ac.list", "r")
    for line in handle:
        d = line.strip().split("\t")
        if len(d) > 2:
            ac = d[2]
            gene = d[1]
            ac_gene[ac] = gene
    return ac_gene


def summary_goa(options, fullgeneset, ac_gene):
    '''
    Generate a file with anotations of all types for genes in PFAM
    '''
    from sys import stdout
    print "# Summarizing GOA database...\n"
    forbidden_ec = set(["ND", "NR", "IKR", "IRD"])
    not_ec = set(["IKR", "IRD"])
    i = 0
    fullgeneset = get_fullgeneset(options)
    if not os.path.exists(options.dbdir+"summary_gene_association.goa_uniprot") or options.force:
        print "# Preparing intermediate files from \"gene_association.goa_uniprot\"...\n"
        infgo = open(options.dbdir+"gene_association.goa_uniprot", 'r')
        out = open(options.dbdir+"summary_gene_association.goa_uniprot_temp", 'w')
        outnot = open(options.dbdir+"summary_gene_association.goa_uniprot_not", 'w')
        for line in infgo:
            d = line.strip().split('\t')
            if len(d) > 10 and d[8] == 'F' and d[0][0:3] == 'Uni':
                if d[6] not in forbidden_ec:
                    if d[10][(d[10].rfind('|')+1):len(d[10])] in fullgeneset:
                        out.write("%s\t%s\t%s\t%s\n" %
                        (d[1], d[4], d[6], d[10][(d[10].rfind('|')+1):len(d[10])]))
                elif d[6] in not_ec:
                    if d[10][(d[10].rfind('|')+1):len(d[10])] in fullgeneset:
                        outnot.write("%s\t%s\t%s\t%s\n" %
                        (d[1], d[4], d[6], d[10][(d[10].rfind('|')+1):len(d[10])]))
                i = i + 1
                if i >= 100:
                    stdout.write("\r%s          " % d[1])
                    stdout.flush()
                    i = 0
        del fullgeneset
        infgo.close()
        out.close()
        outnot.close()
        stdout.write("\r")
        stdout.flush()
        print "# Cleaning not-annotations (IKR and IRD evidence codes) \n"
        goa_not_handle = open(options.dbdir+"summary_gene_association.goa_uniprot_not", 'r')
        goa_not_dict = {}
        for line in goa_not_handle:
            d = line.strip().split('\t')
            try:
                goa_not_dict[d[3]].add(d[1])
            except:
                goa_not_dict[d[3]] = set()
                goa_not_dict[d[3]].add(d[1])
        goa_not_handle.close()
        goa_handle = open(options.dbdir+"summary_gene_association.goa_uniprot", 'w')
        goa_temp_handle = open(options.dbdir+"summary_gene_association.goa_uniprot_temp", 'r')
        counter_not = 0
        for line in goa_temp_handle:
            d = line.strip().split('\t')
            if not d[3] in goa_not_dict:
                goa_handle.write(line)
            elif not d[1] in goa_not_dict[d[3]]:
                goa_handle.write(line)
            else:
                counter_not += 1
        goa_handle.close()
        goa_temp_handle.close()
        if os.path.exists(options.dbdir+"summary_gene_association.goa_uniprot_temp"):
            os.remove(options.dbdir+"summary_gene_association.goa_uniprot_temp")
        print counter_not, "not-annotations were removed.\n"
    print ""


def get_annot_gene_set(options):
    '''
    Load a list of all pfam annotated genes.
    '''
    print "# Storing annotated genes list...\n"
    annot_gene_set = set()
    infgo = open(options.dbdir+"summary_gene_association.goa_uniprot", 'r')
    for line in infgo:
        d = line.strip().split()
        annot_gene_set.add(d[3])
    infgo.close()
    return annot_gene_set


def write_annot_gene_list(options, annot_gene_set):
    '''
    Store a list with all annotated genes in Pfam.
    '''
    annot_gene_list = list(annot_gene_set)
    annot_gene_list.sort()
    handle = open(options.dbdir+"annot_gene.list", "w")
    for gene in annot_gene_list:
        handle.write("%s\n" % gene)
    handle.close()
    del annot_gene_list


def write_gene_sp(options, annot_gene_set):
    '''
    Recover the species for each gene in Pfam.
    '''
    #Preparing gene's Species information
    sp_gene = {}
    if not os.path.exists(options.dbdir+"summary_uniprot_trembl.dat") or options.force:
        print "# Preparing gene's Species (from uniprot_trembl.dat) information...\n"
        uniprot_trembl = open(options.dbdir+"uniprot_trembl.dat", 'r')
        x = {}
        x["ID"] = ""
        x["OX"] = ""
        while 1:
            d = uniprot_trembl.readline().split()
            if len(d) > 1:
                if d[0] == "ID" and d[1] in annot_gene_set:
                    x["ID"] = d[1]
                elif d[0] == "OX" and x["ID"] in annot_gene_set:
                    x["OX"] = d[1][d[1].rfind('=')+1:len(d[1])-1]
                    try:
                        sp_gene[x["OX"]].add(x["ID"])
                    except:
                        sp_gene[x["OX"]] = set()
                        sp_gene[x["OX"]].add(x["ID"])
                    x = {}
                    x["ID"] = ""
                    x["OX"] = ""
                else:
                    continue
            elif d == [] or d == [""]:
                break
        uniprot_trembl.close()
    if not os.path.exists(options.dbdir+"summary_uniprot_sprot.dat") or options.force:
        print "# Preparing gene's Species (from uniprot_sprot.dat) information...\n"
        uniprot_sprot = open(options.dbdir+"uniprot_sprot.dat", 'r')
        x = {}
        x["ID"] = ""
        x["OX"] = ""
        while 1:
            d = uniprot_sprot.readline().split()
            if len(d) > 1:
                if d[0] == "ID" and d[1] in annot_gene_set:
                    x["ID"] = d[1]
                elif d[0] == "OX" and x["ID"] in annot_gene_set:
                    x["OX"] = d[1][d[1].rfind('=')+1:len(d[1])-1]
                    try:
                        sp_gene[x["OX"]].add(x["ID"])
                    except:
                        sp_gene[x["OX"]] = set()
                        sp_gene[x["OX"]].add(x["ID"])
                    x = {}
                    x["ID"] = ""
                    x["OX"] = ""
                else:
                    continue
            elif d == [] or d == [""]:
                break
        uniprot_sprot.close()
    # Store gene_sp
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

    sp_conversion = {}
    wrong_values = (set(desc_anc.values()) - set(desc_anc))
    for item in set(desc_anc):
        if desc_anc[item] in wrong_values:
            sp_conversion[desc_anc[item]] = item


    delnodes = set()
    if os.path.exists(options.dbdir+"delnodes.dmp"):
        with open(options.dbdir+"delnodes.dmp") as handle_delnodes:
            for line in handle_delnodes:
                d = line.strip().split()
                delnodes.add(d[0])

    merged_nodes = {}
    if os.path.exists(options.dbdir+"merged.dmp"):
        with open(options.dbdir+"merged.dmp") as handle_merged:
            for line in handle_merged:
                d = line.strip().split()
                merged_nodes[d[0]] = d[2]

    handle = open(options.dbdir+"gene_sp.list", "w")
    for sp in sp_gene:
        for geneid in sp_gene[sp]:
            if sp in desc_anc:
                handle.write("%s\t%s\n" % (geneid, sp))
            elif sp in merged_nodes:
                handle.write(geneid+"\t"+merged_nodes[sp]+"\n")
            elif sp in delnodes:
                handle.write(geneid+"\t-\n")
            elif sp in sp_conversion:
                handle.write(geneid+"\t"+sp_conversion[sp]+"\n")
            else:
                handle.write(geneid+"\t-\n")
    handle.close()
    del sp_gene


def write_pf_noannot(options, annot_gene_set, fullfamilyset):
    '''
    Store a list of Protein Families without genes with annotation.
    '''
    #list of untractable families
    print "# Checking for untractable families...\n"
    no_annot = set()
    for pf in fullfamilyset:
        handle = open(options.dbdir+"align/gene_list/"+pf.upper()+".gene", "r")
        gene_set = set()
        for line in handle:
            d1 = line.strip()
            d2 = d1[0:d1.find("/")]
            if d2 in annot_gene_set:
                gene_set.add(d1)
        if len(gene_set) == 0:
            no_annot.add(pf)
        handle.close()
    gene_set = set()
    handle = open(options.dbdir+"pf_noannot.list", "w")
    for pf in no_annot:
        handle.write(pf+"\n")
    handle.close()


def rm_no_annot_files(options):
    '''
    Remove families without annotations from the stored database.
    '''
    noannot_list = set()

    with open(options.dbdir+"pf_noannot.list") as handle:
        for line in handle:
            noannot_list.add(line.strip().split()[0])
    for item in noannot_list:
        os.remove(options.dbdir+"align/"+item+".fasta")
        os.remove(options.dbdir+"align/gene_list/"+item+".gene")


def is_number(s):
    '''
    Check if input is a valid number.
    '''
    try:
        int(s)
        return True
    except ValueError:
        return False


def summ_tax1(options):
    '''
    Prepare Taxonomy information from "taxonomy.txt".
    '''
    if not os.path.exists(options.dbdir+"summary_taxonomy.txt") or options.force:
        print "# Preparing Taxonomy information - \"taxonomy.txt\"...\n"
        taxonomy = open(options.dbdir+"taxonomy.txt", 'r')
        outf = open(options.dbdir+"summary_taxonomy.txt", 'w')
        while True:
            d = taxonomy.readline().split("\t")
            if len(d) < 2:
                break
            elif d[0] != "'0'":
                outf.write(d[0].replace(" ", "_").replace("\'", "")+"\t")
                if not is_number(d[5]):
                    outf.write("-\t")
                else:
                    outf.write(d[5].replace(" ", "_").replace("\'", "")+"\t")
                outf.write(d[2].replace(" ", "_").replace("\'", ""))
                outf.write(d[1].replace(" ", "_").replace("\'", "")+"\n")
        taxonomy.close()
        outf.close()


def summ_tax2(options):
    '''
    Prepare Taxonomy information from "ncbi_taxonomy.obo".
    '''
    if not os.path.exists(options.dbdir+"summary_ncbi_taxonomy.obo") or options.force:
        print "# Preparing Taxonomy information - \"ncbi_taxonomy.obo\"...\n"
        ncbi_taxonomy = open(options.dbdir+"ncbi_taxonomy.obo", 'r')
        outf = open(options.dbdir+"summary_ncbi_taxonomy.obo", 'w')
        x = {}
        x["ID"] = ""
        x["IS_A"] = ""
        for line in ncbi_taxonomy:
            d = line.replace(": NCBITaxon:", "\t").split()
            if len(d) > 1:
                if d[0] == "id":
                    x["ID"] = d[1]
                elif d[0] == "is_a":
                    x["IS_A"] = d[1]
                    if is_number(x["ID"]) and is_number(x["IS_A"]):
                        outf.write(x["ID"]+"\t"+x["IS_A"])
                        outf.write("\n")
                    x = {}
                    x["ID"] = ""
                    x["IS_A"] = ""
            else:
                continue
        outf.close()
        ncbi_taxonomy.close()


def obo2_functionontology(options):
    '''
    Convert "gene_ontology.1_2.obo" to function flat file format.
    '''

    from obo2flat import import_obogo
    from obo2flat import create_child_nodes
    from obo2flat import deep_search
    from obo2flat import Stack
    from obo2flat import go_node

    options.input = os.path.abspath(options.dbdir+"gene_ontology.1_2.obo")
    options.output = os.path.abspath(options.dbdir+"/function.ontology")

    if options.force and os.path.exists(options.output):
        os.remove(options.output)

    print "# Building function.ontology...\n"
    #import
    go_tree = import_obogo(options.input)
    go_tree = create_child_nodes(go_tree)

    #export
    outfile = open(options.output, "w")
    outfile.write("!autogenerated-by:     Sifter-T's OBO2Flat script\n")
    if "saved-by" in go_tree:
        outfile.write("!saved-by:             "+go_tree["saved-by"]+"\n")
    else:
        outfile.write("!saved-by:             -\n")
    outfile.write("!date:                 "+go_tree["date"]+"\n")
    outfile.write("!version: "+go_tree["version"]+"\n")
    outfile.write("!type: % OBO_REL:is_a is_a\n")
    outfile.write("$Gene_Ontology ; GO:0003673\n")

    #recursion
    stack = Stack()
    stack.push("GO:0003674")
    last_h = 0

    deep_search(stack, last_h, go_tree, outfile)

#    os.system("python obo2flat.py -f -i "+os.path.abspath(options.dbdir).replace(" ","\ ")+
#              "gene_ontology.1_2.obo -o "+os.path.abspath(options.dbdir).replace(" ","\ ")+
#              "/function.ontology")


def write_go_names(options):
    '''
    Gene ontology "id -> term" conversion file.
    '''
    handle = open(options.dbdir+"gene_ontology.1_2.obo", "r")
    handle_out = open(options.dbdir+"go_names.txt", "w")
    go_names = dict()
    go_id = str()
    name = str()
    namespace = str()
    is_obsolete = False
    for line in handle:
        if "[Term]" in line or "[Typedef]" in line:
            if not is_obsolete and "molecular" in namespace and "function" in namespace:
                go_names[go_id] = name
            go_id = str()
            name = str()
            namespace = str()
            is_obsolete = False
        elif line[0:4] == "id: ":
            go_id = line[line.index(" ")+1:len(line)-1]
        elif "name: " in line:
            name = line[line.index(" ")+1:len(line)-1]
        elif "namespace: " in line:
            namespace = line[line.index(" ")+1:len(line)-1]
        elif "is_obsolete: " in line:
            is_obsolete = True
    if not is_obsolete and "molecular" in namespace and "function" in namespace:
        go_names[go_id] = name
    handle.close()
    for go_id in go_names:
        handle_out.write(str(int(go_id[3:]))+"\t"+go_names[go_id]+"\n")
    handle_out.close()
    del is_obsolete
    del namespace
    del go_names
    del go_id
    del name


def prepare_hmm(options):
    '''
    Use hmmer3 to prepare Pfam-A hmm alignments.
    '''
    print "# Preparing PFam's HMM alignments...\n"
    os.system("hmmpress "+os.path.abspath(options.dbdir).replace(" ", "\ ")+"/Pfam-A.hmm")


def compDecomp(compObj, srcName, dstName):
    source = file(srcName, "r")
    dest = file(dstName, "w")
    block = source.read(2048)
    while block:
        cBlock = compObj.compress(block)
        dest.write(cBlock)
        block = source.read(2048)
    cBlock = compObj.flush()
    dest.write(cBlock)
    source.close()
    dest.close()


def compress_source_files(dbdir, unnecessary_files):
    '''
    Compress sorce database files - worker.
    '''
    while True:
        try:
            dbfile = q.get(block=True, timeout=0.1)
        except Empty:
            break
        else:
            if os.path.exists(dbdir+dbfile) and not os.path.exists(dbdir+"source/"+dbfile+".gz"):
                print "# Compressing file \""+dbfile+"\"..."
                compObj1 = zlib.compressobj()
                compDecomp(compObj1, dbdir+dbfile, dbdir+"source/"+dbfile+".gz")
                if dbfile in unnecessary_files and os.path.exists(dbdir+dbfile) and os.path.exists(dbdir+"source/"+dbfile+".gz"):
                    os.remove(dbdir+dbfile)
            q.task_done()


def compress_source_files_multi(options):
    '''
    Compress sorce database files - multithreading
    '''
    print "# Compressing post-processing unneeded database source files...\n"

    compressible_files = ["gene_association.goa_uniprot", "taxonomy.txt",
                          "ncbi_taxonomy.obo", "delnodes.dmp", "merged.dmp",
                          "uniprot_sprot.dat", "uniprot_trembl.dat",
                          "gene_ontology.1_2.obo", "Pfam-A.hmm.dat",
                          "Pfam-A.hmm", "Pfam-A.full"]

    unnecessary_files = ["gene_association.goa_uniprot", "taxonomy.txt",
                          "ncbi_taxonomy.obo", "delnodes.dmp", "merged.dmp",
                          "uniprot_sprot.dat", "uniprot_trembl.dat",
                          "gene_ontology.1_2.obo", "Pfam-A.full"]

    dbdir = options.dbdir
    file_size = list()
    for dbfile in compressible_files:
        if os.path.exists(dbdir+dbfile):
            file_size.append((os.path.getsize(dbdir+dbfile), dbfile))
    file_size.sort()
    file_size.reverse()
    for item in file_size:
        if os.path.exists(dbdir+item[1]) \
        and not os.path.exists(dbdir+"source/"+item[1]+".gz"):
            q.put(item[1])

    num_proc = int(max(round(multiprocessing.cpu_count()/4.0),1))

    if not os.path.exists(dbdir+"source/"):
        os.mkdir(dbdir+"source/")

    for i in range(num_proc):
        proc = Process(target=compress_source_files, name='%i' % (i+1), args=(dbdir, unnecessary_files))
        proc.start()

    sleep(num_proc*0.2)
    q.join()
    sleep(num_proc*0.05)
    if proc.is_alive() and q.empty():
        sleep(num_proc*0.2)
        if proc.is_alive() and q.empty():
            proc.terminate()


def _main():
    '''
    Main function for standalone usage.
    '''
    #Defines usage and help
    usage = "\n     %prog -d DIR"
    description = "Database preparation for SIFTER pipeline usage."
    #Defines input variables
    parser = OptionParser(usage=usage, version="%prog 0.8", description=description)
    parser.add_option("-d", "--directory",
        dest="dbdir",
        help="Full path for the database directory. " \
             "The needed files are: "                                    \
              "gene_association.goa_uniprot, "                           \
              "taxonomy.txt, "                                           \
              "ncbi_taxonomy.obo, "                                      \
              "delnodes.dmp, "                                           \
              "merged.dmp, "                                             \
              "uniprot_sprot.dat, "                                      \
              "uniprot_trembl.dat, "                                     \
              "gene_ontology.1_2.obo, "                                  \
              "Pfam-A.full, "                                            \
              "Pfam-A.hmm, "                                             \
              "Pfam-A.hmm.dat",
        metavar="/DIRECTORY/DIR/")
    parser.add_option("-f", "--force",
        help="(Optional) In case of new databases, force database preparation."\
             " (default False)",
        action="store_true",
        dest="force",
        default=False)
    parser.add_option("-c", "--compress",
        help="(Optional) Compress database source files."\
             " (default False)",
        action="store_true",
        dest="compress",
        default=False)
    parser.add_option("-t", "--threads",
        dest="threads",
        default=int(round(multiprocessing.cpu_count()/2.0)),
        help="(Optional) Number of parallel CPU workers to use for " \
             "multithreads (default all).",
        metavar="NUM",
        type="int")
    (options, args) = parser.parse_args()

    if len(args) > 0:
        print "\n# Extra arguments. Wrong usage. Exiting... \n"
        sys.exit(1)

    print "\n###  --------------------------------------------------  ###"
    print "###                  Database Preparation.               ###"
    print "###  --------------------------------------------------  ###\n"
    print "# Checking parameters and files...\n"

    ### Exception care ###
    if not options.dbdir:
        print "Not all required parameters were specified. Type \"-h\" for help. \nExiting..."
        sys.exit(1)

    if options.dbdir:
        if not os.path.isdir(options.dbdir):
            print options.dbdir+" is not an existing directory."
            sys.exit()

    if not os.access(options.dbdir, os.R_OK) or not os.access(options.dbdir, os.W_OK):
        print "The actual user does not have READING and/or WRITE access to " \
              "the specified directory. \nExiting..."
        sys.exit(1)

    options.dbdir = os.path.abspath(options.dbdir)
    if options.dbdir[-1] != "/":
        options.dbdir = options.dbdir+"/"

    required_files = ["gene_association.goa_uniprot", "taxonomy.txt",
                      "ncbi_taxonomy.obo", "delnodes.dmp", "merged.dmp",
                      "uniprot_sprot.dat", "uniprot_trembl.dat",
                      "gene_ontology.1_2.obo", "Pfam-A.hmm.dat",
                      "Pfam-A.hmm", "Pfam-A.full"]

    for item in required_files:
        if not os.path.exists(options.dbdir+item):
            print "Not all source files are on %s folder. Please, confirm if"  \
                  " the following files are there:\n\t"                        \
                  "gene_association.goa_uniprot\n\t"                           \
                  "taxonomy.txt\n\t"                                           \
                  "ncbi_taxonomy.obo\n\t"                                      \
                  "delnodes.dmp\n\t"                                           \
                  "merged.dmp\n\t"                                             \
                  "uniprot_sprot.dat\n\t"                                      \
                  "uniprot_trembl.dat\n\t"                                     \
                  "gene_ontology.1_2.obo\n\t"                                  \
                  "Pfam-A.full\n\t"                                            \
                  "Pfam-A.hmm\n\t"                                             \
                  "Pfam-A.hmm.dat \n\n"                                        \
                  "Detailed information in \"README.txt\" file. "              \
                  "\nExiting..." % options.dbdir
            sys.exit(1)

    for item in required_files:
        if not os.access(options.dbdir+item, os.R_OK):
            print "The actual user does not have READING access to all "       \
                  "required files. Type \"-h\" for help. \nExiting..."
            sys.exit()


    remove_former_db(options)
    fullfamilyset = recover_stockholm(options)
    write_pfamlist(options, fullfamilyset)
    stockfasta_multi(options, fullfamilyset)
    get_acession_numbers(options)
    fullgeneset = set()
    ac_gene = {}
    summary_goa(options, fullgeneset, ac_gene)
    del ac_gene
    del fullgeneset
    annot_gene_set = get_annot_gene_set(options)
    write_annot_gene_list(options, annot_gene_set)
    write_pf_noannot(options, annot_gene_set, fullfamilyset)
    rm_no_annot_files(options)
    del fullfamilyset
    summ_tax1(options)
    summ_tax2(options)
    write_gene_sp(options, annot_gene_set)
    del annot_gene_set
    obo2_functionontology(options)
    write_go_names(options)
    prepare_hmm(options)

    if os.path.exists(options.dbdir+"pfam_gene_ac.list"):
        os.remove(options.dbdir+"pfam_gene_ac.list")
    if os.path.exists(options.dbdir+"gene.list"):
        os.remove(options.dbdir+"gene.list")
    if os.path.exists(options.dbdir+"annot_gene.list"):
        os.remove(options.dbdir+"annot_gene.list")

    if options.compress:
        compress_source_files_multi(options)

    print "# Cleaning up and exiting..."
    sys.exit()

if __name__ == '__main__':
    _main()

