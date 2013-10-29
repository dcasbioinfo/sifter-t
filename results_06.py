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
#     Almeida-e-Silva D.C. and Vêncio R.Z.N. Sifter-T: A functional         #
#     framework for large-scale probabilistic protein domain annotation.    #
#     (In preparation...)                                                   #
#                                                                           #
########## ########## ########### ########## ########## ########## ##########

"""
 * Summarizes results from Sifter-T's workflow.

    get_querynames_conversion(options)
    get_go_term_conversion(options)
    get_query_pfam__query_pfam_len(options)
    get_no_pfam(options, query_pfam, querynames_conversion)
    get_noannot_goa(options)
    get_noannot_ec(options)
    get_noannot_all(options, noannot_goa, noannot_ec, query_pfam)
    get_useful_pfam(options)
    get_useful_pfam_func(options, useful_pfam)
    get_useful_pfam_query_prob(options, useful_pfam)
    get_query_len_pfam_prob(options, useful_pfam_query_prob)
    get_memsize()
    get_disk_used(options)
    get_input_type(options)
    write_report_alt_txt(options, useful_pfam, query_pfam, query_pfam_len,
                         querynames_conversion, useful_pfam_query_prob,
                         useful_pfam_func, go_term_conversion)
    write_report_txt(options, useful_pfam, querynames_conversion, 
                    useful_pfam_query_prob, useful_pfam_func, 
                    go_term_conversion, input_type, memsize, disk_used,
                    nopfam, noannot_goa, noannot_ec, noannot_all,
                    query_len_pfam_prob)
    create_report_tab(options, useful_pfam, querynames_conversion, 
                     useful_pfam_query_prob, useful_pfam_func, 
                     go_term_conversion, input_type, memsize, disk_used,
                     nopfam, noannot_goa, noannot_ec, noannot_all, 
                     query_len_pfam_prob)
    _main()

"""

from optparse import OptionParser
import pickle
import time
import os
import sys


def get_querynames_conversion(options):
    '''
    Load gene names conversion. name[query] = original_name
    '''
    if options.type == "aa" or options.type == "nt":
        querynames_conversion = dict()
        handle = open(options.outdir+"input_names.txt","r")
        for line in handle:
            d = line.strip().split("\t")
            querynames_conversion[d[0]] = d[1]
        handle.close()
        return querynames_conversion


def get_go_term_conversion(options):
    '''
    Load GO_id -> "full term name" conversion table.
    '''
    go_term_conversion = dict()
    handle = open(options.dbdir+"go_names.txt","r")
    for line in handle:
        d = line.strip().split("\t")
        go_term_conversion[d[0]] = d[1]
    handle.close()
    return go_term_conversion


def get_query_pfam__query_pfam_len(options):
    '''
    Loads complete set of query names and their aminoacid sequence lenght
    '''
    if options.type == "aa" or options.type == "nt":
        query_pfam = dict()
        query_pfam_len = dict()
    
        handle = open(options.outdir+"query.pfam","r")
        for line in handle:
            d = line.strip().split()
            try:
                query_pfam[d[3]].add(d[0])
                query_pfam_len[d[3]].add(d[0]+"|"+d[1]+"-"+d[2])
            except:
                query_pfam[d[3]] = set()
                query_pfam[d[3]].add(d[0])
                query_pfam_len[d[3]] = set()
                query_pfam_len[d[3]].add(d[0]+"|"+d[1]+"-"+d[2])
        handle.close()
        return query_pfam, query_pfam_len


def get_no_pfam(options, query_pfam, querynames_conversion):
    '''
    Return set with input genes without associated Protein Family
    '''
    if options.type == "aa" or options.type == "nt":
        query_with_pfam = set()
        for fam in query_pfam:
            for item in query_pfam[fam]:
                if "_" in item:
                    query_with_pfam.add(item[0:item.find("_")])
                else:
                    query_with_pfam.add(item)
    
        nopfam = set(querynames_conversion)
        nopfam = nopfam - query_with_pfam
        return nopfam


def get_noannot_goa(options):
    '''
    Loads domains in genes associated to Protein Families without a single
    annotation to propagate.
    '''
    noannot_goa = set()
    handle = open(options.outdir+"without_annotations_goa.txt","r")
    for line in handle:
        d = line.strip().split()
        noannot_goa.add(d[0])
    
    handle.close()
    return noannot_goa


def get_noannot_ec(options):
    '''
    Loads domains in genes associated to Protein Families without annotations
    using the selected evidence codes, but have annotations on the remaining
    evidence codes.
    '''
    noannot_ec = set()
    handle = open(options.outdir+"without_annotations_ec.txt","r")
    for line in handle:
        d = line.strip().split()
        noannot_ec.add(d[0])
    handle.close()
    return noannot_ec


def get_noannot_all(options, noannot_goa, noannot_ec, query_pfam):
    '''
    Return a set of genes with domains that can not be treated by Sifter-T.
    '''
    if options.type == "aa" or options.type == "nt":
        noannot_all = dict()
        for fam in (noannot_goa | noannot_ec):
            noannot_all[fam] = set()
            for item in query_pfam[fam]:
                if "_" in item:
                    noannot_all[fam].add(item[0:item.find("_")])
                else:
                    noannot_all[fam].add(item)
        return noannot_all


def get_useful_pfam(options):
    '''
    Load a set of Protein Families that can be treated.
    '''
    useful_pfam = set()
    handle = open(options.outdir+"useful_pfam.txt","r")
    for line in handle:
        d = line.strip().split()
        useful_pfam.add(d[0])
    handle.close()
    return useful_pfam


def get_useful_pfam_func(options, useful_pfam):
    '''
    Loads the list of functions (according to GO) in each Protein Family
    '''
    useful_pfam_func = dict()
    for fam in useful_pfam:
        handle = open(options.outdir+fam.upper()+"/infer-"+fam.lower()+".fx","r")
        temp = handle.readline()
        temp = handle.readline()
        useful_pfam_func[fam] = temp.strip().split()
        handle.close()
    return useful_pfam_func


def get_useful_pfam_query_prob(options, useful_pfam):
    '''
    Loads the list of probabilities for all query sequences in all Protein 
    Families.
    '''
    useful_pfam_query_prob = dict()
    for fam in useful_pfam:
        useful_pfam_query_prob[fam] = dict()
        handle = open(options.outdir+fam.upper()+"/"+fam.upper()+".rdata","r")
        for line in handle:
            d = line.strip().split()
            if d[0][0:5] == 'query':
                useful_pfam_query_prob[fam][d[0]] = d[1:-1]
            elif d[0][0:4] != 'Node' and options.type == 'pf':
                useful_pfam_query_prob[fam][d[0]] = d[1:-1]
    return useful_pfam_query_prob


def get_query_len_pfam_prob(options, useful_pfam_query_prob):
    '''
    Loads Sifter2.0 annotations and return for each domain in each gene:
        Reading frame
        Start position - aminoacid
        Stop position - aminoacid
        Protein Family code
        List of probabilities
    '''
    if options.type == "aa" or options.type == "nt":
        query_len_pfam_prob = dict()
        for fam in useful_pfam_query_prob:
            for query in useful_pfam_query_prob[fam]:
                d = query.split("|")
                if len(d) == 2:
                    if "_" in d[0]:
                        name = d[0][0:d[0].find("_")]
                        strand = d[0][d[0].find("_")+1:]
                    else:
                        name = d[0]
                        strand = " "
                    start = int(d[1][0:d[1].find("-")])
                    stop = int(d[1][d[1].find("-")+1:])
                    pfam = fam
                try:
                    query_len_pfam_prob[name].add((strand, start, stop, 
                        pfam, tuple(useful_pfam_query_prob[fam][query])))
                except:
                    query_len_pfam_prob[name] = set()
                    query_len_pfam_prob[name].add((strand, start, stop, 
                        pfam, tuple(useful_pfam_query_prob[fam][query])))
        return query_len_pfam_prob


def get_memsize():
    '''
    Return system memory size.
    '''
    handle = open("/proc/meminfo", "r")
    memsize = int(handle.readlines()[0].strip().split()[1])
    handle.close()
    return memsize


def get_disk_used(options):
    '''
    Return the disk space used for the Output directory.
    '''
    fileset = set()
    for item in os.walk(options.outdir):
        for file_item in item[2]:
            fileset.add(item[0]+"/"+file_item)

    disk_used = sum([os.path.getsize(f) for f in fileset if os.path.isfile(f)])
    return disk_used


def get_input_type(options):
    '''
    Return the full name for each input type code.
    '''
    if options.type == "aa":
        return "Aminoacid"
    elif options.type == "nt":
        return "Nucleotide"
    elif options.type == "pf":
        return "Protein Family"


def write_report_alt_txt(options, useful_pfam, query_pfam, query_pfam_len,
                         querynames_conversion, useful_pfam_query_prob,
                         useful_pfam_func, go_term_conversion):
    '''
    Create an alternative report file with 2 blocks:
        1 - Header
        2 - Annotated genes

    The annotated genes shows the follloing structure for aminoacid input:
    [protein family code]
        [gene name]|[start position]-[stop position] 
    		PROB_MAX	GO	GO_DESCRIPTION
    		PROB_MID	GO	GO_DESCRIPTION
    		PROB_MIN	GO	GO_DESCRIPTION

    And this structure for nucleotide input:
    [protein family code]
        [gene name]|[reading frame]|[start position]-[stop position] 
    		PROB_MAX	GO	GO_DESCRIPTION
    		PROB_MID	GO	GO_DESCRIPTION
    		PROB_MIN	GO	GO_DESCRIPTION

    Protein Family input does not need an alternative report.
    '''
    if options.type == "aa" or options.type == "nt":
        handle = open(options.outdir+"REPORT_alt.txt","w")
        handle.write("----- The following genes/protein families were " \
                     "annotated using Sifter-T. -----\n")
        for fam in useful_pfam:
            handle.write(fam+"\n")
            for query in query_pfam[fam]:
                for query2 in query_pfam_len[fam]:
                    d = query2.split("|")
                    if len(d) == 2 and d[0] == query:
                        if "_" in query:
                            name = query[0:query.find("_")]
                            frame = query[query.find("_")+1:]+"|"
                        else:
                            name = query
                            frame = ""
                        handle.write("\t"+querynames_conversion[name]+"" \
                                     "|"+frame+d[1]+"\n")
                        temp_list = list()
                        for i in range(len(useful_pfam_query_prob[fam][query2])):
                            if float(useful_pfam_query_prob[fam][query2][i]) >= options.scut:
                                temp_list.append((float(useful_pfam_query_prob[fam][query2][i]),
                                                  useful_pfam_func[fam][i]))
                        temp_list.sort()
                        temp_list.reverse()
                        for item in temp_list:
                            if float(item[0]) > float(1):
                                handle.write("\t\t"+str("0.999999999999")+"\t")
                            else:
                                handle.write("\t\t"+str(item[0])+"\t")
                            handle.write(item[1]+"\t"+go_term_conversion[item[1]]+"\n")
        handle.close()


def write_report_txt(options, useful_pfam, querynames_conversion, 
                     useful_pfam_query_prob, useful_pfam_func, 
                     go_term_conversion, input_type, memsize, disk_used,
                     nopfam, noannot_goa, noannot_ec, noannot_all,
                     query_len_pfam_prob):
    '''
    Create a report file with 4 blocks:
        1 - Header
        2 - Genes in families without a single annotation - Sifter can not
            annotate those genes
        3 - Families with annotations, but none among the evidence codes 
            selected. Sifter can annotate if you include other evidence codes.
        4 - Annotated genes
  
    The annotated genes shows the follloing structure for aminoacid input:
    [gene name]
    	[start position]-[stop position] - [protein family code]
    		PROB_MAX	GO	GO_DESCRIPTION
    		PROB_MID	GO	GO_DESCRIPTION
    		PROB_MIN	GO	GO_DESCRIPTION

    And this structure for nucleotide input:
    [gene name]
    	[reading frame] - [start pos]-[stop pos] - [protein family code]
    		PROB_MAX	GO	GO_DESCRIPTION
    		PROB_MID	GO	GO_DESCRIPTION
    		PROB_MIN	GO	GO_DESCRIPTION

    The report for Protein Families as inputs uses the same basic structure as 
    the REPORT_alt.txt for other types of input:
    [protein family code]
        [gene name]
    		PROB_MAX	GO	GO_DESCRIPTION
    		PROB_MID	GO	GO_DESCRIPTION
    		PROB_MIN	GO	GO_DESCRIPTION

    '''
    handle = open(options.outdir+"REPORT.txt","w")
    handle.write("# ----- ----- ----- ----- ----- ----- ----- ----- ----- -----\n")
    handle.write("# Sifter-T Workflow (2013)\n")
    handle.write("# Author: Danillo C. Almeida-e-Silva (dcas.bioinfo@gmail.com)\n")
    handle.write("# Date: "+time.asctime( time.localtime(time.time()) )+"\n") 
    handle.write("# Parameters: \n") 
    handle.write("# \tType: \t\t\t\t"+str(input_type)+"\n") 
    handle.write("# \tInput file: \t\t\t"+str(os.path.abspath(options.file))+"\n") 
    handle.write("# \tDatabase dir: \t\t\t"+str(os.path.abspath(options.dbdir))+"\n") 
    handle.write("# \tOutput dir: \t\t\t"+str(os.path.abspath(options.outdir))+"\n") 
    handle.write("# \tSifter2.0 dir: \t\t\t"+str(os.path.abspath(options.sdir))+"\n") 
    handle.write("# \tSifter-T dir: \t\t\t"+str(os.path.abspath(os.getcwd()))+"\n") 
    handle.write("# \tEvidence Code(s) selected: \t"+str(options.experimental)+"\n") 
    handle.write("# \tSpecie's annotations removed: \t"+str(options.species)+"\n") 
    handle.write("# \tSpecies tree branche's annotations removed: \t"+str(options.branch)+"\n") 
    if options.type == "nt":
        handle.write("# \tTranslation Table: \t\t"+str(options.translation)+"\n") 
    if options.reconciliation:
        handle.write("# \tReconciliation: \t\tYes\n")
    else:
        handle.write("# \tReconciliation: \t\tNo\n") 
    if options.input_species != "0":
        handle.write("# \tInput Species: \t\t\t"+options.input_species+"\n")
    handle.write("# \tSifter2.0 Cutoff: \t\t"+str(options.scut)+"\n") 
    if options.type != "pf":
        handle.write("# \tPFamScan's e-value Cutoff: \t"+str(options.pcut)+"\n") 
    handle.write("# \tThreads: \t\t\t"+str(options.threads)+"\n") 
    handle.write("# \tMemory: \t\t\t"+str(memsize/1024)+" MB\n") 
    handle.write("# \tDisk Used (Output dir): \t"+str(disk_used/1024)+" MB\n") 
    handle.write("# ----- ----- ----- ----- ----- ----- ----- ----- ----- -----\n")
    if options.type == "pf":

        handle.write("\n")
        for fam in useful_pfam:
            handle.write(fam+"\n")
            for query in useful_pfam_query_prob[fam]:
                name = query
                handle.write("\t"+name+"\n")
                temp_list = list()
                for i in range(len(useful_pfam_query_prob[fam][query])):
                    if float(useful_pfam_query_prob[fam][query][i]) >= options.scut:
                        temp_list.append((float(useful_pfam_query_prob[fam][query][i]), useful_pfam_func[fam][i]))
                temp_list.sort()
                temp_list.reverse()
                for item in temp_list:
                    if float(item[0]) > float(1):
                        handle.write("\t\t"+str("0.999999999999")+"\t")
                    else:
                        handle.write("\t\t"+str(item[0])+"\t")
                    handle.write(item[1]+"\t"+go_term_conversion[item[1]]+"\n")
            handle.write("\n")
        handle.close()
    
    elif options.type == "aa" or options.type == "nt":    
        if len(nopfam) == 0:
            handle.write("\n# ----- There is no gene without associated PFam." \
                         " -----\n\n")
        else:
            handle.write("\n# ----- "+str(len(nopfam))+" genes do not have" \
                         " associated PFam. -----\n")
            n = 1
            for query in nopfam:
                handle.write(str(n)+"\t"+querynames_conversion[query]+"\n")
                n = n+1
        handle.write("\n")
    
        if len(noannot_goa) == 0:
            handle.write("\n# ----- There is no gene/protein family " \
                "non-annotated due to lack of GOA annotations. -----\n\n")
        else:
            # num of genes
            n = 0
            for fam in noannot_goa:
                n = n + len(noannot_all[fam])
            handle.write("\n# ----- "+str(len(noannot_goa))+" identified " \
                "protein families do not have genes with GOA annotations. " \
                "Certain domains in "+str(n)+" genes cannot be treated. -----\n")
            for fam in noannot_goa:
                handle.write("\n"+fam+"\n")
                for query in noannot_all[fam]:
                    handle.write("\t"+querynames_conversion[query]+"\n")
    
    
        if len(noannot_ec) == 0:
            handle.write("\n# ----- There is no gene/protein family " \
                "non-annotated due to incomplete evidence codes " \
                "selection. -----\n\n")
        else:
            handle.write("\n# ----- "+str(len(noannot_ec))+" genes were not" \
                " annotated due to incomplete evidence code selection. -----\n")
            for fam in noannot_ec:
                handle.write("\n"+fam+"\n")
                for query in noannot_all[fam]:
                    handle.write("\t"+querynames_conversion[query]+"\n")
    
        handle.write("\n# ----- The following genes/protein families were" \
                     " annotated using Sifter-T. -----\n")
        query_list = list()
        for item in query_len_pfam_prob:
            query_list.append(int(item[5:]))
        query_list.sort()
        query_list2 = list()
        for item in query_list:
            query_list2.append("query"+str(item))
    
        for query in query_list2:
            handle.write(querynames_conversion[query]+"\n")
            order_list = list(query_len_pfam_prob[query])
            order_list.sort()
            for item in order_list:
                if item[0] == " ":
                    handle.write("\t"+str(item[1])+"-"+str(item[2])+" - "+item[3]+"\n")
                else:
                    handle.write("\t"+item[0]+" - "+str(item[1])+"-"+str(item[2])+" - "+item[3]+"\n")
                temp_list = list()
                fam = item[3]
                for i in range(len(item[4])):
                    temp_list.append((float(item[4][i]), useful_pfam_func[fam][i]))
                temp_list.sort()
                temp_list.reverse()
                for item2 in temp_list:
                    if item2[0] >= options.scut:
                        if float(item2[0]) > float(1):
                            handle.write("\t\t0.999999999999\t")
                        else:
                            handle.write("\t\t"+str(item2[0])+"\t")
                        handle.write(item2[1]+"\t"+go_term_conversion[item2[1]]+"\n")
            handle.write("\n")    
        handle.close()


def create_report_tab(options, useful_pfam, querynames_conversion, 
                     useful_pfam_query_prob, useful_pfam_func, 
                     go_term_conversion, input_type, memsize, disk_used, 
                     nopfam, noannot_goa, noannot_ec, noannot_all,
                     query_len_pfam_prob):
    '''
    Create a "Tab Separated Values" file with the follloing structure:

    # Collumns: <Gene Name> <reading frame> <aminoacid start position> 
    <aminoacid stop position> <PFam code> <probability> <GO Code> <GO Name>
    
    Gene Name - are returned in the same order as input.

    Reading Frame - Uses the Watson/Crick notation
        c - Crick - sense strand
        w - Watson- antisense strand
        1, 2, 3 - The reading frame

    Start/Stop position - For nucleotides or aminoacid inputs, it shows the 
    begining and the end of an aminoacid region (result of reading frame 
    translation or not) with identified Protein Family.
        
    PFam Code - The code for the Protein Family identified at this region.

    Probability - Sifter2.0 returns a list of probabilities for a list of
    functions. Sifter-T shows only the functions (according to GO) with
    probabilities above the cutoff. The probability can be read as "The 
    probability of the identified domain have the function [name of the
    function] given the Protein Family structure and annotations. Sometimes  
    you can find probabilities of "0.999999999999". It means that this 
    particular Protein Family have only one GO number annotated along the tree,
    so the probability of the identified domain have the function X, given that
    this is the only function on the tree, is almost 1.

    GO Code - The Gene Ontology code for the annotated function.
    
    GO Name - The Gene Ontology full name for the annotated GO code.

    '''
    handle = open(options.outdir+"REPORT.tab","w")
    handle.write("# ----- ----- ----- ----- ----- ----- ----- ----- ----- -----\n")
    handle.write("# Sifter-T Workflow (2013)\n")
    handle.write("# Author: Danillo C. Almeida-e-Silva (dcas.bioinfo@gmail.com)\n")
    handle.write("# Date: "+time.asctime( time.localtime(time.time()) )+"\n") 
    handle.write("# Parameters: \n") 
    handle.write("# \tType: \t\t\t\t"+str(input_type)+"\n") 
    handle.write("# \tInput file: \t\t\t"+str(os.path.abspath(options.file))+"\n") 
    handle.write("# \tDatabase dir: \t\t\t"+str(os.path.abspath(options.dbdir))+"\n") 
    handle.write("# \tOutput dir: \t\t\t"+str(os.path.abspath(options.outdir))+"\n") 
    handle.write("# \tSifter2.0 dir: \t\t\t"+str(os.path.abspath(options.sdir))+"\n") 
    handle.write("# \tSifter-T dir: \t\t\t"+str(os.path.abspath(os.getcwd()))+"\n") 
    handle.write("# \tEvidence Code(s) selected: \t"+str(options.experimental)+"\n") 
    handle.write("# \tSpecie's annotations removed: \t"+str(options.species)+"\n") 
    handle.write("# \tSpecies tree branche's annotations removed: \t"+str(options.branch)+"\n") 
    if options.type == "nt":
        handle.write("# \tTranslation Table: \t\t"+str(options.translation)+"\n") 
    if options.reconciliation:
        handle.write("# \tReconciliation: \t\tYes\n")
    else:
        handle.write("# \tReconciliation: \t\tNo\n") 
    if options.input_species != "0":
        handle.write("# \tInput Species: \t\t\t"+options.input_species+"\n")
    handle.write("# \tSifter2.0 Cutoff: \t\t"+str(options.scut)+"\n") 
    if options.type != "pf":
        handle.write("# \tPFamScan's e-value Cutoff: \t"+str(options.pcut)+"\n") 
    handle.write("# \tThreads: \t\t\t"+str(options.threads)+"\n") 
    handle.write("# \tMemory: \t\t\t"+str(memsize/1024)+" MB\n") 
    handle.write("# \tDisk Used (Output dir): \t"+str(disk_used/1024)+" MB\n") 
    handle.write("# ----- ----- ----- ----- ----- ----- ----- ----- ----- -----\n")

    if options.type == "aa" or options.type == "nt":
        handle.write("# <Gene Name> <reading frame> <aminoacid start position>"\
                      " <aminoacid stop position> <PFam code> <probability>"   \
                      " <GO Code> <GO Name>\n\n")
        for query in nopfam:
            handle.write(querynames_conversion[query].replace("\t"," ")+"\t-\t-\t-\t-\t-\t-\t-\n")
        for fam in noannot_goa:
            for query in noannot_all[fam]:
                handle.write(querynames_conversion[query].replace("\t"," ")+"\t-\t-\t-\t"+fam+"\t-\t-\t-\n")
        for fam in noannot_ec:
            for query in noannot_all[fam]:
                handle.write(querynames_conversion[query].replace("\t"," ")+"\t-\t-\t-\t"+fam+"\t-\t-\t-\n")
        query_list = list()
        for item in query_len_pfam_prob:
            query_list.append(int(item[5:]))
        query_list.sort()
        query_list2 = list()
        for item in query_list:
            query_list2.append("query"+str(item))

        for query in query_list2:
            order_list = list(query_len_pfam_prob[query])
            order_list.sort()
            for item in order_list:
                temp_list = list()
                fam = item[3]
                for i in range(len(item[4])):
                    temp_list.append((float(item[4][i]), useful_pfam_func[fam][i]))
                temp_list.sort()
                temp_list.reverse()
                for item2 in temp_list:
                    if item2[0] >= options.scut:
                        if float(item2[0]) > float(1):
                            handle.write(querynames_conversion[query]+"\t" \
                                ""+item[0].replace(" ","-")+"\t" \
                                ""+str(item[1])+"\t"+str(item[2])+"\t" \
                                ""+item[3]+"\t"+"0.999999999999\t"+item2[1]+"" \
                                "\t"+go_term_conversion[item2[1]]+"\n")
                        else:
                            handle.write(querynames_conversion[query]+"\t" \
                                ""+item[0].replace(" ","-")+"\t"
                                ""+str(item[1])+"\t"+str(item[2])+"\t" \
                                ""+item[3]+"\t"+str(item2[0])+"\t" \
                                ""+item2[1]+"\t" \
                                ""+go_term_conversion[item2[1]]+"\n")
        handle.close()

    elif options.type == "pf":
        handle.write("# <Gene Name> <PFam> <probability> <GO Code> <GO Name>\n\n")
        for fam in useful_pfam:
            for query in useful_pfam_query_prob[fam]:
                name = query
                temp_list = list()
                for i in range(len(useful_pfam_query_prob[fam][query])):
                    if float(useful_pfam_query_prob[fam][query][i]) >= options.scut:
                        temp_list.append((float(useful_pfam_query_prob[fam][query][i]),
                                          useful_pfam_func[fam][i]))
                temp_list.sort()
                temp_list.reverse()
                for item in temp_list:
                    if float(item[0]) > float(1):
                        handle.write(name+"\t"+fam+"\t"+str("0.999999999999")+"\t")
                    else:
                        handle.write(name+"\t"+fam+"\t"+str(item[0])+"\t")
                    handle.write(item[1]+"\t"+go_term_conversion[item[1]]+"\n")
        handle.close()



def _main():
    ''' 
    Main function for standalone usage.
    '''
    #Defines usage and help
    usage = "\n 	%prog -d DIR"
    description = "Summarizes results from Sifter-T's workflow."
    #Defines input variables
    parser = OptionParser(usage=usage, version="%prog 0.1.0", description=description)
    parser.add_option("-d", "--directory", 
        dest="dir", 
        help="full path for the working directory.", 
        metavar="/DIRECTORY/DIR/")
    (options, args) = parser.parse_args()
    if len(args) > 0:
        print "\n# Extra arguments. Wrong usage. Exiting... \n"
        sys.exit(1)

    if not os.path.exists(options.dir+"options.pk"):
        print "\"options.pk\" not found.\nExiting..."
        sys.exit(1)

    options = pickle.load(file(options.dir+"options.pk", "r"))


    print "\n###  --------------------------------------------------  ###" 
    print "###               Generating REPORT files                ###" 
    print "###  --------------------------------------------------  ###\n" 
    
    querynames_conversion = get_querynames_conversion(options)
    go_term_conversion = get_go_term_conversion(options)
    if options.type != "pf":
        query_pfam, query_pfam_len = get_query_pfam__query_pfam_len(options)
        nopfam = get_no_pfam(options, query_pfam, querynames_conversion)
    else: 
        nopfam = set()
    noannot_goa = get_noannot_goa(options)
    noannot_ec = get_noannot_ec(options)
    if options.type != "pf":
        noannot_all = get_noannot_all(options, noannot_goa, noannot_ec, query_pfam)
    else:
        noannot_all = set()
    useful_pfam = get_useful_pfam(options)
    useful_pfam_func = get_useful_pfam_func(options, useful_pfam)
    useful_pfam_query_prob = get_useful_pfam_query_prob(options, useful_pfam)
    query_len_pfam_prob = get_query_len_pfam_prob(options, useful_pfam_query_prob)
    memsize = get_memsize()
    disk_used = get_disk_used(options)
    input_type = get_input_type(options)
    if options.type != "pf":
        write_report_alt_txt(options, useful_pfam, query_pfam, query_pfam_len,
                             querynames_conversion, useful_pfam_query_prob,
                             useful_pfam_func, go_term_conversion)
    write_report_txt(options, useful_pfam, querynames_conversion,
                     useful_pfam_query_prob, useful_pfam_func, go_term_conversion,
                     input_type, memsize, disk_used, nopfam, noannot_goa, noannot_ec,
                     noannot_all, query_len_pfam_prob)
    create_report_tab(options, useful_pfam, querynames_conversion,
                      useful_pfam_query_prob, useful_pfam_func, go_term_conversion,
                      input_type, memsize, disk_used, nopfam, noannot_goa, noannot_ec,
                      noannot_all, query_len_pfam_prob)

    print "# The result summary is available at: \n"
    print "    "+options.outdir+"REPORT.txt\n"
    if options.type != "pf":
        print "    "+options.outdir+"REPORT_alt.txt\n"
    print "    "+options.outdir+"REPORT.tab\n"

    sys.exit()

if __name__ == '__main__':
    _main()
    print "# The result summary is available at: \n"
    print "    "+options.outdir+"REPORT.txt\n"
    if options.type != "pf":
        print "    "+options.outdir+"REPORT_alt.txt\n"
    print "    "+options.outdir+"REPORT.tab\n"
