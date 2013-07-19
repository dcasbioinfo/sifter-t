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
 * 
"""

if __name__ == '__main__':
    from optparse import OptionParser
    import pickle
    import os
    import sys
    try: 
        from Bio import SeqIO
    except:
        raise Exception("Can't load BioPython. Exiting...")
    try: 
        import dendropy
    except:
        raise Exception("Can't load Dendropy. Exiting...")
    import multiprocessing
    from multiprocessing import JoinableQueue


def check_databases(options):
    '''
    Check databases.
    '''
    print "# Checking database...\n"
    if not os.path.isdir(options.dbdir):        
        print options.dbdir+" is not an existing directory."
        sys.exit()
    if not os.access(options.dbdir, os.R_OK):
        print "The actual user does not have READING access to %s directory."  \
              " \nExiting..." % options.dbdir
        sys.exit(1)
    options.dbdir = os.path.abspath(options.dbdir)
    if options.dbdir[-1] != "/":
        options.dbdir = options.dbdir+"/"
    if not os.path.exists(options.dbdir+"summary_gene_association.goa_uniprot") or \
       not os.path.exists(options.dbdir+"summary_ncbi_taxonomy.obo") or        \
       not os.path.exists(options.dbdir+"annot_gene.list") or                  \
       not os.path.exists(options.dbdir+"align") or                            \
       not os.path.exists(options.dbdir+"summary_taxonomy.txt") or             \
       not os.path.exists(options.dbdir+"Pfam-A.hmm.h3f") or                   \
       not os.path.exists(options.dbdir+"Pfam-A.hmm") or                       \
       not os.path.exists(options.dbdir+"Pfam-A.hmm.dat") or                   \
       not os.path.exists(options.dbdir+"Pfam-A.hmm.h3i") or                   \
       not os.path.exists(options.dbdir+"Pfam-A.hmm.h3m") or                   \
       not os.path.exists(options.dbdir+"Pfam-A.hmm.h3p") or                   \
       not os.path.exists(options.dbdir+"function.ontology") or                \
       not os.path.exists(options.dbdir+"gene_sp.list") or                     \
       not os.path.exists(options.dbdir+"go_names.txt") or                     \
       not os.path.exists(options.dbdir+"pfam.list") or                        \
       not os.path.exists(options.dbdir+"pfam_gene_ac.list") or                \
       not os.path.exists(options.dbdir+"function.ontology"):
        if not os.path.exists(options.dbdir+"gene_association.goa_uniprot") or \
           not os.path.exists(options.dbdir+"taxonomy.txt") or                 \
           not os.path.exists(options.dbdir+"ncbi_taxonomy.obo") or            \
           not os.path.exists(options.dbdir+"delnodes.dmp") or                 \
           not os.path.exists(options.dbdir+"merged.dmp") or                   \
           not os.path.exists(options.dbdir+"uniprot_sprot.dat") or            \
           not os.path.exists(options.dbdir+"uniprot_trembl.dat") or           \
           not os.path.exists(options.dbdir+"gene_ontology.1_2.obo") or        \
           not os.path.exists(options.dbdir+"Pfam-A.hmm.dat") or               \
           not os.path.exists(options.dbdir+"Pfam-A.hmm") or                   \
           not os.path.exists(options.dbdir+"Pfam-A.full"):
            print "Not all source files are on %s folder. Please, confirm if"  \
                  " the following files are there:\n\tgene_association."       \
                  "goa_uniprot\n\ttaxonomy.txt\n\tuniprot_sprot.dat\n\tuniprot"\
                  "_trembl.dat\n\tgene_ontology.1_2.obo\n\tPfam-A.hmm\n\t"     \
                  "Pfam-A.hmm.dat \n\nDetailed information in \"README.txt\""  \
                  " file. \nExiting..." % options.dbdir
            sys.exit(1)
        else:
            print "The basic database files are on the correct folder. However"\
                  " they were not prepared for Sifter-T use yet."
            prepare = raw_input("Prepare them now? (it can take several "      \
                                "minutes...) [y/n]: ")
            if prepare == "y":
                dbprep(options)
            elif prepare == "n":
                print "\nExiting..."
                sys.exit(1)
            elif prepare != "y" or prepare != "n":
                print "Not a valid option. \nExiting..."
                sys.exit(1)
    return options


def dbprep(options):
    '''
    Prepare databases for Sifter-T usage.
    '''
    from dbprep import remove_former_db
    from dbprep import recover_stockholm
    from dbprep import write_pfamlist
    from dbprep import get_memsize
    from dbprep import stockfastaconvert 
    from dbprep import stockfasta_multi
    from dbprep import get_acession_numbers
    from dbprep import get_fullgeneset
    from dbprep import write_gene_ac
    from dbprep import get_ac_gene
    from dbprep import summary_goa
    from dbprep import get_annot_gene_set
    from dbprep import write_annot_gene_list
    from dbprep import write_gene_sp
    from dbprep import write_pf_noannot
    from dbprep import is_number
    from dbprep import summ_tax1
    from dbprep import summ_tax2
    from dbprep import obo2_functionontology
    from dbprep import write_go_names
    from dbprep import prepare_hmm

    print "\n###  --------------------------------------------------  ###" 
    print "###                  Database Preparation.               ###" 
    print "###  --------------------------------------------------  ###\n" 
    print "# Checking parameters and files...\n"

    ### Exception care ###
    if not options.dbdir:
        print "Not all required parameters were specified. Type \"-h\" for "   \
              "help. \nExiting..."
        sys.exit(1)

    if options.dbdir:
        if not os.path.isdir(options.dbdir):
            print options.dbdir+" is not an existing directory."
            sys.exit()

    if not os.access(options.dbdir, os.R_OK) or \
       not os.access(options.dbdir, os.W_OK):
        print "The actual user does not have READING and/or WRITE access to " \
              "the specified directory. \nExiting..."
        sys.exit(1)

    if not os.path.exists(options.dbdir+"gene_association.goa_uniprot") or \
      not os.path.exists(options.dbdir+"Pfam-A.full") or \
      not os.path.exists(options.dbdir+"taxonomy.txt") or \
      not os.path.exists(options.dbdir+"uniprot_sprot.dat") or \
      not os.path.exists(options.dbdir+"uniprot_trembl.dat"):
        print "Not all required files are on specified folder. Type \"-h\" for"\
              " help. \nExiting..."
        sys.exit()

    if not os.access(options.dbdir+"gene_association.goa_uniprot", os.R_OK) or \
      not os.access(options.dbdir+"Pfam-A.full", os.R_OK) or \
      not os.access(options.dbdir+"taxonomy.txt", os.R_OK) or \
      not os.access(options.dbdir+"uniprot_sprot.dat", os.R_OK) or \
      not os.access(options.dbdir+"uniprot_trembl.dat", os.R_OK):
        print "The actual user does not have READING access to all required " \
              "files. Type \"-h\" for help. \nExiting..."
        sys.exit()

    remove_former_db(options)
    fullfamilyset = recover_stockholm(options)
    write_pfamlist(options, fullfamilyset)
    stockfasta_multi(options, fullfamilyset)
    gene_ac = get_acession_numbers(options)
    fullgeneset = get_fullgeneset(gene_ac)
    write_gene_ac(options, gene_ac)
    ac_gene = get_ac_gene(gene_ac)
    del(gene_ac)
    summary_goa(options, fullgeneset, ac_gene)
    del(ac_gene)
    del(fullgeneset)
    annot_gene_set = get_annot_gene_set(options)
    write_annot_gene_list(options, annot_gene_set)
    write_pf_noannot(options, annot_gene_set, fullfamilyset)
    del(fullfamilyset)
    write_gene_sp(options, annot_gene_set)
    del(annot_gene_set)
    summ_tax1(options)
    summ_tax2(options)
    obo2_functionontology(options)
    write_go_names(options)
    prepare_hmm(options)


def param_chk():
    '''
    Checks main parameters, files and creates "options.pk" for other scripts.
    '''
    from param_chk_00 import required_parameters
    from param_chk_00 import check_type
    from param_chk_00 import check_input_file0
    from param_chk_00 import is_number
    from param_chk_00 import check_input_file1
    from param_chk_00 import check_outdir
    from param_chk_00 import check_sifter_dir
    from param_chk_00 import check_sifter_files
    from param_chk_00 import check_pfam_scan
    from param_chk_00 import check_fasttree
    from param_chk_00 import check_notung
    from param_chk_00 import check_program
    from param_chk_00 import check_hmmer3
    from param_chk_00 import check_mafft
    from param_chk_00 import check_evcodes
    from param_chk_00 import check_translation_table
    from param_chk_00 import check_sifter_cutoff
    from param_chk_00 import check_pfam_cutoff
    from param_chk_00 import check_reconciliation
    from param_chk_00 import check_species
    from param_chk_00 import check_threads

#    from param_chk_00 import check_databases'''

    #Defines usage and help
    usage = "\n 	%prog -t XX -f FILE -o DIR -d DIR -s DIR " \
            "-i CODE [options]"
    description = "This script reads a nucleotide or protein multifasta FILE," \
                  " translates in 6 frames, identifies the PFAM ids associated"\
                  " to each frame, build a tree based on PFAM alignments and " \
                  "evidence codes, run SIFTER to propagate information from "  \
                  "families annotated to families of interest and returns a "  \
                  "list of annotations and his probabilities."
    exphelp = "Include Evidence Codes: Experimental Evidence Codes (EXP, IDA," \
              " IPI, IMP, IGI, IEP); Computational Analysis Evidence Codes "   \
              "(ISS, ISO, ISA, ISM, IGC, RCA, IBD, IBA), Author Statement "    \
              "Evidence Codes (TAS, NAS); Curator Statement Evidence Codes "   \
              "(IC, ND); Automatically-assigned Evidence Codes (IEA). Each "   \
              "code must be included separately. Ex: \"-i EXP -i IDA -i IPI\"."\
              " For all Experimental Evidence Codes, use \"-i ALL_EXP\". "       \
              "For all Computational Analysis Evidence Codes, use \"-i COMP_ALL\"."\
              " For all Evidence Codes, except IEA, use \"-i ABIEA\". "          \
              "For all Evidence Codes, use \"-i ALL\". "                       

    #Defines input variables
    parser = OptionParser(usage=usage, 
        version="%prog 0.8.9", 
        description=description)
    parser.add_option("-t", "--type", 
        dest="type", 
        help="Data type of the input files. Nucleotide (\"nt\"), " \
             "Aminoacid (\"aa\") or Protein Family (\"pf\").", 
        metavar="XX")
    parser.add_option("-f", "--file", 
        dest="file", 
        help="Full path for the FILE containing the data to be analized.", 
        metavar="FILE")
    parser.add_option("-d", "--database", 
        dest="dbdir", 
        help="Full path for the directory containing the database files. " \
             "Detailed information in \"README.txt\" file.", 
        metavar="/DIRECTORY/DIR/")
    parser.add_option("-o", "--outputdirectory", 
        dest="outdir", 
        help="Full path for the output directory.", 
        metavar="/DIRECTORY/DIR/")
    parser.add_option("-s", "--sifterdirectory", 
        dest="sdir", 
        help="Full path for the directory where SIFTER was prepared.", 
        metavar="/DIRECTORY/DIR/")
    parser.add_option("-i", 
        dest="experimental", 
        action="append", 
        default=[], 
        help=exphelp, 
        metavar="CODE")
    parser.add_option("-x", 
        dest="species", 
        action="append", 
        default=[], 
        help="(Optional) Remove specific species annotation from the analisis."\
             " The code must be the TaxonomyID. (More information on " \
             "http://www.ncbi.nlm.nih.gov/taxonomy)", 
        metavar="CODE")
    parser.add_option("-b", 
        dest="branch", 
        action="append", 
        default=[], 
        help="(Optional) Remove full taxon annotation from the analisis, " \
             "given the TaxonomyID code for the taxon's root. " \
             "(More information on http://www.ncbi.nlm.nih.gov/taxonomy)", 
        metavar="CODE")
    parser.add_option("--threads", 
        dest="threads", 
        default=multiprocessing.cpu_count(), 
        help="(Optional) Number of parallel CPU workers to use for " \
             "multithreads (default all).", 
        metavar="NUM", 
        type="int")
    parser.add_option("-r", "--translation-table", 
        dest="translation", 
        default="1", 
        help="(Optional) Translation table to be used (default - 1 - Standard" \
             " Code). For more information, visit " \
             "http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi.", 
        metavar="NUM", 
        type="int")
    parser.add_option("--sifter-cutoff", 
        dest="scut", 
        default="0.001", 
        help="(Optional) Show only sifter results with probability above" \
             " the cutoff (Default: 0.001).", 
        metavar="NUM", 
        type="float")
    parser.add_option("--pfamscan-cutoff", 
        dest="pcut", 
        default="0.00001", 
        help="(Optional) During \"pfam_scan.pl\" step, detect protein families"\
             " above the e-value cutoff (Default: 0.00001).", 
        metavar="NUM", 
        type="float")
    parser.add_option("--input_species", 
        dest="input_species", 
        default="0", 
        help="(Optional) In order to generate an apropriate Species Tree, "    \
             "input the NCBI Taxon Identification for the Species of the "     \
             "input data. (Default: No Species Tree generation)",
        metavar="NUM", 
        type="string")
    parser.add_option("--reconciliation", 
        dest="reconciliation", 
        default=False, 
        action="store_true", 
        help="(Optional) Enables reconciliation for nucleotide or aminoacid "  \
             "sequences as input. Require \"--input_species\". "               \
             "(Default: No reconciliation.)")
    (options, args) = parser.parse_args()
    print "\n###  --------------------------------------------------  ###" 
    print "###         Running SIFTER-T. Initial checking...        ###" 
    print "###  --------------------------------------------------  ###\n" 

    required_parameters(options)
    options = check_databases(options)
    options.file = os.path.abspath(options.file)
    check_type(options)
    check_input_file0(options)
    check_input_file1(options)
    print "# Checking directories...\n"
    options = check_outdir(options)
    options = check_sifter_dir(options)
    options.stdir = os.path.abspath(os.getcwd())+"/"
    check_sifter_files(options)
    check_pfam_scan(options)
    check_fasttree(options)
    check_hmmer3()
    check_mafft()
    check_notung(options)
    options = check_evcodes(options)
    print "# Checking optional parameters...\n"
    check_reconciliation(options)
    check_translation_table(options)
    check_sifter_cutoff(options)
    check_pfam_cutoff(options)
    check_species(options)
    check_threads(options)
    options.force = True

    ## Registering "options" variable
    pickle.dump(options, file(options.outdir+"options.pk", "w"))
    return options


def ntaapfam(options):
    '''
    Protein Family multiple alignment generation.
    '''
    from ntaapfam_01 import nt_prepare_input
    from ntaapfam_01 import aa_prepare_input
    from ntaapfam_01 import pfam_scan_mp
    from ntaapfam_01 import split_query_fasta
    from ntaapfam_01 import pfam_scan_multi
    from ntaapfam_01 import get_input_pfam_genes
    from ntaapfam_01 import get_pfam_genes
    from ntaapfam_01 import get_annot_genes_all
    from ntaapfam_01 import get_no_annot
    from ntaapfam_01 import get_no_annot2
    from ntaapfam_01 import get_sp_gene
    from ntaapfam_01 import get_sp_desc_anc
    from ntaapfam_01 import get_sp_branch_set
    from ntaapfam_01 import get_forbidden_sp_gene
    from ntaapfam_01 import get_input_genes_pfam__handle_pf
    from ntaapfam_01 import write_input_ntaa
    from ntaapfam_01 import write_fasta_pf
    from ntaapfam_01 import write_selected_pfam_genes
    from ntaapfam_01 import multi_write_selected_pfam_genes
    from ntaapfam_01 import clean_useful_pfam
    from ntaapfam_01 import align_sequences
    from ntaapfam_01 import write_pfam_list
    from ntaapfam_01 import clean_annot_genes_all

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
    return None



def recann(options):
    '''
    Annotation recovery and \".pli\" files preparation.
    '''
    from recann_02 import pfam2pli 
    from recann_02 import get_forbidden_genes
    from recann_02 import get_uniprots
    from recann_02 import get_familylist
    from recann_02 import get_families
    from recann_02 import pfam2pli_multi

    print "\n###  --------------------------------------------------  ###" 
    print "###                  Annotation Recovery                 ###" 
    print "###  --------------------------------------------------  ###\n" 
    pfam2pli_multi(options)
    return None



def gentree(options):
    '''
    Gene tree preparation.
    '''
    from gentree_03 import clean_tree 
    from gentree_03 import nseq
    from gentree_03 import pli2tree 
    from gentree_03 import gentree_multi

    print "\n###  --------------------------------------------------  ###" 
    print "###                Building Gene Trees                   ###" 
    print "###  --------------------------------------------------  ###\n" 

    familylist = list()
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
    return None

def reconciliation(options):
    '''
    Species tree creation and Protein Family reconciliation for Protein Families 
    with suitable size.
    '''
    from reconciliation_04 import get_species_set
    from reconciliation_04 import buid_species_tree
    from reconciliation_04 import build_reconciled_tree
    from reconciliation_04 import build_rec_multi

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
    return None


def runsifter(options):
    '''
    This script automate the inference engine execution.
    '''
    from runsifter_05 import run_sifter0 
    from runsifter_05 import run_sifter1 
    from runsifter_05 import sifter_multi

    print "\n###  --------------------------------------------------  ###" 
    print "###           Posterior probability calculation          ###" 
    print "###  --------------------------------------------------  ###\n" 
    
    #Defines sifter arguments related to evidence codes
    scodes = "--with-iea --with-ic --with-iep --with-igc --with-igi " \
             "--with-ipi --with-iss --with-rca --with-tas --with-nas"

    familylist = list()
    handle = open(options.outdir+"useful_pfam.txt","r")
    for line in handle:
        d = line.strip().split()
        familylist.append(d[0])
    handle.close()

    sifter_multi(options, familylist, scodes)
    return None


def results(options):
    '''
    Summarizes results from Sifter-T's workflow.
    '''
    from results_06 import get_querynames_conversion
    from results_06 import get_go_term_conversion
    from results_06 import get_query_pfam__query_pfam_len
    from results_06 import get_no_pfam
    from results_06 import get_noannot_goa
    from results_06 import get_noannot_ec
    from results_06 import get_noannot_all
    from results_06 import get_useful_pfam
    from results_06 import get_useful_pfam_func
    from results_06 import get_useful_pfam_query_prob
    from results_06 import get_query_len_pfam_prob
    from results_06 import get_memsize
    from results_06 import get_disk_used
    from results_06 import get_input_type
    from results_06 import write_report_alt_txt
    from results_06 import write_report_txt
    from results_06 import create_report_tab

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
        noannot_all = get_noannot_all(options, noannot_goa, 
                                      noannot_ec, query_pfam)
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
                     useful_pfam_query_prob, useful_pfam_func, 
                     go_term_conversion, input_type, memsize, disk_used, nopfam,
                     noannot_goa, noannot_ec, noannot_all, query_len_pfam_prob)
    create_report_tab(options, useful_pfam, querynames_conversion,
                      useful_pfam_query_prob, useful_pfam_func, 
                      go_term_conversion, input_type, memsize, disk_used, 
                      nopfam, noannot_goa, noannot_ec, noannot_all, 
                      query_len_pfam_prob)
    return None




def _main():
    '''
    Main function for standalone usage.
    '''
    options = param_chk()
    ntaapfam(options)
    recann(options)
    gentree(options)
    if options.reconciliation:
        reconciliation(options)
    runsifter(options)
    results(options)
    print "# The result summary is available at: \n"
    print "    "+options.outdir+"REPORT.txt\n"
    if options.type != "pf":
        print "    "+options.outdir+"REPORT_alt.txt\n"
    print "    "+options.outdir+"REPORT.tab\n"
    sys.exit()






if __name__ == '__main__':
    q = JoinableQueue()    
    n = JoinableQueue()    

    _main()

