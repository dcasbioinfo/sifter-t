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
 * This script checks main parameters, files and creates "options.pk" for 
 * other scripts.

    required_parameters(options)
    check_type(options)
    check_input_file0(options)
    is_number(num)
    check_input_file1(options)
    check_outdir(options)
    check_sifter_dir(options)
    check_sifter_files(options)
    check_pfam_scan(options)
    check_fasttree(options)
    check_notung(options)
    check_program(name)
    check_hmmer3()
    check_mafft()
    check_evcodes(options)
    check_translation_table(options)
    check_sifter_cutoff(options)
    check_pfam_cutoff(options)
    check_reconciliation(options)
    check_species(options)
    check_threads(options)
    check_databases(options)
    _main()

"""


from optparse import OptionParser
import pickle
import os
import sys
import multiprocessing
from Bio.Alphabet import IUPAC
try: 
    from Bio.Alphabet import IUPAC
except:
    raise Exception("Can't load Biophyton. Exiting...")

try: 
    import dendropy
except:
    raise Exception("Can't load Dendropy. Exiting...")


def required_parameters(options):
    '''
    Check for required parameters.
    '''
    if not options.outdir or \
       not options.dbdir or \
       not options.sdir or \
       not options.experimental or \
       not options.type or \
       not options.file:
        print "Not all required parameters were specified. "                   \
              "Type \"-h\" for help. \nExiting..."
        sys.exit(1)


def check_type(options):
    '''
    Check if "options.type" is a valid input type.
    '''
    options.type = options.type.lower()
    if options.type not in ["nt", "aa", "pf"]:
        print "%s is not a valid Data Type. Type \"-h\" for help."             \
              " \nExiting..." % options.type
        sys.exit(1)


def check_input_file0(options):
    '''
    Check input file existence and reading access.
    '''
    print "# Checking input file...\n"
    if not os.path.isfile(options.file):        
        print options.file+" is not an existing file."
        sys.exit()
    if not os.access(options.file, os.R_OK):
        print "The actual user does not have READING access to %s file."       \
              " \nExiting..." % options.file
        sys.exit(1)
    

def is_number(num):
    '''
    Return True if number.
    '''
    try:
        int(num)
        return True
    except ValueError:
        return False


def check_input_file1(options):
    '''
    Check if type match input file content.
    '''
    if options.type == "pf":
        pfam_set = set()
        with open(options.dbdir+"pfam.list") as handle:
            for line in handle:
                pfam_set.add(line.strip().split()[0])
        handle = open(options.file,"r")
        for line in handle:
            tmp_str = line.strip().split()
            if not tmp_str[0][:2].lower() == "pf":
                print "The file %s is malformed. All protein families must "   \
                      "start with 'PF' or 'pf', followed by five"              \
                      " numbers." % options.file
                sys.exit(1)
            if not is_number(tmp_str[0][2:]):
                print "The file %s is malformed. All protein families must "   \
                      "start with 'PF' or 'pf', followed by five "             \
                      "numbers." % options.file
                sys.exit(1)
            if not tmp_str in pfam_set:
                print "The file %s is not valid. %s is not present on PFAM"    \
                      " database." % (options.file, tmp_str)
                sys.exit(1)
    if options.type == "nt" or options.type == "aa" :
        from Bio import SeqIO
        handle = open(options.file, "rU")
        if options.type == "nt":
            wanted = set(IUPAC.IUPACAmbiguousDNA.letters+IUPAC.IUPACAmbiguousRNA.letters)
            for nuc_rec in SeqIO.parse(handle, "fasta"):
                if not wanted.issuperset(nuc_rec.seq):
                    print "Nucleotide input sequence malformed: %s "           \
                          "\nExiting..." % nuc_rec.id
                    sys.exit(1)
        elif options.type == "aa":
            wanted = set(IUPAC.ExtendedIUPACProtein.letters)
            for nuc_rec in SeqIO.parse(handle, "fasta"):
                if not wanted.issuperset(nuc_rec.seq):
                    print "Aminoacid input sequence malformed: %s "            \
                          "\nExiting..." % nuc_rec.id
                    sys.exit(1)


def check_outdir(options):
    '''
    Check if output directory exists and if the actual user have reading and 
    wrinting access.
    '''
    if not os.path.isdir(options.outdir):        
        print options.outdir+" is not an existing directory.\n"
        create = raw_input("Create directory? [y/n]: ")
        if create == "y":
            try:
                os.makedirs(options.outdir)
            except:
                raise IOError("\nCould not create the specified directory."    \
                              " \nExiting...")
            print ""
        elif create == "n":
            print "\nExiting..."
            sys.exit()
        elif create != "y" or create != "n":
            print "Not a valid option. \nExiting..."
            sys.exit()
    if not os.access(options.outdir, os.R_OK) or \
       not os.access(options.outdir, os.W_OK):
        print "The actual user does not have READING and/or WRITE access to"   \
              " %s directory. \nExiting..." % options.outdir
        sys.exit(1)
    options.outdir = os.path.abspath(options.outdir)
    if options.outdir[-1] != "/":
        options.outdir = options.outdir+"/"
    return options


def check_sifter_dir(options):
    '''
    Check if sifter2.0 directory exists and the actual user have reading
    access.
    '''
    if not os.path.isdir(options.sdir):        
        print options.sdir+" is not an existing directory."
        sys.exit()
    if not os.access(options.sdir, os.R_OK) or \
       not os.access(options.sdir, os.W_OK):
        print "The actual user does not have READING and/or WRITE access to"   \
              " %s directory. \nExiting..." % options.sdir
        sys.exit(1)
    options.sdir = os.path.abspath(options.sdir)
    if options.sdir[-1] != "/":
        options.sdir = options.sdir+"/"
    return options


def check_sifter_files(options):
    '''
    Check if sifter2.0 files are on the sifter directory, and check if sifter 
    is correctly installed.
    '''
    print "# Checking SIFTER consistency...\n"
    #checking SIFTER files
    handle = open("sifter-t_chk/sifter2.0_files.list1","r")
    for line in handle:
        tmp_str = line.strip().split()
        if not os.path.exists(options.sdir+tmp_str[0]):
            print "Sifter 2.0 is not present on %s. \nExiting..." % options.sdir
            sys.exit(1)
    handle.close()
    #checking SIFTER proper instalation
    handle = open("sifter-t_chk/sifter2.0_files.list2","r")
    for line in handle:
        tmp_str = line.strip().split()
        if not os.path.exists(options.sdir+tmp_str[0]):
            print "Sifter 2.0 is not correctly installed on %s. "              \
                  "\nExiting..." % options.sdir
            sys.exit(1)
    handle.close()


def check_pfam_scan(options):
    '''
    Check if standalone PfamScan is present on Sifter-T directory.
    '''
    print "# Checking PFamScan...\n"
    handle = open(options.stdir+"sifter-t_chk/pfam_scan_files.list","r")
    for line in handle:
        tmp_str = line.strip().split()
        if not os.path.exists(tmp_str[0]):
            print "PfamScan is not present on the actual directory. It must be"\
                  " extracted on the same directory as Sifter-T's scripts."    \
                  " \nExiting..."
            sys.exit(1)
    handle.close()


def check_fasttree(options):
    '''
    Check if FastTree is present on Sifter-T directory.
    '''
    print "# Checking FastTree...\n"
    if not os.path.isfile(options.stdir+"FastTree"):
        print "FastTree is not present on the actual directory. This file must"\
              " be on the same directory as Sifter-T's scripts. \nExiting..."
        sys.exit(1)
    if not os.access(options.stdir+"FastTree", os.X_OK):
        print "The actual user does not have EXECUTING access to FastTree. "   \
              "\nExiting..." 
        sys.exit(1)


def check_notung(options):
    '''
    Check if Notung is present on Sifter-T directory.
    '''
    print "# Checking Notung...\n"
    if not os.path.isfile(options.stdir+"Notung.jar"):
        print "Notung.jar is not present on the actual directory. This file "  \
              "must be on the same directory as Sifter-T's scripts. \nExiting..."
        sys.exit(1)
    if not os.access(options.stdir+"Notung.jar", os.R_OK):
        print "The actual user does not have READING access to Notung. "   \
              "\nExiting..." 
        sys.exit(1)


def check_program(name):
    '''
    Check if "name" is on $PATH
    '''
    for directory in os.environ['PATH'].split(':'):
        prog = os.path.join(directory, name)
        if os.path.exists(prog):
            return prog


def check_hmmer3():
    '''
    Check if hmmer3 is installed on the system.
    '''
    print "# Checking HMMER3...\n"
    if not check_program("hmmpress") or \
       not check_program("hmmalign") or \
       not check_program("hmmbuild") or \
       not check_program("hmmconvert") or \
       not check_program("hmmemit") or \
       not check_program("hmmfetch") or \
       not check_program("hmmscan") or \
       not check_program("hmmsearch") or \
       not check_program("phmmer") or \
       not check_program("jackhmmer"):
        print "HMMER3 is not installed or is not on $PATH. \nExiting..." 
        sys.exit(1)


def check_mafft():
    '''
    Check if mafft is installed on the system.
    '''
    print "# Checking MAFFT...\n"
    if not check_program("mafft"):
        print "MAFFT is not installed or is not on $PATH. \nExiting..." 
        sys.exit(1)


def check_evcodes(options):
    '''
    Check if the evidence codes are valid.
    '''
    print "# Checking evidence codes...\n"
    codes = ["TAS", "IDA", "IMP", "IGI", "IPI", "EXP", "ISO", "ISA", "ISM",
             "ISS", "IBA", "IBD", "IEP", "IGC", "RCA", "IC", "NAS", "ND", "IEA"]
    exp_set = set()
    for exp in list(options.experimental):
        if exp not in codes:
            if (exp == "ALL"):
                exp_set = (exp_set | set(codes)) - set(["ALL"])
            elif (exp == "ABIEA"):
                exp_set = exp_set | (set(codes) - set(["IEA", "ABIEA"]))
            elif (exp == "EXP_ALL"):
                exp_set = (exp_set | set(["EXP", "IDA", "IPI", "IMP", "IGI", 
                                          "IEP"]))- set(["EXP_ALL"])
            elif (exp == "COMP_ALL"):
                exp_set = (exp_set | set(["ISS", "ISO", "ISA", "ISM", "IGC", 
                                      "RCA", "IBD", "IBA"])) - set(["COMP_ALL"])
            else:
                print "%s is not a valid Evidence Code. Type \"-h\" for help." \
                      " \nExiting..." % exp
                sys.exit(1)
    options.experimental = list(exp_set)
    return options


def check_translation_table(options):
    '''
    Check if translation table selected is valid.
    '''
    transtable = [1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15, 16, 21, 22, 23]
    if options.translation not in transtable:
        print "%i is not a valid Translation Table. Type \"-h\" for help."     \
              " \nExiting..." % options.translation
        sys.exit(1)


def check_sifter_cutoff(options):
    '''
    Check if Sifter2.0 cutoff selected is valid.
    '''
    if options.scut < 0 or options.scut >= 1:
        print "%s is not a valid Sifter cutoff. Must be between 0 and 1. "     \
              "Type \"-h\" for help. \nExiting..." % str(options.scut)
        sys.exit(1)
    

def check_pfam_cutoff(options):
    '''
    Check if PfamScan cutoff selected is valid.
    '''
    if options.pcut <= 0:
        print "%s is not a valid \"pfam_scan\" e-value cutoff. Must be a "     \
              "positive value. Type \"-h\" for help. " \
              "\nExiting..." % str(options.pcut)
        sys.exit(1)


def check_reconciliation(options):
    '''
    If reconciliation is true, check if "--input_species" is also used.
    '''
    if options.reconciliation and options.input_species == "0":
        print "%s is not a valid \"--input_species\" value while "             \
              "reconciliation is on. Type \"-h\" for help. "                   \
              "\nExiting..." % str(options.input_species)
        sys.exit(1)
    if options.reconciliation and options.type == "pf":
        print "Reconciliation is not yet suitable for Protein Family list as"  \
              " input. Type \"-h\" for help. "                   \
              "\nExiting..." % str(options.input_species)
        sys.exit(1)

def check_species(options):
    '''
    Check if all NCBI_Tax number is wellformed and valid.
    '''
    for taxid in (set(options.species) | set(options.branch) | set([options.input_species])):
        if not is_number(taxid):
            print "The TaxonomyID %s is malformed. It must contain only "      \
                  "numbers." % taxid
            sys.exit(1)
    if len(set(options.species) | set(options.branch) | set([options.input_species])) > 0:
        sp_set = set()
        with open(options.dbdir+"summary_ncbi_taxonomy.obo") as handle:
            for line in handle:
                sp_set.add(line.strip().split()[0])
                sp_set.add(line.strip().split()[1])
        with open(options.dbdir+"summary_taxonomy.txt") as handle:
            for line in handle:
                sp_set.add(line.strip().split()[0])
        for taxid in (set(options.species) | set(options.branch) | set([options.input_species])):
            if taxid not in sp_set:
                print "The TaxonomyID %s is not a valid NCBI Taxonomy number." % taxid
                sys.exit(1)


def check_threads(options):
    '''
    Check if number of threads selected is valid.
    '''
    if options.threads < 0 or options.threads > multiprocessing.cpu_count():
        print "%s is not a valid number of threads. Must be a value between 1" \
              " and %s. Type \"-h\" for help. \nExiting..." \
              "" % (str(options.threads), str(multiprocessing.cpu_count()))
        sys.exit(1)


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
                os.system("python dbprep.py -d "+options.dbdir+" -f")
            elif prepare == "n":
                print "\nExiting..."
                sys.exit(1)
            elif prepare != "y" or prepare != "n":
                print "Not a valid option. \nExiting..."
    return options



def _main():
    '''
    Main function for standalone usage.
    '''
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
              "For all Experimental Evidence Codes, use \"-i ALL_EXP\"."       \
              "For all Computational Analysis Evidence Codes, use \"-i COMP_ALL\"."\
              "For all Evidence Codes, except IEA, use \"-i ABIEA\". "          \
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

    if len(args) > 0:
        print "\n# Extra arguments. Wrong usage. Exiting... \n"
        sys.exit(1)

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


    ## Registering "options" variable
    pickle.dump(options, file(options.outdir+"options.pk", "w"))
    sys.exit()


if __name__ == '__main__':
    _main()

