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
#     Almeida-e-Silva D.C. and Vêncio R.Z.N. (2015) SIFTER-T: A scalable    #
#     and optimized framework for the SIFTER phylogenomic method of         #
#     probabilistic protein domain annotation. BioTechniques, Vol. 58,      #
#     No. 3, March 2015, pp. 140–142                                        #
#                                                                           #
########## ########## ########### ########## ########## ########## ##########

"""
 * Quick and dirty databases download and decompressing.

    _main()

"""

import os, sys
from optparse import OptionParser

def remove_former_db(options):
    '''
    Clean old files in case of "force".
    '''
    required_files = ["gene_association.goa_uniprot", "taxonomy.txt", 
                      "ncbi_taxonomy.obo", "delnodes.dmp", "merged.dmp", 
                      "uniprot_sprot.dat", "uniprot_trembl.dat", 
                      "gene_ontology.1_2.obo", "Pfam-A.hmm.dat", 
                      "Pfam-A.hmm", "Pfam-A.full"]
    gz_files = ["gene_association.goa_uniprot.gz", "taxonomy.txt.gz", 
                      "taxdump.tar.gz", 
                      "uniprot_sprot.dat.gz", "uniprot_trembl.dat.gz", 
                      "go-basic.obo", "Pfam-A.hmm.dat.gz", 
                      "Pfam-A.hmm.gz", "Pfam-A.full.gz"]
    if options.force:
        for dbfile in required_files:
            if os.path.exists(options.dir+dbfile):
                os.remove(options.dir+dbfile)
        for dbfile in gz_files:
            if os.path.exists(options.dir+dbfile):
                os.remove(options.dir+dbfile)


def download_files(options):
    '''
    Download Sifter-T's required database.
    '''    
    os.chdir(options.dir)
    print "# Downloading files...\n"

    #NCBI Taxonomy 1
    print "# Downloading ncbi_taxonomy.obo...\n"
    if not os.path.exists(options.dir+"ncbi_taxonomy.obo"):
        os.system("wget -q -nv http://www.berkeleybop.org/ontologies/obo-all/ncbi_taxonomy/ncbi_taxonomy.obo")

    #NCBI Taxonomy 2
    print "# Downloading taxdump.tar.gz...\n"
    if not os.path.exists(options.dir+"merged.dmp") or not os.path.exists(options.dir+"delnodes.dmp"):
        os.system("wget -q -nv ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz")
        os.system("tar -zxf taxdump.tar.gz merged.dmp delnodes.dmp")
        os.system("rm taxdump.tar.gz")

    #Gene Ontology
    print "# Downloading gene_ontology.1_2.obo...\n"
    if not os.path.exists(options.dir+"gene_ontology.1_2.obo"):
        os.system("wget -q -nv http://purl.obolibrary.org/obo/go/go-basic.obo")
        os.system("mv go-basic.obo gene_ontology.1_2.obo")

    #GOA
    print "# Downloading gene_association.goa_uniprot...\n"
    if not os.path.exists(options.dir+"gene_association.goa_uniprot"):
        os.system("wget -q -nv ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/gene_association.goa_uniprot.gz")
        os.system("gzip -d gene_association.goa_uniprot.gz &")

    #PFAM
    #pfam_site = "ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/" - Major changes on PFam v29+. Using v28 due to compatibility issues.
    pfam_site = "ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam28.0/"

    print "# Downloading Pfam-A.full...\n"
    if not os.path.exists(options.dir+"Pfam-A.full"):
        os.system("wget -q -nv "+pfam_site+"Pfam-A.full.gz")
        os.system("gzip -d Pfam-A.full.gz &")

    print "# Downloading Pfam-A.hmm.dat...\n"
    if not os.path.exists(options.dir+"Pfam-A.hmm.dat"):
        os.system("wget -q -nv "+pfam_site+"Pfam-A.hmm.dat.gz")
        os.system("gzip -d Pfam-A.hmm.dat.gz &")
 
    print "# Downloading Pfam-A.hmm...\n"
    if not os.path.exists(options.dir+"Pfam-A.hmm"):
        os.system("wget -q -nv "+pfam_site+"Pfam-A.hmm.gz")
        os.system("gzip -d Pfam-A.hmm.gz &")

    print "# Downloading uniprot_trembl.dat...\n"
    if not os.path.exists(options.dir+"uniprot_trembl.dat"):
        #uniprot_trembl.dat.gz not present on PFam v28 directory. Using alternative link.
        #os.system("wget -q -nv "+pfam_site+"uniprot_trembl.dat.gz")
        os.system("wget -q -nv ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.dat.gz")
        os.system("gzip -d uniprot_trembl.dat.gz")

    print "# Downloading uniprot_sprot.dat...\n"
    if not os.path.exists(options.dir+"uniprot_sprot.dat"):
        os.system("wget -q -nv "+pfam_site+"uniprot_sprot.dat.gz")
        os.system("gzip -d uniprot_sprot.dat.gz")

    print "# Downloading taxonomy.txt...\n"
    if not os.path.exists(options.dir+"taxonomy.txt"):
        os.system("wget -q -nv "+pfam_site+"database_files/taxonomy.txt.gz")
        os.system("gzip -d taxonomy.txt.gz")



def check_files(options):
    '''
    Check if databases were successfully downloaded.
    '''
    required_files = ["gene_association.goa_uniprot", "taxonomy.txt", 
                      "ncbi_taxonomy.obo", "delnodes.dmp", "merged.dmp", 
                      "uniprot_sprot.dat", "uniprot_trembl.dat", 
                      "gene_ontology.1_2.obo", "Pfam-A.hmm.dat", 
                      "Pfam-A.hmm", "Pfam-A.full"]
    for item in required_files:
        if not os.path.exists(options.dir+item):
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
                  "\nExiting..." % options.dir
            sys.exit(1)

    for item in required_files:
        if not os.access(options.dir+item, os.R_OK):
            print "The actual user does not have READING access to all "       \
                  "required files. Type \"-h\" for help. \nExiting..."
            sys.exit()



def _main():
    '''
    Main function for standalone usage.
    '''
    #Defines usage and help
    usage = "\n 	%prog -d DIR"
    description = "Quick and dirty Sifter-T's databases download and decompressing." 
    
    #Defines input variables
    parser = OptionParser(usage=usage, version="%prog 0.1.0", 
        description=description)
    parser.add_option("-d", "--directory", 
        dest="dir", 
        help="full path for the output directory.", 
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

    print "\n###  --------------------------------------------------  ###" 
    print "###          Database Download and Decompressing         ###" 
    print "###  --------------------------------------------------  ###\n" 

    print "# Checking directory...\n"
    if not os.path.isdir(options.dir):        
        print options.dir+" is not an existing directory."
        sys.exit()
    if not os.access(options.dir, os.R_OK):
        print "The actual user does not have READING access to %s directory."  \
              " \nExiting..." % options.dir
        sys.exit(1)
    options.dir = os.path.abspath(options.dir)
    if options.dir[-1] != "/":
        options.dir = options.dir+"/"    
    
    remove_former_db(options)
   
    download_files(options)
    
    check_files(options)

    print "\n###  ------------------------------------------------------  ###" 
    print "###   Databases Download completed.                          ###" 
    print "###   Now you can prepare the databases with \"dbprep.py\".  ###" 
    print "###  ------------------------------------------------------  ###\n" 
    
    sys.exit()

if __name__ == '__main__':
    _main()




