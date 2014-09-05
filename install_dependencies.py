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
 * Quick and dirty Sifter-T's software dependencies download and install.

    _main()

"""

import os, sys
from optparse import OptionParser


def check_dir(options):
    '''
    Check if it is Sifter-T's directory.
    '''
    siftert_files = ["dbprep.py", "download_db.py", "sifter-t.py", 
                      "gentree_03.py", "ntaapfam_01.py", "obo2flat.py", 
                      "param_chk_00.py", "recann_02.py", 
                      "recann_02.py", "reconciliation_04.py", 
                      "results_06.py", "runsifter_05.py"]
    for sfile in siftert_files:
        if not os.path.exists(options.dir+sfile):
            print "It appears that the selected directory does'n have Sifter-T installed. Type \"-h\" for help."\
                  " \nExiting..."
            sys.exit(1)


def remove_former_files(options):
    '''
    Clean old files in case of "force".
    '''
    dep_files = ["pfam_scan.pl", "FastTree"]
    dep_dirs = ["Bio", "notung", "PfamScan", "sifter"]
    if options.force:
        for depfile in dep_files:
            if os.path.exists(options.dir+depfile):
                os.remove(options.dir+depfile)
        import shutil
        for depdir in dep_dirs:
            if os.path.exists(options.dir+depdir):
                shutil.rmtree(os.path.abspath(options.dir+depdir))


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


def download_files(options):
    '''
    Download and install Sifter-T's software dependencies.
    '''
#    current_dir = options.dir.replace(" ","\ ")
    os.chdir(options.dir)
    #Biopython, Java, Bioperl
    print "# Downloading and installing Biopython, Java and Bioperl...\n"
    os.system("sudo apt-get -qq --force-yes install python-biopython python-pip openjdk-6-jre openjdk-7-jre bioperl wget tar unzip build-essential w3m gzip grep make sed")
    try: 
        from Bio.Alphabet import IUPAC
    except:
        raise Exception("Can't load Biophyton. Exiting...")

    #CPAN packages
    print "# Downloading and installing perl packages...\n"
    check_program("perl")
    os.system("sudo perl -MCPAN -e 'install CPAN'")
    os.system("sudo perl -MCPAN -e 'install Moose'")
    os.system("sudo perl -MCPAN -e 'install Data::Printer'")

    #Dentropy
    print "# Downloading and installing Dendropy...\n"
    check_program("pip")
    os.system("sudo pip install dendropy --quiet")

    try: 
        import dendropy
    except:
        raise Exception("Can't load Dendropy. Exiting...")

    #FastTree
    print "# Downloading and installing FastTree...\n"
    check_program("wget")
    os.system("wget -q -nv http://meta.microbesonline.org/fasttree/FastTree")
    os.system("chmod u+x FastTree")

    #PfamScan
    print "# Downloading and installing PfamScan...\n"
    check_program("tar")
    os.system("wget -q -nv ftp://ftp.sanger.ac.uk/pub/databases/Pfam/Tools/PfamScan.tar.gz") 
    os.system("tar -zxf PfamScan.tar.gz")
    os.system("rm PfamScan.tar.gz")
    os.system("mv PfamScan/Bio ./ ")
    os.system("mv PfamScan/pfam_scan.pl ./ ")

    #Sifter2.0
    print "# Downloading and installing Sifter2.0...\n"
    check_program("make")
    os.system("wget -q -nv http://sifter.berkeley.edu/code/sifter2.0.tar.gz")
    os.system("tar -zxf sifter2.0.tar.gz")
    os.system("rm sifter2.0.tar.gz")
    os.system("make clean -i -B -s -C sifter2.0/")
    os.system("make -i -B -s -C sifter2.0/")
    os.system("mv sifter2.0 sifter")

    #Notung
    print "# Downloading and installing Notung...\n"
    check_program("unzip")
    os.system("wget -q -nv http://lampetra.compbio.cs.cmu.edu/Notung/distributions/Notung-2.6.zip")
    os.system("unzip -qq Notung-2.6.zip")
    os.system("mv Notung-2.6 notung")
    os.system("mv notung/Notung-2.6.jar notung/Notung.jar")
    os.system("rm Notung-2.6.zip")

    check_program("w3m")
    check_program("sed")
    #Hmmer
    if not check_hmmer3() or options.force:
        print "# Downloading and installing Hmmer3...\n"
        os.system('''w3m ftp://selab.janelia.org/pub/software/hmmer3/CURRENT | grep -i linux | grep -i x86_64 | sed "s/ /\\t/g" | cut -f 1 > temp.txt''')
        handle = open("temp.txt", "r")
        files = handle.readlines()[0].strip()
        handle.close()
        os.system("rm temp.txt")
        os.system("wget -q -nv ftp://selab.janelia.org/pub/software/hmmer3/CURRENT/"+files)
        os.system("tar -zxf "+files)
        hmmerdir = files.replace(".tar.gz","")
        os.chdir(hmmerdir)
        os.system("./configure")
        os.system("make -B -s")
        os.system("make check -B -s")
        os.system("sudo make install -B -s")
        os.chdir("..")
        os.system("rm -R "+hmmerdir)
        os.system("rm "+files)

    #MAFFT
    if not check_mafft() or options.force:
        print "# Downloading and installing Mafft...\n"
        os.system('''w3m http://mafft.cbrc.jp/alignment/software/changelog.html | grep "^v" | head -n 1 | sed "s/ /\\t/g" | cut -f 1 | sed "s/v//g" > temp.txt''')
        handle = open("temp.txt", "r")
        version = handle.readlines()[0].strip()
        handle.close()
        os.system("rm temp.txt")
        os.system("wget -q -nv http://mafft.cbrc.jp/alignment/software/mafft-"+version+"-without-extensions-src.tgz")
        os.system("tar -zxf mafft-"+version+"-without-extensions-src.tgz")
        mafftdir = "mafft-"+version+"-without-extensions"
        os.system("make clean -i -B -s -C "+mafftdir+"/core/")
        os.system("make -i -B -s -C "+mafftdir+"/core/")
        os.system("sudo make install -i -B -s -C "+mafftdir+"/core/")
        os.system("rm mafft-"+version+"-without-extensions-src.tgz")
        os.system("sudo rm -R mafft-"+version+"-without-extensions")

   

def _main():
    '''
    Main function for standalone usage.
    '''
    #Defines usage and help
    usage = "\n 	%prog -d DIR"
    description = "Quick and dirty Sifter-T's software dependencies download" \
                  " and install." 
    
    #Defines input variables
    parser = OptionParser(usage=usage, version="%prog 0.1.0", 
        description=description)
    parser.add_option("-d", "--directory", 
        dest="dir", 
        help="(Optional) Full path for Sifter-T's directory.", 
        metavar="/DIRECTORY/DIR/",
        default = os.getcwd())
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
    print "###          Dependencies Download and Install           ###" 
    print "###  --------------------------------------------------  ###\n" 

    print "# Checking directory...\n"
    if not os.path.isdir(options.dir):        
        print options.dir+" is not an existing directory."
        sys.exit()
    if not os.access(options.dir, os.R_OK):
        print "The actual user does not have READING access to %s directory."  \
              " \nExiting..." % options.dir
        sys.exit(1)
    if not os.access(options.dir, os.W_OK):
        print "The actual user does not have READING access to %s directory."  \
              " \nExiting..." % options.dir
        sys.exit(1)
    options.dir = os.path.abspath(options.dir)
    if options.dir[-1] != "/":
        options.dir = options.dir+"/"    
    
    check_dir(options)
   ###############
    remove_former_files(options)
    download_files(options)
    
    print "\n###  ------------------------------------------------------  ###" 
    print "###   Dependencies installation is complete.                 ###" 
    print "###  ------------------------------------------------------  ###\n" 
    sys.exit()

if __name__ == '__main__':
    _main()

