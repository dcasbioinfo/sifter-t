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

import os, sys, urllib
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
    return False

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
        return False
    else:
        return True

def check_mafft():
    '''
    Check if mafft is installed on the system.
    '''
    print "# Checking MAFFT...\n"
    if not check_program("mafft"):
        return False
    else:
        return True


def download_files(options):
    '''
    Download and install Sifter-T's software dependencies.
    '''
#    current_dir = options.dir.replace(" ","\ ")
    os.chdir(options.dir)
    #Biopython, Java, Bioperl
    print "\n# Downloading and installing Biopython, Java and Bioperl...\n"
    if options.pause:
        raw_input("Press Enter to continue...")
    if not check_program("apt-get"):
        print "'apt-get' is not installed or is not on $PATH. \nExiting..." 
        sys.exit(1)
    if not check_program("sudo"):
        print "'sudo' is not installed or is not on $PATH. \nExiting..." 
        sys.exit(1)
    os.system("sudo apt-get -qq --force-yes update")
    os.system("sudo apt-get -qq --force-yes install python python-biopython python-pip openjdk-6-jre openjdk-7-jre openjdk-6-jdk openjdk-7-jdk perl bioperl tar unzip build-essential w3m grep gzip make sed python-dev libpython-dev libevent-dev")
    if not check_program("javac"):
        print "'javac' (package openjdk-6-jdk) was not correctly installed or is not on $PATH. It is required for Sifter installation. \nExiting..." 
        sys.exit(1)
    try: 
        from Bio.Alphabet import IUPAC
    except:
        raise Exception("Can't load Biophyton. Exiting...")

    if not check_program("tar"):
        print "'tar' (package tar) was not correctly installed or is not on $PATH. It is required for automatic installation. \nExiting..." 
        sys.exit(1)

    if not check_program("make"):
        print "'make' (package make) was not correctly installed or is not on $PATH. It is required for automatic installation. \nExiting..." 
        sys.exit(1)

    if not check_program("unzip"):
        print "'unzip' (package unzip) was not correctly installed or is not on $PATH. It is required for automatic installation. \nExiting..." 
        sys.exit(1)

    if not check_program("sed"):
        print "'sed' (package sed) was not correctly installed or is not on $PATH. It is required for automatic installation. \nExiting..." 
        sys.exit(1)

    if not check_program("perl"):
        print "'perl' (package perl) was not correctly installed or is not on $PATH. It is required for automatic installation. \nExiting..." 
        sys.exit(1)

    if not check_program("w3m"):
        print "'w3m' (package w3m) was not correctly installed or is not on $PATH. It is required for automatic HMMER 3 installation. \nExiting..." 
        sys.exit(1)


    #CPAN packages
    print "\n# Downloading and installing perl packages...\n"
    if options.pause:
        raw_input("Press Enter to continue...")
    #os.system("sudo perl -MCPAN -e 'install CPAN'")
    if not os.system("perl -MMoose -e 1"):
        os.system("sudo perl -MCPAN -e 'install Moose'")
    if not os.system("perl -MData::Printer -e 1"):
        os.system("sudo perl -MCPAN -e 'install Data::Printer'")

    #Dentropy
    print "\n# Downloading and installing Dendropy...\n"
    if options.pause:
        raw_input("Press Enter to continue...")
    if not check_program("pip"):
        print "'pip' (package python-pip) was not correctly installed or is not on $PATH. It is required for automatic Dendropy installation. \nExiting..." 
        sys.exit(1)
    os.system("sudo pip install dendropy==3.12.3 --quiet")

    try: 
        import dendropy
    except:
        raise Exception("Can't load Dendropy. Exiting...")

    #FastTree
    print "\n# Downloading and installing FastTree...\n"
    if options.pause:
        raw_input("Press Enter to continue...")
    urllib.urlretrieve("http://meta.microbesonline.org/fasttree/FastTree", "FastTree")
    os.system("chmod u+x FastTree")

    #PfamScan
    print "\n# Downloading and installing PfamScan...\n"
    if options.pause:
        raw_input("Press Enter to continue...")
    urllib.urlretrieve("ftp://ftp.sanger.ac.uk/pub/databases/Pfam/Tools/PfamScan.tar.gz", "PfamScan.tar.gz")
    os.system("tar -zxf PfamScan.tar.gz")
    os.system("rm PfamScan.tar.gz")
    os.system("mv PfamScan/Bio ./ ")
    os.system("mv PfamScan/pfam_scan.pl ./ ")

    #Sifter2.0
    print "\n# Downloading and installing Sifter2.0...\n"
    if options.pause:
        raw_input("Press Enter to continue...")
    if not os.path.exists("sifter/sifter.jar") or options.force:
        urllib.urlretrieve("http://genome.cshlp.org/content/suppl/2011/07/21/gr.104687.109.DC1/sifter2.0.tar.gz", "sifter2.0.tar.gz")
        os.system("tar -zxf sifter2.0.tar.gz")
        os.system("rm sifter2.0.tar.gz")
        os.system("make -i -B -s -C sifter2.0/")
        os.system("mv sifter2.0 sifter")

    #Notung
    print "\n# Downloading and installing Notung...\n"
    if options.pause:
        raw_input("Press Enter to continue...")
    if not os.path.exists("notung/Notung.jar") or options.force:
        urllib.urlretrieve("http://lampetra.compbio.cs.cmu.edu/Notung/distributions/Notung-2.6.zip", "Notung-2.6.zip")
        os.system("unzip -qq Notung-2.6.zip")
        os.system("mv Notung-2.6 notung")
        os.system("mv notung/Notung-2.6.jar notung/Notung.jar")
        os.system("rm Notung-2.6.zip")

    #Hmmer
    if not check_hmmer3() or options.force:
        print "\n# Downloading and installing Hmmer3...\n"
        if options.pause:
            raw_input("Press Enter to continue...")
        os.system('''w3m ftp://selab.janelia.org/pub/software/hmmer3/CURRENT | grep -i linux | grep -i x86_64 | sed "s/ /\\t/g" | cut -f 1 > temp.txt''')
        handle = open("temp.txt", "r")
        files = handle.readlines()[0].strip()
        handle.close()
        os.system("rm temp.txt")
        urllib.urlretrieve("ftp://selab.janelia.org/pub/software/hmmer3/CURRENT/"+files, files)
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
        print "\n# Downloading and installing Mafft...\n"
        if options.pause:
            raw_input("Press Enter to continue...")
        os.system('''w3m http://mafft.cbrc.jp/alignment/software/changelog.html | grep "^v" | grep -v "only" | head -n 1 | sed "s/ /\\t/g" | cut -f 1 | sed "s/v//g" > temp.txt''')
        handle = open("temp.txt", "r")
        version = handle.readlines()[0].strip()
        handle.close()
        os.system("rm temp.txt")
        urllib.urlretrieve("http://mafft.cbrc.jp/alignment/software/mafft-"+version+"-without-extensions-src.tgz", "mafft-"+version+"-without-extensions-src.tgz")
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
    parser.add_option("-p", "--pause", 
        help="(Optional) Pause between installations. Helps track potential error messages. (default False)", 
        action="store_true", 
        dest="pause", 
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
    remove_former_files(options)
    download_files(options)
    
    print "\n###  ------------------------------------------------------  ###" 
    print "###   Dependencies installation is complete.                 ###" 
    print "###  ------------------------------------------------------  ###\n" 
    sys.exit()

if __name__ == '__main__':
    _main()

