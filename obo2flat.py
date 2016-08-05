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
 * This script converts Gene Ontology's obo ontology file to 
 * legacy Molecular Funcion GO Flat File.

    import_obogo(filename)
    create_child_nodes(go_tree)
    deep_search(stack, last_h)
    class Stack  
    class go_node

"""

# old format definitions: ftp://ftp.geneontology.org/go/www/GO.format.go.shtml
# new format definitions: http://www.geneontology.org/GO.format.obo-1_2.shtml

from optparse import OptionParser
import os
import sys

class Stack: 
    '''
    Defines a stack with some extra functions.
    '''
    def __init__(self):  
        self.stack = []  
    def push(self, variable):
        ''' Push object in the stack. '''
        self.stack.append(variable)  
    def last(self):  
        ''' Return last object of the stack. '''
        return(self.stack[len(self.stack)-1])
    def pop(self):  
        ''' Pop object from the stack. '''
        if len(self.stack) == 0:  
            raise Exception("Error", "stack is empty")
        obj = self.stack[-1]  
        del self.stack[-1]  
        return obj  
    def isempty(self):  
        ''' Check if the stack is empty. '''
        if len(self.stack) == 0:  
            return True  
        else:  
            return False  
    def lenght(self):
        ''' Return the lenght of the stack. '''
        return len(self.stack)



class go_node:
    '''
    Defines node structure for each functional term in GO OBO File.
    '''
    def __init__(self): 
        self.go_id = str()
        self.name = str()    
        self.is_a = dict()
        self.namespace = str()
        self.alt_id = set()
        self.child_nodes = set()
        self.is_obsolete = False    
    def empty(self):
        ''' Return go_node to initial state. '''
        self.go_id = str()
        self.name = str()    
        self.is_a = dict()
        self.namespace = str()
        self.alt_id = set()
        self.child_nodes = set()
        self.is_obsolete = False
    def copy(self):
        ''' Return a copy of go_node. '''
        temp = go_node()
        temp.go_id = self.go_id
        temp.name = self.name 
        temp.is_a = self.is_a.copy()
        temp.namespace = self.namespace
        temp.alt_id = self.alt_id.copy()
        temp.child_nodes = self.child_nodes.copy()
        temp.is_obsolete = self.is_obsolete
        return temp



def import_obogo(filename):
    '''
    Import obo file to a tree-like structure.
    '''
    gene_ontology_obo = open(filename,"r")
    term = False
    go_tree = dict()
    temp_node = go_node()
    empty_node = go_node()
    for line in gene_ontology_obo:
        if not term:
            if "date:" in line:
                go_tree["date"] = line[6:len(line)-1]
            if "remark: cvs version:" in line:
                go_tree["version"] = line[21:len(line)-1]
            if "saved-by: " in line:
                go_tree["saved-by"] = line[10:len(line)-1]
            if "[Term]" in line:
                term = True
        if "[Term]" in line or "[Typedef]" in line:
            if not temp_node.is_obsolete \
                    and "molecular" in temp_node.namespace \
                    and "function" in temp_node.namespace:
                go_tree[temp_node.go_id] = temp_node.copy()
            temp_node = empty_node.copy()
        elif line[0:4] == "id: ":
            temp_node.go_id = line[line.index(" ")+1:len(line)-1]
        elif "name: " in line:
            temp_node.name = line[line.index(" ")+1:len(line)-1]
        elif "alt_id: " in line:
            temp_node.alt_id.add(line[line.index(" ")+1:len(line)-1])
        elif "namespace: " in line:
            temp_node.namespace = line[line.index(" ")+1:len(line)-1]
        elif "is_a: " in line:
            temp_node.is_a[line[line.index(" ")+1:line.index("!")-1]] = line[line.index("!")+2:len(line)-1]
        elif "is_obsolete: " in line:
            temp_node.is_obsolete = True
    if not temp_node.is_obsolete \
            and "molecular" in temp_node.namespace \
            and "function" in temp_node.namespace:
        go_tree[temp_node.go_id] = temp_node.copy()
    gene_ontology_obo.close()
    return go_tree


#create child conections
def create_child_nodes(go_tree):
    '''
    Create child nodes on the tree.
    '''
    for item in set(go_tree.keys()):
        if "GO" in item:
            for item2 in set(go_tree.keys()):
                if "GO" in item2 and item in set(go_tree[item2].is_a.keys()):
                    go_tree[item].child_nodes.add(item2)
    return go_tree


def deep_search(stack, last_h, go_tree, outfile):
    '''
    Perform a deep search along the ontology and store in output file.
    '''
    last_object = stack.last()
    h = stack.lenght()
    if h > last_h:
        term_name_id = go_tree[last_object].name.replace(",","\,").replace(":","\:")+" ; "+go_tree[last_object].go_id
        term_alt_ids = ""
        for alt_id in go_tree[last_object].alt_id:
            term_alt_ids = "%s, %s" % (term_alt_ids, alt_id)
        space = (stack.lenght())*" "
        alt_is_a = ""
        for is_a in go_tree[last_object].is_a.keys():
            if is_a != stack.stack[len(stack.stack)-2]:
                alt_is_a = alt_is_a+" % "+go_tree[last_object].is_a[is_a].replace(",","\,")+" ; "+is_a
        outfile.write(space+"%"+term_name_id+term_alt_ids+alt_is_a+"\n")
    for item in go_tree[last_object].child_nodes:
        stack.push(item)
        last_h = h
        deep_search(stack, last_h, go_tree, outfile)
        stack.pop()


def _main():
    parser = OptionParser(usage="\n 	%prog -i INPUT -o OUTPUT", 
        version="%prog 0.0.2", 
        description="This script converts Gene Ontology's obo ontology file" \
                    " to legacy Molecular Funcion GO Flat File.")
    parser.add_option("-i", 
        dest="input", 
        help="Input file. Must be the Ontology File in OBO v1.2 format.", 
        metavar="FILE")
    parser.add_option("-o", 
        dest="output", 
        help="Output file.", 
        metavar="FILE")
    parser.add_option("-f", "--force", 
        help="(Optional) Force file substitution. (default False)", 
        action="store_true", 
        dest="force", 
        default = False)
    (options, args) = parser.parse_args()
    
    if (bool(options.input) + bool(options.output)) != 2:
        print "\n    No input and/or output file selected."\
              " Type \"-h\" for help.\n"
        exit(1)
    
    options.input = os.path.abspath(options.input)
    options.output = os.path.abspath(options.output)
    
    if not os.path.exists(options.input) or not os.path.isfile(options.input):
        print "\n    Input file does not exists.\n"
        exit(1)
    
    if options.force and os.path.exists(options.output):
        os.remove(options.output)
    
    elif os.path.exists(options.output) \
            or os.path.isfile(options.output) \
            and not options.force:
        print "\n    Output file \""+options.output+"\" already exists."
        overwrite = raw_input('    Overwrite? [y/n] ')
        if overwrite.lower() != "y":
            print "    Choose another output file. Exiting...\n"
            exit(1)
        else:
            os.remove(options.output)

    #import
    go_tree = import_obogo(options.input)
    go_tree = create_child_nodes(go_tree)

    #export
    outfile = open(options.output,"w")
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
    sys.exit()


if __name__ == '__main__':
    _main()


