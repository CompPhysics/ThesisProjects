#Turbomole converter 
#Converts basis files in turbomole format to C++ files readable by Fermion Mingle

import glob
import os

def extract_m_basisinfo(filename, printing = False):
    #Extract multiple basises from a turbomole file
    swtch = 0
    f = open(filename)
    
    comments = []
    basis_sets = [] #list containing all basis sets
    basis = []
    basisW = []
    basisE = []
    basisT = []
    contractedtype = []
    activecontracted = -1
    for i in f:
        I = i.split()
        if len(I)==0:
            continue
        if I[0] == "*":
            continue #ignore wildcards
        if I[0] == '#':        
            comments.append("//")
            comments[-1] += i
            #for e in I:
            #    comments[-1] += e + " "
            continue
        #print I
        
        if I[0][0] == "$":
            continue
        if I[0] in ["h", "he", "li", "be", "b", "c", "n", "o", "f", "ne", "na" ,"mg", "al", "si", "p", "s", "cl", "ar", "k", "ca", "sc", "ti", "v", "cr", "mn", "fe", "co", "ni", "cu", "zn", "ga", "ge", "as", "se","br", "kr", "rb", "sr", "y", "zr", "nb", "mo", "tc", "ru", "rh", "pd", "ag", "cd", "in", "sn", "sb", "te", "i", "xe"]:
            #If not first time here, put current basis into list
            if swtch == 1:
                basis_sets.append([bT, basisW, basisE, contractedtype, len(contractedtype)])
                basisW = []
                basisE = []
                contractedtype = []
                activecontracted = -1
            swtch = 1

            #Begin collecting data
            bT = I[1]+"_"+I[0]
            bT = bT.replace("-", "_") #basistype
            bT = bT.replace("(", "") #basistype
            bT = bT.replace(")", "") #basistype                        
            bT = bT.replace(",", "") #basistype                                    
            continue
        if I[0] in "123456789":
            if I[1] == "s":
                #print "Identified an s orbital."
                contractedtype.append(0)
            if I[1] == "p":
                #print "Identified an p orbital."
                contractedtype.append(1)
            if I[1] == "d":
                #print "Identified an d orbital."
                contractedtype.append(2)
            if I[1] == "f":
                #print "Identified an f orbital."
                contractedtype.append(3) 
            #Create contracted
            basis.append([])
            basisE.append([])
            basisW.append([])
            activecontracted += 1   
            #print I
        else:
            try:
                w = float(I[0].replace("D", "E"))
                e = float(I[1].replace("D", "E"))
                basisE[activecontracted].append(e)
                basisW[activecontracted].append(w)
                #basisE[activecontracted].append(float(I[1]))
                #basisW[activecontracted].append(float(I[0]))
            except:
                comments.append("//Ignored the following information from file:")
                for e in I:
                    comments[-1] += e + " "
                #comments[-1] += "Ignored the following information:" + I
                #print "Failed to read weigths and exponents."
                #print I
        #print I   
    if printing:
        print comments
        print contractedtype
        print len(contractedtype)
        print basisE
        print basisW
    f.close()
    #return comments, contractedtype, len(contractedtype), basisE, basisW
    return comments, basis_sets

def extract_basisinfo(filename, printing = False):
    swtch = 0
    f = open(filename)
    
    comments = []
    basis = []
    basisW = []
    basisE = []
    basisT = []
    contractedtype = []
    activecontracted = -1
    for i in f:
        I = i.split()
        if len(I)==0:
            continue
        if I[0] == "*":
            continue #ignore wildcards
        if I[0] == '#':        
            comments.append("//")
            for e in I:
                comments[-1] += e + " "
            continue
        #print I
        
        if I[0][0] == "$":
            swtch += 1
            swtch = swtch %2
            #if swtch == 1:
            #Add new function   
            continue
        if I[0] in ["h", "he", "li", "be", "b", "c", "n", "o", "f", "ne", "na" ,"mg", "al", "si", "p", "s", "cl", "ar", "k", "ca", "sc", "ti", "v", "cr", "mn", "fe", "co", "ni", "cu", "zn", "ga", "ge", "as", "se","br", "kr", "rb", "sr", "y", "zr", "nb", "mo", "tc", "ru", "rh", "pd", "ag", "cd", "in", "sn", "sb", "te", "i", "xe"]:
            bT = I[1]+"_"+I[0]
            bT = bT.replace("-", "_")
            basisT.append(bT)
            print basisT[-1]
            continue
        if I[0] in "123456789":
            if I[1] == "s":
                #print "Identified an s orbital."
                contractedtype.append(0)
            if I[1] == "p":
                #print "Identified an p orbital."
                contractedtype.append(1)
            if I[1] == "d":
                #print "Identified an d orbital."
                contractedtype.append(2)
            if I[1] == "f":
                #print "Identified an f orbital."
                contractedtype.append(3) 
            #Create contracted
            basis.append([])
            basisE.append([])
            basisW.append([])
            activecontracted += 1   
            #print I
        else:
            try:
                w = float(I[0].replace("D", "E"))
                e = float(I[1].replace("D", "E"))
                basisE[activecontracted].append(e)
                basisW[activecontracted].append(w)
            except:
                comments.append("//Ignored the following information from file:")
                for e in I:
                    comments[-1] += e + " "
                #comments[-1] += "Ignored the following information:" + I
                #print "Failed to read weigths and exponents."
                #print I
        #print I          
    if printing:
        print comments
        print contractedtype
        print len(contractedtype)
        print basisE
        print basisW
    f.close()
    return comments, contractedtype, len(contractedtype), basisE, basisW



def createfunction(fname, comments, orbital_types, N_orbitals, exponents, weights):
    endline = ["\n"]
    cppclass  = []
    cppheader = []
    for i in comments:
        cppclass.append(i)
        #cppheader.append(i)        
        
    #cppclass.append(endline)
    #cppheader.append(endline)
    
    #writing to header
    cppheader.append("    void add_"+fname+"(vec3 corePos);")
    

    cppclass.append("void basisbank::add_"+fname+"(vec3 corePos){\n")
    for i in range(N_orbitals):
        if orbital_types[i] == 0:
            #create an s-orbital
            cppclass.append("    bs.add_state();")
            for e in range(len(exponents[i])):
                cppclass.append("    Primitive S%iA%i = bs.turbomolePrimitive(%.8f,%.8f,0,0,0,corePos);" % (e, i,exponents[i][e], weights[i][e]))
                cppclass.append("    bs.add_primitive_to_state(bs.Nstates-1, S%iA%i);" % (e,i))
        if orbital_types[i] == 1:
            #create an p-orbital
            cppclass.append("    bs.add_state();")
            for e in range(len(exponents[i])):
                cppclass.append("    Primitive P%iA%i = bs.turbomolePrimitive(%.8f,%.8f,1,0,0,corePos);" % (e,i, exponents[i][e], weights[i][e]))
                cppclass.append("    bs.add_primitive_to_state(bs.Nstates-1, P%iA%i);" %( e,i))
                cppclass.append("    Primitive P%iB%i = bs.turbomolePrimitive(%.8f,%.8f,0,1,0,corePos);" % (e,i, exponents[i][e], weights[i][e]))
                cppclass.append("    bs.add_primitive_to_state(bs.Nstates-1, P%iB%i);" % (e,i))
                cppclass.append("    Primitive P%iC%i = bs.turbomolePrimitive(%.8f,%.8f,0,0,1,corePos);" % (e,i, exponents[i][e], weights[i][e]))
                cppclass.append("    bs.add_primitive_to_state(bs.Nstates-1, P%iC%i);" % (e,i))
                
        if orbital_types[i] == 2:
            #create an d-orbital
            cppclass.append("    bs.add_state();")
            for e in range(len(exponents[i])):
                cppclass.append("    Primitive D%iA%i = bs.turbomolePrimitive(%.8f,%.8f,2,0,0,corePos);" % (e, i,exponents[i][e], weights[i][e]))
                cppclass.append("    bs.add_primitive_to_state(bs.Nstates-1, D%iA%i);" % (e,i))
                cppclass.append("    Primitive D%iB%i = bs.turbomolePrimitive(%.8f,%.8f,0,2,0,corePos);" % (e,i, exponents[i][e], weights[i][e]))
                cppclass.append("    bs.add_primitive_to_state(bs.Nstates-1, D%iB%i);" %(e,i))
                cppclass.append("    Primitive D%iC%i = bs.turbomolePrimitive(%.8f,%.8f,0,0,2,corePos);" % (e,i, exponents[i][e], weights[i][e]))
                cppclass.append("    bs.add_primitive_to_state(bs.Nstates-1, D%iC%i);" %( e,i))
                
                cppclass.append("    Primitive D%iD%i = bs.turbomolePrimitive(%.8f,%.8f,1,1,0,corePos);" % (e, i,exponents[i][e], weights[i][e]))
                cppclass.append("    bs.add_primitive_to_state(bs.Nstates-1, D%iD%i);" % (e,i))
                cppclass.append("    Primitive D%iE%i = bs.turbomolePrimitive(%.8f,%.8f,0,1,1,corePos);" % (e,i, exponents[i][e], weights[i][e]))
                cppclass.append("    bs.add_primitive_to_state(bs.Nstates-1, D%iE%i);" % (e,i))
                cppclass.append("    Primitive D%iF%i = bs.turbomolePrimitive(%.8f,%.8f,1,0,1,corePos);" % (e, i,exponents[i][e], weights[i][e]))
                cppclass.append("    bs.add_primitive_to_state(bs.Nstates-1, D%iF%i);" %( e,i))
        if orbital_types[i] == 3:
            #create an d-orbital
            cppclass.append("    bs.add_state();")
            for e in range(len(exponents[i])):
                cppclass.append("    Primitive E%iA%i = bs.turbomolePrimitive(%.8f,%.8f,3,0,0,corePos);" % (e, i,exponents[i][e], weights[i][e]))
                cppclass.append("    bs.add_primitive_to_state(bs.Nstates-1, E%iA%i);" % (e,i))
                cppclass.append("    Primitive E%iB%i = bs.turbomolePrimitive(%.8f,%.8f,0,3,0,corePos);" % (e,i, exponents[i][e], weights[i][e]))
                cppclass.append("    bs.add_primitive_to_state(bs.Nstates-1, E%iB%i);" %(e,i))
                cppclass.append("    Primitive E%iC%i = bs.turbomolePrimitive(%.8f,%.8f,0,0,3,corePos);" % (e,i, exponents[i][e], weights[i][e]))
                cppclass.append("    bs.add_primitive_to_state(bs.Nstates-1, E%iC%i);" %( e,i))
                
                cppclass.append("    Primitive E%iD%i = bs.turbomolePrimitive(%.8f,%.8f,1,2,0,corePos);" % (e, i,exponents[i][e], weights[i][e]))
                cppclass.append("    bs.add_primitive_to_state(bs.Nstates-1, E%iD%i);" % (e,i))
                cppclass.append("    Primitive E%iE%i = bs.turbomolePrimitive(%.8f,%.8f,0,1,2,corePos);" % (e,i, exponents[i][e], weights[i][e]))
                cppclass.append("    bs.add_primitive_to_state(bs.Nstates-1, E%iE%i);" % (e,i))
                cppclass.append("    Primitive E%iF%i = bs.turbomolePrimitive(%.8f,%.8f,1,0,2,corePos);" % (e, i,exponents[i][e], weights[i][e]))
                cppclass.append("    bs.add_primitive_to_state(bs.Nstates-1, E%iF%i);" %( e,i))
                
                cppclass.append("    Primitive E%iG%i = bs.turbomolePrimitive(%.8f,%.8f,2,1,0,corePos);" % (e, i,exponents[i][e], weights[i][e]))
                cppclass.append("    bs.add_primitive_to_state(bs.Nstates-1, E%iG%i);" % (e,i))
                cppclass.append("    Primitive E%iH%i = bs.turbomolePrimitive(%.8f,%.8f,0,2,1,corePos);" % (e,i, exponents[i][e], weights[i][e]))
                cppclass.append("    bs.add_primitive_to_state(bs.Nstates-1, E%iH%i);" % (e,i))
                cppclass.append("    Primitive E%iI%i = bs.turbomolePrimitive(%.8f,%.8f,2,0,1,corePos);" % (e, i,exponents[i][e], weights[i][e]))
                cppclass.append("    bs.add_primitive_to_state(bs.Nstates-1, E%iI%i);" %( e,i))               
                
                cppclass.append("    Primitive E%iJ%i = bs.turbomolePrimitive(%.8f,%.8f,1,1,1,corePos);" % (e, i,exponents[i][e], weights[i][e]))
                cppclass.append("    bs.add_primitive_to_state(bs.Nstates-1, E%iJ%i);" %( e,i))   
                
    cppclass.append("}\n")
    
    return cppheader, cppclass

def saveclass(cppclasses, cppheaders):
    cppfile = "//This file is maintained by an external python script and should not be edited manually.\n#include <basisbank.h>\n#include <armadillo>\n#include <basis.h>\n#include <primitive.h>\nusing namespace std;\nusing namespace arma;\nbasisbank::basisbank(basis BS){\n    bs = BS;} \nbasisbank::basisbank(){} \n \n"
    
    for i in cppclasses:
        for e in i:
            cppfile += e+ "\n"
        cppfile += "\n"
    #print cppfile
    
    
    cppheader = "//This file is maintained by an external python script and should not be edited manually.\n#ifndef BASISBANK_H\n#define BASISBANK_H\n#include <armadillo>\n#include <basis.h>\n#include <primitive.h>\nusing namespace std;\nusing namespace arma;\n \nclass basisbank{\npublic:\n    basisbank(basis BS);\n    basisbank();\n    basis bs;\n    string basistype;"
    for i in cppheaders:
        for e in i:
            cppheader += e + "\n"
    cppheader += "};\n"
    cppheader += "#endif // BASISBANK_H"
    #print cppheader
    #Write files
    cpp_header = open("basisbank.h", "w")
    cpp_class = open("basisbank.cpp", "w")
    cpp_header.write(cppheader)
    cpp_class.write(cppfile)
    cpp_header.close()
    cpp_class.close()

    
def convertfolder():
    
    cppclasses = []
    cppheaders = []
    folder  = [] #a list containing all files in folder
    #os.chdir("/")
    #print os.getcwd()
    for filename in os.listdir(os.getcwd()):
        if filename.endswith(".txt"):
            try:
                #comments, orbital_types, N_orbitals, exponents, weigths = extract_basisinfo(filename)
                comments, I = extract_m_basisinfo(filename)
                commented = 0
                #print comments

                
                for i in I:
                    #i = [bT, basisW, basisE, contractedtype, len(contractedtype)]
                    if commented == 0:
                        cppheader, cppclass = createfunction(i[0], comments, i[3], i[4], i[2], i[1])
                        cppclasses.append(cppclass)
                        cppheaders.append(cppheader)
                        commented = 1
                    else:
                        cppheader, cppclass = createfunction(i[0], [""], i[3], i[4], i[2], i[1])
                        cppclasses.append(cppclass)
                        cppheaders.append(cppheader)

                        
            except:
                print "Failed to import file:", filename

    saveclass(cppclasses, cppheaders)

#filename = "STO6G_H.txt"
#comments, orbital_types, N_orbitals, exponents, weigths = extract_basisinfo(filename)
#cppheader, cppclass = createfunction(filename, comments, orbital_types, N_orbitals, exponents, weigths)
convertfolder() #converts all .txt files in current folder
#saveclass([cppclass],[cppheader])
#filename = "sto3g_all.txt"
#comm, info = extract_m_basisinfo(filename)
#for i in info:
#    print i[0]
#print exponents

convertfolder() #converts all .txt files in current folder
    
    
    
    