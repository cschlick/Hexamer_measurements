import sys, os

# Edit this path to the location of the PDBModule folder on your computer
sys.path.append("../")

from PDBFile import PDBFile

def printhelp():
    helpmessage= """
pdbRefactor -help

Usage: pdbRefactor input.pdb type=TYPE  extract=SECTION
Where 'TYPE' is the output type 'phenix' or 'chimera' or 'coot'
And 'SECTION' is optional, can be 'asu' or 'patch'

Also required is a parameter file located either in the working directory, or one or two directories up.
parameter file should be exactly 'pdbparams.txt' and should include three parameters:

monomerlen = integer    # The length of protein residues in a monomer.
asusize = integer       # The number of monomers in an asu. Note that for this program, all monomers MUST be an identical length
numligands = integer    # The number of non-protein ligands per asymmetric unit. (specified with LIG as the residue ID)
endsequence = string    # Several letters of protein sequence to identify the end of the monomer chain. ie: PNAPILST
patchsize = integer     # The number of nearby asymmetric units to extract when using 'extract=patch'. You will probably need to adjust this. For HBV T=4, patchsize=4 extracts a 6-fold.

These parameters sare used to divide up a pdb into symmetric sections, assuming the order they occur in the pdb is Asu1, Asu2, Asu3, etc.
Different orderings will cause this to fail. For example, putting all the ligands together at the end of the file. To generate a sample parameter file run pdbRefactor -initialize in the directory you want the file.
    """
    print(helpmessage)


def argprocess(arg):
    split = arg.split("=")
    for i, item in enumerate(split):
        split[i] = item.strip(' \n\t')
    return split[-1]

args = sys.argv



for arg in args:
    if (arg == "-h") or (arg == "--h") or (arg == "-H") or (arg=="--H") or (arg == "-help") or (arg == "--help") or (arg == "-Help") or (arg=="--Help"):
        printhelp()
        sys.exit(0)

    if (arg == "-initialize") or (arg == "--initialize"):
        file = open("pdbparams.txt","w")
        samplestring = """
monomerlen=142
asusize=4
numligands=2
endsequence=PNAPILST
patchsize=4

"""
        file.write(samplestring)
        file.close()
        print("\nWrote sample parameter file.\n")
        sys.exit(0)
try:
    type = None
    extract = None

    filename = args[1]
    for arg in args:
        if "type=" in arg:
            type = argprocess(arg)
        elif "extract=" in arg:
            extract = argprocess(arg)
    if (type != None) and (extract!=None):
        print("Refactoring pdb to type="+type+" extract="+str(extract))
    elif  (type != None):
        print("Refactoring pdb to type="+type)
    else:
        print("\nNeed to provide at least an output type.\n")
        sys.exit(0)

except:

    print("\nFailed to process input\n")
    printhelp()
    sys.exit(0)




########################################################################
##PARSE
pdb = PDBFile(filepath=filename)
pdb.parametersParse()
############################################################################

def remove_model_tags(filepath):
    #remove model references
    print("Removing MODEL references for Coot....")
    file = open(outputfile,"r")
    lines = file.readlines()
    newlines = []
    for line in lines:
        if ("MODEL" not in line) and ("ENDMDL" not in line) :
            newlines.append(line)
        else:
            print("Removing MODEL entry "+line)
    file.close()
    os.remove(outputfile)
    file = open(outputfile,"w")
    for line in newlines:
        file.write(line)
    file.close()


if (type == 'phenix') or (type == 'Phenix') or (type == 'PHENIX'):
    pdb.convert_to_phenix()

    if (extract != None):
        if (extract == "asu") or (extract == "ASU") or (extract == "Asu"):
            print("Extracting a single asu....")
            asu= pdb.extract_asu()
            pdb.clear()
            for residue in asu.residues:
                for atom in residue.atoms:
                    pdb.atoms.append(atom)
            outputfile = pdb.filepath.replace(".pdb", "_extractedASU_PHENIX.pdb")
            pdb.write(outputfile)

        if (extract == "patch") or (extract == "6fold") or (extract == "Patch"):
            print("Extracting approximately a 6fold patch...")
            patch_groups = pdb.extract_n_asus(pdb.patchsize)
            pdb.clear()
            for group in patch_groups:
                for residue in group.residues:
                    for atom in residue.atoms:
                        pdb.atoms.append(atom)
            outputfile = pdb.filepath.replace(".pdb", "_extracted6fold_PHENIX.pdb")
            pdb.write(outputfile)
    else:
        outputfile = pdb.filepath.replace(".pdb", "_PHENIX.pdb")
        pdb.write(outputfile)





elif (type == 'chimera') or (type == 'Chimera') or (type == 'CHIMERA'):
    pdb.convert_to_chimera()
    if (extract != None):
        if (extract == "asu") or (extract == "ASU") or (extract == "Asu"):
            print("Extracting a single asu....")
            asu= pdb.extract_asu()
            pdb.clear()
            for residue in asu.residues:
                for atom in residue.atoms:
                    pdb.atoms.append(atom)
            outputfile = pdb.filepath.replace(".pdb", "_extractedASU_CHIMERA.pdb")
            pdb.write(outputfile)
        if (extract == "patch") or (extract == "6fold") or (extract == "Patch"):
            print("Extracting approximately a 6fold patch...")
            patch_groups = pdb.extract_n_asus(pdb.patchsize-1)
            pdb.clear()
            for group in patch_groups:
                for residue in group.residues:
                    for atom in residue.atoms:
                        pdb.atoms.append(atom)
            outputfile = pdb.filepath.replace(".pdb", "_extracted6fold_CHIMERA.pdb")
            pdb.write(outputfile)

    else:
        print("Trying to make outputfile")
        outputfile = pdb.filepath.replace(".pdb", "_CHIMERA.pdb")
        pdb.write(outputfile)

elif (type == 'coot') or (type == 'Coot') or (type == 'COOT'):
    pdb.convert_to_chimera()

    if (extract != None):
        if (extract == "asu") or (extract == "ASU") or (extract == "Asu"):
            print("Extracting a single asu....")
            asu= pdb.extract_asu()
            pdb.clear()
            for residue in asu.residues:
                for atom in residue.atoms:
                    pdb.atoms.append(atom)
            outputfile = pdb.filepath.replace(".pdb", "_extractedASU_COOT.pdb")
            pdb.write(outputfile)
            remove_model_tags(outputfile)

        if (extract == "patch") or (extract == "6fold") or (extract == "Patch"):
            print("Extracting approximately a 6fold patch...")
            patch_groups = pdb.extract_n_asus(pdb.patchsize-1)
            pdb.clear()
            for group in patch_groups:
                for residue in group.residues:
                    for atom in residue.atoms:
                        pdb.atoms.append(atom)
            outputfile = pdb.filepath.replace(".pdb", "_extracted6fold_COOT.pdb")
            pdb.write(outputfile)
            remove_model_tags(outputfile)




    else:
        outputfile = pdb.filepath.replace(".pdb", "_COOT.pdb")
        pdb.write(outputfile)
        print("Writing to "+outputfile)
        remove_model_tags(outputfile)



else:
    print("The supplied 'type=' is not recognized.")
    printhelp()
