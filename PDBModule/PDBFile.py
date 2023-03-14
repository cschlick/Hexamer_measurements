import os, sys
from Atom import Atom
from Residue import Residue
from Chain import Chain
from ResidueGroup import ResidueGroup
import numpy as np
import random
from pathlib import Path
from multiprocessing import Pool

class PDBFile():
    def __init__(self,filepath=None,debug=False,build_residues=True):
        self.filepath = filepath
        self.atoms = []
        self.hetatoms = []
        self.models = []
        self.residue_groups = []
        self.ligand_groups = []
        self.monomerlen = None
        self.asulen = None
        self.splitmodels = False
        if self.filepath != None:
            print("\nStart PDB Parsing")
            if os.path.isfile(self.filepath):
                self.filepath = filepath
                print("*************************************************************************")
                print("Found File ",self.filepath)
            elif os.path.isfile(str(os.getcwd())+self.filepath):
                self.filepath= str(os.getcwd())+self.filepath
                print("\n*************************************************************************")
                print("Found File ",str(os.getcwd())+self.filepath)
            else:
                print("Failed. Not valid filepath.")
            self.parse()
        ##############################
        if build_residues == True:
            self.residues = self.build_residues()
            self.chains = self.build_chains()

        if debug == True:
            self.residue_num_distribution()


        #self.trim_residue(143)
        print("End PDB Parsing")
        print("***********************************************************************\n")

    def rebuild(self): # if you added or removed atoms, do this
        self.residues = self.build_residues()
        self.chains = self.build_chains()


    def parse_atom_line(self,line): # parses an atom line, returns an atom object
        record_name = line[:6]
        serial_number = line[6:11]
        atom_name = line[11:17]
        alt_loc =line[16:17]
        residue_name = line[17:20]
        chain_id = line[20:22]
        residue_number = line[22:26]
        icode = line[26:30]
        x_coord = line[30:38]
        y_coord = line[38:46]
        z_coord = line[46:54]
        occupancy = line[54:60]
        bfactor = line[60:66]
        element = line[66:78]
        charge = line[78:]
        #troubleshooting
        #print(record_name+","+serial_number+","+atom_name+","+alt_loc+","+res_name+","+chain_id+","+res_seq+","+icode+","+xcoord+","+ycoord+","+zcoord+","+occupancy+","+bfactor+","+element+","+charge+",")
        atom = Atom(record_name,serial_number,atom_name,alt_loc,residue_name,chain_id,residue_number,icode,x_coord,y_coord,z_coord,occupancy,bfactor,element,charge)

        #print(record_name,serial_number,atom_name,alt_loc,residue_name,chain_id,residue_number,icode,x_coord,y_coord,z_coord,occupancy,bfactor,element,charge)

        # no hydrogens
        if "H" in atom.atom_name:
             return None
        else:
            return atom

    def parse(self):
        print("Parsing File. Will ignore all metadata.")
        file = open(self.filepath,"r")
        lines = file.readlines()
        file.close()

        for line in lines:
            if line[0:4] == "ATOM":
                atom = self.parse_atom_line(line)
                if atom is not None:
                    self.atoms.append(atom)

            if line[0:6] == "HETATM": # treat just like atom for now
                atom = self.parse_atom_line(line)
                self.atoms.append(atom)
        print("Finished parsing atoms, found ",len(self.atoms)," atoms.")
        print("\n")


    def write(self,outputpath,method="Atoms"): #change to chains,residues, etc

        if os.path.isfile(outputpath):
            os.remove(outputpath)
        file = open(outputpath,"a")
        if self.splitmodels == True:
            file.write("MODEL\tMDL_1\n")
        for atom in self.atoms:
            file.write(atom.generate_entry())
            if atom.term == True:
                file.write("TER\n")
            if atom.modelsplit == True:
                file.write("ENDMDL\n")
                name = "MDL_"+str(hex(int(str(random.randint(0,10000000))))[2:])
                file.write("MODEL\t"+name+"\n")

        file.close()
        print("Written atoms to ",outputpath)




    def build_residues(self):
        atomlist = self.atoms
        residues = []
        current_res_num = int(atomlist[0].residue_number)
        current_residue_name = atomlist[0].residue_name
        current_chain_id = atomlist[0].chain_id
        starti = 0
        endi = None
        for i in range(0,len(atomlist)):
            atom = atomlist[i]

            if (int(atom.residue_number) == current_res_num) and (atom.chain_id == current_chain_id):
                if atom.residue_name == current_residue_name:
                    #print("Residue name mismatch")
                    #print("ith res seq num ",atom.residue_number," Current num ",current_res_num)
                    #print(atom.residue_name,atom.serial_number," not matching current name:",current_residue_name)
                    pass
            else:
                endi = i
                resname = atomlist[i-1].residue_name
                residue = Residue(resname,atomlist[starti:endi])
                residues.append(residue)
                starti = i
                current_residue_name = atom.residue_name
                current_res_num= int(atom.residue_number)
                current_chain_id = atom.chain_id
        #handle last one not being triggered
        endi = i+1
        resname = atomlist[-1].residue_name
        residue = Residue(resname,atomlist[starti:endi])
        residues.append(residue)
        l = len(residues)

        print("Found ",l," residues. If icosahedral, that is ",l/60," residues per asu.")
        #print("Per monomer T3: ",round(l/60/3,2)," residues. With 1 drug, ",round((l-60)/60/3,2)," residues.")
        #print("Per monomer T4: ",round(l/60/4,2)," residues. With 1 drug, ",round((l-60)/60/4,2)," residues.")
        print("\n")
        return residues


    def build_chains(self):
        chains = {}
        for residue in self.residues:
            if residue.chain_id not in chains:
                chains[residue.chain_id] = Chain(residue.chain_id)

        for residue in self.residues:
            chains[residue.chain_id].residues.append(residue)

        chainlist = []
        for chain in chains.values():
            chainlist.append(chain.chain_id)
        print("Found ",len(chains.values())," chains: ",chainlist)
        print("\n")
        return chains.values()


    def res_name(self):
        reslist = []
        for atom in self.atoms:
            if atom.residue_name not in reslist:
                reslist.append(atom.residue_name)
        return reslist

    def build_residue_group(self,expected_length=None,debug=False):
        if expected_length == None:
            print("Provide expected group length. Trying with 142")
            expected_length = 142
        print("Attempting to build Residue  Groups.....")
        groups = []
        l = len(self.residues)
        if l % expected_length != 0:
            print("FAILED!")
            print("Number of residues ",l," not divisible by expected length ",expected_length,". Consider triming or removing ligands.")
            self.residue_num_distribution()
        else:
            print("Making ",l/expected_length," groups....")
            for i in range(0,int((l/expected_length))):

                n = int(i*expected_length)

                group = ResidueGroup(i,self.residues[n:n+expected_length])
                if debug == True:
                    num1 = group.residues[0].residue_number
                    num2 = group.residues[-1].residue_number
                    name1 = group.residues[0].residue_name
                    name2 = group.residues[-1].residue_name
                    print("Making groups with residues ",name1,num1," to ",name2,num2)

                groups.append(group)
            print("Success...probably. Made ",len(groups)," new groups.")
        self.residue_groups = groups

    def endsequence_conversion(self,endsequence): # takes an endsequence of eith upper case single letters, lower case single letters, or three leter codes and returns three letter codes
        uppercasedict = {"R": "ARG","H": "HIS","K": "LYS","D": "ASP","E": "GLU","S":"SER","T":"THR","N":"ASN","Q":"GLN","C":"CYS","G":"GLY","P":"PRO","A":"ALA","V":"VAL","I":"ILE","L":"LEU","M":"MET","F":"PHE","Y":"TYR","W":"TRP"}
        lowercasedict = {"r": "ARG","h": "HIS","k": "LYS","d": "ASP","e": "GLU","s":"SER","t":"THR","n":"ASN","q":"GLN","c":"CYS","g":"GLY","p":"PRO","a":"ALA","v":"VAL","i":"ILE","l":"LEU","m":"MET","f":"PHE","y":"TYR","w":"TRP"}
        if (len(endsequence[0]) != 1) and (len(endsequence[0]) != 3):
            print("There is a problem with the endsequence provided. ")
        elif(len(endsequence[0]) == 3):
            print("Assuming endsequence is three letter codes, not converting")
        elif endsequence[0].islower():
            newsequence = []
            for aa in endsequence:
                newsequence.append(lowercasedict[aa])
            return newsequence
        elif endsequence[0].isupper():
            newsequence = []
            for aa in endsequence:
                newsequence.append(uppercasedict[aa])
            return newsequence

        else:
            print("Problem encountered converting single letter sequence to 3 letter sequence")



    def build_residue_group2(self,endsequence):
        print("Attempting to build Residue groups...")

        endsequence = self.endsequence_conversion(endsequence)
        if self.monomerlen == None or self.asulen == None:
            print("Need to set pdb.monomerlen and pdb.asulen first")
            sys.exit(0)
        self.residue_groups = []
        end_seq_len = len(endsequence)
        sequence = []
        i = 1 # residue counter
        for residue in self.residues:
            sequence.append(residue.residue_name)
            if (len(sequence)>end_seq_len):
                if sequence[-end_seq_len:] == endsequence:
                    start = i-self.monomerlen
                    end = i
                    residues_to_add = self.residues[start:end]
                    #print("Adding "+str(len(residues_to_add))+" residues to group, from "+str(start)+" to "+str(end))
                    group = ResidueGroup(i,residues_to_add)
                    self.residue_groups.append(group)
            i +=1
        print("Success...probably. Made ",len(self.residue_groups)," new groups.")


    def build_ligand_groups(self,numligs):
        self.ligand_groups = []
        ligands = []
        for residue in self.residues:
            if residue.residue_name == "LIG":
                ligands.append(residue)
        print("Found "+str(len(ligands))+" ligands.")
        i = 1 # residue counter
        n = 1 # another counter
        for ligand in ligands:
            if n == self.numligands:
                start = i-numligs
                end = i
                ligs_to_add = ligands[start:end]
                #print("Adding "+str(len(ligs_to_add))+" residues to group, from "+str(start)+" to "+str(end))
                group = ResidueGroup(i,ligs_to_add)
                self.ligand_groups.append(group)
                n=1
            else:
                n+=1
            i+=1




    def distance_atoms(self,atom1,atom2):
        dist = np.linalg.norm(atom1.get_coord()-atom2.get_coord())
        return dist

    def distance_points(self,xyz1,xyz2):
        dist = np.linalg.norm(xyz1-xyz2)
        return dist

    def merge_groups(self,number,debug=False): # this will merge each n successive asus into a single asu
        new_groups = []
        n = len(self.residue_groups)/number
        print("Will try to make ",n," groups out of ",len(self.residue_groups))
        for i in range(0,int(n)):
            residues = []
            n = int(i*number)
            group_subset = self.residue_groups[n:n+number]
            for group in group_subset:
                for residue in group.residues:
                    residues.append(residue)

            group = ResidueGroup(i,residues)
            new_groups.append(group)
        print("Resorted into ",len(new_groups),"  groups")
        self.residue_groups= new_groups




    def residue_num_distribution(self):
        print("Assessing Distribution of Residue numbers:")
        reskeys = []
        resdict = {}
        resnums = []
        for residue in self.residues:
            num = residue.residue_number
            resnums.append(num)
            if num not in reskeys:
                reskeys.append(num)
        for key in reskeys:
            resdict[key] = resnums.count(key)
        for key in reskeys:
            print("Residue Number: ",key,"\tOccurance number: ",resdict[key])

    def atom_num_distribution(self):
        print("Assessing Distribution of Atom numbers:")
        reskeys = []
        resdict = {}
        resnums = []
        for atom in self.atoms:
            num = atom.serial_number
            resnums.append(num)
            if num not in reskeys:
                reskeys.append(num)
        for key in reskeys:
            resdict[key] = resnums.count(key)
        for key in reskeys:
            print("Atom Number: ",key,"\tOccurance number: ",resdict[key])



    def trim_residue(self,residue_number): # for now use ints
        number = int(residue_number)
        new_res = []
        for residue in self.residues:
            if int(residue.residue_number) != number:
                new_res.append(residue)

        self.residues = new_res
        new_atoms = []
        for atom in self.atoms:
            if int(atom.residue_number) != number:
                new_atoms.append(atom)

        self.atoms= new_atoms

        self.chains = self.build_chains()


    def parametersParse(self): # look for a parameter file in the current directory, and two directories up
        cwd = Path(os.getcwd())
        if os.path.isfile(str(cwd)+"/pdbparams.txt"):
            print("Found parameter file in working directory.")
            parameterfilename = str(cwd)+"/pdbparams.txt"
        elif os.path.isfile(str(cwd.parent)+"/pdbparams.txt"):
            print("Found parameter file in parent directory")
            parameterfilename = str(cwd.parent)+"/pdbparams.txt"
        elif os.path.isfile(str(cwd.parent.parent)+"/pdbparams.txt"):
            print("Found parameter file two directories up")
            parameterfilename = str(cwd.parent.parent)+"/pdbparams.txt"
        else:
            print("Cound not find 'pdbparams.txt'. Please place it in either current directoy, parent directory, or two directories up.\nUse pdbRefactor -h for more info.")
            sys.exit(0)
        parameterfile = open(parameterfilename,"r")
        lines = parameterfile.readlines()
        parameterfile.close()
        for line in lines:
            if "monomerlen" in line:
                split = line.split("=")
                for i, item in enumerate(split):
                    split[i] = item.strip(' \n\t')
                self.monomerlen = int(split[-1])
        for line in lines:
            if "asusize" in line:
                split = line.split("=")
                for i, item in enumerate(split):
                    split[i] = item.strip(' \n\t')
                self.asusize = int(split[-1])
        for line in lines:
            if "numligands" in line:
                split = line.split("=")
                for i, item in enumerate(split):
                    split[i] = item.strip(' \n\t')
                self.numligands = int(split[-1])

        for line in lines:
            if "endsequence" in line:
                split = line.split("=")
                for i, item in enumerate(split):
                    split[i] = item.strip(' \n\t')
                seq = list(split[-1])
                self.endsequence = seq

        for line in lines:
            if "patchsize" in line:
                split = line.split("=")
                for i, item in enumerate(split):
                    split[i] = item.strip(' \n\t')
                self.patchsize = int(split[-1])
        try:
            self.asulen = (self.monomerlen*self.asusize)+self.numligands
            print("\nParameter parsing results:")
            print("Residue length of monomer: "+str(self.monomerlen))
            print("Number of ligand residues per asu: "+str(self.numligands))
            print("Number of monomers per asu: "+str(self.asusize))
            print("Total calculated number of residues per asu: "+str(self.asulen))
            print("Monomer end sequence: "+str(self.endsequence))
            print("Number of asymmetric units to extract with 'extract=patch': "+str(self.patchsize)+"\n")
        except:
            print("Failed to completely parse parameter file. Check that it is complete. Use 'pdbRefactor -help' for instructions.")
            sys.exit(0)


###############################################################################################################scripting functions
    def renumber_residue_by_group(self):
        i = 1
        for group in self.residue_groups:
            for residue in group.residues:
                residue.residue_number=i
                i=i+1
            i =1

    def rename_chains_by_group(self):
        newchains = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T','U', 'V', 'W', 'X', 'Y', 'Z','a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't','u', 'v', 'w', 'x', 'y', 'z','1','2','3','4','5','6','7','8','9']
        i = 0
        for group in self.residue_groups:
            for residue in group.residues:
                residue.chain_id = newchains[i]
                
            i=i+1

    def rename_chains_by_group_identical_asus(self):
        newchains =['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T','U', 'V', 'W', 'X', 'Y', 'Z','a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't','u', 'v', 'w', 'x', 'y', 'z','1','2','3','4','5','6','7','8','9']
        i = 0
        for group in self.residue_groups:
            for residue in group.residues:
                residue.chain_id = newchains[i]
                
            if i == self.asusize-1:
                i = 0
            else:
                i=i+1


    def renumber_atoms_sequentially(self):
        i = 1
        for atom in self.atoms:
            atom.serial_number=i
            i+=1
    def renumber_atoms_by_group(self):
        i = 1
        for group in self.residue_groups:
            for residue in group.residues:
                for atom in residue.atoms:
                    atom.serial_number = i
                    i=i+1
            i = 1


    def convert_to_chimera(self):
            print("\nRefactoring pdb to Chimera style......\n")

            self.build_residue_group2(self.endsequence)
            self.build_ligand_groups(self.numligands)
            #self.rename_chains_by_group_identical_asus()
            self.renumber_residue_by_group()
            self.renumber_atoms_by_group()
            #rename ligand chain id
            for ligandgroup in self.ligand_groups:
                for residue in ligandgroup.residues:
                    residue.chain_id = "L"
                    
            # renumber ligand residues by chain
            i = 1
            for ligandgroup in self.ligand_groups:
                for residue in ligandgroup.residues:
                    residue.residue_number=i
                    
                    i=i+1
                i =1
            # renumber ligand atoms by chain
            i = 1
            for ligandgroup in self.ligand_groups:
                for residue in ligandgroup.residues:
                    for atom in residue.atoms:
                        atom.serial_number = i
                        i=i+1
                i =1

            ## Add term at end of each group
            for group in self.residue_groups:
                group.residues[-1].atoms[-1].term = True  #set term flag for last atom in last residue of each group
            #group into asus before output
            self.build_residue_group(expected_length=self.asulen)
            self.renumber_atoms_by_group()

            ## Add model at end of each group
            self.splitmodels = True
            for group in self.residue_groups:
                group.residues[-1].atoms[-1].modelsplit = True  #set term flag for last atom in last residue of each group
            print("Passed: convert_to_chimera()")



    def convert_to_phenix(self):
        print("\nRefactoring pdb to PHENIX style......\n")

        self.build_residue_group(expected_length=self.asulen)
        self.renumber_residue_by_group()
        self.rename_chains_by_group()
        self.renumber_atoms_sequentially()


    def clear(self):
        self.residue_groups = []
        self.residues = []
        self.atoms= []

    def com_group(self,group): #returns the center of mass of all atoms in the group
        x = []
        y = []
        z = []
        for residue in group.residues:
            for atom in residue.atoms:
                x.append(float(atom.x_coord))
                y.append(float(atom.y_coord))
                z.append(float(atom.z_coord))
        x = np.asarray(x)
        y = np.asarray(y)
        z = np.asarray(z)

        x_mean = np.mean(x)
        y_mean = np.mean(y)
        z_mean = np.mean(z)
        return np.array([x_mean,y_mean,z_mean])



    def extract_asu(self):
        self.build_residue_group(expected_length=self.asulen)
        return self.residue_groups[0]




    def extract_n_asus(self,n): # will extract n asus close together
        self.build_residue_group(expected_length=self.asulen)
        comdict = {}
        numbers = []
        for asu in self.residue_groups:
            com = self.com_group(asu)
            comdict[asu.number] = com
            numbers.append(asu.number)

        first_number = min(numbers)


        # #######starting search chain
        patch_groups = []
        for group in self.residue_groups:
            if group.number == first_number:
                patch_groups.append(group)

        dist_dict = {}
        com_1 = comdict[first_number]
        del comdict[first_number]

        current_com = com_1
        for iteration in range(n):
            dists = []
            for key,value in comdict.items():
                dist = self.distance_points(current_com,value)
                dist_dict[dist] = key
                dists.append(dist)

            newgroup_number = dist_dict[min(dists)]
            for group in self.residue_groups:
                if group.number == newgroup_number:
                    patch_groups.append(group)

            del comdict[newgroup_number]
            ####update current com to include all asus added so far
            current_residues = []
            for group in patch_groups:
                for residue in group.residues:
                    current_residues.append(residue)
            current_group = ResidueGroup(iteration,current_residues)
            current_com = self.com_group(current_group)

        return patch_groups
