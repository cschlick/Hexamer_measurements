

class ResidueGroup():
    def __init__(self,number,residues): # a class that assumes the atoms in a pdb are ordered with an new asu every n:n+c residues
        self.type = "Group"
        self.number = number
        self.name = number
        self.residues = residues




    def print(self):
        print("Group ID: ",self.number,"\t\t\t...printing top 10 residues.")
        print("\n")
        for residue in self.residues[0:10]:
            residue.print()