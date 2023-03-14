import os



class PDBWriter():
    def __init__(self,outputpath,data):
        self.outputpath = outputpath
        self.data = data
        self.datatype = self.data[0].type


    def write(self): #change to chains,residues, etc
        print("Writing data")
        outputpath = self.outputpath

        if os.path.isfile(outputpath):
            os.remove(outputpath)
        file = open(outputpath,"a")
        print("Data type is ",self.datatype)
        if self.datatype == "ASU":
            for asu in self.data:
                for residue in asu.residues:
                    for atom in residue.atoms:
                        file.write(atom.generate_entry())

        if self.datatype == "Residue":
            for residue in self.data:
                for atom in residue.atoms:
                    file.write(atom.generate_entry())




        file.close()
        print("Written atoms to ",outputpath)

