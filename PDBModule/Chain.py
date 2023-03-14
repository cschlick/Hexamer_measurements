class Chain():
    def __init__(self,chain_id):
        self.type= "Chain"
        self.chain_id = chain_id
        self.residues = []
        #self.update()


    def update(self):
        for residue in self.residues:
            residue.chain_id = self.chain_id
            residue.update()





    def summary(self):
        print("################Summary Chain ID: ",self.chain_id,":\tContains ",len(self.residues)," residues.")
        print("\n")



    def print(self):
        print("Chain ID: ",self.chain_id,"\t\t\t...printing top 10 residues.")
        print("\n")
        for residue in self.residues[0:10]:
            residue.print()