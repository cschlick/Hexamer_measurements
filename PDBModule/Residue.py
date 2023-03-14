class Residue():
    def __init__(self,residue_name,atoms):
        #######FLAGS
        self.__setup = True



        self.name = residue_name
        self.residue_name = self.name
        self.atoms = atoms
        self.residue_number = atoms[0].residue_number
        self.first_atom_number = atoms[0].serial_number
        self.chain_id = atoms[0].chain_id

        self.__setup = False
################GETTER AND SETTERS############################
    @property
    def residue_name(self):
        return self.__residue_name

    @residue_name.setter
    def residue_name(self,residue_name):
        if self.__setup is False:
            #propagate changes down to atoms
            for atom in self.atoms:
                atom.residue_name = residue_name
        self.__residue_name = residue_name

    @property
    def chain_id(self):
        return self.__chain_id

    @chain_id.setter
    def chain_id(self,chain_id):
        if self.__setup is False:
            #propagate changes down to atoms
            for atom in self.atoms:
                atom.chain_id = chain_id
        self.__chain_id = chain_id


    @property
    def residue_number(self):
        return self.__residue_number

    @residue_number.setter
    def residue_number(self,residue_number):
        if self.__setup is False:
            #propagate changes down to atoms
            for atom in self.atoms:
                atom.residue_number = residue_number
        self.__residue_number = residue_number


################END GETTER AND SETTERS############################



    def print(self):

        for atom in self.atoms:
            atom.print()
        print("\n")

    def add_atom(self,atom):
        self.atoms.append(atom)