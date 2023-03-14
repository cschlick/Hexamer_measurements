from Atom import Atom
import pathlib


class PDBIO():
    def __init__(self,filepath): # takes a filepath and returns a pdbfile object
        self.filepath = filepath



    def parse(self):
        print("Parsing File. Will ignore all metadata.")
        file = open(self.filepath,"r")
        lines = file.readlines()
        file.close()
        for line in lines:
            if line[0:4] == "ATOM":
                atom = self.parse_atom_line(line)
                self.atoms.append(atom)

            if line[0:6] == "HETATM": # treat just like atom for now
                atom = self.parse_atom_line(line)
                self.atoms.append(atom)
        print("Finished parsing atoms, found ",len(self.atoms)," atoms.")
        print("\n")




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
        return atom