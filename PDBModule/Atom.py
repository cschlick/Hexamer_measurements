import numpy as np

class Atom():
    def __init__(self,record_name,serial_number,atom_name,alt_loc,res_name,chain_id,res_seq,icode,x_coord,y_coord,z_coord,occupancy,bfactor,element,charge):


        ######### FLAGS AND SWITCHES ######################
        self.__writing_entry = False
        self.modelsplit = False # if set to true, will write endmodel after this residue and startmodel with a random string
        self.term = False # if set to true, will write a term after this residue when pdb writes out
        #skipping alt loc for now
        self.__alt_loc_exists = False

        #########PROPERTIES######################### (All properties are stored as unpadded strings)

        self.record_name = record_name
        self.serial_number = serial_number
        if self.__alt_loc_exists is True:   # handle the presence of alt loc
            self.atom_name = atom_name[:-2]
            self.alt_loc = alt_loc
        else:
            self.atom_name = atom_name
            self.alt_loc = ""

        self.residue_name = res_name
        self.chain_id = chain_id
        self.residue_number = res_seq
        self.icode = icode
        self.x_coord = x_coord
        self.y_coord = y_coord
        self.z_coord = z_coord
        self.occupancy = occupancy
        self.bfactor = bfactor
        self.element = element
        self.charge = charge



    ###################PROPERTIES AND SETTERS############################

    def __common_setter(self,value):
        value = str(value)
        value = value.replace(" ","")  # try removing this for speed
        return value


    @property
    def record_name(self):
        if self.__writing_entry==True:
            return self.__record_name.ljust(6)

        else:
            return self.__record_name

    @record_name.setter
    def record_name(self,record_name):
        record_name = self.__common_setter(record_name)
        self.__record_name = record_name

    @property
    def serial_number(self):
        if self.__writing_entry is True:
            return self.__serial_number.rjust(5)
        else:
            return self.__serial_number

    @serial_number.setter
    def serial_number(self, serial_number):
        try:
            ## If atom number is to large, use hex codes
            serial_number = int(serial_number)
            if serial_number >= 99999:
                serial_number = str(hex(serial_number))[2:]

        except:
            pass
        serial_number = self.__common_setter(serial_number)
        self.__serial_number = serial_number

    @property
    def atom_name(self):
        if self.__writing_entry is True:
            if  len(self.__atom_name) <= 3:
                return "  "+self.__atom_name.ljust(3)

            elif len(self.__atom_name) == 4:
                return " "+self.__atom_name.ljust(4)
            else:
                print("Error, AtomName longer than 4 characters")

        else:
            return self.__atom_name

    @atom_name.setter
    def atom_name(self,atom_name):
        atom_name = self.__common_setter(atom_name)
        self.__atom_name = atom_name

    @property
    def alt_loc(self):
        if self.__writing_entry==True:
            if self.__alt_loc_exists is True:
                return self.__alt_loc.rjust(1)
            else:
                return self.__alt_loc.rjust(1)
        else:
            return self.__alt_loc

    @alt_loc.setter
    def alt_loc(self, alt_loc):
        alt_loc = self.__common_setter(alt_loc)
        self.__alt_loc = alt_loc

    @property
    def residue_name(self):
        if self.__writing_entry==True:
            return self.__residue_name.rjust(3)
        else:
            return self.__residue_name

    @residue_name.setter
    def residue_name(self, residue_name):
        residue_name = self.__common_setter(residue_name)
        self.__residue_name = residue_name

    @property
    def chain_id(self):
        if self.__writing_entry==True:
            return self.__chain_id.rjust(2)
        else:
            return self.__chain_id

    @chain_id.setter
    def chain_id(self, chain_id):
        chain_id = self.__common_setter(chain_id)
        self.__chain_id = chain_id

    @property
    def residue_number(self):
        if self.__writing_entry==True:
            return self.__residue_number.rjust(4)
        else:
            return self.__residue_number

    @residue_number.setter
    def residue_number(self, residue_number):
        residue_number = self.__common_setter(residue_number)
        self.__residue_number = residue_number

    @property
    def icode(self):
        if self.__writing_entry==True:
            return self.__icode.rjust(4)
        else:
            return self.__icode

    @icode.setter
    def icode(self, icode):
        icode = self.__common_setter(icode)
        self.__icode = icode


    @property
    def x_coord(self):
        if self.__writing_entry==True:
            return self.__x_coord.rjust(8)
        else:
            return self.__x_coord

    @x_coord.setter
    def x_coord(self, x_coord):
        x_coord = float(x_coord)
        x_coord = "{:.3f}".format(x_coord) #trailing zeroes to 3 decimal points
        x_coord = self.__common_setter(x_coord)
        self.__x_coord = x_coord

    @property
    def y_coord(self):
        if self.__writing_entry==True:
            return self.__y_coord.rjust(8)
        else:
            return str(self.__y_coord)

    @y_coord.setter
    def y_coord(self, y_coord):
        y_coord = float(y_coord)
        y_coord = "{:.3f}".format(y_coord) #trailing zeroes to 3 decimal points
        y_coord = self.__common_setter(y_coord)
        self.__y_coord = y_coord

    @property
    def z_coord(self):
        if self.__writing_entry==True:
            return self.__z_coord.rjust(8)
        else:
            return self.__z_coord

    @z_coord.setter
    def z_coord(self, z_coord):
        z_coord = float(z_coord)
        z_coord = "{:.3f}".format(z_coord) #trailing zeroes to 3 decimal points
        z_coord = self.__common_setter(z_coord)
        self.__z_coord = z_coord

    @property
    def occupancy(self):
        if self.__writing_entry==True:
            return self.__occupancy.rjust(6)
        else:
            return self.__occupancy

    @occupancy.setter
    def occupancy(self, occupancy):
        occupancy = float(occupancy)
        occupancy = "{:.2f}".format(occupancy) #trailing zeroes to 2 decimal points
        occupancy = self.__common_setter(occupancy)
        self.__occupancy = occupancy

    @property
    def bfactor(self):
        if self.__writing_entry==True:
            return self.__bfactor.rjust(6)
        else:
            return self.__bfactor

    @bfactor.setter
    def bfactor(self, bfactor):
        bfactor = float(bfactor)
        bfactor = "{:.2f}".format(bfactor) #trailing zeroes to 2 decimal points
        bfactor = self.__common_setter(bfactor)
        self.__bfactor = bfactor

    @property
    def element(self):
        if self.__writing_entry==True:
            return self.__element.rjust(12)
        else:
            return self.__element

    @element.setter
    def element(self, element):
        element = self.__common_setter(element)
        self.__element = element

    @property
    def charge(self):
        if self.__writing_entry==True:
            return self.__charge.rjust(2)
        else:
            return self.__charge

    @charge.setter
    def charge(self, charge):
        charge = self.__common_setter(charge)
        self.__charge = charge

    ###################END PROPERTIES AND SETTERS############################


    def get_coord(self):
        x = float(self.x_coord)
        y = float(self.y_coord)
        z = float(self.z_coord)
        return np.array([x,y,z])




    def set_x_coord(self,value):
        self.x_coord = round(float(value),3)
    def set_y_coord(self,value):
        self.y_coord = round(float(value),3)
    def set_z_coord(self,value):
        self.z_coord = round(float(value),3)



    def generate_entry(self):
        self.__writing_entry = True
        entry =  self.record_name+self.serial_number+self.atom_name+self.alt_loc+self.residue_name+self.chain_id+self.residue_number+self.icode+self.x_coord+self.y_coord+self.z_coord+self.occupancy+self.bfactor+self.element+self.charge+""
        self.__writing_entry = False
        return entry
        # skipped alt_loc entry to accomodate hydrogen double digits

    def print(self):
        print(self.generate_entry().strip("\n"))



    def summary(self):
        print("Record Type: "+str(self.record_name))
        print("Atom Number: "+str(self.serial_number))
        print("Atom Name: "+str(self.atom_name))
        print("Alt Loc: "+str(self.alt_loc))
        print("Residue Name: "+str(self.residue_name))
        print("Chain ID: "+str(self.chain_id))
        print("Residue Seq Number: "+str(self.residue_number))
        print("Code Res Insert: "+str(self.icode))
        print("X Coordinate: "+str(self.x_coord))
        print("Y Coordinate: "+str(self.y_coord))
        print("Z Coordinate: "+str(self.z_coord))
        print("Occupancy: "+str(self.occupancy))
        print("TempFactor: "+str(self.bfactor))
        #print("Segment ID: "+str(self.segment_id))
        print("Element Symbol: "+str(self.element))
        print("Charge: "+str(self.charge))