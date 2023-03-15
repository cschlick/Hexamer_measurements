import numpy as np
from Bio.PDB import *
import tqdm
from pathlib import Path
import itertools
import subprocess
import argparse
from collections import defaultdict
import json



def float_dict(d):
    new_dict = {}
    for k,v in d.items():
        if type(v) == dict:
            new_dict[k] = float_dict(v)
        elif type(v)==str:
            new_dict[k] = v
        elif type(v)==list:
            new_dict[k] = [float(e) for e in v]
    return new_dict


def unit_vector(vector):
    return vector / np.linalg.norm(vector)

def angle_between(v1, v2,degrees=False):
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    a = np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))
    if degrees:
        a = np.degrees(a)
    return a




def get_residue_ca_coord(chain,residue_number):
    for r in chain.get_residues():
        res_number = int(r.get_full_id()[-1][1])
        if res_number ==residue_number:
            for a in r.get_atoms():
                if a.name == "CA":
                    coord = a.get_coord()
                    return coord

def measure_structure(s):
    """
    Process a biopython structure object. It assumes s is a T=4 HBV capsid, with 240 chains (60 A, 60 B, 60 C, 60 D)
    1. Get some "handle" coordinates, coordinates useful for clustering or making measurements
    2. Cluster into hexamer groups, uniquely label each chain
    3. Make measurements for distances and angles
    4  Return measurements as a dictionary
    """
    assert len(list(s.get_chains())) == 240
    
    # Get some handle coords, specific atoms or atom range means useful for figuring out quasiequivalent chains
    handle_dict = {
        "capsomer_center": {"residue_numbers": [132]},
        "dimer_center": {"residue_numbers": [61]},
        "interface_1": {"residue_numbers": [33, 34, 35, 36, 37, 38]},
        "interface_2": {"residue_numbers": [12, 13, 14, 15, 16, 17]},
        "helix_5": {"residue_numbers": [120, 121, 122, 123, 124, 125, 126, 127]},
        "tip": {"residue_numbers": list(range(61, 95))},
        "base": {"residue_numbers": list(range(1, 61)) + list(range(95, 142))},
    }

    # collect handle coords
    for handle_label, handle_info in handle_dict.items():

        handle_coords = []
        for c in s.get_chains():
            assert (
                len(list(c.get_residues())) == 149
            ), "Unexpected number of residues"
            residue_coords = []
            for r in c.get_residues():
                res_number = int(r.get_full_id()[-1][1])
                if res_number in handle_info["residue_numbers"]:
                    for a in r.get_atoms():
                        if a.name == "CA":
                            coord = a.get_coord()
                            residue_coords.append(coord)
            handle_coords.append(np.array(residue_coords).mean(axis=0))
        handle_info["handle_coords"] = np.array(handle_coords)
        
    # cluster interfaces
    from sklearn.neighbors import KDTree

    tree = KDTree(handle_dict["interface_1"]["handle_coords"])
    dists, inds = tree.query(handle_dict["interface_2"]["handle_coords"], k=1)
    # assert dists[:,0].max() < dists[:,1].min(), "Failed sanity check"
    interface_partner_map = {
        i: ind[0] for i, ind in enumerate(inds)
    }  # a dict to map a chain index (integer position in s.get_chains()) to it's partner
    interface_groups = list(
        set([frozenset([k, v]) for k, v in interface_partner_map.items()])
    )
    interface_groups = list([list(s) for s in interface_groups])
    assert len(interface_partner_map) == 240, "Failed sanity check"
    
    tree = KDTree(handle_dict["interface_1"]["handle_coords"])
    dists, inds = tree.query(handle_dict["interface_2"]["handle_coords"], k=1)
    interface_partner_map = {
        i: ind[0] for i, ind in enumerate(inds)
    }  # a dict to map a chain index (integer position in s.get_chains()) to it's partner
    
    def is_cap(i, j):
        # returns whether i is the cap in an i,j interface
        # assumes handle dict is defined in scope
        # If helix 5 of i is closer to interface 1 of j than to interface 2, the i is the cap
        helix5_mean_i = handle_dict["helix_5"]["handle_coords"][i]
        interface_1_j = handle_dict["interface_1"]["handle_coords"][j]
        interface_2_j = handle_dict["interface_2"]["handle_coords"][j]
        d_1 = np.linalg.norm(helix5_mean_i - interface_1_j)
        d_2 = np.linalg.norm(helix5_mean_i - interface_2_j)
        if d_1 < d_2:
            return True
        else:
            return False


    chain_to_cap_map = {}
    chain_to_base_map = {}

    for i, j in interface_partner_map.items():

        if is_cap(i, j):
            chain_to_base_map[i] = j
            chain_to_cap_map[j] = i
        elif is_cap(j, i):
            # i is cap
            chain_to_base_map[j] = i
            chain_to_cap_map[i] = j
            
            
    # cluster dimers
    from sklearn.neighbors import KDTree

    tree = KDTree(handle_dict["dimer_center"]["handle_coords"])
    dists, inds = tree.query(handle_dict["dimer_center"]["handle_coords"], k=2)
    dimer_partner_map = {
        i: ind[1] for i, ind in enumerate(inds)
    }  # a dict to map a chain index (integer position in s.get_chains()) to it's dimer partner
    dimer_groups = set([frozenset([k, v]) for k, v in dimer_partner_map.items()])
    dimer_groups = list([list(s) for s in dimer_groups])
    assert len(dimer_partner_map) == 240, "Failed sanity check"
    assert len(dimer_groups) == 120, "Failed sanity check"

    # Result of above section is dimer_partner_map
    
    
    # cluster hexamers
    tree = KDTree(handle_dict["capsomer_center"]["handle_coords"])
    dists, inds = tree.query(handle_dict["capsomer_center"]["handle_coords"], k=7)
    
    
    # optimize a cutoff distance that results in clean clustering of hexamers
    cutoffs = np.linspace(20,50,1000)

    for cutoff1 in cutoffs:
        count1 = sum(dists[:,:6].max(axis=1)<cutoff1)
        if count1==180:
            #print(cutoff1)
            break

    hex_inds = np.where(dists[:,:6].max(axis=1)<cutoff1)[0]
    for cutoff2 in cutoffs:
        found = dists[:,:5].max(axis=1)<cutoff2
        found[hex_inds]=False
        count2 = sum(found)
        if count2==60:
            #print(cutoff2)
            break
    
    pentamer_center_chains = []
    hexamer_center_chains = []

    for i, row in enumerate(dists):
        #assert max(row[:7]) < 70, "Failed Sanity check"
        #assert max(row[:5]) < 40, "Failed Sanity check"

        if max(row[:6]) < cutoff1:
            hexamer_center_chains.append(i)
        elif max(row[:5]) < cutoff2:
            pentamer_center_chains.append(i)

    pentamer_center_chains = set(pentamer_center_chains)
    hexamer_center_chains = set(hexamer_center_chains)
    assert len(pentamer_center_chains) == 60, "Failed sanity check, n pentamer centric chains: "+str(len(pentamer_center_chains))
    assert len(hexamer_center_chains) == 180, "Failed sanity check"

    # group pentamer chains into groups of chain indices that INCLUDE the outside dimers
    pentamer_index_groups = set()
    hexamer_index_groups = set()

    for i, ind in enumerate(inds):
        if i in pentamer_center_chains:
            pentamer_indices = list(inds[i, :5])
            # add partners
            for j in list(pentamer_indices):
                pentamer_indices.append(dimer_partner_map[j])
            pentamer_index_groups.add(frozenset(pentamer_indices))
        elif i in hexamer_center_chains:
            hexamer_indices = list(inds[i, :6])
            # add partners
            for j in list(hexamer_indices):
                hexamer_indices.append(dimer_partner_map[j])
            hexamer_index_groups.add(frozenset(hexamer_indices))
        else:
            assert False, "Failed Sanity check"
    pentamer_index_groups = [list(l) for l in pentamer_index_groups]
    hexamer_index_groups = [list(l) for l in hexamer_index_groups]
    #assert len(pentamer_index_groups) == 12, "Failed sanity check" # pentamers not working but not needed
    assert len(hexamer_index_groups) == 30, "Failed sanity check"
    
    
    # import itertools
    # test = list(
    #     itertools.chain.from_iterable(pentamer_index_groups + hexamer_index_groups)
    # )
    # assert len(test) == (12 * 5 * 2) + (30 * 6 * 2), "Failed sanity check"
    
    # cluster unique asus
    asu_groups = []
    chains = list(s.get_chains())
    for interface_group in interface_groups:
        chain1, chain2 = interface_group
        partner1 = dimer_partner_map[chain1]
        partner2 = dimer_partner_map[chain2]
        asu_group = [partner1, chain1, chain2, partner2]
        assert len(set(asu_group)) == 4, "Failed sanity check"

        chain_ids = [chains[i].id for i in asu_group]
        if set(chain_ids) == set(["A", "B", "C", "D"]):
            asu_groups.append(asu_group)

    # assert len(set([frozenset(g) for g in asu_groups])) == 60
    
    hexamer_chain_ids = []
    for hexamer_group in hexamer_index_groups:
        hexamer_chain_id = [chains[i].id for i in hexamer_group]
        assert (
            hexamer_chain_id.count("A") == 2
        ), "Hexamer has unexpected number of A chains"
        assert (
            hexamer_chain_id.count("B") == 2
        ), "Hexamer has unexpected number of B chains"
        assert (
            hexamer_chain_id.count("C") == 4
        ), "Hexamer has unexpected number of C chains"
        assert (
            hexamer_chain_id.count("D") == 4
        ), "Hexamer has unexpected number of D chains"
        hexamer_chain_ids.append(hexamer_chain_id)
        
    # Label each chain in the hexamer with a unique label
    distances = {"a":[], # B-B' pore
                 "b":[], # C-C' pore
                 "c":[], # D2-D2' pore
                 "d":[], # A-D, A'-D' expansion
                 "e":[], # C-C2', A'-C2 expansion
                 "f":[]} # C2-D, C2'-D' expansion

    angles = {"cd":[],
              "ab":[]}




    for hexamer_group in hexamer_index_groups:
        hexamer_group_labels = [None for chain_i in hexamer_group]
        for chain_i in hexamer_group:
            group_i = hexamer_group.index(chain_i)
            chain = chains[chain_i]
            chain_id = chain.id
            if chain_id == "B":
                if "B" not in hexamer_group_labels:
                    # this is the "anchor" B
                    chain_i_b = chain_i
                    group_i_b = group_i
                    hexamer_group_labels[group_i_b] = "B"

                    # get partner of B, A
                    chain_i_a = dimer_partner_map[chain_i_b]
                    assert (
                        chain_i_a in hexamer_group
                    ), "Failed sanity check: partner index not in hexamer group"
                    hexamer_group_labels[hexamer_group.index(chain_i_a)] = "A"

                    # Get cap of B, C
                    chain_i_c = chain_to_cap_map[chain_i_b]
                    assert (
                        chain_i_c in hexamer_group
                    ), "Failed sanity check: partner index not in hexamer group"
                    hexamer_group_labels[hexamer_group.index(chain_i_c)] = "C"

                    # Get C partner, D
                    chain_i_d = dimer_partner_map[chain_i_c]
                    assert (
                        chain_i_d in hexamer_group
                    ), "Failed sanity check: partner index not in hexamer group"
                    hexamer_group_labels[hexamer_group.index(chain_i_d)] = "D"

                    # Get the cap C, D2
                    chain_i_d_2 = chain_to_cap_map[chain_i_c]
                    assert (
                        chain_i_d in hexamer_group
                    ), "Failed sanity check: partner index not in hexamer group"
                    hexamer_group_labels[hexamer_group.index(chain_i_d_2)] = "D2"

                    # get partner of D2, C2
                    chain_i_c_2 = dimer_partner_map[chain_i_d_2]
                    assert (
                        chain_i_d in hexamer_group
                    ), "Failed sanity check: partner index not in hexamer group"
                    hexamer_group_labels[hexamer_group.index(chain_i_c_2)] = "C2"

                else:
                    pass
                    # this is the "anchor" B
                    chain_i_b = chain_i
                    group_i_b = group_i
                    hexamer_group_labels[group_i_b] = "B'"

                    # get partner of B, A
                    chain_i_a = dimer_partner_map[chain_i_b]
                    assert (
                        chain_i_a in hexamer_group
                    ), "Failed sanity check: partner index not in hexamer group"
                    hexamer_group_labels[hexamer_group.index(chain_i_a)] = "A'"

                    # Get cap of B, C
                    chain_i_c = chain_to_cap_map[chain_i_b]
                    assert (
                        chain_i_c in hexamer_group
                    ), "Failed sanity check: partner index not in hexamer group"
                    hexamer_group_labels[hexamer_group.index(chain_i_c)] = "C'"

                    # Get C partner, D
                    chain_i_d = dimer_partner_map[chain_i_c]
                    assert (
                        chain_i_d in hexamer_group
                    ), "Failed sanity check: partner index not in hexamer group"
                    hexamer_group_labels[hexamer_group.index(chain_i_d)] = "D'"

                    # Get the interface with C, D2
                    chain_i_d_2 = chain_to_cap_map[chain_i_c]
                    assert (
                        chain_i_d in hexamer_group
                    ), "Failed sanity check: partner index not in hexamer group"
                    hexamer_group_labels[hexamer_group.index(chain_i_d_2)] = "D2'"

                    # get partner of D2, C2
                    chain_i_c_2 = dimer_partner_map[chain_i_d_2]
                    assert (
                        chain_i_d in hexamer_group
                    ), "Failed sanity check: partner index not in hexamer group"
                    hexamer_group_labels[hexamer_group.index(chain_i_c_2)] = "C2'"


        hexamer_chain_dict = {}
        hexamer_index_dict = {}
        for chain_i, chain_label in zip(hexamer_group, hexamer_group_labels):
            io = PDBIO()
            chain = chains[chain_i]
            hexamer_chain_dict[chain_label] = chain
            hexamer_index_dict[chain_label] = chain_i
            #io.set_structure(chain)
            #io.save(chain_label + ".pdb")
        
        # # get hexamer center
        # atoms = itertools.chain.from_iterable([chain.get_atoms() for label,chain in hexamer_chain_dict.items()])
        # all_coords = np.vstack([a.get_coord() for a in atoms])
        # hexamer_center = all_coords.mean(axis=0)

        # value "a"
        v1 = get_residue_ca_coord(hexamer_chain_dict["B"],131)
        v2 = get_residue_ca_coord(hexamer_chain_dict["B'"],131)
        d = np.linalg.norm(v1-v2)
        distances["a"].append(d)

        # value b
        v1 = get_residue_ca_coord(hexamer_chain_dict["C"],131)
        v2 = get_residue_ca_coord(hexamer_chain_dict["C'"],131)
        d = np.linalg.norm(v1-v2)
        distances["b"].append(d)

        # value c
        v1 = get_residue_ca_coord(hexamer_chain_dict["D2"],131)
        v2 = get_residue_ca_coord(hexamer_chain_dict["D2'"],131)
        d = np.linalg.norm(v1-v2)
        distances["c"].append(d)

        # value d
        v1 = get_residue_ca_coord(hexamer_chain_dict["A"],131)
        v2 = get_residue_ca_coord(hexamer_chain_dict["D"],131)
        d = np.linalg.norm(v1-v2)
        distances["d"].append(d)

        v1 = get_residue_ca_coord(hexamer_chain_dict["A'"],131)
        v2 = get_residue_ca_coord(hexamer_chain_dict["D'"],131)
        d = np.linalg.norm(v1-v2)
        distances["d"].append(d)

        # value e
        v1 = get_residue_ca_coord(hexamer_chain_dict["A"],131)
        v2 = get_residue_ca_coord(hexamer_chain_dict["C2'"],131)
        d = np.linalg.norm(v1-v2)
        distances["e"].append(d)

        v1 = get_residue_ca_coord(hexamer_chain_dict["A'"],131)
        v2 = get_residue_ca_coord(hexamer_chain_dict["C2"],131)
        d = np.linalg.norm(v1-v2)
        distances["e"].append(d)

        # value f
        v1 = get_residue_ca_coord(hexamer_chain_dict["C2"],131)
        v2 = get_residue_ca_coord(hexamer_chain_dict["D"],131)
        d = np.linalg.norm(v1-v2)
        distances["f"].append(d)

        v1 = get_residue_ca_coord(hexamer_chain_dict["C2'"],131)
        v2 = get_residue_ca_coord(hexamer_chain_dict["D'"],131)
        d = np.linalg.norm(v1-v2)
        distances["f"].append(d)

        # angles
        #group dimers in this hexamer
        dimer1 = ["A","B"]
        dimer2 = ["C","D"]
        dimer3 = ["D2","C2"]
        dimer4 = ["B'","A'"]
        dimer5 = ["C'","D'"]
        dimer6 = ["C2'","D2'"]
        dimers = [dimer1,dimer2,dimer3,dimer4,dimer5,dimer6]
        
        # get hexamer base center
        hexamer_bases = []
        for label1,label2 in dimers:
            chain_i, chain_j = hexamer_index_dict[label1], hexamer_index_dict[label2]
            base_i, base_j = handle_dict["base"]["handle_coords"][chain_i], handle_dict["base"]["handle_coords"][chain_j]
            base = (base_i+base_j)/2
            hexamer_bases.append(base)
        hexamer_center = np.array(hexamer_bases).mean(axis=0)
        
        for label1,label2 in dimers:
            chain_i, chain_j = hexamer_index_dict[label1], hexamer_index_dict[label2]
            base_i, base_j = handle_dict["base"]["handle_coords"][chain_i], handle_dict["base"]["handle_coords"][chain_j]
            tip_i, tip_j = handle_dict["tip"]["handle_coords"][chain_i], handle_dict["tip"]["handle_coords"][chain_j]
            base = (base_i+base_j)/2
            tip = (tip_i+tip_j)/2
            v1 = hexamer_center-base
            v2 = tip-base
            a = angle_between(v1,v2,degrees=True)
            labels = "".join([label1,label2])
            if "C" in labels:
                angles["cd"].append(a)
            elif "B" in labels:
                angles["ab"].append(a)
    data = {"distances":distances,"angles":angles}
    return data



def main():


    return_data = []
    combined_data = {"distances":{"a":[], # B-B' pore
                     "b":[], # C-C' pore
                     "c":[], # D2-D2' pore
                     "d":[], # A-D, A'-D' expansion
                     "e":[], # C-C2', A'-C2 expansion
                     "f":[]}, # C2-D, C2'-D' expansion
                     "angles":{"cd":[],"ab":[]}}
                 
    parser = argparse.ArgumentParser()
    parser.add_argument("data_dir",help="Directory to look for pdb files. Will assume 'CHIMERA' in filename.")
    args = parser.parse_args()
    print(args.data_dir)
    data_dir = Path(args.data_dir)
    files = [file for file in data_dir.glob("*") if "CHIMERA" in file.name and file.suffix==".pdb"]
    print("Files")
    print(files)
    for f in tqdm.tqdm(files):
        try:
            parser = PDBParser(QUIET=True)
            outpath = Path(f.parent,f.name.replace("pdb","")+"measurements.json")
            if not outpath.exists():
                s = parser.get_structure(str(f),str(f))
                ret = measure_structure(s)

                return_data.append(ret)

                # merge data over multiple structures
                for key,value in ret["distances"].items():
                    combined_data["distances"][key]+=value
                for key,value in ret["angles"].items():
                    combined_data["angles"][key]+=value


                with outpath.open("w") as fh:
                    json.dump(float_dict(ret),fh)
        except:
            print("Failed: "+str(f))
            raise

            
if __name__ == "__main__":
    main()