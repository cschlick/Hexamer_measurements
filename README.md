# Hexamer measurements
A place to store scripts for measuring hexamer structural properties for HBV capsid


### Install external depencies with conda
```bash
conda install numpy scipy matplotlib seaborn scikit-learn biopython tqdm
```

### Clone this repository
```bash
git clone https://github.com/cschlick/Hexamer_measurements.git
cd Hexamer_measurements

```

### Convert files
This script uses some old code to re-format a pdb file. It uses a parameter file with information about the capsid to split the file into multiple "models" using the pdb MODEL ENDMDL tags. This is useful for displaying the full model in in viewers like UCSF Chimera, which seem to give display errors with multiple identical chain IDs in the same model. This also happens to be useful for the subsequent script which uses the Biopython PDB package, in which structures need to have either unique chain IDs or be separated into models. Here we convert all pdb files in the ../data folder.
```bash
cd hexamer_measurements
python convert_files.py ../data
```

The result of this is that all *.pdb files are re-formatted to *_CHIMERA.pdb files. This will serve as input for subsequent processing.

### Measure hexamer properties
The second step is to open the capsid structures, cluster the chains into 30 "hexamers" and then measure properties for each of those hexamers. The script is hard-coded to look for *_CHIMERA.pdb files. The measurements are written to a *measurements.json file with the same prefix as the pdb file. 
```bash
python measure_hexamers.py ../data
```

The measurements json file is a dictionary which includes a list of 8 named measurements.:
```python
{"distances":{"a":[], # B-B' pore
                 "b":[], # C-C' pore
                 "c":[], # D2-D2' pore
                 "d":[], # A-D, A'-D' expansion
                 "e":[], # A-C2', A'-C2 expansion
                 "f":[]}, # C2-D, C2'-D' expansion
 "angles":{"cd":[], # CD dimer angle
           "ab":[]}} # ab dimer angle
 ```

![alt text](https://github.com/cschlick/Hexamer_measurements/blob/e0705ed558d9f63ae154127a2d5109a915dff32a/HexamerDiagram.jpg)

Each distance measurement is defined below:
 - a: The distance between residue 131 CA atoms between chain B and B'
 - b: The distance between residue 131 CA atoms between chain C and C'
 - c: The distance between residue 131 CA atoms between chain D2 and D2'
 - d: The distance between residue 131 CA atoms between chain A and D, AND between chain A' and D'
 - e: The distance between residue 131 CA atoms between chain A and C2', AND between chain A' and C2
 - f: The distance between residue 131 CA atoms between chain C2 and D, AND between chain D' and C2'

For angles, three points are defined for each dimer:
1. The center of the hexamer (constant for all dimer in the hexamer)
2. The center of the dimer "base" region (Mean of all coordinates in residue range: 1-61 and 95, 142
3. The center of the dimer "tip" region (Mean of all coordinates in residue range: 61-95

The angle is defined as the angle between two vectors: base to hexamer center, and base to tip. Angles for all C,D chains are combined, as are angles for all A,B chains around a hexamer.


### Plot measurments
The final script collects all the measurements from the .json files, and plots them. Also specify the name/location to write the figure.
```bash
python analyse_and_plot.py ../data test_fig
```

### Run on a non-capsid hexamer
To measure a pdb file that has a single hexamer (not a T=4 capsid), such as might come from an experimental structure, use a modified version of the script.
```bash
python measure_single_hexamer.py ../experimental_data
```
This should write a json file with the sample contents:

```Python
{"distances": {"a": [35.41972732543945], 
               "b": [25.42247772216797],
               "c": [35.75651550292969],
               "d": [79.69721221923828,79.6960678100586],
               "e": [81.47354888916016, 81.47306823730469], 
               "f": [82.27655792236328, 82.2757797241211]}, 
 "angles": {"cd": [99.78011102553393, 100.80055352499853, 99.7804021241633, 100.8006734705682], 
            "ab": [97.28775036421341, 97.28759887629494]}}
```
