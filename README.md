# Hexamer measurements
A place to store scripts for measuring hexamer structural properties for HBV capsid


### Install external depencies with conda
```bash
conda install numpy scipy matplotlib seaborn scikit-learn biopython tqdm
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

The measurements json file is a dictionary which includes a list of all the named measurements specified below:
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




### Plot measurments
The final script collects all the measurements from the .json files, and plots them. Also specify the name/location to write the figure.
```bash
python analyse_and_plot.py ../data test_fig
```



