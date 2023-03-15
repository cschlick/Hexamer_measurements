# Hexamer measurements
A place to store scripts for measuring hexamer structural properties for HBV capsid


### Install external depencies with conda
```bash
conda install numpy scipy matplotlib seaborn scikit-learn biopython tqdm
```

### Convert files
This script uses some old code to re-format a pdb file. It uses a parameter file with information about the capsid to split the file into multiple "models" using the pdb MODEL ENDMDL flags. This is useful for displaying the full model in in viewers like UCSF Chimera, which struggle to interpret multiple identical chain IDs in the same model. This happens to be useful for the subsequent script which uses the Biopython PDB package. Here we convert all pdb files in the ../data folder.
```bash
cd hexamer_measurements
python convert_files.py ../data
```

The result of this is that all *.pdb files are re-formatted to *_CHIMERA.pdb files. This will serve as input for subsequent processing.

### Measure hexamer properties
The second step is to open the capsid structures, cluster the chains into 30 "hexamers" and then measure properties for each of those hexamers. The measurements are written to a *measurements.json file with the same prefix as the pdb file. 
```bash
python measure_hexamers.py ../data
```

### Plot measurments
The final script collects all the measurements from the .json files, and plots them. Also specify the name/location to write the figure.
```bash
python analyse_and_plot.py ../data test_fig
```



