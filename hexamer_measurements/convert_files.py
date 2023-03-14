from pathlib import Path
import tqdm
import subprocess
import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("data_dir",help="directory to look for pdb files.")
    args = parser.parse_args()
    print(args.data_dir)
    data_dir = Path(args.data_dir)
    files = [file for file in data_dir.glob("*") if "CHIMERA" not in file.name and file.suffix==".pdb"]
    print("Files")
    print(files)
    for file in tqdm.tqdm(files):
        try:
            subprocess.run(["python","../PDBModule/pdbRefactor.py",str(file),"type=CHIMERA"])
        except:
            print("File failed:"+str(file))

            
if __name__ == "__main__":
    main()