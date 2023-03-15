import json
from pathlib import Path
import seaborn as sns
import tqdm
import matplotlib.pyplot as plt
import numpy as np
import argparse
from matplotlib.ticker import FormatStrFormatter


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("data_dir",help="directory to look for json files.")
    parser.add_argument("output_figure",help="The path to save the resulting figure.")


    args = parser.parse_args()
    print(args.data_dir)
    data_dir = Path(args.data_dir)
    files = [file for file in data_dir.glob("*") if file.suffix==".json"]
    print("Files")
    print(files)
    
    combined_data = {"distances":{"a":[], # B-B' pore
                 "b":[], # C-C' pore
                 "c":[], # D2-D2' pore
                 "d":[], # A-D, A'-D' expansion
                 "e":[], # A-C2', A'-C2 expansion
                 "f":[]}, # C2-D, C2'-D' expansion
                 "angles":{"cd":[],"ab":[]}}
                 
    for file in tqdm.tqdm(files):
        with file.open("r") as fh:
            data = json.load(fh)

        for key,value in data["distances"].items():
            combined_data["distances"][key]+=value
        for key,value in data["angles"].items():
            combined_data["angles"][key]+=value


    SMALL_SIZE = 16
    MEDIUM_SIZE = 22
    BIGGER_SIZE = 44

    plt.rc('font', size=22)          # controls default text sizes
    plt.rc('axes', titlesize=22)     # fontsize of the axes title
    plt.rc('axes', labelsize=18)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize


    fig, axs = plt.subplots(2,4,figsize=(22,10))
    axs = axs.flatten()



    stat = "density"

    # first plot
    ax = axs[0]
    vals = combined_data["distances"]["a"]
    sns.histplot(vals,ax=ax,stat=stat,linewidth=0)
    ax.set_xlabel("Interdimer Residue 131 CA distance (Å)")
    ax.set_title("BB Pore (Å)\nMean: "+str(round(np.mean(vals),1)))
    ax.set_ylabel("")
    ax.set_ylim(0,0.3)
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax.set_xlim(15,40)




    # second plot
    ax = axs[1]
    vals = combined_data["distances"]["b"]
    sns.histplot(vals,ax=ax,stat=stat,linewidth=0)
    ax.set_xlabel("Interdimer Residue 131 CA distance (Å)")
    ax.set_title("CC Pore (Å)\nMean: "+str(round(np.mean(vals),1)))
    ax.set_ylabel("")
    ax.set_ylim(0,0.3)
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax.set_xlim(15,40)




    # third plot
    ax = axs[2]
    vals = combined_data["distances"]["c"]
    sns.histplot(vals,ax=ax,stat=stat,linewidth=0)
    ax.set_xlabel("Interdimer Residue 131 CA distance (Å)")
    ax.set_title("DD Pore (Å)\nMean: "+str(round(np.mean(vals),1)))
    ax.set_ylabel("")
    ax.set_ylim(0,0.3)
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax.set_xlim(15,40)



    # fourth plot
    ax = axs[3]
    vals = combined_data["distances"]["d"]
    sns.histplot(vals,ax=ax,stat=stat,linewidth=0)
    ax.set_xlabel("Interdimer Residue 131 CA distance (Å)")
    ax.set_title("AD expansion (Å)\nMean: "+str(round(np.mean(vals),1)))
    ax.set_ylabel("")
    ax.set_ylim(0,0.3)
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax.set_xlim(75,95)



    # fifth plot
    ax = axs[4]
    vals = combined_data["distances"]["e"]
    sns.histplot(vals,ax=ax,stat=stat,linewidth=0)
    ax.set_xlabel("Interdimer Residue 131 CA distance (Å)")
    ax.set_title("AC expansion (Å)\nMean: "+str(round(np.mean(vals),1)))
    ax.set_ylabel("")
    ax.set_ylim(0,0.3)
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax.set_xlim(75,95)

    #  sixth plot
    ax = axs[5]
    vals = combined_data["distances"]["f"]
    sns.histplot(vals,ax=ax,stat=stat,linewidth=0)
    ax.set_xlabel("Interdimer Residue 131 CA distance (Å)")
    ax.set_title("CD expansion (Å)\nMean: "+str(round(np.mean(vals),1)))
    ax.set_ylabel("")
    ax.set_ylim(0,0.3)
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax.set_xlim(75,95)

    # seventh plot
    ax = axs[6]
    vals = combined_data["angles"]["ab"]
    sns.histplot(vals,ax=ax,stat=stat,linewidth=0)
    ax.set_xlabel("Spike-tip/hexamer internal angle (°)")
    ax.set_title("AB angle (°)\nMean: "+str(round(np.mean(vals),1)))
    ax.set_ylabel("")
    ax.set_ylim(0,0.15)
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax.set_xlim(85,120)


    # eighth plot
    ax = axs[7]
    vals = combined_data["angles"]["cd"]
    sns.histplot(vals,ax=ax,stat=stat,linewidth=0)
    ax.set_xlabel("Spike-tip/hexamer internal angle (°)")
    ax.set_title("CD angle (°)\nMean: "+str(round(np.mean(vals),1)))
    ax.set_ylabel("")       
    ax.set_ylim(0,0.15)
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax.set_xlim(85,120)

    plt.tight_layout()
    
    outpath= Path(args.output_figure)
    if outpath.exists():
        outpath.unlink()
    plt.savefig(args.output_figure)



            
if __name__ == "__main__":
    main()