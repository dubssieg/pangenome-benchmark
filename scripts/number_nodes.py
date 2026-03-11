from sys import argv
from data_structures import IntPair,Interval
from gfagraphs import Graph
from datetime import datetime
from statistics import mean
from os import listdir

input_path:str = argv[1]

print(
    "Sp.",
    "Cond.",
    "µ",
    "Rep.",
    "nb_nodes_ms",
    "nb_nodes_mc",
    "nb_nodes_pggb",
    sep=','
)

for directory in listdir(input_path):

    base_path:str = f"{input_path}/{directory}/"
    vcf_sim:str = base_path + "variants/LN0#1#0.mspangepop.vcf"
    vcf_real_MC:str = base_path + "variants/LN0#1#0.mc.vcf"
    vcf_real_PGGB:str = base_path + "variants/LN0#1#0.pggb.vcf"
    tsv_pancat_MC:str = base_path + "dist.ms.mc.tsv"
    tsv_pancat_PGGB:str = base_path + "dist.ms.pggb.tsv"
    gfa_sim:str = base_path + "graph.mspangepop.gfa"
    gfa_real_MC:str = base_path + "graph.mc.gfa"
    gfa_real_PGGB:str = base_path + "graph.pggb.gfa"

    nb_nodes_ms:int = 0
    with open(gfa_sim,'r') as reader:
        for line in reader:
            if line.startswith('S'):
                nb_nodes_ms +=1

    nb_nodes_pggb:int = 0
    with open(gfa_real_PGGB,'r') as reader:
        for line in reader:
            if line.startswith('S'):
                nb_nodes_pggb +=1

    nb_nodes_mc:int = 0
    with open(gfa_real_MC,'r') as reader:
        for line in reader:
            if line.startswith('S'):
                nb_nodes_mc +=1

    specie:str = directory.split("_")[0]
    condition:str = directory.split("_")[1]
    mutation_rate:str = "1e7"
    try:
        replicate:str = directory.split("_")[2]
    except IndexError:
        replicate:str = "rep0"

    print(
        specie,
        condition,
        mutation_rate,
        replicate,
        nb_nodes_ms,
        nb_nodes_mc,
        nb_nodes_pggb,
        sep=','
    )