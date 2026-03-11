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
    "nb_variants_ms",
    "nb_variants_mc",
    "nb_variants_pggb",
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

    # For each variant file, create a set of nodes to have node coverage
    nb_var_ms:int = 0
    nb_var_mc:int = 0
    nb_var_pggb:int = 0
    try:
        sim_vcf_set:set = set()
        with open(vcf_sim,'r',encoding='utf-8') as vcf_reader:
            for line in vcf_reader:
                if line.startswith('#'):
                    continue
                else:
                    nb_var_ms += 1
                    for nodes in line.split('\t')[7].split(';')[3][3:].split(','):
                        sim_vcf_set.update([int(x) for x in nodes.replace('<','>')[1:].split('>')[1:-1]])

        real_vcf_set_MC:set = set()
        with open(vcf_real_MC,'r',encoding='utf-8') as vcf_reader:
            for line in vcf_reader:
                if line.startswith('#'):
                    continue
                else:
                    nb_var_mc += 1
                    for nodes in line.split('\t')[7].split(';')[3][3:].split(','):
                        real_vcf_set_MC.update([int(x) for x in nodes.replace('<','>')[1:].split('>')[1:-1]])

        real_vcf_set_PGGB:set = set()
        with open(vcf_real_PGGB,'r',encoding='utf-8') as vcf_reader:
            for line in vcf_reader:
                if line.startswith('#'):
                    continue
                else:
                    nb_var_pggb += 1
                    for nodes in line.split('\t')[7].split(';')[3][3:].split(','):
                        real_vcf_set_PGGB.update([int(x) for x in nodes.replace('<','>')[1:].split('>')[1:-1]])
    except FileNotFoundError:
        sim_vcf_set,real_vcf_set_MC,real_vcf_set_PGGB = set(),set(),set()

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
        nb_var_ms,
        nb_var_mc,
        nb_var_pggb,
        sep=','
    )