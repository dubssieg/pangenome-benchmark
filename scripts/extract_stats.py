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
    "|A|",
    "|B_mc|",
    "|B_pggb|",
    "|A∩B_mc|",
    "|A∩B_pggb|",
    "|Ω|",
    "E_mc",
    "E_pggb",
    "M_mc",
    "M_pggb",
    "S_mc",
    "S_pggb",
    "∆µ_mc",
    "∆µ_pggb",
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

    # Parse file from rs-pancat-compare to get path lengths
    path_lengths:dict[str,int] = dict()
    with open(tsv_pancat_MC,'r',encoding='utf-8') as tsv_reader:
        for line in tsv_reader:
            if line.startswith('##'):
                name,length = line[2:].strip().split('\t')
                path_lengths[name] = int(length)
            elif line.startswith('# Distance:'):
                distance_string = line.strip()[line.index('(')+1:line.index(')')].replace(' ','')
    try:
        equiv_MC,split_MC,merge_MC = [int(x[2:]) for x in distance_string.split(',')[:-1]]
    except NameError:
        equiv_MC,split_MC,merge_MC = 0,0,0
    with open(tsv_pancat_PGGB,'r',encoding='utf-8') as tsv_reader:
        for line in tsv_reader:
            if line.startswith('##'):
                name,length = line[2:].strip().split('\t')
                path_lengths[name] = int(length)
            elif line.startswith('# Distance:'):
                distance_string = line.strip()[line.index('(')+1:line.index(')')].replace(' ','')
    try:
        equiv_PGGB,split_PGGB,merge_PGGB = [int(x[2:]) for x in distance_string.split(',')[:-1]]
    except NameError:
        equiv_PGGB,split_PGGB,merge_PGGB = 0,0,0

    # For each variant file, create a set of nodes to have node coverage
    try:
        sim_vcf_set:set = set()
        with open(vcf_sim,'r',encoding='utf-8') as vcf_reader:
            for line in vcf_reader:
                if line.startswith('#'):
                    continue
                else:
                    for nodes in line.split('\t')[7].split(';')[3][3:].split(','):
                        sim_vcf_set.update([int(x) for x in nodes.replace('<','>')[1:].split('>')[1:-1]])

        real_vcf_set_MC:set = set()
        with open(vcf_real_MC,'r',encoding='utf-8') as vcf_reader:
            for line in vcf_reader:
                if line.startswith('#'):
                    continue
                else:
                    for nodes in line.split('\t')[7].split(';')[3][3:].split(','):
                        real_vcf_set_MC.update([int(x) for x in nodes.replace('<','>')[1:].split('>')[1:-1]])

        real_vcf_set_PGGB:set = set()
        with open(vcf_real_PGGB,'r',encoding='utf-8') as vcf_reader:
            for line in vcf_reader:
                if line.startswith('#'):
                    continue
                else:
                    for nodes in line.split('\t')[7].split(';')[3][3:].split(','):
                        real_vcf_set_PGGB.update([int(x) for x in nodes.replace('<','>')[1:].split('>')[1:-1]])
    except FileNotFoundError:
        sim_vcf_set,real_vcf_set_MC,real_vcf_set_PGGB = set(),set(),set()

    # Load the gfa to transform node coverage to base coverage (intervals)

    sim_graph:Graph = Graph(gfa_file=gfa_sim)
    sim_graph.sequence_offsets()
    sim_intervals:dict[str,list[IntPair]] = dict()
    for node in sim_vcf_set:
        for path_name in sim_graph.segments[str(node)]['PO'].keys():
            for position_a,position_b,_ in sim_graph.segments[str(node)]['PO'][path_name]:
                try:
                    sim_intervals[path_name].append(IntPair(position_a,position_b))
                except:
                    sim_intervals[path_name] = [IntPair(position_a,position_b)]
    for path_name,intervals in sim_intervals.items():
        sim_intervals[path_name] = Interval(intervals)
    sim_basepairs = sum([x.total() for x in sim_intervals.values()])

    real_graph_MC:Graph = Graph(gfa_file=gfa_real_MC)
    real_graph_MC.sequence_offsets()
    real_intervals_MC:dict[str,list[IntPair]] = dict()
    for node in real_vcf_set_MC:
        for path_name in real_graph_MC.segments[str(node)]['PO'].keys():
            for position_a,position_b,_ in real_graph_MC.segments[str(node)]['PO'][path_name]:
                try:
                    real_intervals_MC[path_name].append(IntPair(position_a,position_b))
                except:
                    real_intervals_MC[path_name] = [IntPair(position_a,position_b)]
    for path_name,intervals in real_intervals_MC.items():
        real_intervals_MC[path_name] = Interval(intervals)
    real_basepairs_MC = sum([x.total() for x in real_intervals_MC.values()])

    real_graph_PGGB:Graph = Graph(gfa_file=gfa_real_PGGB)
    real_intervals_PGGB:dict[str,list[IntPair]] = dict()
    real_graph_PGGB.sequence_offsets()
    for node in real_vcf_set_PGGB:
        for path_name in real_graph_PGGB.segments[str(node)]['PO'].keys():
            for position_a,position_b,_ in real_graph_PGGB.segments[str(node)]['PO'][path_name]:
                try:
                    real_intervals_PGGB[path_name].append(IntPair(position_a,position_b))
                except:
                    real_intervals_PGGB[path_name] = [IntPair(position_a,position_b)]
    for path_name,intervals in real_intervals_PGGB.items():
        real_intervals_PGGB[path_name] = Interval(intervals)
    real_basepairs_PGGB = sum([x.total() for x in real_intervals_PGGB.values()])

    # Diff between intervals coverage
    intersect_intervals_MC = {name:sim_intervals[name].to_set().intersection(real_intervals_MC[name].to_set() if name in real_intervals_MC else set()) for name in sim_intervals.keys()}
    intersect_basepairs_MC = sum([len(x) for x in intersect_intervals_MC.values()])

    intersect_intervals_PGGB = {name:sim_intervals[name].to_set().intersection(real_intervals_PGGB[name].to_set() if name in real_intervals_PGGB else set()) for name in sim_intervals.keys()}
    intersect_basepairs_PGGB = sum([len(x) for x in intersect_intervals_PGGB.values()])

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
        sim_basepairs,
        real_basepairs_MC,
        real_basepairs_PGGB,
        intersect_basepairs_MC,
        intersect_basepairs_PGGB,
        sum(path_lengths.values()),
        equiv_MC,
        equiv_PGGB,
        merge_MC,
        merge_PGGB,
        split_MC,
        split_PGGB,
        mean([x['length'] for x in real_graph_MC.segments.values()]) - mean([x['length'] for x in sim_graph.segments.values()]),
        mean([x['length'] for x in real_graph_PGGB.segments.values()]) - mean([x['length'] for x in sim_graph.segments.values()]),
        sep=','
    )