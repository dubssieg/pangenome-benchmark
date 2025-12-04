from sys import argv
from data_structures import IntPair,Interval
from gfagraphs import Graph
from datetime import datetime
from statistics import mean
from os import listdir

# Takes as input two vcf files from vg deconstruct + 1 file from rs-pancat-compare

input_path:str = argv[2]

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
    "P_var_mc",
    "P_var_pggb",
    "R_var_mc",
    "R_var_pggb",
    "E_mc",
    "E_pggb",
    "M_mc",
    "M_pggb",
    "S_mc",
    "S_pggb",
    "Rd_mc",
    "Rd_pggb",
    "Rb_mc",
    "Rb_pggb",
    "∆µ_mc",
    "∆µ_pggb",
    "P_dist_mc",
    "P_dist_pggb",
    "R_dist_mc",
    "R_dist_pggb",
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

    #print(f"[{datetime.now()}] Parsing rs-pancat-compare output...")
    # Parse file from rs-pancat-compare to get path lengths
    path_lengths:dict[str,int] = dict()
    with open(tsv_pancat_MC,'r',encoding='utf-8') as tsv_reader:
        for line in tsv_reader:
            if line.startswith('##'):
                name,length = line[2:].strip().split('\t')
                path_lengths[name] = int(length)
            elif line.startswith('# Distance:'):
                distance_string = line.strip()[line.index('(')+1:line.index(')')].replace(' ','')
                #print(distance_string)
    equiv_MC,split_MC,merge_MC = [int(x[2:]) for x in distance_string.split(',')[:-1]]
    with open(tsv_pancat_PGGB,'r',encoding='utf-8') as tsv_reader:
        for line in tsv_reader:
            if line.startswith('##'):
                name,length = line[2:].strip().split('\t')
                path_lengths[name] = int(length)
            elif line.startswith('# Distance:'):
                distance_string = line.strip()[line.index('(')+1:line.index(')')].replace(' ','')
                #print(distance_string)
    equiv_PGGB,split_PGGB,merge_PGGB = [int(x[2:]) for x in distance_string.split(',')[:-1]]

    #print(f"[{datetime.now()}] Processing variant files...")
    # For each variant file, create a set of nodes to have node coverage
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

    #print(f"[{datetime.now()}] Processing graph files...")
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

    real_intervals_MC:dict[str,list[IntPair]] = dict()
    real_graph_MC:Graph = Graph(gfa_file=gfa_real_MC)
    real_graph_MC.sequence_offsets()
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

    real_intervals_PGGB:dict[str,list[IntPair]] = dict()
    real_graph_PGGB:Graph = Graph(gfa_file=gfa_real_PGGB)
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

    #print(f"[{datetime.now()}] Computing intervals coverage...")
    # Diff between intervals coverage
    intersect_intervals_MC = {name:sim_intervals[name].to_set().intersection(real_intervals_MC[name].to_set()) for name in sim_intervals.keys()}
    intersect_basepairs_MC = sum([len(x) for x in intersect_intervals_MC.values()])

    intersect_intervals_PGGB = {name:sim_intervals[name].to_set().intersection(real_intervals_PGGB[name].to_set()) for name in sim_intervals.keys()}
    intersect_basepairs_PGGB = sum([len(x) for x in intersect_intervals_PGGB.values()])

    #print(f"[{datetime.now()}] Done!")
    # Variant coverage measures
    # Sp. = espèce
    # Cond. = condition expé
    # µ = taux de mutation
    # Rep. = N° de réplicat
    # |A| = simulated basepairs
    # |B_mc| = real basepairs mc
    # |B_pggb| = real basepairs pggb
    # |A∩B_mc| = intersect basepairs mc
    # |A∩B_pggb| = intersect basepairs pggb
    # |Ω| = total graph size
    # P_var_mc = Precision variants mc
    # P_var_pggb = Precision variants pggb
    # R_var_mc = Recall variants mc
    # R_var_pggb = Recall variants pggb
    specie:str = "yeast"
    condition:str = directory[5:directory.index('1')]
    mutation_rate:str = directory.split("_")[0][directory.index('1'):]
    replicate:str = directory.split("_")[1]

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
        intersect_basepairs_MC/real_basepairs_MC,
        intersect_basepairs_PGGB/real_basepairs_PGGB,
        intersect_basepairs_MC/sim_basepairs,
        intersect_basepairs_PGGB/sim_basepairs,
        equiv_MC,
        equiv_PGGB,
        merge_MC,
        merge_PGGB,
        split_MC,
        split_PGGB,
        (merge_MC+split_MC)/(merge_MC+split_MC+equiv_MC),
        (merge_PGGB+split_PGGB)/(merge_PGGB+split_PGGB+equiv_PGGB),
        split_MC/(merge_MC+split_MC),
        split_PGGB/(merge_PGGB+split_PGGB),
        mean([x['length'] for x in real_graph_MC.segments.values()]) - mean([x['length'] for x in sim_graph.segments.values()]),
        mean([x['length'] for x in real_graph_PGGB.segments.values()]) - mean([x['length'] for x in sim_graph.segments.values()]),
        equiv_MC/(equiv_MC+merge_MC),
        equiv_PGGB/(equiv_PGGB+merge_PGGB),
        equiv_MC/(equiv_MC+split_MC),
        equiv_PGGB/(equiv_PGGB+split_PGGB),
        sep=','
    )