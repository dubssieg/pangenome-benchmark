from sys import argv
from data_structures import IntPair,Interval
from gfagraphs import Graph
from datetime import datetime

# Takes as input two vcf files from vg deconstruct + 1 file from rs-pancat-compare
vcf_sim,vcf_real,tsv_pancat,gfa_sim,gfa_real = argv[1:]

print(f"[{datetime.now()}] Parsing rs-pancat-compare output...")
# Parse file from rs-pancat-compare to get path lengths
path_lengths:dict[str,int] = dict()
with open(tsv_pancat,'r',encoding='utf-8') as tsv_reader:
    for line in tsv_reader:
        if line.startswith('##'):
            name,length = line[2:].strip().split('\t')
            path_lengths[name] = int(length)

print(f"[{datetime.now()}] Processing variant files...")
# For each variant file, create a set of nodes to have node coverage
sim_vcf_set:set = set()
with open(vcf_sim,'r',encoding='utf-8') as vcf_reader:
    for line in vcf_reader:
        if line.startswith('#'):
            continue
        else:
            for nodes in line.split('\t')[7].split(';')[3][3:].split(','):
                sim_vcf_set.update([int(x) for x in nodes.replace('<','>')[1:].split('>')[1:-1]])

real_vcf_set:set = set()
with open(vcf_real,'r',encoding='utf-8') as vcf_reader:
    for line in vcf_reader:
        if line.startswith('#'):
            continue
        else:
            for nodes in line.split('\t')[7].split(';')[3][3:].split(','):
                real_vcf_set.update([int(x) for x in nodes.replace('<','>')[1:].split('>')[1:-1]])

print(f"[{datetime.now()}] Processing graph files...")
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

real_graph:Graph = Graph(gfa_file=gfa_real)
real_graph.sequence_offsets()
real_intervals:dict[str,list[IntPair]] = dict()
for node in real_vcf_set:
    for path_name in real_graph.segments[str(node)]['PO'].keys():
        for position_a,position_b,_ in real_graph.segments[str(node)]['PO'][path_name]:
            try:
                real_intervals[path_name].append(IntPair(position_a,position_b))
            except:
                real_intervals[path_name] = [IntPair(position_a,position_b)]
for path_name,intervals in real_intervals.items():
    real_intervals[path_name] = Interval(intervals)
real_basepairs = sum([x.total() for x in real_intervals.values()])

print(f"[{datetime.now()}] Computing intervals coverage...")
# Diff between intervals coverage
intersect_intervals = {name:sim_intervals[name].to_set().intersection(real_intervals[name].to_set()) for name in sim_intervals.keys()}
intersect_basepairs = sum([len(x) for x in intersect_intervals.values()])

# Results
print(f"[{datetime.now()}] Done!")
print(f"|A|={sim_basepairs} |B|={real_basepairs} |A∩B|={intersect_basepairs} |Ω|={sum(path_lengths.values())}")
print(f"Precision: {intersect_basepairs/real_basepairs}")
print(f"Recall: {intersect_basepairs/sim_basepairs}")