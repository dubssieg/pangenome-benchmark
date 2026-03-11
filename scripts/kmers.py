from Bio import SeqIO
from collections import Counter
from sys import argv
from os import listdir

input_path:str = argv[1]
k:int = 21

print(
    "Sp.",
    "Cond.",
    "Rep.",
    "id",
    "n",
    sep=','
)

for directory in listdir(input_path):

    base_path:str = f"{input_path}/{directory}/singlefasta"

    specie:str = directory.split("_")[0]
    condition:str = directory.split("_")[1]
    mutation_rate:str = "1e7"
    try:
        replicate:str = directory.split("_")[2]
    except IndexError:
        replicate:str = "rep0"

    for i,path_name in enumerate(listdir(base_path)):
        for record in SeqIO.parse(open(f"{base_path}/{path_name}"),"fasta"):
            counts = Counter([record.seq[i:i+k] for i in range(len(record)-k)])
            # keep only non-unique kmers
            n:float = float(sum({x: count for x, count in counts.items() if count > 1}.values()))/float(len(record)-k)

            print(
                specie,
                condition,
                replicate,
                record.id,
                n,
                sep=','
            )




