from subprocess import run,PIPE
from sys import argv
from os import listdir

input_path:str = argv[1]

print(
    "Sp.",
    "Cond.",
    "µ",
    "Rep.",
    "mash",
    sep=','
)

for directory in listdir(input_path):

    base_path:str = f"{input_path}/{directory}/singlefasta"
    mash:list[float] = list()

    for i,path_name_A in enumerate(listdir(base_path)):
        for j,path_name_B in enumerate(listdir(base_path)):
            if i>j:
                mash.append(float(run(['mash', 'dist',f'{base_path}/{path_name_A}',f'{base_path}/{path_name_B}'], stdout=PIPE).stdout.decode('utf-8').split()[2]))
            
    mean_mash:float = sum(mash)/len(mash)

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
        mean_mash,
        sep=','
    )