#!/usr/bin/env python3
from sys import argv
from json import dump
from typing import Any
from statistics import mean, stdev
from os import listdir
from copy import deepcopy

input_path:str = argv[1]
sampling: int = 200
all_dicts:list = list()

for directory in listdir(input_path):

    base_path:str = f"{input_path}/{directory}"

    specie:str = directory.split("_")[0]
    condition:str = directory.split("_")[1]
    mutation_rate:str = "1e7"
    try:
        replicate:str = directory.split("_")[2]
    except IndexError:
        replicate:str = "rep0"

    for soft in ['mc','pggb']:
        file_name:str = f"dist.ms.{soft}.tsv"

        edition_results: dict = {}
        merges_results: dict = {}
        splits_results: dict = {}
        ratios_paths: dict = {}
        path_lengths: dict = {}

        edits_results: dict[str, Any] = {
            'specie':specie,
            'condition':condition,
            'mutation-rate':mutation_rate,
            'replicate':replicate,
            'software': soft,
            'ratios': 0,
            'lengths': 0,
            'editions-mean': {},
            'merges-mean': {},
            'splits-mean': {},
            'editions-stdev': {},
            'merges-stdev': {},
            'splits-stdev': {},
        }

        with open(f"{base_path}/{file_name}", 'r', encoding='utf-8') as f:
            # Skip the header
            next(f)
            # Read the path sizes
            while (l := next(f)).startswith('##'):
                path_name, path_length = l[3:].strip().split('\t')
                path_lengths[path_name] = int(path_length)
                ratios_paths[path_name] = sampling/int(path_length)
            # Skip the header
            next(f)
            # Init the results
            for path_name in ratios_paths.keys():
                edition_results[path_name] = [0] * sampling
                merges_results[path_name] = [0] * sampling
                splits_results[path_name] = [0] * sampling
            # Read the data
            for l in f:
                if l.startswith('#'):
                    continue
                else:
                    path_name, position, operation, nodeA, nodeB, breakpointA, breakpointB = l.strip().split(
                        '\t')
                    if operation == 'S':
                        splits_results[path_name][int(
                            int(position)*ratios_paths[path_name])] += 1
                    elif operation == 'M':
                        merges_results[path_name][int(
                            int(position)*ratios_paths[path_name])] += 1
                    edition_results[path_name][int(
                        int(position)*ratios_paths[path_name])] += 1


        edits_results['ratios'] = mean(ratios_paths.values())
        edits_results['lengths'] = mean(path_lengths.values())
        edits_results['editions-mean'] = [mean([edition_results[y][x] for y in path_lengths.keys()]) for x in range(sampling)]
        edits_results['merges-mean'] = [mean([merges_results[y][x] for y in path_lengths.keys()]) for x in range(sampling)]
        edits_results['splits-mean'] = [mean([splits_results[y][x] for y in path_lengths.keys()]) for x in range(sampling)]
        edits_results['editions-stdev'] = [stdev([edition_results[y][x] for y in path_lengths.keys()]) for x in range(sampling)]
        edits_results['merges-stdev'] = [stdev([merges_results[y][x] for y in path_lengths.keys()]) for x in range(sampling)]
        edits_results['splits-stdev'] = [stdev([splits_results[y][x] for y in path_lengths.keys()]) for x in range(sampling)]

        all_dicts.append(deepcopy(edits_results))

# Save the results
with open(f'editions_{input_path}.json', 'w', encoding='utf-8') as f:
    dump(all_dicts, f)
