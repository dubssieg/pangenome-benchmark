#!/bin/bash
#SBATCH --job-name=bmMSpp
#SBATCH --constraint avx2
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G


###################################################################
## Author: Siegfried Dubois
##
## This script exectues a workflow to benchmark pg builders
##
###################################################################


#### Test cleaning

[ -d $1 ] && rm -r $1

#### Envs

ENV_MS="/home/genouest/genscale/sdubois/.conda/envs/mspangepop"
ENV_PGGB="/home/genouest/genscale/sdubois/.conda/envs/pggb_v0.7.4"
ENV_SAMTOOLS="/home/genouest/genscale/sdubois/.conda/envs/samtools"
ENV_CACTUS="apptainer run /projects/genscale/env/cactus_v2.9.9.sif"
ENV_VG="/projects/genscale/env/vg1.61.0"

#### Paths to files
FASTA=.fasta
PGGB_GFA=.pggb.gfa
MC_GFA=.mc.gfa
MS_GFA=.mspangepop.gfa
COMP_PGGB=.ms.pggb.tsv
COMP_MC=.ms.mc.tsv
PGGB_VCF=.pggb.vcf
MC_VCF=.mc.vcf
MS_VCF=.mspangepop.vcf

. /local/env/envconda.sh

######################## Simulating pangenome with MSpangepop ########################

conda activate $ENV_MS
cd MSpangepop/
#./mspangepop run --unlock
./mspangepop local-run 
cd ..
conda deactivate


######################## Index fasta file ########################
mkdir $1
conda activate $ENV_SAMTOOLS
for d in MSpangepop/results/*
do
    d="${d%/}"
    d="${d##*/}"
    mkdir $1/$d
    gzip -d MSpangepop/results/$d/03_graph/chr_1/fasta/*
    cp MSpangepop/results/$d/03_graph/chr_1/fasta/* $1/$d/multifasta$FASTA
    samtools faidx $1/$d/multifasta$FASTA # Index file for PGGB
    cp MSpangepop/results/$d/03_graph/chr_1/*.gfa $1/$d/graph$MS_GFA
done
conda deactivate

######################## Prepping MC files ########################

for d in $1/*
do
    mkdir $d/singlefasta # Creating a folder for MC
    awk -v output_dir="$d/singlefasta" '
    /^>/ {
        seq_id = substr($0, 2);
        if (index(seq_id, " ")) {
            seq_id = substr(seq_id, 1, index(seq_id, " ") - 1); 
        }
        gsub(/[^a-zA-Z0-9_-]/, "_", seq_id);
        file = sprintf("%s/%s.fasta", output_dir, seq_id);
        print $0 > file;
        next;
    }
    { print >> file; }
    ' "$d"/multifasta"$FASTA"
    mkdir $d/.pggb # Creating a folder for PGGB
    mkdir $d/.mc # Creating a folder for MC
done

######################## Creating the PGGB graph ########################

conda activate $ENV_PGGB
for d in $1/*
do
    pggb -i $d"/multifasta"$FASTA -o $d"/.pggb" -n 6 -t 8 -p 90 -s 5k 
    conda deactivate
    mv $d"/.pggb/"*.smooth.final.gfa $d"/graph"$PGGB_GFA
done
######################## Prep MC files ########################

for d in $1/*
do
REPATH=$(cat <<END
from os import listdir
with open("$d/.mc/pipeline.txt","w",encoding='utf-8') as writer:
    for seq in listdir("$d/singlefasta/"):
        with open(f"$d/singlefasta/{seq}") as f:
            header = f.readline().strip()[1:]
            sample,haplotype = header.split('#')[:-1]
        writer.write(f"{sample}.{haplotype}\t$d/singlefasta/{seq}\n")
END
)
FILE="$(python3 -c "$REPATH")"
sort -o $d/.mc/pipeline.txt{,}
done

######################## Creating the MC graph ########################

for d in $1/*
do
# Getting reference name (first line in file before \t)
echo "$(head -n 1 $d/.mc/pipeline.txt)" | cut -d$'\t' -f1 > $d/.mc/tempfile.txt
NAME_REF=`cat $d/.mc/tempfile.txt`
JB=$d/.mc/.js/
TMP_MC=$d/.mc/tmp/
OUT=$TMP_MC"graph"
[ -d $JB ] && rm -r $JB
mkdir $TMP_MC
$ENV_CACTUS cactus-minigraph $JB $d/.mc/pipeline.txt $OUT.gfa --reference $NAME_REF
$ENV_CACTUS cactus-graphmap $JB $d/.mc/pipeline.txt $OUT.gfa $OUT.paf  --reference $NAME_REF --outputFasta $OUT.sv.gfa.fa.gz
$ENV_CACTUS cactus-align $JB $d/.mc/pipeline.txt $OUT.paf $OUT.hal --pangenome --outGFA --outVG --reference $NAME_REF --workDir $TMP_MC
$ENV_CACTUS cactus-graphmap-join $JB --vg $OUT.vg --outDir $TMP_MC --outName "final" --reference $NAME_REF --clip 0 --filter 0
gzip -d $d"/.mc/tmp/final.full.gfa.gz"
done

######################## Convert and VCFs ########################

conda activate $ENV_VG
for d in $1/*
do
    NAME_REF=`cat $d/.mc/tempfile.txt`
    vg convert -g -f -W $d"/.mc/tmp/final.full.gfa" > $d"/graph"$MC_GFA # Get as GFA1.0 MC graph
    vg deconstruct -a $d/graph$MS_GFA -p $NAME_REF > $d/variants$MS_VCF
    vg deconstruct -a $d/graph$PGGB_GFA -p $NAME_REF > $d/variants$PGGB_VCF
    vg deconstruct -a $d/graph$MC_GFA -p $NAME_REF > $d/variants$MC_VCF
done
conda deactivate

######################## Compare the graphs ########################

if [ ! -f rs-pancat-compare ]; then
    wget https://github.com/dubssieg/rs-pancat-compare/releases/download/0.1.5/rs-pancat-compare
    chmod +x rs-pancat-compare
fi
for d in $1/*
do
    ./rs-pancat-compare $d/graph$MS_GFA $d/graph$PGGB_GFA > $d/dist$COMP_PGGB
    ./rs-pancat-compare $d/graph$MS_GFA $d/graph$MC_GFA > $d/dist$COMP_MC
    #[ -d $d"/.pggb" ] && rm -r $d"/.pggb"
    #[ -d $d"/.mc" ] && rm -r $d"/.mc"
done