#!/bin/bash
#SBATCH --job-name=bmMSpp
#SBATCH --constraint avx2
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G


###################################################################
## Author: Siegfried Dubois
##
## This script exectues a workflow to benchmark pg builders
##
###################################################################


#### Envs

ENV_MS="/home/genouest/genscale/sdubois/.conda/envs/mspangepop"
ENV_PGGB="/home/genouest/genscale/sdubois/.conda/envs/pggb_v0.7.4"
ENV_SAMTOOLS="/home/genouest/genscale/sdubois/.conda/envs/samtools"
ENV_CACTUS="apptainer run /projects/genscale/env/cactus_v2.9.9.sif"
ENV_VG="/projects/genscale/env/vg1.61.0"

#### Paths to files
FASTA=$1.fa
PGGB_GFA=$1.pggb.gfa
MC_GFA=$1.mc.gfa
MS_GFA=$1.mspangepop.gfa
COMP_PGGB=$1.ms.pggb.comp.tsv
COMP_MC=$1.ms.mc.comp.tsv

. /local/env/envconda.sh

######################## Simulating pangenome with MSpangepop ########################

conda activate $ENV_MS
cd MSpangepop/
#./mspangepop run --unlock
./mspangepop run 
cd ..
conda deactivate
mv "MSpangepop/results/test_panmictic/03_graph/chr_1/test_panmictic_chr_1_graph.gfa" $MS_GFA

######################## Extract the fasta from the simulated graph ########################

if [ ! -f rs-pancat-paths ]; then
    wget https://github.com/dubssieg/rs-pancat-paths/releases/download/0.1.2/rs-pancat-paths
    chmod +x rs-pancat-paths
fi
./rs-pancat-paths $MS_GFA reconstruct > $FASTA

######################## Index fasta file ########################

conda activate $ENV_SAMTOOLS
samtools faidx $FASTA
conda deactivate

######################## Creating the PGGB graph ########################
conda activate $ENV_PGGB
TMP_PGGB="."$1"pggb_output/"
mkdir $TMP_PGGB
pggb -i $FASTA -o $TMP_PGGB -n 6 -t 8 -p 90 -s 5k 
conda deactivate
mv $TMP_PGGB$FASTA.*.smooth.final.gfa $PGGB_GFA
[ -d $TMP_PGGB ] && rm -r $TMP_PGGB

######################## Prep MC files ########################

mkdir ".seq/"
awk 'BEGIN{RS=">";FS="\n"} NR>1{fnme=".seq/"$1".fa"; print ">" $0 > fnme; close(fnme);}' $FASTA
REPATH=$(cat <<END
from os import listdir
with open($1"_pipeline.txt","w",encoding='utf-8') as writer:
    for seq in listdir(".seq/"):
        writer.write(f"{seq[:-3]}\t.seq/{seq}")
END
)
FILE="$(python3 -c "$REPATH")"

######################## Creating the MC graph ########################

# Getting reference name (first line in file before \t)
echo "$(head -n 1 $1"_pipeline.txt")" | cut -d$'\t' -f1 > $1"_tempfile.txt"
NAME_REF=`cat $1"_tempfile.txt"`
JB=".js/"
TMP_MC="."$1"cactus_output/"
OUT=$TMP_MC"graph"
[ -d $JB ] && rm -r $JB
mkdir $JB $TMP_MC
$ENV_CACTUS cactus-minigraph $JB $1"_pipeline.txt" $OUT.gfa --reference $NAME_REF
$ENV_CACTUS cactus-graphmap $JB $1"_pipeline.txt" $OUT.gfa $OUT.paf  --reference $NAME_REF --outputFasta $OUT.sv.gfa.fa.gz
$ENV_CACTUS cactus-align $JB $1"_pipeline.txt" $OUT.paf $OUT.hal --pangenome --outGFA --outVG --reference $NAME_REF --workDir $TMP_MC
$ENV_CACTUS cactus-graphmap-join $JB --vg $OUT.vg --outDir $TMP_MC --outName "final" --reference $NAME_REF --clip 0 --filter 0
[ -d $JB ] && rm -r $JB
[ -f $1"_tempfile.txt" ] && rm $1"_tempfile.txt"
[ -f $1"_pipeline.txt" ] && rm $1"_pipeline.txt"
GRAPH=$TMP_MC"final.full.gfa"
gzip -d $GRAPH".gz"
conda activate $ENV_VG
vg convert -g -f -W $GRAPH > $MC_GFA
conda deactivate
[ -d $TMP_MC ] && rm -r $TMP_MC


######################## Compare the graphs ########################

if [ ! -f rs-pancat-compare ]; then
    wget https://github.com/dubssieg/rs-pancat-compare/releases/download/0.1.5/rs-pancat-compare
    chmod +x rs-pancat-compare
fi
./rs-pancat-compare $MS_GFA $PGGB_GFA > $COMP_PGGB
./rs-pancat-compare $MS_GFA $MC_GFA > $COMP_MC