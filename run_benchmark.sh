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
COMP_PGGB=.ms.pggb.comp.tsv
COMP_MC=.ms.mc.comp.tsv
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
    mv MSpangepop/results/$d/03_graph/chr_1/fasta/* $1/$d/multifasta$FASTA
    samtools faidx $1/$d/multifasta$FASTA # Index file for PGGB
    mv MSpangepop/results/$d/03_graph/chr_1/*.gfa $1/$d/graph$MS_GFA
done
conda deactivate

######################## Prepping MC files ########################

for d in $1/*
do
    mkdir $d/singlefasta # Creating a folder for MC
    awk -v output_dir="${$d"/singlefasta"}" '
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
    ' "${$d"multifasta"$FASTA}"
    mkdir $d/.pggb # Creating a folder for PGGB
    mkdir $d/.mc # Creating a folder for MC
done

exit 1

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
PPATH=$1"_pipeline.txt"
awk 'BEGIN{RS=">";FS="\n"} NR>1{fnme=".seq/"$1".fa"; print ">" $0 > fnme; close(fnme);}' $FASTA
REPATH=$(cat <<END
from os import listdir
with open("$PPATH","w",encoding='utf-8') as writer:
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
[ -d $TMP_MC ] && rm -r $TMP_MC

######################## Convert and VCFs ########################

conda activate $ENV_VG
vg convert -g -f -W $GRAPH > $MC_GFA # Get as GFA1.0 MC graph
vg deconstruct -a $MS_GFA -p $NAME_REF > $MS_VCF
vg deconstruct -a $PGGB_GFA -p $NAME_REF > $PGGB_VCF
vg deconstruct -a $MC_GFA -p $NAME_REF > $MC_VCF
conda deactivate

######################## Compare the graphs ########################

if [ ! -f rs-pancat-compare ]; then
    wget https://github.com/dubssieg/rs-pancat-compare/releases/download/0.1.5/rs-pancat-compare
    chmod +x rs-pancat-compare
fi
./rs-pancat-compare $MS_GFA $PGGB_GFA > $COMP_PGGB
./rs-pancat-compare $MS_GFA $MC_GFA > $COMP_MC