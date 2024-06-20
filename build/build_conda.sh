#!/usr/bin/env bash
set -e

eval "$(conda shell.bash hook)"

TAG_OR_BRANCH="${TAG_OR_BRANCH:-develop}"

conda create --name fluffy_${TAG_OR_BRANCH} python=3.9 -y

conda activate fluffy_${TAG_OR_BRANCH}

conda install -c conda-forge pip -y

if [ -d fluffy_${TAG_OR_BRANCH} ];
then
    rm -fr fluffy_${TAG_OR_BRANCH}
fi

mkdir fluffy_${TAG_OR_BRANCH}

git clone --branch ${TAG_OR_BRANCH} https://github.com/clinical-genomics-uppsala/fluffy_hematology_wgs fluffy_${TAG_OR_BRANCH}/fluffy

#pip install -r fluffy_${TAG_OR_BRANCH}/fluffy/requirements.txt 
pip install -r /home/jonas/Snakemake/fluffy_hematology_wgs/requirements.txt

conda deactivate

conda pack -n fluffy_${TAG_OR_BRANCH} -o fluffy_${TAG_OR_BRANCH}/env.tar.gz

mkdir -p fluffy_${TAG_OR_BRANCH}/hydra-genetics

git clone https://github.com/snakemake/snakemake-wrappers.git fluffy_${TAG_OR_BRANCH}/snakemake-wrappers

git clone https://github.com/hydra-genetics/prealignment.git fluffy_${TAG_OR_BRANCH}/hydra-genetics/prealignment
git clone https://github.com/hydra-genetics/alignment.git fluffy_${TAG_OR_BRANCH}/hydra-genetics/alignment
git clone https://github.com/hydra-genetics/snv_indels.git fluffy_${TAG_OR_BRANCH}/hydra-genetics/snv_indels
git clone https://github.com/hydra-genetics/annotation.git fluffy_${TAG_OR_BRANCH}/hydra-genetics/annotation
git clone https://github.com/hydra-genetics/filtering.git fluffy_${TAG_OR_BRANCH}/hydra-genetics/filtering
git clone https://github.com/hydra-genetics/qc.git fluffy_${TAG_OR_BRANCH}/hydra-genetics/qc
git clone https://github.com/hydra-genetics/biomarker.git fluffy_${TAG_OR_BRANCH}/hydra-genetics/biomarker
git clone https://github.com/hydra-genetics/fusions.git fluffy_${TAG_OR_BRANCH}/hydra-genetics/fusions
git clone https://github.com/hydra-genetics/cnv_sv.git fluffy_${TAG_OR_BRANCH}/hydra-genetics/cnv_sv
git clone https://github.com/hydra-genetics/compression.git fluffy_${TAG_OR_BRANCH}/hydra-genetics/compression
git clone https://github.com/hydra-genetics/misc.git fluffy_${TAG_OR_BRANCH}/hydra-genetics/misc
git clone https://github.com/hydra-genetics/reports.git fluffy_${TAG_OR_BRANCH}/hydra-genetics/reports
git clone https://github.com/hydra-genetics/parabricks.git fluffy_${TAG_OR_BRANCH}/hydra-genetics/parabricks
git clone https://github.com/hydra-genetics/sentieon.git fluffy_${TAG_OR_BRANCH}/hydra-genetics/sentieon

tar -zcvf fluffy_${TAG_OR_BRANCH}.tar.gz fluffy_${TAG_OR_BRANCH}

if [ -d fluffy_${TAG_OR_BRANCH} ];
then
    rm -fr fluffy_${TAG_OR_BRANCH}
fi
