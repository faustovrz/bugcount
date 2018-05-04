#!/bin/bash
# USAGE:
# analysis.sh wf_consolidated.tab

DATA=$1
#DATA='wf_consolidated.tab'

GLM='./analysis.R'
GLMPREFIX='resistance_analysis'

MM='./blup_heritability.R'
MMPREFIX='blup_analysis'

### Analysis by plant or leaf

for UNIT in plant leaf
do
nohup bash -c "${GLM} ${UNIT} ${GLMPREFIX} ${DATA}" \
        > ${UNIT}_${GLMPREFIX}.log 2>&1 &
nohup bash -c "${MM} ${UNIT} ${MMPREFIX} ${DATA}" \
        > ${UNIT}_${MMPREFIX}.log 2>&1 &
done
