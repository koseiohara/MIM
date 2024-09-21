#!/bin/bash

#PBS -q tqueue
#PBS -N MIM_MODERN
#PBS -j oe
#PBS -l nodes=1:ppn=1

START=$(date "+%s")

cd /mnt/jet11/kosei/mim/mim_modern/src/
ulimit -s unlimited

NOW=$(date "+%Y%m%d_%H%M%S")
DATA="JRA55"
INI=1990
FIN=2020
FLAG=""
RESULT="../output/${DATA}/result_${INI}_${FIN}_${NOW}${FLAG}.txt"

#NAMELIST="../nml/input_JRA55_2000_2000.nml"
NAMELIST="../nml/input_${DATA}_${INI}_${FIN}${FLAG}.nml"
#NAMELIST="../nml/input_${DATA}_test.nml"

cat $0 > ${RESULT}

echo -e "\n" >> ${RESULT}
echo "-----------" >> ${RESULT}
echo -e "\n\n" >> ${RESULT}

cat ${NAMELIST} >> ${RESULT}

echo -e "\n" >> ${RESULT}
echo "-----------" >> ${RESULT}
echo -e "\n\n" >> ${RESULT}


./MIM < ${NAMELIST} >> ${RESULT} 2>&1


END=$(date "+%s")


DIFF_SEC=$(expr ${END} - ${START})
DIFF_MIN=$(expr ${DIFF_SEC} / 60)
DIFF_SEC=$(expr ${DIFF_SEC} - ${DIFF_MIN} \* 60)
DIFF_HR=$(expr ${DIFF_MIN} / 60)
DIFF_MIN=$(expr ${DIFF_MIN} - ${DIFF_HR} \* 60)
echo " " >> ${RESULT}
echo "ELAPSE : ${DIFF_HR}hr ${DIFF_MIN}min ${DIFF_SEC}sec" >> ${RESULT}


