#!/bin/bash

p=$(realpath -- ./fasta)
m=100

if [[ ! -z "${1}" ]]; then
  m=${1}
fi

c=0
rm fof_${m}_kmtricks.txt &> /dev/null
rm fof_${m}.txt &> /dev/null
while read -r exp;
do
  IFS=" " read name cutoff <<< "$exp"
  if [[ -f "${p}/${name}.fasta.gz" ]]; then
    echo "Error: ${p}/${name}.fasta.gz does not exists !"
    break
  fi
  
  echo "D${c} : ${p}/${name}.fasta.gz ! ${cutoff}" >> fof_${m}_kmtricks.txt
  echo "${p}/${name}.fasta.gz ${cutoff}" >> fof_${m}.txt
  let "c++"
  
  if [[ ${c} == ${m} ]]; then
    break 
  fi
done < sra_id_threshold.txt