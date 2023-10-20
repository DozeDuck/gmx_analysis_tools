#!/bin/bash
# Usage: ./renumber_models.sh extream1.pdb

###
### gro_merge.sh - merge rec.gro and lig_GMX.gro, also edit topol.top as well
###
### Usage:
###       gro_merge.sh <receptor.gro>  <ligand_GMX.gro>  <ligand_GMX.itp>
###
### Options:
###     <receptor.gro>		The gro file generated from gmx pdb2gmx
###     <ligand_GMX.gro>		The gro file generated from acpype
###     <ligand_GMX.itp>		The itp file generated fro acpype
###     -h				Show this message

help() {
	sed -rn 's/^### ?//;T;p;' "$0"
}

if [[ $# == 0 ]]  || [[ "$1" == "-h" ]]; then
	help
	exit 1
fi

# input
input=$1

# 初始化计数器
count=1

# 循环遍历文件中的每一行
while IFS= read -r line; do
    # 如果行包含"MODEL 1"，则替换为"MODEL 1"后面接计数器的值，并将计数器加一
    if [[ $line == *MODEL\ \ \ \ \ \ \ \ 1* ]]; then
        echo "$line" | sed "s/MODEL        1/MODEL        $count/"
        ((count++))
    else
        echo "$line"
    fi
done < $input > renumbered_$input
