#!/bin/bash
# Usage: ./gro_merge.sh rec.gro lig_GMX.gro lig_GMX.itp

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

# This part is for merge two gro files
more  $2 | sed -n '3,$p' > test.gro && sed -i '' test.gro 
# a=`sed -n 2p rec.gro`
a=`sed -n 2p $1`
b=`sed -n 2p $2 `
c=$(($a+$b))
# echo $c
more $1 > complex.gro
sed -i "2c $c" complex.gro
tail -1 $1 > last
sed -i '$d' complex.gro
more test.gro >> complex.gro
more last >> complex.gro
rm last 
rm test.gro

# Below is the part for editing topol.top
echo "; Include ligand topology" >> include.dat
# echo "#include \"lig_GMX.itp\"" >> include.dat
echo "#include $3" >> include.dat
echo "#ifdef POSRES_LIG" >> include.dat
echo "#include \"posre_lig.itp\"" >> include.dat
echo "#endif

" >>  include.dat

a=`more $3 | grep -A1 nrexcl | grep -v nrexcl`
B=`echo $a | awk '{print $1}'`
# sed -i '$a lig\ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ 1' topol.top 
echo "lable"
sed -i '$a lolololololololol' topol.top 
sed -i "s/lolololololololol/$B              1/g" topol.top
sed -i '22r include.dat'     topol.top
rm include.dat

