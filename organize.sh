#!/usr/bin/env bash

# Author: Guilherme Duarte Ramos Matos
# Date: March 2022

###### Help Function ########
helpFunction(){
    echo -e "\tUsage: $0 -m molecules_csv"
    exit 1
}

# Assign typed arguments to variables
while getopts "m:" opt
do
    case $opt in
        m ) SOLUTE="$OPTARG";;
        ? ) helpFunction ;;
    esac
done

# Prints helpFuntion in case the number of parameters do not match what
# the script requires
if [ -z "${SOLUTE}" ]
then
    echo "You are misusing this script"
    helpFunction
fi

# define paths
root=$(pwd)
scriptdir=${root}/zzz.scripts

# Define relevant numbers for gromacs
temperature=298.15 #kelvin
pressure=1.01325 #bar


array=(${SOLUTE})
for val in "${array[@]}"
do

    # Create PDB files for each molecule
    echo "Creating PDB files from SMILES"
    pdbdir=${root}/001.molecules
    mkdir -p ${pdbdir}
    cd ${pdbdir}
    python3 ${scriptdir}/create_pdbfiles.py -c ${root}/$val
    # Remove lines containing 'CONECT'
    for elem in *.pdb
    do
        sed '/CONECT/d' $elem > tmp.pdb
        mv tmp.pdb $elem
    done

    echo "Parameterizing each generated molecule"
    for elem in *.pdb
    do

        echo ""
        echo ${elem}
        echo "All molecules are assumed to be neutral"
        molname=${elem%.*}
        # Define 3-letter code
        tmp_code=${molname::3}
        code=$(echo $tmp_code | tr '[:lower:]' '[:upper:]')
        echo "Starting Antechamber with ${molname}, code ${code}"
        # Fixing residue names
        antechamber -i ${molname}.pdb -fi pdb -o ${molname}_renamed.pdb -fo pdb -rn ${code}
        rm ${molname}.pdb
        mv ${molname}_renamed.pdb ${molname}.pdb
        # Creating simulation files
        antechamber -i ${molname}.pdb -fi pdb -o ${molname}.mol2 -fo mol2 -at gaff2 -c gas -rn ${code} -nc 0
        # Fix excess charges
        python3 ${scriptdir}/fix_charges.py -m ${molname}.mol2
        echo "Starting parmchk2"
        parmchk2 -i ${molname}.mol2 -f mol2 -o ${molname}.frcmod
        rm sqm.*

        # Transform Cl and Br in mol2 files in CL and BR.
        # tleap will fail if don't do it
        echo "Checking if there are halogens with inappropriate names"
        python3 ${scriptdir}/fixChlorineBromine.py -f ${molname}.mol2
        if [ -f "${molname}_fixed.mol2" ]
        then
            mv ${molname}.mol2 ${molname}.mol2.bckp
            mv ${molname}_fixed.mol2 ${molname}.mol2
        fi
    done
    cd ${root}
done


# For each template, generate files
echo "Generate simulation files"
for mol1 in ${root}/001.molecules/*.pdb
do
    # Create name of molecule 1
    tmp_mol1=${mol1%.pdb}
    name1=${tmp_mol1##*/}

    # Create directory and copy files
    mkdir -p ${pdbdir}/${name1}
    cd ${pdbdir}/${name1}
    mv ../${name1}.* .
    echo ${mol1}

    #Define 3-letter code for each molecule
    tmp_code1=${name1::3}
    code1=$(echo $tmp_code1 | tr '[:lower:]' '[:upper:]')

        # Create tleap input
    echo "Create tleap input"
    cat <<EOF > tl.in
source leaprc.gaff2

${code1} = loadmol2 ${name1}.mol2
loadamberparams ${name1}.frcmod

saveAmberParm ${code1} ${name1}.prmtop ${name1}.inpcrd
quit
EOF
    echo "Running tleap at $(pwd)"
    tleap -f tl.in

    # Create gromacs gro e top files
    python3 ${scriptdir}/amber_to_gmx.py -p ${name1}.prmtop -c ${name1}.inpcrd

    # top e gro foram feitos. Agora resta deixá-los apropriados para o restante do trabalho



done

cd ${root}
