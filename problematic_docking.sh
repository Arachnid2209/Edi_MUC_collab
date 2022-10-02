### FULL AUTOMATION OF TF DOCKING ###

#	Defining some paths
DOCKING=/Users/maartenvandenancker/Desktop/docking_related/TF_docking
MGLTools=~/MGLTools/MGLToolsPckgs/AutoDockTools/Utilities24/

echo "What is your Transcription factor? (enter the filename without file extension)"
read PDB

echo "What is your Metal? (enter the filename without file extension)"
read METAL


#	This assumes you are starting in a directory where you want to work
#	This directory must have the TF you want to use as a PDB
#	For preparing structures from alphafold structures you need to delete the first and second 
#	last lines
#	This directory must also have the ligand file
#	Ligand file must be generated with beta carbon and gamma sulfur
#	i.e. SMILES would be C[S]-[METAL++] for a divalent metal
#	This directory must also have a modified AD4_parameters.dat file with the extra
#	Atom types you want to use

##CREATING LIST OF CYSTEINES AND MAKING A FOLDER FOR EACH##

grep -e "SG" ${PDB}.pdb | awk '{ print $5 $6 }' > cys_list.txt 
LIST=$(cat cys_list.txt)
for i in $LIST; do mkdir cys${i}; done
for i in $LIST/; do cp ${METAL}.pdb cys${i}/; cp ${PDB}.pdb cys${i}/; done

#	All the preparation files need python2 to run
#	py27 is a conda environment running python 2.7.18 and numpy 1.12.0
source ~/opt/miniconda3/etc/profile.d/conda.sh
conda activate py27
export PYTHONPATH="${PYTHONPATH}:/usr/local/lib/python2.7/site-packages:/usr/lib/python2.7/site-packages"

##PREPARING COVALENT STRUCTURES##
for i in $LIST; do cd cys${i}; chain=${i:0:1}; number=$(echo $i | sed "s/${chain}//g"); python ~/adcovalent/prepareCovalent.py --ligand ${METAL}.pdb --ligindices 1,2 --receptor ${PDB}.pdb --residue ${chain}:CYS${number} --outputfile ligcovalent.pdb; cd ..; done

##PREPARING THE RECEPTOR##
#	My python package installations are a bit messy, so I need to switch the PYTHONPATH variable
export PYTHONPATH=~/MGLTools/MGLToolsPckgs/:$PYTHONPATH
for i in $LIST; do cd cys${i}; $MGLTools/prepare_receptor4.py -r ${PDB}.pdb -A Hydrogens; cd ..; done
#	The ligand also uses prepare_receptor4.py for flexible docking / covalent docking
#	This is the problematic part of the problematic docking script


for i in $LIST; do cd cys${i}; cat ligcovalent.pdb | grep -v "0.000" > fixedligcovalent.pdb; $MGLTools/prepare_receptor4.py -r fixedligcovalent.pdb; cd ..; done
#	Gasteiger charges are stripped from Ag+ because prepare_receptor4.py doesnt recognise
#	the atom type, so I'll just replace the 0.000 charge with 1.000 (because Ag+)
for i in $LIST; do cd cys${i}; cat fixedligcovalent.pdbqt | sed 's/0.000/1.000/g' > gasteiger_lig.pdbqt; cd ..; done

##GENERATING FLEXIBLE PDBQT FILES##
#	Receptor:
for i in $LIST; do cd cys${i}; chain=${i:0:1}; number=$(echo $i | sed "s/${chain}//g"); $MGLTools/prepare_flexreceptor4.py -r ${PDB}.pdbqt -s ${PDB}:${chain}:CYS${number}; cd ..; done
#	Ligand:
for i in $LIST; do cd cys${i}; chain=${i:0:1}; number=$(echo $i | sed "s/${chain}//g"); $MGLTools/prepare_flexreceptor4.py -r gasteiger_lig.pdbqt -s gasteiger_lig:${chain}:CYS${number}; cd ..; done

##GENERATING PARAMETER FILES##
#	Grid parameter files (GPFs) - for autogrid
for i in $LIST; do cd cys${i}; $MGLTools/prepare_gpf4.py -r ${PDB}_rigid.pdbqt -x gasteiger_lig_flex.pdbqt -l gasteiger_lig_flex.pdbqt -y -I 20 -o ${PDB}.gpf; cd ..; done
#	Docking parameter files (DPFs) - for autodock
for i in $LIST; do cd cys${i}; touch empty; $MGLTools/prepare_dpf4.py -r ${PDB}_rigid.pdbqt -x gasteiger_lig_flex.pdbqt -l gasteiger_lig_flex.pdbqt -o dock_protein.dpf -p move='empty'; cd ..; done
#	Changing the unbound model parameter
for i in $LIST; do cd cys${i}; cat dock_protein.dpf | sed 's/unbound_model extended/unbound_energy 0.0/g' > new_dpf.dpf; cd ..; done
#	Telling autogrid and autodock to use new atom type parameter file
for i in $LIST; do cd cys${i}; echo -e "parameter_file AD4_parameters.dat\n$(cat new_dpf.dpf)" > new_dpf.dpf; echo -e "parameter_file AD4_parameters.dat\n$(cat ${PDB}.gpf)" > ${PDB}.gpf; cd ..; done

##AUTODOCK AND AUTOGRID##
#	For some reason autodock4 and autogrid need a python3 environment and all the 
#	programs required to prepare the files are python2. It makes no sense but it is what it is
conda deactivate
#	Autogrid:
for i in $LIST; do cd cys${i}; cp ../AD4_parameters.dat .; autogrid4 -p ${PDB}.gpf -l ${PDB}.glg; cd ..; done
#	Autodock:
for i in $LIST; do cd cys${i}; ~/autodock4 -p new_dpf.dpf -l cys${i}.dlg; cd ..; done

##THRESHOLD STUFF##
#	Maximum free energy is 0.60kcal/mol
#	This is because the free energy of rotation is always 0.60kcal/mol
#	This is because new atom types can not be used while calculating RMSD
#	Finding the free energy of binding
#		This removes the surrounding text, and if there is a plus it removes it

touch docked_coords.pdb
touch docked_resids.txt

for i in $(cat cys_list.txt); do gibbs=$(grep -e "Estimated Free Energy of Binding" cys${i}/cys${i}.dlg | head -n 1 | sed 's/DOCKED\: USER    Estimated Free Energy of Binding    \=   //g' | sed 's/kcal\/mol  \[\=(1)+(2)+(3)-(4)\]//g' | sed 's/+//g'); echo $gibbs; if (( $(echo "$gibbs < 0.60" | bc -l) )); then echo "$i passes" >> docked_resids.txt; grep -A2 'DOCKED: ATOM      3' cys${i}/cys${i}.dlg | head -2 | sed 's/ATOM      4/HETATM    4/g' | sed 's/DOCKED\: //g' >> docked_coords.pdb; else echo "$i fails" >> docked_resids.txt; fi; done

##STRUCTURE PREPARATION##
#	Make the model template
#	Getting the coords of all docked metals + gamma sulfur
#	Putting the coords into the model template and removing existing gamma sulfurs

echo -e "$(grep -v "SG" ${PDB}.pdb | sed '$d' | sed '$d')\n$(cat docked_coords.pdb)\n$(tail -n 2 ${PDB}.pdb)" > docked_model.pdb

for i in $(cat cys_list.txt); do gibbs=$(grep -e "Estimated Free Energy of Binding" cys${i}/cys${i}.dlg | head -n 1 | sed 's/DOCKED\: USER    Estimated Free Energy of Binding    \=   //g' | sed 's/kcal\/mol  \[\=(1)+(2)+(3)-(4)\]//g' | sed 's/+//g'); if (( $(echo "$gibbs < 0.60" | bc -l) )); then echo "$i $gibbs pass" >> stats.txt; else echo "$i $gibbs fail" >> stats.txt; fi; done

