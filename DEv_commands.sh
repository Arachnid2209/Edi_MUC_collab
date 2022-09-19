##MAKING FOLDERS##
mkdir -p DEv/pdb DEv/no_metal_pdb DEv/obminimized_pdbs DEv/prepared_MTs DEv/autodock_files DEv/structures
DEv=$(pwd)/DEv
cd $DEv/pdb

##DOWNLOADING PDBs OF KNOWN METALLOTHIONEIN STRUCTURES##
#	Now you have to make a comma separated list of PDB IDs
#	Mine (batch.txt) looks like this
#	1DMC,1DMD,1DME,1DMF,6GV6,6GV7,6GV8,6GV9,1FMY,1AQQ,1AQR,1AQS,1AOO

curl -O https://www.rcsb.org/scripts/batch_download.sh
batch_download.sh -f batch.txt -p
gunzip *

##REMOVING ANY METALS IN THE STRUCTURES##
#	Non-typical atom types in PDB files start with the line HETATM
#	This for loop takes all lines of the various PDB files except those starting with HETATM
#	and makes a new PDB file with the same name in the $DEv/no_metal_pdb folder

for i in *.pdb; do grep -v "^HETATM" $i > $DEv/no_metal_pdb/$i; done

##ENERGY MINIMIZATION##
#	This energy minimization uses openbabel
#	The forcefield model we are using is the General Amber Force Field (GAFF)
#	The number of iterations has been reduced to 300 from the usual 2500 (for time reasons)
#	The output of the obminimize function has things we are not interested in
#	Between the lines starting with AUTHOR and MASTER is the pdb file
#	We save this sed output as a new energy minimized pdb file

cd $DEv/no_metal_pdb
for i in *.pdb; do obminimize -ff GAFF -n 300 $i | sed -n '/^AUTHOR/,/^MASTER/p' > $DEv/obminimized_pdbs/$i; done

##PREPARING METALLOTHIONEINS FOR AUTODOCK##
#	For this you need to be in a conda environment running some version of python 2
#	I am using a conda environment with Python 2.7.18, numpy 1.16.6
#	You can create this with

conda create -n "py27" python=2.7.18
conda activate py27
conda install -c conda-forge numpy=1.16.6

#	However, the MGLTools python scripts will still not work
#	You will get an issue with the module MolKit - "ImportError: No module named MolKit"
#	cd into the directory where the MGLTools folder is
#	Type the following command

export PYTHONPATH=$(pwd)/MGLTools/MGLToolsPckgs/:$PYTHONPATH

#	Everything should work with MGLTools now

MGLTools=$(pwd)/MGLTools/MGLToolsPckgs/AutoDockTools/Utilities24/

##PREPARING THE METALLOTHIONEIN PDBs FOR AUTODOCK##
#	prepare_receptor4.py generates the .pdbqt file that AutoDock4 uses
#	The -A hydrogens flag repairs this molecule if it is missing hydrogens 

cd $DEv/obminimized_pdbs
for i in *.pdb; do $MGLTools/prepare_receptor4.py -r $i -o $DEv/prepared_MTs/$i -A hydrogens; done

##PREPARING THE METAL ION LIGAND FOR AUTODOCK##



##PREPARING GRID PARAMETER FILES AND DOCKING PARAMETER FILES##

cd $DEv/autodock_files
for i in $DEv/prepared_MTs/*.pdbqt; do $MGLTools/prepare_gpf4.py -l metal.pdbqt -r $i -y; done
for i in $DEv/prepared_MTs/*.pdbqt; do $MGLTools/prepare_dpf4.py -l metal.pdbqt -r $i; done
for i in *.gpf; do echo -e "parameter_file AD4_parameters.dat\n$(cat ${i})" > $i; done
for i in *.dpf; do echo -e "parameter_file AD4_parameters.dat\n$(cat ${i})" > $i; done

##AUTOGRID AND AUTODOCK##
#	AutoDock4 and AutoGrid use python3 instead of python2
#	I have no idea why the software to prepare the files is in python2
#	AutoGrid will generate a bunch of files mapping the noncovalent interactions between the
#	rigid receptor and probe atom
#	Produces a ton of mapping files, and a .glg file (grid log)
#	AutoDock does the actual docking, and produces a .dlg file (docking log)
#	You might need to do output redirection for both .dlg and .glg

conda deactivate
for i in $DEv/autodock_files/*.gpf; do autogrid4 -p $i; done
for i in $DEv/autodock_files/*.dpf; do autodock4 -p $i; done

##TURNING AUTODOCK OUTPUT INTO A PDB FILE##
#	This code pretty much just finds the coordinates for the ligand and adds it on to 
#	the original receptor file.
#	The .dlg file generated earlier has the kinetic data

for i in $DEv/autodock_files/*.dlg;\
do \
filename=$(basename $i .dlg);\
grep "^DOCKED" $i | cut -c9- | cut -c-66 > $DEv/structures/${filename}_docked_ligand.pdb;\
cat receptor.pdb $DEv/structures/${filename}_docked_ligand.pdb | grep -v '^END   ' | grep -v '^END$' > $DEv/structures/${filename}_complex.pdb;\
done


