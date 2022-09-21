from vina import Vina
from meeko import MoleculePreparation
from rdkit import Chem
from pymol import cmd


def create_pdbqt_file_from_pdb(pdb_file, hydrate, output_name):
    mol = Chem.MolFromPDBFile(pdb_file)
    preparator = MoleculePreparation(hydrate=hydrate)
    preparator.prepare(mol)
    preparator.write_pdbqt_file(f'{output_name}.pdbqt')


def create_pdbqt_file_from_mol2(mol2_file, hydrate, output_name):
    mol = Chem.MolFromMol2File(mol2_file)
    preparator = MoleculePreparation(hydrate=hydrate)
    preparator.prepare(mol)
    preparator.write_pdbqt_file(f'{output_name}.pdbqt')


def save_dna_pdb_from_seq(dna_seq, output_name):
    cmd.fnab(dna_seq, name="dna_1", mode="DNA", form="B", dbl_helix=1)
    pdb = cmd.get_pdbstr("dna_1")
    f = open(f"{output_name}.pdb", "w")
    f.write(pdb)
    f.close()

#creates pdbqt file from pdb file of protein ligand
def protein_ligand(pdb_file_name):
    f = open(pdb_file_name, "r")
    pdb_str = f.read()
    cmd.read_pdbstr(pdb_str, oname='protein_ligand')
    cmd.h_add(selection="protein_ligand")
    mol_str = cmd.get_str(format='mol', selection="protein_ligand")
    f = open("prot_ligand.mol", "w")
    f.write(mol_str)
    f.close()
    out_dna_2 = Chem.MolFromMolFile('prot_ligand.mol')
    preparator = MoleculePreparation(hydrate=True)
    preparator.prepare(out_dna_2)
    preparator.write_pdbqt_file("prot_ligand.pdbqt")

def prepare_dna(dna_seq, output_name):
    cmd.fnab(dna_seq, name="dna_2", mode="DNA", form="B", dbl_helix=1)
    cmd.h_add(selection="dna_2")
    mol_str = cmd.get_str(format='mol', selection="dna_2")
    f = open("dna_2.mol", "w")
    f.write(mol_str)
    f.close()
    out_dna_2 = Chem.MolFromMolFile('dna_2.mol')
    preparator = MoleculePreparation(hydrate=True)
    preparator.prepare(out_dna_2)
    preparator.write_pdbqt_file("dna_2.pdbqt")


def mutate_protein(protein_pdb, position, new_aa, output_name):
    cmd.load(protein_pdb)
    cmd.wizard("mutagenesis")
    cmd.do("refresh_wizard")
    cmd.get_wizard().set_mode(new_aa)
    cmd.get_wizard().do_select(f"{position}/")
    cmd.frame(11)
    cmd.get_wizard().apply()
    pdb = cmd.get_pdbstr(protein_pdb[:-4])
    f = open(f"{output_name}.pdb", "w")
    f.write(pdb)
    f.close()


print("Started calculation!")
create_pdbqt_file_from_mol2('tetracycline.mol2', True, 'tetracycline')
v = Vina(sf_name='vina')
v.set_receptor('tetR_ADFR.pdbqt')
print("Receptor set!")
v.set_ligand_from_file('tetracycline.pdbqt')
print("Ligand set!")
v.compute_vina_maps(center=[15.917, -2.405, -1.171], box_size=[40, 40, 40])  # 40.793, 17.327, -12.182
print("Computed Vina maps!")

# Score the current pose
energy = v.score()
print('Score before minimization: %.3f (kcal/mol)' % energy[0])

# Minimized locally the current pose
energy_minimized = v.optimize()
print('Score after minimization : %.3f (kcal/mol)' % energy_minimized[0])
v.write_pose('prot_tetracycl.pdbqt', overwrite=True)

# Dock the ligand
v.dock(exhaustiveness=32, n_poses=20)
v.write_poses('prot_tetracyc_vina_out.pdbqt', n_poses=5, overwrite=True)
