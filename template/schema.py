import json
import MDAnalysis as mda
import sys, os


def writeTmp(AtomGroup):
	AtomGroup.write('tmp.gro')
	U = mda.Universe('tmp.gro', guess_bonds=True)
	os.remove('tmp.gro')
	return U.atoms

U = mda.Universe(sys.argv[1])

protein_grp = U.select_atoms('resname ALA')
protein_grp = writeTmp(protein_grp)

water_grp = U.select_atoms('resname SOL').select_atoms('resnum 10:12')
water_grp = writeTmp(water_grp)

sodium_grp = U.select_atoms('resname NA')
chloride_grp = U.select_atoms('resname CL')

def mapData(AtomGroup, forcefield):

	if hasattr(AtomGroup.atoms, 'bonds'):
		top = {
			'bonds': AtomGroup.atoms.bonds.to_indices().tolist(),
			'angles': AtomGroup.atoms.angles.to_indices().tolist(),
			'dihedrals': AtomGroup.atoms.dihedrals.to_indices().tolist()
		}
	else:
		top = None

	return  {
		'positions': AtomGroup.atoms.positions.tolist(),
		'masses': AtomGroup.atoms.masses.tolist(),
		'names': AtomGroup.atoms.names.tolist(),
		'topology': top,
		'forcefield': forcefield
	}

# molecule example
protein = mapData(protein_grp, 'AMBER99')
water = mapData(water_grp, 'TIP3P')
sodium = mapData(sodium_grp, 'AMBER99')
chloride = mapData(chloride_grp, 'AMBER99')

# simulation box
box = {
	'shape': 'cube',
	'bound': U.atoms.bbox().tolist()
}

data = {
    'molecule': [protein, water, sodium, chloride],

    'box': box,

    'props': {
    	'pH': 2,
    	'solvent_implicit': None,
    }
}

with open('bio.json', 'w') as fp:
    json.dump(data, fp, sort_keys=True, indent=4, separators=(',', ': '))
