import json
import MDAnalysis as mda
import sys, os


def writeTmp(AtomGroup):
	AtomGroup.write('tmp.gro')
	U = mda.Universe('tmp.gro', guess_bonds=True)
	os.remove('tmp.gro')
	return U.atoms

U = mda.Universe(sys.argv[1])

protein_grp = U.select_atoms('protein')
protein_grp = writeTmp(protein_grp)

water_grp = U.select_atoms('resname SOL')
water_grp = writeTmp(water_grp)

sodium_grp = U.select_atoms('resname NA')
chloride_grp = U.select_atoms('resname CL')

def mapData(AtomGroup, forcefield):

	if hasattr(AtomGroup.atoms, 'bonds'):
		bonds = AtomGroup.atoms.bonds.to_indices().tolist(),
	else:
		bonds = None

	return  {
		'types': AtomGroup.atoms.types.tolist(),
		'positions': AtomGroup.atoms.positions.tolist(),
		'masses': AtomGroup.atoms.masses.tolist(),
		'names': AtomGroup.atoms.names.tolist(),
		'bonds': bonds,
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
	'type' : ('periodic', 'periodic', 'periodic'),
	'bound': U.atoms.bbox().tolist()
}

data = {
    'molecule': [protein, water, sodium, chloride],

    'box': box,

    'solvent': {
	'implicit': True,
	'algorithm': 'gen-born',
	'radii_method': 'hct',
	'radii_freq': 1,
	'surf-tension': -1,
	'dielectric': 80,
	'salinity': 0.1
    }
}

with open('protein.json', 'w') as fp:
    json.dump(data, fp, sort_keys=False, indent=4, separators=(',', ': '))
