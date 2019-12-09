Keywords
========
```
{
    "version": "...",

    "units": "nano",

    "molecule": {
        "types": "atomic element (e.g. C) or entity (e.g. CG) name",
        "positions": [
            [
                "x1",
                "y1",
                "z1"
            ],
            [
                "x2",
                "y2",
                "z2"
            ],
            [
                "..."
            ]
        ],
        "masses": [
            "mass1",
            "mass2",
            "..."
        ],
        "bonds (None)": [
            [
                "atom_index1",
                "atom_index2"
            ],
            [
                1,
                2
            ],
            [
                2,
                3
            ],
            [
                "..."
            ]
        ],
        "names (None)": [
		"list of atomic labels",
		"e.g. CA",
		"..."
	],
        "forcefield (None)": "forcefield name (e.g. charmm27)"
    },

    "box": {
        "shape": "cube",
        "type": [
            "periodic",
            "periodic",
            "periodic"
        ],
        "bound": [
            [
                "xmin",
                "ymin",
                "zmin"
            ],
            [
                "xmax",
                "ymax",
                "zmax"
            ]
        ]
    },

    "frac_coords (False)": "boolean variable specifying if coordinates are fractional",

    "solvent (None)": {

        "implicit (None)": {
            "surf_tension": "surface tension",
            "dielectric": "dielectric constant",
            "salinity": "concentration of (virtual) ions"
        },

        "explicit (None)": {
            "model": "water model e.g. spc216",
            "box (sim box)": [
                [
                    "xmin",
                    "ymin",
                    "zmin"
                ],
                [
                    "xmax",
                    "ymax",
                    "zmax"
                ]
            ],
            "shell (0)": "thickness of optional water layer around solute",
            "maxsol (inf)": "max number of water molecules to add if they fit in the box"
        }
    }
}
```
