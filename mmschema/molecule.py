"""
Molecule Object Model
"""

class Molecule(ProtoModel):
    """
    A QCSchema representation of a Molecule. This model contains
    data for symbols, geometry, connectivity, charges, fragmentation, etc while also supporting a wide array of I/O and manipulation capabilities.

    Molecule objects geometry, masses, and charges are truncated to 8, 6, and 4 decimal places respectively to assist with duplicate detection.
    """

    schema_name: constr(strip_whitespace=True, regex=qcschema_molecule_default) = Field(  # type: ignore
        qcschema_molecule_default,
        description=(
            f"The QCSchema specification this model conforms to. Explicitly fixed as " f"{qcschema_molecule_default}."
        ),
    )
    schema_version: int = Field(  # type: ignore
        2, description="The version number of ``schema_name`` that this Molecule model conforms to."
    )
    validated: bool = Field(  # type: ignore
        False,
        description="A boolean indicator (for speed purposes) that the input Molecule data has been previously checked "
        "for schema (data layout and type) and physics (e.g., non-overlapping atoms, feasible "
        "multiplicity) compliance. This should be False in most cases. A ``True`` setting "
        "should only ever be set by the constructor for this class itself or other trusted sources such as "
        "a Fractal Server or previously serialized Molecules.",
    )

    # Required data
    symbols: Array[str] = Field(  # type: ignore
        ...,
        description="An ordered (nat,) array-like object of atomic elemental symbols of shape (nat,). The index of "
        "this attribute sets atomic order for all other per-atom setting like ``real`` and the first "
        "dimension of ``geometry``. Ghost/Virtual atoms must have an entry in this array-like and are "
        "indicated by the matching the 0-indexed indices in ``real`` field.",
    )
    geometry: Array[float] = Field(  # type: ignore
        ...,
        description="An ordered (nat,3) array-like for XYZ atomic coordinates [a0]. "
        "Atom ordering is fixed; that is, a consumer who shuffles atoms must not reattach the input "
        "(pre-shuffling) molecule schema instance to any output (post-shuffling) per-atom results "
        "(e.g., gradient). Index of the first dimension matches the 0-indexed indices of all other "
        "per-atom settings like ``symbols`` and ``real``."
        "\n"
        "Can also accept array-likes which can be mapped to (nat,3) such as a 1-D list of length 3*nat, "
        "or the serialized version of the array in (3*nat,) shape; all forms will be reshaped to "
        "(nat,3) for this attribute.",
    )

    # Molecule data
    name: Optional[str] = Field(  # type: ignore
        None, description="A common or human-readable name to assign to this molecule. Can be arbitrary."
    )
    identifiers: Optional[Identifiers] = Field(  # type: ignore
        None,
        description="An optional dictionary of additional identifiers by which this Molecule can be referenced, "
        "such as INCHI, canonical SMILES, etc. See the :class:``Identifiers`` model for more details.",
    )
    comment: Optional[str] = Field(  # type: ignore
        None,
        description="Additional comments for this Molecule. Intended for pure human/user consumption " "and clarity.",
    )
    molecular_charge: float = Field(0.0, description="The net electrostatic charge of this Molecule.")  # type: ignore
    molecular_multiplicity: int = Field(1, description="The total multiplicity of this Molecule.")  # type: ignore

    # Atom data
    masses: Optional[Array[float]] = Field(  # type: ignore
        None,
        description="An ordered 1-D array-like object of atomic masses [u] of shape (nat,). Index order "
        "matches the 0-indexed indices of all other per-atom settings like ``symbols`` and ``real``. If "
        "this is not provided, the mass of each atom is inferred from their most common isotope. If this "
        "is provided, it must be the same length as ``symbols`` but can accept ``None`` entries for "
        "standard masses to infer from the same index in the ``symbols`` field.",
    )
    real: Optional[Array[bool]] = Field(  # type: ignore
        None,
        description="An ordered 1-D array-like object of shape (nat,) indicating if each atom is real (``True``) or "
        "ghost/virtual (``False``). Index "
        "matches the 0-indexed indices of all other per-atom settings like ``symbols`` and the first "
        "dimension of ``geometry``. If this is not provided, all atoms are assumed to be real (``True``)."
        "If this is provided, the reality or ghostality of every atom must be specified.",
    )
    atom_labels: Optional[Array[str]] = Field(  # type: ignore
        None,
        description="Additional per-atom labels as a 1-D array-like of of strings of shape (nat,). Typical use is in "
        "model conversions, such as Elemental <-> Molpro and not typically something which should be user "
        "assigned. See the ``comments`` field for general human-consumable text to affix to the Molecule.",
    )
    atomic_numbers: Optional[Array[np.int16]] = Field(  # type: ignore
        None,
        description="An optional ordered 1-D array-like object of atomic numbers of shape (nat,). Index "
        "matches the 0-indexed indices of all other per-atom settings like ``symbols`` and ``real``. "
        "Values are inferred from the ``symbols`` list if not explicitly set.",
    )
    mass_numbers: Optional[Array[np.int16]] = Field(  # type: ignore
        None,
        description="An optional ordered 1-D array-like object of atomic *mass* numbers of shape (nat). Index "
        "matches the 0-indexed indices of all other per-atom settings like ``symbols`` and ``real``. "
        "Values are inferred from the most common isotopes of the ``symbols`` list if not explicitly set.",
    )

    # Fragment and connection data
    connectivity: Optional[List[Tuple[int, int, float]]] = Field(  # type: ignore
        None,
        description="The connectivity information between each atom in the ``symbols`` array. Each entry in this "
        "list is a Tuple of ``(atom_index_A, atom_index_B, bond_order)`` where the ``atom_index`` "
        "matches the 0-indexed indices of all other per-atom settings like ``symbols`` and ``real``.",
    )
    fragments: Optional[List[Array[np.int32]]] = Field(  # type: ignore
        None,
        description="An indication of which sets of atoms are fragments within the Molecule. This is a list of shape "
        "(nfr) of 1-D array-like objects of arbitrary length. Each entry in the list indicates a new "
        "fragment. The index "
        "of the list matches the 0-indexed indices of ``fragment_charges`` and "
        "``fragment_multiplicities``. The 1-D array-like objects are sets of atom indices indicating the "
        "atoms which compose the fragment. The atom indices match the 0-indexed indices of all other "
        "per-atom settings like ``symbols`` and ``real``.",
    )
    fragment_charges: Optional[List[float]] = Field(  # type: ignore
        None,
        description="The total charge of each fragment in the ``fragments`` list of shape (nfr,). The index of this "
        "list matches the 0-index indices of ``fragment`` list. Will be filled in based on a set of rules "
        "if not provided (and ``fragments`` are specified).",
    )
    fragment_multiplicities: Optional[List[int]] = Field(  # type: ignore
        None,
        description="The multiplicity of each fragment in the ``fragments`` list of shape (nfr,). The index of this "
        "list matches the 0-index indices of ``fragment`` list. Will be filled in based on a set of "
        "rules if not provided (and ``fragments`` are specified).",
    )

    # Orientation
    fix_com: bool = Field(  # type: ignore
        False,
        description="An indicator which prevents pre-processing the Molecule object to translate the Center-of-Mass "
        "to (0,0,0) in euclidean coordinate space. Will result in a different ``geometry`` than the "
        "one provided if False.",
    )
    fix_orientation: bool = Field(  # type: ignore
        False,
        description="An indicator which prevents pre-processes the Molecule object to orient via the inertia tensor."
        "Will result in a different ``geometry`` than the one provided if False.",
    )
    fix_symmetry: Optional[str] = Field(  # type: ignore
        None, description="Maximal point group symmetry which ``geometry`` should be treated. Lowercase."
    )
    # Extra
    provenance: Provenance = Field(  # type: ignore
        provenance_stamp(__name__),
        description="The provenance information about how this Molecule (and its attributes) were generated, "
        "provided, and manipulated.",
    )
    id: Optional[Any] = Field(  # type: ignore
        None,
        description="A unique identifier for this Molecule object. This field exists primarily for Databases "
        "(e.g. Fractal's Server) to track and lookup this specific object and should virtually "
        "never need to be manually set.",
    )
    extras: Dict[str, Any] = Field(  # type: ignore
        None, description="Extra information to associate with this Molecule."
    )
   
   # Lattice 
    fractional_coords: Optional[bool] = Field(  # type: ignore
        False,
        description="boolean variable specifying if coordinates are fractional e.g. for unit cells",
    )

    cell_number: Optional[Tuple[int, int, int, int int int]] = Field( # type: ignore
        (1,1,1),
        description="The first three values specify dimensions in number of cells along the unit cell axes; the second three values specify the placement of the block relative to the cell containing the original structure."

    cell_origin: Optional[Tuple[float, float, float]] = Field( # type: ignore
        (0,0,0),
        description="Coordinates in unit cell lengths and describe translations of the unit cell along its axes."
    )

    cell_lengths: Optional[Tuple[float, float, float]] = Field( # type: ignore
        (1,1,1),
        description="Lengths of the unit cell in each direction"
    )

    cell_angles: Optional[Tuple[float, float, float]] = Field( # type: ignore
        (90, 90, 90),
        description="Angles of the unit cell vectors. By default, the cell vectors are assumed to be orthogonal."
    )

    cell_offset: Optional[List[Tuple[int, int, int, int]]] = Field(  # type: ignore
        None,
        description="List of cell atom offsets for each axis: (atom_index, int_x, int_y, int_z)"
    )
