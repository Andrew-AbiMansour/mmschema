MOLSSI Molecular Mechanics Schema
=================================
MMSchema seeks to provide a standard definition of molecular mechanics (MM) simulation input/output data. The goal behind  MMSchema is to provide interoperable data structures that could be used in existing workflow systems or MM applications.

The general theme of this project is to explore sustainable and practical ways by which classical molecular mechanics (MM) users can exchange MM simulation data in a manner that is agnostic to the MM application. The MM data we’re interested in includes “core” input simulation parameters and system definition (molecule, solvent, etc.).

MMSchema is divided into 3 classes:
* [System schema](SystemSchema.md): a general definition of an MM system
* [Simulation schema](SimSchema.md): a set of common input MM parameters that define a simulation

Examples
========
* [Protein](data/alanine/protein.md)
* [Protein in aqueous soln](data/alanine/protein_aq.md)
* [Periodic Polymer](data/polyeth/polymer.md)
