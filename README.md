```js
st140216288280784=>start: start compute
io140216288280912=>inputoutput: input: cls, input_data, config
cond140216288281168=>condition: if isinstance(input_data, cls.input())
sub140216288281616=>subroutine: cls.input()(**input_data.dict())
op140216288281424=>operation: program = cls(name=cls.__name__)
op140216288282128=>operation: (_, exec_output) = program.execute(input_data)
cond140216288282576=>condition: if isinstance(exec_output, cls.output())
sub140216288282000=>subroutine: cls.output()(**exec_output.dict())
io140216288281552=>inputoutput: output:  exec_output
e140216288282832=>end: end function return
op140216288282256=>operation: raise TypeError
op140216288281808=>operation: raise TypeError

st140216288280784->io140216288280912
io140216288280912->cond140216288281168
cond140216288281168(yes)->sub140216288281616
sub140216288281616->op140216288281424
op140216288281424->op140216288282128
op140216288282128->cond140216288282576
cond140216288282576(yes)->sub140216288282000
sub140216288282000->io140216288281552
io140216288281552->e140216288282832
cond140216288282576(no)->op140216288282256
op140216288282256->io140216288281552
cond140216288281168(no)->op140216288281808
op140216288281808->op140216288281424
```

```bash
python -m pyflowchart mmic/components/base/base_component.py -f ProgramHarness.compute
npm install diagrams
diagrams flowchart input flowchart.svg
```

Molecular Mechanics Schema
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
