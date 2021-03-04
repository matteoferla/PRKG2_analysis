# PRKG2_analysis
_In silico_ analysis of the cGMP-dependent protein kinase 2 (PRKG2 encoded).

This protein has two crystallised domains, both binding a cyclic nucleotides.
There is a full-length structure of a 30% identity homologue (PDB:5DYK),
namely the cGMP-dependent protein kinase PKG from Plasmodium falciparum. This domain binds ATP and the target peptide.

## Model

Two ligand bound structures of the N-terminal domains were minimised against their electron density with bound ligands.

* CNB-A domain — PDB:5C6C — This has cAMP (`CMP`). Whereas the calcium ion was kept (two contact points), ethylene glycol, cadmium ions and cobalt ions were removed.
* CNB-B domain — PDB:5BV6 — This has cGMP (`35G`). A sodium ion was kept, but calcium and acetate ions were removed.

For the code used see [notes about code - crystal](code_crystal.md)

The full sequence was submitted to I-Tasser. Model 1 and two ligands were chosen:

* peptide substrate (PDB:4Z83)
* ATP (PDB:1ATP)

The three parts were combined manually, see [notes about code - collage](code_collage.md).
The loops were corrected, after a brief energy minimisation.
Afterwards, the model was energy minimised fully.

The multipart approach is totally overkill as the mutations were not in the crystal structure parts.

For C-terminal extension see [notes about code - extension](code_extension.md)