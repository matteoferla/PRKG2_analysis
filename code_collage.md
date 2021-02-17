## Code: merged

> This is a copy paste from a Jupyter notebook. Some functions are in
> [my pyrosetta_scripts repo](https://github.com/matteoferla/pyrosetta_scripts)

The I-Tasser model, the Swissmodel and the two domain were combined manually.

### I-Tasser

The models were similar. Model 1, C-score: -1.10, estimated RMSD = 10.9±4.6Å.
The difference is simply from a loop, which was less spaghetti from Swissmodel...

Of the I-Tasser ligands, the following were kept:

* peptide ligand from PDB:4Z83
* ATP from 1ATP
* cAMP and cGMP in CNB domains (see [code crystal note](code_crystal.md)).

### Merger

Swissmodel 5DYK-threaded model was used a loop donor for 418-449. Although this is an unimportant part.
All non peptides moved to chain X.

PyMOL cmd:
    
    create combo, 5C6C.r or 5BV6.r or (liganded and resi 1-149 and chain A) or loop or (liganded and resi 449-999999 and chain A) or (liganded and chain B) or (liganded and resn MN) or (liganded and resn ATP)
    
### Pre-Loops

Load:

    pose = get_pose(f'combined.pdb', params_filenames=['35G.params','CMP.params', 'ATP.params'])
    scorefxn = pyrosetta.get_fa_scorefxn()


Prevent all from moving too much

    constraint = pyrosetta.rosetta.protocols.constraint_generator.AtomPairConstraintGenerator()
    #chain_A_sele = pyrosetta.rosetta.core.select.residue_selector.ChainSelector('A')
    #constraint.set_residue_selector(chain_A_sele)
    true_sele = pyrosetta.rosetta.core.select.residue_selector.TrueResidueSelector()
    constraint.set_residue_selector(true_sele)
    constraint.set_ca_only(True)
    constraint.set_use_harmonic_function(True)
    constraint.set_max_distance(2.0)
    constraint.set_sd(0.5)
    # add to pose
    add_csts = pyrosetta.rosetta.protocols.constraint_generator.AddConstraints()
    add_csts.add_generator(constraint)
    add_csts.apply(pose)
    # enable in scorefxn
    stm = pyrosetta.rosetta.core.scoring.ScoreTypeManager()
    scorefxn.set_weight(stm.score_type_from_name("atom_pair_constraint"), 2)
    
Initially this was set to x10 and ony to chain A. This resulted in the expulsion of chain A.
(movemap jump is False by default). Adding checkpoints showed it was the cartesian step that was blowing up.
    
Select not crystals

    first = pyrosetta.rosetta.core.select.residue_selector.ResidueSpanSelector(153, 267)
    second = pyrosetta.rosetta.core.select.residue_selector.ResidueSpanSelector(270, 415)
    or_sele = pyrosetta.rosetta.core.select.residue_selector.OrResidueSelector(first, second)
    nor_sele = pyrosetta.rosetta.core.select.residue_selector.NotResidueSelector(or_sele)
    n = nor_sele.apply(pose)
    
Relax non-crystals with fixed bb x3

    relax = pyrosetta.rosetta.protocols.relax.FastRelax(scorefxn, 3)
    movemap = pyrosetta.MoveMap()
    movemap.set_bb(False)
    movemap.set_chi(allow_chi=n)
    relax.apply(pose)
    pose.dump_pdb('combined.preR1.fixedbb.pdb')
    print('Fixed BB done')
    
Relax non-crystals with flexible bb x3

    relax = pyrosetta.rosetta.protocols.relax.FastRelax(scorefxn, 3)
    movemap = pyrosetta.MoveMap()
    movemap.set_bb(allow_bb=n)
    movemap.set_chi(allow_chi=n)
    relax.apply(pose)
    pose.dump_pdb('combined.preR1.flexibb.pdb')
    print('Flexible BB done')
    
Relax cartesianly x3. This actually might be the step that bows up. Skipped.

    scorefxn_cart = pyrosetta.create_score_function('ref2015_cart')
    scorefxn_cart.set_weight(stm.score_type_from_name("atom_pair_constraint"), 10)
    relax = pyrosetta.rosetta.protocols.relax.FastRelax(scorefxn_cart, 3)
    relax.minimize_bond_angles(True)
    relax.minimize_bond_lengths(True)
    relax.apply(pose)
    pose.dump_pdb('combined.preR1.cart.pdb')
    
Final dihedral x3

    relax = pyrosetta.rosetta.protocols.relax.FastRelax(scorefxn, 3)
    relax.apply(pose)
    print('Flexible BB done')
    
Save
    
    pose.dump_scored_pdb(f'combined.r1.pdb', pyrosetta.get_fa_scorefxn())
    
### Loops

Functions to fix loops. As discussed in [my blog](https://blog.matteoferla.com/2020/07/filling-missing-loops-proper-way.html).

    def remove_termini(pose, preceding:int):
        # preceding precedes the discontinuity
        # LOWER_CONNECT N
        # UPPER_CONNECT C
        rm_upper = pyrosetta.rosetta.core.conformation.remove_upper_terminus_type_from_conformation_residue
        rm_lower = pyrosetta.rosetta.core.conformation.remove_lower_terminus_type_from_conformation_residue
        rm_upper(pose.conformation(), preceding)
        rm_lower(pose.conformation(), preceding)
        rm_lower(pose.conformation(), preceding + 1)
        rm_upper(pose.conformation(), preceding + 1)
    
    def fix_loop(pose, start, stop, cutpoint):
        remove_termini(pose, cutpoint)
        # cutpoint: resi before discontinuity
        lm = pyrosetta.rosetta.protocols.loop_modeler.LoopModeler()
        loops = pyrosetta.rosetta.protocols.loops.Loops()
        loop = pyrosetta.rosetta.protocols.loops.Loop(start,
                                                      stop, 
                                                      cutpoint) #cutpoint
        #loop.auto_choose_cutpoint(pose)
        loops.add_loop(loop)
        lm.set_loops(loops)
        # these are enabled by default. here for quick changing.
        lm.enable_centroid_stage()
        lm.enable_fullatom_stage()
        lm.enable_build_stage()
        lm.apply(pose)
        
First remove ligands or it will segfault.

    pyrosetta.rosetta.protocols.constraint_generator.RemoveConstraints().apply(pose)
    original = pose.clone()
    pyrosetta.rosetta.core.pose.remove_nonprotein_residues(pose)
    
Then fix junctions:
    
    fix_loop(pose, start=146, stop=151, cutpoint=149)
    fix_loop(pose, start=267, stop=270, cutpoint=268)
    fix_loop(pose, start=416, stop=450, cutpoint=418)
    fix_loop(pose, start=416, stop=450, cutpoint=448)
    pose.dump_scored_pdb(f'combined.r1.loopsNoLigs.pdb', pyrosetta.get_fa_scorefxn())
    
Finally return the ligands

    original = get_pose(f'combined.r1.pdb', params_filenames=['35G.params','CMP.params'])
    # select ligands
    LIGAND = pyrosetta.rosetta.core.chemical.ResidueProperty.LIGAND
    lig_sele = pyrosetta.rosetta.core.select.residue_selector.ResiduePropertySelector(LIGAND)
    # get residue vector to find start and stop
    rv = pyrosetta.rosetta.core.select.residue_selector.ResidueVector(lig_sele.apply(original))
    pyrosetta.rosetta.core.pose.append_subpose_to_pose(pose, 
                                                       original,
                                                       min(rv), 
                                                       max(rv),
                                                       True)
    pose.dump_scored_pdb(f'combined.r1.loops.pdb', pyrosetta.get_fa_scorefxn())

### Post-Loops

Relax with CA constraint of x2.

    constraint = pyrosetta.rosetta.protocols.constraint_generator.AtomPairConstraintGenerator()
    chain_A_sele = pyrosetta.rosetta.core.select.residue_selector.ChainSelector('A')
    constraint.set_residue_selector(chain_A_sele)
    constraint.set_ca_only(True)
    constraint.set_use_harmonic_function(True)
    constraint.set_max_distance(2.0)
    constraint.set_sd(0.5)
    add_csts = pyrosetta.rosetta.protocols.constraint_generator.AddConstraints()
    add_csts.add_generator(constraint)
    add_csts.apply(pose)
    stm = pyrosetta.rosetta.core.scoring.ScoreTypeManager()
    scorefxn.set_weight(stm.score_type_from_name("atom_pair_constraint"), 2)
    # relax all
    relax = pyrosetta.rosetta.protocols.relax.FastRelax(scorefxn, 5)
    relax.apply(pose)
    pose.dump_scored_pdb(f'combined.relaxed.pdb', pyrosetta.get_fa_scorefxn())
    
And a final inspection:

    import ngview
    nglview.show_rosetta(pose)