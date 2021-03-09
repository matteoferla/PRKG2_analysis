## Extension
Loading stuff
```jupyterpython
import pyrosetta
from pyrosetta_help.init_ops import make_option_string, configure_logger
from pyrosetta_help.common_ops import pose_from_file

logger = configure_logger()
pyrosetta.distributed.maybe_init(extra_options=make_option_string(no_optH=False,
                                                ex1=None,
                                                ex2=None,
                                                mute='all',
                                                ignore_unrecognized_res=True,
                                                load_PDB_components=False,
                                                ignore_waters=False)
                               )
ref = pose_from_file('final.pdb', params_filenames=['35G.params','CMP.params', 'ATP.params', 'NME.params'])
pose = ref.clone()
```                               
Adding extensions. The following were tried:

* N-methylamide (NME)
* Glycine
* Glu-Leu-Leu

Constraints for turn and last residue's side chain to backbone were created e.g.
```jupyterpython
def get_AtomID(chain:str, resi:int, atomname:str) -> pyrosetta.rosetta.core.id.AtomID:
    r = pose.pdb_info().pdb2pose(res=resi, chain=chain)
    assert r != 0, f'{resi}:{chain} is absent'
    residue = pose.residue(r)
    return pyrosetta.rosetta.core.id.AtomID(atomno_in=residue.atom_index(atomname), rsd_in=r)

HarmonicFunc = pyrosetta.rosetta.core.scoring.func.HarmonicFunc
AtomPairConstraint = pyrosetta.rosetta.core.scoring.constraints.AtomPairConstraint
cons = []
cons.append( AtomPairConstraint(get_AtomID('A', 759, 'O'),
                                get_AtomID('A', 762, 'N'),
                                HarmonicFunc(x0_in=3.3, sd_in=0.2)
                                )
            )
cons.append( AtomPairConstraint(get_AtomID('A', 762, 'CG'),
                                get_AtomID('A', 527, 'O'),
                                HarmonicFunc(x0_in=5.3, sd_in=0.2)
                                )
            )
cl = pyrosetta.rosetta.utility.vector1_std_shared_ptr_const_core_scoring_constraints_Constraint_t()
cl.extend(cons)
cs = pyrosetta.rosetta.core.scoring.constraints.ConstraintSet()
cs.add_constraints(cl)
setup = pyrosetta.rosetta.protocols.constraint_movers.ConstraintSetMover()
setup.constraint_set(cs)
# setup.apply(pose)  # if it were not looped
```
Extensions and calculations (10 replicates at constraint different weights)
```
```jupyterpython
reg_scorefxn = pyrosetta.get_fa_scorefxn()
for w in (20, 10, 5, 0):
    for i in range(10):
        pose = ref.clone()
        #D761E
        pyrosetta.rosetta.protocols.simple_moves.MutateResidue(target=761, new_res='GLU').apply(pose)
        #F762L
        pyrosetta.rosetta.protocols.simple_moves.MutateResidue(target=762, new_res='LEU').apply(pose)
        setup.apply(pose)
        from pyrosetta_help.common_ops import get_last_res_in_chain
        #ELLTEEKLITACTLQKRTSRINNPTHYFLFRVL*
        for resn in 'LTEE':
            chm = pyrosetta.rosetta.core.chemical.ChemicalManager.get_instance()
            resiset = chm.residue_type_set( 'fa_standard' )
            res_type = resiset.get_representative_type_name1(resn) #e.g. A
            residue = pyrosetta.rosetta.core.conformation.ResidueFactory.create_residue(res_type)
            r = get_last_res_in_chain(pose, 'A')
            rm_upper = pyrosetta.rosetta.core.conformation.remove_upper_terminus_type_from_conformation_residue
            rm_lower = pyrosetta.rosetta.core.conformation.remove_lower_terminus_type_from_conformation_residue
            rm_upper(pose.conformation(), r)
            rm_lower(pose.conformation(), r)
            pose.append_polymer_residue_after_seqpos(residue, r, True)
            pose.pdb_info().set_resinfo(r+1,
                                        'A',
                                        int(pose.pdb_info().pose2pdb(r).split()[0]) + 1)
        movemap = pyrosetta.MoveMap()
        span_sele = pyrosetta.rosetta.core.select.residue_selector.ResidueSpanSelector(759,764)
        neigh_sele = pyrosetta.rosetta.core.select.residue_selector.NeighborhoodResidueSelector(span_sele, 12, True)
        n = neigh_sele.apply(pose)
        targets = pyrosetta.rosetta.core.select.residue_selector.ResidueVector(n)
        movemap.set_bb(allow_bb=n)
        movemap.set_chi(allow_chi=n)
        movemap.set_jump(True)

        scorefxn = pyrosetta.get_fa_scorefxn()
        stm = pyrosetta.rosetta.core.scoring.ScoreTypeManager()
        scorefxn.set_weight(stm.score_type_from_name("atom_pair_constraint"), w)

        relax = pyrosetta.rosetta.protocols.relax.FastRelax(scorefxn, 15)
        relax.set_movemap_disables_packing_of_fixed_chi_positions(True)
        relax.set_movemap(movemap)
        relax.apply(pose)
        print(f'ELLTEE_{i}_x{w}', reg_scorefxn(pose) - reg_scorefxn(ref))
        pose.dump_pdb(f'ELLTEE_{i}_x{w}.pdb')
```
Results:
```jupyterpython
scores.groupby(['weight']).min()
```

|   weight |   score |
|---------:|--------:|
|        0 | 48.5476 |
|        5 | 37.469  |
|       10 | 34.1854 |
|       20 | 29.5288 |

The scores are without constraints. The reason why the scores are increasingly
worse without constraints is that the optimal solution becomes harder to find given
how deleterious the extension is.
The best unweighted variant features the C-terminal carboxyl of the D761ELLTEE in the pocket of F762.
Getting the top scoring x20 and doing a 5 cycle x10, x5, x0 gradient does not improve the score.

Relatedly, 
15 relax cycles with constraints gives:

DF761EL: +9.5 kcal/mol
DF761ELX[NME] 23.1 kcal/mol
DF761ELLX[NME] 16.2 kcal/mol

The per residue scores were inspected
```jupyterpython
from pyrosetta_help.common_ops import pose2pandas
pyrosetta.rosetta.protocols.constraint_generator.RemoveConstraints().apply(pose)
scorefxn = pyrosetta.get_fa_scorefxn()
scores = pose2pandas(pose)
scores.loc[scores.total_score > 5][['residue', 'total_score']]
rscores = pose2pandas(ref)
rs_map = dict(zip(rscores.residue, rscores.total_score))
scores['wt_score'] = scores.residue.apply(lambda v: rs_map[v] if v in rs_map else float('nan'))
scores['∆∆G'] = scores.total_score - scores.wt_score
scores.loc[scores['∆∆G'] >= 2][['residue', 'total_score']]
```