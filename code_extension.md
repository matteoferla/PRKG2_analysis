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

```jupyterpython
from pyrosetta_help.common_ops import get_last_res_in_chain

for resn in 'ELL':
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
print(pose.chain_sequence(1))
```
Then it was energy minimised
```jupyterpython
movemap = pyrosetta.MoveMap()
# 759 - 762 turn
span_sele = pyrosetta.rosetta.core.select.residue_selector.ResidueSpanSelector(r-5,r)
neigh_sele = pyrosetta.rosetta.core.select.residue_selector.NeighborhoodResidueSelector(span_sele, 12, True)
n = neigh_sele.apply(pose)
targets = pyrosetta.rosetta.core.select.residue_selector.ResidueVector(n)
print(f'{sum(n)} residues to be moved {targets}')
movemap.set_bb(allow_bb=n)
movemap.set_chi(allow_chi=n)
movemap.set_jump(True)
relax = pyrosetta.rosetta.protocols.relax.FastRelax(scorefxn, 15)
relax.set_movemap_disables_packing_of_fixed_chi_positions(True)
relax.set_movemap(movemap)
relax.apply(pose)
```
A few variations were done. Including with constraints for turn e.g.
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
cl = pyrosetta.rosetta.utility.vector1_std_shared_ptr_const_core_scoring_constraints_Constraint_t()
cl.extend(cons)
cs = pyrosetta.rosetta.core.scoring.constraints.ConstraintSet()
cs.add_constraints(cl)
setup = pyrosetta.rosetta.protocols.constraint_movers.ConstraintSetMover()
setup.constraint_set(cs)
setup.apply(pose)

stm = pyrosetta.rosetta.core.scoring.ScoreTypeManager()
scorefxn.set_weight(stm.score_type_from_name("atom_pair_constraint"), 10)
```
or last residue's side chain to backbone
```jupyterpython
cons.append( AtomPairConstraint(get_AtomID('A', 762, 'CZ'),
                                get_AtomID('A', 527, 'O'),
                                HarmonicFunc(x0_in=3.4, sd_in=0.2)
                                )
            )
```
Constraints were ommitted from final results due to incompatibility of extension with turn + F sidechain in core.