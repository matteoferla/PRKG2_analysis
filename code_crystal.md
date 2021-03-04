## Code: CNB domains

> This is a copy paste from a Jupyter notebook. Some functions are in
> [my pyrosetta_scripts repo](https://github.com/matteoferla/pyrosetta_scripts)

The protein were thus cleaned in PyMOL:

* CNB-A domain — PDB:5C6C — This has cAMP (`CMP`). Whereas the calcium ion was kept (two contact points), ethylene glycol, cadmium ions and cobalt ions were removed.
* CNB-B domain — PDB:5BV6 — This has cGMP (`35G`). A sodium ion was kept, but calcium and acetate ions were removed.


The ligands were parameterised as opposed to using the PDB component library
because they are charged and not all structures have the atom names that match.

35G:

```jupyterpython
from rdkit_to_params import Params
params = Params.from_smiles_w_pdbfile(pdb_file='5BV6_clean.pdb',
                              smiles='c1nc2c(n1[C@H]3[C@@H]([C@H]4[C@H](O3)CO[P@@](=O)(O4)[O-])O)N=C(NC2=O)N',
                            name='35G')
params.dump('35G.params')

import nglview
nglview.show_rosetta(params.test())
```
  
CMP:
    
```jupyterpython
from rdkit_to_params import Params
params = Params.from_smiles_w_pdbfile(pdb_file='5C6C_clean.pdb',
                              smiles='c1nc(c2c(n1)n(cn2)[C@H]3[C@@H]([C@H]4[C@H](O3)CO[P@](=O)(O4)[O-])O)N',
                            name='CMP')
params.dump('CMP.params')

import nglview
nglview.show_rosetta(params.test())
```

ATP:

```jupyterpython
from rdkit_to_params import Params
params = Params.from_smiles_w_pdbfile(pdb_file='ATP.pdb',
                              smiles='c1nc(c2c(n1)n(cn2)[C@H]3[C@@H]([C@@H]([C@H](O3)CO[P@@](=O)([O-])O[P@](=O)([O-])OP(=O)([O-])[O-])O)O)N',
                                      name='ATP',
                                     proximityBonding=False)
params.dump('ATP.params')

import nglview
nglview.show_rosetta(params.test())
```
    
The structures were minimised:

```jupyterpython
import pyrosetta
from init_helper import make_option_string, get_logger
logger = get_logger()
pyrosetta.distributed.maybe_init(extra_options=make_option_string(no_optH=False,
                                                ex1=None,
                                                ex2=None,
                                                mute='all',
                                                ignore_unrecognized_res=True,
                                                load_PDB_components=False,
                                                ignore_waters=False)
                               )
```
                             
For each, loaded:
      
```jupyterpython
code = '5C6C' # '5BV6'
pose = get_pose(f'{code}_clean.pdb', params_filenames=['35G.params','CMP.params'])
ED = prep_ED(pose, f'{code}.ccp4')
print(ED.matchPose(pose))
```
    
Check all is good:

```jupyterpython
import nglview
nglview.show_rosetta(pose)
```

Three cycles of cartesian minimisation with `elec_dens_fast` at x30.
    
```jupyterpython
scorefxn_cart = pyrosetta.create_score_function('ref2015_cart')
elec_dens_fast = pyrosetta.rosetta.core.scoring.ScoreType.elec_dens_fast
scorefxn_cart.set_weight(elec_dens_fast, 30)
relax = pyrosetta.rosetta.protocols.relax.FastRelax(scorefxn_cart, 3)
relax.minimize_bond_angles(True)
relax.minimize_bond_lengths(True)
relax.apply(pose)
```
    
Then 15 cycles regular dihedral minimisation with `elec_dens_fast` at x30.

```jupyterpython
scorefxn = pyrosetta.get_fa_scorefxn()
elec_dens_fast = pyrosetta.rosetta.core.scoring.ScoreType.elec_dens_fast
scorefxn.set_weight(elec_dens_fast, 30)
relax = pyrosetta.rosetta.protocols.relax.FastRelax(scorefxn, 15)
relax.apply(pose)
```
    
.

```jupyterpython
pose.dump_scored_pdb(f'{code}.r.pdb', pyrosetta.get_fa_scorefxn())
```
    
The `elec_dens_fast` was not reduced as once combined they will be relaxed again.