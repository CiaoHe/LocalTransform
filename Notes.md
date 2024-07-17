function `def extract_from_reaction(reaction, setting = default_setting):`
```python
reactants_list, products_list, reagents_split_reagents(reaction:dict)# all return list are List[SMILES]
# clean not in product_maps atom-mapping

# 这一步的目的是为了进一步分开一些混在reactants中的reagents
for reactant in reactants_:
    is_reagent, max_num =  extend_atom_tag(reactant, max_num)
    # 如果reactant is_reagent, 直接返回
    # 如果不是，则对atom周围的neighbor_atom打上atom_map, index从max_num+1开始
    
# try to sanitize reactants and products, including remove Hs
...

# Calculate changed atoms
changed_atoms, changed_atom_tags, err = get_changed_atoms(reactants, products)
"""
- prod_atoms, prod_atom_tags = get_tagged_atoms_from_mols(products)
- reac_atoms, reac_atom_tags = get_tagged_atoms_from_mols(reactants)
- # check: len(set(prod_atom_tags)) == len(set(reac_atom_tags)); check: len(prod_atoms) == len(reac_atoms)
- # find Product atoms that are different from reactant atom equivalent, add those reac_atom and corresponding tag
- # find Reactant atoms that do not appear in product (tagged leaving groups)
- return changed_atoms, changed_atom_tags, err
"""

# find key fragments
reactant_fragments, intra_only, dimer_only = get_fragments_for_changed_atoms(reactants, changed_atom_tags, "reactant")
product_fragments, _, _ = get_fragments_for_changed_atoms(reactants, changed_atom_tags, "product")
"""
for mol in mols:
	# Initialize list of replacement symbols
	symbol_replacements = []
	# Build list of atoms to use
	atoms_to_use = []
	for atom in mol.GetAtoms():
		# check self (only tagged atoms)
		atom_to_use.append(atom.GetIdx())
		symbol = get_strict_smarts_for_atom(atom)
		# remove chiral information in smarts
		# add to symbol_replacements if changed in smarts
		if symbol != atom.GetSmarts():
            symbol_replacements.append((atom.GetIdx(), symbol))
    # collecting all leaving groups if RETRO
    ...
    # update symbols with symbol_replacements
    ...
    # reset molAtomMapNumber
 	mol_copy = deepcopy(mol)
    [x.ClearProp('molAtomMapNumber') for x in mol_copy.GetAtoms()] 
"""
# link fragments
rxn_string = reactant_fragments + '>>' + product_fragment
# atom_dict
atom_dict = {str(atom.GetAtomMapNum()): {'charge': atom.GetFormalCharge(), 'Hs': atom.GetNumExplicitHs()} for atom in changed_atoms}
# canonicalize_transform
    # - canonicalize_template(), sorting like
    # - reassign_atom_mapping(transform, atom_dict) ! very complicated
rxn_canonical, replacement_dict = canonicalize_transform(rxn_string, atom_dict) # remapped rxn_string
canonical_template = canonicalize_smarts(reactants_string)+ '>>'+canonicalize_smarts(products_string)
# check canonical_template valid
rxn = AllChem.ReactionFromSmarts(canonical_template)
if rxn.Validate()[1]!=0: gg

edits, H_change, Charge_change, Chiral_change = match_label(
	reactants_smiles, products_smiles, replacement_dict, changed_atom_tags, 
    retro=RETRO, remote=REMOTE,
)
# 1. label_CHS_change(smiles1, smiles2, edit_num, replacement_dict, use_stereo) -> return atom_map_dict1, H_change::[atom_map_num:value], Charge_change, Chiral_change
# 2. label_foward_edit_site(smiles1, smiles2, edit_num) -> formed_bonds, broken_bonds, changed_bonds, remote_bonds
# # edits: ['A': (bond_idxs::List[(bond_start_atom_id, bond_end_atom_id),], bond_maps, bond_temps)]
```
