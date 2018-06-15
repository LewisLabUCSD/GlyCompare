# import external libraries
from glypy.structure.glycan import glycan
from glypy.structure.glycan import MAIN_BRANCH_SYM
from glypy.structure.monosaccharide import depth #,children,parents
from glypy.structure.constants import Stem
import numpy as np
# import glycompare libraries
import gc_extract_motif
import customize_motif_vector
import extract_motif
import motif_class
import __init__
import json_utility
import glycan_profile
import plot_glycan_utilities
import matplotlib.pyplot as plt
from glypy.io import glycoct

def glycan_stats(glycan):
	functions = {'total_branch_points':total_branch_points,...}

###############
## Branchpoints

# get total number of (non-root) branchpoints 
def total_branch_points(glycan,exclude_root=True):
	branches = set(a_glycan.branch_parent_map.values())
	if exclude_root:
		branches.remove(MAIN_BRANCH_SYM)
	return len(branches)

# total n-level (or all) branchpoints (from/to) [any or specific monosaccharides] (mannose/glcNAc/...)
# mannosaccharide: 'all' or 'glc' or ['glc','man'] 
# exclude_sacch: 'glc' or ['glc','man'] 
def n_level_branching(glycan,dir='from',mannosaccharide='all',exclude_sacch=None)
	assert dir in ['to','from']
	assert mannosaccharide == 'all' or set(monosaccharide) < set(Stem)

	# get nodes, depth of nodes, children and parents
	nodes =  [i for i in a_glycan.depth_first_traversal()]
	ndepth =  np.array([i for i in a_glycan.depth_first_traversal(apply_fn=depth)])
	sacch =  np.array([i.stem[0] for i in a_glycan.depth_first_traversal()]) ### bug here, sometimes .stem returns multiple sugars. this implimentation ignores that
	if dir=='from':
		relatives =  [len(i.children) for i in a_glycan.depth_first_traversal()]
	elif dir=='to':
		relatives =  [len(i.parents) for i in a_glycan.depth_first_traversal()]

	n_level_branches = {}
	for i in range(ndepth.max()):
		# identify sugars at depth i
		if monosaccharide=='all':
			idx = np.where(ndepth==i)[0]
		elif exclude_sacch is None:
			idx = np.where(ndepth==i and (sacch in mannosaccharide))[0]
		else:
			idx = np.where(ndepth==i and (sacch in mannosaccharide) and (sacch not in exclude_sacch))[0]
		if idx.size==0:
			continue
		# iterate over relevant sugars
		sum_i = 0
		for j in idx:
			# add a branch fo every value greater than 1
			sum_i += (relatives[j]-1)
		n_level_branches[i]=sum_i

	return n_level_branches


###############
## Depth & Elongation

# get depth
def branch_depth(glycan):
	 glycan.label_branches()

# In [39]: a_glycan.branch_parent_map
# Out[39]: {'a': '-', 'b': '-', 'c': 'b', 'd': 'b', 'e': 'd', 'f': 'd'}

# In [40]: a_glycan.branch_lengths
# Out[40]: defaultdict(int, {'-': 6, 'a': 1, 'b': 6, 'c': 6, 'd': 5, 'e': 5, 'f': 5})

# deepest branch, number of deepest branches

###############
## Sialation

# get number of sialic acids (optional: divide into Neu5Ac & Neu5Gc)

# get number sialic acids capped branches (optional: divide into Neu5Ac & Neu5Gc)

# count polysialation events -b14Gal-a23Sia-(a28Sia)_n s.t. n>1 (optional: human aka just Neu5Ac)

# count oligosialation events -b14Gal-a23Sia-a28Sia  (optional: human aka just Neu5Ac)

# count terminal sialation events -b14Gal-a2(3/6)Sia  (optional: human aka just Neu5Ac)

# number of noncanonical sialic acid caps (caps - (poly+oligo+terminal))

###############
## Metastructures

# branch symetry (options: n-th branchpoint vs each branch point, binary vs % symetry)
	# isolate each branchs at each level
		# topological_equality

# branch motif similairty
	# isolate each branch
		# number of common motifs
		# largest common motif
		# (significant common motif)
		# largest common motif == branch ?
		# glypy.algorithms/subtree_search/common_subgraph/MaximumCommonSubgraphSolver

###############
## Degree of modification

# manosylation percentage

# high mannose glucose