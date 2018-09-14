# import external libraries
from glypy.structure.glycan import glycan
from glypy.structure.glycan import MAIN_BRANCH_SYM
from glypy.structure.monosaccharide import depth #,children,parents
from glypy.structure.constants import Stem
import numpy as np
# import glycompare libraries
import gc_extract_motif
import customize_motif_vec
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

###############
## Linkage stats
# Get Linkages - Based on observations in motifID
def get_motif_linkages_dict(motifID, monosaccharide_glycoct2labels, output = __init__.motif_linkages_dict_addr):
	'''
	This function generates a motif_linkages_dict based on the linkages observed from the motifs in motifID. motifID could
	be got from its respective json file according to __init__.py file.
	'''
	import re
	import json
	from glypy.io import glycoct

	# Add code to verify if motifs are glycan objects or not
	try:
		motifID = [glycoct.loads(motif) for motif in motifID]
	except:
		pass

	motif_linkages_dict = {'simple': {}, 'specific': {}}

	# Get linkages from glycan(s)
	simple_linkages = []
	specific_linkages = []
	for motif in motifID:
		nodes = {}
		for node in motif.iternodes():
			nodes[str(node.id)] = node

		re_nodes = "\((\d+)\)"
		re_linkages = "\((-?\d\+-?\d)\)"
		for link in motif.link_index:
			nodeIDs = re.findall(re_nodes, str(link))
			node1 = monosaccharide_glycoct2labels[str(nodes[nodeIDs[0]])]
			node2 = monosaccharide_glycoct2labels[str(nodes[nodeIDs[1]])]
			linkage = '_' + re.findall(re_linkages, str(link))[0] + '_'
			simple_linkages.append(node1 + '_' + node2)
			specific_linkages.append(node1 + linkage + node2)

			simple_linkages = list(set(simple_linkages))
			specific_linkages = list(set(specific_linkages))

	# Assign index to each linkage
	link_id = 0
	for linkage in simple_linkages:
		motif_linkages_dict['simple'][link_id] = linkage
		link_id += 1

	link_id = 0
	for linkage in specific_linkages:
		motif_linkages_dict['specific'][link_id] = linkage
		link_id += 1

	print('There are ' + str(len(motif_linkages_dict['simple'])) + ' simple linkages')
	print('There are ' + str(len(motif_linkages_dict['specific'])) + ' simple linkages')
	with open(__init__.json_address + 'motif_linkages_dict.json', 'w') as f:
		json.dump(motif_linkages_dict, f)
	return motif_linkages_dict


def get_glycan_linkage_stats(glycan, monosaccharide_glycoct2labels, motif_linkages_dict):
	'''
	Function is to get a vector for a glycan, containing counts for each linkage that it has. Each count is located at the
	respective linkage index according to their IDs in motif_linkages_dict. Both monosaccharide_glycoct2labels and
	motif_linkages_dict could be found as json files according to __init__.py file.
	'''
	import re
	import numpy as np
	from glypy.io import glycoct
	from collections import Counter

	simple_linkage_stats = None
	specific_linkage_stats = None

	# Code to verify if glycan is a glycan object or not
	if type(glycan) == str:
		try:
			glycan = glycoct.loads(glycan)
		except:
			print(glycan + ' is not in GlycoCT format')

	_man1 = glycoct.loads("""
            RES
            1b:b-dman-HEX-1:5
            LIN""")

	if type(glycan) == type(_man1):
		# Get linkages from glycan(s)
		nodes = {}
		for node in glycan.iternodes():
			nodes[str(node.id)] = node

		simple_linkages = []
		specific_linkages = []
		re_nodes = "\((\d+)\)"
		re_linkages = "\((-?\d\+-?\d)\)"
		NoError = True
		for link in glycan.link_index:
			nodeIDs = re.findall(re_nodes, str(link))
			try:
				node1 = monosaccharide_glycoct2labels[str(nodes[nodeIDs[0]])]
			except:
				node1 = str(nodes[nodeIDs[0]])
				print(node1 + ' not in monosaccharide_glycoct2labels')
				NoError = False
				break
			try:
				node2 = monosaccharide_glycoct2labels[str(nodes[nodeIDs[1]])]
			except:
				node2 = str(nodes[nodeIDs[1]])
				print(node2 + ' not in monosaccharide_glycoct2labels')
				NoError = False
				break
			linkage = '_' + re.findall(re_linkages, str(link))[0] + '_'
			simple_linkages.append(node1 + '_' + node2)
			specific_linkages.append(node1 + linkage + node2)

		if NoError:
			# Get stats from linkages
			simple_linkage_freq = Counter(simple_linkages)
			specific_linkage_freq = Counter(specific_linkages)

			# Generate vector based on linkage observations contained in motif_linkages_dict
			simple_linkage_stats = np.zeros((1, len(motif_linkages_dict['simple'])))
			specific_linkage_stats = np.zeros((1, len(motif_linkages_dict['specific'])))

			for link in simple_linkage_freq.keys():
				for key, value in motif_linkages_dict['simple'].items():
					if str(link) == str(value):
						simple_linkage_stats[0][key] = int(simple_linkage_freq[link])

			for link in specific_linkage_freq.keys():
				for key, value in motif_linkages_dict['specific'].items():
					if str(link) == str(value):
						specific_linkage_stats[0][key] = int(specific_linkage_freq[link])
	else:
		print(str(glycan) + ' not analyzed, because is not a glycan object')
	return simple_linkage_stats, specific_linkage_stats


def get_multiple_glycan_linkage_stats(glycans, monosaccharide_glycoct2labels, motif_linkages_dict):
	'''
	Function is to get a matrix for multiple glycans, containing respective vectors with counts for each linkage that
	they have. Each count is located at the	respective linkage index according to their IDs in motif_linkages_dict.
	Both monosaccharide_glycoct2labels and	motif_linkages_dict could be found as json files according to __init__.py file.
	'''
	if isinstance(glycans, list):
		simple_linkages = np.zeros((len(glycans), len(motif_linkages_dict['simple'])))
		specific_linkages = np.zeros((len(glycans), len(motif_linkages_dict['specific'])))
		glycan_number = 0
		not_included = []
		for glycan in glycans:
			simple, specific = get_glycan_linkage_stats(glycan, monosaccharide_glycoct2labels, motif_linkages_dict)
			if (isinstance(simple, np.ndarray)) and (isinstance(specific, np.ndarray)):
				simple_linkages[glycan_number, :] = simple
				specific_linkages[glycan_number, :] = specific
			else:
				not_included.append(glycan_number)
			glycan_number += 1
		if len(not_included) > 0:
			print('Glycans not included in the analysis have the following indexes:')
			print(not_included)
	else:
		print('glycans must be a list')
		simple_linkages = None
		specific_linkages = None
	return simple_linkages, specific_linkages