# break glycoCT
from glypy.io import glycoct, iupac
from glypy.algorithms.subtree_search import subtree_of
from glypy.plot import plot
from glypy.structure.glycan import fragment_to_substructure
from glypy.io.glycoct import dump
from matplotlib import pyplot as plt
import seaborn as sns;

sns.set(color_codes=True)

from json_utility import load_json
import time
import numpy as np
import NBT_glycan_profile
import NBT_data_processing_pip
import NBT_init
import NBT_glycan_motif

from json_utility import *
import pandas as pd
import seaborn
import scipy
from glypy.io import glycoct
import numpy as np


def find_greatest_common_divisor(motif_list, cluster_dex, vec_):
    """
    two list of glycan obj

    :param motif_list:
    :param motif_vec:
    :return:
    """
    motif_vec = vec_
    motif_cluster_obj = [motif_vec[j] for j in motif_list]
    # print(motif_list)
    _row = len(motif_list)
    _col = len(motif_vec)
    _match_vec = np.zeros((_row, _col))

    for i in range(_row):
        try:
            aa = motif_match_motif_vec(motif_cluster_obj[i], motif_vec)
            _match_vec[i, :] = aa
        except:
            print(cluster_dex)
    # print(_match_vec)
    _final_vec = []
    for i in range(_col):
        if sum(_match_vec[:, i]) == _row:
            _final_vec.append(i)
    # print(_final_vec)
    _, l_dex = get_the_highest_degree([motif_vec[i] for i in _final_vec])
    return [motif_vec[_final_vec[i]] for i in l_dex], [_final_vec[i] for i in l_dex]


def find_the_least_common_multiple():
    pass


def get_the_highest_degree(motif_group):
    high_degree = 0
    a_vec = []
    for i in range(len(motif_group)):
        if len(motif_group[i]) > high_degree:
            high_degree = len(motif_group[i])
            a_vec = [i]
        else:
            a_vec.append(i)
    return high_degree, a_vec


def motif_match_motif_vec(a_tree, motif_vec):
    """
    :param a_tree:
    :param a_vec:
    :return:
    """
    # print(a_tree)
    _checked_matrix = np.zeros((1, len(motif_vec)))
    for jdex, j in enumerate(motif_vec):
        #             print(type(j), type(_tree))
        if not subtree_of(j, a_tree, False) is None:
            # _checked_matrix[jdex] = 1
            # print('finished ', )
            _checked_matrix[0, jdex] = 1
    return _checked_matrix



