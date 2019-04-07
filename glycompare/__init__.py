import getopt
import sys

exact_Ture = False
# root_ = "/Users/apple/PycharmProjects/GlyCompare/"
num_processors = 8

print("""
### There is one file must required and two files optional:
    1. abundance_table.xls

    if external profile name rather than index
    2. external_profile_naming.json

    if has mz or hplz name other than glycan_id
    3. glycoprofile_name_to_glycan_id.json

### There are several parameter should be set up first
    1. working_addr : root working dir

    2. project_name: usually same as the folder of root

    3. __init__.num_processors: number of processes needed

    4. __init__.exact_Ture: False for topology; True for exact linkage matching

    5. complex_profile: True is every profile has unique glycans; False is every profile has same name
    6. complex_naming: use the mz or hplc naming rather than glytoucan_id in abundance table.
    7. external_profile_naming: each glycoprofile has complex name rather than index.
        complex_profile and outside_name helps charactorize how complex of the naming is

    8. structure_loader: for loading glytoucan id, could be:
        list of glycan_id
        glycan_dict

    9. glytoucan_db_addr: the addr for glytoucan database if needed

    # meta_name: a gff type table which includes glycan's naming information

""")
# dir_required = ["intermediate_file", "output_plot", "source_data"]
########################################################################################
# """editing log for glypy
# 1.  plypy.structure.monosaccharide
#         class Monosaccharide
#               def __init__(self,....)
#         add line
#                  self.proportion = 0
#
# 2. from glypy.plot.buchheim import buchheim
#        change def first_walk(v, distance=0.6, visited=None)
#        change def second_walk(v, m=0, depth=0, min=None, visited=None)
#                min = second_walk(w, m + v.mod, depth + 0.55, min, visited=visited)
#
# 3. change draw_tree.py
#        DEFAULT_SYMBOL_SCALE_FACTOR = 0.25
#
# 4. whole plot_glycan_utility is wrapped and modified from plot() function from plot.draw_tree
# """

########################################################################################
# Basic editing


# %matplotlib inline

from glypy import plot

# for HMO


# for HMO
# setting up the basic directory
# exact_Ture = True
# json_address = root_ + "generated_json_file/"
# plot_output_address = root_ + "output_plot/"
# motif_plot_address = root_ + "motif_plot/"
# source_address = root_ + "source_data/"





########################################################################################
# Part 1 check extract_motif.py
# for gc_extract motif
"""
def load_glycoct_for_database():
    1. get the glycanID from Glycan_topolog_list
    2. find ID in glytoucan database: /root_address + r'data_dic_finnn.json'
    3. find ID in self-generated local file: /NBT_init.json_address+_code+".glycoct_condensed"
    4. output a dict ID str -> glycoct str stored in: root_address + 'NBT_for_motif_extraction.json'
"""
# input
# glytoucan_database_addr__ = source_address + r'data_dic_finnn.json'
# topology_list_addr = source_address + r'Glycan_topolog_list.txt'
# # # output
# glycan_dict_addr = json_address + 'glycan_dict.json'
#

# <___________Ben___Start from here _________________>
"""
def get_motif_pip(gly_len, prior=True):
    Please set the prior=True to get the data file please run the NBT_GLYCAN_preprocess file
    If prior=False, it will generate glycan motif for all glycan in glytoucan database
    1. load  dict ID -> glycoct str
    2. convert to glypy.glycan obj
    3. extract motif
    4. save glycan_motif_addr
"""
# input
# glytoucan_database_addr__ = source_address + r'data_dic_finnn.json'  # from above
# glycan_dict_addr = json_address + 'glycan_dict.json'
#
# # glytoucan_data_base_addr__ # from above
# # output
# # this one will go into the next step
# glycan_motif_dict_addr = json_address + 'glycan_motif_dict.json'

########################################################################################
# part 2 check customize_motif_vec.py

"""

def get_motif_dict_degree_list_pipe(glycan_dict, output_motif_dic_degree_list_addr):
    merge the substructure of all glycans into motif dict
    :param glycan_dict: degree -> [motif1, motif2, ... ]/NBT_glycan_dict_degree_list_glycoct_for_motif
    store the glycan motif to NBT_motif_dic_degree_list.json
    :return: sorted motif_vec
    """
# input
# glycan_dict_addr  # from part 1

# output
# merged_motif_dict_addr = json_address + "merged_motif_dict.json"

"""def motif_matching_wrapper(motif_dict, glycan_with_motif_dict, id_list, matched_glycan_dict_addr):
    :param motif_dict: degree - >[motif1, motif2, ...]  /NBT_motif_dic_degree_list
    :param glycan_with_motif_dict: degree -> [motif1, motif2, ... ] /NBT_glycan_dict_degree_list_glycoct_for_motif
    :param id_list: all GlytoucanID of the glycans you are analyzing /NBT_fixed_gylcan_name_list
    :param matched_glycan_dict_addr: output_addr /NBT_fixed_gylcan_name_list
    :return: glycan_match_existed_motif degree - >[motif1, motif2, ...]
"""
# input file
# glycan_dict_addr  # from above
# merged_motif_dict_addr  # from above

# output file
# output_matched_dict_addr = json_address + "match_dict.json"

########################################################################################
# part 3
# for plotting glycan, check plot_glycan_utilities.py
"""
# def plot_glycan_list(glycoct_list, idex_list=[], title='Glycans')

# def plot_glycan(tree, title='', at=(0, 0), ax=None, orientation='h', center=False, label=False,
                symbol_nomenclature='cfg', layout='balanced', **kwargs):
# def output_glycan_motif_vec_to_file
"""

#
# def load_motif_vec():
#     return load_json(json_address + "NBT_motif_vec.json")


# def load_motif_vec_obj():
#     return [glycoct.loads(i) for i in load_json(root_address + "NBT_motif_vec.json")]
#
#
# aaa = ['WT',
#        'mgat4A',
#        'mgat4A/mgat4B',
#        'mgat5',
#        'mgat4A/mgat4B/mgat5',
#        'B4GalT1',
#        'B4GalT2',
#        'B4GalT3',
#        'B4GalT4',
#        'B4GalT1/B4GalT2',
#        'B4GalT1/B4GalT3',
#        'B3gnt1',
#        'B3gnt2',
#        'st3gal3',
#        'st3gal4',
#        'st3gal6',
#        'st3gal3/st3gal4',
#        'st3gal4/st3gal6',
#        'KI_ST6GalNAc1/st3gal4/st3gal6',
#        'B3gnt2/mgat4a/mgat4b/mgat5',
#        'st3gal4/st3gal6/mgat4a/mgat4b/mgat5',
#        'KI_ST6GalNAc1/st3gal4/st3gal6/mgat4a/mgat4b/mgat5',
#        'EPO48(mgat3)',
#        'EPO143(mgat4C)',
#        'EPO174(mgat2)',
#        'EPO200(B4galt1/B4galt2/B4galt3)',
#        'EPO275(B3gnt8)',
#        'EPO78(mgat4B)',
#        'EPO104(mgat5B)',
#        'EPO127(mgat1)',
#        'EPO259(mgat2/st3gal4/st3gal6)',
#        'EPO261(mgat2/mgat4A/mgat4B/mgat5)',
#        'EPO263(mgat2/st3gal4/st3gal6/magt4A/mgat4B/mgat5)',
#        'EPO266(fut8)']
# len(aaa)
#
# aaa_re = ['WT',
#           'mgat4A',
#           'mgat4A.mgat4B',
#           'mgat5',
#           'mgat4A.mgat4B.mgat5',
#           'B4GalT1',
#           'B4GalT2',
#           'B4GalT3',
#           'B4GalT4',
#           'B4GalT1.B4GalT2',
#           'B4GalT1.B4GalT3',
#           'B3gnt1',
#           'B3gnt2',
#           'st3gal3',
#           'st3gal4',
#           'st3gal6',
#           'st3gal3.st3gal4',
#           'st3gal4.st3gal6',
#           'KI_ST6GalNAc1.st3gal4.st3gal6',
#           'B3gnt2.mgat4a.mgat4b.mgat5',
#           'st3gal4.st3gal6.mgat4a.mgat4b.mgat5',
#           'KI_ST6GalNAc1.st3gal4.st3gal6.mgat4a.mgat4b.mgat5',
#           'EPO48.mgat3.',
#           'EPO143.mgat4C.',
#           'EPO174.mgat2.',
#           'EPO200.B4galt1.B4galt2.B4galt3.',
#           'EPO275.B3gnt8.',
#           'EPO78.mgat4B.',
#           'EPO104.mgat5B.',
#           'EPO127.mgat1.',
#           'EPO259.mgat2.st3gal4.st3gal6.',
#           'EPO261.mgat2.mgat4A.mgat4B.mgat5.',
#           'EPO263.mgat2.st3gal4.st3gal6.magt4A.mgat4B.mgat5.',
#           'EPO266.fut8.']
#
# ########################################################################################
# # part 4
# # # for getting glycan stats, check glycan_stats.py
# # monosaccharide_glycoct2labels_addr = json_address + 'monosaccharide_glycoct2labels.json'
# # motifID_addr = json_address + ' Unicarbkb_motif_vec_12259.json'
# # motif_linkages_dict_addr = json_address + 'motif_linkages_dict.json'


if __name__ == '__main__':
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hdi:", ["help", "output="])
    except getopt.GetoptError:
        pass