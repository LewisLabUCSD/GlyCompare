
########################################################################################
"""editing log for glypy
1.  plypy.structure.monosaccharide 
        class Monosaccharide
              def __init__(self,....)
        add line
                 self.proportion = 0
       
2. from glypy.plot.buchheim import buchheim
       change def first_walk(v, distance=0.6, visited=None)
       change def second_walk(v, m=0, depth=0, min=None, visited=None)
               min = second_walk(w, m + v.mod, depth + 0.55, min, visited=visited)
       
3. change draw_tree.py
       DEFAULT_SYMBOL_SCALE_FACTOR = 0.25

4. whole plot_glycan_utility is wrapped and modified from plot() function from plot.draw_tree
"""

########################################################################################
# Basic editing




# setting up the basic directory
root_ = "/Users/apple/PycharmProjects/GlyCompare/"
num_processors = 8
json_address = root_ + "generated_json_file/"
motif_plot_address = root_ + "motif_plot/"
manual_curated_address = root_ + "glycan_structure/"
plot_output_address = root_ + "output_plot/"
source_address = root_ + "source_data/"



########################################################################################
# Part 1 check extract_motif.py
# for gc_extract motif
# """
# def load_glycoct_for_database():
#     1. get the glycanID from Glycan_topolog_list
#     2. find ID in glytoucan database: /root_address + r'data_dic_finnn.json'
#     3. find ID in self-generated local file: /NBT_init.json_address+_code+".glycoct_condensed"
#     4. output a dict ID str -> glycoct str stored in: root_address + 'NBT_for_motif_extraction.json'
# """
# # input
glytoucan_data_base_addr__ = json_address + r'data_dic_finnn.json'
topology_list = json_address + r'Glycan_topolog_list.txt'
# # output
glycoct_dict_goto_extraction_addr = json_address + 'NBT_for_motif_extraction.json'
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
glytoucan_data_base_addr__ = json_address + r'data_dic_finnn.json' # from above
glycoct_dict_goto_extraction_addr = json_address + 'NBT_for_motif_extraction.json'

# glytoucan_data_base_addr__ # from above
# output
success_log_addr = json_address + 'NBT_for_motif_log.json'
glycan_list_addr = json_address + 'NBT_glycan_list_127.json'
# this one will go into the next step
glycan_dict_motif_list_addr = json_address + 'NBT_glycan_dict_degree_list_glycoct_for_motif.json'



########################################################################################
# part 2 check customize_motif_vector.py

"""

def get_motif_dict_degree_list_pipe(glycan_dict, output_motif_dic_degree_list_addr):
    merge the substructure of all glycans into motif dict
    :param glycan_dict: degree -> [motif1, motif2, ... ]/NBT_glycan_dict_degree_list_glycoct_for_motif
    store the glycan motif to NBT_motif_dic_degree_list.json
    :return: sorted motif_vec
    """
# input
glycan_dict_motif_list_addr # from part 1

# output
output_motif_dic_degree_list_addr = json_address + "NBT_motif_dic_degree_list.json"




"""def motif_matching_wrapper(motif_dict, glycan_with_motif_dict, id_list, matched_glycan_dict_addr):
    :param motif_dict: degree - >[motif1, motif2, ...]  /NBT_motif_dic_degree_list
    :param glycan_with_motif_dict: degree -> [motif1, motif2, ... ] /NBT_glycan_dict_degree_list_glycoct_for_motif
    :param id_list: all GlytoucanID of the glycans you are analyzing /NBT_fixed_gylcan_name_list
    :param matched_glycan_dict_addr: output_addr /NBT_fixed_gylcan_name_list
    :return: glycan_match_existed_motif degree - >[motif1, motif2, ...]
"""
# input file
glycan_dict_motif_list_addr # from above
output_motif_dic_degree_list_addr # from above
NBT_fixed_gylcan_name_list_addr= json_address + "NBT_fixed_gylcan_name.json"

# output file
output_matched_glycan_addr = json_address + "NBT_glycan_match_existed_motif.json"




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


aaa = ['WT',
       'mgat4A',
       'mgat4A/mgat4B',
       'mgat5',
       'mgat4A/mgat4B/mgat5',
       'B4GalT1',
       'B4GalT2',
       'B4GalT3',
       'B4GalT4',
       'B4GalT1/B4GalT2',
       'B4GalT1/B4GalT3',
       'B3gnt1',
       'B3gnt2',
       'st3gal3',
       'st3gal4',
       'st3gal6',
       'st3gal3/st3gal4',
       'st3gal4/st3gal6',
       'KI_ST6GalNAc1/st3gal4/st3gal6',
       'B3gnt2/mgat4a/mgat4b/mgat5',
       'st3gal4/st3gal6/mgat4a/mgat4b/mgat5',
       'KI_ST6GalNAc1/st3gal4/st3gal6/mgat4a/mgat4b/mgat5',
       'EPO48(mgat3)',
       'EPO143(mgat4C)',
       'EPO174(mgat2)',
       'EPO200(B4galt1/B4galt2/B4galt3)',
       'EPO275(B3gnt8)',
       'EPO78(mgat4B)',
       'EPO104(mgat5B)',
       'EPO127(mgat1)',
       'EPO259(mgat2/st3gal4/st3gal6)',
       'EPO261(mgat2/mgat4A/mgat4B/mgat5)',
       'EPO263(mgat2/st3gal4/st3gal6/magt4A/mgat4B/mgat5)',
       'EPO266(fut8)']
len(aaa)










