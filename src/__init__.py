from glypy.io import glycoct

from json_utility import *

# %matplotlib inline
"""
root address is json address
"""
root_address = "/Users/apple/PycharmProjects/nbt_glycan_profile/generated_json_file/"
json_address = root_address
motif_plot_address = "/Users/apple/PycharmProjects/nbt_glycan_profile/motif_plot/"
manual_curated_address = "/Users/apple/PycharmProjects/nbt_glycan_profile/glycan_structure/"
plot_output_address = "/Users/apple/PycharmProjects/nbt_glycan_profile/output_plot/"
source_address = "/Users/apple/PycharmProjects/nbt_glycan_profile/source_data/"


def load_motif_vec():
    return load_json(root_address + "BNT_motif_vec.json")


# def load_motif_vec_obj():
#     return [glycoct.loads(i) for i in load_json(root_address + "BNT_motif_vec.json")]


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


"""editing log from glypy
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






