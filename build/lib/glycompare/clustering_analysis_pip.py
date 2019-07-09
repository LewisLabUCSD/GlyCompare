import seaborn as sns

import nglycan_alignment
import plot_glycan_utilities

sns.set_palette("RdBu_r", 7)
import time
import numpy as np
import glycan_profile
import customize_motif_vec
import __init__
import motif_class
from  importlib import reload
import matplotlib.pyplot as plt
# sns.palplot(sns.color_palette("RdBu_r", 7))
import json_utility
import pandas as pd
import scipy
from glypy.io import glycoct
import glycan_io
import numpy as np

# profile_name = ['WT',
#                 'mgat4A',
#                 'mgat4A/mgat4B',
#                 'mgat5',
#                 'mgat4A/mgat4B/mgat5',
#                 'B4GalT1',
#                 'B4GalT2',
#                 'B4GalT3',
#                 'B4GalT4',
#                 'B4GalT1/B4GalT2',
#                 'B4GalT1/B4GalT3',
#                 'B3gnt1',
#                 'B3gnt2',
#                 'st3gal3',
#                 'st3gal4',
#                 'st3gal6',
#                 'st3gal3/st3gal4',
#                 'st3gal4/st3gal6',
#                 'KI_ST6GalNAc1/st3gal4/st3gal6',
#                 'B3gnt2/mgat4a/mgat4b/mgat5',
#                 'st3gal4/st3gal6/mgat4a/mgat4b/mgat5',
#                 'KI_ST6GalNAc1/st3gal4/st3gal6/mgat4a/mgat4b/mgat5',
#                 'EPO48(mgat3)',
#                 'EPO143(mgat4C)',
#                 'EPO174(mgat2)',
#                 'EPO200(B4galt1/B4galt2/B4galt3)',
#                 'EPO275(B3gnt8)',
#                 'EPO78(mgat4B)',
#                 'EPO104(mgat5B)',
#                 'EPO127(mgat1)',
#                 'EPO259(mgat2/st3gal4/st3gal6)',
#                 'EPO261(mgat2/mgat4A/mgat4B/mgat5)',
#                 'EPO263(mgat2/st3gal4/st3gal6/magt4A/mgat4B/mgat5)',
#                 'EPO266(fut8)']


def draw_motif_cluster(g, df, name_prefix, color_threshold, fig_size=(10,35)):
    plt.figure(figsize=fig_size)
    plt.title('Hierarchical Clustering Profile', fontdict={'fontsize': 25})
    plt.xlabel('Distance', fontdict={'fontsize': 25})
    plt.rc_context({'lines.linewidth': 4})
    plt.ylabel('KO Gene Name', fontdict={'fontsize': 25})
    den = scipy.cluster.hierarchy.dendrogram(g.dendrogram_row.linkage,
                                             # truncate_mode='lastp',show_contracted=True,p=50,
                                             labels=df.index,
                                             color_threshold=color_threshold, orientation='left', leaf_font_size=10,
                                             distance_sort='descending', leaf_rotation=0)

    #     plt.show()
    plt.savefig(r"./" + name_prefix + 'motif_cluster.png')
    cccluster_dict = {}
    # for i in zip(scipy.cluster.hierarchy.fcluster(g.dendrogram_row.linkage, t=color_threshold, criterion='distance'),df.index.tolist()):
    #     print(i)


    for i, j in zip(scipy.cluster.hierarchy.fcluster(g.dendrogram_row.linkage, t=color_threshold, criterion='distance'),
                    df.index.tolist()):
        if i not in cccluster_dict.keys():
            cccluster_dict[i] = [j]
        else:
            cccluster_dict[i].append(j)
    return cccluster_dict


def draw_glycan_cluster(glyco_motif_cluster, name_prefix, motif_vec, motif_weights_dict, address, plot_rep=True, threshold_list=[0.5, 0.6, 0.7]):
    # vec_ = load_json(NBT_init.root_address + 'NBT_motif_vec.json')
    _count = 0
    for i in range(1, len(glyco_motif_cluster.keys()) + 1):
        list_ = glyco_motif_cluster[i]
        name_list = [str(j)+':'+str(round(motif_weights_dict[j], 4)) for j in list_]
        _count += len(list_)
        plot_glycan_utilities.plot_glycan_list([motif_vec[i] for i in list_], name_list, str(i))
        print(glyco_motif_cluster[i])
        plt.savefig(__init__.plot_output_address + name_prefix + str(i) + '.png')
        ## _temp_node, _name = NBT_motif_match_motifvec.find_greatest_common_divisor(glyco_motif_cluster[i], i, vec_)
        # plot_glycan_utilities.plot_glycan(_temp_node[0], _name[0])
        if plot_rep:
            a_panel = nglycan_alignment.glycan_model()
            for j in glyco_motif_cluster[i]:
                # print(j)
                # plot_glycan_utilities.plot_glycan(vec_[i], title=str(i))
                gly_nglycan_dict = nglycan_alignment.traves_glycan(motif_vec[j], weight=motif_weights_dict[j])
                a_panel.glycan_walk(gly_nglycan_dict)

                # NBT_nglycan_alignment.travel_str_dict(a_panel.panel)

                # plot_glycan_utilities.plot_glycan(a_panel.get_common_representative(0.1), title=0.1)
            glycan_list = a_panel.get_reps(threshold_list=[0.51])
            glycan_io.save_glycan(glycan_list[0], address+name_prefix+str(i)+'.core.glycoct')
            plot_glycan_utilities.plot_glycan(glycan_list[0], title=str(i), addr=address+name_prefix+str(i)+'.core.eps')

        # plt.savefig(__init__.plot_output_address + name_prefix + str(i) + '.png')
        # plot_glycan_utilities.plot_glycan_list(glycan_list, idex_list=[str(k) for k in [0.51, 0.6, 0.7]])
        # plt.savefig()


# print(len(cccluster_dict.keys()))
def draw_profile_cluster(g, df, profile_name, name_prefix, color_threshold, address):
    """
    three profiles assss
    """
    plt.figure(figsize=(15, 15))
    plt.title('Hierarchical Clustering Profile', fontdict={'fontsize': 25})
    plt.xlabel('Distance', fontdict={'fontsize': 25})
    plt.rc_context({'lines.linewidth': 4})
    plt.ylabel('KO Gene Name', fontdict={'fontsize': 25})
    col_list = [str(x) for x in df.columns.tolist()]
    den = scipy.cluster.hierarchy.dendrogram(g.dendrogram_col.linkage, distance_sort='descending',
                                             labels=[', '.join(i) for i in zip(col_list, profile_name)],
                                             color_threshold=color_threshold, orientation='left', leaf_font_size=25.)
    plt.savefig(address + name_prefix + 'profile_cluster.png')
    return {}


"""The following are pipelines"""


def draw_glycan_clustermap(df_ncore, name_prefix, address, metric="braycurtis"):
    # draw clustermap
    g = sns.clustermap(df_ncore, metric=metric)
    g.savefig(address + name_prefix + 'clustermap.png')
    return g
    # draw profiles


    # draw motif clusters
    # draw_glycan_cluster(cccluster_dict, name_prefix)

#
# def motif_with_n_glycan_core_all_motif(motif_, existed_table, weight_dict, color_threshold=0.55):
#     motif_.get_motif_with_core()
#     ncore_dependence_tree = motif_.motif_with_ncore_dependence_tree()
#     dp_tree = motif_class.MotifDpTree(ncore_dependence_tree, motif_.motif_with_core_list, weight_dict)
#     #     # store_json(r"/Users/apple/PycharmProjects/abundance_matrix.json", weight_dict)
#     len(dp_tree.all_nodes)
#     # existed_table = _ggg.table_exist_or_not()
#     df_ncore = existed_table[existed_table.index.isin(dp_tree.all_nodes)]
#     g = sns.clustermap(df_ncore, metric="jaccard")
#     name_prefix = 'all_nodes.'
#     g.savefig(__init__.plot_output_address + name_prefix + 'clustermap.png')
#     draw_profile_cluster(g, df_ncore, profile_name, name_prefix, color_threshold)
#     cccluster_dict = draw_motif_cluster(g, df_ncore, name_prefix, color_threshold=0.35)
#     draw_glycan_cluster(cccluster_dict, name_prefix)
#
#
# def pipe_motif_ana_with_motif_abundance():
#     """for every glycan, use the motif counts(4 max) vec to represent
#         the abundance of the motif in a profile is represented by the sum of weight*abundance
# """
#
#     output_for_motif, glycan_dict_glycoct, glycan_dict = glycan_profile.get_glycan_glycoct()
#     NBT_dict_name_abundance_cross_profile = glycan_profile.load_glycan_profile(glycan_dict)
#     glycan_profile, glycan_profile_merged = glycan_profile.load_glycan_profile_dic()
#     profile_obj_list = glycan_profile.combine_profile_mz_with_motif_abundance(glycan_profile,
#                                                                                   NBT_dict_name_abundance_cross_profile)
#     _ggg = glycan_profile.glycan_profiles(profile_obj_list)
#     _table = _ggg.table_against_wt_relative_abd()
#     # _table
#     _table.to_csv(NBT_init.source_address + r"abundance_matrix.txt")
#     _np_mat = np.array(_table)
#     weight_dict = {}
#     for i in range(len(_np_mat)):
#         weight_dict[i] = list(_np_mat[i])
#     existed_table = _ggg.table_exist_or_not()
#     vec_ = load_json(NBT_init.root_address + 'NBT_motif_dic_degree_list.json')
#     motif_ = motif_class.GlycanMotifLib(vec_)
#     motif_with_n_glycan_core_all_motif(motif_, existed_table, weight_dict)
#     motif_with_n_glycan_core(motif_, existed_table, weight_dict)

# def get_abd_table_pip(match_dict_addr=__init__.json_address + "match_dict.json"):
#     # load CHO paper abundance table
#     mz_abd_table = glycan_profile.load_cho_mz_abundance(cho_addr=__init__.source_address + 'nbt.3280_cho.txt',
#                                                         mz_abd_addr=__init__.source_address + 'glycan_table.xls')
#     # load glycoprofile Mass Spectrum m/z and glycan structure info
#     profile_mz_to_id = glycan_profile.load_glycan_profile_dic()
#     # normalize CHO abundance table
#     norm_mz_abd_dict = glycan_profile.get_norm_mz_abd_table(mz_abd_table,
#                                                             norm_abd_table_dict_addr=__init__.json_address + "norm_mz_abd_dict.json")
#     # load match_dict
#     match_dict = json_utility.load_json(match_dict_addr)
#     # digitalize the glycoprofile
#     glycoprofile_list = glycan_profile.get_glycoprofile_list(profile_mz_to_id, norm_mz_abd_dict, match_dict)
#     # generate table
#     table_generator = glycan_profile.MotifAbdTableGenerator(glycoprofile_list)
#     motif_abd_table = table_generator.table_against_wt_relative_abd()
#     return motif_abd_table

#
# def motif_with_n_glycan_core(motif_lib, existed_table, color_threshold=0.95):
#     weight_dict = motif_class.get_weight_dict(existed_table)
#     dropper = motif_class.NodesDropper(motif_lib, weight_dict)
#     dropper.drop_node()
#     print("", len(dropper.drop_node()))
#     df_ncore = existed_table[existed_table.index.isin(dropper.nodes_kept)]
#     """ with n_glycan_core using jaccard for binary and use braycurtis for float
#     """
#     df_ncore.to_csv(__init__.json_address + r"abundance_matrix.txt")
#     name_prefix = 'dropped'
#     g = draw_glycan_clustermap(df_ncore, color_threshold, name_prefix)
#     draw_profile_cluster(g, df_ncore, profile_name, name_prefix, color_threshold)
#     cccluster_dict = draw_motif_cluster(g, df_ncore, name_prefix, color_threshold=0.23)

#
# def pipe_motif_ana_with_motif_existance():
#     """ for every glycan use the motif existance vec(0/1) to represent the existance of the motif
#         the abundance of the motif in a profile is represented by the sum of weight*existence
#     """
#
#     motif_abd_table = get_abd_table_pip()
#
#     motif_lib = motif_class.MotifLabNGlycan(
#         json_utility.load_json(__init__.merged_motif_dict_addr))  # unicarbkb_motifs_12259.json
#
#     weight_dict = motif_class.get_weight_dict(motif_abd_table)
#     dropper = motif_class.NodesDropper(motif_lib, weight_dict)
#     dropper.drop_node()
#     print("", len(dropper.drop_node()))
#     df_ncore = motif_abd_table[motif_abd_table.index.isin(dropper.nodes_kept)]
#     # draw plot
#     # motif_with_n_glycan_core_all_motif(motif_, _table, weight_dict)
#     """ with n_glycan_core using jaccard for binary and use braycurtis for float
#     """
#     df_ncore.to_csv(__init__.json_address + r"abundance_matrix.txt")
#     name_prefix = 'dropped'
#     g = draw_glycan_clustermap(df_ncore, name_prefix=name_prefix)
#     draw_profile_cluster(g, df_ncore, profile_name, name_prefix, color_threshold=0.95)
#     cccluster_dict = draw_motif_cluster(g, df_ncore, name_prefix, color_threshold=0.23)
#

if __name__ == '__main__':
    # pipe_motif_ana_with_motif_existance()
    pass