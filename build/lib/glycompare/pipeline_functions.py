import os

import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import __init__
import merge_substructure_vec
import extract_substructures
import glycan_io
import process_glycoprofiles
import json_utility
import clustering_analysis
from clustering_analysis import draw_substructure_representative as draw_substructure_representative_pip


import select_motifs


def load_para_keywords(project_name, working_addr, **kwargs):
    """
    generate all necessary intermediate files
    :param project_name:
    :param working_addr:
    :param kwargs: any other parameter might be used
    :return: a comprehensive para dict
    """
    # for i in
    intermediate_dir = os.path.join(working_addr, "intermediate_file/")
    plot_output_dir = os.path.join(working_addr, "output_plot/")
    source_dir = os.path.join(working_addr, "source_data/")
    glycoct_dir = os.path.join(working_addr, 'source_data/glycoct/')

    name_to_id_addr = os.path.join(source_dir, 'glycan_identifier_to_glytoucan_id.json')
    # abundance_table_addr = os.path.join(source_dir, 'abundance_table')
    external_profile_naming_addr = os.path.join(source_dir, 'external_profile_naming.json')

    glycan_dict_addr = os.path.join(intermediate_dir, project_name + '_glycan_dict.json')
    glycan_motif_dict_addr = os.path.join(intermediate_dir, project_name + '_glycan_motif_dict.json')
    motif_dict_addr = os.path.join(intermediate_dir, project_name + "_motif_dict.json")
    matched_dict_addr = os.path.join(intermediate_dir, project_name + "_matched_dict.json")
    motif_abd_table_addr = os.path.join(intermediate_dir, project_name + "_motif_abd_table.csv")
    substructure_abd_table_addr = os.path.join(intermediate_dir, project_name + '_substructure_abd_table.csv')

    glycoprofile_list_addr = os.path.join(intermediate_dir, project_name + "_glycoprofile_list.json")
    simple_profile = False
    simple_naming = False
    external_profile_naming = False

    para_keyword = {'project_name': project_name,
                    'working_addr': working_addr,
                    'glycoct_dir': glycoct_dir,
                    'source_dir': source_dir,
                    'intermediate_dir': intermediate_dir,
                    'plot_output_dir': plot_output_dir,
                    'glycan_dict_addr': glycan_dict_addr,
                    # 'abundance_table_addr': abundance_table_addr,
                    'glycan_motif_dict_addr': glycan_motif_dict_addr,
                    'motif_dict_addr': motif_dict_addr,
                    'substructure_abd_table_addr': substructure_abd_table_addr,
                    'motif_abd_table_addr': motif_abd_table_addr,
                    'matched_dict_addr': matched_dict_addr,
                    'external_profile_naming_addr': external_profile_naming_addr,
                    'simple_profile': simple_profile,
                    'simple_naming': simple_naming,
                    'external_profile_naming': external_profile_naming,
                    'name_to_id_addr': name_to_id_addr,
                    'glycoprofile_list_addr': glycoprofile_list_addr
                    }
    for key, para in kwargs.items():
        para_keyword[key] = para
    return para_keyword


def _create_files(file_dir_list):
    """
    create the files needed
    :param keyword_dict:
    :return:
    """
    # exact_Ture = True
    # intermediate_address = keyword_dict['intermediate_dir']
    # plot_output_address = keyword_dict['plot_output_address']
    # source_address = keyword_dict['source_address']
    # source_address = keyword_dict['source_address']

    """create all dir for files"""
    for i in file_dir_list:
        if not os.path.isdir(i):
            os.mkdir(i)
            print("created", i)


def _find_dir(keywords_dict):
    checked_list = []
    for i in keywords_dict.keys():
        if i.find('_dir') != -1:
            checked_list.append(keywords_dict[i])
    return checked_list


def _check_exist(checked_list):
    """check all dir"""
    print("Check if the required files exist")
    for i in checked_list:
        if not os.path.isdir(i):
            print("files doesn't exist", i)
            return False
    print("files checked")
    return True


def check_init_dir(keywords_dict):
    """
    check if the working directory exists, if new ask to transfer the

    :param keywords_dict:
    :return:
    """
    working_addr = keywords_dict['working_addr']
    if os.path.isdir(working_addr):
        pass
    else:
        # print(working_addr)
        os.mkdir(working_addr)
        print("Successfully created the directory %s " % working_addr)

    check_dir_list = _find_dir(keywords_dict)
    # print(check_dir_list)
    if not _check_exist(check_dir_list):
        _create_files(check_dir_list)
    assert _check_exist(check_dir_list), 'files created unsuccessfully'
    print("Successfully created the directory need, please add the source file is the directory")


    # break
    # __init__.intermediate_address = keywords_dict["intermediate_address"]
    # __init__.plot_output_address = keywords_dict['plot_output_address']
    # __init__.source_address = keywords_dict['source_address']
    # print("set", __init__.plot_output_address)
    # print("set", __init__.json_address)
    # print("set", __init__.source_address)

    # print("set created the directory need")
    # print("set", keywords_dict["intermediate_address"])
    # print("set", keywords_dict['plot_output_address'])
    # print("set", keywords_dict['source_address'])


#
# def check_glycan_motif_dict(a_motif_dic):
#     for i in a_motif_dic:
#         for j in a_motif_dic[i]:
#             assert isinstance(j, glypy.Glycan)


def load_glycans_pip(keywords_dict, data_type, structure_loader):
    # project_name, structure_loader, data_type, glytoucan_db="", glycoct_address=""):
    """
    :param keywords_dict:
    :param data_type: one in [used, glycan_dict, glytoucanid, local_glycoct, mix]
    :param structure_loader: a list of glycan name/customized id, or a glycan_dict

    load glycan from several type of data and extra parameter may required from keyword_dict:
    if used:
        -
    if glycan_dict:
        -
    if glytoucanid:
       :param glytoucan_database_addr
    if local_glycoct:
        :param glycoct_dir
    if mix:
        :param glycoct_dir & glytoucan_database_addr
    :return:
    """
    try:
        check_init_dir(keywords_dict)
        # project_name = kwargs['project_name']
        if data_type == "glycan_dict":
            glycan_dict = structure_loader
            assert glycan_io.check_glycan_dict(glycan_dict), "Wrong structure_loader"

        elif data_type == "glytoucanid":
            glycan_dict = {}
            if 'glytoucan_db_addr' not in keywords_dict.keys():
                assert False, 'need glytoucan_db_addr'
            glytoucan_db = glycan_io.load_glytoucan_database(keywords_dict['glytoucan_db_addr'])
            for i in structure_loader:
                _re = glycan_io.load_glycan_obj_from_glytoucan(i, glytoucan_db)
                if _re:
                    glycan_dict[i] = _re

        elif data_type == "local_glycoct":
            if 'glycoct_dir' not in keywords_dict.keys():
                assert False, 'need glycoct_dir'
            glycoct_dir = keywords_dict['glycoct_dir']
            glycan_dict = {}
            _glycan_dict = glycan_io.load_glycan_obj_from_glycoct_file(glycoct_dir)
            assert type(structure_loader) is list, 'structure_loader should be a list of glycan_id'
            try:
                for j in structure_loader:
                    _j = j
                    glycan_dict[j] = _glycan_dict[j]
            except KeyError:
                print(_j, 'cannt find it in local dir')
            print('end loading glycoct from ', glycoct_dir)

        elif data_type == "used":
            glycan_dict = glycan_io.load_glycan_dict_from_json(keywords_dict['glycan_dict_addr'])

        elif data_type == "mix":

            if 'glycoct_dir' not in keywords_dict.keys() or 'glytoucan_db_addr' not in keywords_dict.keys():
                assert False, 'need glycoct_dir and glytoucan_db_addr'
            # print('isany', kwargs['glycoct_dir'])
            glycoct_dir = keywords_dict['glycoct_dir']
            glycan_dict = {}
            _glycoct_glycan_dict = glycan_io.load_glycan_obj_from_glycoct_file(glycoct_dir)
            print('end loading glycoct from ', glycoct_dir)
            # print(glycan_dict)
            glytoucan_db = glycan_io.load_glytoucan_database(keywords_dict['glytoucan_db_addr'])
            for i in structure_loader:
                if i not in _glycoct_glycan_dict.keys():
                    _re = glycan_io.load_glycan_obj_from_glytoucan(i, glytoucan_db)
                    if _re:
                        glycan_dict[i] = _re
                    else:
                        assert False, 'missing glycan ' + i
                else:
                    glycan_dict[i] = _glycoct_glycan_dict[i]
            print('end loading glytoucan db total', len(glycan_dict.keys()), 'are loaded')
        else:
            assert False, 'the Wrong type'
        json_utility.store_json(keywords_dict['glycan_dict_addr'],
                                glycan_io.glycan_obj_to_glycan_str(glycan_dict))
        print(keywords_dict['glycan_dict_addr'])
        return glycan_dict
    except KeyError:
        print("No such glycan")


def extract_and_merge_substrutures_pip(keywords_dict, linkage_specific, forced=False):
    # project_name = keywords_dict['project_name']
    # working_addr = keywords_dict['working_addr']
    # intermediate_address = keywords_dict['intermediate_address']
    # source_address = keywords_dict['source_address']
    glycan_dict_addr = keywords_dict['glycan_dict_addr']
    # print(glycan_dict_addr)
    glycan_motif_dict_addr = keywords_dict['glycan_motif_dict_addr']
    motif_dict_addr = keywords_dict['motif_dict_addr']
    matched_dict_addr = keywords_dict['matched_dict_addr']

    if os.path.isfile(glycan_dict_addr):
        print('start glycan_dict')
        if forced or not os.path.isfile(matched_dict_addr):
            if forced or not os.path.isfile(glycan_motif_dict_addr):
                    glycan_dict = glycan_io.load_glycan_dict_from_json(glycan_dict_addr)
                    glycan_motif_dic = extract_substructures.extract_substructures_pip(glycan_dict=glycan_dict,
                                                                                       gly_len=23,
                                                                                       output_file=glycan_motif_dict_addr)
                    print('finished parse motif_dic')
            else:
                    glycan_motif_dic = glycan_io.load_glycan_motif_dict_from_json(glycan_motif_dict_addr)
                    glycan_dict = glycan_io.load_glycan_dict_from_json(glycan_dict_addr)
                    print('loaded existed motif_dic')
            if forced or not os.path.isfile(motif_dict_addr):
                print('start merge motif_dict')
                merge_motif_dict = merge_substructure_vec.merge_substructure_dict_pip(
                    glycan_motif_dict=glycan_motif_dic,
                    glycan_dict=glycan_dict,
                    linkage_specific=linkage_specific,
                    output_merged_motif_dict_addr=motif_dict_addr)
                print('finished merge motif_dic')
            else:
                merge_motif_dict = glycan_io.glycan_str_to_glycan_obj(glycan_io.load_json(motif_dict_addr))
                print('loaded merged motif_dic')

            matched_dict = merge_substructure_vec.substructure_matching_wrapper(motif_dict=merge_motif_dict,
                                                                                glycan_motif_dict=glycan_motif_dic,
                                                                                linkage_specific=linkage_specific,
                                                                                matched_glycan_dict_addr=matched_dict_addr)
        else:
            matched_dict = glycan_io.load_json(matched_dict_addr)
        print('finished glycan deconvolution')
    else:
        assert False, 'cannot find the glycan_dict file'
    return matched_dict


def glycoprofile_pip(keywords_dict, abd_table, unique_glycan_identifier_to_glytoucan_id=False,
                     already_glytoucan_id=False,
                     external_profile_naming=False, forced=False, ):
    """
    required file
    :param keywords_dict:
    :param unique_glycan_identifier_to_glytoucan_id:
    :param external_profile_naming:
    :param already_glytoucan_id:
    :param forced:
    :return:
    """
    name_to_id_addr = keywords_dict['name_to_id_addr']

    naming_abd_dict, profile_columns = glycan_io.abd_table_to_dict(abd_table)

    """generating the glycoprofile naming"""

    matched_dict_addr = keywords_dict['matched_dict_addr']
    external_profile_naming_addr = keywords_dict['external_profile_naming_addr']
    # motif_abd_table_addr = keywords_dict['motif_abd_table_addr']
    substructure_abd_table_addr = keywords_dict['substructure_abd_table_addr']
    glycoprofile_list_addr = keywords_dict['glycoprofile_list_addr']

    # if forced or not os.path.isfile(substructure_abd_table_addr):
    if os.path.isfile(matched_dict_addr):
        if unique_glycan_identifier_to_glytoucan_id:
            if already_glytoucan_id:
                """the easiest way just duplicate everything"""
                glycan_identifier_to_glytoucan_id = {}
                naming = list(naming_abd_dict.keys())
                for i in profile_columns:
                    glycan_identifier_to_glytoucan_id[i] = dict(zip(naming, naming))
                    # print(glycan_identifier_to_glytoucan_id)
            else:
                name_to_id_addr = keywords_dict['name_to_id_addr']
                name_to_id = glycan_io.load_json(name_to_id_addr)
                _ = list(name_to_id.keys())
                if type(_[0]) == dict:
                    glycan_identifier_to_glytoucan_id = glycan_io.load_json(name_to_id_addr)
                else:
                    glycan_identifier_to_glytoucan_id = {}
                    for i in profile_columns:
                        glycan_identifier_to_glytoucan_id[i] = name_to_id
        else:
            if os.path.isfile(name_to_id_addr):
                glycan_identifier_to_glytoucan_id = glycan_io.load_glycoprofile_name_to_id(name_to_id_addr)
            else:
                assert False, 'missing one of them' + name_to_id_addr

        if external_profile_naming:
            if os.path.isfile(external_profile_naming_addr):
                profile_name = json_utility.load_json(external_profile_naming_addr)
            else:
                print("no external profile naming found")
                profile_name = []
        else:
            print("no external profile naming found")
            profile_name = []

        match_dict = glycan_io.load_match_dict_from_json(matched_dict_addr)
        glycoprofile_list = process_glycoprofiles.get_glycoprofile_list(glycan_identifier_to_glytoucan_id,
                                                                        naming_abd_dict,
                                                                        match_dict,
                                                                        profile_columns,
                                                                        profile_name,
                                                                        glycoprofile_list_addr,
                                                                        get_existance=True)
        table_generator = process_glycoprofiles.MotifAbdTableGenerator(glycoprofile_list)
        glycoprofile_vector_table = table_generator.table_against_wt_relative_abd()
        glycoprofile_vector_table.to_csv(substructure_abd_table_addr)
    else:
        assert False, 'missing one of them' + \
                      '\n'.join([matched_dict_addr,
                                 external_profile_naming_addr,
                                 ]) + '\n' + \
                      '\n'.join([str(xx) for xx in [os.path.isfile(matched_dict_addr),
                                                    os.path.isfile(external_profile_naming_addr),
                                                    ]])
    # else:
    #     glycoprofile_vector_table = pd.read_csv(substructure_abd_table_addr)
    #     print('loaded substructure_abd_table')
    return glycoprofile_vector_table, glycoprofile_list


def select_motifs_pip(keywords_dict, linkage_specific, only_substructures_start_from_root, select_col=[],
                      core_substructure_index=-1):
    motif_dict_addr = keywords_dict['motif_dict_addr']
    assert os.path.isfile(motif_dict_addr), 'missing ' + motif_dict_addr
    substructure_abd_table_addr = keywords_dict['substructure_abd_table_addr']
    assert os.path.isfile(substructure_abd_table_addr), 'missing' + substructure_abd_table_addr

    substructure_abd_table = pd.read_csv(substructure_abd_table_addr, index_col=0)
    motif_dict = glycan_io.load_motif_dict_from_json(motif_dict_addr)
    # _motif_lab = select_motifs.MotifLabwithCore(motif_dict, glycan_core=select_motifs.nglycan_core, linkage_specific=False)  # unicarbkb_motifs_12259.json
    # _motif_lab.get_dependence_tree_core()

    if not select_col:
        select_col = substructure_abd_table.columns
    if not only_substructures_start_from_root:
        _motif_lab = select_motifs.MotifLab(motif_=motif_dict,
                                            linkage_specific=linkage_specific)
        _motif_lab.get_dependence_tree_all()
        a_node_state = select_motifs.NodesState(dependence_tree=_motif_lab.motif_dep_tree,
                                                motif_weight=select_motifs.get_weight_dict(
                                                    substructure_abd_table[select_col]),
                                                linkage_specific=linkage_specific)
        node_attri, edge_attri, mod_nodes, mod_edges, merged_weights_dict = a_node_state.nodes_dropping_pipe(
            drop_parellel=False, drop_diff_abund=True)

    else:

        _motif_lab = select_motifs.MotifLabwithCore(motif_=motif_dict,
                                                    glycan_core=select_motifs.nglycan_core,
                                                    linkage_specific=linkage_specific)  # unicarbkb_motifs_12259.json
        _motif_lab.get_dependence_tree_core()
        a_node_state = select_motifs.NodesState(dependence_tree=_motif_lab.motif_dep_tree_core,
                                                motif_weight=select_motifs.get_weight_dict(
                                                    substructure_abd_table[select_col]),
                                                linkage_specific=linkage_specific)
        node_attri, edge_attri, mod_nodes, mod_edges, merged_weights_dict = a_node_state.nodes_dropping_pipe(
            drop_parellel=False, drop_diff_abund=True)
        #
        # motif_dict_addr = keywords_dict['motif_dict_addr']
        # assert os.path.isfile(motif_dict_addr), 'missing ' + motif_dict_addr
        # substructure_abd_table_addr = keywords_dict['substructure_abd_table_addr']
        # assert os.path.isfile(substructure_abd_table_addr), 'missing' + substructure_abd_table_addr
        #
        # substructure_abd_table = pd.read_csv(substructure_abd_table_addr, index_col=0)
        # motif_dict = glycan_io.load_motif_dict_from_json(motif_dict_addr)
        # _motif_lab = select_motifs.MotifLabwithCore(motif_dict, glycan_core=select_motifs.nglycan_core, linkage_specific=False)  # unicarbkb_motifs_12259.json
        # _motif_lab.get_dependence_tree_core()
        #
        # _a = select_motifs.NodesState(_motif_lab.motif_dep_tree_core,
        #                                               select_motifs.get_weight_dict(substructure_abd_table[_table_col]),
        #                             linkage_specific=linkage_specific)
        #
        # _,_,mod_nodes,_, merged_weights_dict=_a.nodes_dropping_pipe(drop_parellel=False, drop_diff_abund=True)
        #
        # # _collapsed_edge, _collapsed_node, _collapsed_dege_attri = _a.collapsing_potential_node()
        # #
        if core_substructure_index != -1:
            if core_substructure_index in mod_nodes:
                mod_nodes.remove(core_substructure_index)

    motif_abd_table = substructure_abd_table[select_col][substructure_abd_table.index.isin(mod_nodes)]
    motif_abd_table_addr = keywords_dict['motif_abd_table_addr']
    motif_abd_table.to_csv(motif_abd_table_addr)

    return motif_abd_table, _motif_lab, merged_weights_dict


def clustering_analysis_pip(keywords_dict, motif_abd_table, select_profile_name=[]):
    # df_ncore = deepcopy(motif_abd_table)
    # print(sorted(mod_nodes))
    # print(df_ncore.shape)
    # draw plot
    # motif_with_n_glycan_core_all_motif(motif_, _table, weight_dict)
    """ with n_glycan_core using jaccard for binary and use braycurtis for float
    """

    import matplotlib.pyplot as plt
    if select_profile_name:
        selected_name_list = select_profile_name
        motif_abd_table.columns = selected_name_list

    else:
        selected_name_list = motif_abd_table.colmuns

    # df_ncore=pd.DataFrame(data=preprocessing.scale(df_ncore.transpose()).transpose(), index=df_ncore.index, columns=df_ncore.columns)
    # motif_abd_table.to_csv(os.path.join(keywords_dict['intermediate_dir'],
    #                                     str(len(selected_name_list)) + r"selected_abundance_matrix.txt"))
    # motif_abd_table.colmuns = selected_name_list

    plt.savefig(keywords_dict['plot_output_dir'] + 'pseudo_profile_clustering.svg')
    cluster_grid = clustering_analysis.draw_glycan_clustermap(motif_abd_table=motif_abd_table,
                                                              address=keywords_dict[
                                                                          'plot_output_dir'] + 'pseudo_profile_clustering.svg',
                                                              metric="correlation",
                                                              cmap=sns.color_palette("RdBu_r", 20),
                                                              linewidths=0.01,
                                                              figsize=(15, 15),
                                                              linecolor='black',
                                                              method='complete')
    glycoprofile_cluster_dict = clustering_analysis.draw_profile_cluster(g=cluster_grid,
                                                                         df=motif_abd_table,
                                                                         profile_name=selected_name_list,
                                                                         color_threshold=0.5,
                                                                         address=keywords_dict[
                                                                                     'plot_output_dir'] + 'profile_clustering.svg')
    glyco_motif_cluster_dict = clustering_analysis.draw_motif_cluster(g=cluster_grid,
                                                                      df=motif_abd_table,
                                                                      color_threshold=0.185,
                                                                      address=keywords_dict[
                                                                                  'plot_output_dir'] + 'motif_cluster.svg',
                                                                      fig_size=(6, 20),
                                                                      )

    # raw_abd_zscore_plot_addr = keywords_dict['raw_abd_zscore_plot_addr']
    # cccluster_dict = clustering_analysis.draw_motif_cluster(g, motif_abd_table, name_prefix, color_threshold=0.185, fig_size=(6, 20))

    # plt.savefig(raw_abd_zscore_plot_addr)
    return glycoprofile_cluster_dict, glyco_motif_cluster_dict


    # (glyco_motif_cluster=glyco_motif_cluster_dict,
    #                                                      substructure_vec=a_node_state.,
    #                                                      motif_weights_dict=,
    #                                                      address_dir=,
    #                                                      threshold=0.51,
    #                                                      plot_rep=True)


def parse_abundance_table(glycan_abundance_table, ):
    """
    load glycan table
    load motif vec
    :param glycan_table:
            column glycan
            index glycoprofile
    :param index=0
    :param
    :return:
    """
    sorted_glycan_data_array = np.zeros(glycan_abundance_table.shape)
    filtered_matrix = np.array(glycan_abundance_table)[:, :]
    for i in range(glycan_abundance_table.shape[0]):
        sorted_glycan_data_array[i, :] = filtered_matrix[i, :] / sum(filtered_matrix[i, :])


import pandas as pd


def parse_meta_table(addr, *kwargs):
    """parse the glycan meta table into Glycoprofile and store it"""
    pd.read_csv(addr, sep='\t', header=True)

    pass


def parse_with_name(name_dict):
    pass
