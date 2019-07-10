import os

import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import __init__
import customize_motif_vec
import extract_substructure
import glycan_io
import glycan_profile
import json_utility
import motif_class


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
    glycoct_dir = os.path.join(working_addr, 'glycoct/')

    name_to_id_addr = os.path.join(source_dir, 'glycoprofile_name_to_glycan_id.json')
    # abundance_table_addr = os.path.join(source_dir, 'abundance_table')
    external_profile_naming_addr = os.path.join(source_dir, 'external_profile_naming.json')

    glycan_dict_addr = os.path.join(intermediate_dir, project_name + '_glycan_dict.json')
    glycan_motif_dict_addr = os.path.join(intermediate_dir, project_name + '_glycan_motif_dict.json')
    motif_dict_addr = os.path.join(intermediate_dir, project_name + "_motif_dict.json")
    matched_dict_addr = os.path.join(intermediate_dir, project_name + "_match_dict.json")
    motif_abd_table_addr = os.path.join(intermediate_dir, project_name + "_motif_abd_table.csv")
    substructure_abd_table_addr = os.path.join(intermediate_dir, project_name + '_substructure_abd_table.csv')
    raw_abd_zscore_plot_addr = os.path.join(plot_output_dir, 'raw_abundance_zscore.eps')
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
                    'raw_abd_zscore_plot_addr': raw_abd_zscore_plot_addr,
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


def load_structure_pip(keywords_dict, data_type, structure_loader):
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


def glycan_deconvoluting_pip(keywords_dict, linkage_specific, forced=False):
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
        if forced or not os.path.isfile(glycan_motif_dict_addr):
            glycan_dict = glycan_io.load_glycan_dict_from_json(glycan_dict_addr)
            glycan_motif_dic = extract_substructure.extract_substructure_pip(glycan_dict=glycan_dict, gly_len=23,
                                                                             output_file=glycan_motif_dict_addr)
            print('finished parse motif_dic')
        else:
            glycan_motif_dic = glycan_io.load_glycan_motif_dict_from_json(glycan_motif_dict_addr)
            glycan_dict = glycan_io.load_glycan_dict_from_json(glycan_dict_addr)
            print('loaded existed motif_dic')
        print('start merge motif_dict')
        if forced or not os.path.isfile(motif_dict_addr):
            merge_motif_dict = customize_motif_vec.merge_motif_dict_pipe(glycan_motif_dict=glycan_motif_dic,
                                                                         glycan_dict=glycan_dict,
                                                                         linkage_specific= linkage_specific,
                                                                         output_merged_motif_dict_addr=motif_dict_addr)
            print('finished merge motif_dic')

        else:
            merge_motif_dict = glycan_io.glycan_str_to_glycan_obj(glycan_io.load_json(motif_dict_addr))
            print('loaded merged motif_dic')

        if forced or not os.path.isfile(matched_dict_addr):
            motif_occurance_vector_dict = customize_motif_vec.motif_matching_wrapper(motif_dict=merge_motif_dict,
                                                                                     glycan_motif_dict=glycan_motif_dic,
                                                                                     linkage_specific=linkage_specific,
                                                                                     matched_glycan_dict_addr=matched_dict_addr)
        print('finished glycan deconvolution')
    else:
        print('cannot find the glycan_dict file')


def abd_table_pip(keywords_dict, abd_table, simple_profile=False, simple_naming=False,
                  external_profile_naming=False, forced=False, ):
    """
    required file
    :param keywords_dict:
    :param simple_profile:
    :param external_profile_naming:
    :param simple_naming:
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

    if forced or not os.path.isfile(substructure_abd_table_addr):
        if os.path.isfile(matched_dict_addr):
            if simple_profile:
                if simple_naming:
                    """the easiest way just duplicate everything"""
                    profile_naming_to_id = {}
                    naming = list(naming_abd_dict.keys())
                    for i in profile_columns:
                        profile_naming_to_id[i] = dict(zip(naming, naming))
                        # print(profile_naming_to_id)
                else:
                    name_to_id_addr = keywords_dict['name_to_id_addr']
                    name_to_id = glycan_io.load_json(name_to_id_addr)
                    _ = list(name_to_id.keys())
                    if type(_[0]) == dict:
                        profile_naming_to_id = glycan_io.load_json(name_to_id_addr)
                    else:
                        profile_naming_to_id = {}
                        for i in profile_columns:
                            profile_naming_to_id[i] = name_to_id
            else:
                if os.path.isfile(name_to_id_addr):
                    profile_naming_to_id = glycan_io.load_glycoprofile_name_to_id(name_to_id_addr)
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
            glycoprofile_list = glycan_profile.get_glycoprofile_list(profile_naming_to_id,
                                                                     naming_abd_dict,
                                                                     match_dict,
                                                                     profile_columns,
                                                                     profile_name,
                                                                     glycoprofile_list_addr,
                                                                     get_existance=True)
            table_generator = glycan_profile.MotifAbdTableGenerator(glycoprofile_list)
            substructure_abd_table = table_generator.table_against_wt_relative_abd()
            substructure_abd_table.to_csv(substructure_abd_table_addr)
        else:
            assert False, 'missing one of them' + \
                          '\n'.join([matched_dict_addr,
                                     external_profile_naming_addr,
                                     ]) + '\n' + \
                          '\n'.join([str(xx) for xx in [os.path.isfile(matched_dict_addr),
                                                        os.path.isfile(external_profile_naming_addr),
                                                        ]])
    else:
        substructure_abd_table = pd.read_csv(substructure_abd_table_addr)
        print('loaded substructure_abd_table')
    # return glycoprofile_list


def motif_vec_pip(keywords_dict, linkage_specific, epitope=False):
    motif_dict_addr = keywords_dict['motif_dict_addr']
    assert os.path.isfile(motif_dict_addr), 'missing ' + motif_dict_addr
    substructure_abd_table_addr = keywords_dict['substructure_abd_table_addr']
    assert os.path.isfile(substructure_abd_table_addr), 'missing' + substructure_abd_table_addr

    substructure_abd_table = pd.read_csv(substructure_abd_table_addr)
    motif_dict = glycan_io.load_motif_dict_from_json(motif_dict_addr)

    if epitope:
        _motif_lab = motif_class.MotifLab(motif_dict, linkage_specific)
        _motif_lab.get_dependence_tree_all()
        a_node_state = motif_class.NodesState(dependence_tree=_motif_lab.motif_dep_tree,
                                              motif_weight=motif_class.get_weight_dict(substructure_abd_table),
                                              linkage_specific=linkage_specific)
    else:
        _motif_lab = motif_class.MotifLabwithCore(motif_dict, linkage_specific)  # unicarbkb_motifs_12259.json
        _motif_lab.get_dependence_tree_core()
        a_node_state = motif_class.NodesState(dependence_tree=_motif_lab.motif_dep_tree_core,
                                              motif_weight=motif_class.get_weight_dict(substructure_abd_table),
                                              linkage_specific=linkage_specific)
    return a_node_state

    # return None


def glyco_motif_and_plotting(keywords_dict, a_node_state):
    substructure_abd_table_addr = keywords_dict['substructure_abd_table_addr']
    motif_abd_table_addr = keywords_dict['motif_abd_table_addr']

    node_attri, edge_attri, mod_nodes, mod_edges, merged_weights_dict = a_node_state.nodes_dropping_pipe(
        drop_diff_abund=True)
    # _collapsed_edge, _collapsed_node, _collapsed_dege_attri = _a.collapsing_potential_node()
    substructure_abd_table = pd.read_csv(substructure_abd_table_addr)
    _table_using = substructure_abd_table[substructure_abd_table.index.isin(mod_nodes)]
    _table_using.to_csv(motif_abd_table_addr)
    g = sns.clustermap(_table_using, metric="correlation", cmap=sns.color_palette("RdBu_r", 20), z_score=0,
                       linewidths=0.01, figsize=(18, 9), linecolor='black', method='complete', )

    raw_abd_zscore_plot_addr = keywords_dict['raw_abd_zscore_plot_addr']
    plt.savefig(raw_abd_zscore_plot_addr)


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
