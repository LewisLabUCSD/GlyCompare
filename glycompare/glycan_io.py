import glypy
from glypy import Glycan
from json_utility import load_json, store_json
import __init__
from glypy.io import glycoct, iupac
from pathlib import Path
import os
import pandas as pd
import numpy as np


def load_glycoprofile_name_to_id(addr, naming='mz', _format='json'):
    """ glycan_profile = {'1': {'2244': 'G04483SK',
                                '2605': 'G30460NZ',
                                '2967': 'G17689DH',
                                '3416': 'G54338PJ',
                                '3865': '3865.1',
                                '4226': 'G86696LV',
                                '4587': '4587.1',
                                '5037': 'G49604DB',
                                '5486': '5486.1'}}
    load such format
                                """
    if _format == 'gff':
        # f = open(addr)
        profile_naming_dict = {}
        # for i in f.readline():
        #     profile, name, glycan_id = i.rstrip('\n').split('\t')
        #     if profile not in profile_naming_dict.keys():
        #         profile_naming_dict[profile] = {name: glycan_id}
        #     else:
        #         assert name not in profile_naming_dict[profile].keys()
        #         profile_naming_dict[profile][name] = glycan_id
        f = pd.read_csv(addr, sep='\t', header=0)
        # naming='mz'
        naming_list = f[naming].tolist()
        glycan_id_list = f['glycan_id'].tolist()
        profile_list = f['profile'].tolist()

        for index, profile in enumerate(profile_list):
            if profile not in profile_naming_dict.keys():
                profile_naming_dict[profile] = {naming_list[index]: glycan_id_list[index]}
            else:
                assert naming_list[index] not in profile_naming_dict[profile].keys()
                profile_naming_dict[profile][naming_list[index]] = glycan_id_list[index]
    elif _format == 'json':
        profile_naming_dict = load_json(addr)
    else:
        assert False, "wrong format"
    return profile_naming_dict


def output_glycoprofile_name_to_id(addr, glycoprofile_naming_dict, _formate='json'):
    if _formate == 'gff':
        file = open(addr, 'w')
        file.write('\t'.join(['profile', 'mz', 'glycan_id']) + "\n")
        for i in glycoprofile_naming_dict:
            for j in glycoprofile_naming_dict[i]:
                file.write('\t'.join([i, j, glycoprofile_naming_dict[i][j]]) + "\n")
        file.close()
    elif _formate == "json":
        store_json(addr, glycoprofile_naming_dict)


def load_glycan_dict_from_json(addr):
    """

    :param addr:
    :return: a glycan_dict
    """
    # if prior:
    a = load_json(addr)
    return glycan_str_to_glycan_obj(a)
    # return glycan_dict


def load_motif_dict_from_json(addr):
    """

    :param addr:
    :return: a glycan_dict
    """
    # if prior:
    a = load_json(addr)
    return glycan_str_to_glycan_obj(a)


def load_glycan_motif_dict_from_json(addr):
    """

    :param addr:
    :return: a glycan dict
    """
    a = load_json(addr)
    return glycan_str_to_glycan_obj(a)


def load_match_dict_from_json(addr):
    """
    :param addr:
    :return: a glycan dict
    """
    return load_json(addr)


def load_glycan_obj_from_database(topology_list_addr, output_file="", loader=glycoct):
    """
        each line in topology_list is defined as m/z, glycan_id

        1. get the glycanID from Glycan_topolog_list
        2. find glycan_id in glytoucan database: /root_address + r'data_dic_finnn.json'
        3. find glycan_id in self-generated local file: /__init__.json_address+glycan_id+".glycoct_condensed"
        4. output a dict glycan_id str -> glycan_str str stored in: root_address + 'NBT_for_motif_extraction.json'
        :return: glycan_dict
        """
    x = load_glytoucan_database()
    glycan_str_dict = {}
    glycan_dict = {}
    f = open(topology_list_addr)
    _count = 0
    _loaded = 0
    for i in f.readlines():
        _count += 1
        _name, glycan_id = i.rstrip('\n').split("\t")

        if len(glycan_id) == 8:
            _glycan_str = load_glycan_str_from_glytoucan(glycan_id, x)
        else:
            _glycan_str = load_glycan_str_from_glycoct(glycan_id)
        if _glycan_str == "":
            pass
        else:
            _loaded += 1
            glycan_str_dict[glycan_id] = _glycan_str
            glycan_dict[glycan_id] = loader.loads(_glycan_str)
    print("There are ", _count, "glycan id found; ", _loaded, "glycans loaded")
    if output_file == "":
        pass
    else:
        store_json(output_file, glycan_str_dict)
    return glycan_dict


def output_dict_to_glycoct(dict, addr):
    for i, j in dict.items():
        out_glycan_obj_as_glycoct(j,i,addr)

def load_glycan_obj_from_glycoct_file(dir_address):
    """
        load all Glycan objs from a directory
        :return: glycan_dict
        """
    glycan_dict = {}
    for (dirpath, dirnames, filenames) in os.walk(dir_address):
        for i in filenames:
            if i.find('.glycoct_condensed') != -1:
                glycan_id = i[:i.find('.glycoct_condensed')]
                print(glycan_id)
                temp_i = load_glycan_str_from_glycoct(glycan_id, dir_address)
                glycan_dict[glycan_id] = glycoct.loads(temp_i)
    return glycan_dict


def glycan_str_to_glycan_obj(a_dict_of_glycan_str):
    """
    :param a_dict_of_glycan_str: load glycoct str transform to Glycan
    :return:
    """
    if type(a_dict_of_glycan_str) == list:
        # if a list
        return [glycoct.loads(i) for i in a_dict_of_glycan_str]
    elif type(a_dict_of_glycan_str) == dict:
        a_dict = {}
        for i in a_dict_of_glycan_str.keys():
            if type(a_dict_of_glycan_str[i]) == dict:
                a_dict[i] = {}
                for j in a_dict_of_glycan_str[i].keys():
                    a_dict[i][j] = [glycoct.loads(k) for k in a_dict_of_glycan_str[i][j]]
            elif type(a_dict_of_glycan_str[i]) == list:
                a_dict[i] = [glycoct.loads(k) for k in a_dict_of_glycan_str[i]]
            elif type(a_dict_of_glycan_str[i]) == str:
                a_dict[i] = glycoct.loads(a_dict_of_glycan_str[i])
        return a_dict


def load_table_to_dict(addr):
    return abd_table_to_dict(load_table(addr))


def load_table(addr, rep=None):
    """
    require the table: column is glycan, row is glycoprofile
    :param addr:
    :return:
    """
    _format = addr.split('.')[-1]
    if _format == 'xls':
        _table = pd.read_excel(addr, index_col=0)
    elif _format in ['csv', 'txt']:
        f = open(addr, 'r')
        f.readline()
        a = f.readline()
        f.close()
        if a.find(',') != -1:
            _table = pd.read_csv(addr, index_col=0, sep=',')
        elif a.find('\t') != -1:
            _table = pd.read_csv(addr, index_col=0, sep='\t')
        else:
            assert False, 'Error in load_table_to_dict, sep is not , or tab'
    else:
        print('format is unrecognizable, will try pd.read_csv')
        f = open(addr, 'r')
        f.readline()
        a = f.readline()
        f.close()
        if a.find(',') != -1:
            _table = pd.read_csv(addr, index_col=0, sep=',')
        elif a.find('\t') != -1:
            _table = pd.read_csv(addr, index_col=0, sep='\t')
        else:
            assert False, 'Error in load_table_to_dict, sep is neither , nor tab'
    if rep:
        return _table.replace(rep, 0)
    else:
        return _table


def abd_table_to_dict(abd_table):
    """
    parse mz, hplc or glycan table to dict
    :param abd_table:
    :return: abd_dict
    """
    mz_list = list(abd_table.index)
    profile_list = list(abd_table.columns.values)
    print('mz_list', len(mz_list), 'profile_list', len(profile_list))
    lr_ = abd_table.shape[0]
    lc_ = abd_table.shape[1]
    print(lr_, lc_)

    # filtered_matrix = np.array(mz_abd_table)[:, 1:]
    # normalized_matrix = np.array(filtered_matrix)
    # for i in range(lc_):
    #     normalized_matrix[:, i] = filtered_matrix[:, i] / sum(filtered_matrix[:, i])
    # print(normal_weight_list)
    _matrix = np.array(abd_table)
    mz_abd_dict = {}
    for idex, i in enumerate(mz_list):
        mz_abd_dict[str(i)] = list(_matrix[idex, :][:])
    # store_json(norm_abd_table_dict_addr, mz_abd_dict)

    return mz_abd_dict, [str(i) for i in abd_table.columns.tolist()]


def glycan_obj_to_glycan_str(a_dict_to_glycan_str):
    """
    :param a_dict_to_glycan_str: transform the glypy.Glycan to str
    :return: a sting dict
    """
    if type(a_dict_to_glycan_str) == list:
        # if a list
        return [str(i) for i in a_dict_to_glycan_str]
    elif type(a_dict_to_glycan_str) == dict:
        a_dict = {}
        for i in a_dict_to_glycan_str.keys():
            if type(a_dict_to_glycan_str[i]) == dict:
                a_dict[i] = {}
                for j in a_dict_to_glycan_str[i].keys():
                    a_dict[i][j] = [str(k) for k in a_dict_to_glycan_str[i][j]]
            elif type(a_dict_to_glycan_str[i]) == list:
                a_dict[i] = [str(k) for k in a_dict_to_glycan_str[i]]
            elif isinstance(a_dict_to_glycan_str[i], glypy.Glycan):
                a_dict[i] = str(a_dict_to_glycan_str[i])
        return a_dict


def load_glycan_str_from_glytoucan(glycan_id, glytoucan):
    """
    Given a glytoucan id, check if they have glycoct format structure, if so, return the str, if not return ''
    :param glycan_id: glytoucan id
    :param glytoucan: dict, glytoucan database
    :return:
    """
    try:
        _gly_stru = glytoucan[glycan_id]['structure_']
    except KeyError:
        print("no name: ", glycan_id, "Please manually draw the glycoct format")
        _gly_stru = ''
    return _gly_stru


def load_glycan_obj_from_glytoucan(glycan_id, glytoucan):
    _gly_stu = load_glycan_str_from_glytoucan(glycan_id, glytoucan)
    if _gly_stu == '':
        print('missing structure', glycan_id)
        return None
    else:
        print('adding', glycan_id)
        return glycoct.loads(_gly_stu)


def load_glycan_id_from_meta_data(addr):
    """load glycan id list with head 'glycan_id"""
    glycan_id = pd.read_csv(addr, sep='\t')['glycan_id'].tolist()
    return glycan_id


def check_glycan_dict(glycan_dict):
    """
        check if the glycan dict contains Glycan obj not string
    """
    for i in glycan_dict.keys():
        assert isinstance(glycan_dict[i], Glycan), 'Failed glycan dict check, should be glypy.Glycan'


def check_glycan_motif_dict(a_glycan_motif_dict):
    for i in a_glycan_motif_dict:
        assert type(i) == str, 'glycan name is str'
        for j in a_glycan_motif_dict[i]:
            assert type(j) == str, ('glycan degree are not str', j, type(j))
            assert isinstance(a_glycan_motif_dict[i][j][0], glypy.Glycan)


def check_motif_dict(a_motif_dict):
    assert str(1) in a_motif_dict.keys(), 'a motif_dict without monossar'
    for j in a_motif_dict.keys():
        assert type(j) == str, ('glycan degree are stored in degree, it stores', j, type(j))
        assert isinstance(a_motif_dict[j][0], glypy.Glycan)


def load_glytoucan_database(addr):
    print('loading glytoucan_database from ', addr)
    return load_json(addr)


def load_glycan_obj_from_glycoct(glycan_id, address):
    _gly_stu = load_glycan_str_from_glycoct(glycan_id, address)
    if _gly_stu == '':
        return None
    else:
        return glycoct.loads(_gly_stu)


def load_glycan_str_from_glycoct(glycan_id, address):
    try:
        if Path(os.path.join(address, glycan_id)).exists():
            f = open(os.path.join(address, glycan_id))
        elif Path(os.path.join(address, glycan_id + '.glycoct_condensed')).exists():
            f = open(os.path.join(address, glycan_id + '.glycoct_condensed'))
        else:
            raise FileNotFoundError
        glycan_str = "".join(f.readlines())

    except FileNotFoundError:
        print("This id: ", glycan_id, " not found, in ", address)
        glycan_str = ''
    return glycan_str


def motif_dict_to_motif_vec(motif_dict):
    """convert the motif dictionary to vector"""
    motif_vec = []
    for i in sorted([int(j) for j in list(motif_dict.keys())]):
        print(i, len(motif_dict[str(i)]))
        motif_vec.extend(motif_dict[str(i)])
    print(len(motif_vec))
    return motif_vec


def motif_vec_to_motif_dict(motif_vec):
    """convert the motif vector to dictionary with degree as key"""
    motif_dict = {}
    for i in motif_vec:
        # print(i,))
        if len(i) not in motif_dict.keys():
            motif_dict[len(i)] = [i]
        else:
            motif_dict[len(i)].append(i)
    # print(len(motif_vec))
    return motif_dict


def out_glycan_obj_as_glycoct(a_glycan, glycan_id, dir_addr):
    if dir_addr.find(glycan_id + '.glycoct_condensed') != -1:
        print('file already exist, will not store')
    else:
        glycan_addr = os.path.join(dir_addr, glycan_id + '.glycoct_condensed')
        _w = open(glycan_addr, 'w')
        _w.write(str(a_glycan))
        _w.close()
