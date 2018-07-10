from json_utility import load_json, store_json
import __init__
from glypy.io import glycoct, iupac
from pathlib import Path


def load_glycan_dict_from_json(addr=__init__.glycan_dict_addr, loader=glycoct):
    """
    :param loader: glycoct or iupac
    :return: a glycan_dict
    """
    # if prior:
    glycan_str_dict = load_json(addr)
    glycan_dict = {}
    for i in glycan_str_dict:
        glycan_dict[i] = loader.loads(glycan_dict[i])
    return glycan_dict


def load_glycan_str_from_database(topology_list_addr, output_file=__init__.glycan_dict_addr, loader=glycoct):
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
            _glycan_str = load_glycan_str_from_manual_drawn(glycan_id)
        if _glycan_str == "":
            pass
        else:
            _loaded += 1
            glycan_str_dict[glycan_id] = _glycan_str
            glycan_dict[glycan_id] = loader.loads(_glycan_str)
    print("There are ", _count, "glycan id found; ", _loaded, "glycans loaded")
    store_json(output_file, glycan_str_dict)
    return glycan_dict


def motif_dict_to_motif_vec(motif_dict):
    motif_vec = []
    for i in sorted([int(j) for j in list(motif_dict.keys())]):
        print(i, len(motif_dict[str(i)]))
        motif_vec.extend(motif_dict[str(i)])
    print(len(motif_vec))
    return motif_vec


def glycan_str_to_glycan(a_dict_of_glycan_str):
    """
    :param a_dict_of_glycan_str: load glycoct str transform to Glycan
    :return:
    """
    if type(a_dict_of_glycan_str) == list:
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


def load_glytoucan_database():
    return load_json(__init__.glytoucan_database_addr__)


def load_glycan_str_from_manual_drawn(glycan_id):
    try:
        my_file = Path(__init__.source_address + glycan_id + ".glycoct_condensed")
        if my_file.exists():
            f = open(__init__.source_address + glycan_id + ".glycoct_condensed")
        else:
            if Path(__init__.source_address + glycan_id).exists():
                f = open(__init__.source_address + glycan_id)
            else:
                raise FileNotFoundError
        glycan_str = "".join(f.readlines())

    except FileNotFoundError:
        print("This id: ", glycan_id, " not found")
        glycan_str = ''
    return glycan_str

def motif_dict_to_motif_vec(motif_dict):
    motif_vec = []
    for i in sorted([int(j) for j in list(motif_dict.keys())]):
        print(i, len(motif_dict[str(i)]))
        motif_vec.extend(motif_dict[str(i)])
    print(len(motif_vec))
    return motif_vec