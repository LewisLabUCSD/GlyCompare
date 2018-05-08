from glypy.io import glycoct, iupac
from glypy.algorithms.subtree_search import subtree_of
# from glypy.plot import plot
from glypy.structure.glycan import fragment_to_substructure
# from glypy.io.glycoct import dump
# from matplotlib import pyplot as plt
# from parser_glytoucan_utility import *
import time
import multiprocessing
import json
import gc_init
from gc_init import *


def store_json(address, dic):
    with open(address, 'w') as fp:
        json.dump(dic, fp)


def load_json(address):
    with open(address, 'r') as f:
        return json.load(f)


def extract_motif(glycoct_obj, idex=0):
    # print('start getmotif')
    _frag_motif_list = {}
    # fig, axes = plt.subplots(6,9)
    # fig.set_size_inches(14,6)
    start_time = time.time()
    for i in glycoct_obj.fragments(max_cleavages=5):
        _frag_gly = fragment_to_substructure(i, glycoct_obj)
        # plot(_frag_gly)
        if not len(_frag_gly) in _frag_motif_list.keys():
            _frag_motif_list[len(_frag_gly)] = [glycoct.loads(str(_frag_gly))]
        else:
            _frag_motif_list[len(_frag_gly)].append(glycoct.loads(str(_frag_gly)))
    mid_time = time.time()
    # print('start clean duplicate')
    _frag_motif_list = clean_duplicate(_frag_motif_list)
    end_time = time.time()
    print(idex, len(glycoct_obj), end_time - mid_time, mid_time - start_time)
    # print('finished getmotif')

    return _frag_motif_list


def clean_duplicate(_frag_motif_dict):
    for i in _frag_motif_dict.keys():
        # print(i)
        ldex = 0
        _check_list = _frag_motif_dict[i]
        while ldex < len(_check_list):
            jdex = ldex + 1
            while jdex < len(_check_list):
                if subtree_of(_check_list[ldex], _check_list[jdex]) == 1:
                    del _check_list[jdex]
                else:
                    jdex += 1
            ldex += 1
        _frag_motif_dict[i] = _check_list
    return _frag_motif_dict


#
# def clean_duplicate_with_root(_frag_motif_dict):
#     for i in _frag_motif_dict.keys():
#         # print(i)
#         ldex = 0
#         _check_list = _frag_motif_dict[i]
#         while ldex < len(_check_list):
#             jdex = ldex + 1
#             while jdex < len(_check_list):
#                 if subtree_of(_check_list[ldex], _check_list[jdex], True) == 1:
#                     del _check_list[jdex]
#                 else:
#                     jdex += 1
#             ldex += 1
#         _frag_motif_dict[i] = _check_list
#     return _frag_motif_dict


def extract_motif_with_count(glycoct_obj, idex=0):
    # print('start getmotif')
    _frag_motif_list = {}
    start_time = time.time()
    for i in glycoct_obj.fragments(max_cleavages=5):
        _frag_gly = fragment_to_substructure(i, glycoct_obj)

        if not len(_frag_gly) in _frag_motif_list.keys():
            _frag_motif_list[len(_frag_gly)] = [str(_frag_gly)]
        else:
            _frag_motif_list[len(_frag_gly)].append(str(_frag_gly))
    mid_time = time.time()
    # print('start clean duplicate')
    #     _frag_motif_list = clean_duplicate_with_root(_frag_motif_list)
    end_time = time.time()
    print(idex, len(glycoct_obj), end_time - mid_time, mid_time - start_time)
    # print('finished getmotif')
    return _frag_motif_list


def extract_motif_wrap(idex, _name, glycoct_obj, motif_dic, success_list):
    if not _name in success_list:
        try:
            # print('start', _name, glycoct_obj)

            _temp_motif = extract_motif(glycoct_obj, idex)
            motif_dic[_name] = _temp_motif
            # print('has_motif', motif_dic[_name])

            success_list.append(_name)
            #     for j in motif_dic.keys():
            #         print(j, len(motif_dic[j]))
        except TypeError:
            print(idex, _name, 'has error')
        except KeyboardInterrupt:
            print('break')


def extract_motif_with_count_wrap(idex, _name, glycoct_obj, motif_dic, success_list):
    if not _name in success_list:
        try:
            # print('start', _name, glycoct_obj)

            _temp_motif = extract_motif_with_count(glycoct_obj, idex)
            motif_dic[_name] = _temp_motif
            # print('has_motif', motif_dic[_name])

            success_list.append(_name)
            #     for j in motif_dic.keys():
            #         print(j, len(motif_dic[j]))
        except TypeError:
            print(idex, _name, 'has error')
        except KeyboardInterrupt:
            print('break')


def load_glycoct_for_database():
    """ 1. get the glycanID from Glycan_topolog_list
        2. find ID in glytoucan database: /root_address + r'data_dic_finnn.json'
        3. find ID in self-generated local file: /NBT_init.json_address+_code+".glycoct_condensed"
        4. output a dict ID str -> glycoct str stored in: root_address + 'BNT_for_motif_extraction.json'
        """
    x = load_json(glytoucan_data_base_addr__)

    # store_json(r'/Users/apple/PycharmProjects/data_dic_finnn.json',x)
    def get_drawed_glycan(addre):
        from glypy.io import glycoct,iupac
        f = open(addre)
        _str = "".join(f.readlines())
    #     for i in f.readlines():
        return _str
    output_for_motif = {}

    f = open(topology_list)
    glycan_dict = {}
    _count=0
    for i in f.readlines():
        _count += 1
        _name, _code = i.rstrip('\n').split("\t")
        _name = int(_name)
        if len(_code) == 8:
            try:
                _gly_stru = x[_code]['structure_']
            except KeyError:
                print("no name: ", _code)
                continue
        else:
            _addr = gc_init.manual_curated_address+_code+".glycoct_condensed"
            _gly_stru = get_drawed_glycan(_addr)

        if _name not in glycan_dict.keys():
            glycan_dict[_name] = {_code:glycoct.loads(_gly_stru)}

        else:
            glycan_dict[_name][_code] = glycoct.loads(_gly_stru)
        output_for_motif[_code] = _gly_stru
    print(_count)
    store_json(for_extraction_glycoct_dict_addr, output_for_motif)

    return output_for_motif


def get_motif_pip(gly_len, prior=True):
    """Please set the prior=True to get the data file please run the NBT_GLYCAN_preprocess file
    If prior=False, it will generate glycan motif for all glycan in glytoucan database
    1. load  dict ID -> glycoct str
    2. convert to glypy.glycan obj
    3. extract motif
    4. save glycan_motif_addr
    """
    # root = r'/Users/apple/PycharmProjects/'
    # multiprocess
    manager = multiprocessing.Manager()
    x = load_json(glytoucan_data_base_addr__)
    glycan_list = []
    if prior:
        glycan_dict = load_json(for_extraction_glycoct_dict_addr)
        glycoct_dict={}
        glycan_list = list(glycan_dict.keys())
        for i in glycan_dict:
            glycoct_dict[i] = glycoct.loads(glycan_dict[i])

    else:

        for i in x.keys():
            if type(x[i]) is dict:
                if 'structure_' in x[i].keys():
                    glycan_list.append(i)
                else:
                    print('no structure', i)
            elif type(x[i]) is str:
                glycan_list.append(i)
        print('load ' + str(len(glycan_list)) + " glycan structure")
        store_json(for_extraction_glycoct_dict_addr, glycan_list)
        glycoct_dict = {}
        count_ = 0
        for i in glycan_list:
            count_ += 1
            print(i)
            if count_ % 1000 == 0:
                print(count_)
            if type(x[i]) is str:
                glycoct_dict[i] = glycoct.loads(x[i])
            elif type(x[i]) is dict and 'structure_' in x[i].keys():
                glycoct_dict[i] = glycoct.loads(x[i]['structure_'])
    glytoucan_motif_dic = manager.dict()
    success_list = manager.list()

    end_time = time.time()

    print('start parallel')

    pool = multiprocessing.Pool(processes=num_processors)
    for idex, i in enumerate(glycoct_dict):
        if len(glycoct_dict[i]) > gly_len: continue
        """ using get motif with count wrapper
            Also check exists wrapper
        """
        pool.apply_async(extract_motif_with_count_wrap, args=(idex, i, glycoct_dict[i], glytoucan_motif_dic, success_list))
        # print('finished ', idex)
    print("closing poll")
    pool.close()
    print('joining pool')
    pool.join()

    print('finished pool')
    a_motif_dic = dict(glytoucan_motif_dic)

    print('success_log')
    success_list = list(success_list)
    store_json(success_log_addr, success_list)

    str_motif = {}
    print('store duplicate')

    for i in glycan_list:
        if len(glycoct_dict[i]) > gly_len: continue
        str_motif[i] = {}
        if i not in a_motif_dic.keys(): continue
        for j in a_motif_dic[i].keys():
            str_motif[i][j] = [str(k) for k in a_motif_dic[i][j]]
    store_json(glycan_dict_motif_list_addr, str_motif)


if __name__ == '__main__':
    load_glycoct_for_database()
    get_motif_pip(gly_len=23, prior=True)
