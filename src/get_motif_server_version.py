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
import NBT_init

def store_json(address, dic):
    with open(address, 'w') as fp:
        json.dump(dic, fp)


def load_json(address):
    with open(address, 'r') as f:
        return json.load(f)


def get_motif(glycoct_obj, idex=0):
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


def get_motif_with_count(glycoct_obj, idex=0):
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


def get_motif_wrap(idex, _name, glycoct_obj, motif_dic, success_list):
    if not _name in success_list:
        try:
            # print('start', _name, glycoct_obj)

            _temp_motif = get_motif(glycoct_obj, idex)
            motif_dic[_name] = _temp_motif
            # print('has_motif', motif_dic[_name])

            success_list.append(_name)
            #     for j in motif_dic.keys():
            #         print(j, len(motif_dic[j]))
        except TypeError:
            print(idex, _name, 'has error')
        except KeyboardInterrupt:
            print('break')


def get_motif_with_count_wrap(idex, _name, glycoct_obj, motif_dic, success_list):
    if not _name in success_list:
        try:
            # print('start', _name, glycoct_obj)

            _temp_motif = get_motif_with_count(glycoct_obj, idex)
            motif_dic[_name] = _temp_motif
            # print('has_motif', motif_dic[_name])

            success_list.append(_name)
            #     for j in motif_dic.keys():
            #         print(j, len(motif_dic[j]))
        except TypeError:
            print(idex, _name, 'has error')
        except KeyboardInterrupt:
            print('break')


def get_motif_pip(gly_len, used=False):
    """to get the init file please run the NBT_GLYCAN_preprocess file"""
    # root = r'/Users/apple/PycharmProjects/'
    root = NBT_init.json_address
    # root = "./"
    glytoucan_data_base_addr__ = root + 'BNT_glycan_iupac_dic.json'
    success_log_addr = root + 'BNT_for_motif_log.json'
    glycan_list_addr = root + 'BNT_glycan_list_127.json'
    glycan_motif_addr = root + 'BNT_glycan_with_motif.json'
    x = load_json(glytoucan_data_base_addr__)

    # multiprocess
    manager = multiprocessing.Manager()

    glycan_list = []
    if used:
        glycan_list = load_json(glycan_list_addr)
        # glytoucan_motif_dic = manager.dict(load_json(glycan_motif_addr))
        # success_list = manager.list(load_json(success_log_addr))

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
        store_json(glycan_list_addr, glycan_list)

    glytoucan_motif_dic = manager.dict()
    success_list = manager.list()
    # print(glycan_list[:2])

    end_time = time.time()

    # success_list = manager.list()
    glycoct_dict = {}
    count_ = 0
    for i in glycan_list:
        count_ += 1
        if count_ % 1000 == 0:
            print(count_)
        if type(x[i]) is str:
            glycoct_dict[i] = glycoct.loads(x[i])
        elif type(x[i]) is dict and 'structure_' in x[i].keys():
            glycoct_dict[i] = glycoct.loads(x[i]['structure_'])

    print('start parappel')
    num_processors = 8
    pool = multiprocessing.Pool(processes=num_processors)
    for idex, i in enumerate(glycan_list):
        if len(glycoct_dict[i]) > gly_len: continue
        """ using get motif with count wrapper
            Also check exists wrapper
        """
        pool.apply_async(get_motif_with_count_wrap, args=(idex, i, glycoct_dict[i], glytoucan_motif_dic, success_list))
        # print('finished ', idex)
    print("closing poll")
    pool.close()
    print('joining pool')
    pool.join()

    print('finished pool')
    a_motif_dic = dict(glytoucan_motif_dic)
    # print(a_motif_dic)

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
    store_json(glycan_motif_addr, str_motif)


if __name__ == '__main__':
    get_motif_pip(gly_len=23, used=False)
