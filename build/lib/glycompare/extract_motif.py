import glypy
from glypy.algorithms.subtree_search.inclusion import subtree_of
from glypy.io import glycoct
from glypy.structure.glycan import fragment_to_substructure
from glypy import Glycan
import time
import multiprocessing

import glycan_io
import __init__
from json_utility import load_json, store_json


# def extract_motif(glycoct_obj, idex=0):
#     # print('start getmotif')
#     _frag_motif_list = {}
#     # fig, axes = plt.subplots(6,9)
#     # fig.set_size_inches(14,6)
#     start_time = time.time()
#     for i in glycoct_obj.fragments(max_cleavages=5):
#         _frag_gly = fragment_to_substructure(i, glycoct_obj)
#         # plot(_frag_gly)
#         if not len(_frag_gly) in _frag_motif_list.keys():
#             _frag_motif_list[len(_frag_gly)] = [glycoct.loads(str(_frag_gly))]
#         else:
#             _frag_motif_list[len(_frag_gly)].append(glycoct.loads(str(_frag_gly)))
#     mid_time = time.time()
#     # print('start clean duplicate')
#     _frag_motif_list = clean_duplicate(_frag_motif_list)
#     end_time = time.time()
#     print(idex, len(glycoct_obj), end_time - mid_time, mid_time - start_time)
#     # print('finished getmotif')
#
#     return _frag_motif_list


def clean_duplicate(a_motif_dict, linkage_specific):
    """

    :param a_motif_dict: remove the duplicates in a motif_dict_degree_list
    :return:
    """
    for i in a_motif_dict.keys():
        # print(i)
        ldex = 0
        _check_list = a_motif_dict[i]
        while ldex < len(_check_list):
            jdex = ldex + 1
            while jdex < len(_check_list):
                if subtree_of(_check_list[ldex], _check_list[jdex], exact=linkage_specific) == 1:
                    del _check_list[jdex]
                else:
                    jdex += 1
            ldex += 1
        a_motif_dict[i] = _check_list
    return a_motif_dict



def extract_motif(a_glycan, branch=5):
    """
    :param a_glycan: Glycan obj
    :param branch:
    :return:
    """
    # print('start getmotif')
    extracted_motif_dic = {}
    # print('aa1')
    # print(a_glycan)
    # a_glycan=glycoct.loads(a_glycan)
    # print(isinstance(a_glycan, Glycan))
    # print('ab')
    for i in a_glycan.fragments(max_cleavages=branch):
        # print('aaa')
        _frag_gly = fragment_to_substructure(i, a_glycan)
        # print('aab')

        if not str(len(_frag_gly)) in extracted_motif_dic.keys():
            extracted_motif_dic[str(len(_frag_gly))] = [_frag_gly]
        else:
            extracted_motif_dic[str(len(_frag_gly))].append(_frag_gly)
    # print('ab')
    extracted_motif_dic[str(len(a_glycan))] = [a_glycan]
    # print('ac')
    return extracted_motif_dic


def extract_motif_wrapper(a_name, a_glycan_str, motif_dic):
    """
    :param idex: idex of Glycan in the list
    :param a_name: the ID of the Glycan
    :param a_glycan_str: Glycan obj
    :param motif_dic: dict for storing the extracted motif dict
    :return:
    """

    try:
        print('start', a_name)
        start_time = time.time()
        motif_dic[a_name] = extract_motif(a_glycan_str)
        end_time = time.time()
        print(a_name, len(motif_dic[a_name]), end_time - start_time)
        # print('has_motif', motif_dic[_name])
        #     for j in motif_dic.keys():
        #         print(j, len(motif_dic[j]))
    except TypeError:
        print(a_name, 'has error')
    except KeyboardInterrupt:
        print('break')


def get_motif_pip(glycan_dict, gly_len, output_file, num_processors=__init__.num_processors):
    """Please set the prior=True to get the data file please run the NBT_GLYCAN_preprocess file
    If prior=False, it will generate glycan motif for all glycan in glytoucan database
    1. load  {glyacn_id: glycan_str}
    2. convert to glypy.glycan obj {glyacn_id: glycan_}
    3. extract {glyacn_id: motif_dict}
    4. save glycan_motif_dic
    :param glycan_dict:
    :param gly_len: the max degree of the glycan that can be processed
    :param output_file: store str type of glycan_motif_dict
    """
    # root = r'/Users/apple/PycharmProjects/'
    # multiprocess
    glycan_io.check_glycan_dict(glycan_dict)
    manager = multiprocessing.Manager()
    motif_dic = manager.dict()
    print('start parallel parsing', len(glycan_dict), 'glycans')
    pool = multiprocessing.Pool(processes=num_processors)
    pool_list = []
    for idex, i in enumerate(glycan_dict):
        if len(glycan_dict[i]) > gly_len:
            print(i, 'larger than max')
            continue
        """ using get motif with count wrapper
            Also check exists wrapper
        """
        pool_list.append(pool.apply_async(extract_motif_wrapper, args=(i, glycan_dict[i], motif_dic)))
    result_list = [xx.get() for xx in pool_list]
    # print('finished ', idex)
    # print("closing poll")
    pool.close()
    # print('joining pool')
    pool.join()
    print('finished pool')
    print('glycan_dict', len(motif_dic))
    glycan_motif_dic = dict(motif_dic)

    str_motif = {}
    for i in glycan_dict:
        # if len(glycan_dict[i]) > gly_len: continue
        if i not in glycan_motif_dic.keys():
            print('missing', i)
            continue
        str_motif[i] = {}
        for j in glycan_motif_dic[i]:
            str_motif[i][j] = [str(k) for k in glycan_motif_dic[i][j]]
    if output_file != '':
        store_json(output_file, str_motif)
    return glycan_motif_dic


def main():
    # glycan_dict = glycan_io.load_glycan_obj_from_database(topology_list_addr=__init__.topology_list_addr, output_file=__init__.glycan_dict_addr, loader=glycoct)
    # glycan_motif_dic = get_motif_pip(glycan_dict=glycan_dict, gly_len=23, output_file=__init__.glycan_motif_dict_addr)
    pass


if __name__ == '__main__':
    main()
