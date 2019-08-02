from glypy.algorithms.subtree_search.inclusion import subtree_of
from glypy.structure.glycan import fragment_to_substructure
import time
import multiprocessing

from . import glycan_io
from . import __init__
from . import json_utility


def clean_duplicate(a_substructure_dict, linkage_specific):
    """

    :param a_substructure_dict: remove the duplicates in a substructure_dict_degree_list
    :return:
    """
    for i in a_substructure_dict.keys():
        # print(i)
        ldex = 0
        _check_list = a_substructure_dict[i]
        while ldex < len(_check_list):
            jdex = ldex + 1
            while jdex < len(_check_list):
                if subtree_of(_check_list[ldex], _check_list[jdex], exact=linkage_specific) == 1:
                    del _check_list[jdex]
                else:
                    jdex += 1
            ldex += 1
        a_substructure_dict[i] = _check_list
    return a_substructure_dict



def extract_substructure(a_glycan, branch=5):
    """
    :param a_glycan: Glycan obj
    :param branch:
    :return:
    """
    extracted_substructure_dic = {}

    for i in a_glycan.fragments(max_cleavages=branch):
        # print('aaa')
        _frag_gly = fragment_to_substructure(i, a_glycan)
        # print('aab')

        if not str(len(_frag_gly)) in extracted_substructure_dic.keys():
            extracted_substructure_dic[str(len(_frag_gly))] = [_frag_gly]
        else:
            extracted_substructure_dic[str(len(_frag_gly))].append(_frag_gly)
    # print('ab')
    extracted_substructure_dic[str(len(a_glycan))] = [a_glycan]
    # print('ac')
    return extracted_substructure_dic


def extract_substructure_wrapper(a_name, a_glycan_str, substructure_dic):
    """
    :param idex: idex of Glycan in the list
    :param a_name: the ID of the Glycan
    :param a_glycan_str: Glycan obj
    :param substructure_dic: dict for storing the extracted substructure dict
    :return:
    """

    try:
        print('start', a_name)
        start_time = time.time()
        substructure_dic[a_name] = extract_substructure(a_glycan_str)
        end_time = time.time()
        print(a_name, len(substructure_dic[a_name]), end_time - start_time)
        # print('has_substructure', substructure_dic[_name])
        #     for j in substructure_dic.keys():
        #         print(j, len(substructure_dic[j]))
    except TypeError:
        print(a_name, 'has error')
    except KeyboardInterrupt:
        print('break')


def extract_substructures_pip(glycan_dict, gly_len, output_file, num_processors):
    """Please set the prior=True to get the data file please run the NBT_GLYCAN_preprocess file
    If prior=False, it will generate glycan substructure for all glycan in glytoucan database
    1. load  {glyacn_id: glycan_str}
    2. convert to glypy.glycan obj {glyacn_id: glycan_}
    3. extract {glyacn_id: substructure_dict}
    4. save glycan_substructure_dic
    :param glycan_dict:
    :param gly_len: the max degree of the glycan that can be processed
    :param output_file: store str type of glycan_substructure_dict
    """
    # root = r'/Users/apple/PycharmProjects/'
    # multiprocess
    glycan_io.check_glycan_dict(glycan_dict)
    manager = multiprocessing.Manager()
    substructure_dic = manager.dict()
    print('start parallel parsing', len(glycan_dict), 'glycans')
    pool = multiprocessing.Pool(processes=num_processors)
    pool_list = []
    for idex, i in enumerate(glycan_dict):
        if len(glycan_dict[i]) > gly_len:
            print(i, 'larger than max')
            continue
        """ using get substructure with count wrapper
            Also check exists wrapper
        """
        pool_list.append(pool.apply_async(extract_substructure_wrapper, args=(i, glycan_dict[i], substructure_dic)))
    result_list = [xx.get() for xx in pool_list]
    # print('finished ', idex)
    # print("closing poll")
    pool.close()
    # print('joining pool')
    pool.join()
    print('finished pool')
    print('glycan_dict', len(substructure_dic))
    glycan_substructure_dic = dict(substructure_dic)

    str_substructure = {}
    for i in glycan_dict:
        # if len(glycan_dict[i]) > gly_len: continue
        if i not in glycan_substructure_dic.keys():
            print('missing', i)
            continue
        str_substructure[i] = {}
        for j in glycan_substructure_dic[i]:
            str_substructure[i][j] = [str(k) for k in glycan_substructure_dic[i][j]]
    if output_file != '':
        json_utility.store_json(output_file, str_substructure)
    return glycan_substructure_dic


# def main():
#     # glycan_dict = glycan_io.load_glycan_obj_from_database(topology_list_addr=__init__.topology_list_addr, output_file=__init__.glycan_dict_addr, loader=glycoct)
#     # glycan_substructure_dic = get_substructure_pip(glycan_dict=glycan_dict, gly_len=23, output_file=__init__.glycan_substructure_dict_addr)
#     pass
#
#
# if __name__ == '__main__':
#     main()
