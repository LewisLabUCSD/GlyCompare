from glypy.algorithms.subtree_search.inclusion import subtree_of
from glypy.io import glycoct
import glypy
import multiprocessing
# from . import __init__
from . import glycan_io
from . import json_utility
import pandas as pd
import numpy as np
import json



def _duplicate_cleaning_wrapper(degree, substructure_list, cleaned_substructure_dic, linkage_specific, reverse_dict):
    """
    :param degree: complexity degree
    :param substructure_list: a list of Glycan
    :param cleaned_substructure_dic: stored
    :return:
    """
    ldex = 0
    _check_list = substructure_list
    _ = []
    while ldex < len(_check_list):
        jdex = ldex + 1
        _.append(ldex)
        while jdex < len(_check_list):
            try:
                if not subtree_of(glycoct.loads(reverse_dict[_check_list[jdex]]), glycoct.loads(reverse_dict[_check_list[ldex]]), exact=linkage_specific) is None:
                    del _check_list[jdex]
                # elif not subtree_of(_check_list[ldex], _check_list[jdex]) is None:
                #     _check_list[ldex] = _check_list[ldex]
                #     del _check_list[jdex]
                else:
                    jdex += 1
            except AttributeError:
                print(ldex, jdex)
        ldex += 1
    cleaned_substructure_dic[degree] = _check_list


def merge_glycan_substructure_dict_to_substructure_dict(glycan_substructure_dict, glycan_dict, reference_dict, gly_len, combine_original=True):
    """
    merge glycan_substructure_dict to substructure_dict
    :param glycan_substructure_dict: dict_ID_dict_degree_list
    :param combine_original: Check if you need add the original glycan to the substructure-vec
    :return:
    """
    print("Start merge_glycan_substructure_to_substructure_dict")
    _substructure_dic = {}
    for i in sorted(list(glycan_substructure_dict.keys()), reverse=True):
        # print(count__)
        for j in glycan_substructure_dict[i].keys():
            if j not in _substructure_dic.keys():
                _substructure_dic[j] = glycan_substructure_dict[i][j][:]
            else:
                _substructure_dic[j].extend(glycan_substructure_dict[i][j][:])
    # print(_substructure_dic.keys())
    if combine_original:
        print("combine original")
        for i in glycan_dict.keys():
            if len(glycan_dict[i]) > gly_len:
                continue
            if str(len(glycan_dict[i])) not in _substructure_dic.keys():
                print('add new glycan degree', )
                _substructure_dic[str(len(glycan_dict[i]))] = [reference_dict[glycoct.dumps(glycan_dict[i])]]
            else:
                _substructure_dic[str(len(glycan_dict[i]))].append(reference_dict[glycoct.dumps(glycan_dict[i])])
    return _substructure_dic


def merge_substructure_dict_pip(glycan_substructure_dict, glycan_dict, linkage_specific, num_processors, reference_dict, reverse_dict, output_merged_substructure_glycoct_vec_addr, output_merged_substructure_glycoct_dict_addr=""):
    """
    merge the substructure of all glycans into substructure dict
    :param glycan_substructure_dict: {degree: [substructure1, substructure2, ... ]} /NBT_glycan_dict_degree_list_glycoct_for_substructure
    store the glycan substructure to NBT_substructure_dic_degree_list.json
    :param output_merged_substructure_glycoct_dict_addr: addr
    :return: sorted substructure_vec
    """
    # output_substructure_dic_degree_list_addr = root + "NBT_substructure_dic_degree_list.json"

#     glycan_io.check_glycan_substructure_dict(glycan_substructure_dict)
    gly_len = 25
    glycan_io.check_glycan_dict(glycan_dict)
    substructure_dict = merge_glycan_substructure_dict_to_substructure_dict(glycan_substructure_dict, glycan_dict=glycan_dict, reference_dict = reference_dict, gly_len = gly_len, combine_original=True)
#     print('substructure_dict is merged with len ', check_substructure_dict_length(_substructure_dic))

#     _substructure_degree = list(_substructure_dic.keys())
#     sorted_len_substructure_degree = sorted(_substructure_degree, key=lambda x: len(_substructure_dic[x]), reverse=True)
#     pool = multiprocessing.Pool(processes=num_processors)
#     manager = multiprocessing.Manager()
#     cleaned_substructure_dic = manager.dict()
#     for i in sorted_len_substructure_degree:
#         pool.apply_async(_duplicate_cleaning_wrapper, args=(i, _substructure_dic[i], cleaned_substructure_dic, linkage_specific, reverse_dict))

#     pool.close()
#     pool.join()

#     print('finished removing duplicate')
#     substructure_dict = dict(cleaned_substructure_dic)
#     substructure_dict_to_save = {}
#     for k in substructure_dict.keys():
#         temp = []
#         for g in substructure_dict[k]:
#             if type(glycoct.loads(reverse_dict[g])) == glypy.structure.glycan.Glycan and g not in temp:
#                 temp.append(g)

#         substructure_dict_to_save[k] = temp
    substructure_vec = []
    for i in sorted([int(j) for j in list(substructure_dict.keys())]):
        for g in substructure_dict[str(i)]:
            if type(glycoct.loads(reverse_dict[g])) == glypy.structure.glycan.Glycan:
                substructure_vec.append(g)
#                 substructure_vec.extend(substructure_dict_to_save[str(i)])
    substructure_vec_ = []
    for i in substructure_vec:
        if i not in substructure_vec_:
            substructure_vec_.append(i)
    substructure_vec = substructure_vec_

    print('after the cleaning the substructure vec\'s length is', check_substructure_dict_length(substructure_dict))

    if output_merged_substructure_glycoct_dict_addr != "":
#         json_utility.store_json(output_merged_substructure_glycoct_dict_addr, substructure_dict_to_save)
        with open(output_merged_substructure_glycoct_vec_addr, "w") as f:
            json.dump(substructure_vec, f)

    return substructure_vec



def match_substructure(substructure_vec, _glycan_substructure_dict, linkage_specific, reverse_dict):
    """
    :param substructure_vec: customized substructure vec
    :param glycan_substructure_dict: extracted_substructure_dic returned from the extract_substructure
    :return: match_vec
    """
#     glycan_substructure_dict = {}
#     for i in _glycan_substructure_dict.keys():
#         if type(_glycan_substructure_dict[i][0]) == str:
#             glycan_substructure_dict[i] = [glycoct.loads(k) for k in _glycan_substructure_dict[i]]
#         else:
#             glycan_substructure_dict[i] = _glycan_substructure_dict[i]

#     match_vec = [0] * len(substructure_vec)
    match_vec = []

    for jdex, j in enumerate(substructure_vec):
        _len = len(glycoct.loads(reverse_dict[j]))
        if _len - 1 >= len(_glycan_substructure_dict) or not _glycan_substructure_dict[_len - 1]: continue
        _iter = 0
#         times = _glycan_substructure_dict[_len - 1].count(j)
#         match_vec.extend([jdex for k in range(times)])
        while _iter < len(_glycan_substructure_dict[_len - 1]):
            if not subtree_of(glycoct.loads(reverse_dict[j]), glycoct.loads(reverse_dict[_glycan_substructure_dict[_len - 1][_iter]]), exact=linkage_specific) is None:
                match_vec.append(jdex)
#                 match_vec[jdex] += 1
#                 del _glycan_substructure_dict[str(_len)][_iter]
#             _iter += 1
                del _glycan_substructure_dict[_len - 1][_iter]
            else:
                _iter += 1

    return match_vec


def match_substructure_for_pip(substructure_vec, glycan_substructure_dict, glycan_id, linkage_specific, reverse_dict, idex=0):
    """
    :param substructure_vec: customized substructure vec
    :param glycan_substructure_dict: extracted_substructure_dic returned from the extract_substructure
    :param glycan_id: ID glytoucan or other
    :param match_dict: dict for storing final data
    :param idex: idex in list (mainly for progress check)
    :return:
    """
    # make duplicate
    try:
        match_vec = match_substructure(substructure_vec, glycan_substructure_dict, linkage_specific=linkage_specific, reverse_dict = reverse_dict)
    except:
        print('error', idex)
    print('finished ', idex)
    return glycan_id, match_vec


def substructure_matching_wrapper(substructure_vec, glycan_substructure_dict, linkage_specific, num_processors, reference_dict, reverse_dict, sub_glycoct_addr, matched_dict_addr=""):
    """
    match the glycan_substructure_dict_degree_list to substructure vec
    :param num_processors:
    :param substructure_: {degree: [substructure1, substructure2, ...]}  /NBT_substructure_dic_degree_list
    :param glycan_substructure_dict: {degree: [substructure1, substructure2, ...]} /NBT_glycan_dict_degree_list_glycoct_for_substructure
    :param matched_glycan_dict_addr: output_addr /NBT_fixed_gylcan_name_list
    :return: glycan_match_existed_substructure degree - >[substructure1, substructure2, ...]
    """
#     if type(substructure_) == dict:
#         substructure_vec = glycan_io.substructure_dict_to_substructure_vec(substructure_)
#     elif type(substructure_) == list:
#         substructure_vec = substructure_
#     else:
#         assert False, 'Error: incorrect type(substructure_)'
    print('get substructure vec, the length is ', len(substructure_vec))
    pool = multiprocessing.Pool(processes=num_processors)
    manager = multiprocessing.Manager()
#     match_dict = manager.dict()  # {substructure_name:[scores]}
    pool_list = []
    glycan_substructure_dict_simple = {}
    for key in glycan_substructure_dict:
        max_length = max([int(i) for i in list(glycan_substructure_dict[key].keys())])
        glycan_substructure_dict_simple[key] = [[] for i in range(max_length)]
        for length in glycan_substructure_dict[key]:
#             for sub_ in glycan_substructure_dict[key][length]:
#                 subtree = [reference_dict[glycoct.dumps(i)] for i in glycoct.loads(reverse_dict[sub_]).subtrees()]
#                 glycan_substructure_dict_simple[key][int(length) - 1].extend(subtree)
            glycan_substructure_dict_simple[key][int(length) - 1].extend(glycan_substructure_dict[key][length])
    for idex, i in enumerate(glycan_substructure_dict_simple):
        print('start processing', i)
        pool_list.append(pool.apply_async(match_substructure_for_pip, args=(substructure_vec, glycan_substructure_dict_simple[i], i, linkage_specific, reverse_dict, idex)))
    matrix = [[0 if j not in i.get()[1] else i.get()[1].count(j) for j in range(len(substructure_vec))] for i in pool_list]
    matched_df = pd.DataFrame(np.array(matrix).transpose(), columns = [i.get()[0] for i in pool_list])
    print("closing poll")
    pool.close()
    print('joining pool')
    pool.join()
    print('converting dict')
    nonzero_rows = (matched_df!=0).any(axis=1)
    substructure_vec_nonzero = [substructure_vec[i] for i in range(len(substructure_vec)) if nonzero_rows[i]]
    with open(sub_glycoct_addr, "w") as f:
        json.dump(substructure_vec_nonzero, f)
        
    matched_df_nonzero = matched_df[(matched_df!=0).any(axis=1)]
    matched_df_nonzero = matched_df_nonzero.reset_index(drop=True)
    
    if matched_dict_addr != "":
        matched_df_nonzero.to_csv(matched_dict_addr)
    return matched_df_nonzero


def check_substructure_dict_length(a_dict):
    return sum([len(a_dict[i]) for i in a_dict.keys()])




if __name__ == '__main__':
    # customizing_substructure_vec_pip()
    pass
