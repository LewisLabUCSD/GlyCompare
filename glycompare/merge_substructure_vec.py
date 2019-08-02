from glypy.algorithms.subtree_search.inclusion import subtree_of
from glypy.io import glycoct
import multiprocessing
# from . import __init__
from . import glycan_io
from . import json_utility


# def profile_loader(profile_dict, name):
#     return glycan_profile_obj(profile_dict['glycan'],
#                               profile_dict['m/z'],
#                               profile_dict['weight'],
#                               profile_dict['substructure_vec'],
#                               name)



# def merge_unzero_vec(prof_n, dict_name_abundance_cross_profile, glycan_dict, glycan_profile):
#     """prof_n: start with 0"""
#
#     if prof_n + 1 < 10:
#         _id = 'Gly0' + str(prof_n + 1)
#     else:
#         _id = 'Gly' + str(prof_n + 1)
#     # print(_id)
#     for i in dict_name_abundance_cross_profile.keys():
#         if dict_name_abundance_cross_profile[i][prof_n] > 0.01:
#             #             print(glycan_profile[_id].keys())
#             if i in glycan_profile[_id].keys(): continue
#             if len(glycan_dict[int(i)].keys()) > 1:
#                 print(dict_name_abundance_cross_profile[i][prof_n], '\'' + i + '\'', list(glycan_dict[int(i)].keys()))
#                 continue
#             glycan_profile[_id][i] = list(glycan_dict[int(i)].keys())[0]
#     return glycan_profile[_id]


# merge_unzero_vec(0, dict_name_abundance_cross_profile, glycan_dict, glycan_profile)

# def merge_profile():
#     glycan_profile = load_json(root_address + "NBT_glycan_profile.json")
#     NBT_glycan_match_existed_substructure = load_json(root_address + "NBT_glycan_match_existed_substructure.json")
#     dict_name_abundance_cross_profile = store_json(root_address + r"NBT_dict_name_abundance_cross_profile.json")
#
#     profile_obj_list = []
#     for i in range(1, 38):
#         if i < 10:
#             _id = 'Gly0' + str(i)
#         else:
#             _id = 'Gly' + str(i)
#             #     print(_id)
#         weighted_matrix = np.zeros((len(NBT_glycan_match_existed_substructure["3865.1"])))
#         #     print(weighted_matrix.shape)
#         abundance_ = []
#         # glycan_hit_array_ = []
#         mz_ = []
#         glycan_id_ = []
#         for i in sorted(list(glycan_profile[_id].keys())):
#             _name = glycan_profile[_id][i]
#             mz_.append(i)
#             glycan_id_.append(_name)
#
#             _bundance = dict_name_abundance_cross_profile[int(i)][0]
#             _temp_hit_matrix = np.array(NBT_glycan_match_existed_substructure[_name])
#
#             # glycan_hit_array.append(_temp_hit_matrix)
#             abundance_.append(_bundance)
#             weighted_matrix += _temp_hit_matrix * _bundance
#         profile_obj_list.append(glycan_profile_obj(glycan_id_, mz_, abundance_, weighted_matrix))


def _duplicate_cleaning_wrapper(degree, substructure_list, cleaned_substructure_dic, linkage_specific):
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
                # print(isinstance(_check_list[jdex], glypy.Glycan))
                # print(isinstance(_check_list[ldex], glypy.Glycan))

                if not subtree_of(_check_list[jdex], _check_list[ldex], exact=linkage_specific) is None:
                    del _check_list[jdex]
                # elif not subtree_of(_check_list[ldex], _check_list[jdex]) is None:
                #     _check_list[ldex] = _check_list[ldex]
                #     del _check_list[jdex]
                else:
                    jdex += 1
                # exit()
            except AttributeError:
                print(ldex, jdex)
                # exit()
        # if not find_same:
        ldex += 1
    cleaned_substructure_dic[degree] = _check_list


def merge_glycan_substructure_dict_to_substructure_dict(glycan_substructure_dict, glycan_dict, combine_original=True):
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
            if str(len(glycan_dict[i])) not in _substructure_dic.keys():
                print('add new glycan degree', )
                _substructure_dic[str(len(glycan_dict[i]))] = [glycan_dict[i]]
            else:
                _substructure_dic[str(len(glycan_dict[i]))].append(glycan_dict[i])
    return _substructure_dic


def merge_substructure_dict_pip(glycan_substructure_dict, glycan_dict, linkage_specific, num_processors, output_merged_substructure_glycoct_dict_addr=""):
    """
    merge the substructure of all glycans into substructure dict
    :param glycan_substructure_dict: {degree: [substructure1, substructure2, ... ]} /NBT_glycan_dict_degree_list_glycoct_for_substructure
    store the glycan substructure to NBT_substructure_dic_degree_list.json
    :param output_merged_substructure_glycoct_dict_addr: addr
    :return: sorted substructure_vec
    """
    # output_substructure_dic_degree_list_addr = root + "NBT_substructure_dic_degree_list.json"

    glycan_io.check_glycan_substructure_dict(glycan_substructure_dict)
    # glycan_dict = glycan_io.load_glycan_dict_from_json(glycan_dict_addr)
    glycan_io.check_glycan_dict(glycan_dict)
    _substructure_dic = merge_glycan_substructure_dict_to_substructure_dict(glycan_substructure_dict, glycan_dict=glycan_dict, combine_original=True)
    glycan_io.check_substructure_dict(_substructure_dic)
    # exit()
    print('substructure_dict is merged with len ', check_substructure_dict_length(_substructure_dic))
    # print("get_substructure_dict_degree_list_pipe")
    # num_processors = 8

    _substructure_degree = list(_substructure_dic.keys())
    sorted_len_substructure_degree = sorted(_substructure_degree, key=lambda x: len(_substructure_dic[x]), reverse=True)
    pool = multiprocessing.Pool(processes=num_processors)
    manager = multiprocessing.Manager()
    cleaned_substructure_dic = manager.dict()
    for i in sorted_len_substructure_degree:
        pool.apply_async(_duplicate_cleaning_wrapper, args=(i, _substructure_dic[i], cleaned_substructure_dic, linkage_specific,))

    # print("closing poll")
    pool.close()
    # print('joining pool')
    pool.join()

    print('finished removing duplicate')
    substructure_dict = dict(cleaned_substructure_dic)

    print('after the cleaning the substructure vec\'s length is', check_substructure_dict_length(substructure_dict))

    substructure_glycoct_dict_str = {}
    # print(substructure_dict.keys())
    # for i in sorted_len_substructure_degree:
    #     print(i, len(substructure_dict[i]))
    #     substructure_glycoct_dict_str[i] = [str(j) for j in substructure_dict


    substructure_glycoct_dict_str = glycan_io.glycan_obj_to_glycan_str(glycan_io.substructure_dict_to_substructure_vec(substructure_dict))


    if output_merged_substructure_glycoct_dict_addr != "":
        json_utility.store_json(output_merged_substructure_glycoct_dict_addr, substructure_glycoct_dict_str)

    return substructure_dict


def match_substructure(substructure_vec, _glycan_substructure_dict, linkage_specific):
    """
    :param substructure_vec: customized substructure vec
    :param glycan_substructure_dict: extracted_substructure_dic returned from the extract_substructure
    :return: match_vec
    """
    glycan_substructure_dict = {}
    for i in _glycan_substructure_dict.keys():
        if type(_glycan_substructure_dict[i][0]) == str:
            glycan_substructure_dict[i] = [glycoct.loads(k) for k in _glycan_substructure_dict[i]]
        else:
            glycan_substructure_dict[i] = _glycan_substructure_dict[i]
    match_vec = [0] * len(substructure_vec)

    for jdex, j in enumerate(substructure_vec):
        #             print(type(j), type(_tree))
        _len = len(j)
        if str(_len) not in glycan_substructure_dict.keys(): continue
        _iter = 0
        # print(jdex)
        while _iter < len(glycan_substructure_dict[str(_len)]):
            # print(type(j), type(i))

            if not subtree_of(j, glycan_substructure_dict[str(_len)][_iter], exact=linkage_specific) is None:
                # print('yes')
                match_vec[jdex] += 1
                del glycan_substructure_dict[str(_len)][_iter]
            else:
                _iter += 1

    return match_vec


def match_substructure_for_pip(substructure_vec, glycan_substructure_dict, glycan_id, match_dict, linkage_specific, idex=0):
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
        match_vec = match_substructure(substructure_vec, glycan_substructure_dict, linkage_specific=linkage_specific)
        match_dict[glycan_id] = match_vec[:]
    except:
        print('error', idex)
    print('finished ', idex)


def substructure_matching_wrapper(substructure_, glycan_substructure_dict, linkage_specific, num_processors, matched_dict_addr=""):
    """
    match the glycan_substructure_dict_degree_list to substructure vec
    :param num_processors:
    :param substructure_: {degree: [substructure1, substructure2, ...]}  /NBT_substructure_dic_degree_list
    :param glycan_substructure_dict: {degree: [substructure1, substructure2, ...]} /NBT_glycan_dict_degree_list_glycoct_for_substructure
    :param matched_glycan_dict_addr: output_addr /NBT_fixed_gylcan_name_list
    :return: glycan_match_existed_substructure degree - >[substructure1, substructure2, ...]
    """
    glycan_io.check_glycan_substructure_dict(glycan_substructure_dict)
    glycan_io.check_substructure_dict(substructure_)

    if type(substructure_) == dict:
        substructure_vec = glycan_io.substructure_dict_to_substructure_vec(substructure_)
    elif type(substructure_) == list:
        substructure_vec = substructure_
    else:
        assert False, 'Error: incorrect type(substructure_)'
    # store_json(output_substructure_vec_addr, [str(i) for i in substructure_vec])
    print('get substructure vec, the length is ', len(substructure_vec))
    pool = multiprocessing.Pool(processes=num_processors)
    manager = multiprocessing.Manager()
    match_dict = manager.dict()  # {substructure_name:[scores]}
    for idex, i in enumerate(glycan_substructure_dict):
        print('start processing', i)
        pool.apply_async(match_substructure_for_pip, args=(substructure_vec, glycan_substructure_dict[i], i, match_dict, linkage_specific, idex))
    print("closing poll")
    pool.close()
    print('joining pool')
    pool.join()
    print('converting dict')
    return_matched_dic = dict(match_dict)
    # print(return_matched_dic)
    if matched_dict_addr != "":
        json_utility.store_json(matched_dict_addr, return_matched_dic)

    return return_matched_dic


def check_substructure_dict_length(a_dict):
    return sum([len(a_dict[i]) for i in a_dict.keys()])


# def customizing_substructure_vec_pip():
#     """ merge the substructure of all glycans into substructure dict"""
#     print('start merging')
#     glycan_substructure_dict = glycan_io.glycan_str_to_glycan_obj(
#         load_json(__init__.glycan_substructure_dict_addr))
#
#     merged_substructure_dict = merge_substructure_dict_pipe(glycan_substructure_dict,
#                                               output_merged_substructure_dict_addr=__init__.merged_substructure_dict_addr)
#
#     """ Start substructure matching"""
#     print('start substructure match')
#     merged_substructure_dict = glycan_io.glycan_str_to_glycan_obj(load_json(__init__.merged_substructure_dict_addr))
#     match_dict = substructure_matching_wrapper(merged_substructure_dict,
#                                         glycan_substructure_dict,
#                                         matched_glycan_dict_addr=__init__.output_matched_dict_addr)

    # nbt_table = pd.read_table(input_profile_addr)
    # a_profile = build_profiles(nbt_glycan_dict, nbt_table)


if __name__ == '__main__':
    # customizing_substructure_vec_pip()
    pass
