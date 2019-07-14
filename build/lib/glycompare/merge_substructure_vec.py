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
#                               profile_dict['motif_vec'],
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
#     NBT_glycan_match_existed_motif = load_json(root_address + "NBT_glycan_match_existed_motif.json")
#     dict_name_abundance_cross_profile = store_json(root_address + r"NBT_dict_name_abundance_cross_profile.json")
#
#     profile_obj_list = []
#     for i in range(1, 38):
#         if i < 10:
#             _id = 'Gly0' + str(i)
#         else:
#             _id = 'Gly' + str(i)
#             #     print(_id)
#         weighted_matrix = np.zeros((len(NBT_glycan_match_existed_motif["3865.1"])))
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
#             _temp_hit_matrix = np.array(NBT_glycan_match_existed_motif[_name])
#
#             # glycan_hit_array.append(_temp_hit_matrix)
#             abundance_.append(_bundance)
#             weighted_matrix += _temp_hit_matrix * _bundance
#         profile_obj_list.append(glycan_profile_obj(glycan_id_, mz_, abundance_, weighted_matrix))


def _duplicate_cleaning_wrapper(degree, motif_list, cleaned_motif_dic, linkage_specific):
    """
    :param degree: complexity degree
    :param motif_list: a list of Glycan
    :param cleaned_motif_dic: stored
    :return:
    """
    ldex = 0
    _check_list = motif_list
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
    cleaned_motif_dic[degree] = _check_list


def merge_glycan_motif_dict_to_motif_dict(glycan_motif_dict, glycan_dict, combine_original=True):
    """
    merge glycan_motif_dict to motif_dict
    :param glycan_motif_dict: dict_ID_dict_degree_list
    :param combine_original: Check if you need add the original glycan to the motif-vec
    :return:
    """
    print("Start merge_glycan_motif_to_motif_dict")
    _motif_dic = {}
    for i in sorted(list(glycan_motif_dict.keys()), reverse=True):
        # print(count__)
        for j in glycan_motif_dict[i].keys():
            if j not in _motif_dic.keys():
                _motif_dic[j] = glycan_motif_dict[i][j][:]
            else:
                _motif_dic[j].extend(glycan_motif_dict[i][j][:])
    # print(_motif_dic.keys())
    if combine_original:
        print("combine original")
        for i in glycan_dict.keys():
            if str(len(glycan_dict[i])) not in _motif_dic.keys():
                print('add new glycan degree', )
                _motif_dic[str(len(glycan_dict[i]))] = [glycan_dict[i]]
            else:
                _motif_dic[str(len(glycan_dict[i]))].append(glycan_dict[i])
    return _motif_dic


def merge_substructure_dict_pip(glycan_motif_dict, glycan_dict, linkage_specific, num_processors, output_merged_motif_dict_addr=""):
    """
    merge the substructure of all glycans into motif dict
    :param glycan_motif_dict: {degree: [motif1, motif2, ... ]} /NBT_glycan_dict_degree_list_glycoct_for_motif
    store the glycan motif to NBT_motif_dic_degree_list.json
    :param output_merged_motif_dict_addr: addr
    :return: sorted motif_vec
    """
    # output_motif_dic_degree_list_addr = root + "NBT_motif_dic_degree_list.json"

    glycan_io.check_glycan_motif_dict(glycan_motif_dict)
    # glycan_dict = glycan_io.load_glycan_dict_from_json(glycan_dict_addr)
    glycan_io.check_glycan_dict(glycan_dict)
    _motif_dic = merge_glycan_motif_dict_to_motif_dict(glycan_motif_dict, glycan_dict=glycan_dict, combine_original=True)
    glycan_io.check_motif_dict(_motif_dic)
    # exit()
    print('motif_dict is merged with len ', check_motif_dict_length(_motif_dic))
    # print("get_motif_dict_degree_list_pipe")
    # num_processors = 8

    _motif_degree = list(_motif_dic.keys())
    sorted_len_motif_degree = sorted(_motif_degree, key=lambda x: len(_motif_dic[x]), reverse=True)
    pool = multiprocessing.Pool(processes=num_processors)
    manager = multiprocessing.Manager()
    cleaned_motif_dic = manager.dict()
    for i in sorted_len_motif_degree:
        pool.apply_async(_duplicate_cleaning_wrapper, args=(i, _motif_dic[i], cleaned_motif_dic, linkage_specific,))

    # print("closing poll")
    pool.close()
    # print('joining pool')
    pool.join()

    print('finished removing duplicate')
    motif_dict = dict(cleaned_motif_dic)

    print('after the cleaning the motif vec\'s length is', check_motif_dict_length(motif_dict))

    motif_dict_str = {}
    # print(motif_dict.keys())
    # for i in sorted_len_motif_degree:
    #     print(i, len(motif_dict[i]))
    #     motif_dict_str[i] = [str(j) for j in motif_dict[i]]

    motif_dict_str = glycan_io.glycan_obj_to_glycan_str(motif_dict)
    if output_merged_motif_dict_addr != "":
        json_utility.store_json(output_merged_motif_dict_addr, motif_dict_str)
    return motif_dict


def match_motif(motif_vec, _glycan_motif_dict, linkage_specific):
    """
    :param motif_vec: customized motif vec
    :param glycan_motif_dict: extracted_motif_dic returned from the extract_motif
    :return: match_vec
    """
    glycan_motif_dict = {}
    for i in _glycan_motif_dict.keys():
        if type(_glycan_motif_dict[i][0]) == str:
            glycan_motif_dict[i] = [glycoct.loads(k) for k in _glycan_motif_dict[i]]
        else:
            glycan_motif_dict[i] = _glycan_motif_dict[i]
    match_vec = [0] * len(motif_vec)

    for jdex, j in enumerate(motif_vec):
        #             print(type(j), type(_tree))
        _len = len(j)
        if str(_len) not in glycan_motif_dict.keys(): continue
        _iter = 0
        # print(jdex)
        while _iter < len(glycan_motif_dict[str(_len)]):
            # print(type(j), type(i))

            if not subtree_of(j, glycan_motif_dict[str(_len)][_iter], exact=linkage_specific) is None:
                # print('yes')
                match_vec[jdex] += 1
                del glycan_motif_dict[str(_len)][_iter]
            else:
                _iter += 1

    return match_vec


def match_motif_for_pip(motif_vec, glycan_motif_dict, glycan_id, match_dict, linkage_specific, idex=0):
    """
    :param motif_vec: customized motif vec
    :param glycan_motif_dict: extracted_motif_dic returned from the extract_motif
    :param glycan_id: ID glytoucan or other
    :param match_dict: dict for storing final data
    :param idex: idex in list (mainly for progress check)
    :return:
    """
    # make duplicate
    try:
        match_vec = match_motif(motif_vec, glycan_motif_dict, linkage_specific=linkage_specific)
        match_dict[glycan_id] = match_vec[:]
    except:
        print('error', idex)
    print('finished ', idex)


def substructure_matching_wrapper(motif_dict, glycan_motif_dict, linkage_specific, num_processors, matched_glycan_dict_addr=""):
    """
    match the glycan_motif_dict_degree_list to motif vec
    :param num_processors:
    :param motif_dict: {degree: [motif1, motif2, ...]}  /NBT_motif_dic_degree_list
    :param glycan_motif_dict: {degree: [motif1, motif2, ...]} /NBT_glycan_dict_degree_list_glycoct_for_motif
    :param matched_glycan_dict_addr: output_addr /NBT_fixed_gylcan_name_list
    :return: glycan_match_existed_motif degree - >[motif1, motif2, ...]
    """
    glycan_io.check_glycan_motif_dict(glycan_motif_dict)
    glycan_io.check_motif_dict(motif_dict)

    if type(motif_dict) == dict:
        motif_vec = glycan_io.motif_dict_to_motif_vec(motif_dict)
    else:
        motif_vec = motif_dict
    # store_json(output_motif_vec_addr, [str(i) for i in motif_vec])
    print('get motif vec, the length is ', len(motif_vec))
    pool = multiprocessing.Pool(processes=num_processors)
    manager = multiprocessing.Manager()
    match_dict = manager.dict()  # {motif_name:[scores]}
    for idex, i in enumerate(glycan_motif_dict):
        print('start processing', i)
        pool.apply_async(match_motif_for_pip, args=(motif_vec, glycan_motif_dict[i], i, match_dict, linkage_specific, idex))
    print("closing poll")
    pool.close()
    print('joining pool')
    pool.join()
    print('converting dict')
    return_matched_dic = dict(match_dict)
    # print(return_matched_dic)
    if matched_glycan_dict_addr != "":
        json_utility.store_json(matched_glycan_dict_addr, return_matched_dic)

    return return_matched_dic


def check_motif_dict_length(a_dict):
    return sum([len(a_dict[i]) for i in a_dict.keys()])


# def customizing_motif_vec_pip():
#     """ merge the substructure of all glycans into motif dict"""
#     print('start merging')
#     glycan_motif_dict = glycan_io.glycan_str_to_glycan_obj(
#         load_json(__init__.glycan_motif_dict_addr))
#
#     merged_motif_dict = merge_motif_dict_pipe(glycan_motif_dict,
#                                               output_merged_motif_dict_addr=__init__.merged_motif_dict_addr)
#
#     """ Start motif matching"""
#     print('start motif match')
#     merged_motif_dict = glycan_io.glycan_str_to_glycan_obj(load_json(__init__.merged_motif_dict_addr))
#     match_dict = motif_matching_wrapper(merged_motif_dict,
#                                         glycan_motif_dict,
#                                         matched_glycan_dict_addr=__init__.output_matched_dict_addr)

    # nbt_table = pd.read_table(input_profile_addr)
    # a_profile = build_profiles(nbt_glycan_dict, nbt_table)


if __name__ == '__main__':
    # customizing_motif_vec_pip()
    pass
