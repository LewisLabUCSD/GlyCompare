import time
from glypy.algorithms.subtree_search import subtree_of
from glypy.structure.glycan import fragment_to_substructure
import multiprocessing
from gc_init import *
import gc_init
import json


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


def store_json(address, dic):
    with open(address, 'w') as fp:
        json.dump(dic, fp)


def load_json(address):
    with open(address, 'r') as f:
        return json.load(f)


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
#     glycan_profile = load_json(root_address + "BNT_glycan_profile.json")
#     BNT_glycan_match_existed_motif = load_json(root_address + "BNT_glycan_match_existed_motif.json")
#     dict_name_abundance_cross_profile = store_json(root_address + r"BNT_dict_name_abundance_cross_profile.json")
#
#     profile_obj_list = []
#     for i in range(1, 38):
#         if i < 10:
#             _id = 'Gly0' + str(i)
#         else:
#             _id = 'Gly' + str(i)
#             #     print(_id)
#         weighted_matrix = np.zeros((len(BNT_glycan_match_existed_motif["3865.1"])))
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
#             _temp_hit_matrix = np.array(BNT_glycan_match_existed_motif[_name])
#
#             # glycan_hit_array.append(_temp_hit_matrix)
#             abundance_.append(_bundance)
#             weighted_matrix += _temp_hit_matrix * _bundance
#         profile_obj_list.append(glycan_profile_obj(glycan_id_, mz_, abundance_, weighted_matrix))


def _duplicate_cleaning_wrapper(degree, motif_list, cleaned_motif_dic):
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
                if not subtree_of(_check_list[jdex], _check_list[ldex]) is None:
                    del _check_list[jdex]
                # elif not subtree_of(_check_list[ldex], _check_list[jdex]) is None:
                #     _check_list[ldex] = _check_list[ldex]
                #     del _check_list[jdex]
                else:
                    jdex += 1
            except AttributeError:
                print(ldex, jdex)
        # if not find_same:
        ldex += 1
    cleaned_motif_dic[degree] = _check_list


def merge_glycan_motif_dict_to_motif_dict(glycan_dict, combine_original=False):
    """

    :param glycan_dict: dict_ID_dict_degree_list
    :param combine_original: Check if you need add the original glycan to the motif-vec
    :return:
    """
    print("Start merge_glycan_motif_to_motif_dict")
    _motif_dic = {}
    for i in sorted(list(glycan_dict.keys()), reverse=True):
        # print(count__)
        for j in glycan_dict[i].keys():
            if j not in _motif_dic.keys():
                _motif_dic[j] = glycan_dict[i][j][:]
            else:
                _motif_dic[j].extend(glycan_dict[i][j][:])
    # print(_motif_dic.keys())
    if combine_original:
        print("combine original")
        glycan_dict_glycoct = glycan_str_to_glycan(
            load_json(glycoct_dict_goto_extraction_addr))
        for i in glycan_dict_glycoct.keys():
            print(type(glycan_dict_glycoct[i]), str(len(glycan_dict_glycoct[i])))
            if str(len(glycan_dict_glycoct[i])) not in _motif_dic.keys():
                print(len(glycan_dict_glycoct[i]))
                _motif_dic[str(len(glycan_dict_glycoct[i]))] = [glycan_dict_glycoct[i]]
            else:
                _motif_dic[str(len(glycan_dict_glycoct[i]))].append(glycan_dict_glycoct[i])
    return _motif_dic


def get_motif_dict_degree_list_pipe(glycan_dict, output_motif_dic_degree_list_addr):
    """
    merge the substructure of all glycans into motif dict
    :param glycan_dict: {degree: [motif1, motif2, ... ]} /BNT_glycan_dict_degree_list_glycoct_for_motif
    store the glycan motif to BNT_motif_dic_degree_list.json
    :return: sorted motif_vec
    """
    # output_motif_dic_degree_list_addr = root + "NBT_motif_dic_degree_list.json"

    _motif_dic = merge_glycan_motif_dict_to_motif_dict(glycan_dict, combine_original=False)
    print('check merged motif vec len', check_motif_dict_length(_motif_dic))
    print("get_motif_dict_degree_list_pipe")
    num_processors = 8
    pool = multiprocessing.Pool(processes=num_processors)
    _motif_degree = list(_motif_dic.keys())
    sorted_len_motif_degree = sorted(_motif_degree, key=lambda x: len(_motif_dic[x]), reverse=True)
    manager = multiprocessing.Manager()
    cleaned_motif_dic = manager.dict()
    for i in sorted_len_motif_degree:
        # print(i, len(_motif_dic[i]))
        pool.apply_async(_duplicate_cleaning_wrapper, args=(i, _motif_dic[i], cleaned_motif_dic))

    print("closing poll")
    pool.close()
    print('joining pool')
    pool.join()

    print('finished removing duplicate')
    cleaned_motif_ = dict(cleaned_motif_dic)

    print('after the cleaning the motif vec\'s length is', check_motif_dict_length(cleaned_motif_))

    str_motif = {}
    # print(cleaned_motif_.keys())
    for i in sorted_len_motif_degree:
        print(i, len(cleaned_motif_[i]))
        str_motif[i] = [str(j) for j in cleaned_motif_[i]]
    store_json(output_motif_dic_degree_list_addr, str_motif)
    return cleaned_motif_


def motif_matching(motif_vec, glycan_with_motif_dict, glycan_name, idex, match_dict):
    """

    :param motif_vec: customized motif vec
    :param glycan_with_motif_dict: extracted_motif_dic returned from the extract_motif
    :param glycan_name: ID glytoucan or other
    :param idex: idex in list (mainly for progress check)
    :param match_dict: dict for storing final data
    :return:
    """
    glycan_obj_dict = {}
    for i in glycan_with_motif_dict.keys():
        if type(glycan_with_motif_dict[i][0]) == str:
            glycan_obj_dict[i] = [glycoct.loads(k) for k in glycan_with_motif_dict[i]]
        else:
            glycan_obj_dict[i] = glycan_with_motif_dict[i]
    _existed_matrix = [0] * len(motif_vec)
    # _count_matrix = [0] * len(motif_vec)
    #     assert type(glycan_with_motif_dict[1][0])!=str
    # print(len(motif_vec))
    for jdex, j in enumerate(motif_vec):

        #             print(type(j), type(_tree))
        _len = len(j)
        if str(_len) not in glycan_obj_dict.keys(): continue
        #         print(_len)
        #         assert type(j) != str
        # print(len())
        _iter = 0
        # print(jdex)
        while _iter < len(glycan_obj_dict[str(_len)]):
            # print(type(j), type(i))

            if not subtree_of(j, glycan_obj_dict[str(_len)][_iter]) is None:
                # print('yes')
                _existed_matrix[jdex] += 1
                del glycan_obj_dict[str(_len)][_iter]
            else:
                _iter += 1

    # print('finished ', idex)
    print('finished ', idex)
    # print()
    match_dict[glycan_name] = _existed_matrix[:]


def motif_matching_wrapper(motif_dict, glycan_with_motif_dict, id_list, matched_glycan_dict_addr):
    """
    match the glycan_motif_dict_degree_list to motif vec
    :param motif_dict: {degree: [motif1, motif2, ...]}  /NBT_motif_dic_degree_list
    :param glycan_with_motif_dict: {degree: [motif1, motif2, ...]} /BNT_glycan_dict_degree_list_glycoct_for_motif
    :param id_list: all GlytoucanID of the glycans you are analyzing /NBT_fixed_gylcan_name_list
    :param matched_glycan_dict_addr: output_addr /NBT_fixed_gylcan_name_list
    :return: glycan_match_existed_motif degree - >[motif1, motif2, ...]
    """
    if type(motif_dict) == dict:
        motif_vec = motif_dict_to_vec(motif_dict)
    else:
        motif_vec = motif_dict
    # store_json(output_motif_vec_addr, [str(i) for i in motif_vec])
    print('get motif vec, the length is ', len(motif_vec))
    pool = multiprocessing.Pool(processes=num_processors)
    manager = multiprocessing.Manager()
    match_dict = manager.dict()  # {motif_name:[scores]}
    # print(type(glycan_with_motif_dict['G85809SI']['1'][0]))
    for idex, i in enumerate(id_list):
        print('start processing', i)
        pool.apply_async(motif_matching, args=(motif_vec, glycan_with_motif_dict[i], i, idex, match_dict))
    print("closing poll")
    pool.close()
    print('joining pool')
    pool.join()
    print('converting dict')
    a_motif_dic = dict(match_dict)
    # print(a_motif_dic)
    store_json(matched_glycan_dict_addr, a_motif_dic)
    return a_motif_dic


def motif_dict_to_vec(motif_dict):
    motif_vec = []
    for i in sorted([int(j) for j in list(motif_dict.keys())]):
        print(i, len(motif_dict[str(i)]))
        motif_vec.extend(motif_dict[str(i)])
    print(len(motif_vec))
    return motif_vec


def compare_profile(profile_a, profile_b):
    pass


def build_profiles(glycan_dict, glycan_profile):
    """
    :param glycan_dict: with hit
    :param glycan_profile: 2D matrix with glycan relative abundance
    :return: a dict {profile1: [(glycan1, abundance1), (g2, ab2)],
                     profile2: [(glycan1, abundance1), (g2, ab2)]}
    """
    pass


def check_motif_dict_length(a_dict):
    return sum([len(a_dict[i]) for i in a_dict.keys()])


def customizing_motif_vec_pip():
    """ merge the substructure of all glycans into motif dict"""
    BNT_glycan_dict_degree_list_glycoct_for_motif = glycan_str_to_glycan(
        load_json(glycan_dict_motif_list_addr))
    # output_motif_dic_degree_list_addr = root + "BNT_motif_dic_degree_list.json"
    NBT_motif_dic_degree_list = get_motif_dict_degree_list_pipe(BNT_glycan_dict_degree_list_glycoct_for_motif,
                                                                output_motif_dic_degree_list_addr)
    # NBT_motif_dic_degree_list = glycan_str_to_glycan(load_json(root+"BNT_motif_dic_degree_list.json"))

    """ Start motif matching"""
    # print('start motif match')
    # BNT_glycan_dict_degree_list_glycoct_for_motif = glycan_str_to_glycan(load_json(glycan_dict_motif_list_addr))
    BNT_motif_dic_degree_list = glycan_str_to_glycan(load_json(output_motif_dic_degree_list_addr))

    NBT_fixed_gylcan_name_list = load_json(NBT_fixed_gylcan_name_list_addr)

    # NBT_fixed_gylcan_name_list = load_json(root + "BNT_fixed_gylcan_name.json")
    glycan_match_existed_motif = motif_matching_wrapper(NBT_motif_dic_degree_list,
                                                        BNT_glycan_dict_degree_list_glycoct_for_motif,
                                                        NBT_fixed_gylcan_name_list, output_matched_glycan_addr)

    # glycan_match_existed_motif = motif_matching_wrapper(nbt_motif_dict, nbt_motif_dic_degree_list, nbt_fixed_gylcan_name_list, output_matched_glycan_addr)

    # nbt_table = pd.read_table(input_profile_addr)
    # a_profile = build_profiles(nbt_glycan_dict, nbt_table)


if __name__ == '__main__':
    customizing_motif_vec_pip()
