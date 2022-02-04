from glypy.algorithms.subtree_search.inclusion import subtree_of
from glypy.structure.glycan import fragment_to_substructure
import time
import multiprocessing

from . import glycan_io
from . import __init__
from . import json_utility
from . import pipeline_functions
import re
from glypy.io import wurcs
from glypy.io import glycoct
import json
import os
from datetime import datetime


#### proposed function

def extract_substructure_wurcs_idx(a_glycan,sub_gly, branch=5,linkage_specific=True):
    """
    :param a_glycan: Glycan obj
    :param sub_gly: dictionary of substructures. WURCS substructure index. Values are dictionary of matched glycans with the number of times each glycan appears.
    :param branch:
    :return:
    """
    if linkage_specific:
        gw = wurcs.dumps(a_glycan)
    else:
        # strip linkage information
        gct = glycoct.dumps(a_glycan)
        # replace all linkages with -1
        gct = re.sub('\(\d','(-1',gct)
        a_glycan = glycoct.loads(gct)
        gw = wurcs.dumps(a_glycan)

    # iterate over glycan fragments
    for i in a_glycan.fragments(max_cleavages=branch):
        # print('aaa')
        _frag_gly = wurcs.dumps(fragment_to_substructure(i, a_glycan))

        # print('aab')
        #add substructure 1st
        try:
            sub_gly[_frag_gly][gw] += 1
        except:
            try:
                sub_gly[_frag_gly][gw] = 1
            except:
                sub_gly[_frag_gly] = {}

reference_dict_ = ""
                
def extract_substructure(a_glycan, reference_dict, linkage_specific, branch=5):
    """
    :param a_glycan: Glycan obj
    :param branch:
    :return:
    """
#     global reference_dict_
    extracted_substructure_dic = {}
    ref_up = []
    for i in a_glycan.fragments(max_cleavages=branch):
        # print('aaa')
        _frag_gly = fragment_to_substructure(i, a_glycan)
        ref_ind, reference_dict = reference_get(glycoct.dumps(_frag_gly), reference_dict, linkage_specific = linkage_specific)
        ref_up.append((glycoct.dumps(_frag_gly), ref_ind))
        if not str(len(_frag_gly)) in extracted_substructure_dic.keys():
            extracted_substructure_dic[str(len(_frag_gly))] = [ref_ind]
        else:
            extracted_substructure_dic[str(len(_frag_gly))].append(ref_ind)
    # print('ab')
    ref_ind = reference_get(glycoct.dumps(a_glycan), reference_dict, linkage_specific = linkage_specific)[0]
    ref_up.append((glycoct.dumps(a_glycan), ref_ind))
    reference_dict = ""
    extracted_substructure_dic[str(len(a_glycan))] = [ref_ind]
    
    # print('ac')
    return extracted_substructure_dic, ref_up

def reference_get(glycan, reference_dict, linkage_specific = True):
#     global reference_dict_
#     reference_dict_ = reference_dict
    if glycan in reference_dict:
        return reference_dict[glycan], reference_dict
    else:
        if linkage_specific:
            reference_dict[glycan] = "L" + str(len(reference_dict))
        else:
            reference_dict[glycan] = "S" + str(len(reference_dict))
        return reference_dict[glycan], reference_dict
    
def reference_update(reference_dict, reference_dict_addr):
    old_reference_dict = json.load(open(reference_dict_addr, "r"))
    count = len(old_reference_dict) - len(reference_dict)
    os.remove(reference_dict_addr)
    num = str(len(reference_dict))
    time = "_".join(str(datetime.now()).split(".")[0].split(" "))
    time = time.replace(":", "-")
    dir_path = "/".join(reference_dict_addr.split("/")[:-1])
    target = reference_dict_addr.split("/")[-1].split(".json")[0]
    with open(dir_path + "/" + "_".join(target.split("_")[:4]) + "_" + num + "_" + time + ".json", 'w') as outfile:
        json.dump(reference_dict, outfile)
    print("Renewed reference dict, " + str(count) + " new motifs are added")

def extract_substructure_wrapper(a_name, a_glycan_str, substructure_dic, reference_dict, linkage_specific):
    """
    :param idex: idex of Glycan in the list
    :param a_name: the ID of the Glycan
    :param a_glycan_str: Glycan obj
    :param substructure_dic: dict for storing the extracted substructure dict
    :return:
    """

#     try:
    print('start', a_name)
    start_time = time.time()
    substructure_dic[a_name], ref_up = extract_substructure(a_glycan_str, reference_dict, linkage_specific, branch = 5)
    end_time = time.time()
    print(a_name, len(substructure_dic[a_name]), end_time - start_time)
    return ref_up
        # print('has_substructure', substructure_dic[_name])
        #     for j in substructure_dic.keys():
        #         print(j, len(substructure_dic[j]))
#     except TypeError:
#         print(a_name, 'has error')
#     except AttributeError:
#         print(a_name, 'has error')
#     except KeyboardInterrupt:
#         print('break')




def extract_substructures_pip(glycan_dict, gly_len, output_file, num_processors, reference_dict, linkage_specific, reference_dict_addr):
    """Please set the prior=True to get the data file please run the NBT_GLYCAN_preprocess file
    If prior=False, it will generate glycan substructure for all glycan in glytoucan database
    1. load  {glyacn_id: glycan_str}
    2. convert to glypy.glycan obj {glyacn_id: glycan_}
    3. extract {glyacn_id: substructure_dict}
    4. save glycan_substructure_dic
    :param num_processors:
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
        pool_list.append(pool.apply_async(extract_substructure_wrapper, args=(i, glycan_dict[i], substructure_dic, reference_dict, linkage_specific)))
    result_list = [xx.get() for xx in pool_list]
#     print(result_list)
    for gly in result_list:
        for sub in gly:
            if sub[0] not in reference_dict:
                reference_dict[sub[0]] = sub[1]
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
    reference_update(reference_dict, reference_dict_addr)
    return glycan_substructure_dic


# def main():
#     # glycan_dict = glycan_io.load_glycan_obj_from_database(topology_list_addr=__init__.topology_list_addr, output_file=__init__.glycan_dict_addr, loader=glycoct)
#     # glycan_substructure_dic = get_substructure_pip(glycan_dict=glycan_dict, gly_len=23, output_file=__init__.glycan_substructure_dict_addr)
#     pass
#
#
# if __name__ == '__main__':
#     main()
