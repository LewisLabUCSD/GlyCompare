import pandas as pd
import numpy as np
from matplotlib import pyplot as plt

from . import json_utility
from . import plot_glycan_utilities


def load_glycoprofile_mz_glycan_map():
    pass


def load_glycoprofile_naming_map(addr):
    f = open(addr)
    profile_naming_dict = {}
    for i in f.readline():
        profile, name, glycan_id = i.rstrip('\n').split('\t')
        if profile not in profile_naming_dict.keys():
            profile_naming_dict[profile] = {name: glycan_id}
        else:
            assert name not in profile_naming_dict[profile].keys()
            profile_naming_dict[profile][name] = glycan_id
    return profile_naming_dict


def output_glycoprofile_naming_map(glycoprofile_naming_dict, addr):
    file = open(addr, 'wr')
    for i in glycoprofile_naming_dict:
        for j in glycoprofile_naming_dict[i]:
            file.write('\t'.join([i, j, glycoprofile_naming_dict[i][j]]) + "\n")
    file.close()


def get_profile_str(glycan_dict, ez_vec):
    temp_dict = {}
    for i in ez_vec:
        if i not in glycan_dict.keys():
            print(i, 'not in glycan dict')
        else:
            temp_dict[i] = list(glycan_dict[i].keys())
    print("{")
    for i in sorted(list(temp_dict.keys())):
        print("\"" + str(i) + "\"" + ":\"" + " ".join(temp_dict[i]) + "\",")
    print("}")


def _plot_glycan_profile(a_profile, glycan_dict):
    _a = len(a_profile.keys())
    _len = 5
    _r = divmod(_a, 5)[0] + 1 if divmod(_a, 5)[1] != 0 else (divmod(_a, 5)[0])
    fig, axes = plt.subplots(_r, 5, squeeze=False)
    fig_size = (divmod(_a, 5)[0] + 1) * 3 if divmod(_a, 5)[1] != 0 else (divmod(_a, 5)[0]) * 3
    fig.set_size_inches(16, fig_size)
    #     print(axes)
    #     fig.set_size_inches(16,5)
    _count = 0
    for i in sorted(a_profile.keys()):
        #             print(i)
        #             print(divmod(_count, _a))
        _x, _y = divmod(_count, _len)
        plot_glycan_utilities.plot_glycan(glycan_dict[a_profile[i]], ax=axes[_x][_y], center=True, title=str(i))
        _count += 1


def merge_unzero_vec(prof_n, dict_name_abundance_cross_profile, glycan_dict, glycan_profile):
    """prof_n: start with 0 but the glycan profile start with 1
        skip the m/z with multiple glycans if the mass is small than 0.01
    """

    undoubt_list = ['3504']

    _id = str(prof_n + 1)
    _temp_dict = dict(glycan_profile[_id])
    for i in sorted(list(dict_name_abundance_cross_profile.keys())):
        if dict_name_abundance_cross_profile[i][prof_n] > 0.02:
            # print(glycan_profile[_id].keys())
            if i in glycan_profile[_id].keys():
                # print('\"' + i + '\": \"' + glycan_profile[_id][i] + '\", ')
                continue
            elif i in undoubt_list:
                _temp_dict[i] = list(glycan_dict[int(i)].keys())[0]
                # print('\"' + i + '\": \"' + glycan_profile[_id][i] + '\",')
                continue
            print(dict_name_abundance_cross_profile[i][prof_n], '\"' + i + '\":', list(glycan_dict[int(i)].keys()), ",")
        # continue
        #             glycan_profile[_id][i] = list(glycan_dict[int(i)].keys())[0]
        elif dict_name_abundance_cross_profile[i][prof_n] > 0.01:
            if i in glycan_profile[_id].keys():
                # print('\"' + i + '\": \"' + glycan_profile[_id][i] + '\",')
                continue
            elif i in undoubt_list:
                _temp_dict[i] = list(glycan_dict[int(i)].keys())[0]
                # print('\"' + i + '\": \"' + glycan_profile[_id][i] + '\",')
                continue
            if len(glycan_dict[int(i)].keys()) == 1:
                _temp_dict[i] = list(glycan_dict[int(i)].keys())[0]
                # else:
                #     if i in glycan_profile[_id].keys():
                #         print('\"' + i + '\": \"' + glycan_profile[_id][i] + '\",')
    return _temp_dict


def check_profile_naming_to_id(profile_naming_to_id):
    """
    :param profile_naming_to_id:
    :return:
    """

    key_list = list(profile_naming_to_id.keys())
    if type(key_list[0]) == str:
        return profile_naming_to_id
    elif type(key_list[0]) == int:
        return_dict = {}
        for i in profile_naming_to_id.keys():
            return_dict[str(i)] = profile_naming_to_id[i]
        return return_dict


def check_profile_naming_order(profile_naming_order):
    assert type(profile_naming_order) == list, 'check_profile_naming_order is not a list'
    if type(profile_naming_order[0]) == str:
        return profile_naming_order
    else:
        return [str(i) for i in profile_naming_order]


def check_external_profile_name(external_profile_name):
    assert type(external_profile_name) is dict, 'check_profile_naming_order is not a list'
    key_list = list(external_profile_name.keys())
    if type(key_list[0]) == str:
        return external_profile_name
    elif type(key_list[0]) == int:
        return_dict = {}
        for i in external_profile_name.keys():
            return_dict[str(i)] = external_profile_name[i]
        return return_dict


"""change the name of the glycan and glycan id"""


def get_glycoprofile_list(profile_naming_to_id, norm_mz_abd_dict, match_df, profile_name_order, external_profile_name,
                          glyprofile_list_addr,
                          absolute=False,
                          get_existance=True):
    """combine glycan m/z with substructure existances, convert the counts 1+/0 into 1/0 in each glycan
    Glycan substructures are modified with
    :param profile_name_order:
    :param get_existance:
    :param glyprofile_list_addr:
    :param profile_naming_to_id:
    :param norm_mz_abd_dict:
    :param match_df:
    return profile_obj_list
    """

    profile_naming_to_id = check_profile_naming_to_id(profile_naming_to_id)
    profile_name_order = check_profile_naming_order(profile_name_order)
    # print(profile_naming_to_id)
    if not external_profile_name:
        _ = list(profile_naming_to_id.keys())
        external_profile_name = dict(zip(_, _))
    else:
        external_profile_name = check_external_profile_name(external_profile_name)

    glycoprofile_list = []

    glycan_abd_dict = {}
    _num = len(profile_naming_to_id.keys())
    for i in match_df.columns:
        glycan_abd_dict[i] = _num * [0]

    for pro_idex, pro in enumerate(profile_name_order):
        # pro_id = pro
        weighted_vec = np.zeros(match_df.shape[0])
        abundance_list = []
        mz_list = []
        glycan_id_list = []
        match_mtrix = []

        for i in sorted(list(profile_naming_to_id[pro].keys())):
            glycan_id = profile_naming_to_id[pro][i]
            if glycan_id in match_df.columns:  
                mz_list.append(i)
                glycan_id_list.append(glycan_id)
                # print()
                _bundance = norm_mz_abd_dict[i][pro_idex]
                glycan_abd_dict[glycan_id][pro_idex] = _bundance
                _existance_list = []
                for _count in match_df[glycan_id]:
                    if _count >= 1:
                        if get_existance:
                            _existance_list.append(1)
                        else:
                            _existance_list.append(_count)
                    elif _count == 0:
                        _existance_list.append(0)
                    else:
                        assert False, 'wired in combine_profile_mz_with_substructure_existance'

                _temp_hit_matrix = np.array(_existance_list)

                abundance_list.append(_bundance)
                match_mtrix.append(_temp_hit_matrix)
        # print(abundance_list)

        for idex, i in enumerate(abundance_list):
                # print(abundance_list
                weighted_vec += match_mtrix[idex] * i



        glycoprofile_list.append(
            Glycoprofile(glycan_id_list, mz_list, abundance_list, weighted_vec, hit_matrix=match_mtrix,
                         name=pro, profile_name=external_profile_name[pro]))

    glycoprofile_output_list = []
    # print([round(i, 3) for i in merged_profile_dict[3]['substructure_vec'][:20]])
    for idex, i in enumerate(glycoprofile_list):
        glycoprofile_output_list.append(str(i.get_dict()))
    # print(glycoprofile_output_list[0])
    json_utility.store_json(glyprofile_list_addr, glycoprofile_output_list)
    # json_utility.store_json(addr_root + r"glycoprofile_list.json", glycoprofile_output_list)
    return glycoprofile_list


class Glycoprofile():
    """we will have the profile with topology
    dict self.profile
    list self.substructure_list
    profile contains glycan name and str structure, and substructure information
    we will have a
    """

    def __init__(self, glycan_id_list, mz_id_list, norm_abd, match_vec_weighted, hit_matrix=[], name=0,
                 profile_name=''):
        self.glycan_id_list = glycan_id_list
        self.mz_id_list = mz_id_list

        self.match_vec_weighted = match_vec_weighted
        self.match_vec_exist = []
        for i in match_vec_weighted:
            if i < 0.001:
                self.match_vec_exist.append(0)
            else:
                self.match_vec_exist.append(1)
        # print("len match_vec_exist", len(self.match_vec_exist))
        self.relative_abundance = norm_abd
        self.hit_matrix = hit_matrix
        self.name = name
        self.profile_name = profile_name

    def get_dict(self):
        rt_ = {'profile_name': self.profile_name,
               'm/z_list': self.mz_id_list,
               'glycan_id_list': self.glycan_id_list,
               'relative_abundance': list(self.relative_abundance),
               'match_vec_exist': self.match_vec_exist,
               'match_vec_weighted': list(self.match_vec_weighted),
               # 'hit_matrix': self.hit_matrix
               }
        return rt_

    def get_table(self):
        return pd.DataFrame({'m/z_list': self.mz_id_list,
                             'glycan_id_list': self.glycan_id_list,
                             'relative_abundance': list(self.relative_abundance),
                             'match_vec_weighted': list(self.match_vec_weighted)})



class substructureAbdTableGenerator:
    def __init__(self, glycoprofile_list):
        self.raw_table = glycoprofile_list[:]

    def table_btwn_two(self, n_a, n_b):
        #     change_f = np.array(a_p.weighted_substructure_vec)/np.array(b_p.weighted_substructure_vec)
        #     change_abs = np.array(a_p.weighted_substructure_vec)/np.array(b_p.weighted_substructure_vec)
        a_p = self.raw_table[n_a]
        b_p = self.raw_table[n_b]
        d = {a_p.name: a_p.match_vec_weighted, b_p.name: b_p.match_vec_weighted}
        table = pd.DataFrame(data=d)
        #     table = table[(table[a_p.name]+table[b_p.name])!=0]
        """find out the abundance of the N-glycan core and use it to balance the weight """
        table['r' + b_p.name] = table[b_p.name]
        table['fc'] = table['r' + b_p.name] / table[a_p.name]
        table['absc'] = table['r' + b_p.name] - table[a_p.name]
        return table

    # def get_table_with_target_index(self):
    #     _table = self.table_against_wt_relative_abd()
    #     for i in range(list(_table.index()):
    #         for ji in ra(_table.col)
    #         df.loc[row_indexer,column_indexer]

    def table_against_wt_relative_abd(self):
        """just generate raw table"""
        #     change_f = np.array(a_p.weighted_substructure_vec)/np.array(b_p.weighted_substructure_vec)
        #     change_abs = np.array(a_p.weighted_substructure_vec)/np.array(b_p.weighted_substructure_vec)
        a_p = self.raw_table[0]
        d = {a_p.name: a_p.match_vec_weighted}
        wt_table = pd.DataFrame(data=d)
        for idex, i in enumerate(self.raw_table):
            # if idex == 0:
            #     continue
            b_p = i
            # _weight = 1 / b_p.weighted_substructure_vec[7]
            wt_table[b_p.name] = np.array(b_p.match_vec_weighted)/sum(b_p.relative_abundance)
            #     wt_table = wt_table[(wt_table[a_p.name]+wt_table[b_p.name])!=0]
            """find out the abundance of the N-glycan core and use it to balance the weight """
        # wt_table = pd.DataFrame(data=d)
        return wt_table

    def table_absolute_abd(self):
        """just generate raw table"""
        #     change_f = np.array(a_p.weighted_substructure_vec)/np.array(b_p.weighted_substructure_vec)
        #     change_abs = np.array(a_p.weighted_substructure_vec)/np.array(b_p.weighted_substructure_vec)
        a_p = self.raw_table[0]
        d = {a_p.name: a_p.match_vec_weighted}
        wt_table = pd.DataFrame(data=d)
        for idex, i in enumerate(self.raw_table):
            b_p = i
            wt_table[b_p.name] = np.array(b_p.match_vec_weighted)
            #     wt_table = wt_table[(wt_table[a_p.name]+wt_table[b_p.name])!=0]
            """find out the abundance of the N-glycan core and use it to balance the weight """
        # wt_table = pd.DataFrame(data=d)
        return wt_table

    def table_existance(self):
        #     change_f = np.array(a_p.weighted_substructure_vec)/np.array(b_p.weighted_substructure_vec)
        #     change_abs = np.array(a_p.weighted_substructure_vec)/np.array(b_p.weighted_substructure_vec)
        a_p = self.raw_table[0]

        d = {a_p.name: a_p.match_vec_exist}
        existance_table = pd.DataFrame(data=d)
        for idex, i in enumerate(self.raw_table):
            # if idex == 0:
            #     continue
            b_p = i
            # _weight = 1 / b_p.weighted_substructure_vec[7]
            existance_table[b_p.name] = np.array(b_p.match_vec_exist)
            #     wt_table = wt_table[(wt_table[a_p.name]+wt_table[b_p.name])!=0]
            """find out the abundance of the N-glycan core and use it to balance the weight """
        # wt_table = pd.DataFrame(data=d)
        return existance_table

    def table_against_wt_fc(self):
        #     change_f = np.array(a_p.weighted_substructure_vec)/np.array(b_p.weighted_substructure_vec)
        #     change_abs = np.array(a_p.weighted_substructure_vec)/np.array(b_p.weighted_substructure_vec)

        wt_table = self.table_against_wt_relative_abd()
        wt_table += 0.001
        for idex, i in enumerate(wt_table.columns):
            if idex == 0:
                continue
            wt_table[i] = wt_table[i] / wt_table[wt_table.columns[0]]
            #     wt_table = wt_table[(wt_table[a_p.name]+wt_table[b_p.name])!=0]
            """find out the abundance of the N-glycan core and use it to balance the weight """
        # wt_table = pd.DataFrame(data=d)
        del wt_table[wt_table.columns[0]]

        return wt_table

    def table_against_wt_abs_val(self):
        #     change_f = np.array(a_p.weighted_substructure_vec)/np.array(b_p.weighted_substructure_vec)
        #     change_abs = np.array(a_p.weighted_substructure_vec)/np.array(b_p.weighted_substructure_vec)

        wt_table = self.table_against_wt_relative_abd()
        for idex, i in enumerate(wt_table.columns):
            if idex == 0:
                continue
            wt_table[i] = wt_table[i] - wt_table[wt_table.columns[0]]
            # wt_table = wt_table[(wt_table[a_p.name]+wt_table[b_p.name])!=0]
            """find out the abundance of the N-glycan core and use it to balance the weight """
        # wt_table = pd.DataFrame(data=d)
        del wt_table[wt_table.columns[0]]
        return wt_table

    def compare_two(self, a_n, b_n):
        wt_table = self.table_against_wt_relative_abd()
        l_a = wt_table.columns[a_n]
        l_b = wt_table.columns[b_n]
        col_a = wt_table[l_a] + 0.001
        col_b = wt_table[l_b] + 0.001

    def table_exist_or_not(self):
        wt_table = pd.DataFrame()
        for idex, i in enumerate(self.raw_table):
            _temp_vec = []
            for j in list(i.match_vec_weighted):
                if j == 0:
                    _temp_vec.append(0)
                else:
                    _temp_vec.append(1)
            wt_table[i.name] = _temp_vec
            #     wt_table = wt_table[(wt_table[a_p.name]+wt_table[b_p.name])!=0]
            """find out the abundance of the N-glycan core and use it to balance the weight """
        # wt_table = pd.DataFrame(data=d)
        # g = sns.clustermap(existed_table, metric="yule")
        return wt_table


if __name__ == '__main__':
    pass
