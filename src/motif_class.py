# break glycoCT
import time
from glypy.algorithms.subtree_search import subtree_of
from glypy.structure.glycan import fragment_to_substructure
import seaborn as sns
from scipy.spatial import distance
from glypy.io import glycoct
import __init__
import numpy as np
from scipy import stats

sns.set(color_codes=True)
# glyco_motif_list={}
# glycoct_list = []
profile_name = ['WT',
                'mgat4A',
                'mgat4A/mgat4B',
                'mgat5',
                'mgat4A/mgat4B/mgat5',
                'B4GalT1',
                'B4GalT2',
                'B4GalT3',
                'B4GalT4',
                'B4GalT1/B4GalT2',
                'B4GalT1/B4GalT3',
                'B3gnt1',
                'B3gnt2',
                'st3gal3',
                'st3gal4',
                'st3gal6',
                'st3gal3/st3gal4',
                'st3gal4/st3gal6',
                'KI_ST6GalNAc1/st3gal4/st3gal6',
                'B3gnt2/mgat4a/mgat4b/mgat5',
                'st3gal4/st3gal6/mgat4a/mgat4b/mgat5',
                'KI_ST6GalNAc1/st3gal4/st3gal6/mgat4a/mgat4b/mgat5',
                'EPO48(mgat3)',
                'EPO143(mgat4C)',
                'EPO174(mgat2)',
                'EPO200(B4galt1/B4galt2/B4galt3)',
                'EPO275(B3gnt8)',
                'EPO78(mgat4B)',
                'EPO104(mgat5B)',
                'EPO127(mgat1)',
                'EPO259(mgat2/st3gal4/st3gal6)',
                'EPO261(mgat2/mgat4A/mgat4B/mgat5)',
                'EPO263(mgat2/st3gal4/st3gal6/magt4A/mgat4B/mgat5)',
                'EPO266(fut8)']


# len(aaa)

#
# def get_motif(glycoct_obj, idex=0):
#     # print('start getmotif')
#     _frag_motif_list = {}
#     # fig, axes = plt.subplots(6,9)
#     # fig.set_size_inches(14,6)
#     start_time = time.time()
#     for i in glycoct_obj.fragments(max_cleavages=len(glycoct_obj)):
#         _frag_gly = fragment_to_substructure(i, glycoct_obj)
#
#         # plot(_frag_gly)
#         if not len(_frag_gly) in _frag_motif_list.keys():
#             _frag_motif_list[len(_frag_gly)] = [glycoct.loads(_frag_gly)]
#         else:
#             _frag_motif_list[len(_frag_gly)].append(glycoct.loads(_frag_gly))
#     mid_time = time.time()
#     # print('start clean duplicate')
#     _frag_motif_list = clean_duplicate(_frag_motif_list)
#     end_time = time.time()
#     print(idex, len(glycoct_obj), end_time - mid_time, mid_time - start_time)
#     # print('finished getmotif')
#
#     return _frag_motif_list


def clean_duplicate(_frag_motif_list):
    for i in _frag_motif_list.keys():
        # print(i)
        ldex = 0
        _check_list = _frag_motif_list[i]
        while ldex < len(_check_list):
            jdex = ldex + 1
            while jdex < len(_check_list):
                if subtree_of(_check_list[ldex], _check_list[jdex], __init__.exact_Ture) == 1 and subtree_of(
                        _check_list[jdex],
                        _check_list[ldex]) == 1:
                    del _check_list[jdex]
                else:
                    jdex += 1
            # if not find_same:
            ldex += 1
        _frag_motif_list[i] = _check_list
    return _frag_motif_list


# glytoucan_data_base = load_json(r'/Users/apple/PycharmProjects/glycan_within_all_lectins.json')
# # for i in glytoucan_data_base.keys():
# #     glycoct_list.append(glycoct.loads(glytoucan_data_base[i]['GlycoCT']))
# #     if len(glycoct.loads(glytoucan_data_base[i]['GlycoCT']))==10:
# #         print(i)
# # # print(set([len(i) for i in glycoct_list]))
# ten_glycan = glycoct.loads(glytoucan_data_base['G28566CQ']['GlycoCT'])
# print(get_motif(ten_glycan))

class MotifLab():
    _man1 = glycoct.loads("""
        RES
        1b:b-dman-HEX-1:5
        LIN""")

    def __init__(self, motif_):
        if type(motif_) == dict:
            print(type(list(motif_.keys())[0]))
            dict_keys = sorted([int(i) for i in motif_.keys()])
            self.motif_vec = []
            for i in dict_keys:
                for j in motif_[str(i)]:
                    if isinstance(j, type(self._man1)):
                        self.motif_vec.append(j)
                    else:
                        self.motif_vec.append(glycoct.loads(j))
            self.motif_dict = {}
            for idex, i in enumerate(self.motif_vec):
                if len(i) not in self.motif_dict.keys():
                    self.motif_dict[len(i)] = [idex]
                else:
                    self.motif_dict[len(i)].append(idex)
        elif type(motif_) == list:
            if isinstance(motif_[0], type(self._man1)):
                self.motif_vec = motif_
            else:
                self.motif_vec = [glycoct.loads(i) for i in motif_]
            self.motif_dict = {}
            for idex, i in enumerate(self.motif_vec):
                if len(i) not in self.motif_dict.keys():
                    self.motif_dict[len(i)] = [idex]
                else:
                    self.motif_dict[len(i)].append(idex)
        else:
            assert False, "should be either list or dict"
        self.motif_list = [i for i in range(len(self.motif_vec))]
        self.motif_dep_tree = {}

    def dep_tree_to_edge_list(self, dep_tree):
        """

        :param dep_tree: motif_with_ncore_dependence_tree, motif_dependence_tree, motif_single_dependence_tree
        :return: edge_list
        """
        edge_list = []
        for i in dep_tree:
            for k in dep_tree[i]:
                edge_list.append((i, k))
        return edge_list

    def build_dependence_tree(self, a_motif_dict):
        """ connect motif to all parents"""
        print('start building dependence_tree')
        edge_list = []
        _dep_tree = {}
        # self.motif_dep_tree = {}
        for i in sorted(list(a_motif_dict.keys())):
            print(i)
            if i - 1 not in a_motif_dict.keys():
                for j in a_motif_dict[i]:
                    _dep_tree[j] = []
                continue
            for j in a_motif_dict[i]:
                """
                    motif j in i degree/motif in i-1 degree
                """
                _dep_tree[j] = []
                for k in a_motif_dict[i - 1]:
                    if subtree_of(self.motif_vec[k], self.motif_vec[j], exact=__init__.exact_Ture) == 1:
                        _dep_tree[k].append(j)
                        edge_list.append((k, j))
        return _dep_tree, edge_list

    def get_dependence_tree_all(self):
        """
        get the dep tree for all node
        :return: dep_tree, edge_list
        """
        if self.motif_dep_tree == {}:
            _dep_tree, _edge_list = self.build_dependence_tree(self.motif_dict)
            self.motif_dep_tree = _dep_tree
            return _dep_tree, _edge_list
        else:
            return self.motif_dep_tree, self._get_edge(self.motif_dep_tree)

    def _get_edge(self, dep_tree):
        edge_list = []
        for i in dep_tree:
            edge_list.extend([(i, j) for j in dep_tree[i]])
        return edge_list


class MotifLabNGlycan(MotifLab):
    """
    store vec
    """
    nglycan_core = glycoct.loads(
        "RES\n1b:x-dglc-HEX-1:5\n2s:n-acetyl\n3b:b-dglc-HEX-1:5\n4s:n-acetyl\n5b:b-dman-HEX-1:5\n6b:a-dman-HEX-1:5"
        "\n7b:a-dman-HEX-1:5\nLIN\n1:1d(2+1)2n\n2:1o(4+1)3d\n3:3d(2+1)4n\n4:3o(4+1)5d\n5:5o(3+1)6d\n6:5o(6+1)7d\n ")
    _with_sia_core = glycoct.loads("""RES
        1b:b-dglc-HEX-1:5
        2s:n-acetyl
        3b:b-dgal-HEX-1:5
        4b:a-dgro-dgal-NON-2:6|1:a|2:keto|3:d
        5s:n-acetyl
        LIN
        1:1d(2+1)2n
        2:1o(4+1)3d
        3:3o(3+2)4d
        4:4d(5+1)5n""")
    _no_sia_core = glycoct.loads("""RES
        1b:b-dglc-HEX-1:5
        2s:n-acetyl
        3b:b-dgal-HEX-1:5
        LIN
        1:1d(2+1)2n
        2:1o(4+1)3d""")
    _man2 = glycoct.loads("""
        RES
        1b:a-dman-HEX-1:5
        LIN""")

    def __init__(self, motif_):
        """
        self.motif_dict stores the id of the self.motif_vec
        :param motif_: motif vec or motif dict_degree_list:
        """
        assert motif_, "motif vector is empty"
        MotifLab.__init__(self, motif_)
        self.motif_dict_with_core = {}
        self.motif_dep_tree_ncore = {}
        self.motif_with_core_list = []
        self.extract_motif_with_core()

    #     self.motif_vec_sia_ept = []
    #     self.motif_vec_gala_ept = []
    #
    # def create_epitope_vec(self):
    #     print("start motif with sia")
    #     if not self.motif_vec_sia_ept:
    #         for i in sorted(list(self.motif_dict.keys())):
    #             # if (i) > 9:
    #             #     break
    #             # print("len", i)
    #             for j in self.motif_dict[i]:
    #                 """
    #                 motif j in i degree/motif in i-1 degree
    #                 """
    #                 if subtree_of(self._man1, self.motif_vec[j], exact=__init__.exact_Ture) is not None or subtree_of(
    #                         self._man2, self.motif_vec[
    #                             j], exact=__init__.exact_Ture) is not None:
    #                     continue
    #                 if subtree_of(self._with_sia_core, self.motif_vec[j], exact=__init__.exact_Ture) is not None:
    #                     if len(self.motif_vec[j]) % 2 == 1:
    #                         self.motif_vec_sia_ept.append(j)
    #                         # self.gala_ept_vec.append(j)
    #                 elif subtree_of(self._no_sia_core, self.motif_vec[j], exact=__init__.exact_Ture) is not None:
    #                     if len(self.motif_vec[j]) % 2 == 0:
    #                         self.motif_vec_gala_ept.append(j)
    #
    #         print("Finish sia match ", len(self.motif_vec_sia_ept),
    #               " motifs are find with sia core ", len(self.motif_vec_gala_ept), " motifs are find with no sia core ")
    #     else:
    #         print("Finish sia match ", len(self.motif_vec_sia_ept),
    #               " motifs are find with sia core ", len(self.motif_vec_sia_ept), " motifs are find with no sia core ")

    def extract_motif_with_core(self):
        """ store the result in self.motif_with_core_list
        and return the count"""
        # count = []
        print("start motif_with core")
        if self.motif_dict_with_core == {}:
            for i in sorted(list(self.motif_dict.keys())):
                if len(self.nglycan_core) > i:
                    continue
                print("len", i)
                self.motif_dict_with_core[i] = []
                for j in self.motif_dict[i]:
                    """
                    motif j in i degree/motif in i-1 degree
                    """
                    if subtree_of(self.nglycan_core, self.motif_vec[j], exact=__init__.exact_Ture) == 1:
                        self.motif_dict_with_core[i].append(j)
                        self.motif_with_core_list.append(j)
            print("Finish the n-glycan match ", len(self.motif_with_core_list),
                  " motifs are matched to the n-glycan core")
        else:
            print("Finish the n-glycan match ", len(self.motif_with_core_list),
                  " motifs are matched to the n-glycan core")

    def get_dependence_tree_ncore(self):
        """
        get the dep tree for ncore's node
        :return: dep_tree, edge_list
        """
        if self.motif_dep_tree_ncore == {}:
            _dep_tree, _edge_list = self.build_dependence_tree(self.motif_dict_with_core)
            self.motif_dep_tree_ncore = _dep_tree
            return _dep_tree, _edge_list
        else:
            return self.motif_dep_tree_ncore, self._get_edge(self.motif_dep_tree_ncore)


def get_weight_dict(motif_abd_table):
    _np_mat = np.array(motif_abd_table)
    weight_dict = {}
    for i in range(len(_np_mat)):
        weight_dict[i] = list(_np_mat[i])
    return weight_dict


class NodesState():
    threshold = 200

    def __init__(self, dependence_tree, motif_weight):
        self.dep_tree = dependence_tree
        self.nodes_sta = []
        # self.nodes = self._get_node(self.dep_tree)
        self.edges = self._get_edge(dependence_tree)
        self.nodes = self._get_node(dependence_tree)

        self.motif_weight = motif_weight
        self.normalized_motif_weight = {}
        self._normalized_weight()
        self.parents_dic = {}
        for i in dependence_tree:
            self.parents_dic[i] = {}
        self.nodes_kept = []
        self._out_degree_list = []
        self._in_degree_list = []
        self.flat_paired_diff = []
        self.flat_normed_paired_diff = []
        self.nodes_kept = []
        self.nodes_stat_value = []

    def _get_node(self, dep_tree):
        """
        get nodes from dep_tree
        :param dep_tree:
        :return:
        """
        node_list = list(dep_tree.keys())
        for i in dep_tree:
            node_list.extend(dep_tree[i])
        return sorted(list(set(node_list)))

    def _get_edge(self, dep_tree):
        """
        get edge from dep_tree
        :param dep_tree:
        :return:
        """
        edge_list = []
        for i in dep_tree:
            edge_list.extend([(i, j) for j in dep_tree[i]])
        return edge_list

    # wrapper
    # def get_edge_ttest_dis(self):
    #     _list = []
    #     for i, j in self.edge:
    #         _list.append(self.get_value(i, j, method=self.one_vs_rest_t))
    #     return _list

    # wrapper
    def _no_num(self, _list, _num):
        _len = len(_list)
        count_ = 0
        while _num in _list:
            _list.remove(_num)
            count_ += 1
        print("there are ", count_, " removed from ", _len)
        return _list

    def get_edge_ttest_dis(self):
        _list = []
        for i, j in self.edges:
            _list.append(self.get_value_unnormed(i, j, method=self.one_vs_rest_t))
        return _list

    # wrapper
    def get_edge_corr_dis(self):
        _list = []
        for i, j in self.edges:
            _list.append(self.get_value_unnormed(i, j, method=self.get_corr))
        return _list

    def _normalized_weight(self):
        for i in self.motif_weight.keys():
            _max = max(self.motif_weight[i])
            self.normalized_motif_weight[i] = [j / _max for j in self.motif_weight[i]]

    def get_vector(self, i):
        return self.normalized_motif_weight[i]

    def get_value_unnormed(self, i, j, method):
        # print()
        return method(vec_a=self.motif_weight[i], vec_b=self.motif_weight[j])

    #
    # def get_value_normed(self, i, j, method):
    #     # print()
    #     return method(vec_a=self.normalized_motif_weight[j], vec_b=self.normalized_motif_weight[i])

    def get_corr(self, vec_a, vec_b):
        return 1 - distance.braycurtis(vec_a, vec_b)

    def one_vs_rest_t(self, vec_a, vec_b):
        diff_vec = [vec_a[i] - vec_b[i] for i in range(len(vec_a))]
        # for i in range(len(vec_a)):
        _min = min(diff_vec)
        _temp_vec = diff_vec[:]
        _temp_vec.remove(_min)
        neg_log_p_min = self.get_neg_log_p_ttest(_min, _temp_vec)
        _max = max(diff_vec)
        _temp_vec = diff_vec[:]
        _temp_vec.remove(_max)
        neg_log_p_max = self.get_neg_log_p_ttest(_max, _temp_vec)

        return_neg_log_p = max(neg_log_p_min, neg_log_p_max)
        if return_neg_log_p > self.threshold:
            # print(_min, diff_vec)
            return_neg_log_p = self.threshold
        return return_neg_log_p

    def all_one_vs_rest_t(self, vec_a, vec_b):
        diff_vec = [vec_a[i] - vec_b[i] for i in range(len(vec_a))]
        # print(diff_vec)
        # for i in range(len(vec_a)):
        _return_list = []
        for i in diff_vec:
            _temp_vec = diff_vec[:]
            _temp_vec.remove(i)
            # print(i, _temp_vec)
            if np.var(_temp_vec) == 0:
                return []
            neg_log_p = self.get_neg_log_p_ttest(i, _temp_vec)
            if neg_log_p > self.threshold:
                # print(_min, diff_vec)
                neg_log_p = self.threshold
            _return_list.append(neg_log_p)
        return _return_list

    # def get_weight_list(self):
    #     for i in self.normalized_motif_weight:
    #         for j in self.normalized_motif_weight[i]:
    #             _temp_vec = diff_vec[:]
    #             _temp_vec.remove(i)

    def get_edge_all_ttest(self):
        _list = []
        self.flat_normed_paired_diff = []
        self.flat_paired_diff = []
        for i, j in self.edges:
            _temp = self.get_value_unnormed(i, j, method=self.all_one_vs_rest_t)
            if _temp:
                _list.extend(_temp)
                self.flat_paired_diff.extend(
                    [self.motif_weight[i][k] - self.motif_weight[j][k] for k in range(len(self.motif_weight[j]))])
                self.flat_normed_paired_diff.extend(
                    [self.normalized_motif_weight[i][k] - self.normalized_motif_weight[j][k] for k in
                     range(len(self.normalized_motif_weight[j]))])
        return _list

    def get_neg_log_p_ttest(self, _ele, _vec):
        _vec = np.array(_vec)
        if _vec.var() == 0:
            if _ele - _vec.mean() == 0:
                return 0
            else:
                return self.threshold
        else:
            # print(_vec.mean(), _vec.var())
            tt = (_ele - _vec.mean()) / np.sqrt(_vec.var() / len(_vec))
            scipy_tt = stats.ttest_1samp(_vec, _ele)
            # print(type(scipy_tt.pvalue))
            p = scipy_tt.pvalue
        # p = stats.t.sf(np.abs(tt), len(_vec) - 1)*2
        return -np.log(p)

    # def get_neg_log_p_ttest_sci(self, _ele, _vec):
    #     _vec = np.array(_vec)
    #     if _vec.var() == 0:
    #         if _ele - _vec.mean() == 0:
    #             return 0
    #         else:
    #             return self.threshold
    #     else:
    #         # print(_vec.mean(), _vec.var())
    #         tt = (_ele - _vec.mean()) / np.sqrt(_vec.var() / len(_vec))
    #
    #     p = stats.t.sf(np.abs(tt), len(_vec) - 1)*2
    #     return -np.log(p)


    # def one_vs_rest_oneside_t(self, vec_a, vec_b):
    #     diff_vec = [vec_a[i] - vec_b[i] for i in range(len(vec_a))]
    #     print(diff_vec)
    #     # for i in diff_vec:
    #     #     if i < 0:
    #     #         print()
    #             # assert False
    #     # for  i in range(len(vec_a)):
    #     _min = min(diff_vec)
    #     _temp_vec = diff_vec[:]
    #     _temp_vec.remove(_min)
    #     neg_log_p_min = self.get_neg_log_p_ttest(_min, _temp_vec)
    #
    #     if neg_log_p_min > 100:
    #         # print(_min, diff_vec)
    #         neg_log_p_min = 100
    #     assert 0 <= neg_log_p_min <= 100
    #     # elif p_min == float('inf'):
    #     #     print(_min, diff_vec)
    #     # if -np.log(p_min) > 1000 or -np.log(p_min) < -1000:
    #     #     print(_min, diff_vec)
    #     return neg_log_p_min

    def get_node_value(self, redo=True):
        if redo:
            self.nodes_stat_value = []
            self._in_degree_list = []
            self._out_degree_list = []
            self.nodes_sta = []
            for i in sorted(list(self.parents_dic.keys())):
                # print(i)
                for j in self.dep_tree[i]:
                    self.parents_dic[j][i] = (self.get_value_unnormed(i, j, self.get_corr),
                                              self.get_value_unnormed(i, j, self.one_vs_rest_t))
                self._out_degree_list.append(len(self.dep_tree[i]))

            # can be split
            for i in sorted(list(self.parents_dic.keys())):

                self._in_degree_list.append(len(self.parents_dic[i]))
                _list = []
                for j in self.parents_dic[i].keys():
                    _list.append(self.parents_dic[i][j])
                if _list:
                    # self.nodes_kept.append(max(_list))
                    _list = sorted(_list, key=lambda x: x[0])
                    self.nodes_stat_value.append(_list[0])
                else:
                    self.nodes_stat_value.append((-1, -1))

            return zip(*self.nodes_stat_value)

    def get_node_sta(self):
        self.nodes_sta = []
        for i in sorted(list(self.parents_dic.keys())):
            _mean = np.mean(self.motif_weight[i])
            _var = np.var(self.motif_weight[i])
            self.nodes_sta.append((_mean, _var))
        return zip(*self.nodes_sta)

    def get_edge_node_degree(self, edges=[]):
        out_degree = {}
        in_degree = {}
        out_degree_list = []
        in_degree_list = []
        for i in self.nodes:
            out_degree[i] = 0
            in_degree[i] = 0
        if not edges:
            edges = self.edges
        for i, j in edges:
            out_degree[i] += 1
            in_degree[j] += 1
        for i, j in edges:
            out_degree_list.append(out_degree[i])
            in_degree_list.append(in_degree[j])
        return out_degree_list, in_degree_list

    def get_nodes_degree(self, edges=[]):
        out_degree = {}
        in_degree = {}
        out_degree_list = []
        in_degree_list = []
        for i in self.nodes:
            out_degree[i] = 0
            in_degree[i] = 0
        if not edges:
            edges = self.edges
        for i, j in edges:
            out_degree[i] += 1
            in_degree[j] += 1

        return out_degree, in_degree

    def nodes_dropper(self):
        """drop nodes based on features"""


class NodesDropper:
    def __init__(self, dependence_tree, motif_weight):
        self.dep_tree = dependence_tree
        self.all_nodes = self._get_node(dependence_tree)

        self.parents_vec = {}
        for i in dependence_tree:
            self.parents_vec[i] = {}

        self.motif_weight = motif_weight
        self.normalized_motif_weight = {}
        self._normalized_weight()
        self.heavy_dependency = {}
        self.most_dependent_child = {}
        self.nodes_kept = []

    def _get_node(self, dep_tree):
        node_list = list(dep_tree.keys())

        for i in dep_tree:
            node_list.extend(dep_tree[i])
        return sorted(list(set(node_list)))
        # self.gala_ept_vec = a_glycan_motif_lib.gala_ept_vec[:]
        # self.sia_ept_vec = a_glycan_motif_lib.sia_ept_vec[:]
        # self.sia_gala_ept_vec = self.gala_ept_vec[:]
        # self.sia_gala_ept_vec.extend(self.sia_ept_vec[:])

    def _normalized_weight(self):
        for i in self.motif_weight.keys():
            _max = max(self.motif_weight[i])
            self.normalized_motif_weight[i] = [j / _max for j in self.motif_weight[i]]

    def _z_score(self, vec_a, vec_b):
        diff_vec = [vec_a[i] - vec_b[i] for i in range(len(vec_a))]
        # for i in range(len(vec_a)):
        _min = min(diff_vec)
        _temp_vec = diff_vec[:]
        _temp_vec.remove(_min)
        p_min = self.ttest_wrapper(_min, _temp_vec)
        _max = max(diff_vec)
        _temp_vec = diff_vec[:]
        _temp_vec.remove(_max)
        p_max = self.ttest_wrapper(_max, _temp_vec)
        return min(p_min, p_max)

    def ttest_wrapper(self, _ele, _vec):
        _vec = np.array(_vec)
        tt = (_ele - _vec.mean()) / np.sqrt(_vec.var() / len(_vec))
        p = stats.t.sf(np.abs(tt), len(_vec) - 1) * 2
        return p

    # def get_sia_gal_vec(self):
    #     rt_lst = []
    #     rt_lst.extend(self.gala_ept_vec[:])
    #     rt_lst.extend(self.sia_ept_vec[:])
    #     return rt_lst

    # drop_list = add_sia_gal(a_dp_tree)
    #
    # def get_drop_node_with_sia(self):
    #     rt_lst = self.drop_node(distance.correlation)
    #     rt_lst.extend(self.sia_gala_ept_vec)
    #     return rt_lst

    def compare_abundance(self):
        pass

        # def drop_node_with_parents(self):
        #     for i in self.parents_vec.keys():
        #         find_ = False
        #         for j in self.parents_vec[i].keys():
        #             if self.parents_vec[i][j] > 0.999:
        #                 find_ = True
        #                 if j not in self.heavy_dependency.keys():
        #                     self.heavy_dependency[j] = [i]
        #                 else:
        #                     self.heavy_dependency[j].append(i)
        #         if find_: continue
        # self.after_drop_node.append(i)

    # def _normalized_weight(self):
    #     for i in self.motif_weight.keys():
    #         _max = max(self.motif_weight[i])
    #         self.normalized_motif_weight[i] = [j / _max for j in self.motif_weight[i]]

    def drop_node(self, method=distance.braycurtis, redo=True):
        if self.nodes_kept != [] and redo is False:
            return self.nodes_kept
        else:
            self.nodes_kept = []
            for i in self.parents_vec.keys():
                # print(i)
                for j in self.dep_tree[i]:
                    self.parents_vec[j][i] = 1 - method(self.normalized_motif_weight[j],
                                                        self.normalized_motif_weight[i])
                    # sns.clustermap
            for i in self.parents_vec.keys():
                find_ = False
                _max = 0
                _max_parent = ''
                for j in self.parents_vec[i].keys():

                    if self.parents_vec[i][j] > 0.995:
                        # if _max_parent == '':
                        #     print(i,j,"find", _max,self.parents_vec[i][j],_max < self.parents_vec[i][j])
                        if _max < self.parents_vec[i][j]:
                            _max = self.parents_vec[i][j]
                            _max_parent = j
                            # print(j,_max)
                        # elif
                        # if self.parents_vec[i][j] < 0.999999:
                        #     print(self.parents_vec[i][j], i, j)
                        # elif self.parents_vec[i][j] ==1:
                        #     print('100', i,j)
                        find_ = True

                if find_:
                    # if _max_parent == '':
                    #     print('wtf',i,_max_parent,_max)
                    if _max_parent not in self.heavy_dependency.keys():
                        # if _max_parent == "":
                        #     print(_max)
                        self.heavy_dependency[_max_parent] = [i]
                    else:
                        self.heavy_dependency[_max_parent].append(i)

                else:
                    self.nodes_kept.append(i)
            return self.nodes_kept

    def drop_node_with_t_test(self, method=distance.braycurtis, redo=True):
        if self.nodes_kept != [] and redo is False:
            return self.nodes_kept
        else:
            self.nodes_kept = []
            for i in self.parents_vec.keys():
                # print(i)
                for j in self.dep_tree[i]:
                    self.parents_vec[j][i] = 1 - method(self.normalized_motif_weight[j],
                                                        self.normalized_motif_weight[i])
                    # sns.clustermap
            for i in self.parents_vec.keys():
                find_ = False
                _max = 0
                _max_parent = ''
                for j in self.parents_vec[i].keys():
                    if self.parents_vec[i][j] > 0.995:
                        # if _max_parent == '':
                        #     print(i,j,"find", _max,self.parents_vec[i][j],_max < self.parents_vec[i][j])
                        if _max < self.parents_vec[i][j]:
                            _max = self.parents_vec[i][j]
                            _max_parent = j
                            # print(j,_max)
                        # elif
                        # if self.parents_vec[i][j] < 0.999999:
                        #     print(self.parents_vec[i][j], i, j)
                        # elif self.parents_vec[i][j] ==1:
                        #     print('100', i,j)
                        find_ = True
                if find_:
                    # if _max_parent == '':
                    #     print('wtf',i,_max_parent,_max)
                    if _max_parent not in self.heavy_dependency.keys():
                        # if _max_parent == "":
                        #     print(_max)
                        self.heavy_dependency[_max_parent] = [i]
                    else:
                        self.heavy_dependency[_max_parent].append(i)
                else:
                    self.nodes_kept.append(i)
            return self.nodes_kept
            #
            # def get_the_most_dependent_node(self):
            #     """heavy dependency parents child_lst
            #         parents_vec child -> parents
            #     """
            #     for j, i_list in self.heavy_dependency.items():
            #         g = sorted(zip(i_list, [self.parents_vec[i][j] for i in i_list]), key=lambda x: x[1])[0][0]
            #         if g not in self.most_depedent_child.keys():
            #             self.most_depedent_child[j] = [g]
            #         else:
            #             self.most_depedent_child[j].append(g)
            # self.most_depedent_child
            # def generate_tree(self):
            # def draw_dependency_with_abundance_with_parents_vec(self):
            # def most_common_strcutre(self, a_vec, motif_vec):
            #     NBT_motif_match_motifvec.find_common_structure_in_cluster(a_vec, )
