# break glycoCT
from glypy.algorithms.subtree_search import subtree_of
from glypy.structure.glycan import fragment_to_substructure
import seaborn as sns
from scipy.spatial import distance
sns.set(color_codes=True)
import time
import seaborn as sns
sns.set(color_codes=True)
from glypy.io import glycoct





# glyco_motif_list={}
# glycoct_list = []
aaa = ['WT',
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


def get_motif(glycoct_obj, idex=0):
    # print('start getmotif')
    _frag_motif_list = {}
    # fig, axes = plt.subplots(6,9)
    # fig.set_size_inches(14,6)
    start_time = time.time()
    for i in glycoct_obj.fragments(max_cleavages=len(glycoct_obj)):
        _frag_gly = fragment_to_substructure(i, glycoct_obj)

        # plot(_frag_gly)
        if not len(_frag_gly) in _frag_motif_list.keys():
            _frag_motif_list[len(_frag_gly)] = [glycoct.loads(_frag_gly)]
        else:
            _frag_motif_list[len(_frag_gly)].append(glycoct.loads(_frag_gly))
    mid_time = time.time()
    # print('start clean duplicate')
    _frag_motif_list = clean_duplicate(_frag_motif_list)
    end_time = time.time()
    print(idex, len(glycoct_obj), end_time - mid_time, mid_time - start_time)
    # print('finished getmotif')

    return _frag_motif_list


def clean_duplicate(_frag_motif_list):
    for i in _frag_motif_list.keys():
        # print(i)
        ldex = 0
        _check_list = _frag_motif_list[i]
        while ldex < len(_check_list):
            jdex = ldex + 1
            while jdex < len(_check_list):
                if subtree_of(_check_list[ldex], _check_list[jdex]) == 1 and subtree_of(_check_list[jdex],
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

class GlycanMotifLib:
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
    _man1 = glycoct.loads("""
        RES
        1b:b-dman-HEX-1:5
        LIN""")
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
        if type(motif_) == dict:
            print(type(list(motif_.keys())[0]))
            dict_keys = sorted([int(i) for i in motif_.keys()])
            self.motif_vec = []
            for i in dict_keys:
                for j in motif_[str(i)]:
                    if isinstance(j, type(self._man2)):
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
            if type(motif_[0])==type(self._man2):
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
        self.motif_with_core_dict = {}
        self.sia_ept_vec = []
        self.gala_ept_vec = []
        self.motif_with_core_list = []
        self.motif_list = [i for i in range(len(self.motif_vec))]
        self.motif_ncore_dep_tree = {}
        self.motif_dep_tree = {}
        self.motif_single_dep_tree = {}
        self.extract_motif_with_core()

    def get_motif_index(self):
        pass

    def create_epitope_vec(self):
        print("start motif with sia")
        if not self.sia_ept_vec:
            for i in sorted(list(self.motif_dict.keys())):
                # if (i) > 9:
                #     break
                # print("len", i)
                for j in self.motif_dict[i]:
                    """
                    motif j in i degree/motif in i-1 degree
                    """
                    if subtree_of(self._man1, self.motif_vec[j]) is not None or subtree_of(self._man2, self.motif_vec[j]) is not None:
                        continue
                    if subtree_of(self._with_sia_core, self.motif_vec[j]) is not None:
                        if len(self.motif_vec[j]) % 2 == 1:

                            self.sia_ept_vec.append(j)
                        # self.gala_ept_vec.append(j)
                    elif subtree_of(self._no_sia_core, self.motif_vec[j]) is not None:
                        if len(self.motif_vec[j]) % 2 == 0:
                            self.gala_ept_vec.append(j)

            print("Finish sia match ", len(self.sia_ept_vec),
                  " motifs are find with sia core ", len(self.gala_ept_vec),  " motifs are find with no sia core ")
        else:
            print("Finish sia match ", len(self.motif_with_core_list),
                  " motifs are find with sia core ", len(self.gala_ept_vec),  " motifs are find with no sia core ")

    def extract_motif_with_core(self):
        """ store the result in self.motif_with_core_list
        and return the count"""
        # count = []
        print("start motif_with core")
        if self.motif_with_core_dict == {}:
            for i in sorted(list(self.motif_dict.keys())):
                if len(self.nglycan_core) > i:
                    continue
                print("len", i)
                self.motif_with_core_dict[i] = []
                for j in self.motif_dict[i]:
                    """
                    motif j in i degree/motif in i-1 degree
                    """
                    if subtree_of(self.nglycan_core, self.motif_vec[j]) == 1:
                        self.motif_with_core_dict[i].append(j)
                        self.motif_with_core_list.append(j)
            print("Finish the n-glycan match ", len(self.motif_with_core_list),
                  " motifs are matched to the n-glycan core")
        else:
            print("Finish the n-glycan match ", len(self.motif_with_core_list),
                  " motifs are matched to the n-glycan core")

    def motif_with_ncore_dependence_tree(self):
        """ just connect motif to all parent"""
        print('start building ncore_dependence_tree')
        edge_list = []
        if self.motif_ncore_dep_tree == {}:
            for i in sorted(list(self.motif_with_core_dict.keys())):
                print(i)
                if len(self.nglycan_core) == i:
                    for j in self.motif_with_core_dict[i]:
                        self.motif_ncore_dep_tree[j] = []
                    continue
                for j in self.motif_with_core_dict[i]:
                    """
                    motif j in i degree/motif in i-1 degree
                    """
                    self.motif_ncore_dep_tree[j] = []
                    for k in self.motif_with_core_dict[i - 1]:
                        if subtree_of(self.motif_vec[k], self.motif_vec[j]) == 1:
                            self.motif_ncore_dep_tree[k].append(j)
                            edge_list.append((k, j))
        return self.motif_ncore_dep_tree, edge_list


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

    def motif_dependence_tree(self):
        """ connect motif to all parents"""
        print('start building dependence_tree')
        edge_list = []
        if self.motif_dep_tree == {}:
            # self.motif_dep_tree = {}
            id_list = []
            for i in sorted(list(self.motif_dict.keys())):
                print(i)
                if 1 == i:
                    for j in self.motif_dict[i]:
                        self.motif_dep_tree[j] = []
                    continue
                for j in self.motif_dict[i]:
                    """
                    motif j in i degree/motif in i-1 degree
                    """
                    self.motif_dep_tree[j] = []
                    for k in self.motif_dict[i - 1]:
                        if subtree_of(self.motif_vec[k], self.motif_vec[j]) == 1:
                            self.motif_dep_tree[k].append(j)
                            edge_list.append((k, j))
        return self.motif_dep_tree, edge_list

    def motif_single_dependence_tree(self):
        """ just connect motif to one parent"""
        edge_list = []
        if self.motif_single_dep_tree == {}:
            for i in sorted(list(self.motif_dict.keys())):
                if 1 == i:
                    for j in self.motif_dict[i]:
                        self.motif_single_dep_tree[j] = []
                    continue
                for j in self.motif_dict[i]:
                    """
                    motif j in i degree/motif in i-1 degree
                    """
                    self.motif_single_dep_tree[j] = []
                    for k in self.motif_dict[i - 1]:
                        if subtree_of(self.motif_vec[k], self.motif_vec[j]) == 1:
                            self.motif_single_dep_tree[k].append(j)
                            break
        return self.motif_single_dep_tree, edge_list


class MotifDpTree:
    def __init__(self, a_glycan_motif_lib, motif_weight):
        self.dep_tree = a_glycan_motif_lib.motif_with_ncore_dependence_tree()
        self.node_len = len(a_glycan_motif_lib.motif_with_core_list)
        self.all_nodes = a_glycan_motif_lib.motif_with_core_list
        self.parents_vec = {}
        for i in a_glycan_motif_lib.motif_with_core_list:
            self.parents_vec[i] = {}
        #### self.generate_tree()
        self.motif_weight = motif_weight
        self.normalized_motif_weight = {}
        self._normalized_weight()
        self.heavy_dependency = {}
        self.most_dependent_child = {}
        self.gala_ept_vec = a_glycan_motif_lib.gala_ept_vec[:]
        self.sia_ept_vec = a_glycan_motif_lib.sia_ept_vec[:]
        self.sia_gala_ept_vec = self.gala_ept_vec[:]
        self.sia_gala_ept_vec.extend(self.sia_ept_vec[:])


    def get_sia_gal_vec(self):
        rt_lst = []
        rt_lst.extend(self.gala_ept_vec[:])
        rt_lst.extend(self.sia_ept_vec[:])
        return rt_lst

    # drop_list = add_sia_gal(a_dp_tree)

    def get_drop_node_with_sia(self):
        rt_lst = self.drop_node(distance.correlation)
        rt_lst.extend(self.sia_gala_ept_vec)
        return rt_lst


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

    def _normalized_weight(self):
        for i in self.motif_weight.keys():
            _max = max(self.motif_weight[i])
            self.normalized_motif_weight[i] = [j/_max for j in self.motif_weight[i]]
            # _array[i,:] = _array[i,:]/_max

    def drop_node(self, method=distance.correlation):
        node_kept = []
        for i in self.parents_vec.keys():
            # print(i)
            for j in self.dep_tree[i]:
                self.parents_vec[j][i] = 1 - method(self.normalized_motif_weight[j], self.normalized_motif_weight[i])
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
                node_kept.append(i)
        return node_kept
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



