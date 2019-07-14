import glypy
# from glypy.io import glycoct

# monosaccharides = glypy.monosaccharides
mono_dic = {'M1': "bdMan",
            'M2_1': "adMan",
            'M2_2': "adMan",
            'M3_1':"adMan",
            'M3_2':"adMan",
            "G": "GlcNAc",
            "Glc_1": "GlcNAc",
            "Glc_2": "GlcNAc",
            "F": "Fuc",
            "R": "Gal",
            "S": "Neu5Ac"
            }


class glycan_model():
#     Sia_ = glycoct.loads("""RES
# 1b:a-dgro-dgal-NON-2:6|1:a|2:keto|3:d
# 2s:n-acetyl
# LIN
# 1:1d(5+1)2n""")
    def __init__(self):
        """ a plain model"""
        # self.motif_vec = motif_vec
        # self.motif_weight_dict = motif_weight
        self.motif_weight_list = []
        self.normed_motif_weight_list = []
        self.total_weight = 0
        self.panel = None

        # self.panel = a_gly_dic
        # if self.panel is None:
        #     self.total_count = 0
        # elif isinstance(self.panel, glycan_mono):
        #     self.total_count = 1

    def _model_builder(self):
        pass

    def glycan_walk(self, a_gly_dic):
        """build glycan with 4 branches
            each node have count, """
        if self.panel is None:
            self.panel = a_gly_dic
            self.motif_weight_list = [a_gly_dic.count]
            self.total_weight = a_gly_dic.count
        else:
            assert a_gly_dic.name == self.panel.name, 'root not match'
            # _weight = self.motif_weight_dict[motif_index]
            self.panel.count += a_gly_dic.count
            for i in a_gly_dic.child:
                _found = False
                for j in self.panel.child:
                    if j.name == i.name:
                        # print(i.name, 'add',i.count)
                        j.count += i.count
                        _found = True
                        self._recur_glycan_walk(i, j)
                if not _found:
                    self.panel.child.append(i)
            self.motif_weight_list.append(a_gly_dic.count)
            self.total_weight += a_gly_dic.count

        # self.total_count += 1

    def _recur_glycan_walk(self, i_gly_dic, panel_dic):
        assert i_gly_dic.name == panel_dic.name, 'recur root not match'
        for i in i_gly_dic.child:
            _found = False
            for j in panel_dic.child:
                if j.name == i.name:
                    # print(i.name, 'add',i.count)
                    j.count += i.count
                    _found = True
                    self._recur_glycan_walk(i, j)
            if not _found:
                panel_dic.child.append(i)

    def show_statistics(self):
        pass

    def plot_statistics(self):
        pass

    def travel_str_dict(self):
        return self._travel_str_dict(self.panel)

    def _travel_str_dict(self, a_mono):
        return a_mono.name + ',' + str(a_mono.count / self.total_weight) + '[' + ','.join(
            [self._travel_str_dict(i) for i in a_mono.child]) + ']'

    def get_common_representative(self, threshold):
        """ go through glycans
        G,1.0[G,1.0[M1,1.0[M2_1,1.0[Glc_1,1.0[R,0.9[S,0.6[]]],Glc_2,1.0[R,0.45[S,0.25[]]]],
        M2_2,1.0[Glc_1,1.0[R,0.5[S,0.2[]]],Glc_2,0.95[R,0.15[]]]]],F,0.5[]]'"""
        root_mono = glypy.monosaccharides["GlcNAc"]
        for i in self.panel.child:
            # print(i.name, i.count)
            if i.count / self.total_weight < threshold: continue
            root_mono.add_monosaccharide(self._helper_get_common(i, threshold))
        return glypy.Glycan(root=root_mono)

    def get_reps(self, threshold):
        return_list = []
        # for _thre in threshold_list:
        root_mono = glypy.monosaccharides["GlcNAc"]
        for i in self.panel.child:
            # print(i.name, i.count)
            if i.count / self.total_weight < threshold: continue
            root_mono.add_monosaccharide(self._helper_get_common(i, threshold))
        return_list.append(glypy.Glycan(root=root_mono))
        return return_list

    def _helper_get_common(self, a_glycan_mono, threshold):
        # print(mono_dic[a_glycan_mono.name])
        # print(glypy.monosaccharides(mono_dic[a_glycan_mono.name]))
        _mono = glypy.monosaccharides[mono_dic[a_glycan_mono.name]]
        for i in a_glycan_mono.child:
            # print(i.name, i.count, self.total_weight, i.count / self.total_weight, threshold)
            if i.count / self.total_weight < threshold:
                continue
            _mono.add_monosaccharide(self._helper_get_common(i, threshold))
        return _mono


# mano = ['M1', 'M2']
# glunac = ["G", "Glc_1", "Glc_2"]
# fuc = ["F"]
# gal = ["R"]
# sia = ["S"]


class glycan_mono():
    def __init__(self, name, child=[], count=1):
        self.child = child
        self.name = name
        self.count = count
        self.std_count = 0


        # def add_child(self, a_glycan_mono):
        #     self.child.append(a_glycan_mono)


# def glycan_mono_copy():
#     """only copy child"""


# _mono=_gly.root.children()[0][1]#.children()#[0][1].children()
ROOT = 1
ROOTG = 2
FU = 3
M1 = 4
M2 = 5
M3 = 7
PASM = 6


def _travel_pasm(mono, weight):
    if mono.children():
        _, _child = mono.children()[0]
        _branch = []
        if str(_child).find('n-acetyl') != -1 and str(_child).find('glc') != -1:
            # _branch['G'] = _travel_pasm(_child)
            _temp_mono = glycan_mono('G', _travel_pasm(_child, weight=weight), count=weight)
            _branch.append(_temp_mono)
        elif str(_child).find('dgro-dgal') != -1:
            # _branch['S'] = _travel_pasm(_child)
            _temp_mono = glycan_mono('S', _travel_pasm(_child, weight=weight), count=weight)
            _branch.append(_temp_mono)
        elif str(_child).find('gal') != -1:
            # _branch['R'] = _travel_pasm(_child)
            _temp_mono = glycan_mono('R', _travel_pasm(_child, weight=weight), count=weight)
            _branch.append(_temp_mono)
        else:
            assert False, 'Error in _travel_pasm'

        return _branch
    else:
        return []


def travel_str_dict(a_mono):
    return a_mono.name + ',' + str(a_mono.count) + '[' + ','.join([travel_str_dict(i) for i in a_mono.child]) + ']'


def traves_glycan(a_gly, weight):
    _root = a_gly.root
    _state = ROOT
    return glycan_mono('G', _re_travel_glycan(_root, _state, weight), count=weight)


def get_depth_man(a_dict):
    """get the depth of a glycan"""
    # if len(a_dict.keys()) ==
    if a_dict.child:
        return get_depth_gal(a_dict.child[0])
    else:
        return 0


def get_depth_gal(a_dict):
    if a_dict.child:
        # list(a_dict.keys())[0]
        return 1 + get_depth_gal(a_dict.child[0])
    else:
        return 1


def map_glycan(a_g_dict, b_g_dict):
    """combine two glycan together, if one node doesn't have branch, copy the other"""
    pass
    # for i in a_g_dict:


def _re_travel_glycan(mono, state, weight):
    """return a list of glycan_mono as child"""
    if mono.children():
        if state == M2:
            _count = 0
            _man_count = 0
            _branch = []
            for _, _child in mono.children():
                _count += 1
                if str(_child).find('glc') != -1:
                    # _branch['Glc_' + str(_count)] = _travel_pasm(_child)  # which is man2 passed through
                    _temp_mono = glycan_mono('Glc_' + str(_count), _travel_pasm(_child, weight), count=weight)
                    _branch.append(_temp_mono)
                elif str(_child).find('man') != -1:

                    _man_count += 1
                    # _branch['M2_' + str(_count)] = _re_travel_glycan(_child, state=M2)  # which is man2 passed through
                    _temp_mono = glycan_mono('M3_' + str(_man_count), count=weight)
                    _branch.append(_temp_mono)

            if len(mono.children()) == 2:
                d1 = get_depth_gal(_branch[0])
                d2 = get_depth_gal(_branch[1])
                if d1 > d2:
                    pass
                else:
                    d1_dict = _branch[0].child[:]
                    _branch[0].child = _branch[1].child[:]
                    _branch[1].child = d1_dict[:]
            return _branch

        elif state == ROOT:
            if len(mono.children()) == 2:
                _branch = []
                for _, _child in mono.children():
                    if str(_child).find('a-lgal') != -1:
                        # _branch['F'] = {}
                        _temp_mono = glycan_mono('F', count=weight)
                        _branch.append(_temp_mono)
                    elif str(_child).find('n-acetyl') != -1 and str(_child).find('dglc') != -1:
                        #                         _branch['G'] = {}
                        # _branch['G'] = _re_travel_glycan(_child, state=ROOTG)
                        _temp_mono = glycan_mono('G', _re_travel_glycan(_child, state=ROOTG, weight=weight), count=weight)
                        _branch.append(_temp_mono)
            elif len(mono.children()) == 1:
                _branch = []
                _, _child = mono.children()[0]
                _temp_mono = glycan_mono('G', _re_travel_glycan(_child, state=ROOTG, weight=weight), count=weight)
                _branch.append(_temp_mono)
            else:
                assert False, 'Wrong in root'
            return _branch

        elif state == ROOTG:
            assert len(mono.children()) == 1, 'Wrong length in RGal'
            _branch = []
            _, _child = mono.children()[0]
            assert str(_child).find('man') != -1, 'Wrong man in RGal'
            # _branch['M1'] = _re_travel_glycan(_child, state=M1)
            _temp_mono = glycan_mono('M1', _re_travel_glycan(_child, state=M1, weight=weight), count=weight)
            _branch.append(_temp_mono)
            return _branch

        elif state == M1:
            # print('m1')
            _count = 0
            _branch = []
            for _, _child in mono.children():
                _count += 1
                if str(_child).find('man') != -1:
                    # _branch['M2_' + str(_count)] = _re_travel_glycan(_child, state=M2)  # which is man2 passed through
                    _temp_mono = glycan_mono('M2_' + str(_count), _re_travel_glycan(_child, state=M2, weight=weight), count=weight)
                    _branch.append(_temp_mono)
                else:
                    assert False, 'find other mono in M1'
            if len(mono.children()) == 2:
                d1 = get_depth_man(_branch[0])
                # print('two branches')
                d2 = get_depth_man(_branch[1])
                # print(d1, d2)
                if d1 > d2:
                    pass
                else:
                    d1_dict = _branch[0].child[:]
                    _branch[0].child = _branch[1].child[:]
                    _branch[1].child = d1_dict[:]
            return _branch

    else:
        return []
