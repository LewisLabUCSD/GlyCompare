from glypy.structure.glycan_composition import GlycanComposition, FrozenGlycanComposition
from glypy.algorithms.subtree_search.inclusion import subtree_of
from glypy.io import glycoct
import __init__

nglycan_mgat1 = """
RES
1b:x-dglc-HEX-1:5
2s:n-acetyl
3b:b-dglc-HEX-1:5
4s:n-acetyl
5b:b-dman-HEX-1:5
6b:a-dman-HEX-1:5
7b:b-dglc-HEX-1:5
8s:n-acetyl
9b:a-dman-HEX-1:5
LIN
1:1d(2+1)2n
2:1o(4+1)3d
3:3d(2+1)4n
4:3o(4+1)5d
5:5o(3+1)6d
6:6o(2+1)7d
7:7d(2+1)8n
8:5o(6+1)9d
"""

tri_branch_core = """
RES
1b:x-dglc-HEX-1:5
2s:n-acetyl
3b:b-dglc-HEX-1:5
4s:n-acetyl
5b:b-dman-HEX-1:5
6b:a-dman-HEX-1:5
7b:b-dglc-HEX-1:5
8s:n-acetyl
9b:b-dglc-HEX-1:5
10s:n-acetyl
11b:a-dman-HEX-1:5
12b:b-dglc-HEX-1:5
13s:n-acetyl
14b:a-lgal-HEX-1:5|6:d
LIN
1:1d(2+1)2n
2:1o(4+1)3d
3:3d(2+1)4n
4:3o(4+1)5d
5:5o(3+1)6d
6:6o(2+1)7d
7:7d(2+1)8n
8:6o(4+1)9d
9:9d(2+1)10n
10:5o(6+1)11d
11:11o(2+1)12d
12:12d(2+1)13n
13:1o(6+1)14d
"""

nglycan_core = """
RES
1b:x-dglc-HEX-1:5
2s:n-acetyl
3b:b-dglc-HEX-1:5
4s:n-acetyl
5b:b-dman-HEX-1:5
6b:a-dman-HEX-1:5
7b:a-dman-HEX-1:5
LIN
1:1d(2+1)2n
2:1o(4+1)3d
3:3d(2+1)4n
4:3o(4+1)5d
5:5o(3+1)6d
6:5o(6+1)7d """

nglycan_m5 = """
RES
1b:x-dglc-HEX-1:5
2s:n-acetyl
3b:b-dglc-HEX-1:5
4s:n-acetyl
5b:b-dman-HEX-1:5
6b:a-dman-HEX-1:5
7b:a-dman-HEX-1:5
8b:a-dman-HEX-1:5
9b:a-dman-HEX-1:5
LIN
1:1d(2+1)2n
2:1o(4+1)3d
3:3d(2+1)4n
4:3o(4+1)5d
5:5o(3+1)6d
6:5o(6+1)7d
7:7o(3+1)8d
8:7o(6+1)9d"""

abbrev_dict = {'Glc2NAc': 'A',
               'Man': 'M',
               'Fuc': 'F',
               'Neu5Ac': 'S',
               'Gal': 'G',
               'Glu': 'g'}


class nglycan_composition():
    def __init__(self, a_glycan):
        self.shorthand_dict = {"A": 0,
                               "M": 0,
                               "F": 0,
                               "S": 0,
                               "G": 0,
                               "g": 0}

        if type(a_glycan) == str:
            self._glycan = glycoct.loads(a_glycan)
        else:
            self._glycan = a_glycan

        a_composition = GlycanComposition.from_glycan(a_glycan)
        items_dict = {}
        for i, j in list(a_composition.items()):
            items_dict[str(i)] = j
        # print(items_dict)
        for i in abbrev_dict:
            try:
                self.shorthand_dict[abbrev_dict[i]] = items_dict[i]
            except:
                pass
        self.possible = self._is_possible()

    def _is_possible(self):
        """ Check if the glycan structure exists biologically"""
        if self.shorthand_dict['A'] >= 3:
            nglycan_3GlcNAc = glycoct.loads(nglycan_mgat1)
            if subtree_of(nglycan_3GlcNAc, self._glycan, exact=__init__.exact_Ture) != 1:
                return False

        if self.shorthand_dict['M'] >= 5:
            nglycan_5Man = glycoct.loads(nglycan_m5)
            if subtree_of(nglycan_5Man, self._glycan, exact=__init__.exact_Ture) != 1:
                return False
        elif self.shorthand_dict['M'] >= 3:
            nglycan_3Man = glycoct.loads(nglycan_core)
            if subtree_of(nglycan_3Man, self._glycan, exact=__init__.exact_Ture) != 1:
                return False

        # should be disabled for b4Galt
        if self.shorthand_dict['G'] != self.shorthand_dict['A'] - 2:
            return False

        # nglycan_5GlcNAc = glycoct.loads(tri_branch_core)
        # if subtree_of(nglycan_5GlcNAc, self._glycan, exact=False) == 1:
        #     if subtree_of(nglycan_5GlcNAc, self._glycan, exact=True) != 1:
        #         return False

        return True

    def shorthand(self):
        if self.shorthand_dict['M'] > 4 and self.shorthand_dict['g']:
            return "M" + str(self.shorthand_dict['M'] + self.shorthand_dict['g'])
        if self.shorthand_dict['A'] - 2 < 0:
            assert False, 'GlcNAc is greater than 2'
        elif self.shorthand_dict['A'] - 2 == 0:
            return_str = ''
        else:
            return_str = 'A' + str(self.shorthand_dict['A'] - 2)

        if self.shorthand_dict['F'] == 0:
            pass
        elif self.shorthand_dict['F'] == 1:
            return_str += "F"
        else:
            return_str += "F" + str(self.shorthand_dict['F'])

        if self.shorthand_dict['G'] == 0:
            pass
        elif self.shorthand_dict['G'] == 1:
            return_str += "G"
        else:
            return_str += "G" + str(self.shorthand_dict['G'])

        if self.shorthand_dict['S'] == 0:
            pass
        elif self.shorthand_dict['S'] == 1:
            return_str += "S"
        else:
            return_str += "S" + str(self.shorthand_dict['S'])
        return return_str

    def mono_composition(self):
        HexNAc = self.shorthand_dict['A']
        dHex = self.shorthand_dict['F']
        Hex = self.shorthand_dict['M'] + self.shorthand_dict['G'] + self.shorthand_dict['g']
        NeuAc = self.shorthand_dict['S']
        return {'HexNAc': HexNAc,
                'dHex': dHex,
                'Hex': Hex,
                'NeuAc': NeuAc}

    def composition(self):
        _composition = self._glycan.total_composition()
        return_str = ''
        for i in sorted(_composition.keys()):
            return_str += i + str(_composition[i])
        return return_str
