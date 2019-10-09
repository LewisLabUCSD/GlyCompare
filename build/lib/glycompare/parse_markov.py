
import re
import glypy

glypy2mk = {'GlcNAc': 'bNG',
            'adMan': 'aM',
            'bdMan': 'bM',
            'Fuc': 'aF',
            'Neu5Ac': 'aNN',
            'Gal': 'bA',}
glypy2pos = {'GlcNAc': 1,
          'adMan': 1,
          'bdMan': 1,
          'Fuc': 1,
          'Neu5Ac': 2,
          'Gal': 1}

mk2glypy = {'aF': 'Fuc',
 'aM': 'adMan',
 'aNN': 'Neu5Ac',
 'bA': 'Gal',
 'bM': 'bdMan',
 'bNG': 'GlcNAc'}

monosaccharides_dic = glypy.monosaccharides

def translate_mkov2glypy(temp_str, test_mode=False):
    mono_re = re.compile(r"(?P<position>[2-6])(?P<mono_mkv>[a-z][A-Z]+)")
    braket_degree=0

    temp_str = temp_str[::-1]
    if temp_str.find('nsA;)')==-1:
        temp_str = temp_str+'('
        temp_str = "nsA;)"+ temp_str[4:]
    glycan_list = []
    degree_list = []
    group_list = []
    while len(temp_str) != 0:
        if test_mode:
            print(temp_str[::-1])
        if temp_str[:5] == "nsA;)":
            temp_str=temp_str[5:]
        # elif temp_str[:4] == "nsA;":
        #     temp_str=temp_str[4:]
        elif temp_str[:2] == "NG":
            glycan_list.append((glypy.monosaccharides["GlcNac"], -1))
            degree_list.append(braket_degree)
            temp_str=temp_str[2:]
        elif temp_str[0] == ")":
            braket_degree+=1
    #         group_category = max(degree_list)+1
            temp_str=temp_str[1:]
        elif temp_str[0] =="(":

    #         group_category=
            temp_str=temp_str[1:]
            #合并同类项
            i=len(degree_list)-1
#             print(i)
#             print(len(degree_list), len(degree_list))
            while degree_list[i-1] == braket_degree and i>1:
                if test_mode:
                    print(glycan_list[i-1][0])
                    print(glycan_list[i][0])
    #             glycan_list[i-1][0]
    #             if glycan_list[i][0]
                glycan_list[i-1][0].add_monosaccharide(glycan_list[i][0],
                                                       position=glycan_list[i][1],
                                                       child_position=glycan_list[i][2])
                degree_list.pop()
                glycan_list.pop()
                i-=1
            glycan_list[i-1][0].add_monosaccharide(glycan_list[i][0],
                                                   position=glycan_list[i][1],
                                                   child_position=glycan_list[i][2])
            degree_list.pop()
            glycan_list.pop()
            braket_degree-=1
        elif mono_re.match(temp_str):
            _match = mono_re.match(temp_str)
            temp_str = temp_str[_match.end():]
            glycan_list.append((glypy.monosaccharides[mk2glypy[_match.group('mono_mkv')]],
                                int(_match.group('position')),
                                int(glypy2pos[mk2glypy[_match.group('mono_mkv')]])))
            degree_list.append(braket_degree)
        else:
            print("?", temp_str[:5])
            break
    # plot_glycan_utilities.plot_glycan(glypy.Glycan(root=glycan_list[0][0]))

    return glypy.Glycan(root=glycan_list[0][0])