import pandas as pd
import numpy as np
from matplotlib import pyplot as plt

from json_utility import *
import plot_glycan_utilities


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


#
# def load_glycan_profile_dic():
#     glycan_profile = {'1': {'2244': 'G04483SK',
#                             '2605': 'G30460NZ',
#                             '2967': 'G17689DH',
#                             '3416': 'G54338PJ',
#                             '3865': '3865.1',
#                             '4226': 'G86696LV',
#                             '4587': '4587.1',
#                             '5037': 'G49604DB',
#                             '5486': '5486.1'}}
#     glycan_profile_merged = {'1': {'2244': 'G04483SK',
#                                    '2326': 'G00176HZ',
#                                    '2530': 'G79457WN',
#                                    '2605': 'G30460NZ',
#                                    '2734': 'G00536FZ',
#                                    '2939': 'G37597FW',
#                                    '2967': 'G17689DH',
#                                    '3416': 'G54338PJ',
#                                    '3865': '3865.1',
#                                    '4226': 'G86696LV',
#                                    '4587': '4587.1',
#                                    '5037': 'G49604DB',
#                                    '5486': '5486.1',
#                                    '3055': '3055.1',
#                                    '3777': 'G76812VG',
#                                    '3504': 'G10292TC'
#                                    }}
#     glycan_profile_merged['2'] = {
#         '2244': 'G04483SK',
#         '2326': 'G00176HZ',
#         '2489': 'G10691MJ',
#         '2530': 'G79457WN',
#         '2605': 'G30460NZ',
#         '2734': 'G00536FZ',
#         '2939': 'G37597FW',
#         '2967': 'G17689DH',
#         '3416': 'G54338PJ',
#         '3661': '3661.1',
#         '3865': '3865.1',
#         '4226': 'G86696LV',
#         '4587': '4587.1',
#         '5037': 'G49604DB',
#         '3143': 'G20924UR',
#         '3055': '3055.1',
#         '2693': 'G07568IR',
#         '3777': 'G76812VG',
#         '3504': 'G76812VG',
#     }
#     glycan_profile['2'] = {
#         "2244": "G04483SK",
#         "2605": "G30460NZ",
#         "2967": "G17689DH",
#         "3416": "G54338PJ",
#         "3865": "3865.1",
#         "4226": "G86696LV",
#         "4587": "4587.1",
#         "5037": "G49604DB",
#     }
#     # mgat4A/mgat4B
#     glycan_profile['3'] = {
#         "2244": "G04483SK",
#         "2693": "G07568IR",
#         "3055": "G88127MB",
#         "3416": "3416.3",
#         "3777": "G76812VG",
#         "4226": "G56516KW",
#         "4675": "4675.1",
#     }
#     glycan_profile_merged['3'] = {
#         '2244': 'G04483SK',
#         '2693': 'G07568IR',
#         '3055': 'G88127MB',
#         '3416': '3416.3',
#         '3777': 'G76812VG',
#         '4226': 'G56516KW',
#         '4675': '4675.1',
#         '2605': 'G30460NZ',
#         '2967': 'G17689DH',
#         '3865': '3865.1'
#     }
#     # mgat5
#     glycan_profile['4'] = {
#         "2244": "G04483SK",
#         "2605": "G30460NZ",
#         "2967": "G17689DH",
#         "3416": "G54338PJ",
#         "3777": "G76812VG",
#         "4226": "G56516KW",
#     }
#     glycan_profile_merged['4'] = {
#         "2244": "G04483SK",
#         "2605": "G30460NZ",
#         "2967": "G17689DH",
#         "3416": "G54338PJ",
#         "3777": "G76812VG",
#         "4226": "G56516KW",
#         '3055': '3055.1',
#         '2693': 'G07568IR',
#         '1836': 'G80858MF'
#     }
#     glycan_profile['5'] = {
#         "2244": "G04483SK",
#         "2605": "G30460NZ",
#         "2967": "G17689DH",
#         "3416": "3416.1",
#     }
#     glycan_profile_merged['5'] = {'2040': 'G58667NI',
#                                   '2244': 'G04483SK',
#                                   '2605': 'G30460NZ',
#                                   '2967': 'G17689DH',
#                                   '3416': '3416.1'}
#
#     glycan_profile['6'] = {
#         "1591": "G07483YN",
#         "1836": "G80858MF",
#         "2081": "G80393PG",
#         "2326": "G00176HZ",
#         "2530": "G79457WN",
#         "2892": "G79412GP",
#         "3096": "G40242TG",  # ?
#         "3457": "3457.1",
#         "3777": "G76812VG",
#         "4022": "G80264ZA",
#         "4226": "G86696LV",
#         "4587": "4587.1",
#         "5037": "G49604DB",
#     }
#     glycan_profile_merged['6'] = {'1591': 'G07483YN',
#                                   '1836': 'G80858MF',
#                                   '2081': 'G80393PG',
#                                   '2326': 'G00176HZ',
#                                   '2489': 'G10691MJ',
#                                   '2530': 'G79457WN',
#                                   '2734': 'G00536FZ',
#                                   '2892': 'G79412GP',
#                                   '3096': 'G40242TG',
#                                   '3457': '3457.1',
#                                   '3777': 'G76812VG',
#                                   '4022': 'G80264ZA',
#                                   '4226': 'G86696LV',
#                                   '4587': '4587.1',
#                                   '5037': 'G49604DB',
#                                   '2040': 'G58667NI',
#                                   '2285': 'G03445UI',
#                                   '2244': 'G04483SK',
#                                   '2401': 'G23295TF',
#                                   '2605': 'G30460NZ',
#                                   '2967': 'G17689DH',
#                                   '3212': '3212.1',
#                                   '3661': '3661.1'}
#     # B4GalT2
#     glycan_profile['7'] = {
#         "2244": "G04483SK",
#         "2605": "G30460NZ",
#         "2967": "G17689DH",
#         "3416": "G54338PJ",
#         "3777": "G76812VG",
#         "3865": "3865.1",
#         "4226": "G86696LV",
#         "4587": "4587.1",
#         "4675": "G09280JF",
#         "5037": "G49604DB",
#         "5486": "5486.1",
#     }
#     # glycan_profile_merged['7'] = {'2081': 'G80393PG',
#     #                                   '2244': 'G04483SK',
#     #                                   '2605': 'G30460NZ',
#     #                                   '2967': 'G17689DH',
#     #                                   '3416': 'G54338PJ',
#     #                                   '3504': 'G10292TC',
#     #                                   '3777': 'G76812VG',
#     #                                   '3865': '3865.1',
#     #                                   '4226': 'G86696LV',
#     #                                   '4314': 'G16873YG',
#     #                                   '4587': '4587.1',
#     #                                   '4675': 'G09280JF',
#     #                                   '5037': 'G49604DB',
#     #                                   '5486': '5486.1',
#     #                                   '3055': '3055.1'}
#     # B4GalT3
#     glycan_profile['8'] = {
#         "2605": "G30460NZ",
#         "2967": "G17689DH",
#         "3416": "G54338PJ",
#         "3777": "G76812VG",
#         "3865": "3865.1",
#         "4226": "G86696LV",
#         "4587": "4587.1",
#         "4675": "G09280JF",
#         "5037": "G49604DB",
#         "5486": "5486.1",
#     }
#     glycan_profile_merged['8'] = {'2605': 'G30460NZ',
#                                   '2967': 'G17689DH',
#                                   "3055": '3055.1',
#                                   '3416': 'G54338PJ',
#                                   '3504': 'G10292TC',
#                                   '3661': '3661.1',
#                                   '3777': 'G76812VG',
#                                   '3865': '3865.1',
#                                   '3953': '3953.1',
#                                   '4226': 'G86696LV',
#                                   "4314": 'G16873YG',
#                                   '4587': '4587.1',
#                                   '4675': 'G09280JF',
#                                   '5037': 'G49604DB',
#                                   '5486': '5486.1'}
#     # B4GalT4
#     glycan_profile['9'] = {
#         "2967": "G17689DH",
#         "3416": "G54338PJ",
#         "3777": "G76812VG",
#         "4226": "G86696LV",
#         "4587": "4587.1",
#         "5037": "G49604DB",
#         "5486": "5486.1",
#     }
#     glycan_profile_merged['9'] = {'2967': 'G17689DH',
#                                   '3416': 'G54338PJ',
#                                   '3777': 'G76812VG',
#                                   "3865": '3865.1',
#                                   '4022': 'G80264ZA',
#                                   '4226': 'G86696LV',
#                                   '4587': '4587.1',
#                                   '5037': 'G49604DB',
#                                   '5486': '5486.1'}
#     # B4GalT1/B4GalT2
#     glycan_profile['10'] = {
#         "1591": "G07483YN",
#         "1836": "G80858MF",
#         "2081": "G80393PG",
#         "2326": "G00176HZ",
#         "2530": "G79457WN",
#         "2734": "G00536FZ",
#         "2892": "G79412GP",
#         "3096": "G40242TG",
#         "3457": "3457.1",
#         "3661": "3661.1",
#         "3865": "3865.1",
#         "4022": "G80264ZA",
#         "4226": "G86696LV",
#         "4587": "4587.1",
#         "5037": "G49604DB",
#     }
#     glycan_profile_merged['10'] = {'1591': 'G07483YN',
#                                    '1836': 'G80858MF',
#                                    "2040": 'G58667NI',
#                                    '2081': 'G80393PG',
#                                    '2285': 'G03445UI',
#                                    '2326': 'G00176HZ',
#                                    '2530': 'G79457WN',
#                                    '2646': '2646.1',
#                                    '2734': 'G00536FZ',
#                                    '2892': 'G79412GP',
#                                    '3096': 'G40242TG',
#                                    '3212': '3212.1',
#                                    '3457': '3457.1',
#                                    '3661': '3661.1',
#                                    '3865': '3865.1',
#                                    '4022': 'G80264ZA',
#                                    '4226': 'G86696LV',
#                                    '4587': '4587.1',
#                                    '5037': 'G49604DB'}
#     # B4GalT1/B4GalT2
#     glycan_profile['11'] = {
#         "1591": "G07483YN",
#         "1836": "G80858MF",
#         "2081": "G80393PG",
#         "2326": "G00176HZ",
#         "2646": "2646.1",
#         "2892": "G79412GP",
#         "3212": "3212.1",
#         "3457": "3457.1",
#         "4022": "G80264ZA",
#     }
#     glycan_profile_merged['11'] = {'1591': 'G07483YN',
#                                    '1836': 'G80858MF',
#                                    '2040': 'G58667NI',
#                                    '2081': 'G80393PG',
#                                    '2285': 'G03445UI',
#                                    '2326': 'G00176HZ',
#                                    '2530': 'G79457WN',
#                                    '2646': '2646.1',
#                                    '2892': 'G79412GP',
#                                    '3212': '3212.1',
#                                    '3457': '3457.1',
#                                    '4022': 'G80264ZA'}
#     # B3gnt1
#     glycan_profile['12'] = {'1836': 'G80858MF',
#                             '2040': 'G58667NI',
#                             '2081': 'G80393PG',
#                             '2244': 'G04483SK',
#                             '2285': 'G03445UI',
#                             '2489': 'G10691MJ',
#                             '2605': 'G30460NZ',
#                             '2693': 'G07568IR',
#                             '2967': 'G17689DH',
#                             '3055': '3055.1',
#                             '3416': 'G54338PJ',
#                             '3504': 'G10292TC',
#                             '3777': 'G76812VG',
#                             '3865': '3865.1',
#                             '3953': '3953.1',
#                             '4226': 'G86696LV',
#                             '4314': 'G16873YG',
#                             '4587': '4587.1',
#                             '5037': 'G49604DB'}
#     # B3gnt2
#     glycan_profile['13'] = {
#         "2605": "G30460NZ",
#         "3055": "3055.1",
#         "3416": "G54338PJ",
#         "3777": "G76812VG",
#         "3865": "3865.1",
#         "4226": "G86696LV",
#         "4587": "4587.1",
#     }
#     # st3gal3
#     glycan_profile['14'] = {
#         "2244": "G04483SK",
#         "2605": "G30460NZ",
#         "3055": "3055.1",
#         "3416": "G54338PJ",
#         "3865": "3865.1",
#         "4226": "G86696LV",
#         "4587": "4587.1",
#         "5037": "G49604DB",
#     }
#     # st3gal3
#     glycan_profile['14'] = {
#         "2244": "G04483SK",
#         "2605": "G30460NZ",
#         "3055": "3055.1",
#         "3416": "G54338PJ",
#         "3865": "3865.1",
#         "4226": "G86696LV",
#         "4587": "4587.1",
#         "5037": "G49604DB",
#     }
#     # st3gal4
#     glycan_profile['15'] = {
#         "2244": "G04483SK",
#         "2489": "G10691MJ",
#         "2693": "G07568IR",
#         "2939": "G37597FW",
#         "3143": "G20924UR",
#         "3504": "G10292TC",
#         "3592": "G75308SV",
#         "3953": "3953.1",
#         "4402": "4402.1",
#     }
#     # st3gal6
#     glycan_profile['16'] = {
#         "1591": "G07483YN",
#         "1836": "G80858MF",
#         "2040": "G58667NI",
#         "2081": "G80393PG",
#         "2244": "G04483SK",
#         "2285": "G03445UI",
#         "2489": "G10691MJ",
#         "2693": "G07568IR",
#         "3055": "3055.1",
#         "3143": "G20924UR",
#         "3416": "G54338PJ",
#         "3504": "G10292TC",
#         "3592": "G75308SV",
#         "3865": "3865.1",
#         "4226": "G86696LV",
#         "4314": "G16873YG",
#         "4675": "G09280JF",
#         "5037": "G49604DB",
#     }
#     # st3gal3/4
#     glycan_profile['17'] = {
#         "2040": "G58667NI",
#         "2244": "G04483SK",
#         "2489": "G10691MJ",
#         "2693": "G07568IR",
#         "2939": "G37597FW",
#         "3143": "G20924UR",
#         "3504": "G10292TC",
#         "3592": "G75308SV",
#         "3953": "3953.1",
#         "4041": "4041.1",
#         "4402": "4402.1",
#         "4490": "4490.1",
#         "4851": "4851.1",
#     }
#     # st3gal4/6
#     glycan_profile['18'] = {
#         "2040": "G58667NI",
#         "2244": "G04483SK",
#         "2489": "G10691MJ",
#         "2693": "G07568IR",
#         "2939": "G37597FW",
#         "3143": "G20924UR",
#         "3388": "G88966ZO",
#         "3592": "G75308SV",
#         "4041": "4041.1",
#     }
#     # st3gal4/6
#     glycan_profile['19'] = {
#         "1836": "G80858MF",
#         "2081": "G80393PG",
#         "2605": "G30460NZ",
#         "2967": "G17689DH",
#         "3055": "3055.1",
#         "3416": "G54338PJ",
#         "3504": "G10292TC",
#         "3865": "3865.1",
#         "4226": "G86696LV",
#         "4675": "G09280JF",
#         "5037": "G49604DB",
#         "5486": "5486.1",
#     }
#     # B3gnt2/mgat4a/mgat4b/mgat5
#     glycan_profile['20'] = {
#         "1836": "G80858MF",
#         "2244": "G04483SK",
#         "2605": "G30460NZ",
#         "2967": "G17689DH",
#     }
#     # st3gal4/6/mgat4A/4B/5
#     glycan_profile['21'] = {
#         "1836": "G80858MF",
#         "2040": "G58667NI",
#         "2244": "G04483SK",
#         "2605": "G30460NZ",
#         "2693": "2693.2",
#         "3143": "G54953LX",
#         "3592": "3592.1",
#         "4041": "G90130AG",
#     }
#     # st3gal4/6/mgat4A/4B/5
#     glycan_profile['22'] = {
#         "1836": "G80858MF",
#         "2401": "G23295TF",
#         "2605": "G30460NZ",
#         "2967": "G17689DH",
#     }
#     # mgat3
#     glycan_profile['23'] = {
#         "2244": "G04483SK",
#         "2605": "G30460NZ",
#         "2967": "G17689DH",
#         "3416": "G54338PJ",
#         "3777": "G76812VG",
#         "4226": "G86696LV",
#         "4587": "4587.1",
#         "5037": "G49604DB",
#         "5486": "5486.1",
#     }
#     # mgat3
#     glycan_profile['24'] = {
#         "2605": "G30460NZ",
#         "2967": "G17689DH",
#         "3416": "G54338PJ",
#         "3777": "G76812VG",
#         "4226": "G86696LV",
#         "4587": "4587.1",
#         "5037": "G49604DB",
#         "5486": "5486.1",
#     }
#     # mgat2
#     glycan_profile['25'] = {
#         "1591": "G07483YN",
#         "1795": "G12398HZ",
#         "2156": "G39439UR",
#         "2401": "2401.1",
#         "2605": "2605.1",
#         "2967": "2967.1",
#         "3416": "3416.2",
#     }
#     # B4galt1/B4galt2/B4galt3
#     glycan_profile['26'] = {
#         "1836": "G80858MF",
#         "2081": "G80393PG",
#         "2326": "G00176HZ",
#         "2646": "2646.1",
#         "2892": "G79412GP",
#     }
#     # B3gnt8
#     glycan_profile['27'] = {
#         "2244": "G04483SK",
#         "2605": "G30460NZ",
#         "2967": "G17689DH",
#         "3416": "G54338PJ",
#         "3777": "G76812VG",
#         "4226": "G86696LV",
#         "4587": "4587.1",
#         "5037": "G49604DB",
#         "5486": "5486.1",
#     }
#     # mgat4B
#     glycan_profile['28'] = {
#         "2244": "G04483SK",
#         "2605": "G30460NZ",
#         "2967": "G17689DH",
#         "3055": "3055.1",
#         "3416": "G54338PJ",
#         "3777": "G76812VG",
#         "4226": "G56516KW",
#         "4675": "4675.1",
#     }
#     # mgat5B
#     glycan_profile['29'] = {
#         "2605": "G30460NZ",
#         "2967": "G17689DH",
#         "3416": "G54338PJ",
#         "3777": "G76812VG",
#         "4226": "G86696LV",
#         "4587": "4587.1",
#         "5037": "G49604DB",
#     }
#     # mgat1
#     glycan_profile['30'] = {
#         "1375": "G60415BS",
#         "1580": "G49721VX",
#         "1754": "1754.1",
#     }
#     # mgat2/st3gal4/st3gal6
#     glycan_profile['31'] = {
#         "1591": "G07483YN",
#         "1795": "G12398HZ",
#         "2156": "G39439UR",
#         "2244": "G52428MJ",
#         "2605": "2605.1",
#         "2693": "2693.1",
#         "2967": "2967.1",
#     }
#     # mgat2/mgat4A/mgat4B/mgat5
#     glycan_profile['32'] = {
#         "1795": "G12398HZ",
#         "2156": "G39439UR",
#         "2605": "2605.2",
#     }
#     # mgat2/mgat4A/mgat4B/mgat5
#     glycan_profile['33'] = {
#         "1795": "G12398HZ",
#         "2244": "2244.1",
#         "2605": "2605.2",
#         "2693": "2693.3",
#         "3055": "3055.2",
#         "3143": "3143.1",
#     }
#     # mgat2/mgat4A/mgat4B/mgat5
#     glycan_profile['34'] = {
#         "2792": "G39764AC",
#         "3242": "G39813YP",
#         "3603": "G05098FE",
#         "4052": "G99891PR",
#         "4413": "G05203UQ",
#         "4862": "G85809SI",
#         "5312": "5312.1",
#     }
#     return glycan_profile
#


# def load_mz_glycan_dict(addr=__init__.json_address + "NBT_mz_dict_glycan_glycoct.json"):
#     """
#
#     :param addr:
#     :return:
#     """
#     return load_json(addr)

#
# def load_cho_mz_abundance(cho_addr=__init__.source_address + 'nbt.3280_cho.txt',
#                           mz_abd_addr=__init__.source_address + 'glycan_table.xls'):
#     """
#
#     :param cho_addr:
#     :param mz_abd_addr:
#     :return:
#     """
#     mz_glycan_dict = load_mz_glycan_dict()
#     glycan_profile_table = pd.read_table(cho_addr)
#     glycan_profile_table = glycan_profile_table.fillna(0)
#     mz_list = [int(i) for i in mz_glycan_dict.keys()]
#     mz_abd_table = glycan_profile_table[glycan_profile_table['m/z exp'].isin(mz_list)]
#
#     del mz_abd_table['Composition']
#     del mz_abd_table['m/z calc']
#     _col = mz_abd_table.columns.values
#     # _col =
#     # _col[0] = 'm/z'
#     # _col = list(range(len(_col)))
#     _col[0] = 'm/z'
#     # print(_col)
#     for i in range(1, len(_col)):
#         _col[i] = i
#     mz_abd_table.columns = _col
#     mz_abd_table = mz_abd_table.set_index('m/z')
#     mz_abd_table.to_excel(mz_abd_addr)
#     return mz_abd_table
#

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


#
# def combine_profile_mz_with_motif_abundance(a_glycan_profile, NBT_dict_name_abundance_cross_profile):
#     """combine glycan m/z with motif hit, the maximum will be 4-8
#     return profile_obj_list
#     """
#     print("start combine")
#     NBT_glycan_match_existed_motif = load_json(__init__.json_address + r"NBT_glycan_match_existed_motif.json")
#     profile_obj_list = []
#     for _pro in range(1, len(a_glycan_profile) + 1):
#         if _pro < 10:
#             _id = '' + str(_pro)
#         else:
#             _id = 'Gly' + str(_pro)
#             # print(_id)
#         weighted_matrix = np.zeros((len(NBT_glycan_match_existed_motif["3865.1"])))
#         # print(weighted_matrix.shape)
#         abundance_ = []
#         # glycan_hit_array_ = []
#         mz_ = []
#         glycan_id_ = []
#         hit_matrix_ = []
#         for i in sorted(list(a_glycan_profile[_id].keys())):
#             # print(_pro, i)
#             _name = a_glycan_profile[_id][i]
#             mz_.append(i)
#             glycan_id_.append(_name)
#
#             _bundance = NBT_dict_name_abundance_cross_profile[i][_pro - 1]
#             # print(_pro, _bundance)
#             _temp_hit_matrix = np.array(NBT_glycan_match_existed_motif[_name])
#
#             abundance_.append(_bundance)
#             weighted_matrix += _temp_hit_matrix * _bundance
#             hit_matrix_.append(list(_temp_hit_matrix))
#         profile_obj_list.append(
#             glycan_profile_obj(glycan_id_, mz_, abundance_, weighted_matrix, hit_matrix_, name=str(_pro - 1)))
#     merged_profile_dict = {}
#
#     # print([round(i, 3) for i in merged_profile_dict[3]['motif_vec'][:20]])
#     for idex, i in enumerate(profile_obj_list):
#         merged_profile_dict[idex] = i.get_dict()
#     store_json(__init__.json_address + r"NBT_merged_profile_dict_merged.json", merged_profile_dict)
#     return profile_obj_list

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


def get_glycoprofile_list(profile_naming_to_id, norm_mz_abd_dict, match_dict, profile_name_order, external_profile_name,
                          glyprofile_list_addr,
                          get_existance=True):
    """combine glycan m/z with motif existances, convert the counts 1+/0 into 1/0 in each glycan
    Glycan motifs are modified with
    :param profile_name_order:
    :param get_existance:
    :param glyprofile_list_addr:
    :param profile_naming_to_id:
    :param norm_mz_abd_dict:
    :param match_dict:
    return profile_obj_list
    """

    profile_naming_to_id = check_profile_naming_to_id(profile_naming_to_id)
    profile_name_order = check_profile_naming_order(profile_name_order)

    if not external_profile_name:
        _ = list(profile_naming_to_id.keys())
        external_profile_name = dict(zip(_, _))
    else:
        external_profile_name = check_external_profile_name(external_profile_name)

    glycoprofile_list = []

    glycan_abd_dict = {}
    _num = len(profile_naming_to_id.keys())
    for i in match_dict.keys():
        glycan_abd_dict[i] = _num * [0]

    for pro_idex, pro in enumerate(profile_name_order):
        # pro_id = pro
        weighted_vec = np.zeros(len(match_dict[list(match_dict.keys())[0]]))
        abundance_list = []
        mz_list = []
        glycan_id_list = []
        match_mtrix = []

        for i in sorted(list(profile_naming_to_id[pro].keys())):
            glycan_id = profile_naming_to_id[pro][i]
            mz_list.append(i)
            glycan_id_list.append(glycan_id)

            _bundance = norm_mz_abd_dict[i][pro_idex]
            glycan_abd_dict[glycan_id][pro_idex] = _bundance
            _existance_list = []
            for _count in match_dict[glycan_id]:
                if _count >= 1:
                    if get_existance:
                        _existance_list.append(1)
                    else:
                        _existance_list.append(_count)
                elif _count == 0:
                    _existance_list.append(0)
                else:
                    assert False, 'wired in combine_profile_mz_with_motif_existance'

            _temp_hit_matrix = np.array(_existance_list)

            abundance_list.append(_bundance)
            match_mtrix.append(_temp_hit_matrix)
        # print(abundance_list)
        for idex, i in enumerate(abundance_list):
            # print(abundance_list
            weighted_vec += match_mtrix[idex] * i / sum(abundance_list)

        glycoprofile_list.append(
            Glycoprofile(glycan_id_list, mz_list, abundance_list, weighted_vec, hit_matrix=match_mtrix,
                         name=pro, profile_name=external_profile_name[pro]))

    glycoprofile_output_list = []
    # print([round(i, 3) for i in merged_profile_dict[3]['motif_vec'][:20]])
    for idex, i in enumerate(glycoprofile_list):
        glycoprofile_output_list.append(i.get_dict())
    store_json(glyprofile_list_addr, glycoprofile_output_list)
    # store_json(addr_root + r"glycoprofile_list.json", glycoprofile_output_list)
    return glycoprofile_list


class Glycoprofile():
    """we will have the profile with topology
    dict self.profile
    list self.motif_list
    profile contains glycan name and str structure, and motif information
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


# def combine_profile_mz(a_glycan_profile, NBT_dict_name_abundance_cross_profile):
#     """combine glycan m/z with motif existances, convert the counts into 1/0 in each glycan
#     return profile_obj_list
#     """
#     NBT_glycan_match_existed_motif = load_json(__init__.json_address + "NBT_glycan_match_count_motif.json")
#     profile_obj_list = []
#     for _pro in range(1, len(a_glycan_profile) + 1):
#         if _pro < 10:
#             _id = '' + str(_pro)
#         else:
#             _id = 'Gly' + str(_pro)
#             # print(_id)
#         weighted_matrix = np.zeros((len(NBT_glycan_match_existed_motif["3865.1"])))
#         # print(weighted_matrix.shape)
#         abundance_ = []
#         # glycan_hit_array_ = []
#         mz_ = []
#         glycan_id_ = []
#         hit_matrix_ = []
#         for i in sorted(list(a_glycan_profile[_id].keys())):
#             # print(_pro, i)
#             _name = a_glycan_profile[_id][i]
#             mz_.append(i)
#             glycan_id_.append(_name)
#
#             _bundance = NBT_dict_name_abundance_cross_profile[i][_pro - 1]
#             # print(_pro, _bundance)
#             _existance_list = []
#             for _count in NBT_glycan_match_existed_motif[_name]:
#                 if _count >= 1:
#                     _existance_list.append(1)
#                 elif _count == 0:
#                     _existance_list.append(0)
#                 else:
#                     assert False, 'wired in combine_profile_mz_with_motif_existance'
#
#             _temp_hit_matrix = np.array(_existance_list)
#
#             abundance_.append(_bundance)
#             weighted_matrix += _temp_hit_matrix * _bundance
#             hit_matrix_.append(list(_temp_hit_matrix))
#         profile_obj_list.append(
#             glycan_profile_obj(glycan_id_, mz_, abundance_, weighted_matrix, hit_matrix_, name=str(_pro - 1)))
#     merged_profile_dict = {}
#
#     # print([round(i, 3) for i in merged_profile_dict[3]['motif_vec'][:20]])
#     for idex, i in enumerate(profile_obj_list):
#         merged_profile_dict[idex] = i.get_dict()
#     store_json(__init__.json_address + r"NBT_merged_profile_dict_merged.json", merged_profile_dict)
#     return profile_obj_list


class MotifAbdTableGenerator:
    def __init__(self, glycoprofile_list):
        self.raw_table = glycoprofile_list[:]

    def table_btwn_two(self, n_a, n_b):
        #     change_f = np.array(a_p.weighted_motif_vec)/np.array(b_p.weighted_motif_vec)
        #     change_abs = np.array(a_p.weighted_motif_vec)/np.array(b_p.weighted_motif_vec)
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
        #     change_f = np.array(a_p.weighted_motif_vec)/np.array(b_p.weighted_motif_vec)
        #     change_abs = np.array(a_p.weighted_motif_vec)/np.array(b_p.weighted_motif_vec)
        a_p = self.raw_table[0]
        d = {a_p.name: a_p.match_vec_weighted}
        wt_table = pd.DataFrame(data=d)
        for idex, i in enumerate(self.raw_table):
            # if idex == 0:
            #     continue
            b_p = i
            # _weight = 1 / b_p.weighted_motif_vec[7]
            wt_table[b_p.name] = np.array(b_p.match_vec_weighted)
            #     wt_table = wt_table[(wt_table[a_p.name]+wt_table[b_p.name])!=0]
            """find out the abundance of the N-glycan core and use it to balance the weight """
        # wt_table = pd.DataFrame(data=d)
        return wt_table

    def table_existance(self):
        #     change_f = np.array(a_p.weighted_motif_vec)/np.array(b_p.weighted_motif_vec)
        #     change_abs = np.array(a_p.weighted_motif_vec)/np.array(b_p.weighted_motif_vec)
        a_p = self.raw_table[0]

        d = {a_p.name: a_p.match_vec_exist}
        existance_table = pd.DataFrame(data=d)
        for idex, i in enumerate(self.raw_table):
            # if idex == 0:
            #     continue
            b_p = i
            # _weight = 1 / b_p.weighted_motif_vec[7]
            existance_table[b_p.name] = np.array(b_p.match_vec_exist)
            #     wt_table = wt_table[(wt_table[a_p.name]+wt_table[b_p.name])!=0]
            """find out the abundance of the N-glycan core and use it to balance the weight """
        # wt_table = pd.DataFrame(data=d)
        return existance_table

    def table_against_wt_fc(self):
        #     change_f = np.array(a_p.weighted_motif_vec)/np.array(b_p.weighted_motif_vec)
        #     change_abs = np.array(a_p.weighted_motif_vec)/np.array(b_p.weighted_motif_vec)

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
        #     change_f = np.array(a_p.weighted_motif_vec)/np.array(b_p.weighted_motif_vec)
        #     change_abs = np.array(a_p.weighted_motif_vec)/np.array(b_p.weighted_motif_vec)

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
