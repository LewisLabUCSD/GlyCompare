# glycan_dict['5486.1']



mono2markov = {'RES 1b:b-dglc-HEX-1:5 2s:n-acetyl LIN 1:1d(2+1)2n':'bNG',
               'RES 1b:x-dglc-HEX-1:5 2s:n-acetyl LIN 1:1d(2+1)2n':'bNG',
               'RES 1b:x-dman-HEX-1:5':'bm',
               'RES 1b:b-dman-HEX-1:5':'bM',
               'RES 1b:a-dman-HEX-1:5':'aM',
               'RES 1b:b-dgal-HEX-1:5':'bA',
               'RES 1b:x-dgal-HEX-1:5':'bA',
               'RES 1b:a-dgro-dgal-NON-2:6|1:a|2:keto|3:d 2s:n-acetyl LIN 1:1d(5+1)2n':'aNN',
               'RES 1b:a-lgal-HEX-1:5|6:d':'aF',
              }

# for i in glycan_dict['5486.1'].index:
#     print(i)
#     print(mono2markov[str(i)])

def _glycoct2linear(a_mono):
    children = a_mono.children()
    children = sorted(children, key=lambda x:x[0], reverse=True)
#     print(children)
    _str=''
    if len(children) == 1:
        index, mono = children[0]
        return  str(index) + mono2markov[str(mono)] + _glycoct2linear(mono)
    elif len(children) > 1:
        for index, mono in list(children)[:-1]:
            _str += ')'+str(index)+mono2markov[str(mono)] + _glycoct2linear(mono) + '('
        index, mono = children[-1]
#         print(index, mono)
        _str += str(index) + mono2markov[str(mono)] + _glycoct2linear(mono)
        return _str
    else:
        return ''

def glycoct2linear(a_glycan):
    a = 'nsA;NG' + _glycoct2linear(a_glycan.root)

    return a[::-1]
# plot_glycan_utilities.plot_glycan(glycan_dict['5486.1'])
# print(glycoct2markov(glycan_dict['5486.1']))