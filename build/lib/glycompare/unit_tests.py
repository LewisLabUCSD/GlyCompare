import pandas as pd
import json
import os

def sub_abd_validation(keywords_dict, glycan_abd, get_existance, absolute):
    print("Start Test sub_abd_validation")
    assert os.path.isfile(keywords_dict['glycan_substructure_occurance_dict_addr']), "Test sub_abd_validation failed: Invalid glycan substructure occurance dict path " + keywords_dict['glycan_substructure_occurance_dict_addr']
    assert os.path.isfile(keywords_dict['substructure_abd_table_addr']), "Test sub_abd_validation failed: Invalid substructure abundance table path " + keywords_dict['substructure_abd_table_addr']
    
    sub_occur = pd.read_csv(keywords_dict['glycan_substructure_occurance_dict_addr'], index_col = 0)
#     sub_occur = json.load(open(keywords_dict['glycan_substructure_occurance_dict_addr'], "r"))
    sub_abd = pd.read_csv(keywords_dict['substructure_abd_table_addr'], index_col = 0)
  
    print("Generating standard substructure_abundance_table")
    sub_abd_new = pd.DataFrame()
    existance = get_existance
    for pro in glycan_abd.columns.tolist():
        weight_sub_vec = [0.0 for i in range(sub_occur.shape[0])]
        for gly in glycan_abd.index.tolist():
            g_abd = glycan_abd[pro][gly]
            if not existance:
                sub_vec = sub_occur[gly]
            else:
                sub_vec = [1 if i > 0 else 0 for i in sub_occur[gly].tolist()]
            temp_weight_sub_vec = [float(sub_vec[i] * g_abd) for i in range(len(sub_vec))]
            weight_sub_vec = [weight_sub_vec[i] + temp_weight_sub_vec[i] for i in range(len(weight_sub_vec))]
        sub_abd_new[str(pro)] = weight_sub_vec
        
    if not absolute:
        for col_ind, col in enumerate(sub_abd_new.columns):
            sub_abd_new[col] = sub_abd_new[col] / glycan_abd.sum()[col_ind]
    
    
    print("Comparing sub_abd with standard sub_abd")
    assert sub_abd.shape == sub_abd_new.shape, "Test sub_abd_validation failed: Wrong shape of substructure abundance table. Get " + str(sub_abd.shape) + ". Should be " + str(sub_abd_new.shape)

    errors = []
    for i in list(sub_abd.index):
        for j in list(sub_abd.columns):
            if abs(sub_abd[j][i] - sub_abd_new[j][i]) > 1e-10:
                errors.append((j, i))

    assert not errors, "Test sub_abd_validation failed: " + str(len(errors)) + " Mismatched entries found: \n" + "\n".join([str(i)[1:len(str(i)) - 1] for i in errors])
    print("Passed Test sub_abd_validation")

