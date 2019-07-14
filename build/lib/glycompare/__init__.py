import getopt
import sys
import pipeline_functions
import merge_substructure_vec


# root_ = "/Users/apple/PycharmProjects/GlyCompare/"
num_processors = 8

print("""
### An abundance table (1.1) is necessary to run the full GlyCompare pipeline:
    The abundance_table have column: sample (glycoprofile), row: glycan (the glycan identifier can be glytoucan_id, or custimized: m/z, hplc)
    1.1. abundance_table.xls

    Glycompare defaults to index names for each sample (e.g. 1,2,3). If you want to specify a name for each sample use:
    1.2. (optional) external_profile_naming.json

    If the glycan identifiers in abundance_table.xls are not glytoucan_id (if uses m/z or costumized). See 2.4.
    1.3. (optional) glycan_identifier_to_glytoucan_id.json


### There are several parameter should be set up first
    2.1. working_addr : root working dir

    2.2. project_name: usually same as the folder of root

    2.3. __init__.num_processors: number of processes needed

    2.4. complex_naming: 
        True if they use the mz, hplc naming or self-created naming rather than glytoucan_id in abundance table.
        
    2.5. complex_profile: 
        False, if glycan identifiers in abundance_table.xls is not glytoucan_id.
        
        True, if a glycan identifier, for example m/z, in abundance_table.xls has isomers, the m/z-structure matching 
        must be specified in glycan_identifier_to_glytoucan_id.json (1.3)
        
    #    glycan_profile['5'] = {
    #         "2244": "G04483SK",
    #         "2605": "G30460NZ",
    #         "2967": "G17689DH",
    #         "3416": "3416.1",
    #     }
    #    glycan_profile['6'] = {
    #         "1591": "G07483YN",
    #         "1836": "G80858MF",
    #         "2081": "G80393PG",
    #         "2326": "G00176HZ",
    #         "2530": "G79457WN",
    #         "2892": "G79412GP",
    #         "3096": "G40242TG", 
    #         "3457": "3457.1",
    #         "3777": "G76812VG",
    #         "4022": "G80264ZA",
    #         "4226": "G86696LV",
    #         "4587": "4587.1",
    #         "5037": "G49604DB",
    #     }
        False, if there is no ambiguity between abundance table index and glycan structure ID.

    
    2.6. external_profile_naming: as 1.2.
        True, 
        False

    2.7. structure_loader: for loading glytoucan id, could be:
        list of glycan_id
        glycan_dict

    2.8. glytoucan_db_addr: the addr for glytoucan database if needed

    # meta_name: a gff type table which includes glycan's naming information

""")


if __name__ == '__main__':
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hdi:", ["help", "output="])
    except getopt.GetoptError:
        pass