import getopt
import sys


# root_ = "/Users/apple/PycharmProjects/GlyCompare/"
num_processors = 8

print("""
### There is one file must required and two files optional:
    1. abundance_table.xls

    if external profile name rather than index
    2. external_profile_naming.json

    if has mz or hplz name other than glycan_id
    3. glycoprofile_name_to_glycan_id.json

### There are several parameter should be set up first
    1. working_addr : root working dir

    2. project_name: usually same as the folder of root

    3. __init__.num_processors: number of processes needed

    4. __init__.exact_Ture: False for topology; True for exact linkage matching

    5. complex_profile: True is every profile has unique glycans; False is every profile has same name
    6. complex_naming: use the mz or hplc naming rather than glytoucan_id in abundance table.
    7. external_profile_naming: each glycoprofile has complex name rather than index.
        complex_profile and outside_name helps charactorize how complex of the naming is

    8. structure_loader: for loading glytoucan id, could be:
        list of glycan_id
        glycan_dict

    9. glytoucan_db_addr: the addr for glytoucan database if needed

    # meta_name: a gff type table which includes glycan's naming information

""")


if __name__ == '__main__':
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hdi:", ["help", "output="])
    except getopt.GetoptError:
        pass