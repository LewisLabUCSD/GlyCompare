import getopt
import sys
from . import pipeline_functions
from . import merge_substructure_vec
from . import glycan_io
from . import json_utility
from . import nglycan_alignment
from . import plot_glycan_utilities
from . import clustering_analysis
from . import extract_substructures
from . import process_glycoprofiles
from . import test_glycan_io

# root_ = "/Users/apple/PycharmProjects/GlyCompare/"
# num_processors = 8

print("""
Thanks for using the GlyCompare v1.0, 
Please check our github for the latest update.


Bokan & Ben,
08/1/2019

""")


if __name__ == '__main__':
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hdi:", ["help", "output="])
    except getopt.GetoptError:
        pass