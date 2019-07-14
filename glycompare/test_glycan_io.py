# import unittest
# from glypy.io import glycoct, iupac
#
# # import __init__
# from . import glycan_io
# from . import pipeline_functions
# import os
#
# """assertEqual(a, b)	a == b
# assertNotEqual(a, b)	a != b
# assertTrue(x)	bool(x) is True
# assertFalse(x)	bool(x) is False
# assertIs(a, b)	a is b	3.1
# assertIsNot(a, b)	a is not b	3.1
# assertIsNone(x)	x is None	3.1
# assertIsNotNone(x)	x is not None	3.1
# assertIn(a, b)	a in b	3.1
# assertNotIn(a, b)	a not in b	3.1
# assertIsInstance(a, b)	isinstance(a, b)	3.2
# assertNotIsInstance(a, b)	not isinstance(a, b)	3.2"""
#
# nglycan_core = """
# RES
# 1b:x-dglc-HEX-1:5
# 2s:n-acetyl
# 3b:b-dglc-HEX-1:5
# 4s:n-acetyl
# 5b:b-dman-HEX-1:5
# 6b:a-dman-HEX-1:5
# 7b:a-dman-HEX-1:5
# LIN
# 1:1d(2+1)2n
# 2:1o(4+1)3d
# 3:3d(2+1)4n
# 4:3o(4+1)5d
# 5:5o(3+1)6d
# 6:5o(6+1)7d """
#
# """root_ = "/Users/apple/PycharmProjects/GlyCompare/"
# num_processors = 8
# # exact_Ture = True
# json_address = root_ + "generated_json_file/"
# plot_output_address = root_ + "output_plot/"
# motif_plot_address = root_ + "motif_plot/"
# source_address = root_ + "source_data/"
# """
#
# root_ = root_
# # json_address = os.path.join(root_, "generated_json_file/")
# # plot_output_address = os.path.join(root_, "output_plot/")
# # source_address = os.path.join(root_, "source_data/")
#
#
# # test_glycan_dict = {'nglycan_core': glycoct.loads(nglycan_core)}
#
#
# class test_glycan_io(unittest.TestCase):
#     a_glycan_dict = {'nglycan_core': glycoct.loads(nglycan_core)}
#
#     def test_glycan_str_to_glycan(self):
#         pass
#
#     def test_glycan_to_glycan_str(self):
#         self.assertDictEqual({'nglycan_core': str(self.a_glycan_dict['nglycan_core'])},
#                              glycan_io.glycan_obj_to_glycan_str(self.a_glycan_dict))
#
#     def test_load_structure_pip(self):
#         b_glycan_dict = pipeline_functions.load_glycans_pip(project_name="test",
#                                                             working_dir=root_,
#                                                             data_type='glycan_dict',
#                                                             structure_loader=self.a_glycan_dict)
#         c_glycan_dict = pipeline_functions.load_glycans_pip(project_name="test",
#                                                             working_dir=root_,
#                                                             data_type='used',
#                                                             structure_loader=self.a_glycan_dict)
#         d_glycan_dict = pipeline_functions.load_glycans_pip(project_name="test",
#                                                             working_dir=root_,
#                                                             data_type='used',
#                                                             structure_loader=self.a_glycan_dict)
#         self.assertDictEqual(self.a_glycan_dict,
#                              b_glycan_dict, "Should be glypy glycan")
#         self.assertSequenceEqual(self.a_glycan_dict['nglycan_core'],
#                                  c_glycan_dict['nglycan_core'], "Should be glypy glycan")
#
#     def test_sum_tuple(self):
#         self.assertEqual(sum((1, 2, 3)), 6, "Should be 6")
#
#
#
#
# if __name__ == '__main__':
#     # unittest.main()
#     runner = unittest.TextTestRunner()
#     runner.run(test_glycan_io())
