from ndex.networkn import NdexGraph
from ndex.client import Ndex
import networkx as nx

def test_create_from_edge_list():
    G = NdexGraph()
    edge_list = [('A', 'B'), ('B', 'C')]

    G.create_from_edge_list(edge_list, interaction=['A-B', 'B-C'])

    G.set_name('create_from_edge_list')

    network_id = G.upload_to("http://dev.ndexbio.org", "scratch", "scratch")
    print(network_id)

    #ndex = Ndex("http://dev.ndexbio.org", "scratch", "scratch")
    #ndex.make_network_public(network_id)
    #ndex.make_network_private(network_id)

#def test_cartesian():
 #   G = NdexGraph(server="http://test.ndexbio.org",
  #                username='scratch', password='scratch',
   #               uuid='aa6e7426-3f14-11e6-a7fa-028f28cd6a5b')
    #G.write_to('cartesian2.cx')

def test_layout():
    nx_G = nx.complete_graph(11)
    # nx_G.remove_node(0)
    # nx_G.remove_node(1)
    G = NdexGraph()
    G.add_new_node('1')
    G.add_new_node('2')
    G.add_new_node('3')
    G.add_new_node('4')
    G.add_new_node('5')
    G.add_new_node('6')

    G.add_edge_between(1, 3)
    G.add_edge_between(1, 4)
    G.add_edge_between(2, 3)
    G.add_edge_between(2, 4)

    G.add_edge_between(3, 5)
    G.add_edge_between(3, 6)
    G.add_edge_between(4, 5)
    G.add_edge_between(4, 6, interaction='AAABBBCCC') #attr_dict={'interaction':'testing'})

    initial_pos = {
        1: (0.0, 1.0),
        2: (0.0, 0.0),
        # 3: (0.5, 0.5),
        # 4: (0.5, 0.5),
        5: (1.0, 1.0),
        6: (1.0, 0.0),
    }

    fixed = [1,2,5,6]

    # G.add_new_node('3')
    # G.add_new_node('4')
    # G.add_new_node('5')
    G.pos = nx.spring_layout(G.to_undirected(), pos=initial_pos, fixed=fixed)
    print(G.pos)
    G.set_name('spring_layout undirected')
    #G.upload_to('http://dev.ndexbio.org','scratch','scratch')






def scratch_test():
    G = NdexGraph()
    G.add_new_node('1')
    G.add_new_node('2')
    G.add_edge_between(1,2)

    G.set_name('scratch-test1 - jing')
    print(G.graph())
    #G.upload_to('http://dev.ndexbio.org', 'scratch', 'scratch')


if __name__ == "__main__":
    test_layout()


