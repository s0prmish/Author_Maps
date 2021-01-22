import unittest
from backend.authormaps import networkcreation
import graphviz
test_data={("Ilya","Marlo"):3,("Pragya","Dhruv"):4,("Marlo","Dhruv"):2,("Ilya","Dhruv"):1,("Pragya","Ilya"):7}

class Test_preparation(unittest.TestCase):
    def test_something(self):
        self.assertEqual(True, True)

class Test_auth_network():

    def test_init(self):
        test_graph=networkcreation.auth_network(test_data)
        assert test_graph!=None
        assert test_graph.data ==test_data
        assert type(test_graph.graph) == graphviz.Graph

    def test_build_network(self):
        tg = networkcreation.auth_network(test_data)
        assert "Marlo -- Dhruv [label=2 penwidth=3.4285714285714284]" in tg.graph.source
        #Non annotated test! Maybe prone to randomness?????????????
        tgnoannotation= networkcreation.auth_network(test_data,enable_annotations=False)
        nonannotext="""
	Ilya -- Marlo
	Pragya -- Dhruv
	Marlo -- Dhruv
	Ilya -- Dhruv
	Pragya -- Ilya
}"""
        assert nonannotext in tgnoannotation.graph.source



if __name__ == '__main__':
    unittest.main()
