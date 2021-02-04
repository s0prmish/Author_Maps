import unittest
from backend.authormaps import author_network
import graphviz
from os import path
import os

test_data = {("Ilya", "Marlo"): 3, ("Pragya", "Dhruv"): 4, ("Marlo", "Dhruv"): 2, ("Ilya", "Dhruv"): 1,
             ("Pragya", "Ilya"): 7}


class TestAuthorNetwork:

    def test_init(self):
        test_graph = author_network.AuthorNetwork(test_data)
        assert test_graph
        assert test_graph.data == test_data
        assert type(test_graph.graph) == graphviz.Graph

    def test_build_network(self):
        tg = author_network.AuthorNetwork(test_data)
        print(tg.graph.source)
        assert 'Ilya -- Marlo [label=3 color="0 0 0.4" penwidth=3.0]' in tg.graph.source
        # Non annotated test! Maybe prone to randomness?
        tgnoannotation = author_network.AuthorNetwork(test_data, enable_annotations=False)
        nonannotext = """
	Ilya -- Marlo
	Pragya -- Dhruv
	Marlo -- Dhruv
	Ilya -- Dhruv
	Pragya -- Ilya
}"""
        assert nonannotext in tgnoannotation.graph.source

    def test_save_image(self):
        tg = author_network.AuthorNetwork(test_data)
        savepath = "Authorgraph.pdf"
        if path.exists(savepath):
            os.remove(savepath)
        assert not path.exists(savepath)

        tg.save_graph("pdf", view=False)

        assert path.exists(savepath)
        if path.exists(savepath):
            os.remove(savepath)


if __name__ == '__main__':
    unittest.main()
