from typing import List,  Dict, Tuple
from itertools import combinations
import graphviz
import os

# Task 2: Task 2 - Finding the Shared Work for All-Pairs of Authors (1 pt)


def get_publications_for_pair(authors: Tuple[str, str]) -> List[str]:
    """
   Returns the union of two (as tuple) given authors publications
   :param authors: Tuple of two authors
   :return: A list of publications that both authors are involved in
   """
    authors1_publications = None# get_list_of_publications(authors[0])
    authors2_publications = None# get_list_of_publications(authors[1])
    # return intersection between authors1_publications  and authors2_publications
    return list(set(authors1_publications) & set(authors2_publications))


def get_shared_publications_for_authors_pairs(authors: List[str]) -> Dict[Tuple[str, str], List[str]]:
    """
    from a given list of authors, it collects information about all shared publications between them.
    It returns a dictionary where the key is the tuple of authors and the value is a list of publications they share.
    Warning!: The key pairs are only stored in one direction (eg (A,B):[1,2,3] but not (B,A):[1,2,3] in addition)
    :param authors: a list of authors for which the shared publications should be checked
    :return: a dictionary with author tuples as key and the list of shared publications as value (only one key per pair)
    """
    result = {}
    for combination in combinations(authors, 2):
        publications = get_publications_for_pair(combination)
        result[combination] = publications

    # check that there is no redundancy
    # combination should lead to no redundancy.. permutations would but i think this saves more memory
    return result


# Task 3
def count_shared_publications(authors: List[str]) -> Dict[Tuple[str, str], int]:
    """
    For a given list of authors, this fuction returns a dictionary with the number of shared publications
    between all combinations of given authors. Only one tuple per pair will be created.
    :param authors:
    :return: dictionary where keys are made
    """
    shared_publications = get_shared_publications_for_authors_pairs(authors)
    result = {k: len(v) for k, v in shared_publications.items()}
    return result


# Task 4 - Build a Network (2 pts)
# Build an undirected graph in which nodes represent authors (i.e. names) and edges represent at least
# #1 shared publication between the pair of nodes

# Annotate the edges with the number of shared publications. The line thickness should correspond to the number of
# shared publications i.e. the more publications shared between the pair of authors, the thicker the edge

class auth_network():
    def __init__(self, shared_publications_counts: Dict[Tuple[str, str], int], enable_annotations: bool = True):
        """
        Creates a authormap as grphviz graph
        :param shared_publications_counts: Dictionary of tuples author names as key and the number of shared
        publications as values (int).
        :param enable_annotations: Switch for the annotation of the Network with the labels of the edges and a change
        in line thickness.
        """
        self.data = shared_publications_counts
        self.graph = self.build_network(enable_annotations)


    def _get_different_authors(self):
        """
        Helper function to extract all different authors from self.data (the dictionary that was given
        upon initialisation.
        :return: List of different authors.
        """
        result = set()
        for pair in self.data.keys():
            result.add(pair[0])
            result.add(pair[1])
        return list(result)

    def build_network(self, enable_annotations: bool = True) -> graphviz.Graph:
        """
        Creates a undircrected graph from the given data. Also allows to turn off the annotation
        (numbers next to the edges and varying line thickness). The graphviz.Graph object is then returned.
        :param enable_annotations:
        :return: Constructed Graph as graphviz.Graph
        """
        graphobject = graphviz.Graph("Authornetwork")

        # adding nodes
        for author in self._get_different_authors():
            graphobject.node(author, author)

        # adding edges
        for authors,count in self.data.items():
            #print(authors[0])
            if enable_annotations:
                graphobject.edge(authors[0],authors[1], label=str(count), penwidth=str(count))
            else:
                graphobject.edge((authors[0],authors[1]))

        return graphobject

    # Add additional functionality so that one can visualize the network
    # def show_graph(graph: graphviz.Graph)
    def visualize_as_string(self) -> str:
        """
        Visualises the stored graph as a string with all the attributes of graphviz
        :return: String with all the nodes, edges and their attributes
        """
        return self.graph.source

    # graph: graphviz.Graph,
    # Allow one to be able to export the network in several formats including png, jpg, svg, and pdf
    def save_graph(self,format: str = "png",view=False):
        #where save the graph? Cache?
        folder=""
        viable_formats=["png", "jpg", "svg", "pdf"]

        if format.startswith("."):
            format=format[1:]
        if format not in  viable_formats:
            return False
        filename="Authorgraph"
        print(filename)
        self.graph.render( filename,format=format,view=view,cleanup=True)#os.path.join(folder,
        return True

test_data={("Ilya","Marlo"):3,("Pragya","Dhruv"):4,("Marlo","Dhruv"):2,("Ilya","Dhruv"):1,("Pragya","Ilya"):7}

testobj=auth_network(test_data)
print("Created")
print(testobj.visualize_as_string())
testobj.save_graph("pdf")










