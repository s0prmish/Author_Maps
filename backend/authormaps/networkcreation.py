from typing import Dict, Tuple

import graphviz
import math


# Add wrapper for two author names

# Task 4 - Build a Network (2 pts)
# Build an undirected graph in which nodes represent authors (i.e. names) and edges represent at least
# #1 shared publication between the pair of nodes

# Annotate the edges with the number of shared publications. The line thickness should correspond to the number of
# shared publications i.e. the more publications shared between the pair of authors, the thicker the edge

class AuthorNetwork:
    def __init__(self, shared_publications_counts: Dict[Tuple[str, str], int], enable_annotations: bool = True,
                 max_penwidth=7):
        """
        Creates a author map as graphviz graph
        :param shared_publications_counts: Dictionary of tuples author names as key and the number of shared
        publications as values (int).
        :param enable_annotations: Switch for the annotation of the Network with the labels of the edges and a change
        in line thickness.
        :param max_penwidth: Maximum width of the edges in the graph. this will be the scale for all entries.
        """
        self.data = shared_publications_counts
        self.graph = self.build_network(enable_annotations, max_penwidth=max_penwidth)

    def build_network(self, enable_annotations: bool = True, max_penwidth=7) -> graphviz.Graph:
        """
        Creates a undirected graph from the given data. Also allows to turn off the annotation
        (numbers next to the edges and varying line thickness). The graphviz.Graph object is then returned.
        :param enable_annotations:
        :param max_penwidth: Maximum width of the edges in the graph. this will be the scale for all entries.
        :return: Constructed Graph as graphviz.Graph
        """
        graphobject = graphviz.Graph("Authornetwork")

        highest_number = max(self.data.values())

        # find top 10% of edges
        top_ten_counts = math.ceil(len(self.data.keys()) / 10)
        sorted_edges = [k for k, v in sorted(self.data.items(), key=lambda item: item[1], reverse=True)]
        top_edges = sorted_edges[:top_ten_counts]

        # adding edges
        for authors, count in self.data.items():
            # print(authors[0])
            if enable_annotations:

                if authors in top_edges:
                    graphobject.edge(authors[0], authors[1], label=str(count),
                                     penwidth=str((count / highest_number) * max_penwidth), color="#d7191c")
                else:
                    # print(highest_number, count, max_penwidth)
                    graphobject.edge(authors[0], authors[1], label=str(count),
                                     penwidth=str((count / highest_number) * max_penwidth), color="0 0 0.4")
            else:
                graphobject.edge(authors[0], authors[1])

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
    def save_graph(self, output_format: str = "png", view=False):
        # where save the graph? Cache?

        viable_formats = ["png", "jpg", "svg", "pdf"]

        if output_format.startswith("."):
            output_format = output_format[1:]
        if output_format not in viable_formats:
            return False
        filename = "Authorgraph"
        print(filename)
        self.graph.render(filename, format=output_format, view=view, cleanup=True)  # os.path.join(folder,
        print("Rendered")
        return True


test_data = {("Ilya", "Marlo"): 3, ("Pragya", "Dhruv"): 4, ("Marlo", "Dhruv"): 2, ("Ilya", "Dhruv"): 1,
             ("Pragya", "Ilya"): 7}

#testobj = AuthorNetwork(test_data)  # ,enable_annotations=False)
#print("Created")
#print(testobj.visualize_as_string())
#testobj.save_graph("pdf", view=True)
