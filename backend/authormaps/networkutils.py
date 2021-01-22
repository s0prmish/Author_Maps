from typing import List, Dict, Tuple
from itertools import combinations


# Task 2: Task 2 - Finding the Shared Work for All-Pairs of Authors (1 pt)


def get_publications_for_pair(authors: Tuple[str, str]) -> List[str]:
    """
   Returns the union of two (as tuple) given authors publications
   :param authors: Tuple of two authors
   :return: A list of publications that both authors are involved in
   """
    authors1_publications = None  # get_list_of_publications(authors[0])
    authors2_publications = None  # get_list_of_publications(authors[1])
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
    For a given list of authors, this function returns a dictionary with the number of shared publications
    between all combinations of given authors. Only one tuple per pair will be created.
    :param authors:
    :return: dictionary where keys are made
    """
    shared_publications = get_shared_publications_for_authors_pairs(authors)
    result = {k: len(v) for k, v in shared_publications.items()}
    return result
