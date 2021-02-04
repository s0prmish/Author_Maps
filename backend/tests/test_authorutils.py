import unittest
from backend.authormaps import networkutils

test_data = {("Ilya", "Marlo"): 3, ("Pragya", "Dhruv"): 4, ("Marlo", "Dhruv"): 2, ("Ilya", "Dhruv"): 1,
             ("Pragya", "Ilya"): 7}


class AuthorDataFactoryDummy():
    def __init__(self, authorname: str):
        self.author = authorname

    def get_list_of_publications(self):
        name = self.author
        if name == "Test1":
            return ["A", "B", "C"]
        elif name == "Test2":
            return ["C"]
        elif name == "Test3":
            return ["A", "B"]
        elif name == "Test4":
            return ["A", "B", "C", "D"]
        else:
            return []


class TestAuthorNetwork:
    def test_get_publication_for_pair(self):
        t1 = "Test1"
        t2 = "Test2"
        t3 = "Test3"
        t4 = "Test4"
        t5 = "Test5"

        # One intersection
        ret_12 = networkutils.get_publications_for_pair((t1, t2), author_data_factory=AuthorDataFactoryDummy)
        assert type(ret_12) == list
        assert len(ret_12) == 1
        assert ret_12[0] == "C"

        # Two intersection
        ret_13 = networkutils.get_publications_for_pair((t1, t3), author_data_factory=AuthorDataFactoryDummy)
        assert type(ret_13) == list
        assert len(ret_13) == 2
        assert "A" in ret_13 and "B" in ret_13

        # Three intersection
        ret_14 = networkutils.get_publications_for_pair((t1, t4), author_data_factory=AuthorDataFactoryDummy)
        assert type(ret_14) == list
        assert len(ret_14) == 3
        assert "A" in ret_14 and "B" in ret_14 and "C" in ret_14

        # No intersection
        ret_15 = networkutils.get_publications_for_pair((t1, t5), author_data_factory=AuthorDataFactoryDummy)
        assert type(ret_15) == list
        assert len(ret_15) == 0

        # print(ret_12,ret_13,ret_14,ret_15)

    def test_get_shared_publications_for_authors_pairs(self):
        cd = AuthorDataFactoryDummy
        t1 = "Test1"
        t2 = "Test2"
        t3 = "Test3"
        t4 = "Test4"
        t5 = "Test5"
        tl = [t1, t2, t3, t4, t5]
        test_dict = networkutils.get_shared_publications_for_authors_pairs(tl, author_data_factory=cd)
        print(test_dict)
        assert type(test_dict) == dict

        # check keys
        ks = test_dict.keys()
        assert (t1, t2) in ks
        assert (t2, t3) not in ks

        # check contents
        content = test_dict[(t1, t3)]
        assert type(content) == list
        assert content == ["A", "B"] or content == ["B", "A"]


if __name__ == '__main__':
    unittest.main()
