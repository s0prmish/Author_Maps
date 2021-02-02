from authormaps.authorinfo import AuthorData
import pytest
import Bio
import os

test_data1=AuthorData("zcksdj")
test_data2=AuthorData("Svjetlana Miocinovic")
test_data3=AuthorData("William joyce")

class TestAuthorData:
    """ Unit tests for AuthorData class"""

    def test_get_list_of_publications(self):
        assert type(test_data1.get_list_of_publications())==list
        assert type(test_data2.get_list_of_publications()) == list
        assert type(test_data3.get_list_of_publications()) == list

        assert len(test_data1.get_list_of_publications())==0
        assert len(test_data2.get_list_of_publications()) >= 1
        assert len(test_data3.get_list_of_publications()) >= 1

    def test_get_list_of_coauthors(self):
        if test_data1.get_list_of_coauthors==[]:
            assert os.path.isdir(test_data1.AUTHOR_DIR)

        if test_data2.get_list_of_coauthors!=[]:
            assert os.path.isdir(test_data2.AUTHOR_DIR)

        if test_data3.get_list_of_coauthors!=[]:
            assert os.path.isdir(test_data3.AUTHOR_DIR)


