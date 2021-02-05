from authormaps.authorinfo import AuthorData
import pytest
import os
import shutil




class TestAuthorData:
    """ Unit tests for AuthorData class"""

    @pytest.fixture(autouse=True)
    def run_before_and_after_tests(self):
        self.test_data1 = AuthorData("christiano ronaldo")
        self.test_data2 = AuthorData("Simmons Adam M")
        self.test_data3 = AuthorData("Nasioudis Andreas")

        yield
        shutil.rmtree(self.test_data1.AUTHOR_DIR, ignore_errors=True)
        shutil.rmtree(self.test_data2.AUTHOR_DIR, ignore_errors=True)
        shutil.rmtree(self.test_data3.AUTHOR_DIR, ignore_errors=True)

    def test_get_list_of_publications(self):

        assert type(self.test_data1.get_list_of_publications()) == list
        assert type(self.test_data2.get_list_of_publications()) == list
        assert type(self.test_data3.get_list_of_publications()) == list

        assert len(self.test_data1.get_list_of_publications()) == 0
        assert len(self.test_data2.get_list_of_publications()) >= 1
        assert len(self.test_data3.get_list_of_publications()) >= 1

    def test_get_list_of_coauthors(self):

        if self.test_data1.get_list_of_coauthors == []:
            assert os.path.isdir(self.test_data1.AUTHOR_DIR)

        if self.test_data2.get_list_of_coauthors != []:
            assert os.path.isdir(self.test_data2.AUTHOR_DIR)

        if self.test_data3.get_list_of_coauthors != []:
            assert os.path.isdir(self.test_data3.AUTHOR_DIR)

    def test_get_list_of_coauthors_from_list_of_publications(self):
        assert type(self.test_data1.get_list_of_coauthors_from_list_of_publications())==list
        assert type(self.test_data2.get_list_of_coauthors_from_list_of_publications()) == list
        assert type(self.test_data3.get_list_of_coauthors_from_list_of_publications()) == list
