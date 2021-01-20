from Bio import Entrez
from typing import List
import os
from pathlib import Path
import json

Entrez.email = 'rathodhruv007@gmail.com'

home_dir = str(Path.home())
PROJECT_DIR = os.path.join(home_dir, ".project", "group5")
DATA_DIR = os.path.join(PROJECT_DIR,"data")


class AuthorData:
    def __init__(self, author_name=None,IDs_to_retrive=20):
        self.author = author_name
        self.AUTHOR_DIR = os.path.join(DATA_DIR, self.author)
        os.makedirs(self.AUTHOR_DIR, exist_ok=True)
        self.retmax=IDs_to_retrive

    def get_list_of_publications(self) -> List['PMID']:

        type_ = "[Author - Full]"
        term = self.author + type_
        handle = Entrez.esearch(db="pubmed", term=term,retmax=self.retmax)
        record = Entrez.read(handle)
        idlist = record['IdList']
        print('Number of ids requested : {} '.format(self.retmax))
        print('Number of ids available on NCBI : {}'.format(len(idlist)))
        return idlist

    def __cache_publications(self) -> None:
        ids = self.get_list_of_publications()
        for i in ids:
            if str(i)+".json" in os.listdir(self.AUTHOR_DIR):
#                 print("passing")
                pass
            else:
#                 print("downloading")
                summ_handle_3 = Entrez.esummary(db="pubmed", id=i, version="2.0", retmode="json")
                data = json.load(summ_handle_3)
                path = os.path.join(self.AUTHOR_DIR, f"{i}.json")
                with open(path, 'w') as outfile:
                    json.dump(data, outfile)

    def get_list_of_coauthors_from_list_of_publications(self, publications: List['PMID']) -> List[str]:
        self.__cache_publications()
        out = {}
        for id_ in publications:
            handle = Entrez.esummary(db="pubmed", id=id_)
            record = Entrez.read(handle)
            out[id_] = record[0]['AuthorList']
        final = []
        for i in out.values():
            set_1 = set(final)
            set_2 = set(i)
            list_2_items_not_in_list_1 = list(set_2 - set_1)
            final = final + list_2_items_not_in_list_1
        return final

    def get_list_of_coauthors(self) -> List[str]:
        co_authors = self.get_list_of_coauthors_from_list_of_publications(self.get_list_of_publications())
        if co_authors == []:
            os.rmdir(self.AUTHOR_DIR)
            return "No such author found"
        else:
            print(co_authors)
            return co_authors



d=AuthorData('Svjetlana Miocinovic',5)
d.get_list_of_coauthors()

#for author not found
# d=AuthorData('xcfgsfg',50)
# d.get_list_of_coauthors()