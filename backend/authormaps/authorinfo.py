from Bio import Entrez
from typing import List
import os
from pathlib import Path
import json
from Bio import Medline

from authormaps.startup import DATA_DIR

Entrez.email = 'rathodhruv007@gmail.com'

class AuthorData:
    """Downloads Author data as json and provides co-author information."""
    def __init__(self, author_name=None,IDs_to_retrive=9999):

        self.author = author_name
        self.AUTHOR_DIR = os.path.join(DATA_DIR, self.author)
        os.makedirs(self.AUTHOR_DIR, exist_ok=True)
        self.retmax = 9999

    def get_list_of_publications(self) -> List['PMID']:
        """
        Creates the list of PMIDs os the queried author name. Maximum of 9999 records
        can be retrived.

        :return: List of PMIDs
        """
        type_ = "[Author]"
        term = self.author + type_
        handle = Entrez.esearch(db="pubmed", term=term,retmax=self.retmax)
        record = Entrez.read(handle)
        idlist = record['IdList']
        # print('Number of ids requested : {} '.format(self.retmax))
        # print('Number of ids available on NCBI : {}'.format(len(idlist)))
        return idlist

    def __cache_publications(self) -> None:
        """
        Creates the cache folder in user home directory  Adds all the data to folder named as author's name
        :return: None
        """
        ids = self.get_list_of_publications()
        for i in ids:
            if str(i) + ".json" in os.listdir(self.AUTHOR_DIR):
                pass
            else:
                print('downloading')
                record_list = []

                for i in ids:
                    h = Entrez.efetch(db="pubmed", id=i, rettype="medline", retmode="text")
                    records = Medline.parse(h)
                    record_list.extend(list(records))

                for i in record_list:
                    fp = os.path.join(self.AUTHOR_DIR, str(i['PMID']) + ".json")
                    with open(fp, 'w') as fout:
                        json.dump(i, fout)

    def get_list_of_coauthors_from_list_of_publications(self, publications: List['PMID']) -> List[str]:
        """
        Creates the list of co-authors. Deleted duplicates.
        :param publications:List of PMIDs

        :return: List of co-authors
        """
        self.__cache_publications()
        final = []
        for file in os.listdir(self.AUTHOR_DIR):
            path=os.path.join(self.AUTHOR_DIR,file)
            with open(path) as json_file:
                data = json.load(json_file)
                for i in data['FAU']:
                    i=i.replace(',','')
                    final.append(i)
        final = list(dict.fromkeys(final))
        return final

    def get_list_of_coauthors(self) -> List[str]:
        """
        Creates the list of co-authors from author name
        :return: List of co-authors
        """
        co_authors = self.get_list_of_coauthors_from_list_of_publications(self.get_list_of_publications())
        if co_authors == []:
            os.rmdir(self.AUTHOR_DIR)
            return "No such author found"
        else:
            print(co_authors)
            return co_authors

d=AuthorData('Nasioudis Andreas')
d.get_list_of_coauthors()

#for author not found
# d=AuthorData('xcfgsfg')
# d.get_list_of_coauthors()
# print(DATA_DIR)