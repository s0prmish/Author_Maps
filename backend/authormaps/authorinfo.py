from Bio import Entrez
from typing import List
import os

import json
from Bio import Medline

from authormaps.startup import DATA_DIR

Entrez.email = 'rathodhruv007@gmail.com'


class AuthorData:
    """Downloads Author data as json and provides co-author information."""

    def __init__(self, author_name=None):

        self.author = author_name
        self.AUTHOR_DIR = os.path.join(DATA_DIR, self.author)
        os.makedirs(self.AUTHOR_DIR, exist_ok=True)
        self.retmax = 9999

    def get_list_of_publications(self) -> List[str]:
        """
        Creates the list of PMIDs os the queried author name. Maximum of 9999 records
        can be retrived.It also checks for updates at the NCBI side. Downloads even
        if one publication is added of given author.

        :return: List of PMIDs
        """
        type_ = "[Author]"
        term = self.author + type_
        idlist = []
        if os.listdir(self.AUTHOR_DIR)!=[]:
            for i in os.listdir(self.AUTHOR_DIR):
                idlist.append(i.split('.')[0])
        else:
            handle = Entrez.esearch(db="pubmed", term=term, retmax=self.retmax)
            record = Entrez.read(handle)
            for id in record['IdList']:
                idlist.append(id)
        return idlist

    def __cache_publications(self) -> None:
        """
        Adds all the data to folder named as author's name in the cache directory.
        :return: None
        """

        ids = self.get_list_of_publications()
        record_list = []
        for i in (ids):
            if str(i) + ".json" in os.listdir(self.AUTHOR_DIR):
                pass
            else:
                h = Entrez.efetch(db="pubmed", id=i, rettype="medline", retmode="text")
                records = Medline.parse(h)
                record_list.extend(list(records))

        for i in record_list:
            fp = os.path.join(self.AUTHOR_DIR, str(i['PMID']) + ".json")
            with open(fp, 'w') as fout:
                json.dump(i, fout)

    def get_list_of_coauthors_from_list_of_publications(self) -> List[str]:
        """
        Creates the list of co-authors. Deletes duplicates.

        :return: List of co-authors
        """
        self.__cache_publications()
        final = []
        for file in os.listdir(self.AUTHOR_DIR):
            path = os.path.join(self.AUTHOR_DIR, file)
            with open(path) as json_file:
                data = json.load(json_file)
                for i in data['FAU']:
                    i = i.replace(',', '')
                    final.append(i)
        final = list(dict.fromkeys(final))
        return final

    def get_list_of_coauthors(self) -> List[str]:
        """
        Creates the list of co-authors from author name
        :return: List of co-authors
        """
        co_authors = self.get_list_of_coauthors_from_list_of_publications()
        if not co_authors:
            os.rmdir(self.AUTHOR_DIR)
        return co_authors

