# Group 05 - Project 02 - Author Maps

The package allows you to generate a co-author network which details how strongly connected the queried author is to his/her co-authors.

It supports the interaction using either CLI or Web UI.

## Installation

The package can be installed as a regular python package using `pip`:
```
cd group05/backend
pip install .
```
Also, the package requires to install GraphViz software. This can be done with one of 2 approaches:
* a standalone installation (see https://www.graphviz.org/download/)
* conda installation:
  ```
  conda install graphviz
  conda install python-graphviz
  ```

Potential problem: it can happen that the GraphViz executables are not found in PATH. The solution is to add the GraphViz path to the PATH environment variable. For example, on Windows you might need to go to `Control Panel > System and Security > System > Advanced System Settings > Environment Variables > Path > Edit > New` and add the location of GraphViz (e.g. `C:\Users\UserName\Anaconda3\Library\bin\graphviz`).

## CLI
The package supports 6 CLI commands:
* `get_publications AUTHOR_NAME` - get list of PMIDs of publications for a given author;
* `get_coauthors AUTHOR_NAME` - get list of coauthors names from author's publications;
* `get_common_publications AUTHOR_NAME_1 AUTHOR_NAME_2` - get list of PMID's of publications shared by 2 given authors;
* `get_common_publications_for_all_coauthors AUTHOR_NAME` - get PMID's of common publications for every author's coauthor pair;
* `get_common_publications_counts_for_all_coauthors AUTHOR_NAME` - get count of common publications for every author's coauthor pair;
* `generate_graph AUTHOR_NAME OUTPUT_FILE [FORMAT]` - generate graph of coworkers for a given author.

### CLI usage examples
```
> python backend/authormaps/cli.py get_publications "Birkenbihl Colin"
['33344750', '33285634', '32843907', '31993440']
```

```
> python backend/authormaps/cli.py get_coauthors "Birkenbihl Colin"
['Golriz Khatami Sepehr', 'Robinson Christine', 'Birkenbihl Colin', 'Domingo-Fernandez Daniel', 'Hoyt Charles Tapley', 'Hofmann-Apitius Martin', 'Emon Mohammad Asif', 'Vrooman Henri',
'Westwood Sarah', 'Lovestone Simon', 'Frohlich Holger', 'Shi Liu', 'Nevado-Holgado Alejo', 'Westman Eric', 'Salimi Yasamin']
```

```
> python backend/authormaps/cli.py get_common_publications "Birkenbihl Colin" "Domingo-Fernandez Daniel"
['33344750', '31993440']
```

```
> python backend/authormaps/cli.py get_common_publications_for_all_coauthors "Birkenbihl Colin"
{('Golriz Khatami Sepehr', 'Robinson Christine'): ['31993440'], ('Golriz Khatami Sepehr', 'Birkenbihl Colin'): ['31993440'], ('Golriz Khatami Sepehr', 'Domingo-Fernandez Daniel'): ['31993440'], ('Golriz Khatami Sepehr', 'Hoyt Charles Tapley'): ['31993440'], ('Golriz Khatami Sepehr', 'Hofmann-Apitius Martin'): ['31993440', '32073441'], ('Robinson Christine', 'Birkenbihl Colin'): ['31993440'], ('Robinson Christine', 'Domingo-Fernandez Daniel'): ['31993440'], ('Robinson Christine', 'Hoyt Charles Tapley'): ['31993440'], ('Robinson Christine', 'Hofmann-Apitius Martin'): ['31993440'], ('Birkenbihl Colin', 'Domingo-Fernandez Daniel'): ['31993440', '33344750'], ('Birkenbihl Colin', 'Hoyt Charles Tapley'): ['31993440'], ('Birkenbihl Colin', 'Hofmann-Apitius Martin'): ['31993440', '33344750', '33285634', '32843907'], ('Birkenbihl Colin', 'Emon Mohammad Asif'): ['32843907'], ('Birkenbihl Colin', 'Vrooman Henri'): ['32843907'], ('Birkenbihl Colin', 'Westwood Sarah'): ['33285634', '32843907'], ('Birkenbihl Colin', 'Lovestone Simon'): ['33344750', '33285634', '32843907'], ('Birkenbihl Colin', 'Frohlich Holger'): ['33344750', '32843907'], ('Birkenbihl Colin', 'Shi Liu'): ['33285634'], ('Birkenbihl Colin', 'Nevado-Holgado Alejo'): ['33285634'], ('Birkenbihl Colin', 'Westman Eric'): ['33285634'], ('Birkenbihl Colin', 'Salimi Yasamin'): ['33344750'], ('Domingo-Fernandez Daniel', 'Hoyt Charles Tapley'): ['31092193', '32411185', '32211383', '30854225', '31993440', '29873705', '30564458', '30768158', '31225582', '30576488', '32503412', '31824580'], ('Domingo-Fernandez Daniel', 'Hofmann-Apitius Martin'): ['31092193', '31260040', '32211383', '28651363', '32976572', '31225582', '30576488', '32503412', '33154531', '33344750', '30564458', '33264280', '30854225', '30042519', '32925069', '32411185', '33367476', '29873705', '31993440', '31824580'], ('Domingo-Fernandez Daniel', 'Emon Mohammad Asif'): ['28651363', '30042519', '32503412', '33154531'], ('Domingo-Fernandez Daniel', 'Vrooman Henri'): ['33154531'], ('Domingo-Fernandez Daniel', 'Lovestone Simon'): ['33344750'], ('Domingo-Fernandez Daniel', 'Frohlich Holger'): ['33344750', '32411185', '31824580', '30042519', '33154531'], ('Domingo-Fernandez Daniel', 'Salimi Yasamin'): ['33344750'], ('Hoyt Charles Tapley', 'Hofmann-Apitius Martin'): ['31092193', '32411185', '32211383', '30854225', '31993440', '31604427', '29949955', '29873705', '30564458', '32750869', '31225582', '30576488', '32503412', '31824580'], ('Hoyt Charles Tapley', 'Emon Mohammad Asif'): ['32503412'], ('Hoyt Charles Tapley', 'Frohlich Holger'): ['32411185', '31824580'], ('Hofmann-Apitius Martin', 'Emon Mohammad Asif'): ['32843907', '28651363', '28035920', '32620927', '31730697', '30042519', '32503412', '33154531'], ('Hofmann-Apitius Martin', 'Vrooman Henri'): ['28731430', '32843907', '32620927', '31730697', '33154531'], ('Hofmann-Apitius Martin', 'Westwood Sarah'): ['33285634', '32843907'], ('Hofmann-Apitius Martin', 'Lovestone Simon'): ['32843907', '33344750', '33285634', '25754460'], ('Hofmann-Apitius Martin', 'Frohlich Holger'): ['33344750', '32411185', '32843907', '32620927', '31824580', '31730697', '30042519', '33154531'], ('Hofmann-Apitius Martin', 'Shi Liu'): ['33285634'], ('Hofmann-Apitius Martin', 'Nevado-Holgado Alejo'): ['33285634'], ('Hofmann-Apitius Martin', 'Westman Eric'): ['33285634'], ('Hofmann-Apitius Martin', 'Salimi Yasamin'): ['33344750'], ('Emon Mohammad Asif', 'Vrooman Henri'): ['31730697', '32620927', '33154531', '32843907'], ('Emon Mohammad Asif', 'Westwood Sarah'): ['32843907'], ('Emon Mohammad Asif', 'Lovestone Simon'): ['32843907'], ('Emon Mohammad Asif', 'Frohlich Holger'): ['32843907', '32620927', '31730697', '30042519', '33154531'], ('Vrooman Henri', 'Westwood Sarah'): ['32843907'], ('Vrooman Henri', 'Lovestone Simon'): ['32843907'], ('Vrooman Henri', 'Frohlich Holger'): ['31730697', '32620927', '33154531', '32843907'], ('Westwood Sarah', 'Lovestone Simon'): ['29562526', '30853464', '32831200', '32843907', '31495601', '31078433', '28441961', '27031486', '27239491', '31985466', '31890857', '26635716', '31047856', '30618716', '33285634'], ('Westwood Sarah', 'Frohlich Holger'): ['32843907'], ('Westwood Sarah', 'Shi Liu'): ['29562526', '32831200', '31495601', '31985466', '33285634'], ('Westwood Sarah', 'Nevado-Holgado Alejo'): ['28441961', '31078433', '32831200', '31495601', '31985466', '31890857', '31047856', '30618716', '33285634'], ('Westwood Sarah', 'Westman Eric'): ['33285634'], ('Lovestone Simon', 'Frohlich Holger'): ['27079753', '33344750', '32843907'], ('Lovestone Simon', 'Shi Liu'): ['29562526', '32831200', '31495601', '31985466', '33285634'], ('Lovestone Simon', 'Nevado-Holgado Alejo'): ['27872687', '30108252', '28441961', '31078433', '32831200', '31985466', '31890857', '27903560', '30775436', '30618716', '33285634', '28704727', '33303834', '32252806', '31072055', '27911300', '29954452', '31495601', '31047856'], ('Lovestone Simon', 'Westman Eric'): ['20447732', '31636452', '25370484', '20603455', '21971470', '26726737', '19683363', '33491853', '28098162', '27694991', '20800095', '23028511', '33285634', '26560730', '25012867', '31848485', '24121966', '24768341', '20693633', '25141298', '28394034', '22205954', '26440606', '21157852', '26401776', '30072774', '21763442', '25607358', '20847553', '22246242', '28157104', '29455029', '26973105', '26284520', '25649652', '24326516', '19906260', '20847449', '33007638', '21811624', '25071554', '31551603', '30705627', '23047370', '23541334'], ('Lovestone Simon', 'Salimi Yasamin'): ['33344750'], ('Frohlich Holger', 'Salimi Yasamin'): ['33344750'], ('Shi Liu', 'Nevado-Holgado Alejo'): ['31495601', '33285634', '32831200', '31985466'], ('Shi Liu', 'Westman Eric'): ['33285634'], ('Nevado-Holgado Alejo', 'Westman Eric'): ['33285634']}
```

```
> python backend/authormaps/cli.py get_common_publications_counts_for_all_coauthors "Birkenbihl Colin"
{('Golriz Khatami Sepehr', 'Robinson Christine'): 1, ('Golriz Khatami Sepehr', 'Birkenbihl Colin'): 1, ('Golriz Khatami Sepehr', 'Domingo-Fernandez Daniel'): 1, ('Golriz Khatami Sepehr', 'Hoyt Charles Tapley'): 1, ('Golriz Khatami Sepehr', 'Hofmann-Apitius Martin'): 2, ('Robinson Christine', 'Birkenbihl Colin'): 1, ('Robinson Christine', 'Domingo-Fernandez Daniel'): 1, ('Robinson Christine', 'Hoyt Charles Tapley'): 1, ('Robinson Christine', 'Hofmann-Apitius Martin'): 1, ('Birkenbihl Colin', 'Domingo-Fernandez Daniel'): 2, ('Birkenbihl Colin', 'Hoyt Charles Tapley'): 1, ('Birkenbihl Colin', 'Hofmann-Apitius Martin'): 4, ('Birkenbihl Colin', 'Emon Mohammad Asif'): 1, ('Birkenbihl Colin', 'Vrooman Henri'): 1, ('Birkenbihl Colin', 'Westwood Sarah'): 2, ('Birkenbihl Colin', 'Lovestone Simon'): 3, ('Birkenbihl Colin', 'Frohlich Holger'): 2, ('Birkenbihl Colin', 'Shi Liu'): 1, ('Birkenbihl Colin', 'Nevado-Holgado Alejo'): 1, ('Birkenbihl Colin', 'Westman Eric'): 1, ('Birkenbihl Colin', 'Salimi Yasamin'): 1, ('Domingo-Fernandez Daniel', 'Hoyt Charles Tapley'): 12, ('Domingo-Fernandez Daniel', 'Hofmann-Apitius Martin'): 20, ('Domingo-Fernandez Daniel', 'Emon Mohammad Asif'): 4, ('Domingo-Fernandez Daniel', 'Vrooman Henri'): 1, ('Domingo-Fernandez Daniel', 'Lovestone Simon'): 1, ('Domingo-Fernandez Daniel', 'Frohlich Holger'): 5, ('Domingo-Fernandez Daniel', 'Salimi Yasamin'): 1, ('Hoyt Charles Tapley', 'Hofmann-Apitius Martin'): 14, ('Hoyt Charles Tapley','Emon Mohammad Asif'): 1, ('Hoyt Charles Tapley', 'Frohlich Holger'): 2, ('Hofmann-Apitius Martin', 'Emon Mohammad Asif'): 8, ('Hofmann-Apitius Martin', 'Vrooman Henri'): 5, ('Hofmann-Apitius Martin', 'Westwood Sarah'): 2, ('Hofmann-Apitius Martin', 'Lovestone Simon'): 4, ('Hofmann-Apitius Martin', 'Frohlich Holger'): 8, ('Hofmann-Apitius Martin', 'Shi Liu'): 1, ('Hofmann-Apitius Martin', 'Nevado-Holgado Alejo'): 1, ('Hofmann-Apitius Martin', 'Westman Eric'): 1, ('Hofmann-Apitius Martin', 'Salimi Yasamin'): 1, ('Emon Mohammad Asif', 'Vrooman Henri'): 4, ('Emon Mohammad Asif', 'Westwood Sarah'): 1, ('Emon Mohammad Asif', 'Lovestone Simon'): 1, ('Emon Mohammad Asif', 'Frohlich Holger'): 5, ('Vrooman Henri', 'Westwood Sarah'): 1, ('Vrooman Henri', 'Lovestone Simon'): 1, ('Vrooman Henri', 'Frohlich Holger'): 4, ('Westwood Sarah', 'Lovestone Simon'): 15, ('Westwood Sarah', 'Frohlich Holger'): 1, ('Westwood Sarah', 'Shi Liu'): 5, ('Westwood Sarah', 'Nevado-Holgado Alejo'): 9, ('Westwood Sarah', 'Westman Eric'): 1, ('Lovestone Simon', 'Frohlich Holger'): 3, ('Lovestone Simon', 'Shi Liu'): 5, ('Lovestone Simon', 'Nevado-Holgado Alejo'): 19, ('Lovestone Simon', 'Westman Eric'): 45, ('Lovestone Simon', 'Salimi Yasamin'): 1, ('Frohlich Holger', 'Salimi Yasamin'): 1, ('Shi Liu','Nevado-Holgado Alejo'): 4, ('Shi Liu', 'Westman Eric'): 1, ('Nevado-Holgado Alejo', 'Westman Eric'): 1}
```

```
> python backend/authormaps/cli.py generate_graph "Birkenbihl Colin" "graph.png" -f png
```

## Web UI
You can run Web UI by running `python group05/web/run.py`

By default, the application will be deployed to the port `5000`.

* You can navigate between the home and about page via the buttons on top. 
* In order to find out how many co-author's does an author have, you can input the author's name in the textbox and press enter.
* Once that's done, you will get a tabulated list of all the co-authors that qxist for the queried author.
* If you want to visualize a co-author network: select any two co-authors and click on Submit. You will then see the selected co-author's nodes colored in green, and the others in white.
* If you want to save this image, you can select one of the options given (pdf, png, svg, jpg) and click the "Download" button.
