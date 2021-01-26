from flask import Flask, flash, request, redirect, url_for, render_template, session
from authormaps.authorinfo import AuthorData
from authormaps.networkcreation import AuthorNetwork
from authormaps.networkutils import count_shared_publications
from pathlib import Path
import os,io,base64
import pandas as pd
from graphviz import Graph


home_dir = str(Path.home())
PROJECT_DIR = os.path.join(home_dir, ".project", "group5")
DATA_DIR = os.path.join(PROJECT_DIR,"data")

app = Flask(__name__)
app.secret_key = "someSecretKey"

@app.route("/")
def home():
    return render_template('home.html',msg='',li=[],next_msg='',info='')

@app.route("/about")
def upload():
    return render_template('about.html')

@app.route('/info', methods=['GET', 'POST'])
def display_data():#home_info():
    if request.method == 'POST':
        seq=request.form.get("author")
        if type(seq) == str: #isinstance(seq, str):
            data = AuthorData(seq)   #make changes here
            li = data.get_list_of_coauthors()

            if type(li) == list:
                li = data.get_list_of_coauthors()
                info="The total number of Co-author's for {} is : {} . The co-author list is given below:".format(seq,len(li))
                co_msg="Would you like to see Co-author Network?"
                return render_template('home.html', msg=info, res=sorted(li),next_msg=co_msg)

            else:
                return render_template('home.html', msg="The author {} does not have any associated co-author's".format(seq), res="",next_msg='')

        else:
            return render_template('home.html',msg="Please enter it in the specified format",res="",next_msg='')   #I don't understand why this isn't working



@app.route('/choose', methods=['GET', 'POST'])
def coauthor_map():
    if request.method == 'POST':
        # print(request.form['check'])
        selected = request.form.getlist('check')
        any_selected = bool(selected)
        if any_selected == True:
            # some check boxes are checked
            if len(selected)==2:
                # display the graph
                dic=count_shared_publications(selected)
                # testobj = AuthorNetwork(dic)
                # graph=testobj.visualize_as_string()
                # test_data = {("Ilya", "Marlo"): 3, ("Pragya", "Dhruv"): 4, ("Marlo", "Dhruv"): 2, ("Ilya", "Dhruv"): 1,
                #              ("Pragya", "Ilya"): 7}
                testobj = AuthorNetwork(dic)  # ,enable_annotations=False)
                graph=testobj.build_network()
                # print(type(graph))
                # testobj.save_graph("pdf", view=True)
                # print("works")
                chart_output = graph.pipe(format='png')
                chart_output = base64.b64encode(chart_output).decode('utf-8')
                # print("should have printed")
                return render_template('plot.html', op=chart_output,info="The Co-author Network ")
            else:
                # display msg to select only 2 checkboxes
                return render_template('home.html', info="Please select only any 2 co-author's")
        else:
            # no checkbox is checked, give error to check again
            return render_template('home.html', info="Please select any 2 co-author's before clicking on submit.")



if __name__ == '__main__':
    Link = 'http://127.0.0.1:5000'
    print(f"{Link} ")
    app.run(debug=True)



#data = AuthorData('William Joyce')
# Svjetlana Miocinovic
# Andrew Miller