from flask import Flask, flash, request, redirect, url_for, render_template, session
from authormaps.authorinfo import AuthorData
from authormaps.author_network import AuthorNetwork
from authormaps.networkutils import count_shared_publications
import os, io, base64
import urllib

app = Flask(__name__)
app.secret_key = "someSecretKey"


@app.route("/")
def home():
    session.clear()
    return render_template('home.html', msg='', li=[], next_msg='', info='',newinfo='')


@app.route("/about")
def upload():
    return render_template('about.html')


@app.route('/info', methods=['GET', 'POST'])
def display_data():

    try:
        if request.method == 'POST':
            seq = request.form.get("author")

            if type(seq) == str:
                session["author_name"] = seq
                data = AuthorData(seq)
                li = data.get_list_of_coauthors()

                if type(li) == list:
                    li = data.get_list_of_coauthors()
                    info = "The total number of Co-author's for {} is : {} . The co-author list is given below:".format(seq,len(li))
                    co_msg = "Would you like to see Co-author Network?"
                    return render_template('home.html', msg=info, res=sorted(li), next_msg=co_msg)

                else:
                    return render_template('home.html',
                                           msg="The author {} does not have any associated co-author's".format(seq), res="",
                                           next_msg='')

            else:
                return render_template('home.html', msg="Please enter it in the specified format", res="",
                                       next_msg='')  # I don't understand why this isn't working

    except urllib.error.URLError:
        return render_template('home.html', msg="Please enter the name before clicking on Submit", res="",
                               next_msg='')


@app.route('/choose', methods=['GET', 'POST'])
def coauthor_map():

    if request.method == 'POST':
        selected = request.form.getlist('check')
        any_selected = bool(selected)

        if any_selected == True:    # some check boxes are checked

            if len(selected) == 2:   # display the graph

                data = AuthorData(session["author_name"])
                coauthors = data.get_list_of_coauthors()
                shared_publications = count_shared_publications(coauthors)
                new = []
                for i, k in shared_publications.items():
                    new.append({"a1": i[0], "a2": i[1], "count": k})

                testobj = AuthorNetwork(shared_publications,highlight_authors=[selected[0], selected[1]])

                if testobj is not None:

                    session["author1"] = selected[0]
                    session["author2"] = selected[1]
                    session["shared_pub"] = new

                    graph = testobj.graph
                    chart_output = graph.pipe(format='png')
                    chart_output = base64.b64encode(chart_output).decode('utf-8')
                    return render_template('plot.html', op=chart_output,
                                           info="The Co-author Network for {} and {}".format(selected[0], selected[1]))
                else:
                    return render_template('home.html',
                                           info="The co-author's {} and {} do not have any relationship".format(
                                               selected[0], selected[1]))

            else:     # display msg to select only 2 checkboxes

                return render_template('home.html', info="Please select only any 2 co-author's")

        else:    # no checkbox is checked, give error to check again

            return render_template('home.html', info="Please select any 2 co-author's before clicking on submit.")


@app.route('/saveoptions', methods=['GET', 'POST'])
def save_rendered_img():

    try:
        if request.method == 'POST':
            option = request.form['radopt']

            print(option)
            shared_publications = {}

            for i in session["shared_pub"]:
                j = list(i.values())
                shared_publications[(j[0], j[1])] = j[2]

            testobj = AuthorNetwork(shared_publications, highlight_authors=[session["author1"], session["author2"]])

            if option=='jpg':
                print("selected option is jpg")
                testobj.save_graph("jpg")


            elif option=='pdf':
                print("selected option is pdf")
                testobj.save_graph("pdf")
    #
            elif option=='svg':
                print("selected option is svg")
                testobj.save_graph("svg")

            elif option=='png':
                print("selected option is png")
                testobj.save_graph("png")

            else:
                return render_template('plot.html', savemsg="Please select any file format before clicking on submit.")

            print("Saved")
            return render_template('plot.html', savemsg="File has been saved")

    except KeyError:
        return render_template('plot.html', savemsg="Please select any file format before clicking on submit.")







if __name__ == '__main__':
    # Link = 'http://127.0.0.1:5000'
    # print(f"{Link} ")
    app.run(debug=True)

# data = AuthorData('William Joyce')
# Svjetlana Miocinovic
# Andrew Miller

# Nasioudis Andreas
# Simmons Adam M
