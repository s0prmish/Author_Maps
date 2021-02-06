import os
import io
import base64
import graphviz
from pathlib import Path
from authormaps.author_network import AuthorNetwork
from authormaps.authorinfo import AuthorData
from authormaps.networkutils import count_shared_publications
from flask import Flask, request, render_template, session, send_file

home_dir = str(Path.home())
PROJECT_DIR = os.path.join(home_dir, ".project", "group5")
DATA_DIR = os.path.join(PROJECT_DIR, "data")

app = Flask(__name__)
app.secret_key = "someSecretKey"

MIMETYPE = {
    'jpg': 'image/jpeg',
    'pdf': 'application/pdf',
    'svg': 'image/svg+xml',
    'png': 'image/png'
}


@app.route("/")
def home():
    session.clear()
    return render_template('home.html', msg='', li=[], next_msg='', info='',newinfo='')


@app.route("/about")
def upload():
    return render_template('about.html')


@app.route('/info', methods=['GET', 'POST'])
def display_data():

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

                testobj = AuthorNetwork(shared_publications,highlight_authors=[selected[0], selected[1]])

                if testobj is not None:
                    graph = testobj.graph
                    session["graph_source"] = graph.source
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

    if request.method == 'POST':
        option = request.form['radopt']
        graph = graphviz.Source(session["graph_source"])
        return send_file(
            io.BytesIO(graph.pipe(format=option)),
            mimetype=MIMETYPE[option],
            as_attachment=True,
            attachment_filename=f'AuthorGraph.{option}'
        )


if __name__ == '__main__':
    Link = 'http://127.0.0.1:5000'
    print(f"{Link} ")
    app.run(debug=True)

# data = AuthorData('William Joyce')
# Svjetlana Miocinovic
# Andrew Miller

# Nasioudis Andreas
# Simmons Adam M
