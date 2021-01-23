from flask import Flask, flash, request, redirect, url_for, render_template, session
from authormaps.authorinfo import AuthorData
from pathlib import Path
import os
import pandas as pd


home_dir = str(Path.home())
PROJECT_DIR = os.path.join(home_dir, ".project", "group5")
DATA_DIR = os.path.join(PROJECT_DIR,"data")
ALLOWED_EXTENSIONS = {'png', 'jpg', 'svg', 'pdf'}

app = Flask(__name__)
app.secret_key = "someSecretKey"

@app.route("/")
def home():
    return render_template('home.html',msg='',li=[],next_msg='')

@app.route("/about")
def upload():
    return render_template('about.html')

@app.route('/info', methods=['GET', 'POST'])
def display_data():#home_info():
    if request.method == 'POST':
        seq=request.form.get("author")
        if type(seq) == str: #isinstance(seq, str):
            data = AuthorData(seq,5)   #make changes here
            li = data.get_list_of_coauthors()

            if type(li) == list:
                li = data.get_list_of_coauthors()
                co_msg="Would you like to see Co-author Network?"
                return render_template('home.html', msg="The Co-author list for {} is:\n".format(seq), res=sorted(li),next_msg=co_msg)

            else:
                return render_template('home.html', msg="The author {} does not have any associated co-author's".format(seq), res="",next_msg='')#df.to_html())

        else:
            return render_template('home.html',msg="Please enter it in the specified format",res="",next_msg='')   #I don't understand why this isn't working
#
# @app.route('/choice', methods=['GET', 'POST'])
# def display_data():#home_info():
#     if request.method == 'POST':
#         author_names=request.form.get("author")




def allowed_file(filename):
    return '.' in filename and \
           filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS

if __name__ == '__main__':
    app.run(debug=True)



#data = AuthorData('William Joyce')
# Svjetlana Miocinovic