from flask import Flask, render_template


home_dir = str(Path.home())
PROJECT_DIR = os.path.join(home_dir, ".ex05", "bschultz")
DATA_DIR = os.path.join(PROJECT_DIR, "data")
UPLOAD_DIR = os.path.join(PROJECT_DIR, "uploads")
os.makedirs(UPLOAD_DIR,exist_ok=True)
TEMP_DIR= os.path.dirname(__file__)
UPLOAD_FOLDER = UPLOAD_DIR
ALLOWED_EXTENSIONS = {'png' ,'jpg' ,'svg' ,'pdf'}

OUTPUT_DIR = os.path.join(PROJECT_DIR,"output")
os.makedirs(OUTPUT_DIR,exist_ok=True)

app = Flask(__name__)
app.secret_key = "someSecretKey"
app.config['MAX_CONTENT_PATH'] = 10 * 1024 * 1024  # Max 10MB

@app.route("/")
def home():
    return render_template('home.html')

@app.route("/about")
def upload():
    return render_template('about.html')

def allowed_file(filename):
    return '.' in filename and \
           filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS

if __name__ == '__main__':
    app.run(debug=True)

