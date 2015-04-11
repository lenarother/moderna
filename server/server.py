
"""

Primitive Flask server for use locally and within Docker containers.

usage: python server.py

"""

from flask import Flask, request, session, g, redirect, url_for, \
     abort, render_template, flash
from werkzeug import secure_filename
import os
import moderna


DEBUG = True
SECRET_KEY = 'development key'
USERNAME = 'admin'
PASSWORD = 'default'
UPLOAD_FOLDER = 'uploads'
ALLOWED_EXTENSIONS = set(['txt', 'pdf', 'png', 'jpg', 'jpeg', 'gif'])


# create our little application :)
app = Flask(__name__)
app.config.from_object(__name__)
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER

# app.config.from_envvar('FLASKR_SETTINGS', silent=True)

@app.route('/', methods=['GET', 'POST'])
def main():
    log = ''
    model = None
    if request.method == 'POST':
        template = request.files['template']
        alignment = request.files['alignment']
        chain = 'A'
        if template:
            template_fn = secure_filename(template.filename)
            template_fn = os.path.join(app.config['UPLOAD_FOLDER'], template_fn)
            template.save(template_fn)
        if alignment:
            alignment_fn = secure_filename(alignment.filename)
            alignment_fn = os.path.join(app.config['UPLOAD_FOLDER'], alignment_fn)
            alignment.save(alignment_fn)
        if template and alignment:
            # run moderna
            t = moderna.load_template(template_fn, chain)
            a = moderna.load_alignment(alignment_fn)
            m = moderna.create_model(t, a)
            moderna.write_model(m, 'static/model.pdb')
            moderna.write_logfile()
            log = open('moderna.log').read()
        elif template:
            t = moderna.load_template(template_fn, chain)
            moderna.get_sequence(t)
            moderna.examine_structure(t)
            moderna.write_logfile()
            log = open('moderna.log').read()
            # return redirect(url_for('uploaded_file', filename=filename))
        if model:
            model = 'static/model.pdb'
    return render_template('main.html', log=log, model=model)
    


if __name__ == '__main__':
    app.run(host='0.0.0.0')
