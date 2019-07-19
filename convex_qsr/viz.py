import os
import webbrowser

from flask import Flask, send_from_directory, render_template


def show_viz(data_dir):
    app = Flask(__name__)
    port = 16180
    abs_data_dir = os.path.abspath(data_dir)

    @app.route('/api/<path:filename>')
    def send_file(filename):
        return send_from_directory(abs_data_dir, filename)

    @app.route('/')
    def index():
        return render_template('index.html')

    webbrowser.open('http://localhost:' + str(port), new=2)
    print('Serving visualization on port', port, '...')
    print('Using data from directory', abs_data_dir, '...')
    app.run(host='0.0.0.0', port=port)
