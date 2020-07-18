from flask import Flask, redirect, url_for, render_template, request
import os

from engine_v2 import *


app = Flask(__name__)

@app.route("/")
def home():
    return render_template("index.html", content="Testing")

@app.route("/learn")
def learn():
    return render_template("learn_more.html", content="Testing")


@app.route("/ptg", methods=["POST","GET"])
def sequence():
    if request.method == "POST":
        PTG_input = request.form["sequence_spacers"]
        PTG_name = request.form["PTG_name"]
        PTG_input = PTG_input.split('|')
        PTG_structure = []
        PTG_index=1
        for element in PTG_input:
            element_list=[]
            element_list.append(PTG_name+'_'+str(PTG_index))
            PTG_index+=1
            for e in element.split(';'):
                element_list.append(e)

            PTG_structure.append(element_list)
        print(PTG_structure)
        out,ftrs = runall(PTG_structure)
        return render_template("primer_list.html", out=out)

    else:
        return render_template("sequence.html")

@app.route("/peg", methods=["POST","GET"])
def peg_generation():
    if request.method == "POST":
        PEG_sequence = request.form["sequence"]
        PEG_edit = request.form["edits"]

        peg_result = pegbldr(PEG_sequence, PEG_edit)

        return peg_result[0][2]
    else:
        return render_template("peg_generation.html")

@app.route("/sucess")
def success():
    return "Sucess"


if __name__=='__main__':
    app.run(debug=True, host='0.0.0.0', port=5000) ##host 0.0.0.0 is important for docker container
