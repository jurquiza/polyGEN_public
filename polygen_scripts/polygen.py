from flask import Flask, redirect, url_for, render_template, request
import os

from engine_v2_2 import *

session = {}  #here you can store variables that will be passed around routes 
session['msg'] = None

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
        runall_args = {}
        PTG_input = request.form["sequence_spacers"]
        session['PTG_transfer'] = PTG_input
        PTG_name = request.form["PTG_name"]
        session['PTG_oligo'] = request.form["oligo_prefix"]
        runall_args['tm_range'] = [int(request.form['min_temp'][:2]), int(request.form['max_temp'][:2])]
        if request.form['max_len']:
            runall_args['max_ann_len'] = int(request.form['max_len'])
        if request.form['bb_ovrhng']:
            runall_args['bb_overlaps'] = request.form['bb_ovrhng'].split(';')
        if request.form['add_ovrhng']:
            runall_args['additional_overhangs'] = request.form['add_ovrhng'].split(';')
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

        out,ftrs,session['msg'] = runall(PTG_structure, **runall_args)
        if session['msg'] == 'comb_error':
            return render_template("sequence.html", PTG_transfer=session.get('PTG_transfer', None), session=session)
        
        return render_template("primer_list.html", out=out, session=session)

    else:

        return render_template("sequence.html",PTG_transfer=session.get('PTG_transfer', None))

@app.route("/peg", methods=["POST","GET"])
def peg_generation():
    if request.method == "POST":
        PEG_sequence = request.form["sequence"]
        PEG_edits = request.form["edits"]
        PEG_edits = PEG_edits.split('|')

        pegs_list = []
        for edt in PEG_edits:
            pegs_list.append(edt.split(';'))

        edits_list = pegbldr(PEG_sequence, pegs_list)

        session['PTG_transfer'] = ''
        for c,peg in enumerate(edits_list):
            session['PTG_transfer'] += str(peg[1]) + ';' +  str(peg[2])
            if c < len(edits_list)-1:
                session['PTG_transfer'] += '|'

        return redirect(url_for('sequence'))
    else:
        return render_template("peg_generation.html")

@app.route("/sucess")
def success():
    return "Sucess"


if __name__=='__main__':
    app.run(debug=True, host='0.0.0.0', port=5000) ##host 0.0.0.0 is important for docker container
