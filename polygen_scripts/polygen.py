from flask import Flask, redirect, url_for, render_template, request, Response
import os

from engine_v2_2 import *
    
session = {}		# dictionary to pass variables around routes


app = Flask(__name__)
app.secret_key = 'ee839687a282e6493054c86e00e86925dc04de931acf4aed'

@app.route("/")
def home():
    return render_template("index.html", content="Testing")


@app.route("/learn")
def learn():
    return render_template("learn_more.html", content="Testing")


@app.route("/ptg", methods=["POST","GET"])
def sequence():
    session['msg'] = None
    session['clr'] = {'sequence_spacers': '#FFFFFF', 'link': '#FFFFFF', 'poltype_input': '#FFFFFF'}
    
    if request.method == "POST":
        if request.form['submit_button'] == 'submit':
            runall_args = {}
            runall_args['poltype_run'] = request.form["poltype_input"]
            session['PTG_name'] = request.form["PTG_name"]
            session['oligo_prefix'] = request.form["oligo_prefix"]
            session['oligo_index'] = request.form["oligo_index"]
            session['PTG_transfer'] = request.form["sequence_spacers"]
            session['bb_ovrhng'] = request.form["bb_ovrhng"]
            session['add_ovrhng'] = request.form["add_ovrhng"]
        
            runall_args['tm_range'] = [int(request.form["min_temp"][:2]), int(request.form["max_temp"][:2])]
            if session['bb_ovrhng']:
                runall_args['bb_overlaps'] = session["bb_ovrhng"].split(';')
            if session['add_ovrhng']:
                runall_args['additional_overhangs'] = session["add_ovrhng"].split(';')
            
            PTG_input = session['PTG_transfer'].split('|')
            PTG_structure = []
            PTG_index=1
            for element in PTG_input:
                element_list=[]
                element_list.append([session['oligo_prefix'] if session['oligo_prefix'] else 'o'][0]+'_'+str(PTG_index))
                PTG_index+=1
                for e in element.split(';'):
                    element_list.append(e)
                PTG_structure.append(element_list)  

            session['out'],ftrs,session['msg'],session['primers'] = runall(PTG_structure, **runall_args)
            if session['msg'] == 'comb_error':
                return render_template("sequence.html", PTG_transfer=session.get('PTG_transfer', None), session=session)
        
            if not session['PTG_name']:
                session['PTG_name'] = 'PTG'
            if not session['oligo_prefix']:
                session['oligo_prefix'] = 'o'
            if not session['oligo_index']:
                session['oligo_index'] = '0'
        
            return render_template("primer_list.html", session=session)
            
        elif request.form['submit_button'] == 'reset':
            session['PTG_name'] = session['oligo_prefix'] = session['oligo_index'] = session['bb_ovrhng'] = session['add_ovrhng'] = session['PTG_transfer'] = ''
            return render_template("sequence.html",PTG_transfer=session.get('PTG_transfer', None), session=session)
            
    else:
        
        return render_template("sequence.html",PTG_transfer=session.get('PTG_transfer', None), session=session)


@app.route("/peg", methods=["POST","GET"])
def peg_generation():
    session['msg'] = None
    session['clr'] = {'sequence': '#FFFFFF', 'edits': '#FFFFFF'}

    if request.method == "POST":
        session['PEG_sequence'] = request.form["sequence"]
        session['PEG_edits'] = request.form["edits"]
        PEG_edits = session['PEG_edits'].split('|')
        PEG_mode = request.form["mode"]

        pegs_list = []
        for edt in PEG_edits:
            pegs_list.append(edt.split(';'))

        edits_list = pegbldr(session['PEG_sequence'], pegs_list, PEG_mode)

        session['PTG_transfer'] = ''
        for c,peg in enumerate(edits_list):
            session['PTG_transfer'] += str(peg[1]) + ';' +  str(peg[2])
            if c < len(edits_list)-1:
                session['PTG_transfer'] += '|'

        return redirect(url_for('sequence'))
    else:
        return render_template("peg_generation.html", PEG_transfer=session.get('PEG_sequence', None), session=session)
      
        
@app.route("/primer_list", methods=["POST","GET"])
def serve_primers():
    csv = 'List of oligos:\n'
    csv += 'oligo_ID,sequence\n'
    positions = len(session['oligo_index'])
    collapsed_index = int(session['oligo_index'])
    oligo_ids = []
    for primer in session['primers']:
        oligo_ids.append(session['oligo_prefix']+format(collapsed_index,'0'+str(positions)))
        csv += session['oligo_prefix']+format(collapsed_index,'0'+str(positions))+','+primer+'\n'
        collapsed_index += 1
    
    csv += '\nTable of fragments:\n'
    csv += 'fragment_id,fragment_type,forward_primer,Tm_forw,reverse_primer,Tm_rev\n'
    for c,fragment in enumerate(session['out']):
        csv += session['PTG_name']+'_f'+str(c)+','+fragment.type+','+oligo_ids[c*2]+','+str(np.round(fragment.primer_forward_tm,1))+','+oligo_ids[c*2+1]+','+str(np.round(fragment.primer_reverse_tm, 1))+'\n'
    return Response(
        csv,
        mimetype="text/csv",
        headers={"Content-disposition":
                 "attachment; filename=%s_oligos.csv"%session['PTG_name']})
                 

@app.route("/impressum")
def impress():
    return render_template("impressum.html", content="Testing")


@app.route("/success")
def success():
    return "Success"


@app.errorhandler(InvalidUsage)
def handle_invalid_usage(error):
    session['msg'] = error.to_dict()['message']
    session['clr'][error.to_dict()['box']] = '#FF0000'
    return render_template(error.to_dict()['pge'], PEG_transfer=session.get('PEG_sequence'), PTG_transfer=session.get('PTG_transfer', None), session=session)


if __name__=='__main__':
    app.run(debug=True, host='0.0.0.0', port=5000) ##host 0.0.0.0 is important for docker container
